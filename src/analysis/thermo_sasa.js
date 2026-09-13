/**
 * thermo_sasa.js — Stage-3 real ΔSASA burial for the ΔH/ΔS thermo path.
 *
 * Pain (open debt from Stage-6 calibration §13): the solvent term was dead —
 * as-run ΔSASA = 0 (no contact-count series ever passed) → −TΔS_solv 0.00
 * with a ±50% swing of ±0.00. The LCPO module
 * (src/physics/solvation/lcpo_sasa.js) existed but was never wired into the
 * holo−apo thermo legs.
 *
 * Method (stated, additive, zero deps beyond lcpo_sasa.js):
 *   ΔSASA = Δ_lig + Δ_prot, with
 *   - Δ_lig (ligand-side burial, the dominant term): per holo frame, the free
 *     ligand LCPO areas (ligand-only subsystem — in-regime for the small
 *     molecule) multiplied by the cross-burial survival fraction from protein
 *     neighbours. For ligand atom a with free LCPO area A_free,a and cross
 *     caps C_pa from protein neighbours (LCPO capAreaDeriv, independent-
 *     overlap deconvolution):
 *       buried_a = A_free,a · (1 − Π_p max(0, 1 − C_pa / S_a)),
 *     Δ_lig = ⟨Σ_a buried_a⟩ over subsampled holo frames.
 *   - Δ_prot (protein reorganization): ⟨S_prot⟩ stripped-protein full-LCPO on
 *     subsampled apo frames minus the same on holo frames (same regime both
 *     legs, self-burial cancels; typically small, ±tens of Å²).
 *   SEs are SEMs over the evaluated frame means; dsasa SE = √(SE_lig²+SE_prot²)
 *   (legs treated as independent — documented approximation).
 *
 * Why NOT full-complex LCPO totals (measured, 4W52 CG BNZ, seeds 101/1101):
 * atomic LCPO P3/P4 triple/crowding fits are OUT of regime at Cα-bead scale
 * (bead Rp 3.4 Å vs 3.8 Å spacing — everything overlaps everything):
 *   - bound-ligand LCPO areas (419.8) EXCEED free-ligand areas (340.8) by
 *     79 Å² (anti-burial inversion — protein neighbours inflate P3/P4 terms);
 *   - full SEP (apo-stripped + free − complex) = −49.0 ± 3.9 Å² (complex
 *     LARGER than its parts — not subadditive, unphysical);
 *   - literal apo−holo totals = −15.0 Å² (same-n legs; apo keeps the pinned
 *     ligand per the existing protocol, so no free reference is sampled).
 * The cross-burial route above is monotonic by construction (caps ≥ 0,
 * survival fractions ∈ [0,1]) and lands at the benzene burial scale:
 * Δ_lig 158.9 ± 2.6 + Δ_prot 6.4 ± 3.6 = 165.3 ± 4.4 Å² (stride-10, 100
 * frames/leg; stride-5 gives 159.0 — subsampling-insensitive).
 *
 * Cost (documented): per evaluated holo frame one ligand-only LCPO
 * (O(nL²)) + O(nProt·nSel) cross distances; per stripped-protein frame one
 * full LCPO (O(nProt²), forces off). CG 4W52 (164+6) ≈ 0.3 ms/frame →
 * 100 frames/leg is sub-second. Default stride 10 keeps even 10k-frame
 * recordings tractable; the chunked variant keeps the UI responsive.
 *
 * Units: Å for areas/distances; kcal/mol conversions use SASA_GAMMA in
 * thermodynamics.js (this module returns areas only).
 */

import {
  lcpoSasa, capAreaDeriv, vdwRadiusFor, PROBE_RADIUS,
} from "../physics/solvation/lcpo_sasa.js";

/** Default frame subsample stride for SASA evaluation (every Nth frame). */
export const THERMO_SASA_STRIDE = 10;

/** Default evaluated frames per async chunk (keeps slices < ~50 ms). */
export const THERMO_SASA_CHUNK = 25;

/**
 * Extended CG bead radius (contact radius + water probe) for cross caps.
 * @param {object} [ff] live CG ForceField (uses .sigmaRep, default 4.0 Å)
 * @returns {number} Rp in Å (default 3.4)
 */
export function cgBeadExtendedRadius(ff) {
  return (ff?.sigmaRep ?? 4.0) / 2 + PROBE_RADIUS;
}

/**
 * Element symbols of the SELECTED ligand subset in concatenation order.
 * @param {Array} ligands full molecule list (parseLigands order)
 * @param {number[]|null} molIdx resolved molecule indices (null ⇒ all)
 * @returns {string[]} uppercase element per selected ligand atom
 */
export function selectedLigandElements(ligands, molIdx) {
  const mols = Array.isArray(ligands) ? ligands : [];
  const keep = molIdx == null ? mols.map((_, i) => i) : [...molIdx];
  const out = [];
  for (const i of keep) {
    for (const at of mols[i]?.atoms ?? []) out.push(String(at?.element ?? "C").toUpperCase());
  }
  return out;
}

/**
 * Extended radii for ligand atoms (vdW + probe, per element).
 * @param {string[]} elements uppercase element symbols
 * @returns {Float64Array} Rp per atom (Å)
 */
export function ligandExtendedRadii(elements) {
  const out = new Float64Array(elements?.length ?? 0);
  for (let a = 0; a < out.length; a++) out[a] = vdwRadiusFor(elements[a]) + PROBE_RADIUS;
  return out;
}

/**
 * Per-atom isolated-sphere areas 4πRp².
 * @param {Float64Array|number[]} Rp extended radii
 * @returns {Float64Array}
 */
export function isolatedAreas(Rp) {
  const out = new Float64Array(Rp.length);
  for (let i = 0; i < Rp.length; i++) out[i] = 4 * Math.PI * Rp[i] * Rp[i];
  return out;
}

/**
 * Ligand-side burial for ONE holo frame (free LCPO × cross survival).
 * Deterministic: no RNG, pure function of coordinates.
 * @param {ArrayLike} frame flat coords (protein block + ligand block)
 * @param {object} o
 * @param {number} o.nProt protein particle count (protein occupies 0..nProt−1)
 * @param {number} [o.ligStart] ligand-block offset (default o.nProt; heavy
 *   mode passes ff.ligandStart since hetero atoms sit in between)
 * @param {number[]|null} [o.holoSel] ligand-block-relative indices of the
 *   SELECTED subset (null ⇒ whole block). Order must match o.selElements.
 * @param {string[]} o.selElements element per selected ligand atom
 * @param {Float64Array} [o.RpL] precomputed ligand extended radii (else derived)
 * @param {Float64Array} [o.SL] precomputed isolated areas (else derived)
 * @param {Float64Array|number[]|null} [o.protRp] per-protein-atom extended
 *   radii (uniform number ⇒ broadcast; null ⇒ element-derived C default —
 *   prefer the explicit CG bead radius from cgBeadExtendedRadius)
 * @returns {{burial:number, free:number}} buried + free ligand SASA (Å²)
 */
export function frameLigandBurial(frame, o) {
  const nProt = o.nProt;
  const ligStart = o.ligStart ?? nProt;
  const selEls = o.selElements ?? [];
  const nSel = selEls.length;
  if (!nSel || !frame?.length) return { burial: 0, free: 0 };
  const RpL = o.RpL ?? ligandExtendedRadii(selEls);
  const SL = o.SL ?? isolatedAreas(RpL);
  const holoSel = o.holoSel ?? null;
  // selected ligand coordinates (strided read, no full-frame copy)
  const lp = new Float64Array(nSel * 3);
  for (let k = 0; k < nSel; k++) {
    const a = holoSel ? holoSel[k] : k;
    lp[3 * k] = frame[3 * (ligStart + a)];
    lp[3 * k + 1] = frame[3 * (ligStart + a) + 1];
    lp[3 * k + 2] = frame[3 * (ligStart + a) + 2];
  }
  const free = lcpoSasa(lp, [...selEls], { includeForces: false });
  const protRp = o.protRp ?? null;
  const uniformRp = typeof protRp === "number" ? protRp : null;
  let burial = 0;
  for (let k = 0; k < nSel; k++) {
    let prod = 1;
    for (let i = 0; i < nProt; i++) {
      const rpi = uniformRp ?? (protRp ? protRp[i] : vdwRadiusFor("C") + PROBE_RADIUS);
      const jx = lp[3 * k], jy = lp[3 * k + 1], jz = lp[3 * k + 2];
      const d = Math.hypot(frame[3 * i] - jx, frame[3 * i + 1] - jy, frame[3 * i + 2] - jz);
      const c = capAreaDeriv(RpL[k], rpi, d).A;
      if (c > 0) prod *= Math.max(0, 1 - c / Math.max(1e-9, SL[k]));
      if (prod === 0) break;
    }
    burial += free.areas[k] * (1 - prod);
  }
  return { burial, free: free.total };
}

/**
 * Stripped-protein total LCPO for one frame (protein block only).
 * @param {ArrayLike} frame flat coords (protein first)
 * @param {number} nProt protein particle count
 * @param {object} [o]
 * @param {Float64Array|number[]|null} [o.protRadii] explicit extended radii
 *   (uniform CG bead radius preferred; null ⇒ element-derived)
 * @param {string[]|null} [o.protElements] per-atom elements (default all C)
 * @returns {number} total SASA (Å²)
 */
export function frameProteinSasa(frame, nProt, o = {}) {
  if (!nProt || !frame?.length) return 0;
  const sub = frame.slice(0, nProt * 3);
  const els = o.protElements ?? new Array(nProt).fill("C");
  const opts = { includeForces: false };
  if (o.protRadii) opts.radii = o.protRadii;
  return lcpoSasa(sub, els, opts).total;
}

function semOf(series) {
  const n = series.length;
  if (n < 2) return 0;
  const m = series.reduce((s, v) => s + v, 0) / n;
  const sd = Math.sqrt(series.reduce((s, v) => s + (v - m) ** 2, 0) / (n - 1));
  return sd / Math.sqrt(n);
}

function meanOf(series) {
  if (!series.length) return 0;
  return series.reduce((s, v) => s + v, 0) / series.length;
}

/**
 * Real ΔSASA burial from recorded holo + relaxed apo frames (sync).
 * Subsamples every `stride`-th frame (documented cost above). Deterministic:
 * identical frames give bit-identical numbers (no RNG anywhere).
 * @param {ArrayLike[]} holoFrames recorded holo frames (protein + ligand block)
 * @param {ArrayLike[]} apoFrames relaxed apo frames (protein block required;
 *   any ligand block present is ignored — apo legs retain the pinned ligand
 *   per the existing protocol, so only the protein block is read)
 * @param {object} o
 * @param {number} o.nProt protein particle count
 * @param {number} [o.ligStart] ligand-block offset in holo frames
 *   (default nProt; heavy mode: ff.ligandStart)
 * @param {number[]|null} [o.holoSel] selected-subset ligand indices
 *   (Stage-2 picker; null ⇒ whole block)
 * @param {string[]} [o.selElements] element per selected ligand atom
 *   (default all C — BNZ-safe; EPE callers should pass real elements)
 * @param {number|null} [o.protRadius] uniform protein extended radius
 *   (CG: cgBeadExtendedRadius(ff); heavy: null + protElements for
 *   element-derived radii)
 * @param {string[]|null} [o.protElements] per-protein-atom elements
 * @param {number} [o.stride] subsample stride (default THERMO_SASA_STRIDE)
 * @returns {{dsasa:number, se:number, dLig:number, dLigSE:number,
 *   dProt:number, dProtSE:number, ligFreeMean:number,
 *   stride:number, nHoloEval:number, nApoEval:number,
 *   nHoloFrames:number, nApoFrames:number, method:string}}
 */
export function sasaBurial(holoFrames, apoFrames, o = {}) {
  const nProt = o.nProt ?? 0;
  const stride = Math.max(1, Math.round(o.stride ?? THERMO_SASA_STRIDE));
  const selEls = o.selElements ?? [];
  const RpL = ligandExtendedRadii(selEls);
  const SL = isolatedAreas(RpL);
  const protRp = ("protRadius" in o)
    ? o.protRadius
    : (o.protRadii ?? null);
  const protRadiiForLCPO = (typeof protRp === "number")
    ? new Float64Array(nProt).fill(protRp)
    : (o.protRadii ?? null);
  const shared = {
    nProt, ligStart: o.ligStart ?? nProt, holoSel: o.holoSel ?? null,
    selElements: selEls, RpL, SL,
    protRp: (typeof protRp === "number") ? protRp : (o.protRadii ?? null),
  };
  const ligSeries = [], freeSeries = [];
  const hProtSeries = [], aProtSeries = [];
  const H = holoFrames ?? [], A = apoFrames ?? [];
  for (let f = 0; f < H.length; f += stride) {
    if (!H[f]?.length) continue;
    const r = frameLigandBurial(H[f], shared);
    ligSeries.push(r.burial);
    freeSeries.push(r.free);
    if (nProt) hProtSeries.push(frameProteinSasa(H[f], nProt, { protRadii: protRadiiForLCPO, protElements: o.protElements ?? null }));
  }
  for (let f = 0; f < A.length; f += stride) {
    if (!A[f]?.length) continue;
    if (nProt) aProtSeries.push(frameProteinSasa(A[f], nProt, { protRadii: protRadiiForLCPO, protElements: o.protElements ?? null }));
  }
  const dLig = meanOf(ligSeries), dLigSE = semOf(ligSeries);
  const dProt = meanOf(aProtSeries) - meanOf(hProtSeries);
  const dProtSE = Math.sqrt(semOf(aProtSeries) ** 2 + semOf(hProtSeries) ** 2);
  const dsasa = dLig + dProt;
  return {
    dsasa, se: Math.sqrt(dLigSE ** 2 + dProtSE ** 2),
    dLig, dLigSE, dProt, dProtSE,
    ligFreeMean: meanOf(freeSeries),
    stride,
    nHoloEval: ligSeries.length, nApoEval: aProtSeries.length,
    nHoloFrames: H.length, nApoFrames: A.length,
    method: "lcpo-cross-burial",
  };
}

/**
 * Chunked async variant of sasaBurial (Stage-5 pattern: setTimeout slices,
 * progress callback, UI stays responsive). Resolves the identical object
 * shape as sasaBurial (same evaluation order ⇒ same frame means; the final
 * reductions run once over the full series either way).
 * @param {ArrayLike[]} holoFrames
 * @param {ArrayLike[]} apoFrames
 * @param {object} [o] same as sasaBurial plus:
 * @param {number} [o.chunkFrames] evaluated frames per slice (default 25)
 * @param {(done:number, total:number) => void} [o.onProgress]
 * @returns {Promise<ReturnType<sasaBurial>>}
 */
export function sasaBurialChunked(holoFrames, apoFrames, o = {}) {
  const chunk = Math.max(1, Math.round(o.chunkFrames ?? THERMO_SASA_CHUNK));
  const { chunkFrames: _cf, onProgress: _op, ...rest } = o;
  const onProgress = o.onProgress ?? null;
  return new Promise((resolve) => {
    // Evaluate index lists up-front so chunking cannot change the series.
    const H = holoFrames ?? [], A = apoFrames ?? [];
    const stride = Math.max(1, Math.round(o.stride ?? THERMO_SASA_STRIDE));
    const hIdx = [];
    for (let f = 0; f < H.length; f += stride) hIdx.push(f);
    const aIdx = [];
    for (let f = 0; f < A.length; f += stride) aIdx.push(f);
    const total = hIdx.length + aIdx.length;
    const nProt = o.nProt ?? 0;
    const selEls = o.selElements ?? [];
    const RpL = ligandExtendedRadii(selEls);
    const SL = isolatedAreas(RpL);
    const protRpOpt = ("protRadius" in o) ? o.protRadius : (o.protRadii ?? null);
    const protRadiiForLCPO = (typeof protRpOpt === "number")
      ? new Float64Array(nProt).fill(protRpOpt)
      : (o.protRadii ?? null);
    const shared = {
      nProt, ligStart: o.ligStart ?? nProt, holoSel: o.holoSel ?? null,
      selElements: selEls, RpL, SL,
      protRp: (typeof protRpOpt === "number") ? protRpOpt : (o.protRadii ?? null),
    };
    const ligSeries = [], freeSeries = [], hProtSeries = [], aProtSeries = [];
    let cursor = 0;
    const jobs = [
      ...hIdx.map((f) => ({ leg: "h", f })),
      ...aIdx.map((f) => ({ leg: "a", f })),
    ];
    const step = () => {
      const n = Math.min(chunk, jobs.length - cursor);
      for (let k = 0; k < n; k++) {
        const { leg, f } = jobs[cursor++];
        if (leg === "h") {
          if (!H[f]?.length) continue;
          const r = frameLigandBurial(H[f], shared);
          ligSeries.push(r.burial);
          freeSeries.push(r.free);
          if (nProt) hProtSeries.push(frameProteinSasa(H[f], nProt, { protRadii: protRadiiForLCPO, protElements: o.protElements ?? null }));
        } else {
          if (!A[f]?.length) continue;
          if (nProt) aProtSeries.push(frameProteinSasa(A[f], nProt, { protRadii: protRadiiForLCPO, protElements: o.protElements ?? null }));
        }
      }
      try { if (onProgress) onProgress(Math.min(cursor, jobs.length), jobs.length); } catch (_) {}
      if (cursor < jobs.length) { setTimeout(step, 0); return; }
      const dLig = meanOf(ligSeries), dLigSE = semOf(ligSeries);
      const dProt = meanOf(aProtSeries) - meanOf(hProtSeries);
      const dProtSE = Math.sqrt(semOf(aProtSeries) ** 2 + semOf(hProtSeries) ** 2);
      void total;
      resolve({
        dsasa: dLig + dProt, se: Math.sqrt(dLigSE ** 2 + dProtSE ** 2),
        dLig, dLigSE, dProt, dProtSE,
        ligFreeMean: meanOf(freeSeries),
        stride,
        nHoloEval: ligSeries.length, nApoEval: aProtSeries.length,
        nHoloFrames: H.length, nApoFrames: A.length,
        method: "lcpo-cross-burial",
      });
    };
    if (!jobs.length) {
      resolve({
        dsasa: 0, se: 0, dLig: 0, dLigSE: 0, dProt: 0, dProtSE: 0,
        ligFreeMean: 0, stride, nHoloEval: 0, nApoEval: 0,
        nHoloFrames: H.length, nApoFrames: A.length, method: "lcpo-cross-burial",
      });
      return;
    }
    setTimeout(step, 0);
  });
}
