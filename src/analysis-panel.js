/**
 * analysis-panel.js — panel 6: trajectory analysis report + PMF CSV download
 * (item 5 modularization: moved verbatim from main.js; only the imports
 * changed — `ui` / `state` / `recorder` now come from ui.js).
 *
 * Phase 4: opt-in translational workflows (alanine scanning, DCCM heatmap,
 * SMD unbinding ensemble, cryptic-pocket tracking) as buttons that degrade
 * gracefully when no system / trajectory / ligand is present. The original
 * Analyze flow is untouched.
 */

import { analyzeTrajectory, pmfCsv } from "./analysis.js?v=10";
import { downloadText } from "./recorder.js?v=10";
import { ui, state, viewer, recorder } from "./ui.js?v=10";
import { scanPocket, pocketResidues, formatMutationTable } from "./analysis/alanine_scanning.js?v=10";
import { computeDCCM, renderDCCMHeatmap, topCorrelations, dccmPick, highlightCorrelatedPair } from "./analysis/dccm.js?v=10";
import { trackPocketVolume, detectCryptic } from "./analysis/cryptic_pockets.js?v=10";
import { runPullingEnsemble, jarzynskiFreeEnergy, koffSurrogate } from "./analysis/unbinding_smd.js?v=10";
import { computeThermodynamics, formatThermoTable } from "./analysis/thermodynamics.js?v=10";
import {
  THERMO_LIG_AUTO, THERMO_POCKET_RCUT,
  resolveThermoLigand, selectedLigandAtomIndices, selectedLigandCom,
  pocketFromCom, sliceSelectedPositions, refreshThermoLigOptions,
} from "./analysis/thermo_ligand.js?v=10";
import {
  THERMO_SASA_STRIDE, sasaBurialChunked, cgBeadExtendedRadius,
  selectedLigandElements,
} from "./analysis/thermo_sasa.js?v=10";
import { ForceField } from "./forcefield.js?v=10";
import { LangevinIntegrator } from "./integrator.js?v=10";
import { findRotatableBonds } from "./analysis/rotbonds.js?v=10";

/**
 * Count rotatable bonds in the selected thermo ligand subset (Stage-5 display
 * only — concatenated subset graph, same rule as the headless
 * scripts/validate_flexlig.mjs EPE path; never throws, headless-safe).
 * EPE → 4, BNZ → 0, all (BNZ+EPE, disconnected) → 4.
 * @param {Array} subset resolved ligand subset (parseLigands molecules)
 * @returns {number} rotatable-bond count (0 on any failure)
 */
function selectedRotatableCount(subset) {
  try {
    const atoms = [];
    const bonds = [];
    let off = 0;
    for (const mol of subset ?? []) {
      for (const a of mol?.atoms ?? []) atoms.push(a);
      for (const b of mol?.bonds ?? []) bonds.push([b[0] + off, b[1] + off]);
      off += mol?.atoms?.length ?? 0;
    }
    if (!atoms.length || !bonds.length) return 0;
    return findRotatableBonds(atoms, bonds).count ?? 0;
  } catch (_) { return 0; }
}

if (ui.anaBtn) ui.anaBtn.addEventListener("click", () => {
  if (recorder.count === 0) {
    ui.analysisOut.textContent = "⚠ Nothing recorded yet — press ● Rec, run the sim, then Stop.";
    return;
  }
  ui.analysisOut.textContent = "Analyzing… (RMSIP may take a few seconds)";
  // run off the animation frame so the status text paints first
  setTimeout(() => {
    try {
      const rep = analyzeTrajectory({
        frames: recorder.frames,
        times: recorder.times,
        ref: state.ff.ref,
        nProt: state.ff.nProt,
        n: state.ff.n,
        beads: state.sel.beads,
        ff: state.ff,
        funnel: state.funnel,
      });
      ui.analysisOut.textContent = rep.lines.join("\n");
    } catch (err) {
      ui.analysisOut.textContent = "⚠ " + err.message;
    }
  }, 20);
});

if (ui.anaPmfBtn) ui.anaPmfBtn.addEventListener("click", () => {
  if (!state.funnel) { ui.analysisOut.textContent = "⚠ No active funnel (load a ligand + build the system)."; return; }
  try {
    downloadText(pmfCsv(state.funnel), `pmf_${(state.contactsFn || "enm").replace(/\W+/g, "_")}.csv`);
  } catch (err) {
    ui.analysisOut.textContent = "⚠ " + err.message;
  }
});

/* ------------------------------------------------------------------ */
/*  Phase 4 — translational workflows (opt-in, gracefully degrading)   */
/* ------------------------------------------------------------------ */

/** Current live system as a Phase-4 ScanSystem, or null when not built. */
function currentScanSystem() {
  if (!state.ff || !state.sel) return null;
  return {
    mode: state.heavyMode ? "heavy" : "cg",
    sel: state.sel,
    ff: state.ff,
    par: undefined, // rebuilt from the live ff (parFromCgFF / gamma)
    ligands: state.ligands ?? [],
  };
}

function fmtE(v) {
  return Number.isFinite(v) ? v.toFixed(2) : "NaN";
}

function setDccmCaption(txt) {
  if (ui.dccmCaption) ui.dccmCaption.textContent = txt;
  else {
    const el = typeof document !== "undefined" ? document.getElementById("dccmCaption") : null;
    if (el) el.textContent = txt;
  }
}

let _dccmHasData = false;
let _lastDccmEmpty = 0;

/** Three-state (P2) no-data: actionable empty DCCM with caption + legend key. */
export function drawDccmEmpty(reason) {
  _dccmHasData = false;
  const cv = (ui.dccmCanvas) || (typeof document !== "undefined" ? document.getElementById("dccmCanvas") : null);
  if (!cv) return;
  try {
    const w = cv.clientWidth || 280, h = cv.clientHeight || 180;
    if (w === 0 || h === 0) return;
    const dpr = Math.min(2, (typeof window !== "undefined" && window.devicePixelRatio) || 1);
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
      cv.width = Math.round(w * dpr);
      cv.height = Math.round(h * dpr);
    }
    const ctx = cv.getContext("2d");
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = "#475569";
    ctx.font = "11px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(reason || "Record ≥ 3 frames, then DCCM Heatmap", w / 2, h / 2 - 6);
    ctx.font = "9.5px sans-serif";
    ctx.fillStyle = "#334155";
    ctx.fillText("cross-correlation −1 (blue) … +1 (amber)", w / 2, h / 2 + 10);
  } catch (_) { /* headless: caption only */ }
  setDccmCaption("DCCM: record ≥ 3 frames, then compute. Legend −1…+1.");
}

/** Steady-state (P2): ≤1 Hz empty redraw so the canvas never reads as dead. */
export function dccmTick(now) {
  if (_dccmHasData) return;
  const t = now ?? ((typeof performance !== "undefined" && performance.now()) || Date.now());
  if (t - _lastDccmEmpty < 1000) return;
  _lastDccmEmpty = t;
  drawDccmEmpty();
}

// Paint the empty state once at startup (no-data, before any trajectory).
try { drawDccmEmpty(); } catch (_) {}

// ---- Thermodynamics ΔH/ΔS (Loop-2 S5, Stage-5 async, Stage-2 ligand picker) --
// Holo leg = recorded trajectory (needs ≥30 frames); apo leg = internal
// relaxation of a binding-off clone (same ENM, no ligand coupling).
// Stage-5: the apo relaxation runs chunked (setTimeout slices) so the click
// handler never blocks the UI; progress goes to #thermoCaption, thermoBtn is
// disabled during the run, and a generation counter cancels a stale run on
// rebuild or on a fresh click. Final table render is identical
// (formatThermoTable); defaults/numbers unchanged.
// Stage-2 (ligand picker): pocket COM + binding-energy attribution come from
// the SELECTED ligand only (#thermoLig, default auto → BNZ cavity). The
// Loop-2 record path used BNZ+EPE (all HETATM): surface EPE pulls the
// all-mol COM off the benzene cavity and inflates |ΔH| by −3.18 (row A
// −6.82 vs row B BNZ-only −3.64; calibration_4w52.mjs, BINDING_LOOP2_DONE
// §13). History is NOT rewritten — "all" reproduces the record bit-for-bit;
// new runs default to the buffer-free BNZ cavity. Pocket rule unchanged
// (8 Å). See src/analysis/thermo_ligand.js (BNZ-first fallback rule).
/** Apo-relaxation steps per UI slice (Stage-5: keeps each slice < ~50 ms). */
export const THERMO_CHUNK_STEPS = 150;
let _thermoGen = 0;
/**
 * Cancel any in-flight chunked thermo run (e.g. on system rebuild).
 * Additive Stage-5 hook; safe to call when idle or headless.
 */
export function invalidateThermo() {
  _thermoGen++;
  try { if (ui.thermoBtn) ui.thermoBtn.disabled = false; } catch (_) {}
  try { if (ui.thermoCancelBtn) ui.thermoCancelBtn.disabled = true; } catch (_) {}
}
/**
 * Toggle the thermo run/cancel button pair (Stage-7 cancel UX: Cancel is
 * enabled only while a chunked run is in flight, disabled otherwise).
 * Additive + headless-safe (falls back to getElementById when ui.js captured
 * null, e.g. headless import).
 * @param {boolean} running true while the apo/SASA legs are in flight
 */
function setThermoRunning(running) {
  try { if (ui.thermoBtn) ui.thermoBtn.disabled = !!running; } catch (_) {}
  try {
    const cb = ui.thermoCancelBtn
      ?? (typeof document !== "undefined" ? document.getElementById("thermoCancelBtn") : null);
    if (cb) cb.disabled = !running;
  } catch (_) {}
}
// Idle state at load: Cancel disabled until a run starts (HTML also ships disabled).
try {
  const cb0 = typeof document !== "undefined" ? document.getElementById("thermoCancelBtn") : null;
  if (cb0) cb0.disabled = true;
} catch (_) {}
if (ui.thermoCancelBtn) ui.thermoCancelBtn.addEventListener("click", () => {
  invalidateThermo(); // bump generation → stale stepChunk/SASA continuations return early
  setThermoRunning(false); // belt-and-braces idle state (invalidateThermo already resets)
  setThermoCaption("cancelled — thermo run cancelled by user.");
});
function setThermoCaption(txt) {
  try {
    if (ui.thermoCaption) ui.thermoCaption.textContent = txt;
    else if (typeof document !== "undefined") {
      const el = document.getElementById("thermoCaption");
      if (el) el.textContent = txt;
    }
  } catch (_) {}
}
/**
 * Refresh the thermo ligand picker from the live ligand list (Stage-2).
 * Additive: called on buildSystem (see src/main.js) + lazily on click;
 * in-memory default auto (persist NOT required); headless-safe.
 * @returns {string} effective picker value ("auto" default)
 */
export function refreshThermoLigPicker() {
  try {
    const el = (ui.thermoLig)
      || (typeof document !== "undefined" ? document.getElementById("thermoLig") : null);
    return refreshThermoLigOptions(el, state.ligands ?? [], el?.value ?? THERMO_LIG_AUTO);
  } catch (_) { return THERMO_LIG_AUTO; }
}
try { refreshThermoLigPicker(); } catch (_) {}
/** Read the picker value with a headless-safe fallback (default auto). */
function thermoLigValue() {
  try {
    const v = ui.thermoLig?.value
      ?? (typeof document !== "undefined" ? document.getElementById("thermoLig")?.value : null);
    return String(v ?? THERMO_LIG_AUTO) || THERMO_LIG_AUTO;
  } catch (_) { return THERMO_LIG_AUTO; }
}
if (ui.thermoBtn) ui.thermoBtn.addEventListener("click", () => {
  const myGen = ++_thermoGen;
  try {
    if (!state.ff) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
    if (recorder.frames.length < 30) {
      ui.analysisOut.textContent = `⚠ Need ≥30 recorded holo frames for the Schlitter covariance (have ${recorder.frames.length}). ● Rec, run ~100 ps, Stop, retry.`;
      return;
    }
    const ff = state.ff;
    const nProt = ff.nProt;
    // Stage-2: pocket COM from the SELECTED ligand only (default auto→BNZ).
    // Record path ("all") keeps the old all-mol COM bit-identically.
    let thermoSel = thermoLigValue();
    try { thermoSel = refreshThermoLigPicker() && thermoLigValue(); } catch (_) {}
    const resolved = resolveThermoLigand(state.ligands ?? [], thermoSel);
    const selAtomIdx = selectedLigandAtomIndices(state.ligands ?? [], resolved.molIdx);
    // Stage-5 display only: rotatable-bond count for the selected ligand subset
    // (EPE → 4, BNZ → 0) surfaced in the ligand-note line; no physics change.
    const thermoRotCount = selectedRotatableCount(resolved.subset);
    let lcom = [0, 0, 0];
    if (selAtomIdx.length > 0 && ff.nLigAtoms > 0) {
      lcom = selectedLigandCom(ff.ref, nProt, selAtomIdx);
    } else if (ff.nLigAtoms > 0) {
      for (let a = 0; a < ff.nLigAtoms; a++) {
        lcom[0] += ff.ref[3 * (nProt + a)] / ff.nLigAtoms;
        lcom[1] += ff.ref[3 * (nProt + a) + 1] / ff.nLigAtoms;
        lcom[2] += ff.ref[3 * (nProt + a) + 2] / ff.nLigAtoms;
      }
    } else lcom = centroidOf(ff.ref, nProt);
    const pocketIdx = pocketFromCom(ff.ref, nProt, lcom, THERMO_POCKET_RCUT);
    // holo frames from recorder; energies from the SELECTED ligand only:
    // "all" (or full-coverage single-ligand) reuses the BindLog channel
    // bit-identically; a proper subset is recomputed per recorded frame with
    // a selected-only FF copying the live binding flags (charges/hbMode).
    const holoFrames = recorder.frames.map((f) => Float32Array.from(f));
    const totalLigAtoms = (state.ligands ?? []).reduce((s, m) => s + (m?.atoms?.length ?? 0), 0);
    const coversAll = selAtomIdx.length === 0
      || selAtomIdx.length === ff.nLigAtoms
      || (totalLigAtoms > 0 && selAtomIdx.length === totalLigAtoms)
      || resolved.mode === "all";
    const holoEnergies = [];
    let holoEnergyNote = "";
    if (coversAll) {
      if (state.bindLog && state.bindLog.nEvents > 0) {
        // collect the last 7-term energy snapshot per captured frame time
        const perTime = new Map();
        for (let i = 0; i < state.bindLog.nEvents; i++) {
          if (state.bindLog.evType[i] !== 0) continue;
          const t = state.bindLog.evTime[i];
          if (!perTime.has(t)) perTime.set(t, new Array(7).fill(0));
          perTime.get(t)[state.bindLog.evA[i] % 7] = state.bindLog.evX[i];
        }
        for (const row of perTime.values()) holoEnergies.push(row);
      }
    } else {
      try {
        const ffSel = new ForceField(state.sel, {
          rc: ff.rc ?? 10, gamma: ff.gamma ?? 2.0,
          binding: {
            on: true, holo: ff.holoOn ?? true,
            charges: ff.chargesOn === true, hbMode: ff.hbMode ?? "off",
          },
        }, resolved.subset);
        ffSel.trackTerms = true;
        for (const frame of holoFrames) {
          const sliced = sliceSelectedPositions(frame, nProt, selAtomIdx);
          ffSel.compute(sliced);
          holoEnergies.push([ffSel.bindLJU, ffSel.bindCoulU, ffSel.bindHBU, ffSel.desolvU, 0, 0, 0]);
        }
        holoEnergyNote = `\n(ligand ${resolved.label} — ΔH recomputed selected-only)`;
      } catch (err) {
        holoEnergyNote = `\n(note: selected-ligand recompute failed (${err?.message ?? err}) — ΔH skipped)`;
      }
    }
    // apo leg: internal relaxation with binding off (same count as before;
    // Stage-2: built with the SELECTED subset — "all" matches the old call).
    const nApo = Math.min(2000, Math.max(600, holoFrames.length * 4));
    let ffApo = null, integApo = null;
    try {
      ffApo = new ForceField(state.sel, { rc: 10, gamma: 2.0, binding: { on: false } }, resolved.subset.length ? resolved.subset : (coversAll ? (state.ligands ?? []) : resolved.subset));
      integApo = new LangevinIntegrator(ffApo.ref, ffApo, 110.0);
      integApo.setTemperature(300); integApo.setFriction(8.0);
    } catch (err) {
      ui.analysisOut.textContent = "⚠ " + (err && err.message ? err.message : String(err));
      return;
    }
    const apoFrames = [];
    setThermoRunning(true); // thermoBtn off, Cancel on for the chunked run
    if (ui.analysisOut) ui.analysisOut.textContent = `Relaxing apo leg… 0/${nApo} steps (UI stays responsive).`;
    setThermoCaption(`relaxing apo… 0/${nApo} steps`);
    let done = 0;
    const stepChunk = () => {
      if (myGen !== _thermoGen) return; // cancelled by rebuild / newer run
      try {
        const n = Math.min(THERMO_CHUNK_STEPS, nApo - done);
        for (let k = 0; k < n; k++) {
          integApo.step();
          if (done % 2 === 0) apoFrames.push(Float32Array.from(integApo.pos));
          done++;
        }
        setThermoCaption(`relaxing apo… ${done}/${nApo} steps`);
        if (done < nApo) { setTimeout(stepChunk, 0); return; }
        // Stage-3: real LCPO solvent burial over the recorded holo + relaxed
        // apo frames (chunked like the apo leg — progress to #thermoCaption,
        // stale generations discarded). Protein blocks are read from
        // both legs; the SELECTED ligand subset (Stage-2 picker) supplies the
        // cross-burial + free-ligand reference (see thermo_sasa.js method).
        // Layout guard: CG frames pack protein + concatenated ligands from
        // nProt; heavy frames offset the ligand block by ff.ligandStart. When
        // the recorded tail block disagrees with the picker concatenation
        // order, SASA degrades honestly to the whole tail block + note.
        const ligStartHolo = ff.ligandStart ?? nProt;
        const tailLen = holoFrames.length ? holoFrames[0].length / 3 - ligStartHolo : 0;
        let sasaSel = selAtomIdx;
        let sasaEls = selectedLigandElements(state.ligands ?? [], resolved.molIdx);
        let sasaNote = "";
        if (!Number.isInteger(tailLen) || tailLen < 0 || tailLen !== totalLigAtoms
          || selAtomIdx.some((a) => a < 0 || a >= tailLen)) {
          sasaSel = null; // whole tail block fallback
          const tailAtoms = Array.isArray(ff.ligandAtoms) ? ff.ligandAtoms : [];
          sasaEls = (tailAtoms.length === tailLen)
            ? tailAtoms.map((a) => String(a?.element ?? "C").toUpperCase())
            : new Array(Math.max(0, tailLen)).fill("C");
          sasaNote = "\n(note: ligand layout differs from picker order — SASA over the full ligand block)";
        }
        const heavySasa = !!state.heavyMode;
        const sasaOpts = {
          nProt, ligStart: ligStartHolo, holoSel: sasaSel, selElements: sasaEls,
          protRadius: heavySasa ? null : cgBeadExtendedRadius(ff),
          protElements: (heavySasa && Array.isArray(ff.atoms))
            ? ff.atoms.slice(0, nProt).map((a) => String(a?.element ?? "C").toUpperCase())
            : null,
          stride: THERMO_SASA_STRIDE,
          chunkFrames: 25,
          onProgress: (d, t) => setThermoCaption(`solvent SASA… ${d}/${t} frames`),
        };
        setThermoRunning(true); // keep Cancel armed across the SASA leg
        setThermoCaption(`solvent SASA… 0 (stride ${THERMO_SASA_STRIDE})`);
        sasaBurialChunked(holoFrames, apoFrames, sasaOpts).then((sasa) => {
          if (myGen !== _thermoGen) return; // rebuilt / newer run
          const res = computeThermodynamics({
            holoFrames, apoFrames, holoEnergies,
            pocketIdx, nProt, mass: 110, T: 300,
            sasa: { dsasa: sasa.dsasa, se: sasa.se, dLig: sasa.dLig, dProt: sasa.dProt, method: sasa.method, stride: sasa.stride, nHoloEval: sasa.nHoloEval, nApoEval: sasa.nApoEval },
          });
          ui.analysisOut.textContent = formatThermoTable(res) +
            `\n(ligand ${resolved.label} [${thermoRotCount} rotatable]; pocket ${pocketIdx.length} residues @8Å of selected COM)` +
            (holoEnergyNote ||
              (holoEnergies.length ? "" : "\n(note: BindLog capture was off — ΔH from recorded-frame recomputation skipped; run with BindLog on for the component split)")) +
            sasaNote;
          setThermoCaption(`ΔH/ΔS done — ${pocketIdx.length} pocket residues, ${holoFrames.length} holo frames (${resolved.label}), ΔSASA ${sasa.dsasa.toFixed(1)} ± ${sasa.se.toFixed(1)} Å².`);
        }).catch((err) => {
          if (myGen !== _thermoGen) return;
          // SASA failure must not strand the thermo report: fall back to the
          // legacy (ΔSASA = 0) path with an honest note.
          try {
            const res = computeThermodynamics({
              holoFrames, apoFrames, holoEnergies,
              pocketIdx, nProt, mass: 110, T: 300,
            });
            ui.analysisOut.textContent = formatThermoTable(res) +
              `\n(ligand ${resolved.label} [${thermoRotCount} rotatable]; pocket ${pocketIdx.length} residues @8Å of selected COM)` +
              `\n(note: real-SASA failed (${err?.message ?? err}) — solvent term is the legacy ΔSASA = 0)`;
            setThermoCaption("thermo done (legacy solvent term — SASA failed).");
          } catch (err2) {
            try { ui.analysisOut.textContent = "⚠ " + (err2 && err2.message ? err2.message : String(err2)); } catch (_) {}
            setThermoCaption("thermo failed — see report above.");
          }
        }).finally(() => {
          if (myGen === _thermoGen) {
            setThermoRunning(false);
          }
        });
        return;
      } catch (err) {
        if (myGen !== _thermoGen) return;
        try { ui.analysisOut.textContent = "⚠ " + (err && err.message ? err.message : String(err)); } catch (_) {}
        setThermoCaption("thermo failed — see report above.");
      } finally {
        if (myGen === _thermoGen) {
          setThermoRunning(false);
        }
      }
    };
    setTimeout(stepChunk, 0);
  } catch (err) {
    try { ui.analysisOut.textContent = "⚠ " + (err && err.message ? err.message : String(err)); } catch (_) {}
    try { if (myGen === _thermoGen) setThermoRunning(false); } catch (_) {}
  }
});
function centroidOf(ref, n) {
  let x = 0, y = 0, z = 0;
  for (let i = 0; i < n; i++) { x += ref[3 * i]; y += ref[3 * i + 1]; z += ref[3 * i + 2]; }
  return [x / n, y / n, z / n];
}

// ---- Alanine scanning ----------------------------------------------
if (ui.alaScanBtn) ui.alaScanBtn.addEventListener("click", () => {
  const sys = currentScanSystem();
  if (!sys) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
  if ((sys.ff.n - (sys.ff.ligandStart ?? sys.ff.nProt)) <= 0) {
    ui.analysisOut.textContent = "⚠ Ala-scan needs a holo system (load/place a ligand, then Build System).";
    return;
  }
  let ids;
  const txt = ((ui.alaRes && ui.alaRes.value) || "").trim();
  try {
    if (txt) {
      ids = txt.split(/[,\s;]+/).filter(Boolean).map((s) => (/^-?\d+$/.test(s) && sys.mode === "cg" ? Number(s) : s));
    } else {
      ids = pocketResidues(sys, { rCut: 6.0, maxN: 6 }).map((p) => p.resId);
    }
  } catch (err) { ui.analysisOut.textContent = "⚠ " + err.message; return; }
  if (!ids.length) { ui.analysisOut.textContent = "⚠ No pocket residues found (ligand has no protein neighbours ≤ 6 Å)."; return; }
  ids = ids.slice(0, 8); // UI responsiveness cap; script API has no cap
  ui.analysisOut.textContent = `Alanine scanning ${ids.length} residue(s)… (short relaxations)`;
  setTimeout(() => {
    try {
      const { rows, wtHolo, wtApo } = scanPocket(sys, ids, { relaxSteps: sys.mode === "heavy" ? 25 : 80 });
      ui.analysisOut.textContent =
        `Ala-scan (${sys.mode}): WT holo ${fmtE(wtHolo)} · WT apo ${fmtE(wtApo)} kcal/mol (ranking only — |ΔΔG| < 0.05 ≈ noise, no CG hotspot claim)\n` +
        formatMutationTable(rows);
    } catch (err) {
      ui.analysisOut.textContent = "⚠ " + err.message;
    }
  }, 20);
});

// ---- DCCM heatmap ----------------------------------------------------
if (ui.dccmBtn) ui.dccmBtn.addEventListener("click", () => {
  const sys = currentScanSystem();
  if (!sys) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
  if (recorder.count < 3) {
    ui.analysisOut.textContent = "⚠ DCCM needs ≥ 3 recorded frames — press ● Rec, run, then Stop.";
    return;
  }
  ui.analysisOut.textContent = "Computing DCCM…";
  setTimeout(() => {
    try {
      const dcc = computeDCCM(recorder.frames, { nProt: sys.ff.nProt, ref: sys.ff.ref });
      const cv = ui.dccmCanvas;
      if (cv) {
        cv.style.display = "block";
        renderDCCMHeatmap(cv, dcc.matrix, { n: dcc.n, size: 220 });
        _dccmHasData = true;
        setDccmCaption(`DCCM live: ${dcc.n}×${dcc.n} from ${dcc.nFrames} frames — click heatmap for pair. Legend −1…+1.`);
        if (!cv._dccmClick) {
          cv._dccmClick = true;
          cv.addEventListener("click", (e) => {
            const p = dccmPick(cv, e.clientX, e.clientY);
            if (p && viewer) {
              highlightCorrelatedPair(viewer, p.i, p.j);
              ui.analysisOut.textContent += `\n◉ pair (${p.i}, ${p.j}) C=${dcc.matrix[p.i * dcc.n + p.j].toFixed(2)} highlighted in viewer.`;
            }
          });
        }
      }
      const top = topCorrelations(dcc, 6);
      let meanAbs = 0;
      for (let i = 0; i < dcc.matrix.length; i++) meanAbs += Math.abs(dcc.matrix[i]);
      meanAbs /= dcc.matrix.length;
      ui.analysisOut.textContent =
        `DCCM: ${dcc.n}×${dcc.n} from ${dcc.nFrames} frames · mean|C| = ${meanAbs.toFixed(3)}\n` +
        `top pairs: ${top.map((t) => `(${t.i},${t.j}) ${t.c >= 0 ? "+" : ""}${t.c.toFixed(2)}`).join(" · ")}\n` +
        `click the heatmap to highlight a pair in the viewer.`;
    } catch (err) {
      ui.analysisOut.textContent = "⚠ " + err.message;
    }
  }, 20);
});

// ---- SMD unbinding ensemble ------------------------------------------
if (ui.smdBtn) ui.smdBtn.addEventListener("click", () => {
  const sys = currentScanSystem();
  if (!sys) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
  if ((sys.ff.n - (sys.ff.ligandStart ?? sys.ff.nProt)) <= 0) {
    ui.analysisOut.textContent = "⚠ SMD needs a holo system (load/place a ligand, then Build System).";
    return;
  }
  ui.analysisOut.textContent = "SMD: 4 pulls…";
  setTimeout(() => {
    try {
      const pos0 = (state.integ && state.integ.pos.length === sys.ff.n * 3)
        ? Float64Array.from(state.integ.pos) : sys.ff.ref;
      const ens = runPullingEnsemble({
        ff: sys.ff, pos0, nPulls: 4, mode: "velocity", seed: 42,
        pull: { k: 5.0, v: 4.0, nSteps: 2000, dt: 0.0015, sRupture: 6.0, T: Number(ui.temp?.value || 300) },
      });
      const je = jarzynskiFreeEnergy(ens.works, { T: Number(ui.temp?.value || 300) });
      const koff = koffSurrogate(ens.pulls, je.dF, { T: Number(ui.temp?.value || 300) });
      ui.analysisOut.textContent =
        `SMD ×4 (constant-velocity, k=2, v=2 Å/ps):\n` +
        `works = [${Array.from(ens.works).map(fmtE).join(", ")}] kcal/mol\n` +
        `Jarzynski ΔF = ${fmtE(je.dF)}${je.se !== null ? ` ± ${je.se.toFixed(2)}` : ""} · ⟨W⟩ = ${fmtE(je.meanWork)} · dissipated = ${fmtE(je.dissipated)}\n` +
        `rupture ⟨F⟩ = ${fmtE(koff.meanRuptureForce)} kcal/mol/Å · koff-score = ${fmtE(koff.meanScore)} kT (${koff.note})`;
    } catch (err) {
      ui.analysisOut.textContent = "⚠ " + err.message;
    }
  }, 20);
});

// ---- Cryptic pockets ---------------------------------------------------
if (ui.crypticBtn) ui.crypticBtn.addEventListener("click", () => {
  const sys = currentScanSystem();
  if (!sys) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
  try {
    // Pocket set: funnel pocket, else protein within 8 Å of the ligand COM.
    let pocket = (state.funnel && state.funnel.active && state.funnel.pocket)
      ? [...state.funnel.pocket] : null;
    if (!pocket) {
      pocket = pocketResidues(sys, { rCut: 8.0, maxN: 40 }).map((p) =>
        (sys.mode === "cg" ? p.resId : sys.sel.atoms.findIndex((a) =>
          a.isProtein && a.atomName === "CA" && `${a.chain}|${a.resSeq}` === p.resId)));
      pocket = pocket.filter((i) => i >= 0);
    }
    if (!pocket.length) { ui.analysisOut.textContent = "⚠ No pocket residues (needs a ligand)."; return; }
    const frames = recorder.count >= 2 ? recorder.frames
      : [(state.integ && state.integ.pos.length === sys.ff.n * 3) ? state.integ.pos : sys.ff.ref];
    const track = trackPocketVolume(frames, { pocketIndices: pocket });
    const det = detectCryptic(track, {
      kSigma: 1.5, frames: frames.length >= 2 ? frames : null,
      pocketIndices: pocket, nProt: sys.ff.nProt,
    });
    const vols = Array.from(track.volumes);
    const topRes = (det.heatmap ? [...det.heatmap.perResidue].sort((a, b) => b.score - a.score).slice(0, 5) : []);
    ui.analysisOut.textContent =
      `Cryptic pockets (${pocket.length} residues, ${track.nFrames} frame(s)):\n` +
      `V range ${fmtE(Math.min(...vols))}–${fmtE(Math.max(...vols))} Å³ · ⟨V⟩ ${fmtE(det.mean)} · open threshold ${fmtE(det.threshold)}\n` +
      `open ${(det.openFrac * 100).toFixed(1)}% · ${det.events.length} event(s)` +
      (det.events.length ? `: ${det.events.slice(0, 4).map((e) => `#${e.start}–${e.end} peak ${fmtE(e.peak)}`).join("; ")}` : "") +
      (topRes.length ? `\nprobe-accessible: ${topRes.map((r) => `#${r.index} ${(r.score * 100).toFixed(0)}%`).join(" · ")}` : "");
  } catch (err) {
    ui.analysisOut.textContent = "⚠ " + err.message;
  }
});
