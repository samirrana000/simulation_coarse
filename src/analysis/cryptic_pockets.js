/**
 * cryptic_pockets.js — transient (cryptic) pocket detection + well-tempered
 * metadynamics along pocket collective variables.
 *
 * Collective variables (Å, Å³; positions in Å):
 *   CV1 — pocket size: radius of gyration Rg of the pocket Cα set (smooth,
 *           biasable) with the convex-hull volume tracked as the observable:
 *             Rg² = (1/m) Σ_k |rk − r̄|²,
 *             V_hull = (1/6) Σ_facets a·(b×c)   (outward-oriented facets).
 *           The hull is exact brute-force for m ≤ 64 pocket atoms (O(m⁴) with
 *           a tiny constant; pocket sets are small) with deterministic stride
 *           subsampling above that, and an ellipsoid fallback V = 4π/3·Rg³
 *           when the set is degenerate (collinear/coplanar).
 *   CV2 — gating distance d = |cᴬ − cᴮ|, the COM separation of two loop
 *           selections flanking the pocket mouth.
 *
 * Well-tempered metadynamics (Barducci–Bussi–Parrinello 2008, cf. funnel.js):
 *   V(s,t) = Σ_k w_k exp(−|s−s_k|²/2σ²),
 *   w_k    = w0 · exp(−V(s_k,t) / (kB·T·(γ−1))),
 *   F(s)   = −(γ/(γ−1)) · V(s)   (Tiwary–Parrinello unbiased estimator form).
 * Bias forces use the analytic gradients of Rg and of the COM distance, so
 * the bias translates pocket/loop atoms without torque artifacts.
 *
 * Detection: a frame is "open" when V_hull > ⟨V⟩ + kσ·σV (default kσ = 1.5);
 * contiguous open frames form events with peak/mean volumes. A
 * probe-accessibility voxel heatmap (probe radius 1.4 Å) quantifies which
 * pocket residues line the transient cavity.
 *
 * DOM-free; runs under plain Node.
 */

import { KB_KCAL } from "../units.js?v=10";

/* ------------------------------------------------------------------ */
/*  Geometry primitives                                                */
/* ------------------------------------------------------------------ */

/**
 * Center of mass of a particle subset.
 * @param {Float64Array} pos  length 3n
 * @param {Array<number>|Int32Array|Uint32Array} indices
 * @param {Float64Array} [masses] per-particle masses (uniform if omitted)
 * @returns {[number, number, number]}
 */
export function subsetCOM(pos, indices, masses = null) {
  let x = 0, y = 0, z = 0, M = 0;
  for (const k of indices) {
    const m = masses ? masses[k] : 1;
    x += m * pos[3 * k]; y += m * pos[3 * k + 1]; z += m * pos[3 * k + 2];
    M += m;
  }
  const inv = 1 / Math.max(1e-30, M);
  return [x * inv, y * inv, z * inv];
}

/**
 * Radius of gyration of a particle subset (Å).
 * @param {Float64Array} pos
 * @param {Array<number>|Int32Array|Uint32Array} indices
 * @returns {number}
 */
export function radiusOfGyration(pos, indices) {
  const [cx, cy, cz] = subsetCOM(pos, indices);
  let s = 0;
  for (const k of indices) {
    const dx = pos[3 * k] - cx, dy = pos[3 * k + 1] - cy, dz = pos[3 * k + 2] - cz;
    s += dx * dx + dy * dy + dz * dz;
  }
  return Math.sqrt(s / Math.max(1, indices.length));
}

/**
 * Analytic gradient of Rg w.r.t. every coordinate (flat 3n vector, zeros
 * outside the subset): ∂Rg/∂rk = (rk − r̄)/(m·Rg).
 * @param {Float64Array} pos
 * @param {Array<number>|Int32Array|Uint32Array} indices
 * @param {number} n  total particle count (gradient length 3n)
 * @returns {{grad: Float64Array, rg: number}}
 */
export function rgGradient(pos, indices, n) {
  const grad = new Float64Array(3 * n);
  const rg = radiusOfGyration(pos, indices);
  if (rg < 1e-9) return { grad, rg };
  const [cx, cy, cz] = subsetCOM(pos, indices);
  const sc = 1 / (indices.length * rg);
  for (const k of indices) {
    grad[3 * k] = (pos[3 * k] - cx) * sc;
    grad[3 * k + 1] = (pos[3 * k + 1] - cy) * sc;
    grad[3 * k + 2] = (pos[3 * k + 2] - cz) * sc;
  }
  return { grad, rg };
}

/**
 * Convex-hull volume of a point set (Å³) by brute-force facet enumeration.
 * Deterministic stride-subsamples to 64 points when larger. Falls back to the
 * ellipsoid estimate V = 4π/3·Rg³ for degenerate sets.
 *
 * @param {Float64Array} pts  length 3m
 * @returns {{volume: number, nFacets: number, approx: boolean}}
 */
export function convexHullVolume(pts) {
  const m = Math.floor(pts.length / 3);
  if (m < 4) {
    const rg = m > 0 ? radiusOfGyration(pts, [...Array(m).keys()]) : 0;
    return { volume: (4 / 3) * Math.PI * rg ** 3, nFacets: 0, approx: true };
  }
  // Deterministic stride subsample for large sets (keeps O(m⁴) bounded).
  let P = pts, mm = m;
  if (m > 64) {
    const stride = m / 64;
    const sub = new Float64Array(64 * 3);
    for (let k = 0; k < 64; k++) {
      const src = Math.min(m - 1, Math.floor(k * stride));
      sub[3 * k] = pts[3 * src]; sub[3 * k + 1] = pts[3 * src + 1]; sub[3 * k + 2] = pts[3 * src + 2];
    }
    P = sub; mm = 64;
  }
  const tol = 1e-9;
  let vol = 0, nFacets = 0;
  for (let i = 0; i < mm; i++) {
    const ax = P[3 * i], ay = P[3 * i + 1], az = P[3 * i + 2];
    for (let j = i + 1; j < mm; j++) {
      const abx = P[3 * j] - ax, aby = P[3 * j + 1] - ay, abz = P[3 * j + 2] - az;
      for (let k = j + 1; k < mm; k++) {
        const acx = P[3 * k] - ax, acy = P[3 * k + 1] - ay, acz = P[3 * k + 2] - az;
        // normal n = ab × ac
        const nx = aby * acz - abz * acy;
        const ny = abz * acx - abx * acz;
        const nz = abx * acy - aby * acx;
        if (nx * nx + ny * ny + nz * nz < 1e-24) continue; // collinear
        const d = nx * ax + ny * ay + nz * az;
        let pos = false, neg = false;
        for (let q = 0; q < mm; q++) {
          const s = P[3 * q] * nx + P[3 * q + 1] * ny + P[3 * q + 2] * nz - d;
          if (s > tol) pos = true;
          else if (s < -tol) neg = true;
          if (pos && neg) break;
        }
        if (pos && neg) continue; // interior triangle, not a hull facet
        // Outward orientation: hull lies on the negative side ⇒ (i,j,k) as
        // ordered is outward; else flip to (i,k,j).
        let b = j, c = k;
        if (pos && !neg) { b = k; c = j; }
        const bx = P[3 * b], by = P[3 * b + 1], bz = P[3 * b + 2];
        const cx = P[3 * c], cy = P[3 * c + 1], cz = P[3 * c + 2];
        vol += (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) / 6;
        nFacets++;
      }
    }
  }
  vol = Math.abs(vol);
  if (nFacets === 0 || vol <= 1e-9) {
    const idx = [...Array(mm).keys()];
    const rg = radiusOfGyration(P, idx);
    return { volume: (4 / 3) * Math.PI * rg ** 3, nFacets, approx: true };
  }
  return { volume: vol, nFacets, approx: false };
}

/**
 * Pocket volume + Rg for one configuration.
 * @param {Float64Array} pos  length 3n
 * @param {Array<number>} pocketIndices
 * @returns {{volume: number, rg: number, approx: boolean, nFacets: number}}
 */
export function pocketVolumeAt(pos, pocketIndices) {
  const m = pocketIndices.length;
  const pts = new Float64Array(m * 3);
  for (let k = 0; k < m; k++) {
    pts[3 * k] = pos[3 * pocketIndices[k]];
    pts[3 * k + 1] = pos[3 * pocketIndices[k] + 1];
    pts[3 * k + 2] = pos[3 * pocketIndices[k] + 2];
  }
  const rg = radiusOfGyration(pts, [...Array(m).keys()]);
  const { volume, nFacets, approx } = convexHullVolume(pts);
  return { volume, rg, approx, nFacets };
}

/* ------------------------------------------------------------------ */
/*  Trajectory tracking                                                */
/* ------------------------------------------------------------------ */

/**
 * Track pocket size observables along frames.
 *
 * @param {Float64Array|Array<Float64Array|Float32Array>} data
 *   single snapshot or array of frames (length 3n each)
 * @param {object} opts
 * @param {Array<number>} opts.pocketIndices  pocket particle indices (required)
 * @param {Array<number>} [opts.loopA]  first gating-loop selection
 * @param {Array<number>} [opts.loopB]  second gating-loop selection
 * @returns {{nFrames: number, volumes: Float64Array, rgs: Float64Array,
 *   loopDists: Float64Array, approxFlags: Uint8Array}}
 */
export function trackPocketVolume(data, opts = {}) {
  const pocketIndices = opts.pocketIndices;
  if (!pocketIndices || pocketIndices.length < 1) {
    throw new Error("trackPocketVolume: opts.pocketIndices (≥ 1 index) is required.");
  }
  const frames = (data instanceof Float64Array || data instanceof Float32Array) ? [data] : data;
  if (!frames || frames.length < 1) throw new Error("trackPocketVolume: no frames given.");
  const nF = frames.length;
  const volumes = new Float64Array(nF);
  const rgs = new Float64Array(nF);
  const loopDists = new Float64Array(nF);
  const approxFlags = new Uint8Array(nF);
  const hasLoops = opts.loopA && opts.loopB && opts.loopA.length > 0 && opts.loopB.length > 0;
  for (let f = 0; f < nF; f++) {
    const pos = frames[f];
    const { volume, rg, approx } = pocketVolumeAt(pos, pocketIndices);
    volumes[f] = volume;
    rgs[f] = rg;
    approxFlags[f] = approx ? 1 : 0;
    if (hasLoops) {
      const ca = subsetCOM(pos, opts.loopA), cb = subsetCOM(pos, opts.loopB);
      loopDists[f] = Math.hypot(ca[0] - cb[0], ca[1] - cb[1], ca[2] - cb[2]);
    } else {
      loopDists[f] = NaN;
    }
  }
  return { nFrames: nF, volumes, rgs, loopDists, approxFlags };
}

/* ------------------------------------------------------------------ */
/*  Probe-accessibility heatmap                                        */
/* ------------------------------------------------------------------ */

/**
 * Voxel probe-accessibility map of the pocket region.
 *
 * A bounding box over all frames' pocket atoms (+pad) is voxelized; a voxel
 * is accessible when no protein atom centre lies within rAtom + rProbe. The
 * accumulated counts (÷ nFrames) form the heatmap; per-residue scores average
 * the accessible fraction of voxels within 4 Å of the residue's mean position.
 *
 * @param {Array<Float64Array|Float32Array>} frames
 * @param {Array<number>} pocketIndices
 * @param {object} [opts]
 * @param {number} [opts.nProt]  protein particle count (clash set 0..nProt−1;
 *   default: pocket bounding atoms only — pass nProt for full occlusion)
 * @param {number} [opts.spacing=1.0]  Å
 * @param {number} [opts.rProbe=1.4]   Å (water-sized probe)
 * @param {number} [opts.rAtom=1.9]    Å (united-atom clash radius)
 * @param {number} [opts.pad=3.0]      Å box padding
 * @returns {{origin: number[], spacing: number, dims: number[],
 *   frac: Float32Array, perResidue: Array<{index:number, score:number}>}}
 */
export function probeAccessibility(frames, pocketIndices, opts = {}) {
  const spacing = opts.spacing ?? 1.0;
  const rProbe = opts.rProbe ?? 1.4;
  const rAtom = opts.rAtom ?? 1.9;
  const pad = opts.pad ?? 3.0;
  const clash2 = (rAtom + rProbe) ** 2;
  const nF = frames.length;
  // Bounding box over pocket atoms across all frames.
  let mnx = Infinity, mny = Infinity, mnz = Infinity, mxx = -Infinity, mxy = -Infinity, mxz = -Infinity;
  for (const pos of frames) {
    for (const k of pocketIndices) {
      const x = pos[3 * k], y = pos[3 * k + 1], z = pos[3 * k + 2];
      if (x < mnx) mnx = x; if (y < mny) mny = y; if (z < mnz) mnz = z;
      if (x > mxx) mxx = x; if (y > mxy) mxy = y; if (z > mxz) mxz = z;
    }
  }
  mnx -= pad; mny -= pad; mnz -= pad; mxx += pad; mxy += pad; mxz += pad;
  let nx = Math.max(2, Math.ceil((mxx - mnx) / spacing) + 1);
  let ny = Math.max(2, Math.ceil((mxy - mny) / spacing) + 1);
  let nz = Math.max(2, Math.ceil((mxz - mnz) / spacing) + 1);
  // Memory guard: coarsen spacing until ≤ 200k voxels.
  let sp = spacing;
  while (nx * ny * nz > 200000) {
    sp *= 1.5;
    nx = Math.max(2, Math.ceil((mxx - mnx) / sp) + 1);
    ny = Math.max(2, Math.ceil((mxy - mny) / sp) + 1);
    nz = Math.max(2, Math.ceil((mxz - mnz) / sp) + 1);
  }
  const counts = new Float64Array(nx * ny * nz);
  const nProt = opts.nProt ?? 0;
  for (const pos of frames) {
    const nAll = Math.floor(pos.length / 3);
    const nClash = nProt > 0 ? Math.min(nProt, nAll) : nAll;
    for (let ix = 0; ix < nx; ix++) {
      const x = mnx + ix * sp;
      for (let iy = 0; iy < ny; iy++) {
        const y = mny + iy * sp;
        for (let iz = 0; iz < nz; iz++) {
          const z = mnz + iz * sp;
          let clash = false;
          for (let k = 0; k < nClash; k++) {
            const dx = pos[3 * k] - x, dy = pos[3 * k + 1] - y, dz = pos[3 * k + 2] - z;
            if (dx * dx + dy * dy + dz * dz < clash2) { clash = true; break; }
          }
          if (!clash) counts[ix * ny * nz + iy * nz + iz]++;
        }
      }
    }
  }
  const frac = new Float64Array(counts.length);
  for (let i = 0; i < counts.length; i++) frac[i] = counts[i] / nF;
  // Per-residue scores: mean accessible fraction within 4 Å of mean position.
  const perResidue = pocketIndices.map((k) => {
    let mx = 0, my = 0, mz = 0;
    for (const pos of frames) { mx += pos[3 * k] / nF; my += pos[3 * k + 1] / nF; mz += pos[3 * k + 2] / nF; }
    let s = 0, c = 0;
    for (let ix = 0; ix < nx; ix++) {
      const x = mnx + ix * sp;
      for (let iy = 0; iy < ny; iy++) {
        const y = mny + iy * sp;
        for (let iz = 0; iz < nz; iz++) {
          const z = mnz + iz * sp;
          const dx = x - mx, dy = y - my, dz = z - mz;
          if (dx * dx + dy * dy + dz * dz <= 16) { s += frac[ix * ny * nz + iy * nz + iz]; c++; }
        }
      }
    }
    return { index: k, score: c ? s / c : 0 };
  });
  return { origin: [mnx, mny, mnz], spacing: sp, dims: [nx, ny, nz], frac, perResidue };
}

/* ------------------------------------------------------------------ */
/*  Cryptic-event detection                                            */
/* ------------------------------------------------------------------ */

/**
 * Detect transiently open (cryptic) pocket states from a volume track.
 *
 * Open ⇔ V > ⟨V⟩ + kSigma·σV. Optionally builds the probe-accessibility
 * heatmap when frames are supplied.
 *
 * @param {{nFrames:number, volumes:Float64Array}} track  from trackPocketVolume
 * @param {object} [opts]
 * @param {number} [opts.kSigma=1.5]
 * @param {Array} [opts.frames]  frames for the heatmap (with pocketIndices)
 * @param {Array<number>} [opts.pocketIndices]
 * @param {number} [opts.nProt]
 * @returns {{mean: number, sd: number, threshold: number, openFrac: number,
 *   nOpen: number, events: Array<{start:number,end:number,peak:number,mean:number}>,
 *   labels: Uint8Array, heatmap: object|null}}
 */
export function detectCryptic(track, opts = {}) {
  const kSigma = opts.kSigma ?? 1.5;
  const vols = Array.from(track.volumes).filter(Number.isFinite);
  if (!vols.length) throw new Error("detectCryptic: track contains no finite volumes.");
  const mean = vols.reduce((a, v) => a + v, 0) / vols.length;
  const sd = Math.sqrt(vols.reduce((a, v) => a + (v - mean) ** 2, 0) / vols.length);
  const threshold = mean + kSigma * sd;
  const labels = new Uint8Array(track.nFrames);
  const events = [];
  let cur = null;
  for (let f = 0; f < track.nFrames; f++) {
    const open = Number.isFinite(track.volumes[f]) && track.volumes[f] > threshold;
    labels[f] = open ? 1 : 0;
    if (open && !cur) cur = { start: f, end: f, peak: track.volumes[f], sum: track.volumes[f], n: 1 };
    else if (open && cur) {
      cur.end = f;
      cur.sum += track.volumes[f]; cur.n++;
      if (track.volumes[f] > cur.peak) cur.peak = track.volumes[f];
    } else if (!open && cur) {
      events.push({ start: cur.start, end: cur.end, peak: cur.peak, mean: cur.sum / cur.n });
      cur = null;
    }
  }
  if (cur) events.push({ start: cur.start, end: cur.end, peak: cur.peak, mean: cur.sum / cur.n });
  const nOpen = labels.reduce((a, v) => a + v, 0);
  let heatmap = null;
  if (opts.frames && opts.pocketIndices) {
    heatmap = probeAccessibility(opts.frames, opts.pocketIndices, { nProt: opts.nProt });
  }
  return {
    mean, sd, threshold,
    openFrac: nOpen / track.nFrames, nOpen, events, labels, heatmap,
  };
}

/* ------------------------------------------------------------------ */
/*  Well-tempered metadynamics on (Rg, loopDist)                       */
/* ------------------------------------------------------------------ */

/**
 * 2-D well-tempered MetaD bias over s = (Rg, loopDist).
 *
 * Hills are stored explicitly (cryptic runs deposit hundreds, not millions)
 * so bias and analytic gradient are exact Gaussian sums:
 *   V(s)   = Σ_k w_k exp(−[(rg−rg_k)²/2σrg² + (d−d_k)²/2σd²]),
 *   ∂V/∂rg = Σ_k −(rg−rg_k)/σrg² · g_k   (analogous for d).
 * Reported free energy uses the Tiwary–Parrinello scaling F = −γ/(γ−1)·V.
 */
export class PocketMetaD {
  /**
   * @param {object} [opts]
   * @param {number} [opts.sigmaRg=0.3]  Å
   * @param {number} [opts.sigmaD=0.5]   Å
   * @param {number} [opts.w0=0.05]      kcal/mol initial hill height
   * @param {number} [opts.gamma=6.0]    bias factor (1 ⇒ plain metadynamics)
   * @param {number} [opts.T=300]        K
   * @param {number} [opts.depositEvery=50]  calls between hills
   */
  constructor(opts = {}) {
    this.sigmaRg = opts.sigmaRg ?? 0.3;
    this.sigmaD = opts.sigmaD ?? 0.5;
    this.w0 = opts.w0 ?? 0.05;
    this.gamma = opts.gamma ?? 6.0;
    this.T = opts.T ?? 300;
    this.depositEvery = opts.depositEvery ?? 50;
    this.hills = []; // {rg, d, w}
    this.calls = 0;
    this._last = { rg: 0, d: 0, V: 0 };
  }

  /** Well-tempered hill height at s given current bias V(s). */
  hillHeight(V) {
    if (this.gamma === 1) return this.w0;
    return this.w0 * Math.exp(-V / (KB_KCAL * this.T * (this.gamma - 1)));
  }

  /**
   * Deposit a hill centred at (rg, d).
   * @param {number} rg
   * @param {number} d
   */
  deposit(rg, d) {
    if (!Number.isFinite(rg) || !Number.isFinite(d)) return;
    const { V } = this.biasAt(rg, d);
    this.hills.push({ rg, d, w: this.hillHeight(V) });
  }

  /**
   * Bias value and analytic gradient at (rg, d).
   * @returns {{V: number, dVdrg: number, dVdd: number}}
   */
  biasAt(rg, d) {
    let V = 0, dVdrg = 0, dVdd = 0;
    const sRg2 = this.sigmaRg * this.sigmaRg, sD2 = this.sigmaD * this.sigmaD;
    for (const h of this.hills) {
      const dr = rg - h.rg, dd = d - h.d;
      const g = h.w * Math.exp(-(dr * dr / (2 * sRg2) + dd * dd / (2 * sD2)));
      V += g;
      dVdrg += (-dr / sRg2) * g;
      dVdd += (-dd / sD2) * g;
    }
    return { V, dVdrg, dVdd };
  }

  /**
   * Apply bias forces to a force accumulator (adds in place) and return the
   * bias energy. Deposits a hill every depositEvery calls.
   *
   * Chain rule: F_k = −(dV/dRg)·∇Rg_k − (dV/dd)·∇d_k with ∇Rg from
   * rgGradient() and ∇d the COM-separation unit vector split over loop atoms.
   *
   * @param {Float64Array} pos  positions, length 3n
   * @param {Float64Array} f    force accumulator, length 3n
   * @param {Array<number>} pocketIndices
   * @param {Array<number>} loopA
   * @param {Array<number>} loopB
   * @returns {{U: number, rg: number, loopDist: number, nHills: number}}
   */
  addForces(pos, f, pocketIndices, loopA, loopB) {
    const n = Math.floor(pos.length / 3);
    const { grad, rg } = rgGradient(pos, pocketIndices, n);
    const ca = subsetCOM(pos, loopA), cb = subsetCOM(pos, loopB);
    const dx = ca[0] - cb[0], dy = ca[1] - cb[1], dz = ca[2] - cb[2];
    const loopDist = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
    const ux = dx / loopDist, uy = dy / loopDist, uz = dz / loopDist;
    const { V, dVdrg, dVdd } = this.biasAt(rg, loopDist);
    for (let i = 0; i < 3 * n; i++) f[i] += -dVdrg * grad[i];
    const fa = -dVdd / loopA.length, fb = dVdd / loopB.length;
    for (const k of loopA) { f[3 * k] += fa * ux; f[3 * k + 1] += fa * uy; f[3 * k + 2] += fa * uz; }
    for (const k of loopB) { f[3 * k] += fb * ux; f[3 * k + 1] += fb * uy; f[3 * k + 2] += fb * uz; }
    this.calls++;
    if (this.calls % this.depositEvery === 0) this.deposit(rg, loopDist);
    this._last = { rg, d: loopDist, V };
    return { U: V, rg, loopDist, nHills: this.hills.length };
  }

  /** Unbiased free-energy estimate F(s) = −γ/(γ−1)·V(s) on a grid. */
  freeEnergyGrid(rgVals, dVals) {
    const scale = this.gamma > 1 ? this.gamma / (this.gamma - 1) : 1;
    return dVals.map((d) => rgVals.map((rg) => -scale * this.biasAt(rg, d).V));
  }

  /** Reset hills and schedule. */
  reset() {
    this.hills.length = 0;
    this.calls = 0;
    this._last = { rg: 0, d: 0, V: 0 };
  }
}
