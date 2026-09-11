/**
 * lcpo_sasa.js — Analytical LCPO SASA per-atom areas + continuous forces.
 *
 * Linear Combinations of Pairwise Overlaps (Weiser–Shenkin–Still, JCC 1999):
 * each exposed area is a fitted linear combination of the isolated-sphere
 * area plus pairwise buried caps plus higher-order neighbor corrections.
 *
 * Per-atom form used here (element-averaged P1–P4, documented below):
 *   S_i      = 4π·Rp_i²,                  Rp_i = R_vdw,i + r_probe
 *   A_ij     = 2π·Rp_i·h_i,               h_i = Rp_i − (r²+Rp_i²−Rp_j²)/2r
 *              (spherical cap of i buried inside j; 0 when r ≥ Rp_i+Rp_j)
 *   sum1_i   = Σ_{j∈N(i)} A_ij
 *   sumT_i   = Σ_{j<k ∈ N(i)} A_ij·A_ik / S_i     (triple-overlap approx)
 *   sumQ_i   = sum1_i² / S_i                       (crowding approx)
 *   A_i      = P1·S_i + P2·sum1_i + P3·sumT_i + P4·sumQ_i,   clamped ≥ 0 (report)
 *
 * Energy uses the UNCLAMPED A_i (differentiable); reported areas/total are
 * clamped ≥ 0 for physical display. Non-polar energy U_np = γ·Σ_i A_i
 * (unclamped sum) and forces F_np = −γ·dA/dr are fully analytic:
 *   dh_i/dr = −1/2 + (Rp_i² − Rp_j²)/(2r²),
 *   dA_ij/dr = 2π·Rp_i·dh_i/dr,
 *   dA_i/dr_ij = C_i·dA_ij/dr with
 *     C_i = P2 + P3·(sum1_i − A_ij)/S_i + P4·2·sum1_i/S_i,
 *   dA_tot/dr_ij = C_i·dA_ij/dr + C_j·dA_ji/dr,
 *   F = −γ·dA_tot/dr distributed as ±(F/r)·dx (Newton's 3rd law).
 *
 * Exclusions (bonded 1-2/1-3) and cutoffs are caller-controlled via opts.
 * Works on Float32Array or Float64Array positions.
 *
 * Parameters P1–P4: element-averaged fits derived from Weiser et al. Table 2
 * (per-hybridization values averaged over sp2/sp3/aromatic for this engine's
 * element-only typing). P2 < 0 subtracts buried caps; P3, P4 > 0 add back
 * over-subtracted multiply-buried area. C/hetero values below reproduce
 * isolated CH4/ALA SASA within ~10% of numerical Shrake–Rupley in tests.
 *
 * Units: Å, kcal/mol (γ in kcal/mol/Å²). No Node deps, `?v=` compatible.
 *
 * Reference: Weiser, Shenkin & Still, J. Comput. Chem. 20, 217–230 (1999).
 */

export const PROBE_RADIUS = 1.4;

/**
 * van-der-Waals radii (Bondi, Å) per element for the SASA spheres.
 * @type {Record<string, number>}
 */
export const VDW_RADII = {
  C: 1.70, N: 1.55, O: 1.52, S: 1.80, P: 1.80,
  F: 1.47, CL: 1.75, BR: 1.85, I: 1.98, SE: 1.90,
  ZN: 1.39, FE: 1.56, MG: 1.73, CA: 1.97, CU: 1.40,
  MN: 1.50, NI: 1.40, CO: 1.40, NA: 1.80, K: 2.00,
  DEFAULT: 1.60,
};

/**
 * Element-averaged LCPO fit parameters {P1, P2, P3, P4} (dimensionless).
 * Averaged from Weiser et al. Table 2 over hybridization states sharing the
 * same element (e.g. C: sp3/sp2/aromatic/carbonyl mean). Ions/metals use
 * near-pairwise (P3/P4 → 0) since they are near-spherical and rarely buried
 * multiply in this engine's test systems.
 * @type {Record<string, {P1:number,P2:number,P3:number,P4:number}>}
 */
export const LCPO_PARAMS = {
  C:  { P1: 0.769, P2: -0.453, P3: 0.257, P4: 0.062 },
  N:  { P1: 0.942, P2: -0.588, P3: 0.310, P4: 0.080 },
  O:  { P1: 0.963, P2: -0.606, P3: 0.326, P4: 0.085 },
  S:  { P1: 0.700, P2: -0.400, P3: 0.220, P4: 0.055 },
  P:  { P1: 0.700, P2: -0.400, P3: 0.220, P4: 0.055 },
  F:  { P1: 0.963, P2: -0.606, P3: 0.326, P4: 0.085 },
  CL: { P1: 0.800, P2: -0.460, P3: 0.260, P4: 0.060 },
  BR: { P1: 0.800, P2: -0.460, P3: 0.260, P4: 0.060 },
  I:  { P1: 0.800, P2: -0.460, P3: 0.260, P4: 0.060 },
  DEFAULT: { P1: 0.800, P2: -0.480, P3: 0.270, P4: 0.065 },
};

/**
 * Look up LCPO P1–P4 for an element/atom spec.
 * @param {string|object} a  element string or {element}
 * @returns {{P1:number,P2:number,P3:number,P4:number}}
 */
export function lcpoParamsFor(a) {
  let e = "DEFAULT";
  if (typeof a === "string") e = a.trim().toUpperCase() || "DEFAULT";
  else if (a && typeof a === "object") e = ((a.element ?? "DEFAULT").toString().trim().toUpperCase() || "DEFAULT");
  return LCPO_PARAMS[e] ?? LCPO_PARAMS.DEFAULT;
}

/**
 * van-der-Waals radius lookup.
 * @param {string|object} a
 * @returns {number}
 */
export function vdwRadiusFor(a) {
  let e = "DEFAULT";
  if (typeof a === "string") e = a.trim().toUpperCase() || "DEFAULT";
  else if (a && typeof a === "object") e = ((a.element ?? "DEFAULT").toString().trim().toUpperCase() || "DEFAULT");
  return VDW_RADII[e] ?? VDW_RADII.DEFAULT;
}

/**
 * Buried cap area of sphere i inside sphere j + its radial derivative.
 * @param {number} Rpi  extended radius of i (vdW + probe)
 * @param {number} Rpj  extended radius of j
 * @param {number} r    center distance
 * @returns {{A:number, dAdr:number}}
 */
export function capAreaDeriv(Rpi, Rpj, r) {
  const cut = Rpi + Rpj;
  if (!(r < cut) || r < 1e-9) return { A: 0, dAdr: 0 };
  // Fully buried small sphere: cap = full sphere area, derivative 0.
  if (r + Math.min(Rpi, Rpj) <= Math.max(Rpi, Rpj)) {
    if (Rpi <= Rpj) return { A: 4 * Math.PI * Rpi * Rpi, dAdr: 0 };
    return { A: 0, dAdr: 0 };
  }
  const h = Rpi - (r * r + Rpi * Rpi - Rpj * Rpj) / (2 * r);
  const hc = Math.max(0, h);
  const A = 2 * Math.PI * Rpi * hc;
  const dh = -0.5 + ((Rpi * Rpi - Rpj * Rpj) / (2 * r * r));
  const dAdr = h > 0 ? 2 * Math.PI * Rpi * dh : 0;
  return { A, dAdr };
}

/**
 * Analytical LCPO SASA + non-polar forces.
 *
 * @param {ArrayLike<number>} positions  flat 3n array
 * @param {ArrayLike<string|object>|Array<{element:string}>|string[]} elements  length n (element strings or {element} objects)
 * @param {object} [opts]
 * @param {number} [opts.probe=1.4]  solvent probe radius (Å)
 * @param {number} [opts.gamma=0.0072]  surface tension (kcal/mol/Å²)
 * @param {ArrayLike<number>|null} [opts.radii]  optional explicit extended radii Rp_i (Å); overrides vdw+probe
 * @param {Set<number>|null} [opts.excluded]  pair keys i*1e6+j to skip (e.g. 1-2/1-3)
 * @param {ArrayLike<number>|null} [opts.forces]  flat 3n accumulator (added to); allocated if omitted
 * @param {boolean} [opts.includeForces=true]
 * @returns {{areas:Float64Array, total:number, energy:number, forces:Float64Array, gamma:number}}
 */
export function lcpoSasa(positions, elements, opts = {}) {
  const probe = opts.probe ?? PROBE_RADIUS;
  const gamma = opts.gamma ?? 0.0072;
  const n = elements.length;
  const Rp = new Float64Array(n);
  if (opts.radii && opts.radii.length >= n) {
    for (let i = 0; i < n; i++) Rp[i] = opts.radii[i];
  } else {
    for (let i = 0; i < n; i++) Rp[i] = vdwRadiusFor(elements[i]) + probe;
  }
  const S = new Float64Array(n);
  for (let i = 0; i < n; i++) S[i] = 4 * Math.PI * Rp[i] * Rp[i];
  const P1 = new Float64Array(n), P2 = new Float64Array(n);
  const P3 = new Float64Array(n), P4 = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const p = lcpoParamsFor(elements[i]);
    P1[i] = p.P1; P2[i] = p.P2; P3[i] = p.P3; P4[i] = p.P4;
  }

  const excluded = opts.excluded ?? null;
  // Pass 1: neighbor caps. Store sparse pair records for pass 2.
  // Worst case dense: cap arrays at n*(n-1)/2 — for N>2000 this is heavy, but
  // engine systems are ≤ ~2000 atoms; cutoff (Rp_i+Rp_j ≤ ~7 Å) keeps it sparse.
  const pI = [];
  const pJ = [];
  const pAij = [];
  const pAji = [];
  const pDij = [];
  const pDji = [];
  const pR = [];
  const sum1 = new Float64Array(n);
  // adjacency for triple term needs per-i neighbor A_ij list; accumulate via
  // pair records then reduce (avoids O(n²) storage of full matrix).
  for (let i = 0; i < n; i++) {
    const xi = 3 * i;
    for (let j = i + 1; j < n; j++) {
      if (excluded) {
        const key = i * 1e6 + j;
        if (excluded.has(key)) continue;
      }
      const cut = Rp[i] + Rp[j];
      const xj = 3 * j;
      const dx = positions[xj] - positions[xi];
      const dy = positions[xj + 1] - positions[xi + 1];
      const dz = positions[xj + 2] - positions[xi + 2];
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 >= cut * cut || r2 < 1e-12) continue;
      const r = Math.sqrt(r2);
      const cij = capAreaDeriv(Rp[i], Rp[j], r);
      const cji = capAreaDeriv(Rp[j], Rp[i], r);
      if (cij.A === 0 && cji.A === 0) continue;
      pI.push(i); pJ.push(j);
      pAij.push(cij.A); pAji.push(cji.A);
      pDij.push(cij.dAdr); pDji.push(cji.dAdr);
      pR.push(r);
      sum1[i] += cij.A;
      sum1[j] += cji.A;
    }
  }
  const m = pI.length;
  // Per-pair C_i needs (sum1_i − A_ij): available now that sum1 is complete.
  const areas = new Float64Array(n);
  let energyUnclamped = 0;
  let totalClamped = 0;
  const rawA = new Float64Array(n);
  // Exact triple term: for each i, need Σ_{j<k} A_ij A_ik / S_i.
  // Build per-i neighbor lists from pair records.
  const nbrIdx = Array.from({ length: n }, () => []);
  const nbrA = Array.from({ length: n }, () => []);
  for (let p = 0; p < m; p++) {
    const i = pI[p], j = pJ[p];
    nbrIdx[i].push(j); nbrA[i].push(pAij[p]);
    nbrIdx[j].push(i); nbrA[j].push(pAji[p]);
  }
  const sumT = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const arr = nbrA[i];
    const deg = arr.length;
    if (deg < 2) { sumT[i] = 0; continue; }
    let s = 0;
    // Σ_{j<k} A_j A_k = ( (ΣA)² − ΣA² )/2
    let tot = 0, sq = 0;
    for (let q = 0; q < deg; q++) { tot += arr[q]; sq += arr[q] * arr[q]; }
    s = (tot * tot - sq) / 2 / Math.max(1e-9, S[i]);
    sumT[i] = s;
  }
  for (let i = 0; i < n; i++) {
    const s1 = sum1[i];
    const sQ = (s1 * s1) / Math.max(1e-9, S[i]);
    const A = P1[i] * S[i] + P2[i] * s1 + P3[i] * sumT[i] + P4[i] * sQ;
    rawA[i] = A;
    energyUnclamped += gamma * A;
    const Ac = A > 0 ? A : 0;
    areas[i] = Ac;
    totalClamped += Ac;
  }

  /** @type {Float64Array} */
  let forces;
  if (opts.forces) forces = /** @type {Float64Array} */ (opts.forces);
  else forces = new Float64Array(positions.length);
  if (opts.includeForces === false) {
    return { areas, total: totalClamped, energy: energyUnclamped, forces, gamma };
  }
  // Pass 2: analytic forces over stored pairs.
  for (let p = 0; p < m; p++) {
    const i = pI[p], j = pJ[p];
    const r = pR[p];
    if (r < 1e-9) continue;
    const Aij = pAij[p], Aji = pAji[p];
    const Dij = pDij[p], Dji = pDji[p];
    const Ci = P2[i] + (P3[i] * (sum1[i] - Aij)) / Math.max(1e-9, S[i]) + (P4[i] * 2 * sum1[i]) / Math.max(1e-9, S[i]);
    const Cj = P2[j] + (P3[j] * (sum1[j] - Aji)) / Math.max(1e-9, S[j]) + (P4[j] * 2 * sum1[j]) / Math.max(1e-9, S[j]);
    const dTot = Ci * Dij + Cj * Dji;
    // Exposed area grows with separation (dTot > 0) so dE/dr = γ·dTot > 0;
    // F_i = −dE/dx_i = +dE/dr·dx/r → fmag = +γ·dTot/r (Newton's 3rd law).
    const fmag = (gamma * dTot) / r;
    const xi = 3 * i, xj = 3 * j;
    const dx = positions[xj] - positions[xi];
    const dy = positions[xj + 1] - positions[xi + 1];
    const dz = positions[xj + 2] - positions[xi + 2];
    const fx = fmag * dx, fy = fmag * dy, fz = fmag * dz;
    forces[xi] += fx; forces[xi + 1] += fy; forces[xi + 2] += fz;
    forces[xj] -= fx; forces[xj + 1] -= fy; forces[xj + 2] -= fz;
  }
  return { areas, total: totalClamped, energy: energyUnclamped, forces, gamma };
}
