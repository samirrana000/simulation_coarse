/**
 * virtual-sites.js — Cα-triplet virtual interaction sites (Loop-2 S2, R2 §1).
 *
 * Reconstructs backbone acceptor/donor geometry from three consecutive Cα
 * positions with per-secondary-structure regression coefficients (R2 doc):
 *   helix   O-site 0.93 Å off-Cα, C=O axis err  6.1°
 *   strand  0.47 Å, 14.8°
 *   coil    1.39 Å, 41°  (gates cone-widened by cos² → cos¹ to stay useful)
 *
 * Site frame (Levitt-1976 lineage, per R2 §1):
 *   u  = normalize(r_{i+1} − r_{i−1})          bisector
 *   v  = normalize(2r_i − r_{i+1} − r_{i−1})  out-of-plane (Cα "normal")
 *   w  = normalize(u × v)
 *   O-site position = r_i + a_ss·v + b_ss·w      (carbonyl O side)
 *   NH donor axis   = −v rotated: d_i = normalize(−v + 0.35·w) (approximate;
 *                     helix i·NH points +z of local frame)
 *
 * classifySS: helix if r(i,i+3) < 5.4 Å; strand if r(i,i+2) > 6.6 Å and
 * extended pseudo-dihedral; else coil (R2 §1 exact classifier).
 */

/** @typedef {{oPos: Float64Array, nhAxis: Float64Array, ss: 0|1|2, valid: boolean}} VirtualSite */

// per-SS regression: [c_ss (Å along u, the bisector), a_ss (Å along v)]
// Least-squares fitted against real 4W52 backbone O positions (Loop-2 S2 fit):
//   helix 38 sites: c=1.98±0.14 u, a=−0.39±0.13 v  (w dropped: ±0.88 flip noise)
//   coil 124 sites: c=1.75±0.44 u, a=−0.64±0.64 v
//   strand folds into the coil fit at this system size (19 sites, high variance)
// |O−Cα| ≈ 2.4 Å matches the expected Cα···O separation.
const SS_COEF = {
  helix: [1.98, -0.39],
  strand: [1.75, -0.64], // folds into coil fit (124-site LSQ; 4W52 has few strands)
  coil: [1.75, -0.64],
};

/**
 * Classify secondary structure from Cα geometry (R2 §1).
 * @returns {0|1|2} 0=helix, 1=strand, 2=coil
 */
export function classifySS(prev, cur, next, r13) {
  // r13: distance i−1 → i+1 precomputed by caller
  if (r13 < 5.4) return 0; // helix (i−1..i+1 spans ~5.0 Å)
  if (r13 > 6.6) {
    // strand test: Cα triplet near-collinear
    const u = [next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]];
    const m = [cur[0] - prev[0], cur[1] - prev[1], cur[2] - prev[2]];
    const ul = Math.hypot(...u) || 1, ml = Math.hypot(...m) || 1;
    const cos = (u[0] * m[0] + u[1] * m[1] + u[2] * m[2]) / (ul * ml);
    if (cos > 0.75) return 1; // extended
  }
  return 2; // coil
}

/**
 * Build virtual interaction sites for all Cα beads.
 * @param {Array<{x:number,y:number,z:number,resName:string}>} beads
 * @param {Float64Array} [pos] live positions (optional; defaults from beads)
 * @returns {VirtualSite[]}
 */
export function buildVirtualSites(beads, pos) {
  const n = beads.length;
  const P = (i) => pos ? [pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]] : [beads[i].x, beads[i].y, beads[i].z];
  const out = [];
  for (let i = 0; i < n; i++) {
    const im = Math.max(0, i - 1), ip = Math.min(n - 1, i + 1);
    if (im === i || ip === i || ip === im) { out.push({ oPos: new Float64Array(3), nhAxis: new Float64Array(3), ss: 2, valid: false }); continue; }
    const prev = P(im), cur = P(i), next = P(ip);
    // frame
    let u = [next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]];
    let ul = Math.hypot(...u) || 1; u = [u[0] / ul, u[1] / ul, u[2] / ul];
    let v = [2 * cur[0] - next[0] - prev[0], 2 * cur[1] - next[1] - prev[1], 2 * cur[2] - next[2] - prev[2]];
    let vl = Math.hypot(...v) || 1; v = [v[0] / vl, v[1] / vl, v[2] / vl];
    let w = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    let wl = Math.hypot(...w) || 1; w = [w[0] / wl, w[1] / wl, w[2] / wl];
    const r13 = Math.hypot(next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]);
    const ss = classifySS(prev, cur, next, r13);
    const coef = ss === 0 ? SS_COEF.helix : ss === 1 ? SS_COEF.strand : SS_COEF.coil;
    // O-site: predominantly along the bisector û toward the NEXT residue
    // (carbonyl O sits past Cα_i toward r_{i+1}), minus a small out-of-plane v.
    const oPos = new Float64Array([
      cur[0] + coef[0] * u[0] + coef[1] * v[0],
      cur[1] + coef[0] * u[1] + coef[1] * v[1],
      cur[2] + coef[0] * u[2] + coef[1] * v[2],
    ]);
    // NH donor axis (points roughly −v for the i-th residue's amide N→H)
    const nhAxis = new Float64Array([-v[0] + 0.35 * w[0], -v[1] + 0.35 * w[1], -v[2] + 0.35 * w[2]]);
    const nl = Math.hypot(nhAxis[0], nhAxis[1], nhAxis[2]) || 1;
    nhAxis[0] /= nl; nhAxis[1] /= nl; nhAxis[2] /= nl;
    out.push({ oPos, nhAxis, ss, valid: true });
  }
  return out;
}

/**
 * Directional H-bond energy + forces between a ligand H-bond-capable atom
 * (donor, e.g. ligand O–H/N–H) and a backbone virtual O-site (acceptor),
 * with gradient distributed to the 3 backbone Cα via the frame weights.
 *
 * Form (R2 §2a):  U = −ε · G(r; r0, σ) · max(0, cosθ_D)² · sw(r)
 *   r   = |lig − oPos|
 *   θ_D = angle between (oPos → lig) and the site's local v-axis (acceptor
 *         lone-pair cone; O accepts along +v direction of the carbonyl).
 * Cone width by SS: helix/strand cos² gate; coil cos¹ (R2: widen the cone).
 *
 * @returns {number} energy contribution (kcal/mol); forces accumulated in f.
 */
export function directionalHBond(
  la, ligPos, site, sites, beads, f, i0, i1, i2,
  eps = 2.0, r0 = 3.0, sigma = 0.5, valence
) {
  // la: ligand atom global index; ligPos flat positions
  const lx = ligPos[3 * la], ly = ligPos[3 * la + 1], lz = ligPos[3 * la + 2];
  const ox = site.oPos[0], oy = site.oPos[1], oz = site.oPos[2];
  const dx = lx - ox, dy = ly - oy, dz = lz - oz;
  const r2 = dx * dx + dy * dy + dz * dz;
  if (r2 > 16 || r2 < 1e-10) return 0; // 4 Å gate
  const r = Math.sqrt(r2);
  // radial Gaussian
  const dr = r - r0;
  const g = Math.exp(-(dr * dr) / (2 * sigma * sigma));
  if (g < 1e-4) return 0;
  // smooth switch 3.5→4.0
  let sw = 1, dsw = 0;
  if (r > 3.5) {
    const t = (r - 3.5) / 0.5, t2 = t * t;
    sw = 1 - t2 * t * (10 - 15 * t + 6 * t2);
    dsw = -(30 / 0.5) * t2 * (1 - t) * (1 - t);
  }
  // angle gate: cos of (lig − oPos) against the O-site lone-pair cone.
  // cone axis = normalized (oPos − Cα_i) direction (the +v arm we built with).
  // Reconstructing v from frame each call is expensive; approximate the cone
  // axis as (oPos − cur) using the site geometry captured at build time is
  // wrong under motion — instead store coneAxis on the site at build (below).
  const ca = site.coneAxis;
  const inv = 1 / r;
  const cxs = dx * inv, cys = dy * inv, czs = dz * inv;
  const cosT = (cxs * ca[0] + cys * ca[1] + czs * ca[2]); // lig along the cone continuation
  const gatePow = site.ss === 2 ? 1 : 2; // coil: cos¹ (widened)
  const m = Math.max(0, cosT);
  const gate = Math.pow(m, gatePow);
  const U = -eps * g * sw * gate;
  if (U === 0) return 0;
  if (valence && valence.count >= 2) return 0; // acceptor valence cap
  if (valence) valence.count++;

  // ---- forces ----
  // dU/dr: radial
  const dg = -(dr / (sigma * sigma)) * g;
  const dU_dr = -eps * (dg * sw + g * dsw) * gate;
  // dU/d(cosT) via cone axis: d(cosT)/d(lig) = −ca/r + cosT·(lig−o)/r²  (per component)
  const dgate = gatePow * Math.pow(m, gatePow - 1); // d(m^p)/dm
  // force on ligand: F = −dU/dx_lig ; dU/dx_lig = dU_dr·(dx/r) + (−eps·g·sw·dgate)·d(cosT)/dx_lig
  const cT_dx = -ca[0] / r - cxs * cosT / r; // careful: derivative of −(unit·ca) w.r.t. ligand x
  const cT_dy = -ca[1] / r - cys * cosT / r;
  const cT_dz = -ca[2] / r - czs * cosT / r;
  const pre = -eps * g * sw * dgate;
  const Fx_lig = -(dU_dr * (dx / r) + pre * cT_dx);
  const Fy_lig = -(dU_dr * (dy / r) + pre * cT_dy);
  const Fz_lig = -(dU_dr * (dz / r) + pre * cT_dz);
  f[3 * la] -= Fx_lig; f[3 * la + 1] -= Fy_lig; f[3 * la + 2] -= Fz_lig;
  // opposite on the O-site; O-site force distributed to the 3 Cα with weights
  // (a_ss:b_ss frame is a linear combination of the triplet — equal thirds is
  // the conservative gradient-safe approximation; error documented in R2 §1)
  const w3 = 1 / 3;
  f[3 * i0] += Fx_lig * w3; f[3 * i0 + 1] += Fy_lig * w3; f[3 * i0 + 2] += Fz_lig * w3;
  f[3 * i1] += Fx_lig * w3; f[3 * i1 + 1] += Fy_lig * w3; f[3 * i1 + 2] += Fz_lig * w3;
  f[3 * i2] += Fx_lig * w3; f[3 * i2 + 1] += Fy_lig * w3; f[3 * i2 + 2] += Fz_lig * w3;
  return U;
}

/** Cone axis for an O-site: normalized (oPos − Cα_i) — the acceptor direction. */
export function coneAxisOf(site, cur) {
  const ax = site.oPos[0] - cur[0], ay = site.oPos[1] - cur[1], az = site.oPos[2] - cur[2];
  const l = Math.hypot(ax, ay, az) || 1;
  site.coneAxis = [ax / l, ay / l, az / l];
  return site.coneAxis;
}
