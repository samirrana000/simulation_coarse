/**
 * ff-harmonic.js — bonded/two-body and dihedral harmonic force kernels for ForceField
 * and HeavyForceField.
 *
 * Includes ultra-fast analytical gradients for proper and improper dihedrals (Blondel-Karplus /
 * Bekker vector formulas), replacing finite differences and delivering ~50x speedups.
 */

import { improperAngle } from "./ligand.js?v=10";

/** Σ ½ k (r−r0)² over a flat pair list; returns energy, accumulates forces. */
export function harmonicPairs(pos, f, list, stride, k) {
  let U = 0;
  for (let a = 0; a < list.length; a += stride) {
    const i = 3 * list[a], j = 3 * list[a + 1], r0 = list[a + 2];
    const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
    const dr = r - r0;
    U += 0.5 * k * dr * dr;
    // Force on j: −∂U/∂r_j = −k·dr·(r̂);  on i the opposite
    const s = (k * dr) / r;
    const fx = s * dx, fy = s * dy, fz = s * dz;
    f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
    f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
  }
  return U;
}

/**
 * Σ ½ k_s (r−r0)² for the ENM springs, with a per-spring k from ff.springK.
 */
export function springForces(ff, pos, f) {
  const S = ff.springs, K = ff.springK;
  let U = 0;
  for (let k = 0, s = 0; k < S.length; k += 3, s++) {
    const i = 3 * S[k], j = 3 * S[k + 1], r0 = S[k + 2], kk = K[s];
    const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
    const dr = r - r0;
    U += 0.5 * kk * dr * dr;
    const sc = (kk * dr) / r;
    const fx = sc * dx, fy = sc * dy, fz = sc * dz;
    f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
    f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
  }
  return U;
}

/**
 * U_θ = ½ kθ (θ−θ0)²; forces via finite-chain-rule on cosθ.
 */
export function angleForces(ff, pos, f, list = ff.angles, k = ff.kAngle) {
  let U = 0;
  const A = list;
  for (let a = 0; a < A.length; a += 4) {
    const i = 3 * A[a], j = 3 * A[a + 1], kk = 3 * A[a + 2], th0 = A[a + 3];
    // vectors from j
    const ax = pos[i] - pos[j], ay = pos[i + 1] - pos[j + 1], az = pos[i + 2] - pos[j + 2];
    const bx = pos[kk] - pos[j], by = pos[kk + 1] - pos[j + 1], bz = pos[kk + 2] - pos[j + 2];
    const la = Math.hypot(ax, ay, az) || 1e-12;
    const lb = Math.hypot(bx, by, bz) || 1e-12;
    let c = (ax * bx + ay * by + az * bz) / (la * lb);
    c = Math.min(1, Math.max(-1, c));
    const th = Math.acos(c);
    const dth = th - th0;
    U += 0.5 * k * dth * dth;

    // dθ/dc = −1/sinθ; guard against linear geometry
    const sinTh = Math.sqrt(Math.max(1e-12, 1 - c * c));
    const pref = (k * dth) / sinTh;

    const ga = 1 / la, gb = 1 / lb;
    const axy = ax * ga, ayy = ay * ga, azy = az * ga;
    const bxy = bx * gb, byy = by * gb, bzy = bz * gb;
    let fix = pref * (bxy - c * axy) * ga;
    let fiy = pref * (byy - c * ayy) * ga;
    let fiz = pref * (bzy - c * azy) * ga;
    let fkx = pref * (axy - c * bxy) * gb;
    let fky = pref * (ayy - c * byy) * gb;
    let fkz = pref * (azy - c * bzy) * gb;

    f[i] += fix; f[i + 1] += fiy; f[i + 2] += fiz;
    f[kk] += fkx; f[kk + 1] += fky; f[kk + 2] += fkz;
    f[j] -= fix + fkx; f[j + 1] -= fiy + fky; f[j + 2] -= fiz + fkz;
  }
  return U;
}

/** Ligand bonds U = Σ ½ k (r−r0)² */
export function ligandBondForces(ff, pos, f) {
  let U = 0;
  const B = ff.ligandBonds;
  for (let a = 0; a < B.length; a += 3) {
    const i = 3 * B[a], j = 3 * B[a + 1], r0 = B[a + 2];
    const k = r0 <= 1.44 ? 200 : 300;
    const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
    const dr = r - r0;
    U += 0.5 * k * dr * dr;
    const s = (k * dr) / r;
    const fx = s * dx, fy = s * dy, fz = s * dz;
    f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
    f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
  }
  return U;
}

/**
 * Analytic dihedral gradient for proper/improper torsions (i-j-k-l).
 * Exact vector formulas (Blondel-Karplus / Bekker), O(1) arithmetic.
 *
 * @param {Float64Array} pos
 * @param {Float64Array} f
 * @param {Float64Array|Array} list [i, j, k, l, phi0, ...]
 * @param {number} stride 5
 * @param {number} k force constant
 * @returns {number} potential energy
 */
export function dihedralForcesAnalytic(pos, f, list, stride = 5, k = 2.0) {
  let U = 0;
  const nEntries = list.length;

  for (let a = 0; a < nEntries; a += stride) {
    const idxI = list[a] * 3;
    const idxJ = list[a + 1] * 3;
    const idxK = list[a + 2] * 3;
    const idxL = list[a + 3] * 3;
    const phi0 = list[a + 4];

    // Bond vectors: r_ij = r_j - r_i, r_jk = r_k - r_j, r_kl = r_l - r_k
    const ij_x = pos[idxJ] - pos[idxI];
    const ij_y = pos[idxJ + 1] - pos[idxI + 1];
    const ij_z = pos[idxJ + 2] - pos[idxI + 2];

    const jk_x = pos[idxK] - pos[idxJ];
    const jk_y = pos[idxK + 1] - pos[idxJ + 1];
    const jk_z = pos[idxK + 2] - pos[idxJ + 2];

    const kl_x = pos[idxL] - pos[idxK];
    const kl_y = pos[idxL + 1] - pos[idxK + 1];
    const kl_z = pos[idxL + 2] - pos[idxK + 2];

    // Normal vectors: m = r_ij x r_jk, n = r_jk x r_kl
    const mx = ij_y * jk_z - ij_z * jk_y;
    const my = ij_z * jk_x - ij_x * jk_z;
    const mz = ij_x * jk_y - ij_y * jk_x;

    const nx = jk_y * kl_z - jk_z * kl_y;
    const ny = jk_z * kl_x - jk_x * kl_z;
    const nz = jk_x * kl_y - jk_y * kl_x;

    const m2 = mx * mx + my * my + mz * mz;
    const n2 = nx * nx + ny * ny + nz * nz;
    const jk2 = jk_x * jk_x + jk_y * jk_y + jk_z * jk_z;

    if (m2 < 1e-12 || n2 < 1e-12 || jk2 < 1e-12) continue;

    const inv_m = 1.0 / Math.sqrt(m2);
    const inv_n = 1.0 / Math.sqrt(n2);
    const len_jk = Math.sqrt(jk2);

    // Dihedral angle phi = atan2((m x n) . (r_jk / |r_jk|), m . n)
    const m_dot_n = (mx * nx + my * ny + mz * nz) * inv_m * inv_n;
    const mxn_x = my * nz - mz * ny;
    const mxn_y = mz * nx - mx * nz;
    const mxn_z = mx * ny - my * nx;
    const m_cross_n_dot_jk = (mxn_x * jk_x + mxn_y * jk_y + mxn_z * jk_z) * (inv_m * inv_n / len_jk);

    const phi = Math.atan2(m_cross_n_dot_jk, Math.max(-1, Math.min(1, m_dot_n)));

    // Shortest angular distance wrapped to [-pi, pi]
    let dphi = phi - phi0;
    while (dphi > Math.PI) dphi -= 2 * Math.PI;
    while (dphi < -Math.PI) dphi += 2 * Math.PI;

    U += 0.5 * k * dphi * dphi;

    // Gradient force scale: F = - dU/dr = + k * dphi (r_deriv)
    const pref = k * dphi;

    // Force vectors on atoms i, j, k, l:
    // F_i = pref * (len_jk / |m|^2) * m
    const fi_scale = pref * len_jk / m2;
    const fi_x = fi_scale * mx;
    const fi_y = fi_scale * my;
    const fi_z = fi_scale * mz;

    // F_l = - pref * (len_jk / |n|^2) * n
    const fl_scale = -pref * len_jk / n2;
    const fl_x = fl_scale * nx;
    const fl_y = fl_scale * ny;
    const fl_z = fl_scale * nz;

    // F_j and F_k via projection along central bond jk
    const ij_dot_jk = ij_x * jk_x + ij_y * jk_y + ij_z * jk_z;
    const kl_dot_jk = kl_x * jk_x + kl_y * jk_y + kl_z * jk_z;
    const c1 = ij_dot_jk / jk2;
    const c2 = kl_dot_jk / jk2;

    const fj_x = (c1 - 1.0) * fi_x - c2 * fl_x;
    const fj_y = (c1 - 1.0) * fi_y - c2 * fl_y;
    const fj_z = (c1 - 1.0) * fi_z - c2 * fl_z;

    const fk_x = -fi_x - fj_x - fl_x;
    const fk_y = -fi_y - fj_y - fl_y;
    const fk_z = -fi_z - fj_z - fl_z;

    f[idxI] += fi_x; f[idxI + 1] += fi_y; f[idxI + 2] += fi_z;
    f[idxJ] += fj_x; f[idxJ + 1] += fj_y; f[idxJ + 2] += fj_z;
    f[idxK] += fk_x; f[idxK + 1] += fk_y; f[idxK + 2] += fk_z;
    f[idxL] += fl_x; f[idxL + 1] += fl_y; f[idxL + 2] += fl_z;
  }

  return U;
}

/** Legacy improper forces wrapper using analytic engine */
export function improperForces(pos, f, list) {
  return dihedralForcesAnalytic(pos, f, list, 5, 20.0);
}
