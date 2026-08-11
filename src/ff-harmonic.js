/**
 * ff-harmonic.js — bonded/two-body harmonic force kernels for the ForceField
 * class in forcefield.js (item 5 modularization).
 *
 * Each function is the verbatim body of the original ForceField method with
 * the instance receiver (`this`) replaced by an explicit `ff` parameter.
 * forcefield.js keeps thin wrapper methods of the same names/signatures, so
 * every existing call site (including the test suite, which reaches into
 * ff._angleForces / ff._ligandBondForces etc.) is untouched.
 */

import { improperAngle } from "./ligand.js?v=8";

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
 * Lets an ML contact prior stiffen/weaken individual residue pairs without
 * rebuilding the whole pair list (e.g. map a model's contact probability
 * onto the ENM stiffness of that pair).
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
 * Shares one kernel between the protein backbone angles and the ligand
 * angle set (the latter is stiffer: k = 40 kcal/mol/rad², see ligand.js).
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
    const pref = (k * dth) / sinTh; // −dU/dθ · dθ/dc → pref = kθ·Δθ / sinθ

    // ∂c/∂r_i etc. (standard 3-body angle gradient; Allen & Tildesley App.)
    const ga = 1 / la, gb = 1 / lb;
    const axy = ax * ga, ayy = ay * ga, azy = az * ga; // unit a
    const bxy = bx * gb, byy = by * gb, bzy = bz * gb; // unit b
    // force on i:  pref * ∂c/∂r_i = pref * (b̂ − c·â)/la
    let fix = pref * (bxy - c * axy) * ga;
    let fiy = pref * (byy - c * ayy) * ga;
    let fiz = pref * (bzy - c * azy) * ga;
    // force on k:  pref * (â − c·b̂)/lb
    let fkx = pref * (axy - c * bxy) * gb;
    let fky = pref * (ayy - c * byy) * gb;
    let fkz = pref * (azy - c * bzy) * gb;

    f[i] += fix; f[i + 1] += fiy; f[i + 2] += fiz;
    f[kk] += fkx; f[kk + 1] += fky; f[kk + 2] += fkz;
    f[j] -= fix + fkx; f[j + 1] -= fiy + fky; f[j + 2] -= fiz + fkz; // Newton's 3rd law
  }
  return U;
}

/**
 * Ligand bonds U = Σ ½ k (r−r0)². The list stores [i,j,r0] only — the force
 * constant is recovered from r0: aromatic ring bonds are clamped into the
 * resonance window r0 ≤ 1.44 Å ⇒ k = 200 kcal/mol/Å², all others (measured
 * r0 ≈ 1.5 Å) ⇒ k = 300 (see constants in ligand.js). Same harmonic kernel
 * as the protein bonds, but with per-bond k.
 */
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
    // Force on j: −∂U/∂r_j = −k·dr·(r̂);  on i the opposite
    const s = (k * dr) / r;
    const fx = s * dx, fy = s * dy, fz = s * dz;
    f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
    f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
  }
  return U;
}

/**
 * Improper (out-of-plane) term U = Σ ½ k (φ−φ0)² with k = 20 kcal/mol/rad²
 * for [i, j, k, l, φ0] (central atom j, planar φ0 = native angle). Forces are
 * taken by central finite differences of φ (step h = 1e-5 Å) because φ is
 * defined through an absolute atan2 in ligand.improperAngle, whose analytic
 * sign is easy to get wrong; FD on the same function is self-consistent by
 * construction (verified against total-U finite differences in the tests).
 */
export function improperForces(pos, f, list) {
  const k = 20;
  const h = 1e-5;
  let U = 0;
  for (let a = 0; a < list.length; a += 5) {
    const i = list[a], j = list[a + 1], kk = list[a + 2], l = list[a + 3], phi0 = list[a + 4];
    const phi = improperAngle(pos, i, j, kk, l);
    U += 0.5 * k * (phi - phi0) * (phi - phi0);
    // −dU/dx_m = −k·(φ−φ0)·dφ/dx_m ; pos is restored after each sweep
    const pref = -k * (phi - phi0);
    for (const m of [i, j, kk, l]) {
      const c = 3 * m;
      for (let ax = 0; ax < 3; ax++) {
        const ci = c + ax, save = pos[ci];
        pos[ci] = save + h;
        const phiP = improperAngle(pos, i, j, kk, l);
        pos[ci] = save - h;
        const phiM = improperAngle(pos, i, j, kk, l);
        pos[ci] = save;
        f[ci] += pref * (phiP - phiM) / (2 * h);
      }
    }
  }
  return U;
}
