/**
 * heavy/kernels.js — the harmonic bond and angle energy/force kernels over the
 * flat typed-array lists built by heavy/topology.js.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10).
 *
 *   harmonicFlat / angleFlat              uniform stiffness (legacy path)
 *   harmonicFlatPerK / angleFlatPerK      per-entry stiffness (AMBER ff14SB
 *                                         opt-in path, see useAmber14)
 *
 * These four functions are the most numerically sensitive code in heavy mode:
 * they are called every integration step, they accumulate into a shared force
 * buffer, and their SUMMATION ORDER is part of the result. They were moved
 * verbatim; no expression, no operand order, no accumulator was altered.
 * tests/test_heavy_golden.js pins their output bit-exactly across 6
 * configurations x 12 poses, which is what catches a regrouping like
 * 0.5*k*dr*dr -> 0.5*k*(dr*dr) that a 1e-3 finite-difference tolerance cannot
 * see. Do not "simplify" them.
 *
 * Zero DOM globals. Node-importable.
 */
export function harmonicFlat(pos, f, list, stride, k) {
  let U = 0;
  for (let a = 0; a < list.length; a += stride) {
    const i = 3 * list[a], j = 3 * list[a + 1], r0 = list[a + 2];
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

export function angleFlat(pos, f, list, stride, k) {
  let U = 0;
  for (let a = 0; a < list.length; a += stride) {
    const i = 3 * list[a], j = 3 * list[a + 1], kk = 3 * list[a + 2], th0 = list[a + 3];
    const ax = pos[i] - pos[j], ay = pos[i + 1] - pos[j + 1], az = pos[i + 2] - pos[j + 2];
    const bx = pos[kk] - pos[j], by = pos[kk + 1] - pos[j + 1], bz = pos[kk + 2] - pos[j + 2];
    const la = Math.hypot(ax, ay, az) || 1e-12;
    const lb = Math.hypot(bx, by, bz) || 1e-12;
    let c = (ax * bx + ay * by + az * bz) / (la * lb);
    c = Math.min(1, Math.max(-1, c));
    const th = Math.acos(c);
    const dth = th - th0;
    U += 0.5 * k * dth * dth;
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

/**
 * Per-bond harmonic energy/forces with per-entry k (AMBER14 opt-in path).
 * @param {ArrayLike<number>} pos
 * @param {Float64Array} f
 * @param {ArrayLike<number>} list  flat [i,j,r0...]
 * @param {number} stride
 * @param {ArrayLike<number>} kArr  per-bond stiffness
 * @returns {number}
 */
export function harmonicFlatPerK(pos, f, list, stride, kArr) {
  let U = 0;
  for (let a = 0, b = 0; a < list.length; a += stride, b++) {
    const i = 3 * list[a], j = 3 * list[a + 1], r0 = list[a + 2];
    const k = kArr[b] ?? 200;
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
 * Per-angle bending with per-entry k (AMBER14 opt-in path).
 * @param {ArrayLike<number>} pos
 * @param {Float64Array} f
 * @param {ArrayLike<number>} list  flat [i,j,k,th0...]
 * @param {number} stride
 * @param {ArrayLike<number>} kArr
 * @returns {number}
 */
export function angleFlatPerK(pos, f, list, stride, kArr) {
  let U = 0;
  for (let a = 0, b = 0; a < list.length; a += stride, b++) {
    const k = kArr[b] ?? 40;
    const i = 3 * list[a], j = 3 * list[a + 1], kk = 3 * list[a + 2], th0 = list[a + 3];
    const ax = pos[i] - pos[j], ay = pos[i + 1] - pos[j + 1], az = pos[i + 2] - pos[j + 2];
    const bx = pos[kk] - pos[j], by = pos[kk + 1] - pos[j + 1], bz = pos[kk + 2] - pos[j + 2];
    const la = Math.hypot(ax, ay, az) || 1e-12;
    const lb = Math.hypot(bx, by, bz) || 1e-12;
    let c = (ax * bx + ay * by + az * bz) / (la * lb);
    c = Math.min(1, Math.max(-1, c));
    const th = Math.acos(c);
    const dth = th - th0;
    U += 0.5 * k * dth * dth;
    const sinTh = Math.sqrt(Math.max(1e-12, 1 - c * c));
    const pref = (k * dth) / sinTh;
    const ga = 1 / la, gb = 1 / lb;
    const axy = ax * ga, ayy = ay * ga, azy = az * ga;
    const bxy = bx * gb, byy = by * gb, bzy = bz * gb;
    const fix = pref * (bxy - c * axy) * ga;
    const fiy = pref * (byy - c * ayy) * ga;
    const fiz = pref * (bzy - c * azy) * ga;
    const fkx = pref * (axy - c * bxy) * gb;
    const fky = pref * (ayy - c * byy) * gb;
    const fkz = pref * (azy - c * bzy) * gb;
    f[i] += fix; f[i + 1] += fiy; f[i + 2] += fiz;
    f[kk] += fkx; f[kk + 1] += fky; f[kk + 2] += fkz;
    f[j] -= fix + fkx; f[j + 1] -= fiy + fky; f[j + 2] -= fiz + fkz;
  }
  return U;
}
