/**
 * metals.js — Transition-metal coordination polyhedra enforcer (Phase 2: Chemistry).
 *
 * Point-charge metal ions (Zn²⁺, Fe²⁺/³⁺, Mg²⁺, Ca²⁺, Mn²⁺, ...) collapse into
 * carboxylates or drift out of their site under pairwise-only potentials.
 * This module adds the missing many-body term: for each metal, N/O/S donors
 * within the coordination radius are classified into a coordination polyhedron
 * by coordination number, and a cross-angle harmonic restraint
 *
 *   U_angle = Σ_pairs ½ k_θ (θ − θ_ideal(geometry))²
 *
 * keeps donor–metal–donor angles at the ideal polyhedron values, while a
 * radial harmonic term U_bond = Σ ½ k_r (r − r0)² pins the shell radius.
 * Together they prevent both collapse (donors piling onto one side / into the
 * ion) and dissociation (donors leaving the shell).
 *
 * Geometry assignment (CN = coordination number):
 *   2 linear (180°) · 3 trigonal-planar (120°) · 4 tetrahedral (109.47°) ·
 *   5 trigonal-bipyramidal (1×180° axial + 3×120° equatorial + 6×90°) ·
 *   6 octahedral (3×180° trans + 12×90° cis) · 7 capped-octahedral ·
 *   8 square-antiprismatic · else irregular (radial only).
 * Trans/axial pairs are chosen by greedy disjoint matching on the largest
 * current angles, so the restraint follows the native arrangement instead of
 * fighting it.
 *
 * Ideal M–donor distances follow CSD surveys (Harding 1999/2004; Barber &
 * Clark metal-site statistics, rounded): Zn–N/O 2.00, Zn–S 2.30, Fe–N/O 2.05,
 * Fe–S 2.30, Mg–O 2.10, Ca–O 2.40, Mn–N/O 2.15 Å.
 *
 * Detection radii default to src/ff-params.js METAL_ELEMENT coordR (imported
 * for consistency with heavy.js buildMetalCoordination); an explicit import is
 * used rather than a copy so the two can never drift apart.
 *
 * Pure ES module, no npm deps. Positions are flat 3n array-likes (engine
 * convention); elements parallel array of uppercase symbols.
 */

import { METAL_ELEMENT, METAL_ELEMENT_DEFAULT } from "../ff-params.js?v=10";

/** Metals handled by the enforcer (formal charges for documentation). */
export const COORDINATION_METALS = Object.freeze({
  ZN: { charge: 2 }, FE: { charge: 2 }, MG: { charge: 2 },
  CA: { charge: 2 }, MN: { charge: 2 }, CU: { charge: 2 },
  NI: { charge: 2 }, CO: { charge: 2 },
});

/** Ideal M–donor distances (Å) per metal → donor element. */
export const IDEAL_METAL_DIST = Object.freeze({
  ZN: { N: 2.00, O: 2.00, S: 2.30 },
  FE: { N: 2.05, O: 2.05, S: 2.30 },
  MG: { N: 2.10, O: 2.10, S: 2.50 },
  CA: { N: 2.50, O: 2.40, S: 2.80 },
  MN: { N: 2.15, O: 2.15, S: 2.45 },
  CU: { N: 2.00, O: 2.00, S: 2.30 },
  NI: { N: 2.05, O: 2.05, S: 2.40 },
  CO: { N: 2.05, O: 2.05, S: 2.35 },
});

const DONOR_ELEMENTS = new Set(["N", "O", "S"]);
const TETRA_ANGLE = Math.acos(-1 / 3); // 109.471°

/**
 * Classify a coordination polyhedron by coordination number.
 * @param {number} cn
 * @returns {string} geometry name
 */
export function classifyGeometry(cn) {
  switch (cn) {
    case 0: return "none";
    case 1: return "monodentate";
    case 2: return "linear";
    case 3: return "trigonal-planar";
    case 4: return "tetrahedral";
    case 5: return "trigonal-bipyramidal";
    case 6: return "octahedral";
    case 7: return "capped-octahedral";
    case 8: return "square-antiprismatic";
    default: return "irregular";
  }
}

/**
 * Detect N/O/S donors around one metal.
 *
 * @param {ArrayLike<number>} pos  flat 3n positions (Å)
 * @param {number} metalIndex
 * @param {Array<string>} elements  parallel element symbols (uppercase)
 * @param {object} [opts]
 * @param {number} [opts.coordR]  detection radius (default: METAL_ELEMENT table)
 * @param {number} [opts.maxDonors]  default: METAL_ELEMENT coordN
 * @returns {{ donorIndices:number[], distances:number[], coordinationNumber:number, geometry:string }}
 */
export function detectCoordination(pos, metalIndex, elements, opts = {}) {
  const el = String(elements[metalIndex] ?? "").toUpperCase();
  const table = METAL_ELEMENT[el] ?? METAL_ELEMENT_DEFAULT;
  const coordR = opts.coordR ?? table.coordR;
  const maxDonors = opts.maxDonors ?? table.coordN;
  const mx = pos[3 * metalIndex], my = pos[3 * metalIndex + 1], mz = pos[3 * metalIndex + 2];
  const found = [];
  for (let j = 0; j < elements.length; j++) {
    if (j === metalIndex) continue;
    if (!DONOR_ELEMENTS.has(String(elements[j] ?? "").toUpperCase())) continue;
    const dx = pos[3 * j] - mx, dy = pos[3 * j + 1] - my, dz = pos[3 * j + 2] - mz;
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
    if (r < coordR) found.push({ j, r });
  }
  found.sort((a, b) => a.r - b.r);
  const take = found.slice(0, maxDonors);
  const donorIndices = take.map((d) => d.j);
  const distances = take.map((d) => d.r);
  return {
    donorIndices, distances,
    coordinationNumber: donorIndices.length,
    geometry: classifyGeometry(donorIndices.length),
  };
}

/**
 * Assign ideal cross-angles for a donor set of a given geometry.
 * Trans/axial slots go to the widest currently-open donor pairs (greedy
 * disjoint matching), so an octahedron keeps its native trans axis.
 *
 * @param {ArrayLike<number>} pos  flat 3n
 * @param {number} metalIndex
 * @param {Array<number>} donors  donor atom indices
 * @param {string} geometry  from classifyGeometry()
 * @returns {Array<{a:number,b:number,ideal:number}>}  donor pairs + ideal angle (rad)
 */
export function assignIdealAngles(pos, metalIndex, donors, geometry) {
  const m = donors.length;
  const pairs = [];
  const angleNow = (a, b) => {
    const ax = pos[3 * a] - pos[3 * metalIndex];
    const ay = pos[3 * a + 1] - pos[3 * metalIndex + 1];
    const az = pos[3 * a + 2] - pos[3 * metalIndex + 2];
    const bx = pos[3 * b] - pos[3 * metalIndex];
    const by = pos[3 * b + 1] - pos[3 * metalIndex + 1];
    const bz = pos[3 * b + 2] - pos[3 * metalIndex + 2];
    const la = Math.hypot(ax, ay, az) || 1e-12;
    const lb = Math.hypot(bx, by, bz) || 1e-12;
    return Math.acos(Math.min(1, Math.max(-1, (ax * bx + ay * by + az * bz) / (la * lb))));
  };
  const allPairs = [];
  for (let x = 0; x < m; x++) {
    for (let y = x + 1; y < m; y++) {
      allPairs.push({ x, y, th: angleNow(donors[x], donors[y]) });
    }
  }
  // Greedy disjoint matching, widest first.
  const used = new Set();
  const pickDisjoint = (count) => {
    const picked = [];
    const sorted = [...allPairs].sort((p, q) => q.th - p.th);
    for (const p of sorted) {
      if (picked.length >= count) break;
      if (used.has(p.x) || used.has(p.y)) continue;
      used.add(p.x); used.add(p.y);
      picked.push(p);
    }
    return picked;
  };
  const D2R = Math.PI / 180;
  const emit = (x, y, deg) => pairs.push({ a: donors[x], b: donors[y], ideal: deg * D2R });

  if (geometry === "linear" && m === 2) {
    emit(0, 1, 180);
  } else if (geometry === "trigonal-planar" && m === 3) {
    emit(0, 1, 120); emit(0, 2, 120); emit(1, 2, 120);
  } else if (geometry === "tetrahedral" && m === 4) {
    for (const p of allPairs) emit(p.x, p.y, 109.471);
  } else if (geometry === "trigonal-bipyramidal" && m === 5) {
    const [axial] = pickDisjoint(1);
    const eqSet = new Set([0, 1, 2, 3, 4]);
    if (axial) {
      eqSet.delete(axial.x); eqSet.delete(axial.y);
      emit(axial.x, axial.y, 180);
    }
    const eq = [...eqSet];
    for (const p of allPairs) {
      const inEqX = eq.includes(p.x), inEqY = eq.includes(p.y);
      if (inEqX && inEqY) emit(p.x, p.y, 120);
      else if (inEqX !== inEqY) emit(p.x, p.y, 90);
    }
  } else if (geometry === "octahedral" && m === 6) {
    const trans = pickDisjoint(3);
    const transKey = new Set(trans.map((p) => `${Math.min(p.x, p.y)}-${Math.max(p.x, p.y)}`));
    for (const p of allPairs) {
      const key = `${Math.min(p.x, p.y)}-${Math.max(p.x, p.y)}`;
      emit(p.x, p.y, transKey.has(key) ? 180 : 90);
    }
  } else {
    // capped / antiprismatic / irregular: distribute evenly on sphere —
    // ideal = mean pairwise angle for the shell (keeps donors spread out).
    let mean = 0;
    for (const p of allPairs) mean += p.th;
    mean = allPairs.length ? mean / allPairs.length : Math.PI / 2;
    for (const p of allPairs) pairs.push({ a: donors[p.x], b: donors[p.y], ideal: mean });
  }
  return pairs;
}

/**
 * Multi-body coordination enforcer: radial shell springs + cross-angle
 * harmonic restraints. Accumulates into opts.forces (or a fresh buffer).
 *
 * @param {ArrayLike<number>} pos  flat 3n positions (Å, read-only)
 * @param {Array<{index:number, element?:string, donors?:number[]}>} metals
 * @param {Array<string>} elements  parallel element symbols
 * @param {object} [opts]
 * @param {Float64Array} [opts.forces]  3n accumulation buffer (allocated if absent)
 * @param {number} [opts.kRadial=40.0]  kcal/mol/Å² (matches heavy.js METAL_K)
 * @param {number} [opts.kAngle=20.0]  kcal/mol/rad²
 * @param {boolean} [opts.ideal=false]  r0 from IDEAL_METAL_DIST table instead of native
 * @param {number} [opts.coordR]  donor detection radius override
 * @returns {{ energy:number, forces:Float64Array, details:Array }}
 *   details: [{ metal, element, donors, distances, r0, coordinationNumber, geometry, nAngles }]
 */
export function enforceCoordination(pos, metals, elements, opts = {}) {
  const n = elements.length;
  const f = opts.forces ?? new Float64Array(n * 3);
  const kR = opts.kRadial ?? 40.0;
  const kA = opts.kAngle ?? 20.0;
  let U = 0;
  const details = [];

  for (const m of metals) {
    const mi = m.index;
    const el = String(m.element ?? elements[mi] ?? "").toUpperCase();
    let donors = m.donors ?? null;
    let distances = null;
    let geometry;
    if (!donors) {
      const det = detectCoordination(pos, mi, elements, opts);
      donors = det.donorIndices;
      distances = det.distances;
      geometry = det.geometry;
    } else {
      distances = donors.map((j) => Math.hypot(
        pos[3 * j] - pos[3 * mi], pos[3 * j + 1] - pos[3 * mi + 1], pos[3 * j + 2] - pos[3 * mi + 2]));
      geometry = classifyGeometry(donors.length);
    }
    const ideal = IDEAL_METAL_DIST[el] ?? IDEAL_METAL_DIST.ZN;
    const r0 = donors.map((j, k) => {
      if (opts.ideal) return ideal[String(elements[j] ?? "").toUpperCase()] ?? 2.1;
      return distances[k]; // pin native shell radius (call at native geometry)
    });

    // --- radial shell springs (prevent dissociation / over-collapse) ---
    for (let k = 0; k < donors.length; k++) {
      const j = donors[k];
      const i3 = 3 * mi, j3 = 3 * j;
      const dx = pos[j3] - pos[i3], dy = pos[j3 + 1] - pos[i3 + 1], dz = pos[j3 + 2] - pos[i3 + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
      const dr = r - r0[k];
      U += 0.5 * kR * dr * dr;
      const s = (kR * dr) / r;
      const fx = s * dx, fy = s * dy, fz = s * dz;
      f[i3] += fx; f[i3 + 1] += fy; f[i3 + 2] += fz;
      f[j3] -= fx; f[j3 + 1] -= fy; f[j3 + 2] -= fz;
    }

    // --- cross-angle polyhedron restraints (prevent one-sided collapse) ---
    const pairs = (geometry === "none" || geometry === "monodentate")
      ? [] : assignIdealAngles(pos, mi, donors, geometry);
    for (const { a, b, ideal: th0 } of pairs) {
      const a3 = 3 * a, m3 = 3 * mi, b3 = 3 * b;
      const ax = pos[a3] - pos[m3], ay = pos[a3 + 1] - pos[m3 + 1], az = pos[a3 + 2] - pos[m3 + 2];
      const bx = pos[b3] - pos[m3], by = pos[b3 + 1] - pos[m3 + 1], bz = pos[b3 + 2] - pos[m3 + 2];
      const la = Math.hypot(ax, ay, az) || 1e-12;
      const lb = Math.hypot(bx, by, bz) || 1e-12;
      let c = (ax * bx + ay * by + az * bz) / (la * lb);
      c = Math.min(1, Math.max(-1, c));
      const th = Math.acos(c);
      const dth = th - th0;
      U += 0.5 * kA * dth * dth;
      const sinTh = Math.sqrt(Math.max(1e-12, 1 - c * c));
      const pref = (kA * dth) / sinTh;
      const ga = 1 / la, gb = 1 / lb;
      const nax = ax * ga, nay = ay * ga, naz = az * ga;
      const nbx = bx * gb, nby = by * gb, nbz = bz * gb;
      const fax = pref * (nbx - c * nax) * ga;
      const fay = pref * (nby - c * nay) * ga;
      const faz = pref * (nbz - c * naz) * ga;
      const fbx = pref * (nax - c * nbx) * gb;
      const fby = pref * (nay - c * nby) * gb;
      const fbz = pref * (naz - c * nbz) * gb;
      f[a3] += fax; f[a3 + 1] += fay; f[a3 + 2] += faz;
      f[b3] += fbx; f[b3 + 1] += fby; f[b3 + 2] += fbz;
      f[m3] -= fax + fbx; f[m3 + 1] -= fay + fby; f[m3 + 2] -= faz + fbz;
    }
    details.push({
      metal: mi, element: el, donors: [...donors], distances: [...distances], r0: [...r0],
      coordinationNumber: donors.length, geometry, nAngles: pairs.length,
    });
  }
  return { energy: U, forces: f, details };
}
