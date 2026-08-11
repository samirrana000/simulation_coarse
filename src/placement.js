/**
 * placement.js — rigid-body clash-relaxation placement for library ligands
 * (item 2 of the roadmap).
 *
 * Given a parsed library ligand (mol2.js shape) and a world-space target
 * point (from Viewer.screenToWorld), this module produces a sterically
 * acceptable rigid pose:
 *
 *   1. The molecule is re-centered on its own centroid, given a small
 *      random rotation (seeded → deterministic for tests), and translated
 *      so the centroid lands exactly on the target point.
 *   2. If a protein is present, the pose is relaxed with steepest descent
 *      over the 6 rigid degrees of freedom only (3 translations + 3
 *      rotations about the ligand centroid).
 *
 * RIGIDITY: the pose is stored as (rotation matrix R, centroid c) applied to
 * the centered reference geometry, x_a = R·x0_a + c. Rotations are composed
 * with the exact Rodrigues exponential map, so intramolecular distances are
 * preserved to float64 roundoff no matter how large a step is taken — the
 * small-angle approximation is never used for the pose update itself.
 *
 * Clash objective: for every protein-bead/ligand-atom pair closer than
 * margin·r_e, a quadratic penalty ½k(r − r_e)² with k = 1 (reduced units).
 * The equilibrium pair distance r_e uses the SAME mixing rule as the
 * force field's binding kernel (ff-binding.js): arithmetic-mean σ,
 *   r_e = (σ_protein[i] + σ_ligand[a]) / 2,
 * so a clash-free pose here is also clash-free under the live FF.
 *
 * Gradient (analytic): the force on ligand atom a from a clashing pair is
 * f_a = k(r_e − r)·(x_a − x_i)/r (pushing the atom away from the bead); the
 * net force F = Σ f_a drives translation
 * and the net torque τ = Σ (x_a − c) × f_a about the centroid drives
 * rotation (δU = −F·δx − τ·δω). A simple backtrack line search halves the
 * step whenever the energy increases.
 *
 * Termination: converged as soon as every pair satisfies r ≥ clashMargin·r_e
 * (reported residualClash = max r/r_e over all pairs; ≤ 1 ⇒ clash-free),
 * or maxIters is reached.
 *
 * Cost: O(nProt·m) pair evaluations per iteration with a squared-distance
 * pre-check; scratch buffers are allocated once per call (no per-iteration
 * allocation). Worst case 300 iters × 2000 × 20 pairs ≈ 12M evals — fine.
 *
 * Pure math module: no DOM/browser imports.
 */

import { parseMol2 } from "./mol2.js?v=9";
import { LIG_ELEMENT, LIG_ELEMENT_DEFAULT } from "./ff-params.js?v=9";

// ============================================================================
// Deterministic RNG (mulberry32) — reproducible placement for tests
// ============================================================================

/**
 * mulberry32 PRNG: tiny, fast, seedable; returns uniform doubles in [0, 1).
 * @param {number} seed  32-bit integer seed
 * @returns {() => number}
 */
function mulberry32(seed) {
  let a = (seed >>> 0) || 1;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ============================================================================
// Library parsing (cached per entry id)
// ============================================================================

const _libCache = new Map();   // entry.id -> parsed molecule

/**
 * Parse (and cache) a LIGAND_LIBRARY entry's inline MOL2 text.
 * Returns parseMol2(entry.mol2)[0]; subsequent calls with the same entry id
 * return the SAME object (identity caching — callers must not mutate it).
 * @param {{id: string, mol2: string}} entry  LIGAND_LIBRARY record
 * @returns {{resName: string, chain: string, atoms: Array, bonds: Array}}
 */
export function parseLibraryLigand(entry) {
  let mol = _libCache.get(entry.id);
  if (!mol) {
    mol = parseMol2(entry.mol2)[0];
    _libCache.set(entry.id, mol);
  }
  return mol;
}

// ============================================================================
// Small 3×3 rotation helpers (row-major Float64Array(9))
// ============================================================================

/** R_out = Exp(ω)·R — Rodrigues exponential map, exact for any |ω|. */
function rotExp(ox, oy, oz, R, out) {
  const th = Math.sqrt(ox * ox + oy * oy + oz * oz);
  let kx, ky, kz, A, B;
  if (th < 1e-12) {
    // first-order fallback: Exp(ω) ≈ I + [ω]×
    kx = ox; ky = oy; kz = oz; A = 1; B = 0;
  } else {
    kx = ox / th; ky = oy / th; kz = oz / th;
    A = Math.sin(th); B = 1 - Math.cos(th);
  }
  // E = I·cosθ + sinθ·[k]× + (1−cosθ)·k·kᵀ
  const c = th < 1e-12 ? 1 : Math.cos(th);
  const E0 = c + B * kx * kx, E1 = B * kx * ky - A * kz, E2 = B * kx * kz + A * ky;
  const E3 = B * ky * kx + A * kz, E4 = c + B * ky * ky, E5 = B * ky * kz - A * kx;
  const E6 = B * kz * kx - A * ky, E7 = B * kz * ky + A * kx, E8 = c + B * kz * kz;
  const r0 = R[0], r1 = R[1], r2 = R[2], r3 = R[3], r4 = R[4], r5 = R[5];
  const r6 = R[6], r7 = R[7], r8 = R[8];
  out[0] = E0 * r0 + E1 * r3 + E2 * r6;
  out[1] = E0 * r1 + E1 * r4 + E2 * r7;
  out[2] = E0 * r2 + E1 * r5 + E2 * r8;
  out[3] = E3 * r0 + E4 * r3 + E5 * r6;
  out[4] = E3 * r1 + E4 * r4 + E5 * r7;
  out[5] = E3 * r2 + E4 * r5 + E5 * r8;
  out[6] = E6 * r0 + E7 * r3 + E8 * r6;
  out[7] = E6 * r1 + E7 * r4 + E8 * r7;
  out[8] = E6 * r2 + E7 * r5 + E8 * r8;
}

// ============================================================================
// Clash energy + rigid gradient (translation force / rotation torque)
// ============================================================================

/**
 * Clash energy and rigid-body gradient for a pose given as centroid-relative
 * coords `rel` plus centroid (cx, cy, cz).
 *
 * For each pair with r < margin·r_e: U += 1⁄2k(r − r_e)2 (k = 1); the force on
 * ligand atom a is f_a = k(r − r_e)·(x_a − x_i)/r. Net force F = Σ f_a
 * (translation descent direction) and net torque τ = Σ rel_a × f_a
 * (rotation descent direction) are accumulated. minRatio tracks the min
 * r/r_e over ALL pairs — the worst overlap (nearest pair) — which drives the
 * convergence test and the residual report.
 *
 * @param {Float64Array} rel    ligand coords relative to centroid (3·m)
 * @param {Float64Array} ligSig per-atom σ (m)
 * @param {Float64Array} pp     protein coords (3·nProt)
 * @param {Float64Array} ps     protein per-bead σ (nProt)
 * @param {number} nProt
 * @param {number} cx           centroid x
 * @param {number} cy           centroid y
 * @param {number} cz           centroid z
 * @param {number} margin       count pairs with r < margin·r_e as clashing
 * @param {Float64Array} out    [U, Fx, Fy, Fz, Tx, Ty, Tz, minRatio]
 */
function clashGrad(rel, ligSig, pp, ps, nProt, cx, cy, cz, margin, out) {
  const m = ligSig.length;
  let U = 0, Fx = 0, Fy = 0, Fz = 0, Tx = 0, Ty = 0, Tz = 0, minRatio = Infinity;
  for (let a = 0; a < m; a++) {
    const a3 = 3 * a;
    const rx = rel[a3], ry = rel[a3 + 1], rz = rel[a3 + 2];
    const ax = rx + cx, ay = ry + cy, az = rz + cz;
    const sL = ligSig[a];
    for (let i = 0; i < nProt; i++) {
      const i3 = 3 * i;
      const dx = ax - pp[i3], dy = ay - pp[i3 + 1], dz = az - pp[i3 + 2];
      const r2 = dx * dx + dy * dy + dz * dz;
      const rE = 0.5 * (ps[i] + sL);        // FF mixing rule (ff-binding.js)
      if (r2 < 1e-20) {
        minRatio = 0;                       // coincident pair → worst clash
        continue;
      }
      const ratio = Math.sqrt(r2) / rE;
      if (ratio < minRatio) minRatio = ratio;
      if (r2 >= margin * margin * rE * rE) continue;
      const r = Math.sqrt(r2);
      const dr = rE - r;                    // > 0 when clashing (r < r_e)
      U += 0.5 * dr * dr;                   // 1⁄2k(r − r_e)2, k = 1
      // repulsive force on ligand atom a: f = k(r_e − r)·(x_a − x_i)/r pushes
      // the atom AWAY from the protein bead (down the clash gradient)
      const s = dr / r;
      const fx = s * dx, fy = s * dy, fz = s * dz;
      Fx += fx; Fy += fy; Fz += fz;
      Tx += ry * fz - rz * fy;              // τ += rel_a × f_a
      Ty += rz * fx - rx * fz;
      Tz += rx * fy - ry * fx;
    }
  }
  out[0] = U; out[1] = Fx; out[2] = Fy; out[3] = Fz;
  out[4] = Tx; out[5] = Ty; out[6] = Tz; out[7] = minRatio;
}

/** rel = R·ref for all atoms (pose reconstruction from rigid DOFs). */
function applyRot(R, ref, rel, m) {
  const r0 = R[0], r1 = R[1], r2 = R[2], r3 = R[3], r4 = R[4], r5 = R[5];
  const r6 = R[6], r7 = R[7], r8 = R[8];
  for (let a = 0; a < m; a++) {
    const a3 = 3 * a;
    const x = ref[a3], y = ref[a3 + 1], z = ref[a3 + 2];
    rel[a3]     = r0 * x + r1 * y + r2 * z;
    rel[a3 + 1] = r3 * x + r4 * y + r5 * z;
    rel[a3 + 2] = r6 * x + r7 * y + r8 * z;
  }
}

// ============================================================================
// Rigid clash relaxation (exported core — reusable after protein moves)
// ============================================================================

/**
 * Relax a ligand pose against the protein with steepest descent over the
 * rigid degrees of freedom ONLY (3 translations + 3 exact rotations about
 * the ligand centroid). Intramolecular distances are preserved to float64
 * roundoff by construction.
 *
 * @param {Float64Array} pos  ligand atom coords (3·m), modified IN PLACE
 * @param {{atoms: Array}} mol  parsed molecule (elements → σ via LIG_ELEMENT)
 * @param {{pos: Float64Array, sigma: Float64Array}} protein
 *        protein bead coords (3·nProt) and per-bead σ (from RES_CLASS)
 * @param {object} [opts]
 * @param {number} [opts.maxIters=300]      steepest-descent iteration cap
 * @param {number} [opts.step0=0.05]        initial step (Å per unit force);
 *        halved on energy increase (backtrack line search)
 * @param {number} [opts.clashMargin=0.85]  convergence threshold: stop when
 *        every pair has r ≥ clashMargin·r_e
 * @returns {{pos: Float64Array, iterations: number, residualClash: number,
 *            converged: boolean}}  residualClash = max(r_e/r) over all pairs
 *          (≤ 1 ⇒ clash-free, i.e. every pair at or beyond r_e; the worst
 *          overlap factor); converged ⇒ the clashMargin test passed before
 *          maxIters
 */
export function relaxClash(pos, mol, protein, opts = {}) {
  const maxIters = opts.maxIters ?? 300;
  const step0 = opts.step0 ?? 0.05;
  const clashMargin = opts.clashMargin ?? 0.85;

  const m = mol.atoms.length;
  const pp = protein.pos, ps = protein.sigma;
  const nProt = pp.length / 3;

  // ---- scratch (allocated once per call) -----------------------------------
  const ref = new Float64Array(3 * m);      // centered reference geometry
  const rel = new Float64Array(3 * m);      // R·ref (centroid-relative pose)
  const R = new Float64Array(9);            // cumulative rotation matrix
  const Rn = new Float64Array(9);           // trial rotation
  const grad = new Float64Array(8);         // [U, F(3), T(3), minRatio]

  // per-atom σ from element table (same source as the FF)
  const ligSig = new Float64Array(m);
  for (let a = 0; a < m; a++) {
    ligSig[a] = (LIG_ELEMENT[mol.atoms[a].element] || LIG_ELEMENT_DEFAULT).sigma;
  }

  // Decompose the incoming pose into (ref, R, c): R starts at identity, so
  // the current coordinates ARE the reference, re-centered on the centroid.
  let cx = 0, cy = 0, cz = 0;
  for (let a = 0; a < m; a++) {
    cx += pos[3 * a]; cy += pos[3 * a + 1]; cz += pos[3 * a + 2];
  }
  cx /= m; cy /= m; cz /= m;
  for (let a = 0; a < m; a++) {
    ref[3 * a] = pos[3 * a] - cx;
    ref[3 * a + 1] = pos[3 * a + 1] - cy;
    ref[3 * a + 2] = pos[3 * a + 2] - cz;
  }
  R[0] = R[4] = R[8] = 1;
  rel.set(ref);

  clashGrad(rel, ligSig, pp, ps, nProt, cx, cy, cz, clashMargin, grad);
  let E = grad[0], minRatio = grad[7];
  let iter = 0, converged = minRatio >= clashMargin;

  if (!converged) {
    let step = step0;
    for (iter = 1; iter <= maxIters; iter++) {
      const Fx = grad[1], Fy = grad[2], Fz = grad[3];
      const Tx = grad[4], Ty = grad[5], Tz = grad[6];
      if (Fx === 0 && Fy === 0 && Fz === 0 && Tx === 0 && Ty === 0 && Tz === 0) {
        iter--; break;   // zero gradient, still below margin: stalled
      }

      // trial rigid step: translate by step·F, rotate by ω = step·T with the
      // exact exponential map (pose updated from DOFs → always rigid)
      rotExp(step * Tx, step * Ty, step * Tz, R, Rn);
      applyRot(Rn, ref, rel, m);
      const ncx = cx + step * Fx, ncy = cy + step * Fy, ncz = cz + step * Fz;
      clashGrad(rel, ligSig, pp, ps, nProt, ncx, ncy, ncz, clashMargin, grad);

      if (grad[0] < E) {         // accept trial
        E = grad[0]; minRatio = grad[7];
        R.set(Rn); cx = ncx; cy = ncy; cz = ncz;
        if (minRatio >= clashMargin) { converged = true; break; }
      } else {                   // reject: restore pose, halve the step
        applyRot(R, ref, rel, m);
        clashGrad(rel, ligSig, pp, ps, nProt, cx, cy, cz, clashMargin, grad);
        step *= 0.5;
        if (step < 1e-10) break;   // line search exhausted
      }
    }
  }

  // commit final pose (covers zero-iteration and break paths too)
  for (let a = 0; a < m; a++) {
    pos[3 * a] = rel[3 * a] + cx;
    pos[3 * a + 1] = rel[3 * a + 1] + cy;
    pos[3 * a + 2] = rel[3 * a + 2] + cz;
  }
  // residualClash = max(r_e/r) over all pairs = 1/minRatio — the worst overlap
  // factor: ≤ 1 ⇔ every pair at or beyond r_e (clash-free).
  const residualClash = minRatio > 0 ? 1 / minRatio : Infinity;
  return {
    pos,
    iterations: Math.min(iter, maxIters),
    residualClash,
    converged,
  };
}

// ============================================================================
// Top-level placement
// ============================================================================

/**
 * Place a library ligand at a world-space target point.
 *
 * The molecule is centered on its centroid, given a small random rotation
 * (deterministic under opts.seed), translated so the centroid lands exactly
 * on `target`, and — if opts.protein is given — rigidly relaxed out of
 * clashes via relaxClash().
 *
 * @param {{atoms: Array}} mol  parsed molecule (mol2.js shape)
 * @param {[number, number, number]} target  world point (Å) for the centroid
 * @param {object} [opts]
 * @param {?{pos: Float64Array, sigma: Float64Array}} [opts.protein=null]
 *        protein beads + per-bead σ; null ⇒ no relaxation (no system loaded)
 * @param {number} [opts.seed=1]   PRNG seed for the random rotation
 * @param {number} [opts.maxIters=300]
 * @param {number} [opts.step0=0.05]
 * @param {number} [opts.clashMargin=0.85]
 * @returns {{pos: Float64Array, iterations: number, residualClash: number,
 *            converged: boolean}}
 */
export function placeLigand(mol, target, opts = {}) {
  const m = mol.atoms.length;
  const seed = opts.seed ?? 1;

  // ---- centered reference geometry -----------------------------------------
  let cx = 0, cy = 0, cz = 0;
  for (const at of mol.atoms) { cx += at.x; cy += at.y; cz += at.z; }
  cx /= m; cy /= m; cz /= m;

  // ---- small random rotation (mulberry32, deterministic per seed) ----------
  const rand = mulberry32(seed);
  const ax = (rand() - 0.5) * 1.2;   // ±0.6 rad tumble per axis
  const ay = (rand() - 0.5) * 1.2;
  const az = (rand() - 0.5) * 1.2;
  const cX = Math.cos(ax), sX = Math.sin(ax);
  const cY = Math.cos(ay), sY = Math.sin(ay);
  const cZ = Math.cos(az), sZ = Math.sin(az);

  const pos = new Float64Array(3 * m);
  for (let a = 0; a < m; a++) {
    const x = mol.atoms[a].x - cx, y = mol.atoms[a].y - cy, z = mol.atoms[a].z - cz;
    const y1 = cX * y - sX * z, z1 = sX * y + cX * z;             // Rx
    const x2 = cY * x + sY * z1, z2 = -sY * x + cY * z1;          // Ry
    pos[3 * a]     = cZ * x2 - sZ * y1 + target[0];               // Rz + target
    pos[3 * a + 1] = sZ * x2 + cZ * y1 + target[1];
    pos[3 * a + 2] = z2 + target[2];
  }

  // ---- no protein: return the raw rigid pose --------------------------------
  if (!opts.protein) {
    return { pos, iterations: 0, residualClash: 1.0, converged: true };
  }

  // ---- clash relaxation against the protein ---------------------------------
  return relaxClash(pos, mol, opts.protein, opts);
}
