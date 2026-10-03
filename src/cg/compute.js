/**
 * cg/compute.js — the CG compute() term assembly. Split out of
 * src/forcefield.js; moved verbatim. THIS IS THE NUMERICALLY LOAD-BEARING
 * FILE OF THE SPLIT.
 *
 * WHAT THIS FILE OWNS
 * -------------------
 * The order in which energy terms are accumulated into a single scalar U and
 * into the flat force buffer. Nothing else decides that order, so nothing
 * else can silently change it.
 *
 * WHY THE SPLIT IS SAFE HERE (and why it is stated, not assumed)
 * ---------------------------------------------------------------
 * U is a float64 sum. Reassociating it changes the last bits, so the ORDER
 * below is a physics contract, not a style choice. The split preserved it
 * literally: the same sequence of `U += ...` statements, the same operand
 * expressions, the same inlined holo-floor loop, in the same positions. No
 * term was hoisted into a local, no call was reordered, no expression was
 * factored. tests/golden/4w52_10steps.json is a bit-exact CG trajectory
 * (posHash/velHash over a seeded 10-step run) and pins the result; it is the
 * proof that this file still produces the same bits.
 *
 * The inlined holo-pinned compression floor is deliberately NOT delegated to
 * a kernel module. It is the only energy term that lives directly in the
 * assembly rather than in ff-harmonic.js, and it accumulates U between the
 * holo-spring call and the angle call — moving it would change the
 * summation order. It stays here, byte-identical.
 *
 * ZERO-ALLOC AUDIT — hot loop (called every integration step):
 *   Expected heap allocs per compute(): ~0 in steady state.
 *   - f.fill(0) reuses preallocated Float64Array(ff.forces) — no alloc.
 *   - harmonicPairs / springForces / angleForces / ligandInternal: pure loops on flat
 *     Float64Array views; only scalar locals (no `new` inside loops).
 *   - _repulsion / _binding: uniform-grid cell list uses preallocated Maps
 *     (`ff._grid`, `ff._gridB`) and reused scratch buffers (`_dens`, `_dBdn`,
 *     `_bpA/J/R`, `_repSigma`); `arr.length=0` clears in place, `arr.push(i)` reuses
 *     existing bucket arrays (capacity grows once, then stable). The `for (const [k,arr] of grid)`
 *     iterator does allocate a lightweight iterator object (~O(cells)), but no per-pair arrays.
 *   - No `new Float64Array` or `new Array` inside the hot path after construction.
 *   Measured via bench/alloc.js (heap delta over 1000 steps). If a future change
 *   adds `new Float64Array` inside compute(), hoist it to a scratch buffer (see
 *   HeavyForceField/SasaModel fix in src/heavy/energy.js).
 */

import { HOLO_FLOOR_RMIN, HOLO_FLOOR_K } from "../ff-params.js";

/**
 * compute(pos) → fills ff.forces and returns total potential energy.
 * @param {object} ff   the ForceField
 * @param {Float64Array} pos  Float64Array(3n); forces are −∇U in kcal/mol/Å
 * @returns {number} total potential energy (kcal/mol), also stored on ff.energy
 */
export function computeForces(ff, pos) {
  const f = ff.forces;
  f.fill(0);
  let U = 0;

  // --- 1/3. Harmonic two-body terms (bonds + ENM springs share the kernel)
  U += ff._harmonicPairs(pos, f, ff.bonds, 3, ff.kBond);
  U += ff._springForces(pos, f);
  if (ff.holoSprings.length) U += ff._harmonicPairs(pos, f, ff.holoSprings, 3, ff.holoGamma);
  if (ff.nativeContacts.length) U += ff._harmonicPairs(pos, f, ff.nativeContacts, 3, 1.0);
  // One-sided compression floor on holo-pinned pairs — keeps the funnel or
  // desolvation from ever squeezing the ligand through a pocket wall
  // (holo pairs are excluded from both grid repulsion and the binding pass,
  // so without this term nothing resists r < r_min).
  if (ff.holoSprings.length) {
    // RMIN/KF now come from ff-params.js (HOLO_FLOOR_RMIN/HOLO_FLOOR_K) —
    // they were duplicated as literals here and in respa.js.
    const HS = ff.holoSprings, RMIN = HOLO_FLOOR_RMIN, KF = HOLO_FLOOR_K;
    for (let a = 0; a < HS.length; a += 3) {
      const i = 3 * HS[a], j = 3 * HS[a + 1];
      const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
      if (r >= RMIN) continue;
      const dr = r - RMIN;                      // < 0
      U += 0.5 * KF * dr * dr;
      const s = (KF * dr) / r;                  // dU/dr ÷ r (dr<0 ⇒ repulsive)
      const fx = s * dx, fy = s * dy, fz = s * dz;
      f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
      f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
    }
  }

  // --- 2. Angle bending ---------------------------------------------------
  U += ff._angleForces(pos, f);

  // --- 3b. Ligand united-atom internal terms (no-op without ligands) ------
  U += ff._ligandInternal(pos, f);

  // --- 4. Excluded volume via grid ----------------------------------------
  U += ff._repulsion(pos, f);

  // --- 5. Protein–ligand binding (cross LJ, electrostatics, H-bonds) -------
  // (physical documentation lives on the _binding wrapper in cg/forcefield.js;
  //  the kernel itself is in ff-binding.js)
  U += ff._binding(pos, f);

  // --- 5c. Optional funnel + well-tempered metadynamics bias (off by default)
  if (ff.funnel && ff.funnelOn) U += ff.funnel.addForces(pos, f);

  // --- 6. NaN/∞ guard ------------------------------------------------------
  // A non-finite coordinate or force kills every downstream pair test and
  // (worse) makes the spatial-hash Map grow unboundedly via garbage cell
  // keys. Detect here (zero cost in the normal case: a single isFinite on
  // the accumulated energy plus a fused isfinite-OR loop on the force
  // buffer), sanitize the force buffer, and leave ff.energy = NaN so the
  // app can auto-pause on the next frame.
  if (!Number.isFinite(U)) {
    let bad = 0;
    for (let i = 0; i < f.length; i++) {
      if (!Number.isFinite(f[i])) { f[i] = 0; bad++; }
    }
    ff.nanStrikes = (ff.nanStrikes ?? 0) + 1;
    if (bad === 0) {
      // energy NaN but forces finite — scan positions for the culprit
      for (let i = 0; i < pos.length && bad < 4; i++) {
        if (!Number.isFinite(pos[i])) bad++;
      }
    }
    U = NaN;
  }
  ff.energy = U;
  return U;
}