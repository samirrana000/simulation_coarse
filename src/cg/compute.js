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
 * a kernel module and is deliberately NOT a registry term. It is the only
 * energy term that lives directly in the assembly rather than in ff-harmonic.js,
 * and it accumulates U pair-by-pair into the running total that already carries
 * bonds + springs + holo springs + native contacts. A registry term returns ONE
 * float64, so giving the floor its own accumulator would reassociate the sum
 * and move the last bits. It stays here, byte-identical, and the plan walk in
 * computeForces is SPLIT AROUND IT at HOLO_FLOOR_BARRIER_ORDER.
 *
 * THE TERM REGISTRY (2026-10)
 * ----------------------------
 * The nine remaining `U += …` statements are no longer written here. They are
 * declarative descriptors in src/cg/terms.js, validated by
 * src/physics/terms/registry.js, and this file walks the resulting plan:
 *
 *     const plan = cgPlan();                     // memoised, 1 int compare
 *     for (let t = 0; t < plan.barrier; t++) { … }   // orders <  50
 *     …inlined holo floor, byte-identical…
 *     for (let t = plan.barrier; t < n; t++) { … }   // orders >= 50
 *
 * Each iteration is `U += d.energy.call(ff, pos, f)` — the identical call the
 * inlined statement made, returning the identical scalar, at the identical
 * position. A third-party term added to the registry lands in the same walk
 * with no edit to this file; see docs/PHYSICS_TERMS.md. Dispatch costs one
 * ~2 ns `.call` per term against a 0.167 ms compute (0.011 %), i.e. inside the
 * benchmark's own spread.
 *
 * ZERO-ALLOC AUDIT — hot loop (called every integration step):
 *   Expected heap allocs per compute(): ~0 in steady state.
 *   - f.fill(0) reuses preallocated Float64Array(ff.forces) — no alloc.
 *   - cgPlan() returns a memoised frozen object; the walk allocates nothing.
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
import { cgPlan } from "./terms.js";

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

  // --- Registry plan, segment 1: terms with order < HOLO_FLOOR_BARRIER_ORDER.
  // Bonds (cg.bonds, 10), ENM springs (cg.enm, 20), holo springs (cg.holoSprings,
  // 30) and intra-ligand native contacts (cg.nativeContacts, 40) — descriptors
  // in src/cg/terms.js. Each is `U += <the identical scalar>` at the identical
  // position as the statement it replaced; a disabled term is SKIPPED, so U is
  // not touched, which is what the `if (ff.holoSprings.length)` guard did.
  const plan = cgPlan();
  const terms = plan.terms;
  const bar = plan.barrier;
  const nTerms = terms.length;
  for (let t = 0; t < bar; t++) {
    const d = terms[t];
    const en = d.enabled;
    if (en !== null && !en.call(ff)) continue;
    U += d.energy.call(ff, pos, f);
  }

  // --- 1/4. One-sided compression floor on holo-pinned pairs — keeps the
  // funnel or desolvation from ever squeezing the ligand through a pocket wall
  // (holo pairs are excluded from both grid repulsion and the binding pass,
  // so without this term nothing resists r < r_min).
  // NOT A REGISTRY TERM, and byte-identical to the pre-registry statement:
  // it accumulates U PAIR BY PAIR into the running total, so it cannot return
  // one float64 (see the header note and docs/PHYSICS_TERMS.md). This is why
  // the plan walk above is cut in two here.
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

  // --- Registry plan, segment 2: orders >= HOLO_FLOOR_BARRIER_ORDER.
  // Angle bending (60), ligand united-atom internal FF (70), grid excluded
  // volume (80), protein–ligand binding incl. EEF1-lite desolvation (90) and
  // the optional funnel bias (100). The physical documentation for the binding
  // kernel lives on the _binding wrapper in cg/forcefield.js; the kernel itself
  // is in ff-binding.js.
  for (let t = bar; t < nTerms; t++) {
    const d = terms[t];
    const en = d.enabled;
    if (en !== null && !en.call(ff)) continue;
    U += d.energy.call(ff, pos, f);
  }

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