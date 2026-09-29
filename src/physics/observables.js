/**
 * observables.js — the ONE implementation of the two analysis observables both
 * engines need: instantaneous kinetic temperature and RMSD-to-native.
 *
 * WHY: ForceField (Cα CG) and HeavyForceField (all-atom) each carried their own
 * copy of both, and the copies were not the same function — the CG/heavy
 * asymmetry bug class docs/BINDING_PHYSICS_R1..R7 spent seven rounds on. One
 * physical quantity, two silent copies. Now: one pure kernel (no `this`, no
 * class state, no DOM/browser globals — node-importable) plus an explicit
 * argument per call site. Each call site carries a D1/D2 note below.
 *
 * D1 mass layout — CG assumed a flat per-DOF table (3n) unconditionally and
 *    returned NaN on a per-atom table; heavy auto-detected by length. The
 *    layout is now the explicit `mode` argument; default "dof" == old CG.
 * D2 Boltzmann constant — CG read KB_KCAL from units.js (0.001987204), heavy
 *    from ff-params.js (0.0019872041): a 5.03e-8 relative split (~1.5e-5 K at
 *    300 K) in a HUD number that depended only on which engine was selected.
 *    Unified on units.js, the module that declares itself the leaf of the
 *    dependency graph and the single unit contract — and the one every other
 *    k_B consumer already read (integrator, respa, network, funnel, analysis).
 */

import { KB_KCAL, KCONV } from "../units.js?v=10";

/**
 * Instantaneous kinetic temperature from equipartition.
 *   ke_mech = ½ Σ_{i=0..3n-1} m_i·v_i²                     [Da·Å²/ps²]
 *   T       = (ke_mech / KCONV) / (1.5·n·KB_KCAL)         [K]
 * 3n Cartesian DOF for n particles ⇒ 1.5·n·k_B. KCONV (kcal/mol → Da·Å²/ps²)
 * and KB_KCAL (kcal/mol/K) are the src/units.js contract. The 3n sequential
 * adds and the `m·v·v` grouping are verbatim from both former copies, so the
 * result is bit-identical to whichever of them the call site used to run.
 *
 * @param {Float64Array} vel  flat velocities, 3n (x,y,z interleaved)
 * @param {Float64Array|null} mass mass table; layout per `mode` (divergence D1)
 * @param {number} nAtoms    particle count n
 * @param {"dof"|"atom"} [mode="dof"] "dof" ⇒ mass[3a+j], a flat per-DOF table
 *   of length 3n (what LangevinIntegrator builds, integrator.js:78-91 repeats
 *   each particle mass over x,y,z) — the CG contract, hence the default.
 *   "atom" ⇒ mass[a], a per-atom table of length n (what ff.masses and
 *   RESPAStepper expose), repeated over x,y,z. The call site now says which,
 *   instead of the kernel guessing from array length.
 * @returns {number} kelvin; 0 when `mass` is null/undefined (mass-less system).
 */
export function kineticTemp(vel, mass, nAtoms, mode = "dof") {
  if (!mass) return 0;
  const n3 = nAtoms * 3, perAtom = mode === "atom";
  let ke = 0;
  for (let i = 0; i < n3; i++) ke += (perAtom ? mass[(i / 3) | 0] : mass[i]) * vel[i] * vel[i];
  ke *= 0.5 / KCONV;                       // → kcal/mol
  return ke / (1.5 * nAtoms * KB_KCAL);
}

/**
 * RMSD to a reference conformation over a SELECTED set of atoms.
 *   rmsd = sqrt( Σ_{a∈S} Σ_{k∈{x,y,z}} (pos[3a+k] − ref[3a+k])² / |S| )
 *        = sqrt( Σ_coords d² / count )                      [Å]
 * The divisor is the ATOM count, not the coordinate count — i.e. the summed
 * squared deviation over all 3|S| coordinates is divided by |S|, which is the
 * convention both engines have always used (ForceField.rmsd divided by
 * nProt, HeavyForceField.rmsdAll by n). Preserved verbatim: a 1 Å shift in every
 * axis of one atom of an N-atom selection reads sqrt(3/N) Å, not sqrt(3) Å.
 * Changing the divisor would move the HUD number and every golden value.
 * NO implicit alignment (no Kabsch superposition), by design: the ENM keeps the
 * COM and orientation nearly fixed, so aligned and raw metrics coincide and the
 * raw number stays comparable frame-to-frame and across CG/heavy. Do not
 * "improve" this by aligning — it would move the value the HUD reports.
 * Atoms are visited in ascending index order whatever the mask spelling, so the
 * float summation order — hence the last bit — is mask-independent.
 *
 * @param {Float64Array} pos  current flat coordinates, 3n
 * @param {Float64Array} ref  reference flat coordinates, 3n (ff.ref)
 * The scan covers every atom present in `pos`/`ref` — NOT just the first
 * `count` — because a metric's atoms need not be a prefix (the ligand block
 * starts at ligandStart). `count` is purely the divisor.
 *
 * @param {Uint8Array|Int32Array|number[]|null} mask `null` → every atom, the
 *   all-atom metric. Uint8Array → per-atom flags, atom a included iff mask[a]
 *   !== 0; "protein only" and "ligand only" are two cached flag arrays, i.e.
 *   two call sites, not two copies of the loop below. An index list → explicit
 *   ascending atom indices, normalised to flags so there is still exactly one
 *   summation loop.
 * @param {number} count |S|, the divisor: the number of set flags, or n when
 *   `mask` is null. Callers must keep it equal to the flags actually set; the
 *   cached _maskProt/_maskLig arrays satisfy that by construction.
 * @returns {number} Å.
 */
export function rmsdTo(pos, ref, mask, count) {
  const lim = Math.min((pos.length / 3) | 0, (ref.length / 3) | 0);
  if (mask !== null && !(mask instanceof Uint8Array)) {
    const flags = new Uint8Array(lim);
    for (let k = 0; k < Math.min(lim, mask.length); k++) flags[mask[k]] = 1;
    mask = flags;
  }
  let s = 0;
  for (let a = 0; a < lim; a++) {
    if (mask !== null && !mask[a]) continue;
    const i = 3 * a;
    // Three separate adds, NOT (dx²+dy²+dz²): the old loops summed per
    // component and the golden bit is worth preserving.
    const dx = pos[i] - ref[i]; s += dx * dx;
    const dy = pos[i + 1] - ref[i + 1]; s += dy * dy;
    const dz = pos[i + 2] - ref[i + 2]; s += dz * dz;
  }
  return Math.sqrt(s / count);
}
