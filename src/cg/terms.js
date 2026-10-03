/**
 * cg/terms.js — the CG (Cα ENM) term descriptors + the CG plan cache.
 *
 * WHAT MOVED HERE
 * ---------------
 * Nine of the ten `U += …` statements that used to be written out inside
 * src/cg/compute.js are now descriptors here. What is left in compute.js is the
 * plan walk, the ONE term that cannot move (the inlined holo compression
 * floor — see HOLO_FLOOR_BARRIER_ORDER below), and the NaN/∞ guard, which is
 * not a term.
 *
 * SUMMATION-ORDER SAFETY OF THE MIGRATION (why the goldens did not move)
 * ---------------------------------------------------------------------
 * The original assembly was already `U += <one call returning one float64>`
 * nine times out of ten, so the descriptor for each of those terms is a
 * one-line arrow that performs the IDENTICAL call and returns the IDENTICAL
 * scalar. Concretely, for each term:
 *   • the call expression is byte-identical (same receiver, same arguments);
 *   • the value is `U += <that scalar>` at the same position in the sequence;
 *   • no two terms are merged into one `+=` and no one term is split;
 *   • a conditionally-run term keeps its `if (...)` semantics through the
 *     `enabled` predicate, and the dispatcher SKIPS rather than adding a zero,
 *     so `U` is not touched at all when the term is off — which is what the
 *     `if (ff.holoSprings.length) U += …` did.
 * The only structural difference is dispatch (one `.call` per term instead of
 * an inlined expression), which is arithmetic-free.
 *
 * HOLO_FLOOR_BARRIER_ORDER — why the plan is walked in TWO segments
 * ---------------------------------------------------------------
 * The one-sided holo compression floor is NOT a descriptor and cannot be one.
 * It accumulates `U += 0.5*KF*dr*dr` PER PAIR directly into the running total
 * that already carries bonds + springs + holo springs + native contacts. A
 * descriptor must return ONE float64, so it would have to accumulate into its
 * own local and then be added once — and that reassociation changes the last
 * bits of the CG ligand golden. So the floor stays inlined in compute.js, and
 * the plan is split around it: orders < 50 run before, orders >= 50 run after.
 * A third-party term picks a side by picking an integer order. There is no
 * order 50 slot: 50 is the barrier's own.
 *
 * PER-TERM PROVENANCE (what each descriptor replaced, and why it was safe)
 * ------------------------------------------------------------------------
 *   order  cg.bonds           U += ff._harmonicPairs(pos, f, ff.bonds, 3, ff.kBond)
 *   order  cg.enm             U += ff._springForces(pos, f)
 *   order  cg.holoSprings     if (ff.holoSprings.length) U += ff._harmonicPairs(…, 3, ff.holoGamma)
 *   order  cg.nativeContacts  if (ff.nativeContacts.length) U += ff._harmonicPairs(…, 3, 1.0)
 *   ----    (inlined holo compression floor — NOT a term, see above) ----
 *   order  cg.angles          U += ff._angleForces(pos, f)
 *   order  cg.ligandInternal  U += ff._ligandInternal(pos, f)
 *   order  cg.repulsion       U += ff._repulsion(pos, f)
 *   order  cg.binding         U += ff._binding(pos, f)
 *   order  cg.funnel          if (ff.funnel && ff.funnelOn) U += ff.funnel.addForces(pos, f)
 *
 * Every one of those was already `U += <one call returning one float64>`, so
 * each descriptor is the same call with the same receiver and the same
 * arguments, returning the same scalar, added at the same position. The three
 * conditional ones become an `enabled` predicate, and the dispatcher SKIPS a
 * disabled term rather than adding a zero, which is exactly what the `if`
 * guard did (the running U is never −0.0, but skip-and-do-not-touch is the
 * stronger statement and costs nothing). Nothing was merged, nothing was
 * split, nothing was hoisted into a local.
 *
 * NO PER-TERM ACCUMULATORS ON THE CG FIELD
 * ----------------------------------------
 * Unlike the heavy engine, the CG assembly never published per-term energies
 * on the force field, and this migration does not start: doing so would add
 * nine own instance fields to ForceField, which tests/test_cg_class_shape.js
 * pins to exactly 70 in a fixed order. So a CG descriptor leaves `reports` null
 * and `describeTerms("cg", ff)` reports `energy: null` — a real, documented
 * limitation of the introspection surface, not a silent gap.
 *
 * Zero DOM globals. Node-importable.
 */
import { installBuiltins, listTerms, describeTerms, TERM_GENERATION } from "../physics/terms/registry.js";

/**
 * The order at which the CG plan is cut in two by the inlined holo
 * compression floor. Orders strictly below this run BEFORE the floor; orders
 * at or above it run AFTER. 50 is unclaimed by any built-in on purpose — it is
 * the barrier's own slot.
 */
export const HOLO_FLOOR_BARRIER_ORDER = 50;

/**
 * The CG terms, in the order they contribute to U.
 *
 * Every `energy` is `function (pos, f)` with `this` = the ForceField: the
 * registry's `this`-binding contract (see physics/terms/registry.js), which is
 * also the convention every kernel in this repo already uses. None of them
 * allocates; each returns the single float64 the dispatcher adds.
 */
export const CG_TERMS = [
  {
    // U_bond = Σ ½ k_b (r − r0)² over ff.bonds, r0 = the input Cα–Cα distance
    // (≈ 3.81 Å), k_b = ff.kBond ≈ 100 kcal/mol/Å². Kernel: harmonicPairs in
    // src/ff-harmonic.js via the _harmonicPairs wrapper on ForceField.
    id: "cg.bonds", engine: "cg", order: 10,
    label: "Cα peptide bonds — harmonic, k = ff.kBond",
    energy: function (pos, f) {
      return this._harmonicPairs(pos, f, this.bonds, 3, this.kBond);
    },
  },
  {
    // U_ENM = Σ ½ γ_ij (r − r0_ij)² over the H(Rc − r0) native-contact set, with
    // a PER-SPRING k from ff.springK (setSpringScale / applySeqWeights /
    // useTirionNetwork all rewrite that table). Kernel: springForces in
    // src/ff-harmonic.js. This is the term that holds the native fold.
    id: "cg.enm", engine: "cg", order: 20,
    label: "ENM native-contact springs — harmonic, per-spring k from ff.springK",
    energy: function (pos, f) {
      return this._springForces(pos, f);
    },
  },
  {
    // U = Σ ½ γ_lig (r − r0)² over the crystallographic native pose, γ_lig =
    // ff.holoGamma = 0.5. Excluded from the grid repulsion and from the binding
    // pair pass, so this spring alone governs the native-contact distance.
    id: "cg.holoSprings", engine: "cg", order: 30,
    label: "Holo (native-pose) protein–ligand contact springs, k = ff.holoGamma",
    // was: `if (ff.holoSprings.length) U += ff._harmonicPairs(…)`. The
    // dispatcher skips the term entirely when this is false, so U is untouched.
    enabled: function () { return this.holoSprings.length > 0; },
    energy: function (pos, f) {
      return this._harmonicPairs(pos, f, this.holoSprings, 3, this.holoGamma);
    },
  },
  {
    // Intra-ligand native contacts, k = 1.0. Empty for a protein-only system,
    // which is why this is a predicate rather than an unconditional call.
    id: "cg.nativeContacts", engine: "cg", order: 40,
    label: "Intra-ligand native contacts — harmonic, k = 1.0",
    enabled: function () { return this.nativeContacts.length > 0; },
    energy: function (pos, f) {
      return this._harmonicPairs(pos, f, this.nativeContacts, 3, 1.0);
    },
  },
  {
    // U_θ = Σ ½ k_θ (θ − θ0)² over the backbone pseudo-angle i−1, i, i+1 within
    // each contiguous segment, θ0 from the input coordinates. Kernel:
    // angleForces in src/ff-harmonic.js. Runs AFTER the inlined holo floor.
    id: "cg.angles", engine: "cg", order: 60,
    label: "Backbone pseudo-angle bending, k = ff.kAngle",
    energy: function (pos, f) {
      return this._angleForces(pos, f);
    },
  },
  {
    // Ligand bonds + angles + out-of-plane impropers from src/ligand.js,
    // assembled by _ligandInternal on the force field. Returns 0 for a
    // protein-only system — which is why this is unconditional: the original
    // `U += ff._ligandInternal(…)` was too, and adding a zero is not the same
    // statement as skipping, so it stays.
    id: "cg.ligandInternal", engine: "cg", order: 70,
    label: "Ligand united-atom internal terms (bonds, angles, impropers)",
    energy: function (pos, f) {
      return this._ligandInternal(pos, f);
    },
  },
  {
    // U_rep = ε[(r_e/r)¹² − 2(r_e/r)⁶ + 1] for r < r_e = 2^{1/6}σ, per-pair
    // σ = arithmetic mean of the two bead sizes. Repulsive tail only: implicit
    // solvent, beads may not collapse. Spatial-hash grid kernel in
    // src/ff-repulsion.js; the grid maps and scratch buffers are preallocated
    // (src/cg/grid.js), so this term is zero-alloc in steady state.
    id: "cg.repulsion", engine: "cg", order: 80,
    label: "Excluded volume — repulsive-only 12-6 LJ over the spatial hash",
    energy: function (pos, f) {
      return this._repulsion(pos, f);
    },
  },
  {
    // Cross 12-6 LJ (Lorentz–Berthelot) + screened electrostatics
    // (εr = 4 + 76·tanh(r/8)) + a directional H-bond well, each smoothstep-
    // switched off at bindRcut ≈ 9 Å, PLUS the two-pass EEF1-lite burial /
    // desolvation pass. Full physics documentation is on the _binding wrapper in
    // src/cg/forcefield.js; the kernel is src/ff-binding.js. This is the single
    // most expensive CG term and the one test_cg_ligand_golden.js exists for.
    id: "cg.binding", engine: "cg", order: 90,
    label: "Protein–ligand binding: cross LJ + screened Coulomb + H-bond + EEF1-lite desolvation",
    energy: function (pos, f) {
      return this._binding(pos, f);
    },
  },
  {
    // External collective-variable bias: flat inside the bound state plus
    // well-tempered metadynamics hills for the PMF. Attached by the CALLER via
    // ForceField.setFunnel() — this module never imports src/funnel.js, which is
    // why the predicate has to test both `funnel` and `funnelOn`.
    id: "cg.funnel", engine: "cg", order: 100,
    label: "Funnel + well-tempered metadynamics bias (opt-in; off by default)",
    enabled: function () { return !!(this.funnel && this.funnelOn); },
    energy: function (pos, f) {
      return this.funnel.addForces(pos, f);
    },
  },
];

installBuiltins("cg", CG_TERMS);

let _plan = null;
let _planGen = -1;

/**
 * The CG plan, memoised against TERM_GENERATION.
 *
 * The cache is a MODULE variable rather than a field on the force field on
 * purpose: tests/test_cg_class_shape.js pins ForceField to exactly 70 own
 * instance fields in a fixed order, so a plan cache on the instance would be a
 * visible shape change. It also does not need to be per-field — the registry
 * is process-global, so the plan is too.
 *
 * Per compute() this is one integer compare against a module-local `let` (a
 * monomorphic load) and, on a hit, nothing else.
 *
 * @returns {{terms: ReadonlyArray<object>, barrier: number}} `terms` ordered by
 *   `order`; `barrier` is the index of the first term with
 *   order >= HOLO_FLOOR_BARRIER_ORDER, i.e. where the inlined holo floor sits.
 */
export function cgPlan() {
  if (_planGen === TERM_GENERATION) return _plan;
  const terms = listTerms("cg");
  let barrier = terms.length;
  for (let i = 0; i < terms.length; i++) {
    if (terms[i].order >= HOLO_FLOOR_BARRIER_ORDER) { barrier = i; break; }
  }
  _plan = Object.freeze({ terms, barrier });
  _planGen = TERM_GENERATION;
  return _plan;
}

/**
 * The plan as plain data, for docs and for the analysis modules that want to
 * report which terms are live. NOT a hot-path function: it allocates a record
 * per term. `ff` is optional; without it the `active` flag cannot be resolved.
 * @param {object|null} ff
 * @returns {Array<{id:string,label:string,order:number,active:boolean}>}
 */
export function cgTermTable(ff = null) {
  return describeTerms("cg", ff).map((t) => ({
    id: t.id, label: t.label, order: t.order, active: t.active,
  }));
}
