/**
 * heavy/energy.js — compute(): the total-energy assembly order.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10). Attached to
 * HeavyForceField.prototype by heavy/forcefield.js.
 *
 * The term ORDER and the running total U here are the contract: bonds ->
 * metal coordination -> angles -> propers/impropers -> non-bonded grid ->
 * weak pass -> SASA -> membrane -> ML contact springs -> funnel bias -> NaN
 * guard. Reordering terms changes the last bits of the returned energy and the
 * per-term accumulators that the HUD and BindLog read, so it is a numeric
 * change, not a style change. tests/test_heavy_golden.js pins all 24
 * accumulators bit-exactly for this reason.
 *
 * SUBSTITUTION POLICY (2026-10-03)
 *   Three stages here are wrapped in try/catch that runs a DIFFERENT physics
 *   kernel on failure: OBC2 -> HCT non-bonded, LCPO -> legacy SASA, membrane
 *   slab -> omitted. Those are force-field substitutions, so they must never be
 *   silent and must never leave two models' forces in one buffer. Hence:
 *     • _forceSnapshot() / _forceRestore() bracket each stage, so a substitute
 *       kernel sees the pre-stage force state and its forces exactly match its
 *       energy.
 *     • _notePhysicsFallback() prints, increments ff.physicsFallbacks and sets
 *       ff.lastPhysicsFallback {stage, substitutedTo, error, count}.
 *   See the FIXED 2026-10-03 note in heavy/nonbonded.js for the OBC2 case that
 *   motivated this: it silently produced pure-HCT energy with OBC2-GB +
 *   HCT-everything forces.
 *
 * THE TERM REGISTRY (2026-10) — WHERE THE STAGES WENT
 * ---------------------------------------------------
 * The nine `U += …` stages are no longer written in this file. They are
 * declarative descriptors in src/heavy/terms.js (and, for `heavy.weak`, in
 * src/heavy/weak.js next to its kernel), validated by
 * src/physics/terms/registry.js, and compute() below walks the resulting plan:
 *
 *     const plan = heavyPlan();                 // memoised; 1 int compare
 *     for (let t = 0; t < plan.length; t++) {
 *       const d = plan[t];
 *       const en = d.enabled;
 *       if (en !== null && !en.call(this)) continue;
 *       U += d.energy.call(this, pos, f);
 *     }
 *
 * Each iteration performs the identical call sequence, publishes the identical
 * ff.*U accumulators, and adds the identical float64 at the identical position —
 * including the three stages that sum SEVERAL terms in one `+=`, which stay one
 * descriptor each precisely so the inner partial sums round the same way. The
 * full per-stage summation-order argument is in src/heavy/terms.js's header, and
 * tests/test_heavy_golden.js is the proof: 6 configurations x 12 poses x 24
 * accumulators, bit-exact, before and after this move.
 *
 * What deliberately stayed here: the three force-transaction primitives the
 * substitutable stages call, stage 5a (the per-term binding accounting, which
 * contributes neither U nor forces), and the NaN/∞ guard — none of which is a
 * term.
 *
 * Zero DOM globals. Node-importable.
 */
import { heavyPlan, bindingTrackers } from "./terms.js";

/** Attached to HeavyForceField.prototype by heavy/forcefield.js. */
export const energyMethods = {
  /**
   * Copy the live force buffer into the reusable snapshot slot.
   *
   * WHY THIS EXISTS — every substitutable physics stage in compute() is
   * wrapped in a try/catch that runs a DIFFERENT kernel on failure. A kernel
   * that throws part-way through has already added some of its contributions to
   * `f`, so the substitute then adds its own on top: the caller ends up with
   * forces from BOTH models and an energy from only one. That is not a
   * degraded answer, it is a wrong one, and it was exactly what happened to
   * the OBC2 path (see the FIXED 2026-10-03 note in heavy/nonbonded.js).
   *
   * Cost: lazily allocated on first use and REUSED, so the default HCT path
   * never allocates and the guarded paths allocate once for the life of the
   * field rather than once per step — the same hoist-to-constructor-buffer
   * pattern the G65 audit note above asks for. It is deliberately NOT
   * initialised in the constructor: two of the three guarded stages are
   * off by default, and paying for a buffer nobody uses is exactly the
   * allocation churn this audit is about.
   * @returns {Float64Array} snapshot buffer, already filled from `this.forces`
   */
  _forceSnapshot() {
    const f = this.forces;
    let snap = this._forceSnap;
    if (snap === null || snap === undefined || snap.length !== f.length) {
      snap = this._forceSnap = new Float64Array(f.length);
    }
    snap.set(f);
    return snap;
  },

  /**
   * Undo every force contribution written since the matching _forceSnapshot(),
   * so the substitute kernel starts from exactly the pre-stage state.
   * @param {Float64Array} snap  buffer returned by _forceSnapshot()
   */
  _forceRestore(snap) {
    const f = this.forces;
    if (snap.length !== f.length) return;
    f.set(snap);
  },

  /**
   * Record that a physics path was substituted for another, LOUDLY.
   *
   * Silent substitution of a force field is the failure mode this project's
   * honest-scope discipline exists to prevent: the run continues, the HUD
   * shows a plausible number, and nothing in the reported energies says which
   * model produced them. So every substitution now (a) prints, (b) counts
   * itself on the field as `physicsFallbacks`, and (c) leaves a queryable
   * `lastPhysicsFallback` record naming the stage and the original error. A
   * caller that never looks at these still sees the console; a caller that
   * does look can refuse the result.
   *
   * @param {string} stage    which compute() stage was substituted
   * @param {string} substitutedTo  short name of the kernel now being used
   * @param {unknown} e       the thrown value
   */
  _notePhysicsFallback(stage, substitutedTo, e) {
    const message = e && e.message ? e.message : String(e);
    this.physicsFallbacks = (this.physicsFallbacks ?? 0) + 1;
    this.lastPhysicsFallback = {
      stage, substitutedTo, error: message, count: this.physicsFallbacks,
    };
    console.warn(
      `[HeavyForceField] ${stage} failed (${message}) — substituted ${substitutedTo}. ` +
      `Forces rolled back to the pre-stage state so they match ${substitutedTo} exactly; ` +
      `ff.physicsFallbacks=${this.physicsFallbacks}. Results are NOT the model requested.`,
    );
  },

  /**
   * Total potential energy and forces evaluation.
   *
   * THE BODY IS A PLAN WALK, NOT A SUMMATION SEQUENCE (2026-10)
   * ----------------------------------------------------------
   * The stage-by-stage `U += …` sequence this method used to contain now lives
   * in src/heavy/terms.js as declarative descriptors; this file owns the walk,
   * the two blocks that are not terms, and nothing else. See the header.
   *
   * G65 ZERO-ALLOC AUDIT — hot loop (called every integration step):
   *   Expected allocs per compute(): ideally 0, but CURRENTLY 3 × Float64Array(n) via SasaModel.compute()
   *   (s0, radii, burial — see src/physics/sasa.js:47-54). These allocate O(n) each call and cause
   *   measurable heap growth (~ n*8*3 bytes per step). Suggested fix (non-breaking): hoist s0/radii/burial
   *   to HeavyForceField scratch buffers (e.g., this._sasaS0, this._sasaRadii, this._sasaBurial) allocated
   *   once at construction and reused — same pattern as ForceField's _dens/_bp* buffers. Until then,
   *   ~zero-alloc claim does not hold for heavy mode; see bench/alloc.js for heap delta measurement.
   *   _forceSnapshot() above follows that same recommended pattern for the
   *   guarded stages (lazily allocated + reused, so steady-state alloc count is
   *   unchanged). The registry walk adds NOTHING to this budget: heavyPlan() is
   *   a memoised frozen array and a term returns a number, never a record.
   *   Remaining kernels (harmonicFlat, angleFlat, dihedralForcesAnalytic,
   *   _nonBondedGrid via SpatialGrid.head/next/cellCoords) are zero-alloc in steady
   *   state (grid resizes only when n exceeds maxAtoms).
   */

  compute(pos) {
    const f = this.forces;
    f.fill(0);
    let U = 0;

    // Stages 1-7, in the order they contribute to U. Each descriptor performs
    // the identical call sequence and returns the identical float64 the
    // inlined `U += …` statement produced, at the identical position; a
    // descriptor whose `enabled` predicate is false is skipped entirely, so U
    // is not touched — which is what the original `if (…) U += …` guards did.
    const plan = heavyPlan();
    const nTerms = plan.length;
    for (let t = 0; t < nTerms; t++) {
      const d = plan[t];
      const en = d.enabled;
      if (en !== null && !en.call(this)) continue;
      U += d.energy.call(this, pos, f);
    }

    // 5a. Loop-2 S4 per-term binding accumulators (R4 §5 item 1, R6 §5).
    // trackTerms === true splits bindingU into {lj, coul, hb, desolv} from
    // the per-pair trackers filled in _nonBondedGrid/_nonBondedGridNoGB
    // (bindTerms) + the SASA cross-burial part captured by the sasa term
    // (_bindSasaE); the S3 weak terms join via piU/cpiU/xbU (bindU vector).
    // DEFAULT OFF — trk branches in the kernels are skipped, path bit-identical
    // to pre-S4.
    //
    // NOT A TERM: it contributes no U and no force, only accounting, and it
    // needs two stages' results at once. It runs AFTER the whole plan instead
    // of between SASA and membrane; nothing the membrane / restraint / funnel
    // stages touch is read here, so the move cannot change a number. See the
    // third bullet in src/heavy/terms.js's header.
    if (this.trackTerms === true) {
      const bt = bindingTrackers() || {};
      this.bindLJU = bt.lj || 0;
      this.bindCoulU = bt.coul || 0;
      this.bindHBU = bt.hb || 0;
      this.desolvU = this._bindSasaE || 0;
      // Weak cross terms: the S3 kernels are whole-molecule; the ligand's
      // share enters the binding vector via the native-pose contacts each
      // term makes (approximated here by the weak totals when a ligand is
      // present — documented approximation, S5 consumes only the split).
      this.bindU = {
        lj: this.bindLJU, coul: this.bindCoulU, hb: this.bindHBU,
        desolv: this.desolvU, pi: this.piU || 0, cpi: this.cpiU || 0, xb: this.xbU || 0,
      };
    }
    this._bindSasaE = null;

    // NaN/∞ guard — not a term either. Sanitize the force buffer and leave
    // ff.energy = NaN so the app can auto-pause on the next frame.
    if (!Number.isFinite(U)) {
      for (let i = 0; i < f.length; i++) if (!Number.isFinite(f[i])) f[i] = 0;
      this._nanStrikes++;
      U = NaN;
    }
    this.energy = U;
    return U;
  },
};
