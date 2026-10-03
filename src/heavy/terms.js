/**
 * heavy/terms.js — the heavy-atom term descriptors + the heavy plan cache.
 *
 * WHAT MOVED HERE
 * ---------------
 * Eight of the nine `U += …` stages of HeavyForceField.compute() are now
 * declarative descriptors here (the ninth, `heavy.weak`, lives with its kernel
 * in heavy/weak.js — see the note there). src/heavy/energy.js keeps the plan
 * walk, the binding-accounting block (stage 5a, which contributes neither U nor
 * forces), the NaN/∞ guard, and the three force-transaction primitives
 * (_forceSnapshot / _forceRestore / _notePhysicsFallback) the substitutable
 * stages call.
 *
 * SUMMATION-ORDER SAFETY OF THE MIGRATION (why the heavy golden did not move)
 * --------------------------------------------------------------------------
 * The heavy assembly is where the grouping rule bites. Two stages summed
 * SEVERAL terms in a single `+=`:
 *
 *     U += this.bondU + this.coordU;                        (stage 1)
 *     U += this.improperU + this.properU;                    (stage 3)
 *     U += nb.lj + nb.elec + nb.gb + nb.hbond;               (stage 4)
 *
 * Those are ONE descriptor each, and each returns that exact expression, so
 * the inner partial sums round identically and the outer `U +=` happens once,
 * in the same place. Splitting any of them into N descriptors would insert
 * N−1 extra roundings — that is the whole reason the registry interface is
 * "one float64 per term" and not "one float64 per physics component".
 *
 * For every other stage the descriptor performs the identical call sequence,
 * publishes the identical `ff.*U` accumulators in the identical order, and
 * returns the identical scalar, and the dispatcher adds it at the identical
 * position. Three details that are easy to get wrong and were checked here:
 *
 *   • `heavy.membrane` and `heavy.weak` keep their `else { … = 0 }` branches
 *     INSIDE the descriptor rather than expressing them as an `enabled`
 *     predicate. The heavy golden hashes `membraneU` / `weakU` / `piU` / `cpiU`
 *     / `xbU` for every configuration, including the four where those terms are
 *     off, so the zeroing must still happen — a skipped term would leave the
 *     previous step's value on the field.
 *   • `heavy.sasa` reads the non-bonded stage's `bindE` from `_nb` below, which
 *     the non-bonded descriptor assigns unconditionally. This replaces the
 *     loop-local `nb` the two stages used to share; it allocates nothing
 *     (`_nonBondedGrid` already returns a fresh record) and cannot be stale,
 *     because `heavy.nonbonded` has no `enabled` predicate and therefore always
 *     runs before it.
 *   • The stage-5a binding accounting block stayed in energy.js and now runs
 *     AFTER the whole plan. It writes only `ff.bindLJU / bindCoulU / bindHBU /
 *     desolvU / bindU / _bindSasaE` and reads only `nb.bindTerms` and
 *     `_bindSasaE`; no term reads any of those, so running it after the
 *     membrane / restraint / funnel stages instead of between SASA and membrane
 *     cannot change a number. Every one of the 24 golden accumulators is read
 *     after compute() returns, so only the final values matter.
 *
 * Zero DOM globals. Node-importable.
 */
import { installBuiltins, listTerms, TERM_GENERATION } from "../physics/terms/registry.js";
import { springForces, dihedralForcesAnalytic } from "../ff-harmonic.js";
import { lcpoSasa } from "../physics/solvation/lcpo_sasa.js";
import { membraneEnergyForces } from "../physics/solvation/membrane_slab.js";
import { enforceCoordination } from "../chem/metals.js";
import { harmonicFlat, angleFlat, harmonicFlatPerK, angleFlatPerK } from "./kernels.js";
import { weakTerm } from "./weak.js";

/**
 * The non-bonded stage's result record, published by `heavy.nonbonded` and
 * consumed by `heavy.sasa` — the two stages that the pre-registry loop shared
 * a loop-local `nb` between. Module-private rather than an instance field
 * because nothing outside this file reads it and adding a HeavyForceField
 * instance field is a visible shape change. See the header's second bullet.
 */
let _nb = null;

/**
 * The non-bonded stage's per-term binding trackers (S4 / `trackTerms`), or null
 * when tracking is off. Read once by the stage-5a accounting block in
 * heavy/energy.js — the same value the pre-registry loop read off its local
 * `nb`. Not a hot-path function: the block runs only when `trackTerms === true`.
 * @returns {{lj:number,coul:number,hb:number}|null}
 */
export function bindingTrackers() {
  return _nb === null ? null : _nb.bindTerms;
}

/**
 * The heavy terms, in the order they contribute to U.
 *
 * `energy` is `function (pos, f)` with `this` = the HeavyForceField, matching
 * every kernel in this engine. None allocates per call.
 */
export const HEAVY_TERMS = [
  {
    // Stage 1. Bonds (uniform or AMBER14 per-bond k) PLUS metal coordination
    // (legacy k=40 distance springs, upgraded to radial+angular by
    // enforceCoordination for metals with a detected geometry). The two are ONE
    // term because the original summed them in one `+=`.
    id: "heavy.covalent", engine: "heavy", order: 10,
    label: "Covalent bonds (uniform or AMBER14 per-bond k) + metal coordination",
    reports: ["bondU", "coordU"],
    report: function () { return this.bondU + this.coordU; },
    energy: function (pos, f) {
      // AMBER14 per-bond k when opted in
      if (this.useAmber14 && this._bondK && this._bondK.length === this.bonds.length / 3) {
        this.bondU = harmonicFlatPerK(pos, f, this.bonds, 3, this._bondK);
      } else {
        this.bondU = harmonicFlat(pos, f, this.bonds, 3, this.kBond);
      }
      // Metal coordination: legacy k=40 distance springs for metals WITHOUT a
      // detected coordination geometry. When par.metalAngles === true (R3 §1e),
      // metals WITH a geometry are handled by enforceCoordination (radial
      // k=40 + cross-angle k=20) instead of these springs.
      const me = this._metalEnforce;
      this.coordU = 0;
      if (this.coord.length > 0) {
        this.coordU = harmonicFlat(pos, f, this.coord, 3, this.metalK);
      }
      if (me) {
        const sub = me.metals.filter((m) => me.hasGeometry.has(m.index));
        if (sub.length > 0) {
          const res = enforceCoordination(pos, sub, me.elements, {
            forces: f, kRadial: this.metalK, kAngle: 20.0, ideal: true,
          });
          this.coordU += res.energy;
        }
      }
      return this.bondU + this.coordU;
    },
  },
  {
    id: "heavy.angles", engine: "heavy", order: 20,
    label: "Angle bending (uniform or AMBER14 per-angle k)",
    reports: ["angleU"],
    report: function () { return this.angleU; },
    energy: function (pos, f) {
      if (this.useAmber14 && this._angleK && this._angleK.length === this.angles.length / 4) {
        this.angleU = angleFlatPerK(pos, f, this.angles, 4, this._angleK);
      } else {
        this.angleU = angleFlat(pos, f, this.angles, 4, this.kAngle);
      }
      return this.angleU;
    },
  },
  {
    id: "heavy.dihedrals", engine: "heavy", order: 30,
    label: "Proper + improper dihedrals (analytic Blondel–Karplus / Bekker gradients)",
    reports: ["improperU", "properU"],
    report: function () { return this.improperU + this.properU; },
    energy: function (pos, f) {
      this.improperU = dihedralForcesAnalytic(pos, f, this.impropers, 5, this.kImproper);
      this.properU = dihedralForcesAnalytic(pos, f, this.propers, 5, this.kProper);
      return this.improperU + this.properU;
    },
  },
  {
    // Stage 4. Spatial-grid LJ + Generalized Born + screened Coulomb + H-bond.
    // gbModel "obc2" routes the GB reaction field through gb_obc2.js with OBC-II
    // radii (LJ/Coulomb/H-bond stay on the grid kernel); default "hct" preserves
    // the legacy GeneralizedBorn path.
    //
    // The try/catch is a physics substitution, not an error swallow: OBC2 and
    // HCT are DIFFERENT force fields, so falling back changes the answer. The
    // snapshot makes the substitution self-consistent (the substitute sees the
    // pre-OBC2 force state, so forces are pure HCT to match the pure-HCT
    // energy) and _notePhysicsFallback makes it visible on the field.
    id: "heavy.nonbonded", engine: "heavy", order: 40,
    label: "Spatial-grid non-bonded: LJ + GB/screened Coulomb + directional H-bond",
    reports: ["repU", "elecU", "gbU", "hbondU"],
    report: function () { return this.repU + this.elecU + this.gbU + this.hbondU; },
    energy: function (pos, f) {
      let nb;
      if (this.gbModel === "obc2") {
        const snap = this._forceSnapshot();
        try {
          nb = this._nonBondedGridOBC2(pos, f);
        } catch (e) {
          this._forceRestore(snap);
          this._notePhysicsFallback("OBC2 non-bonded", "HCT non-bonded", e);
          nb = this._nonBondedGrid(pos, f);
        }
      } else {
        nb = this._nonBondedGrid(pos, f);
      }
      this.elecU = nb.elec;
      this.repU = nb.lj;
      this.gbU = nb.gb;
      this.hbondU = nb.hbond;
      _nb = nb;
      return nb.lj + nb.elec + nb.gb + nb.hbond;
    },
  },
  weakTerm,
  {
    // Stage 5. Hydrophobic SASA burial; sasaModel "lcpo" routes through
    // lcpo_sasa.js, anything else through the legacy SasaModel.
    id: "heavy.sasa", engine: "heavy", order: 60,
    label: "Hydrophobic SASA burial (LCPO or legacy per-atom SASA)",
    reports: ["sasaU", "bindingU"],
    report: function () { return this.sasaU; },
    energy: function (pos, f) {
      if (this.sasaModel === "lcpo") {
        const snap = this._forceSnapshot();
        try {
          const res = lcpoSasa(pos, this._lcpoElements, {
            gamma: this.sasa.gamma, excluded: this._excluded, forces: f,
          });
          this.sasaU = res.energy;
          this.bindingU = _nb.bindE;
          if (this.trackTerms === true) this._bindSasaE = 0; // LCPO: bindE carries no SASA part
          return this.sasaU;
        } catch (e) {
          // Same defect class as the OBC2 catch above: lcpoSasa accumulates into
          // `f` per pair, so a mid-pass throw would leave burial forces behind for
          // the legacy SASA kernel to double up on. Roll back first, then
          // substitute. (Unreachable today — lcpo_sasa.js has no throw statement —
          // but the invariant is what makes that fact unimportant.)
          this._forceRestore(snap);
          this._notePhysicsFallback("LCPO SASA", "legacy SASA", e);
          const sasaRes = this.sasa.compute(pos, f, this._elem, this.n, this.nProt, this.ligandStart);
          this.sasaU = sasaRes.energy;
          this.bindingU = _nb.bindE + sasaRes.bindSasaE;
          if (this.trackTerms === true) this._bindSasaE = sasaRes.bindSasaE;
          return this.sasaU;
        }
      }
      const sasaRes = this.sasa.compute(pos, f, this._elem, this.n, this.nProt, this.ligandStart);
      this.sasaU = sasaRes.energy;
      this.bindingU = _nb.bindE + sasaRes.bindSasaE;
      if (this.trackTerms === true) this._bindSasaE = sasaRes.bindSasaE;
      return this.sasaU;
    },
  },
  {
    // Stage 5b. Implicit membrane slab, opt-in via par.membrane = {on:true,…}.
    id: "heavy.membrane", engine: "heavy", order: 70,
    label: "Implicit membrane slab (opt-in via par.membrane.on)",
    reports: ["membraneU"],
    report: function () { return this.membraneU; },
    energy: function (pos, f) {
      if (!this.membraneOpts?.on) {
        this.membraneU = 0;
        return 0;
      }
      const snap = this._forceSnapshot();
      try {
        const radii = new Float64Array(this.n);
        for (let i = 0; i < this.n; i++) radii[i] = 2.0;
        const mem = membraneEnergyForces(pos, this._charges, radii, this._lcpoElements, f, {
          thickness: this.membraneOpts.thickness ?? 15,
          width: this.membraneOpts.width ?? 2,
          epsWater: this.membraneOpts.epsWater ?? this.gbEpsOut ?? 78.5,
          epsMem: this.membraneOpts.epsMem ?? 2.0,
          zCenter: this.membraneOpts.zCenter ?? 0,
        });
        this.membraneU = mem.energy;
        return mem.energy;
      } catch (e) {
        // "skipped" is only true if the slab's forces are removed too — a
        // mid-loop throw in membraneEnergyForces leaves per-atom z-forces in `f`
        // that no energy term accounts for. Roll back, then skip.
        this._forceRestore(snap);
        this._notePhysicsFallback("membrane slab", "no membrane term", e);
        this.membraneU = 0;
        return 0;
      }
    },
  },
  {
    id: "heavy.restraints", engine: "heavy", order: 80,
    label: "ML contact restraints (active springs list)",
    reports: ["springU"],
    report: function () { return this.springU; },
    energy: function (pos, f) {
      this.springU = this.springs.length ? springForces(this, pos, f) : 0;
      return this.springU;
    },
  },
  {
    id: "heavy.funnel", engine: "heavy", order: 90,
    label: "Funnel bias (opt-in; off by default)",
    enabled: function () { return !!(this.funnel && this.funnelOn); },
    energy: function (pos, f) {
      return this.funnel.addForces(pos, f);
    },
  },
];

installBuiltins("heavy", HEAVY_TERMS);

let _plan = null;
let _planGen = -1;

/**
 * The heavy plan, memoised against TERM_GENERATION. One integer compare per
 * compute() against a module-local `let`; nothing else on a hit. See the CG twin
 * in src/cg/terms.js for why this is a module variable and not an instance field.
 * @returns {ReadonlyArray<object>} terms ordered by `order`
 */
export function heavyPlan() {
  if (_planGen === TERM_GENERATION) return _plan;
  _plan = listTerms("heavy");
  _planGen = TERM_GENERATION;
  return _plan;
}
