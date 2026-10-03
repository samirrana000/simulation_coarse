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
 * Zero DOM globals. Node-importable.
 */
import { springForces, dihedralForcesAnalytic } from "../ff-harmonic.js";
import { lcpoSasa } from "../physics/solvation/lcpo_sasa.js";
import { membraneEnergyForces } from "../physics/solvation/membrane_slab.js";
import { enforceCoordination } from "../chem/metals.js";
import { harmonicFlat, angleFlat, harmonicFlatPerK, angleFlatPerK } from "./kernels.js";
import { METAL_K } from "./params.js";

/** Attached to HeavyForceField.prototype by heavy/forcefield.js. */
export const energyMethods = {
  /**
   * Total potential energy and forces evaluation.
   *
   * G65 ZERO-ALLOC AUDIT — hot loop (called every integration step):
   *   Expected allocs per compute(): ideally 0, but CURRENTLY 3 × Float64Array(n) via SasaModel.compute()
   *   (s0, radii, burial — see src/physics/sasa.js:47-54). These allocate O(n) each call and cause
   *   measurable heap growth (~ n*8*3 bytes per step). Suggested fix (non-breaking): hoist s0/radii/burial
   *   to HeavyForceField scratch buffers (e.g., this._sasaS0, this._sasaRadii, this._sasaBurial) allocated
   *   once at construction and reused — same pattern as ForceField's _dens/_bp* buffers. Until then,
   *   ~zero-alloc claim does not hold for heavy mode; see bench/alloc.js for heap delta measurement.
   *   Remaining kernels (harmonicFlat, angleFlat, dihedralForcesAnalytic, _nonBondedGrid via
   *   SpatialGrid.head/next/cellCoords) are zero-alloc in steady state (grid resizes only when n exceeds maxAtoms).
   */

  compute(pos) {
    const f = this.forces;
    f.fill(0);
    let U = 0;

    // 1. Covalent bonds + metal coordination (AMBER14 per-bond k when opted in)
    if (this.useAmber14 && this._bondK && this._bondK.length === this.bonds.length / 3) {
      this.bondU = harmonicFlatPerK(pos, f, this.bonds, 3, this._bondK);
    } else {
      this.bondU = harmonicFlat(pos, f, this.bonds, 3, this.kBond);
    }
    // Metal coordination: legacy k=40 distance springs for metals WITHOUT a
    // detected coordination geometry. When par.metalAngles === true (R3 §1e),
    // metals WITH a geometry are handled below by enforceCoordination (radial
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
    U += this.bondU + this.coordU;

    // 2. Angles (AMBER14 per-angle k when opted in)
    if (this.useAmber14 && this._angleK && this._angleK.length === this.angles.length / 4) {
      this.angleU = angleFlatPerK(pos, f, this.angles, 4, this._angleK);
    } else {
      this.angleU = angleFlat(pos, f, this.angles, 4, this.kAngle);
    }
    U += this.angleU;

    // 3. Fast Analytic Proper & Improper Dihedrals
    this.improperU = dihedralForcesAnalytic(pos, f, this.impropers, 5, this.kImproper);
    this.properU = dihedralForcesAnalytic(pos, f, this.propers, 5, this.kProper);
    U += this.improperU + this.properU;

    // 4. Fast Spatial-Grid Non-Bonded (LJ + Generalized Born + Screened Coulomb + H-bonds)
    // gbModel "obc2" routes the GB reaction field through gb_obc2.js with
    // OBC-II radii (LJ/Coulomb/H-bond stay on the grid kernel); default "hct"
    // preserves the legacy GeneralizedBorn path so existing tests pass.
    let nb;
    if (this.gbModel === "obc2") {
      try {
        nb = this._nonBondedGridOBC2(pos, f);
      } catch (e) {
        console.warn(`[HeavyForceField] OBC2 path failed (${e.message}) — HCT fallback`);
        nb = this._nonBondedGrid(pos, f);
      }
    } else {
      nb = this._nonBondedGrid(pos, f);
    }
    this.elecU = nb.elec;
    this.repU = nb.lj;
    this.gbU = nb.gb;
    this.hbondU = nb.hbond;
    U += nb.lj + nb.elec + nb.gb + nb.hbond;

    // 4b. Weak interactions (opt-in par.weak === "on", Loop-2 S3 / R3 §1a–c):
    // π-stack ring-ring, cation-π cation-ring, halogen σ-hole X···acceptor.
    // Runs nonbonded-adjacent (after the grid pass); pairs already excluded
    // by _excluded (1-2/1-3/intra-ligand) or same-ring are skipped.
    if (this.weakOn) {
      const w = this._weakInteractions(pos, f);
      this.weakU = w.pi + w.cpi + w.xb;
      this.piU = w.pi; this.cpiU = w.cpi; this.xbU = w.xb;
      U += this.weakU;
    } else {
      this.weakU = 0; this.piU = 0; this.cpiU = 0; this.xbU = 0;
    }

    // 5. Hydrophobic SASA burial ("lcpo" routes through lcpo_sasa.js)
    if (this.sasaModel === "lcpo") {
      try {
        const res = lcpoSasa(pos, this._lcpoElements, { gamma: this.sasa.gamma, excluded: this._excluded, forces: f });
        this.sasaU = res.energy;
        this.bindingU = nb.bindE;
        if (this.trackTerms === true) this._bindSasaE = 0; // LCPO: bindE carries no SASA part
        U += this.sasaU;
      } catch (e) {
        console.warn(`[HeavyForceField] LCPO SASA failed (${e.message}) — legacy SASA fallback`);
        const sasaRes = this.sasa.compute(pos, f, this._elem, this.n, this.nProt, this.ligandStart);
        this.sasaU = sasaRes.energy;
        this.bindingU = nb.bindE + sasaRes.bindSasaE;
        if (this.trackTerms === true) this._bindSasaE = sasaRes.bindSasaE;
        U += this.sasaU;
      }
    } else {
      const sasaRes = this.sasa.compute(pos, f, this._elem, this.n, this.nProt, this.ligandStart);
      this.sasaU = sasaRes.energy;
      this.bindingU = nb.bindE + sasaRes.bindSasaE;
      if (this.trackTerms === true) this._bindSasaE = sasaRes.bindSasaE;
      U += this.sasaU;
    }

    // 5a. Loop-2 S4 per-term binding accumulators (R4 §5 item 1, R6 §5).
    // trackTerms === true splits bindingU into {lj, coul, hb, desolv} from
    // the per-pair trackers filled in _nonBondedGrid/_nonBondedGridNoGB
    // (bindTerms) + the SASA cross-burial part captured above (bindSasaE);
    // the S3 weak terms join via piU/cpiU/xbU (bindU vector). DEFAULT OFF —
    // trk branches in the kernels are skipped, path bit-identical to pre-S4.
    if (this.trackTerms === true) {
      const bt = nb.bindTerms || {};
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

    // 5b. Implicit membrane slab (opt-in via par.membrane = {on:true,...})
    if (this.membraneOpts?.on) {
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
        U += mem.energy;
      } catch (e) {
        console.warn(`[HeavyForceField] membrane slab failed (${e.message}) — skipped`);
        this.membraneU = 0;
      }
    } else {
      this.membraneU = 0;
    }

    // 6. ML contact restraints (if active)
    this.springU = this.springs.length ? springForces(this, pos, f) : 0;
    U += this.springU;

    // 7. Funnel bias (if active)
    if (this.funnel && this.funnelOn) U += this.funnel.addForces(pos, f);

    if (!Number.isFinite(U)) {
      for (let i = 0; i < f.length; i++) if (!Number.isFinite(f[i])) f[i] = 0;
      this._nanStrikes++;
      U = NaN;
    }
    this.energy = U;
    return U;
  },
};
