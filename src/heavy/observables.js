/**
 * heavy/observables.js — heavy-side fold metrics.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10). Attached to
 * HeavyForceField.prototype by heavy/forcefield.js.
 *
 *   kineticTemp  T = ke_kcal/(1.5 n k_B). Units, formula and the mass-layout
 *                 contract are documented once, in physics/observables.js; this
 *                 is a delegating shim, not a second implementation.
 *   rmsd          protein heavy atoms only (comparable with CG's rmsd, and
 *                 ligand undocking does not read as fold instability)
 *   rmsdLig       the external-ligand block [ligandStart, n), no alignment
 *   rmsdAll       the legacy all-atom value, kept for callers/tests that need it
 *
 * Zero DOM globals. Node-importable.
 */
import { kineticTemp, rmsdTo } from "../physics/observables.js";

/** Attached to HeavyForceField.prototype by heavy/forcefield.js. */
export const observableMethods = {
  /**
   * Instantaneous kinetic temperature, T = ke_kcal/(1.5·n·k_B) — units, formula
   * and the mass-layout contract are documented once, in physics/observables.js.
   * D1: heavy used to infer the layout from the array length and fall back to
   * this.masses; CG assumed flat 3n. The resolution now happens HERE and the
   * layout is passed as an explicit `mode` — the kernel never guesses. Values
   * are unchanged on both engines.
   */
  kineticTemp(vel, mass) {
    const m = (mass && mass.length === this.n * 3) ? mass : (this.masses || mass);
    return kineticTemp(vel, m, this.n, m && m.length === this.n * 3 ? "dof" : "atom");
  },

  /**
   * RMSD to native over the protein heavy atoms only (fold-stability metric).
   * Mirrors ForceField.rmsd (protein-only): ligand/hetero drift is excluded so
   * HUD RMSD is comparable across CG/heavy and ligand undocking does not read as
   * fold instability. No alignment — ENM keeps the COM/orientation nearly fixed.
   * Use rmsdLig() for the ligand part, rmsdAll() for the legacy all-atom value.
   * Formula: rmsdTo() in physics/observables.js.
   */
  rmsd(pos) {
    const nP = this.nProt;
    if (!Number.isFinite(nP) || nP <= 0) return this.rmsdAll(pos);
    return rmsdTo(pos, this.ref, this._maskProt, nP);
  },

  /**
   * Ligand-only RMSD over the external-ligand block [ligandStart, n), no
   * alignment. Returns 0 when no ligand is present. Hetero/cofactor atoms
   * (protein..ligandStart) are excluded — same slice the HUD ligRMSD uses.
   */
  rmsdLig(pos) {
    const nL = this.nLigAtoms;
    if (!Number.isFinite(nL) || nL <= 0) return 0;
    return rmsdTo(pos, this.ref, this._maskLig, nL);
  },

  /**
   * Legacy all-atom RMSD (protein + hetero + ligand). Preserved for backward
   * compatibility with callers/tests that need the old heavy rmsd number.
   */
  rmsdAll(pos) { return rmsdTo(pos, this.ref, null, this.n); }
};
