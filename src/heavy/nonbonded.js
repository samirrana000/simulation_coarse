/**
 * heavy/nonbonded.js — the spatial-grid non-bonded kernels: LJ + Generalized
 * Born + screened Coulomb + directional H-bond, with and without GB, plus the
 * OBC-II variant.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10). These were the three biggest
 * methods of HeavyForceField and are attached to its prototype by
 * heavy/forcefield.js; see that file for why the descriptor is exactly a class
 * method descriptor.
 *
 *   _nonBondedGrid       the hot kernel. Every heavy compute() call runs this
 *                        (or the OBC2 variant) exactly once. Accumulator and
 *                        operand order are part of the result.
 *   _nonBondedGridNoGB   short-range-only variant used by the OBC2 path, so GB
 *                        is not double-counted
 *   _nonBondedGridOBC2   OBC-II Born radii + reaction field from
 *                        physics/solvation/gb_obc2.js, short range delegated
 *                        to _nonBondedGridNoGB
 *
 * FIXED 2026-10-03 — the temporal-dead-zone defect and the hybrid it produced
 *   _nonBondedGridNoGB used to read `ligStart` in its `trk` initializer, five
 *   lines BEFORE the `const nProt = ..., ligStart = ...` declaration that
 *   introduces it. With ff.trackTerms !== true the initializer short-circuits
 *   and never evaluates `ligStart`, so the ReferenceError was latent; with
 *   par.gbModel === "obc2" AND ff.trackTerms === true it threw, and compute()'s
 *   catch then ran the HCT fallback ON TOP of forces that _nonBondedGridOBC2
 *   had already merged — so the field reported pure-HCT energies with
 *   OBC2-GB + HCT-everything forces. A force field that is neither model, and
 *   neither energy nor forces agreed with it.
 *   Three changes, all required — reordering the declaration alone would only
 *   have silenced the throw and left the hybrid reachable for any FUTURE throw
 *   inside the short-range pass:
 *     1. `nProt` / `ligStart` are declared before their first use.
 *     2. _nonBondedGridOBC2 is TRANSACTIONAL: the OBC-II reaction field is
 *        accumulated into the scratch buffer `gbF` and merged into `f` only
 *        after _nonBondedGridNoGB has returned. A throw inside the short-range
 *        pass therefore leaves the OBC2 contribution unmerged.
 *     3. `this._obc2Radii` is published after the merge, so a failed pass does
 *        not leave Born radii on the field that no force term ever used.
 *   The energy.js catch additionally snapshots/restores `f` around each
 *   substitutable stage, so even a throw from inside the short-range loop
 *   cannot leave partial contributions behind; and every substitution now
 *   increments ff.physicsFallbacks / sets ff.lastPhysicsFallback instead of
 *   only writing to the console.
 *   tests/test_heavy_golden.js was regenerated from the fixed tree (the `full`
 *   configuration is the one that changed); see the commit message. The golden
 *   is a bit-exact net, NOT a tolerance, and was not loosened.
 *
 * Zero DOM globals. Node-importable.
 */
import { COULOMB_CONST } from "../physics/gb.js";
import { GB_RADII } from "../physics/charges.js";
import {
  computeBornRadii as computeOBC2Radii,
  gbEnergyForces as gbOBC2Forces,
  debyeKappa,
} from "../physics/solvation/gb_obc2.js";
import { R_CUT, switchFunc, switchDeriv } from "./params.js";
import { pairKey } from "./pairs.js";

/**
 * The three grid kernels, attached to HeavyForceField.prototype by
 * heavy/forcefield.js. Object-literal method shorthand: dynamic `this`, exactly
 * like a class method.
 */
export const nonBondedKernels = {
  /**
   * OBC2 non-bonded path: LJ + screened Coulomb + H-bonds on the spatial
   * grid (same as _nonBondedGrid) but GB reaction field from gb_obc2.js with
   * freshly computed OBC-II Born radii. Exclusions/1-4 scales honored.
   * @param {ArrayLike<number>} pos
   * @param {Float64Array} f
   */
  _nonBondedGridOBC2(pos, f) {
    const intrinsic = new Float64Array(this.n);
    for (let i = 0; i < this.n; i++) {
      const el = this.atoms[i]?.element ?? "C";
      intrinsic[i] = GB_RADII[el] ?? GB_RADII.DEFAULT ?? 1.6;
    }
    const born = computeOBC2Radii(pos, intrinsic, {});
    // The GB pair energy/forces go into a private scratch buffer, never
    // straight into the caller's `f`, so this kernel is transactional in `f`.
    const gbF = new Float64Array(this.n * 3);
    const kappa = debyeKappa(this.gbSaltM, 300, this.gbEpsOut);
    const gbRes = gbOBC2Forces(pos, this._charges, born, gbF, {
      epsIn: this.gbEpsIn, epsOut: this.gbEpsOut, kappa,
      excluded: this._excluded, scale14: this._scale14,
    });
    // TRANSACTIONAL ORDER: the short-range pass runs FIRST and the OBC-II
    // reaction field is merged into the shared force buffer only after it has
    // returned. _nonBondedGridNoGB writes into `f` pair by pair, so if it
    // throws the caller must not find a half-merged GB contribution already
    // sitting in `f` — energy.js's catch restores from a snapshot, and this
    // ordering is what makes the snapshot's restore exact rather than
    // best-effort. See the FIXED 2026-10-03 note in the header.
    // LJ + H-bond + solute-dielectric Coulomb, no GB double-count: the short-
    // range grid loop _nonBondedGridNoGB computes those three terms and skips
    // the HCT pairInteraction() entirely.
    const base = this._nonBondedGridNoGB(pos, f);
    for (let i = 0; i < f.length; i++) f[i] += gbF[i];
    this._obc2Radii = born; // published only once the whole pass has succeeded
    return {
      lj: base.lj, elec: base.elec, gb: gbRes.gbEnergy, hbond: base.hbond,
      // `* 0`: the OBC-II reaction field is NOT attributed to the protein-ligand
      // binding cross term. Unchanged by the 2026-10-03 fix and deliberately
      // left alone — turning it on is a physics DECISION (how much of the
      // reaction field counts as "binding"), not a defect repair, and flipping
      // it would silently change bindingU/desolvU for every OBC2 user.
      // Recorded here because it is a genuine asymmetry with the HCT path,
      // whose bindE at _nonBondedGrid DOES include gbRes.energy.
      bindE: base.bindE + gbRes.gbEnergy * 0,
      bindTerms: base.bindTerms, // S4 per-term trackers pass through (null when off)
    };
  },

  /**
   * Grid LJ + Coulomb(solute dielectric, NO GB reaction field) + H-bond.
   * Short-range half of the OBC2 branch, mirroring _nonBondedGrid minus
   * pairInteraction(). The OBC-II reaction field is added by
   * _nonBondedGridOBC2, AFTER this method returns — do not move that merge
   * above this call (see the transactional-order note there).
   * @param {ArrayLike<number>} pos
   * @param {Float64Array} f
   */
  _nonBondedGridNoGB(pos, f) {
    // The full loop is written out here rather than delegated: it has to skip
    // the HCT `gb.pairInteraction()` term entirely (no double-counting with the
    // OBC-II reaction field that _nonBondedGridOBC2 merges afterwards) and to
    // charge Coulomb to the solute dielectric, which pairInteraction() folds
    // together with the HCT descreening.
    const n = this.n;
    const nProt = this.nProt, ligStart = this.ligandStart;
    let ljTot = 0, elecTot = 0, hbondTot = 0, bindTot = 0;
    // Loop-2 S4: per-term trackers (skipped when trackTerms is off)
    const trk = this.trackTerms === true && ligStart > 0;
    let bLJ = 0, bCoul = 0, bHB = 0;
    this.grid.build(pos, n);
    const isDonor = this._hbClassification.isDonor;
    const isAcceptor = this._hbClassification.isAcceptor;
    this.grid.forEachPair(pos, n, R_CUT, (i, j, dx, dy, dz, r2, r) => {
      const k = pairKey(i, j);
      if (this._excluded.has(k)) return;
      const s14 = this._scale14.get(k) ?? 1.0;
      const ei = this._elem[i], ej = this._elem[j];
      const s = 0.5 * (ei.sigma + ej.sigma);
      const eps = Math.sqrt(ei.eps * ej.eps);
      const sr = s / r, sr6 = sr * sr * sr * sr * sr * sr;
      const ljE = 4 * eps * (sr6 * sr6 - sr6);
      const ljF = 4 * eps * (12 * sr6 * sr6 - 6 * sr6) / r;
      let hbE = 0, hbFx = 0, hbFy = 0, hbFz = 0;
      if ((isDonor[i] && isAcceptor[j]) || (isDonor[j] && isAcceptor[i])) {
        const hbRes = this.hbond.evaluatePair(i, j, dx, dy, dz, r);
        hbE = hbRes.energy; hbFx = hbRes.fx; hbFy = hbRes.fy; hbFz = hbRes.fz;
      }
      const S = switchFunc(r), dS = switchDeriv(r);
      const totRadialF = s14 * (-S * ljF + dS * ljE);
      ljTot += s14 * S * ljE;
      hbondTot += s14 * S * hbE;
      const xi = 3 * i, xj = 3 * j;
      const fx = (totRadialF * dx / r) + (s14 * S * hbFx);
      const fy = (totRadialF * dy / r) + (s14 * S * hbFy);
      const fz = (totRadialF * dz / r) + (s14 * S * hbFz);
      // Screened Coulomb (solute dielectric, no GB): 332 q_i q_j/(epsIn r)
      const qi = ei.q, qj = ej.q;
      if (qi !== 0 && qj !== 0) {
        const uC = ((COULOMB_CONST / this.gbEpsIn) * qi * qj * s14 * S) / r;
        elecTot += uC;
        // U = C·qq·S(r)/(eps·r); dU/dr = C·qq·(dS/r − S/r²)/eps; F_i = +dU/dr·dx/r
        const dUdrC = ((COULOMB_CONST / this.gbEpsIn) * qi * qj * s14 * dS) / r - uC / r;
        const fmagC = dUdrC / r;
        f[xi] += fmagC * dx; f[xi + 1] += fmagC * dy; f[xi + 2] += fmagC * dz;
        f[xj] -= fmagC * dx; f[xj + 1] -= fmagC * dy; f[xj + 2] -= fmagC * dz;
      }
      f[xi] += fx; f[xi + 1] += fy; f[xi + 2] += fz;
      f[xj] -= fx; f[xj + 1] -= fy; f[xj + 2] -= fz;
      if (i < nProt && j >= ligStart && ligStart > 0) {
        bindTot += s14 * S * (ljE + hbE);
        if (trk) {
          bLJ += s14 * S * ljE;
          bCoul += s14 * S * (((COULOMB_CONST / this.gbEpsIn) * qi * qj) / r);
          bHB += s14 * S * hbE;
        }
      }
    });
    return {
      lj: ljTot, elec: elecTot, gb: 0, hbond: hbondTot, bindE: bindTot,
      bindTerms: trk ? { lj: bLJ, coul: bCoul, hb: bHB } : null,
    };
  },

  /**
   * Fast O(N) Spatial Grid Non-Bonded Kernel.
   */
  _nonBondedGrid(pos, f) {
    const n = this.n;
    const nProt = this.nProt;
    const ligStart = this.ligandStart;
    let ljTot = 0, elecTot = 0, gbTot = 0, hbondTot = 0, bindTot = 0;
    // Loop-2 S4: per-term trackers (skipped when trackTerms is off — zero overhead)
    const trk = this.trackTerms === true && ligStart > 0;
    let bLJ = 0, bCoul = 0, bHB = 0;

    // Build spatial hash
    this.grid.build(pos, n);

    const isDonor = this._hbClassification.isDonor;
    const isAcceptor = this._hbClassification.isAcceptor;

    this.grid.forEachPair(pos, n, R_CUT, (i, j, dx, dy, dz, r2, r) => {
      const k = pairKey(i, j);
      if (this._excluded.has(k)) return;
      const s14 = this._scale14.get(k) ?? 1.0;

      const ei = this._elem[i];
      const ej = this._elem[j];
      const qi = ei.q;
      const qj = ej.q;

      // LJ
      const s = 0.5 * (ei.sigma + ej.sigma);
      const eps = Math.sqrt(ei.eps * ej.eps);
      const sr = s / r, sr6 = sr * sr * sr * sr * sr * sr;
      const ljE = 4 * eps * (sr6 * sr6 - sr6);
      const ljF = 4 * eps * (12 * sr6 * sr6 - 6 * sr6) / r;

      // Generalized Born + Screened Coulomb
      const gbRes = this.gb.pairInteraction(i, j, dx, dy, dz, r, qi, qj, this._bornRadii[i], this._bornRadii[j], s14);

      // Directional H-Bond
      let hbE = 0, hbFx = 0, hbFy = 0, hbFz = 0;
      if ((isDonor[i] && isAcceptor[j]) || (isDonor[j] && isAcceptor[i])) {
        const hbRes = this.hbond.evaluatePair(i, j, dx, dy, dz, r);
        hbE = hbRes.energy;
        hbFx = hbRes.fx;
        hbFy = hbRes.fy;
        hbFz = hbRes.fz;
      }

      // Smooth cutoff switch
      const S = switchFunc(r);
      const dS = switchDeriv(r);

      const totE = s14 * S * (ljE + hbE) + gbRes.energy;
      const totRadialF = s14 * (-S * ljF + dS * ljE);

      ljTot += s14 * S * ljE;
      elecTot += gbRes.coulombE;
      gbTot += gbRes.gbE;
      hbondTot += s14 * S * hbE;

      const xi = 3 * i, xj = 3 * j;
      const fx = (totRadialF * dx / r) + gbRes.fx + (s14 * S * hbFx);
      const fy = (totRadialF * dy / r) + gbRes.fy + (s14 * S * hbFy);
      const fz = (totRadialF * dz / r) + gbRes.fz + (s14 * S * hbFz);

      f[xi] += fx; f[xi + 1] += fy; f[xi + 2] += fz;
      f[xj] -= fx; f[xj + 1] -= fy; f[xj + 2] -= fz;

      if (i < nProt && j >= ligStart && ligStart > 0) {
        bindTot += totE;
        // S4 tracker: coul = full electrostatic pair term (Coulomb + GB
        // reaction field) so lj+coul+hb+desolv sums exactly to bindingU.
        if (trk) { bLJ += s14 * S * ljE; bCoul += gbRes.energy; bHB += s14 * S * hbE; }
      }
    });

    return {
      lj: ljTot, elec: elecTot, gb: gbTot, hbond: hbondTot, bindE: bindTot,
      bindTerms: trk ? { lj: bLJ, coul: bCoul, hb: bHB } : null,
    };
  },
};
