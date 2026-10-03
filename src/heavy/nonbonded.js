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
 * KNOWN PRE-EXISTING DEFECT, NOT INTRODUCED BY THE SPLIT and NOT FIXED HERE
 *   (a fix is a physics change, which this refactor is not allowed to be):
 *   _nonBondedGridNoGB reads `ligStart` in its `trk` initializer, five lines
 *   BEFORE the `const nProt = ..., ligStart = ...` declaration that introduces
 *   it. With ff.trackTerms !== true the initializer short-circuits and never
 *   evaluates `ligStart`, so the temporal-dead-zone ReferenceError is latent.
 *   With par.gbModel === "obc2" AND ff.trackTerms === true it throws, compute()
 *   catches it and falls back to _nonBondedGrid — but the OBC2 forces were
 *   already merged into the force buffer before the throw, so the returned
 *   forces are OBC2-GB plus HCT-everything while the returned energy is pure
 *   HCT. tests/test_heavy_golden.js pins those numbers on purpose: fixing the
 *   TDZ makes that test go red, which is the correct outcome for a physics
 *   change and cannot be mistaken for a refactor.
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
    this._obc2Radii = born;
    // GB pair energy/forces accumulated on a scratch buffer, then merged, so
    // the grid LJ/Coulomb pass below stays identical to the HCT path.
    const gbF = new Float64Array(this.n * 3);
    const kappa = debyeKappa(this.gbSaltM, 300, this.gbEpsOut);
    const gbRes = gbOBC2Forces(pos, this._charges, born, gbF, {
      epsIn: this.gbEpsIn, epsOut: this.gbEpsOut, kappa,
      excluded: this._excluded, scale14: this._scale14,
    });
    for (let i = 0; i < f.length; i++) f[i] += gbF[i];
    // LJ + H-bond (no GB double-count): reuse grid loop for short-range only.
    // To avoid duplicating the full kernel, call the legacy grid then subtract
    // its HCT GB contribution and add OBC2 instead.
    const base = this._nonBondedGridNoGB(pos, f);
    return {
      lj: base.lj, elec: base.elec, gb: gbRes.gbEnergy, hbond: base.hbond,
      bindE: base.bindE + gbRes.gbEnergy * 0,
      bindTerms: base.bindTerms, // S4 per-term trackers pass through (null when off)
    };
  },

  /**
   * Grid LJ + Coulomb(screened via GB pair coulomb part) + H-bond without GB.
   * Helper for the OBC2 branch; mirrors _nonBondedGrid minus pairInteraction GB.
   * @param {ArrayLike<number>} pos
   * @param {Float64Array} f
   */
  _nonBondedGridNoGB(pos, f) {
    // Delegate to the standard grid but with GB charges zeroed is wasteful;
    // instead run the standard grid and rely on _nonBondedGridOBC2 to have
    // already added OBC2 GB: here we compute LJ/H-bond/Coulomb-only by calling
    // _nonBondedGrid on a probe that skips GB via zero Born radii trick is not
    // clean, so implement the short-range loop directly (no GB term).
    const n = this.n;
    let ljTot = 0, elecTot = 0, hbondTot = 0, bindTot = 0;
    // Loop-2 S4: per-term trackers (skipped when trackTerms is off)
    const trk = this.trackTerms === true && ligStart > 0;
    let bLJ = 0, bCoul = 0, bHB = 0;
    this.grid.build(pos, n);
    const isDonor = this._hbClassification.isDonor;
    const isAcceptor = this._hbClassification.isAcceptor;
    const nProt = this.nProt, ligStart = this.ligandStart;
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
