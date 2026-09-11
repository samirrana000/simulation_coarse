/**
 * respa.js — r-RESPA multiple-time-step (MTS) Langevin integrator.
 *
 * Reversible reference-system propagator (Tuckerman, Berne & Martyna,
 * J. Chem. Phys. 97, 1990 (1992)): split the force into a fast part
 * (stiff bonded terms: bonds, angles, torsions, metal coordination) evaluated
 * every inner step h = 1 fs, and a slow part (soft non-bonded terms: LJ,
 * Coulomb, GB-OBC reaction field, SASA, membrane, restraint springs, funnel
 * bias) evaluated every outer step H = nInner * h (2-4 fs).
 *
 * Outer propagator per step H (impulse-MTS + BAOAB-compatible thermostat):
 *   1. B_slow : v += (H/2) * KCONV * F_slow / m      (cached slow impulse)
 *   2. nInner x inner velocity-Verlet on F_fast:
 *        v += (h/2) KCONV F_fast/m; x += h v;
 *        F_fast = bondedFn(x);      v += (h/2) KCONV F_fast/m
 *   3. O      : exact Ornstein-Uhlenbeck thermostat at the OUTER step
 *        v <- c1 v + c2 sqrt(kB T KCONV / m) xi, c1 = exp(-zeta H)
 *        (same stationary Maxwell-Boltzmann variance as integrator.js BAOAB)
 *   4. F_slow = nonbondedFn(x)                       (one slow eval per H)
 *   5. B_slow : v += (H/2) * KCONV * F_slow_new / m
 *
 * Resonance note: impulse-MTS is stable while H << pi/omega_fast; with
 * omega_fast ~ 100 ps^-1 (light-atom bonds) H = 4 fs keeps H*omega ~ 0.4.
 * The outer OU thermostat additionally damps resonance growth, which is why
 * the O step lives outside the inner loop (BAOAB-compatible placement).
 *
 * Units: A, ps, kcal/mol, Da (KCONV/KB_KCAL contract, see src/units.js).
 * Zero new npm deps. Vanilla ES module. All slow-path failures must be
 * handled by the caller (parity fallback to single-step BAOAB).
 *
 * @module physics/integrators/respa
 */

import { KB_KCAL, KCONV } from "../../units.js?v=10";
import {
  harmonicPairs,
  springForces,
  angleForces,
  dihedralForcesAnalytic,
} from "../../ff-harmonic.js?v=10";
import { lcpoSasa } from "../solvation/lcpo_sasa.js?v=10";
import { membraneEnergyForces } from "../solvation/membrane_slab.js?v=10";

/** Per-particle |dv| clamp per kick (A/ps). Mirrors integrator.js _kick. */
export const RESPA_DVMAX = 2.0;
/** Default inner step: 1 fs bonded. */
export const RESPA_DT_INNER = 0.001;
/** Default outer step: 4 fs non-bonded. */
export const RESPA_DT_OUTER = 0.004;

/**
 * r-RESPA stepper. Works on caller-owned Float64Array buffers; never
 * allocates in the stepping loop after construction (scratch reuse).
 */
export class RESPAStepper {
  /**
   * @param {object} opts
   * @param {Float64Array} opts.masses per-particle masses (Da), length n
   * @param {number} [opts.temperature=300] bath T (K)
   * @param {number} [opts.friction=5.0] Langevin zeta (ps^-1)
   * @param {number} [opts.dtInner=0.001] inner step h (ps, bonded)
   * @param {number} [opts.dtOuter=0.004] outer step H (ps, non-bonded)
   */
  constructor(opts = {}) {
    const masses = opts.masses;
    if (!masses || !masses.length) throw new Error("RESPAStepper: opts.masses (per-particle, Da) required");
    this.n = masses.length;
    this.n3 = this.n * 3;
    this.mass = new Float64Array(this.n);
    for (let i = 0; i < this.n; i++) this.mass[i] = masses[i] > 0 ? masses[i] : 12.0;
    this.invMass = new Float64Array(this.n);
    this.thermal = new Float64Array(this.n);
    this.T = opts.temperature ?? 300;
    this.zeta = opts.friction ?? 5.0;
    this.dtInner = opts.dtInner ?? RESPA_DT_INNER;
    this.dtOuter = opts.dtOuter ?? RESPA_DT_OUTER;
    this.nInner = Math.max(1, Math.round(this.dtOuter / this.dtInner));
    // Snap H to an integer multiple of h so inner/outer stay commensurate.
    this.dtOuter = this.nInner * this.dtInner;
    this._rebuildThermal();

    // Scratch (zero-alloc stepping): cached slow/fast forces + RNG buffer.
    this._slowF = new Float64Array(this.n3);
    this._fastF = new Float64Array(this.n3);
    this._slowReady = false;
    this._fastReady = false;
    this._rnd = new Float64Array(this.n3);
    this._haveSpare = false;
    this._spare = 0;

    this.fastU = 0;
    this.slowU = 0;
    this.time = 0;
    this.steps = 0; // outer steps completed
    this.innerSteps = 0;
  }

  /** Build from any force field exposing per-particle .masses (CG or heavy). */
  static fromForceField(ff, opts = {}) {
    if (!ff || !ff.masses || !ff.n) throw new Error("RESPAStepper.fromForceField: ff with .masses/.n required");
    const m = new Float64Array(ff.n);
    if (ff.masses.length === ff.n) {
      m.set(ff.masses);
    } else {
      for (let i = 0; i < ff.n; i++) m[i] = ff.masses[3 * i];
    }
    return new RESPAStepper({ masses: m, ...opts });
  }

  /** Set bath temperature (K); refreshes the per-particle thermal scale. */
  setTemperature(T) { this.T = Math.max(1, T); this._rebuildThermal(); }
  /** Set friction zeta (ps^-1). */
  setFriction(z) { this.zeta = Math.max(0.1, z); }

  /**
   * Set the inner/outer steps (ps). Outer is snapped to nInner * inner.
   * @param {number} dtInner bonded step (ps)
   * @param {number} dtOuter non-bonded step (ps)
   */
  setSteps(dtInner, dtOuter) {
    this.dtInner = dtInner;
    this.dtOuter = dtOuter;
    this.nInner = Math.max(1, Math.round(dtOuter / dtInner));
    this.dtOuter = this.nInner * this.dtInner;
    this.invalidateCache();
  }

  /** Drop cached forces (call after externally moving positions). */
  invalidateCache() { this._slowReady = false; this._fastReady = false; }

  _rebuildThermal() {
    for (let i = 0; i < this.n; i++) {
      this.invMass[i] = 1 / this.mass[i];
      this.thermal[i] = Math.sqrt(KB_KCAL * this.T * KCONV * this.invMass[i]);
    }
  }

  /**
   * Advance one outer step H.
   *
   * @param {Float64Array} positions flat 3n, mutated in place
   * @param {Float64Array} velocities flat 3n, mutated in place
   * @param {Float64Array} forces flat 3n scratch, overwritten with F_fast + F_slow
   * @param {(pos: Float64Array, fOut: Float64Array) => number} bondedFn
   *   fast forces; must OVERWRITE fOut (zero-fill first) and return energy
   * @param {(pos: Float64Array, fOut: Float64Array) => number} nonbondedFn
   *   slow forces; same contract
   * @returns {{fastU: number, slowU: number, time: number}} energies (kcal/mol)
   */
  step(positions, velocities, forces, bondedFn, nonbondedFn) {
    const H = this.dtOuter, h = this.dtInner;
    const { n } = this;

    if (!this._slowReady) {
      this._slowF.fill(0);
      this.slowU = nonbondedFn(positions, this._slowF);
      this._slowReady = true;
    }
    if (!this._fastReady) {
      this._fastF.fill(0);
      this.fastU = bondedFn(positions, this._fastF);
      this._fastReady = true;
    }

    // 1. Outer half-kick with cached slow forces.
    this._kick(velocities, this._slowF, H / 2);

    // 2. Inner velocity-Verlet loop on fast forces.
    for (let k = 0; k < this.nInner; k++) {
      this._kick(velocities, this._fastF, h / 2);
      for (let i = 0; i < this.n3; i++) positions[i] += h * velocities[i];
      this._fastF.fill(0);
      this.fastU = bondedFn(positions, this._fastF);
      this._kick(velocities, this._fastF, h / 2);
      this.innerSteps++;
    }
    this._fastReady = true;

    // 3. Outer exact-OU thermostat (BAOAB-compatible O step at stride H).
    const c1 = Math.exp(-this.zeta * H);
    const c2 = Math.sqrt(Math.max(0, 1 - c1 * c1));
    this._fillGaussian(this.n3);
    for (let p = 0; p < n; p++) {
      const sig = c2 * this.thermal[p];
      velocities[3 * p] = c1 * velocities[3 * p] + sig * this._rnd[3 * p];
      velocities[3 * p + 1] = c1 * velocities[3 * p + 1] + sig * this._rnd[3 * p + 1];
      velocities[3 * p + 2] = c1 * velocities[3 * p + 2] + sig * this._rnd[3 * p + 2];
    }

    // 4-5. Fresh slow forces + closing half-kick.
    this._slowF.fill(0);
    this.slowU = nonbondedFn(positions, this._slowF);
    this._slowReady = true;
    this._kick(velocities, this._slowF, H / 2);

    // Observer buffer: total forces (callers read ff-style forces here).
    for (let i = 0; i < this.n3; i++) forces[i] = this._fastF[i] + this._slowF[i];

    // NaN guard (same policy as integrator.js: sanitize, caller auto-pauses).
    if (!Number.isFinite(this.fastU + this.slowU)) {
      for (let i = 0; i < this.n3; i++) if (!Number.isFinite(forces[i])) forces[i] = 0;
    }

    this.time += H;
    this.steps++;
    return { fastU: this.fastU, slowU: this.slowU, time: this.time };
  }

  /**
   * Half-kick with per-particle |dv| clamp (integrator.js _kick parity).
   * @param {Float64Array} vel flat 3n
   * @param {Float64Array} F flat 3n forces (kcal/mol/A)
   * @param {number} hdt half step (ps)
   */
  _kick(vel, F, hdt) {
    const { invMass } = this;
    for (let p = 0; p < this.n; p++) {
      const im = KCONV * hdt * invMass[p];
      const i3 = 3 * p;
      let dvx = im * F[i3], dvy = im * F[i3 + 1], dvz = im * F[i3 + 2];
      if (!Number.isFinite(dvx)) dvx = 0;
      if (!Number.isFinite(dvy)) dvy = 0;
      if (!Number.isFinite(dvz)) dvz = 0;
      const m2 = dvx * dvx + dvy * dvy + dvz * dvz;
      if (m2 > RESPA_DVMAX * RESPA_DVMAX) {
        const sc = RESPA_DVMAX / Math.sqrt(m2);
        dvx *= sc; dvy *= sc; dvz *= sc;
      }
      vel[i3] += dvx; vel[i3 + 1] += dvy; vel[i3 + 2] += dvz;
    }
  }

  /** Fill this._rnd[0..k) with N(0,1) via Box-Muller (paired deviates). */
  _fillGaussian(k) {
    let i = 0;
    if (this._haveSpare) { this._rnd[0] = this._spare; this._haveSpare = false; i = 1; }
    while (i < k) {
      let u = 0, v = 0, s = 0;
      do { u = Math.random() * 2 - 1; v = Math.random() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
      const mul = Math.sqrt((-2 * Math.log(s)) / s);
      this._rnd[i] = u * mul;
      if (i + 1 < k) this._rnd[i + 1] = v * mul;
      else { this._spare = v * mul; this._haveSpare = true; }
      i += 2;
    }
  }
}

// ---------------------------------------------------------------------------
// Force-splitting adapters (exact mirrors of ForceField/HeavyForceField terms)
// ---------------------------------------------------------------------------

/**
 * Split a coarse-grained ForceField into RESPA fast/slow callbacks.
 * Fast (inner, 1 fs): peptide bonds, ENM springs, holo + native contacts,
 *   holo compression floor, backbone angles, ligand internal terms.
 * Slow (outer, 2-4 fs): excluded-volume repulsion, protein-ligand binding,
 *   funnel/metadynamics bias. Mirrors ForceField.compute() sections 1-5c;
 *   keep in sync with src/forcefield.js if that kernel changes.
 * @param {object} ff ForceField instance (ff.heavy falsy)
 * @returns {{bondedFn: Function, nonbondedFn: Function}}
 */
export function splitCoarsegrained(ff) {
  if (!ff || ff.heavy) throw new Error("splitCoarsegrained: coarse-grained ForceField required");
  const bondedFn = (pos, fOut) => {
    fOut.fill(0);
    let U = 0;
    U += ff._harmonicPairs(pos, fOut, ff.bonds, 3, ff.kBond);
    U += ff._springForces(pos, fOut);
    if (ff.holoSprings.length) U += ff._harmonicPairs(pos, fOut, ff.holoSprings, 3, ff.holoGamma);
    if (ff.nativeContacts.length) U += ff._harmonicPairs(pos, fOut, ff.nativeContacts, 3, 1.0);
    if (ff.holoSprings.length) {
      const HS = ff.holoSprings, RMIN = 2.6, KF = 8.0;
      for (let a = 0; a < HS.length; a += 3) {
        const i = 3 * HS[a], j = 3 * HS[a + 1];
        const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
        if (r >= RMIN) continue;
        const dr = r - RMIN;
        U += 0.5 * KF * dr * dr;
        const s = (KF * dr) / r;
        fOut[i] += s * dx; fOut[i + 1] += s * dy; fOut[i + 2] += s * dz;
        fOut[j] -= s * dx; fOut[j + 1] -= s * dy; fOut[j + 2] -= s * dz;
      }
    }
    U += ff._angleForces(pos, fOut);
    U += ff._ligandInternal(pos, fOut);
    return U;
  };
  const nonbondedFn = (pos, fOut) => {
    fOut.fill(0);
    let U = 0;
    U += ff._repulsion(pos, fOut);
    U += ff._binding(pos, fOut);
    if (ff.funnel && ff.funnelOn) U += ff.funnel.addForces(pos, fOut);
    return U;
  };
  return { bondedFn, nonbondedFn };
}

/**
 * Split a HeavyForceField into RESPA fast/slow callbacks.
 * Fast (inner, 1 fs): covalent bonds, metal coordination, angles, proper +
 *   improper torsions (uniform-k kernels bit-identical to heavy.js
 *   harmonicFlat/angleFlat/dihedralForcesAnalytic calls in compute()).
 * Slow (outer, 2-4 fs): grid non-bonded (HCT or OBC2 branch, same call as
 *   compute()), SASA (legacy or LCPO branch), membrane slab, ML springs,
 *   funnel bias. Mirrors HeavyForceField.compute() sections 1-7.
 * Throws for useAmber14 (per-bond/per-angle stiffness tables are private to
 * heavy.js) so the caller can parity-fallback to single-step BAOAB.
 * @param {object} ff HeavyForceField instance (ff.heavy truthy)
 * @returns {{bondedFn: Function, nonbondedFn: Function}}
 */
export function splitHeavy(ff) {
  if (!ff || !ff.heavy) throw new Error("splitHeavy: HeavyForceField required");
  if (ff.useAmber14) {
    throw new Error("splitHeavy: useAmber14 per-bond stiffness has no split kernel — use single-step BAOAB");
  }
  const bondedFn = (pos, fOut) => {
    fOut.fill(0);
    let U = 0;
    U += harmonicPairs(pos, fOut, ff.bonds, 3, ff.kBond);
    U += harmonicPairs(pos, fOut, ff.coord, 3, ff.metalK);
    U += angleForces(ff, pos, fOut, ff.angles, ff.kAngle);
    U += dihedralForcesAnalytic(pos, fOut, ff.impropers, 5, ff.kImproper);
    U += dihedralForcesAnalytic(pos, fOut, ff.propers, 5, ff.kProper);
    return U;
  };
  const nonbondedFn = (pos, fOut) => {
    fOut.fill(0);
    let U = 0;
    let nb;
    if (ff.gbModel === "obc2") {
      try {
        nb = ff._nonBondedGridOBC2(pos, fOut);
      } catch (e) {
        console.warn(`[splitHeavy] OBC2 path failed (${e.message}) — HCT fallback`);
        nb = ff._nonBondedGrid(pos, fOut);
      }
    } else {
      nb = ff._nonBondedGrid(pos, fOut);
    }
    U += nb.lj + nb.elec + nb.gb + nb.hbond;
    if (ff.sasaModel === "lcpo") {
      try {
        const res = lcpoSasa(pos, ff._lcpoElements, {
          gamma: ff.sasa.gamma, excluded: ff._excluded, forces: fOut,
        });
        U += res.energy;
      } catch (e) {
        console.warn(`[splitHeavy] LCPO SASA failed (${e.message}) — legacy SASA fallback`);
        U += ff.sasa.compute(pos, fOut, ff._elem, ff.n, ff.nProt, ff.ligandStart).energy;
      }
    } else {
      U += ff.sasa.compute(pos, fOut, ff._elem, ff.n, ff.nProt, ff.ligandStart).energy;
    }
    if (ff.membraneOpts?.on) {
      try {
        const radii = new Float64Array(ff.n).fill(2.0);
        U += membraneEnergyForces(pos, ff._charges, radii, ff._lcpoElements, fOut, {
          thickness: ff.membraneOpts.thickness ?? 15,
          width: ff.membraneOpts.width ?? 2,
          epsWater: ff.membraneOpts.epsWater ?? ff.gbEpsOut ?? 78.5,
          epsMem: ff.membraneOpts.epsMem ?? 2.0,
          zCenter: ff.membraneOpts.zCenter ?? 0,
        }).energy;
      } catch (e) {
        console.warn(`[splitHeavy] membrane slab failed (${e.message}) — skipped`);
      }
    }
    // ML contact restraints (per-spring k) + funnel live on the slow clock.
    // Exact: same springForces(ff, pos, f) call as HeavyForceField.compute().
    if (ff.springs && ff.springs.length) U += springForces(ff, pos, fOut);
    if (ff.funnel && ff.funnelOn) U += ff.funnel.addForces(pos, fOut);
    return U;
  };
  return { bondedFn, nonbondedFn };
}

/**
 * Auto-split any supported force field (CG or heavy) for RESPAStepper.
 * @param {object} ff ForceField or HeavyForceField instance
 * @returns {{bondedFn: Function, nonbondedFn: Function}}
 */
export function splitForceField(ff) {
  if (!ff) throw new Error("splitForceField: force field required");
  return ff.heavy ? splitHeavy(ff) : splitCoarsegrained(ff);
}
