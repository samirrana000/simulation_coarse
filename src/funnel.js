/**
 * funnel.js — binding funnel + well-tempered metadynamics for reconstructing
 * the protein–ligand binding free-energy profile (PMF) along a 1-D collective
 * variable: the distance r between the ligand centre of mass and the (fixed)
 * centre of mass of the protein pocket.
 *
 * Funnel bias (style inspired by funnel-MD / CG-FMD): a collective-variable
 * restraint that pulls the ligand into the pocket at large separations but is
 * EXACTLY FLAT inside the bound state (r ≤ rFlat), so the bound-state
 * equilibrium is undisturbed:
 *     U_wall(r) = ½ kf (r − rFlat)²   for r > rFlat,
 *              = 0                     otherwise,
 * with dU_wall/dr = kf (r − rFlat). Because only the ligand COM is coupled and
 * the force is split equally over the ligand atoms, the restraint translates
 * the ligand as a rigid body — no torque, so molecular orientation is free
 * inside the pocket.
 *
 * Well-tempered metadynamics (Barducci, Bussi, Parrinello 2008): every
 * hillStride force calls a Gaussian hill is deposited on the bias grid along
 * r. The hill height is damped with the accumulated bias,
 *     w = w0 · exp(−V(r) / (k_B · T · (γ − 1))),    γ = biasFactor ≥ 1,
 * so the bias grows like the free energy and never over-fills the basin (γ = 1
 * recovers plain metadynamics, w = w0). The bias grids store both the value
 * V(r_k) and its exact r-derivative dV/dr at each bin centre, so the bias
 * force −(dV/dr)·û follows the analytic gradient of the deposited hill sum.
 *
 * Reconstruction: the negative bias is the free energy profile,
 *     F(r) ≈ −V_bias(r) + const,
 * shifted so the minimum inside the bound-state window is 0, and the binding
 * free energy is read off the bias difference between the bound state and a
 * far-apart reference: ΔG_bind ≈ V_bias(r_far) − V_bias(r_bound).
 *
 * The system is a unified particle set (units Å, ps, kcal/mol, Da): protein
 * Cα beads occupy global indices 0..nProt−1, ligand heavy atoms nProt..n−1.
 * A module with no ligand atoms degrades to a safe no-op (active === false).
 */

import { KB_KCAL } from "./forcefield.js?v=8";

export class Funnel {
  /**
   * @param {object} opts
   * @param {number} opts.nProt      protein Cα bead count (indices 0..nProt−1)
   * @param {number} opts.n          total particle count (protein + ligand atoms)
   * @param {Float64Array} opts.ref  native coordinates, length 3n (never mutated)
   * @param {number} [opts.rPocket=8.0]   pocket radius (Å): protein beads within
   *   this distance of the native ligand COM define the pocket
   * @param {number} [opts.rFlat=5.0]     funnel wall onset (Å); CV is flat inside
   * @param {number} [opts.kf=2.0]        funnel wall spring, kcal/mol/Å²
   * @param {number} [opts.sigma=0.3]     metadynamics hill width (Å)
   * @param {number} [opts.w0=0.02]       initial hill height (kcal/mol)
   * @param {number} [opts.biasFactor=6.0] well-tempered γ (1 ⇒ plain metadynamics)
   * @param {number} [opts.hillStride=20]  deposit one hill every this-many
   *   addForces calls
   * @param {number} [opts.rMax=24.0]     CV range covered by the PMF grid (Å)
   * @param {number} [opts.bins=96]       PMF grid resolution
   */
  constructor({ nProt, n, ref, rPocket = 8.0, rFlat = 5.0, kf = 2.0, sigma = 0.3,
    w0 = 0.02, biasFactor = 6.0, hillStride = 20, rMax = 24.0, bins = 96 }) {
    this.active = true;
    this.nProt = nProt;
    this.n = n;
    this.nLig = n - nProt;       // ligand heavy atoms = global indices nProt..n−1
    this.rPocket = rPocket;
    this.rFlat = rFlat;
    this.kf = kf;
    this.sigma = sigma;
    this.w0 = w0;
    this.biasFactor = biasFactor;
    this.hillStride = hillStride;
    this.rMax = rMax;
    this.bins = bins;
    this.ref = ref;

    this.T = 300;                 // K, WTM temperature (setTemperature)
    this.on = true;               // bias/funnel enable switch (setOn)
    this.calls = 0;               // addForces call counter (hill schedule)
    this._lastCV = 0;
    this._nHills = 0;             // deposited hills (estimateDG needs ≥ 1)
    this.pmf = null;
    this.pmfShift = 0;
    this._scratch = new Float64Array(3); // reused ligand-COM buffer (GC-free)

    this._rStep = rMax / bins;    // grid spacing (bin k centre at (k+0.5)·rStep)

    // No ligand atoms ⇒ nothing for the funnel to act on: every method below
    // becomes a safe no-op so callers never need to branch.
    if (this.nLig <= 0) {
      this.active = false;
      return;
    }

    // Native ligand COM — the centre of the bound pose from the reference.
    let lx = 0, ly = 0, lz = 0;
    for (let a = 0; a < this.nLig; a++) {
      const c = 3 * (nProt + a);
      lx += ref[c]; ly += ref[c + 1]; lz += ref[c + 2];
    }
    this.ligCOM0 = [lx / this.nLig, ly / this.nLig, lz / this.nLig];

    // Pocket: every protein bead whose NATIVE distance to the ligand COM is
    // within rPocket. Its COM is fixed at construction — the funnel and the
    // PMF are defined against an immobile reference pocket.
    const pocketIdx = [];
    for (let i = 0; i < nProt; i++) {
      const c = 3 * i;
      const dx = ref[c] - this.ligCOM0[0], dy = ref[c + 1] - this.ligCOM0[1], dz = ref[c + 2] - this.ligCOM0[2];
      if (Math.sqrt(dx * dx + dy * dy + dz * dz) <= rPocket) pocketIdx.push(i);
    }
    // Degenerate pocket fallback: the 6 nearest protein beads (all if fewer).
    if (pocketIdx.length === 0) {
      const byDist = [];
      for (let i = 0; i < nProt; i++) {
        const c = 3 * i;
        const dx = ref[c] - this.ligCOM0[0], dy = ref[c + 1] - this.ligCOM0[1], dz = ref[c + 2] - this.ligCOM0[2];
        byDist.push([dx * dx + dy * dy + dz * dz, i]);
      }
      byDist.sort((a, b) => a[0] - b[0]);
      for (let k = 0; k < Math.min(6, nProt); k++) pocketIdx.push(byDist[k][1]);
    }
    this.pocket = pocketIdx;
    let px = 0, py = 0, pz = 0;
    for (const i of pocketIdx) {
      const c = 3 * i;
      px += ref[c]; py += ref[c + 1]; pz += ref[c + 2];
    }
    this.pocketCOM = [px / pocketIdx.length, py / pocketIdx.length, pz / pocketIdx.length];

    this.cv0 = this._dist3(this.ligCOM0, this.pocketCOM); // native (bound) CV

    // Well-tempered metadynamics bias grids over r ∈ [0, rMax].
    //   _bias[k]      — deposited bias V at bin centre r_k
    //   _biasForce[k] — analytic dV/dr at r_k (derivative of the hill sum)
    this._bias = new Float64Array(bins);
    this._biasForce = new Float64Array(bins);
  }

  /** |a − b| for two {0,1,2} triplets. */
  _dist3(a, b) {
    const dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  /** Instantaneous ligand COM (mean over ligand global indices nProt..n−1). */
  _ligMean(pos) {
    const nProt = this.nProt, nLig = this.nLig, s = this._scratch;
    let x = 0, y = 0, z = 0;
    for (let a = 0; a < nLig; a++) {
      const c = 3 * (nProt + a);
      x += pos[c]; y += pos[c + 1]; z += pos[c + 2];
    }
    s[0] = x / nLig; s[1] = y / nLig; s[2] = z / nLig;
    return s;
  }

  /**
   * Collective variable: distance from the instantaneous ligand COM to the
   * fixed pocket COM.
   * @param {Float64Array} pos  positions, length 3n
   * @returns {number} r in Å
   */
  cv(pos) {
    if (!this.active) return 0;
    return this._dist3(this._ligMean(pos), this.pocketCOM);
  }

  /**
   * Linear interpolation of a bias grid at coordinate r (bins hold centre
   * samples r_k = (k+0.5)·rStep; out-of-range r clamps to the edge bin).
   */
  _interp(arr, r) {
    const nb = this.bins;
    const x = r / this._rStep - 0.5;
    if (x < 0) return arr[0];
    const k0 = Math.floor(x);
    if (k0 >= nb - 1) return arr[nb - 1];
    const f = x - k0;
    return arr[k0] + f * (arr[k0 + 1] - arr[k0]);
  }

  /**
   * Apply the funnel wall + well-tempered bias forces to f and return their
   * potential energy U_bias (kcal/mol). Self-deposits one hill every
   * hillStride calls.
   *
   * Sign convention: with û = (ligMean − pocketCOM)/r the radial outward unit
   * vector, a radial potential U(r) exerts on the ligand COM the force
   * −(dU/dr)·û. Splitting that force equally over the nLig atoms translates
   * the COM only — the correct action of a COM restraint (no torque, the
   * ligand never gets spun around).
   *
   * @param {Float64Array} pos  positions, length 3n
   * @param {Float64Array} f    force accumulator, length 3n (kcal/mol/Å);
   *   funnel/bias forces are added in place
   * @returns {number} U_wall + V_bias (kcal/mol); 0 when disabled/inactive
   */
  addForces(pos, f) {
    if (!this.active || !this.on) return 0;
    const m = this._ligMean(pos);
    const r = this._dist3(m, this.pocketCOM);
    this._lastCV = r;
    const nLig = this.nLig;
    // Radial unit vector pocketCOM → ligand COM (û = ∇r).
    const ux = (m[0] - this.pocketCOM[0]) / r;
    const uy = (m[1] - this.pocketCOM[1]) / r;
    const uz = (m[2] - this.pocketCOM[2]) / r;

    let U = 0;

    // Funnel wall — harmonic pull toward the pocket beyond rFlat, exactly
    // flat (no force) inside the bound state so the bound equilibrium is
    // undisturbed. dU/dr = kf·(r − rFlat) is the radial force magnitude.
    let dUdr = 0;
    if (r > this.rFlat) {
      const dr = r - this.rFlat;
      U += 0.5 * this.kf * dr * dr;
      dUdr = this.kf * dr;
    }

    // Well-tempered metadynamics bias along r: value for the energy, and the
    // interpolated analytic derivative for the radial force.
    U += this._interp(this._bias, r);
    dUdr += this._interp(this._biasForce, r);

    // Force on the ligand COM = −(dU/dr)·û, split equally over its atoms.
    const F = -dUdr / nLig;
    const fx = F * ux, fy = F * uy, fz = F * uz;
    for (let a = 0; a < nLig; a++) {
      const c = 3 * (this.nProt + a);
      f[c] += fx; f[c + 1] += fy; f[c + 2] += fz;
    }

    // Deposition schedule: one hill per hillStride force evaluations.
    this.calls++;
    if (this.calls % this.hillStride === 0) this.deposit(r);

    return U;
  }

  /**
   * Add a well-tempered Gaussian hill of width σ centred at CV = r to the
   * bias grid, including its exact derivative. Public so an external driver
   * can deposit hills on any schedule.
   *
   *   g(r') = w · exp(−(r'−r)² / 2σ²),
   *   w     = w0 · exp(−V(r) / (k_B·T·(γ−1)))   (γ = biasFactor; 1 ⇒ w = w0)
   *   dg/dr'= −(r'−r)/σ² · g(r')
   *
   * @param {number} r  hill centre (CV value, Å)
   */
  deposit(r) {
    if (!this.active) return;
    // Well-tempered damping: hills shrink exponentially where bias already
    // accumulated, so V converges to a scaled copy of the free energy and
    // never keeps growing (plain metadynamics when γ = 1).
    const w = this.biasFactor === 1
      ? this.w0
      : this.w0 * Math.exp(-this._interp(this._bias, r) / (KB_KCAL * this.T * (this.biasFactor - 1)));
    const inv2sig2 = 1 / (2 * this.sigma * this.sigma);
    const invSig2 = 1 / (this.sigma * this.sigma);
    for (let k = 0; k < this.bins; k++) {
      const rk = (k + 0.5) * this._rStep;
      const dr = rk - r;
      const g = w * Math.exp(-dr * dr * inv2sig2);
      this._bias[k] += g;
      this._biasForce[k] += -dr * invSig2 * g; // dg/dr at bin centre
    }
    this._nHills++;
  }

  /**
   * Reconstructed binding PMF: pmf = −V_bias, shifted so its minimum over the
   * bound-state window r ∈ [0, rFlat·1.5] is 0. Returns copies.
   * @returns {{r: Float64Array, pmf: Float64Array}} both of length bins
   */
  getPMF() {
    const r = new Float64Array(this.bins);
    const pmf = new Float64Array(this.bins);
    if (!this.active) {
      for (let k = 0; k < this.bins; k++) r[k] = (k + 0.5) * this._rStep;
      this.pmfShift = 0;
      return { r, pmf };
    }
    const shiftLimit = this.rFlat * 1.5;
    let min = Infinity;
    for (let k = 0; k < this.bins; k++) {
      const rk = (k + 0.5) * this._rStep;
      r[k] = rk;
      const v = -this._bias[k];
      pmf[k] = v;
      if (rk <= shiftLimit && v < min) min = v;
    }
    this.pmfShift = min === Infinity ? 0 : min;
    for (let k = 0; k < this.bins; k++) pmf[k] -= this.pmfShift;
    return { r, pmf };
  }

  /**
   * Binding free energy from the WTM bias difference between the bound state
   * (native CV cv0) and a far reference (cv0 + 6 Å, clamped into the grid):
   *     ΔG_bind ≈ V_bias(r_far) − V_bias(r_bound)
   * (kcal/mol; negative ⇒ binding favorable). Requires ≥ 1 deposited hill,
   * else NaN. Note: the exact well-tempered reconstruction carries the
   * (γ/(γ−1)) scaling factor on the bias difference; this estimator returns
   * the raw bias difference as requested.
   * @returns {number} kcal/mol
   */
  estimateDG() {
    if (!this.active || this._nHills === 0) return NaN;
    const rBound = this.cv0;
    const rFar = Math.min(this.cv0 + 6, this.rMax - this._rStep);
    return this._interp(this._bias, rFar) - this._interp(this._bias, rBound);
  }

  /** Clear hills, bias grids, PMF and the deposition counter. */
  reset() {
    if (this.active) {
      this._bias.fill(0);
      this._biasForce.fill(0);
    }
    this._nHills = 0;
    this.calls = 0;
    this.pmf = null;
    this.pmfShift = 0;
    this._lastCV = 0;
  }

  /** Set the temperature (K) used for the well-tempered height scaling. */
  setTemperature(T) { this.T = T; }

  /** Enable/disable the funnel + bias forces (default on). */
  setOn(on) { this.on = on; }

  /** Most recent CV value computed by addForces(). */
  get lastCV() { return this._lastCV; }
}
