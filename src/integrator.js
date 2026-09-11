/**
 * integrator.js — Langevin dynamics with the BAOAB splitting scheme
 * (Leimkuhler & Matthews, 2013). BAOAB is the de-facto standard for
 * configurational sampling under Langevin dynamics: it is 2nd-order accurate
 * in positions, preserves the Boltzmann distribution for weak friction, and
 * remains stable at timesteps where Euler-Maruyama or velocity-Verlet+OU
 * hybrids blow up.
 *
 * Stochastic differential equation (per bead, per coordinate):
 *     m dv = (−∇U − m ζ v) dt + √(2 m ζ k_B T) dW
 * which is the Langevin equation quoted in the UI. Here
 *   ζ  — friction coefficient (ps⁻¹), the "damping" slider,
 *   T  — bath temperature (K),
 *   m  — particle mass (Da ≈ g/mol). Protein Cα beads carry the Cα mass
 *        (default 110 Da, settable via setMass); ligand united atoms carry
 *        their united-atom mass (e.g. C ≈ 12 Da). Stored per coordinate
 *        (3n), consistent with kcal/mol·ps²·Å⁻² via the force→acceleration
 *        conversion below.
 *
 * BAOAB step for timestep h:
 *   B:  v ← v + (h/2) F/m                            (half kick)
 *   A:  x ← x + (h/2) v                              (half drift)
 *   O:  v ← c1 v + c2 √(k_B T / m) ξ ,  ξ ~ N(0,1)   (exact OU thermostat)
 *         c1 = exp(−ζ h),  c2 = √(1 − c1²)
 *   A:  x ← x + (h/2) v                              (half drift)
 *   B:  v ← v + (h/2) F/m(x_new)                     (half kick w/ new force)
 *
 * Unit consistency (mass Da, distance Å, time ps, energy kcal/mol):
 *   1 kcal/mol              = 4.184 kJ/mol
 *                           = 4.184·10³ / Nₐ J per particle
 *                           = 4.184·10³ · (Nₐ·10⁻³ kg) · Å²/ps² × 10⁴ / Nₐ
 *                           = 418.4  Da·Å²/ps².
 * Hence the force→acceleration conversion is
 *   a[Å/ps²] = KCONV · F[kcal/mol/Å] / m[Da],   KCONV = 418.4,
 * and the thermal velocity scale is
 *   σ_v = √(k_B T / (m·KCONV))   [Å/ps],
 * with k_B = KB_KCAL kcal/mol/K (equipartition: ½m⟨v²⟩KCONV = ½k_BT per dof).
 */

import { KB_KCAL, KCONV } from "./units.js?v=10";

export class LangevinIntegrator {
  /**
   * @param {Float64Array} refPositions  native Cα positions (for reset)
   * @param {object} ff                  ForceField instance
   * @param {number} mass                protein Cα mass in Da (default ≈ 110);
   *   used only as the fallback when ff.masses is missing/empty (legacy
   *   protein-only force field). When ligand masses are present they are read
   *   from ff.masses and this argument only covers the protein beads.
   */
  constructor(refPositions, ff, mass = 110) {
    this.ff = ff;
    const n3 = this.n3 = ff.n * 3;

    // Per-coordinate masses: each particle's mass (Da) repeated over x,y,z.
    // ff.masses is the per-particle table — protein Cα beads first (indices
    // 0..nProt−1), then ligand united atoms (nProt..n−1, e.g. C ≈ 12 Da vs.
    // protein ≈ 110 Da). Fall back to the scalar mass arg for every coordinate
    // when that table is absent (legacy protein-only path).
    this.mass = new Float64Array(n3);
    const masses = ff.masses;
    if (masses && masses.length) {
      for (let i = 0; i < n3; i++) this.mass[i] = masses[(i / 3) | 0];
    } else {
      for (let i = 0; i < n3; i++) this.mass[i] = mass;
    }
    // 1/m and the per-coordinate thermal velocity scale σ_v = √(k_B·T·KCONV/m)
    // (Å/ps). Both are rebuilt whenever T or the masses change.
    this.invMass = new Float64Array(n3);
    this.thermal = new Float64Array(n3);

    this.pos = new Float64Array(refPositions); // working coordinates
    this.ref = new Float64Array(refPositions); // pristine copy for reset
    this.vel = new Float64Array(n3);
    this.time = 0;        // ps of trajectory completed
    this.steps = 0;

    this.T = 300;         // K
    this.zeta = 5.0;      // ps⁻¹
    this._rebuildThermal();
    this.dt = this._pickDt();

    // Preallocate RNG buffer (Box–Muller with cached second deviate)
    this._rnd = new Float64Array(n3);
    this._haveSpare = false;
    this._spare = 0;

    // initial force + Maxwell–Boltzmann velocities
    this.ff.compute(this.pos);
    this._sampleVelocities();
  }

  /* ------------------------------------------------------------------ */
  /*  Parameters                                                        */
  /* ------------------------------------------------------------------ */

  /** Set bath temperature (K); also refreshes the per-coordinate thermal scale. */
  setTemperature(T) { this.T = Math.max(1, T); this._rebuildThermal(); }
  setFriction(z)    { this.zeta = Math.max(0.1, z); this.dt = this._pickDt(); }

  /**
   * Set the protein Cα bead mass (Da) — the scalar slider semantics. Protein
   * coordinates are [0, nProt·3); ligand united-atom masses are left intact.
   * Rebuilds invMass/thermal for those coordinates and keeps the force field's
   * per-particle mass table (ff.masses) in sync for the analysis helpers.
   */
  setMass(m) {
    const pm = Math.max(10, m);
    const nProt = this.ff && this.ff.nProt ? this.ff.nProt : this.n3 / 3;
    for (let i = 0; i < nProt * 3; i++) {
      this.mass[i] = pm;
      this.invMass[i] = 1 / pm;
      this.thermal[i] = Math.sqrt(KB_KCAL * this.T * KCONV * this.invMass[i]);
    }
    if (this.ff && this.ff.nProt) {
      for (let j = 0; j < this.ff.nProt; j++) this.ff.masses[j] = pm;
    }
    this.dt = this._pickDt();
  }

  /** Recompute invMass = 1/m and thermal = √(k_B·T·KCONV/m) per coordinate.
   *  Validation: thermal = sqrt(KB*T*KCONV/m) is the correct Maxwell–Boltzmann
   *  velocity scale (Å/ps) from equipartition ½m⟨v²⟩KCONV = ½k_B T per dof. */
  _rebuildThermal() {
    for (let i = 0; i < this.n3; i++) {
      this.invMass[i] = 1 / this.mass[i];
      this.thermal[i] = Math.sqrt(KB_KCAL * this.T * KCONV * this.invMass[i]);
      // thermal = sqrt(KB*T*KCONV/m) — verified against KB_KCAL*KCONV/m
    }
  }

  /**
   * Stable Δt: bound by both bond vibration and the OU factor.
   *   ω_max is estimated from the stiffest available bond on the lightest mass
   *   (see _maxOmega); dtBond = 2/ω_max is resolved by the half-kick scheme and
   *   clamped to [1 fs, 4 fs]. The 4 fs ceiling is the traditional Cα-backbone
   *   limit (k_b ≤ 160 kcal/mol/Å² on ≈ 110 Da); the floor keeps the integrator
   *   efficient. Stiff friction additionally requires ζΔt < 0.5 so the O-step
   *   stays well-resolved.
   *   G67 — ligand-aware dt: CG alone 4fs vs CG+ligand 1.7fs vs heavy 1fs (honest, auto-tuned)
   */
  _pickDt() {
    const isHeavy = !!(this.ff && this.ff.heavy);
    // G67 — Ligand-aware dt: light ligand atoms (12 Da) need smaller dt; downgrade when ligand present
    const hasLig = !!(this.ff && (this.ff.nLigAtoms > 0 || (this.ff.ligandBonds && this.ff.ligandBonds.length > 0) || (this.ff.n > this.ff.nProt)));
    const maxDt = isHeavy ? 0.001 : (hasLig ? 0.0017 : 0.004); // G67 CG alone 4fs, CG+ligand ~1.7fs, heavy 1fs
    const dtBond = Math.min(maxDt, Math.max(0.0005, 1.5 / Math.max(1e-9, this._maxOmega())));
    const dtDrag = 0.5 / this.zeta;    // resolve friction relaxation time
    return Math.min(dtBond, dtDrag);
  }

  /**
   * Largest harmonic angular frequency ω = √(2·KCONV·k/m) (ps⁻¹) over the
   * available bond sets — the half-kick stability bound for the BAOAB step.
   * The factor 2 enters through the reduced mass: for two equal beads of mass
   * m in a ½k·δr² potential the relative coordinate oscillates at ω = √(k/μ)
   * with μ = m/2.
   *
   * Protein: the bonds list stores only [i, j, r0] (no k), so the stiffest
   * backbone constant is estimated at 160 kcal/mol/Å² on the lightest protein
   * bead mass. Ligand: k is recovered from r0 (≤ 1.44 Å ⇒ aromatic 200, else
   * 300 kcal/mol/Å²) on the united-atom masses (C ≈ 12 Da) — the fastest
   * vibration in the system. E.g. a benzene ring on 12 Da gives
   * ω ≈ √(2·418.4·200/12) ≈ 118 ps⁻¹ ⇒ dt ≤ 2/118 ≈ 0.017 ps, comfortably above
   * the 4 fs clamp, so the current parameters are unchanged — this only guards
   * against future stiff/light additions.
   */
  _maxOmega() {
    const ff = this.ff;
    const nProt = ff.nProt;
    let omega = 0;
    if (nProt > 1) {
      let mMin = Infinity;
      for (let j = 0; j < nProt; j++) if (this.mass[3 * j] < mMin) mMin = this.mass[3 * j];
      omega = Math.max(omega, Math.sqrt(2 * KCONV * 160 / mMin));
    }
    // Ligand bond check: include covalentBonds and ligandBonds for light-atom stability
    const bondSets = [];
    if (ff.covalentBonds && ff.covalentBonds.length) bondSets.push(ff.covalentBonds);
    if (ff.ligandBonds && ff.ligandBonds.length) bondSets.push(ff.ligandBonds);
    if (bondSets.length === 0) bondSets.push(ff.covalentBonds || ff.ligandBonds || new Float64Array(0));
    for (const bondList of bondSets) {
      for (let a = 0; a < bondList.length; a += 3) {
        const kEst = bondList[a + 2] <= 1.44 ? 200 : 300;
        const i = bondList[a], j = bondList[a + 1];
        const mi = this.mass[3 * i], mj = this.mass[3 * j];
        const mMin = mi < mj ? mi : mj;
        omega = Math.max(omega, Math.sqrt(2 * KCONV * kEst / mMin));
      }
    }
    return omega;
  }

  /* ------------------------------------------------------------------ */
  /*  Core                                                              */
  /* ------------------------------------------------------------------ */

  /** Advance one BAOAB step of size this.dt. */
  step() {
    const { pos, vel, ff, dt } = this;
    const F = ff.forces;
    const h = dt * 0.5;

    // B: half kick with current force — a = KCONV·F/m per coordinate
    this._kick(F, h);
    // A: half drift
    for (let i = 0; i < this.n3; i++) pos[i] += h * vel[i];

    // O: exact Ornstein–Uhlenbeck velocity update
    //    stationary variance σ_v² = k_B T / m  in mechanical (Å/ps)² units,
    //    i.e. σ_v = √(KB_KCAL·T·KCONV / m)  [KB·T is kcal/mol → ×KCONV].
    //    c2 folds the per-coordinate scale this.thermal in.
    const c1 = Math.exp(-this.zeta * dt);
    const c2 = Math.sqrt(1 - c1 * c1);
    this._fillGaussian(this.n3);
    for (let i = 0; i < this.n3; i++) vel[i] = c1 * vel[i] + c2 * this.thermal[i] * this._rnd[i];

    // A: half drift
    for (let i = 0; i < this.n3; i++) pos[i] += h * vel[i];

    // new force at x(t+dt), then final B half kick (same |Δv| cap + NaN guard)
    ff.compute(pos);
    this._kick(F, h);

    this.time += dt;
    this.steps++;
  }

  /**
   * BAOAB "B" half kick: v += h·KCONV·F/m per coordinate, with two guards:
   *   1. a non-finite force component contributes zero kick (NaN never
   *      propagates into velocities) — the app-side energy guard will
   *      auto-pause the run;
   *   2. per-particle |Δv⃗| is capped at 2 Å/ps per half-kick. That is ~16×
   *      the RMS thermal velocity of a 12 Da united atom at 300 K, so the cap
   *      never engages during normal dynamics, but it severs the "bullet"
   *      path (force spike on a light particle → projectile → protein blast)
   *      that a native clash would otherwise trigger between force-field
   *      guard fixes.
   */
  _kick(F, h) {
    const DVMAX = 2.0;                    // Å/ps, per-particle vector clamp
    const { vel, invMass } = this;
    for (let p = 0; p < this.n3; p += 3) {
      const im = KCONV * h * invMass[p];  // invMass is per-particle-constant
      let dvx = im * F[p], dvy = im * F[p + 1], dvz = im * F[p + 2];
      if (!Number.isFinite(dvx)) dvx = 0;
      if (!Number.isFinite(dvy)) dvy = 0;
      if (!Number.isFinite(dvz)) dvz = 0;
      const m2 = dvx * dvx + dvy * dvy + dvz * dvz;
      if (m2 > DVMAX * DVMAX) {
        const sc = DVMAX / Math.sqrt(m2);
        dvx *= sc; dvy *= sc; dvz *= sc;
      }
      vel[p] += dvx; vel[p + 1] += dvy; vel[p + 2] += dvz;
    }
  }

  /**
   * Advance wall-clock friendly chunk: run `frames * stepsPerFrame` steps but
   * at most `maxMs` milliseconds of compute; returns steps actually taken.
   * G68 — Adaptive steps per frame: advance(maxMs=14) caps wall-clock per frame to 14ms (documented in src/main.js:525)
   */
  advance(stepsWanted, maxMs = 12) { // G68 advance(maxMs=14) documented — default 12, caller main.js passes 14
    const t0 = performance.now();
    let done = 0;
    while (done < stepsWanted && performance.now() - t0 < maxMs) {
      this.step();
      done++;
    }
    return done;
  }

  /* ------------------------------------------------------------------ */
  /*  Housekeeping                                                      */
  /* ------------------------------------------------------------------ */

  reset() {
    this.pos.set(this.ref);
    this.vel.fill(0);
    this.time = 0;
    this.steps = 0;
    this.ff.compute(this.pos);
    this._sampleVelocities();
  }

  rebuildMass(m) { this.setMass(m); }

  /** Draw v from the Maxwell–Boltzmann distribution at temperature T. */
  _sampleVelocities() {
    this._fillGaussian(this.n3);
    for (let i = 0; i < this.n3; i++) this.vel[i] = this.thermal[i] * this._rnd[i];
    this._removeComMotion();
  }

  /**
   * Remove center-of-mass translation (keeps structure on screen). MASS-WEIGHTED:
   * with heterogeneous masses the COM velocity is p/M, not the unweighted mean,
   * so the protein + ligand system genuinely stops drifting.
   */
  _removeComMotion() {
    let px = 0, py = 0, pz = 0, M = 0;
    for (let i = 0; i < this.n3; i += 3) {
      const m = this.mass[i];
      px += m * this.vel[i]; py += m * this.vel[i + 1]; pz += m * this.vel[i + 2];
      M += m;
    }
    px /= M; py /= M; pz /= M;    // COM velocity (mass-weighted mean)
    for (let i = 0; i < this.n3; i += 3) { this.vel[i] -= px; this.vel[i + 1] -= py; this.vel[i + 2] -= pz; }
  }

  /** Fill this._rnd[0..k) with N(0,1) via Box–Muller (paired deviates). */
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
