/**
 * cg/forcefield.js — the CG ForceField class. Split out of src/forcefield.js;
 * every method body that did real work moved verbatim into the sibling module
 * named in its wrapper's comment. What remains here is the class: a
 * constructor that sequences the build steps, and the thin wrappers that keep
 * the public/protected surface byte-for-byte what it was.
 *
 * WHY A CLASS MOVE AND NOT PROTOTYPE COMPOSITION
 * ----------------------------------------------
 * src/heavy.js used prototype composition with Object.defineProperty and
 * hand-copied class-method descriptors, because it had to keep a 1674-line
 * file's exports stable while relocating bodies out of an ALREADY-COMPOSED
 * class. Here the situation is the opposite: there is exactly one class, it
 * moves as a unit, and nothing composes it from pieces. A plain `class`
 * declaration keeps `Object.getOwnPropertyDescriptor(ForceField.prototype,
 * 'compute')` a real class-method descriptor (non-enumerable, writable,
 * configurable) and makes `for..in` over an instance enumerate exactly the
 * own instance fields and nothing from the prototype — identical to before,
 * by construction rather than by care. Composition here would be strictly
 * more machinery for the same result.
 *
 * CONSTRUCTION ORDER (a contract, not an accident):
 *   initScalars        — rc/gamma/k/stiffness/flags; needs nothing
 *   buildParticles     — ref/n/masses/mask; needs nLigAtoms = 0 from scalars
 *   assignParticleParams — per-bead + per-atom tables; needs chargesOn, hbMode,
 *                         ligandAtoms and nLigAtoms from the two steps above
 *   buildTopology      — bonds/angles/ENM/holo/exclusions/contacts; needs ref,
 *                        rc, gamma, _dist/_pairKey, and the ligSigma table
 *   allocBuffers       — LAST: needs final n and nLigAtoms
 * Every step below does the same work in the same order as the monolithic
 * constructor it replaces, which is why the seeded golden is bit-identical.
 */

import { initScalars, assignParticleParams } from "./params.js";
import { buildParticles } from "./system.js";
import { buildTopology } from "./topology.js";
import { allocBuffers, encodeCell, cellKey, decodeX, decodeY, decodeZ } from "./grid.js";
import { computeForces } from "./compute.js";
import { describePhysics } from "./level.js";
import {
  rebuildHoloSprings, setSpringScale, clearSpringScale, applySeqWeights, useTirionNetwork,
} from "./springs.js";
import {
  harmonicPairs, springForces, angleForces, ligandBondForces, improperForces,
} from "../ff-harmonic.js";
import { repulsion } from "../ff-repulsion.js";
import { binding } from "../ff-binding.js";
import { kineticTemp, rmsdTo } from "../physics/observables.js";

export class ForceField {
  /**
   * @param {object} sel  output of pdb.selectSystem(): {beads, segments}
   * @param {object} par  { rc, gamma }  — ENM cutoff (Å) and spring constant
   * @param {Array}  ligands  output of pdb.parseLigands() (empty ⇒ protein-only)
   */
  constructor(sel, par = {}, ligands = []) {
    const { beads, segments } = sel;
    const nProt = this.nProt = beads.length;

    // 1. scalars + tier resolution (src/cg/params.js)
    initScalars(this, par);
    // 2. reference coords, appended ligand particles, masses, rmsd mask
    //    (src/cg/system.js)
    buildParticles(this, beads, par, ligands);
    // 3. per-residue-class and per-element type tables (src/cg/params.js)
    assignParticleParams(this, beads);
    // 4. connectivity, ENM contact set, exclusion ledger, native contacts
    //    (src/cg/topology.js — includes the holo + seqWeight + Tirion opt-ins)
    buildTopology(this, par, segments);
    // 5. preallocated force/grid/desolvation buffers (src/cg/grid.js)
    allocBuffers(this, par);
  }

  _pairKey(i, j) { return i < j ? i * 1e6 + j : j * 1e6 + i; }

  _dist(p, i, j) {
    const dx = p[3 * j] - p[3 * i], dy = p[3 * j + 1] - p[3 * i + 1], dz = p[3 * j + 2] - p[3 * i + 2];
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  _angle(p, i, j, k) {
    const ax = p[3 * i] - p[3 * j], ay = p[3 * i + 1] - p[3 * j + 1], az = p[3 * i + 2] - p[3 * j + 2];
    const bx = p[3 * k] - p[3 * j], by = p[3 * k + 1] - p[3 * j + 1], bz = p[3 * k + 2] - p[3 * j + 2];
    const la = Math.hypot(ax, ay, az), lb = Math.hypot(bx, by, bz);
    let c = (ax * bx + ay * by + az * bz) / (la * lb);
    c = Math.min(1, Math.max(-1, c));
    return Math.acos(c);
  }

  /**
   * compute(pos) → fills this.forces and returns total potential energy.
   * pos: Float64Array(3n). Forces are −∇U in kcal/mol/Å.
   *
   * The term order, the inlined holo compression floor and the NaN guard all
   * live in src/cg/compute.js — see that file for the zero-alloc audit and
   * for why the summation order is a physics contract.
   */
  compute(pos) { return computeForces(this, pos); }

  /** Attach an optional Funnel instance (constructed by main.js, not imported here). */
  setFunnel(fn) { this.funnel = fn; }

  /**
   * Explicit physics-fidelity descriptor (read-only, headless-safe). Pure
   * introspection for UI/spec wiring and regression tests; no energy effect.
   * Implementation: src/cg/level.js:describePhysics.
   * @returns {{level:string, charges:boolean, hbMode:string, enm:string,
   *   coulombActive:boolean, directionalHBActive:boolean,
   *   isSimplifiedDefault:boolean}}
   */
  describePhysics() { return describePhysics(this); }

  /**
   * Rebuild the holo (native-pose) protein–ligand springs from the CURRENT
   * reference coordinates, keeping the exclusion set in sync. Implementation:
   * src/cg/springs.js:rebuildHoloSprings.
   */
  rebuildHoloSprings() { rebuildHoloSprings(this); }

  /** Σ 1⁄2 k (r−r0)2 over a flat pair list; kernel lives in ff-harmonic.js. */
  _harmonicPairs(pos, f, list, stride, k) { return harmonicPairs(pos, f, list, stride, k); }

  /** ENM spring forces with per-spring k — kernel in ff-harmonic.js. */
  _springForces(pos, f) { return springForces(this, pos, f); }

  /**
   * Apply an ML contact prior to the ENM spring constants.
   * @param {Array} contacts  [[i, j, p], ...] residue pairs (0-indexed) with p ∈ [0,1].
   * @param {number} alpha    global scale (k_per_pair = gamma · (1 + alpha·p)).
   * Marks springScaleActive so the UI knows a custom map is live.
   * Implementation: src/cg/springs.js:setSpringScale.
   */
  setSpringScale(contacts, alpha = 1) { setSpringScale(this, contacts, alpha); }

  /** Restore uniform ENM stiffness (undo an ML contact-prior scaling). */
  clearSpringScale() { clearSpringScale(this); }

  /**
   * Sequence-dependent ENM (Bahar-style) — opt-in stub.
   * Rebuilds springK as K_seq(i,j) = gamma * (1 + 0.2*(w_i + w_j)/2). Preserves
   * the uniform topology — only the stiffness is chemistry-weighted. See
   * src/ff-params.js:SEQ_WEIGHT, docs/CG_HEAVY.md, tests/test_enm_seq.js.
   * Implementation: src/cg/springs.js:applySeqWeights.
   * @param {Array} beads  optional override (defaults to constructor beads if stored)
   */
  applySeqWeights(beads = null) { applySeqWeights(this, beads); }

  /**
   * Opt-in Tirion distance-weighted ENM + SS dihedral basins (Phase 1).
   * Rebuilds springs/springK as γ_ij = γ0·(R0/r0_ij)^6 over the same
   * H(Rc−r0) topology and stores backbone pseudo-dihedral basins on
   * ff.tirionDihedrals / ff.tirionDihedralK / ff.tirionSS. Falls back to the
   * existing uniform network when the module is unavailable.
   * Implementation: src/cg/springs.js:useTirionNetwork.
   * @param {object} [opts]  { gamma0, R0, cutoff, segments }
   */
  useTirionNetwork(opts = {}) { useTirionNetwork(this, opts); }

  /** U_θ = ½ kθ (θ−θ0)² — angle-bending kernel lives in ff-harmonic.js. */
  _angleForces(pos, f, list = this.angles, k = this.kAngle) { return angleForces(this, pos, f, list, k); }

  /**
   * Ligand united-atom internal energy (bonds, angles, improper planarity).
   * No-op (returns 0) when no ligands are present.
   */
  _ligandInternal(pos, f) {
    if (this.ligandAtoms.length === 0) return 0;
    let U = 0;
    U += this._ligandBondForces(pos, f);
    U += this._angleForces(pos, f, this.ligandAngles, 40);
    U += this._improperForces(pos, f, this.ligandImpropers);
    return U;
  }

  /** Ligand bonds U = Σ ½ k (r−r0)² — kernel in ff-harmonic.js. */
  _ligandBondForces(pos, f) { return ligandBondForces(this, pos, f); }

  /** Improper (out-of-plane) term — FD kernel lives in ff-harmonic.js. */
  _improperForces(pos, f, list) { return improperForces(pos, f, list); }

  /**
   * Protein–ligand binding potential — pair pass over every protein bead
   * (0..nProt−1) vs every ligand atom (global la = nProt + a, table a).
   * Three short-range cross terms, each cut off smoothly at bindRcut:
   *
   *   U_LJ = 4ε[(σ/r)¹² − (σ/r)⁶]·sw            cross 12-6 (attractive well)
   *   U_el = COULOMB_CONST·q1·q2/(εr·r)·sw     screened electrostatics,
   *                                              εr(r) = 4 + 76·tanh(r/8)
   *   U_HB = −epsHB·exp(−(r−r0)²/2w²)·sw          H-bond well, r0=3.2, w=0.6
   *
   * with Lorentz–Berthelot combining rules σ = (σi+σa)/2, ε = √(εi·εa).
   *
   * SWITCHING: the smoothstep is identically 1 for r ≤ rsw and falls to 0 at
   * rc with a vanishing derivative, so both energy and force are continuous:
   *     t   = (r − rsw)/(rc − rsw),   Δ = rc − rsw
   *     s   = 1 − t³(10 − 15t + 6t²)          [cubic smoothstep, s(0)=1, s(1)=0]
   *     s′  = ds/dr = −(30/Δ)·t²(1−t)²
   * (ds/dt = −30t²(1−t)² by direct differentiation; dt/dr = 1/Δ by chain rule.)
   *
   * SIGN CONVENTION: the force on protein bead i is −∂U/∂pos_i. With
   * dx = pos_i − pos_la pointing from the ligand atom to the bead,
   * ∂r/∂pos_i = dx/r, so f_i += −(dU/dr)·(dx/r) and f_la −= that (Newton's
   * 3rd law). Each term accumulates its dU/dr into one scalar and a single
   * vector update is applied — no sign ambiguity.
   *
   * EEF1-LITE BURIAL (Pass 2): the same grid scan additionally accumulates a
   * soft protein-occupancy count for each ligand atom a,
   *     n_a = Σ_j g(r_aj),   g(r) = exp(−(r−r0)²/(2σ²)),   r0=4.5 Å, σ=1.8 Å
   * (g has decayed to ≈ 0.04 at r = 9 Å = bindRcut, so no second cutoff is
   * needed). The burial fraction
   *     B_a = 1 − exp(−n_a/3) ∈ [0,1]     (0 exposed, 1 fully buried)
   * turns the per-atom hydrophobic transfer free energy ΔG_a into the burial
   * cost U_desolv = Σ_a ΔG_a·B_a (ΔG_a < 0 ⇒ burial lowers the binding energy).
   * The recorded (a, j, r) pairs from Pass 1 feed the force without re-scanning:
   *     dU_desolv/dr = ΔG_a·(dB/dn)·(dg/dr)
   *     dB/dn = exp(−n_a/3)/3 ,   dg/dr = −g(r)·(r−r0)/σ²
   * applied with the same sign convention as the pair terms above. Passing
   * through a burial maximum at r0 gives the correct physics: a fully exposed
   * atom (no stored pairs) feels no desolvation force, an approaching atom is
   * pulled into the protein (favorable ΔG·B), and one pushed inside below r0
   * is pushed back out. The ± sign is verified numerically against total-U
   * finite differences in the test suite.
   *
   * Two-pass structure:
   *   Pass 1 — pair terms (LJ / electrostatics / H-bond) computed exactly as
   *            before, plus dens[a] accumulation and (a, j, r) recording.
   *   Pass 2 — burial fraction/energy per ligand atom, then desolvation forces
   *            over the recorded pairs.
   */
  _binding(pos, f) { return binding(this, pos, f); }

  /** Repulsive-only 12-6 LJ excluded-volume pass — kernel in ff-repulsion.js. */
  _repulsion(pos, f) { return repulsion(this, pos, f); }

  /* --- numeric spatial-hash cell keys (no string allocs) -------------- */
  /** Codec + non-finite sentinel live in src/cg/grid.js. */
  _encodeCell(cx, cy, cz) { return encodeCell(cx, cy, cz); }
  _cellKey(x, y, z, cell) { return cellKey(x, y, z, cell); }
  _decodeX(k) { return decodeX(k); }
  _decodeY(k) { return decodeY(k); }
  _decodeZ(k) { return decodeZ(k); }

  /* ------------------------------------------------------------------ */
  /*  Analysis helpers (HUD)                                             */
  /* ------------------------------------------------------------------ */

  /**
   * Instantaneous kinetic temperature, T = ke_kcal/(1.5·n·k_B) — units, formula
   * and the mass-layout contract are documented once, in physics/observables.js.
   * CG is always the flat per-DOF table (3n) LangevinIntegrator builds, so the
   * layout is named, not inferred. D1: a per-atom table still reads past its
   * end and yields NaN here — the historical CG behaviour, kept verbatim.
   */
  kineticTemp(vel, mass) { return kineticTemp(vel, mass, this.n, "dof"); }

  /**
   * RMSD to native over the protein Cα beads only (ligand coords are rigid
   * internal DOF and excluded from the fold-space metric). No alignment; ENM
   * keeps the COM/orientation nearly fixed. See rmsdTo() for the formula.
   */
  rmsd(pos) { return rmsdTo(pos, this.ref, this._maskProt, this.nProt); }
}