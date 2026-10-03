/**
 * heavy/forcefield.js — the HeavyForceField class: construction, state, and the
 * four kernel groups attached to its prototype.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10). This file owns the class itself:
 * the constructor (index blocks, masses, reference coords, charges, element
 * parameters, exclusion/1-4 sets, ligand bookkeeping, the opt-in adapter
 * defaults), _buildAmberStiffness, refreshOBC2Radii, the spring-scale setters,
 * and the prototype wiring below.
 *
 * THE FOUR ATTACHED GROUPS, AND WHY THEY ARE ATTACHED RATHER THAN INHERITED
 *   energyMethods       compute()                      -> heavy/energy.js
 *   nonBondedKernels    the three grid kernels         -> heavy/nonbonded.js
 *   weakKernels         weak lists + weak pass         -> heavy/weak.js
 *   observableMethods   kineticTemp / rmsd*            -> heavy/observables.js
 *   physicsLevelMethods describePhysics()              -> heavy/level.js
 *
 * Composition over inheritance is deliberate: HeavyForceField has exactly one
 * base class in the whole repo (none), a subclass would invert the dependency
 * (a kernel module would have to import the class), and a subclass is a new
 * public type that a caller could observe. attachMethods() below therefore
 * defines each method with the EXACT property descriptor a `class` method has
 * — non-enumerable, writable, configurable — so method lookup order, `this`
 * binding, and any reflection over the instance or its prototype are
 * byte-identical to before the split. It is not `Object.assign`, whose
 * enumerable own properties would be visible to for..in and
 * Object.getOwnPropertyNames(...).enumerable.
 */
import { resolveHeavyElementParams } from "../ff-params.js";
import { assignCharges, GB_RADII } from "../physics/charges.js";
import { GeneralizedBorn } from "../physics/gb.js";
import { SasaModel } from "../physics/sasa.js";
import { DirectionalHBond } from "../physics/hbond.js";
import { SpatialGrid } from "../spatial-grid.js";
import { getBondParams, getAngleParams, getNonbondedParams } from "../physics/forcefield/amber14sb.js";
import { computeBornRadii as computeOBC2Radii } from "../physics/solvation/gb_obc2.js";
import { buildTopology, buildMetalCoordination, buildLists } from "./topology.js";
import { R_CUT, METAL_K, heavyMass } from "./params.js";
import { pairKey } from "./pairs.js";
import { resolveHeavyPhysicsLevel, physicsLevelMethods } from "./level.js";
import { energyMethods } from "./energy.js";
import { nonBondedKernels } from "./nonbonded.js";
import { weakKernels } from "./weak.js";
import { observableMethods } from "./observables.js";

/**
 * Define each named method on a prototype with a class-method descriptor
 * (non-enumerable, writable, configurable). Object.assign is deliberately not
 * used: it would make every attached method an enumerable own property.
 * @param {object} proto target prototype
 * @param {object} methods object literal of methods
 */
function attachMethods(proto, methods) {
  for (const name of Object.keys(methods)) {
    Object.defineProperty(proto, name, {
      value: methods[name],
      writable: true,
      configurable: true,
      enumerable: false,
    });
  }
}

/**
 * HeavyForceField — High performance all-atom heavy force field.
 *
 * FP5: optional 4th arg `opts.topo` accepts a prebuilt topology
 * ({bonds, angles, propers, impropers} as in buildTopology, e.g. from the
 * chunked buildTopologyChunked) and skips the internal buildTopology call.
 * Absent → internal build, bit-identical to before (all existing 3-arg
 * callers untouched).
 *
 * Revolution 3 / Issue 1: queryable physics-level mirror (no kernel change).
 * Pre-fix the constructor ignored par.physicsLevel / par.binding.charges /
 * par.binding.hbMode, so UI L1/L2 silently no-opped in heavy except the weak
 * flag. The heavy kernels already run full physics (AMBER charges + GB +
 * directional H-bonds always on), so these flags are stored queryably only
 * (this.physicsLevel / chargesOn / hbMode + describePhysics()) and never
 * gate energy/force math — U/forces are bit-identical across tiers.
 */

export class HeavyForceField {
  constructor(system, par = {}, ligands = [], opts = {}) {
    const atoms = system.atoms;
    this.n = atoms.length;
    this.nProt = atoms.filter((a) => a.isProtein).length;
    this.nLigAtoms = atoms.filter((a) => a.isLigand).length;
    this.nHetero = this.n - this.nProt - this.nLigAtoms;
    this.ligandStart = this.nProt + this.nHetero;
    this.heteroAtoms = atoms.slice(this.nProt, this.ligandStart);
    this.gamma = par.gamma ?? 1.0;
    this.heavy = true;
    // Revolution 3 / Issue 1: store the requested tier + mirror the CG
    // binding flags queryably. Explicit par.binding.charges/hbMode win over
    // the tier (same rule as ForceField/resolvePhysicsLevel); unknown levels
    // fall back to L0. Queryable only — no kernel math reads these flags.
    {
      const _phys = resolveHeavyPhysicsLevel(par);
      this.physicsLevel = _phys.level;
      this.chargesOn = _phys.charges;
      this.hbMode = _phys.hbMode;
    }

    // Fast Spatial Grid for O(N) neighbor searches
    // G66 — Verlet skin 2Å, rebuild every 10 steps, 20% cut — aspirational target; currently rebuilds every step via SpatialGrid.build() with R_CUT=8.5Å (skin not yet implemented)
    this.grid = new SpatialGrid(R_CUT, this.n);

    // Biophysics modules
    this.gb = new GeneralizedBorn({ epsIn: 4.0, epsOut: 78.5, saltM: 0.15, temperature: par.temp ?? 300 });
    this.sasa = new SasaModel({ gamma: 0.0072 });
    this.hbond = new DirectionalHBond({ epsHB: 2.5 });
    this._hbClassification = this.hbond.classifyAtoms(atoms);

    // Residue -> CA map
    this._resCa = [];
    const caIdx = new Map();
    for (let i = 0; i < this.n; i++) {
      const a = atoms[i];
      if (a.isProtein && a.atomName === "CA" && !caIdx.has(a.chain + "|" + a.resSeq)) {
        caIdx.set(a.chain + "|" + a.resSeq, i);
      }
    }
    const seenRes = new Set();
    for (let i = 0; i < this.n; i++) {
      const a = atoms[i];
      if (!a.isProtein) continue;
      const key = a.chain + "|" + a.resSeq;
      if (seenRes.has(key)) continue;
      seenRes.add(key);
      this._resCa.push(caIdx.get(key));
    }

    // Assign physically grounded partial charges & LJ parameters
    const charges = assignCharges(atoms);
    this._elem = new Array(this.n);
    this._charges = charges;
    this._bornRadii = this.gb.computeBornRadii(atoms);

    for (let i = 0; i < this.n; i++) {
      const el = atoms[i].element;
      // M7: the element rule is `METAL_ELEMENT[el] ?? resolveElementParams(el)`,
      // now a named function in the canonical parameter module instead of a
      // ternary here. Same precedence, same values, byte for byte — only the
      // number of places a per-element parameter can be typed is reduced to 1.
      this._elem[i] = { ...resolveHeavyElementParams(el), q: charges[i] };
    }

    // Masses (Da)
    this.masses = new Float64Array(this.n);
    for (let i = 0; i < this.n; i++) this.masses[i] = heavyMass(atoms[i].element);

    this.ref = new Float64Array(this.n * 3);
    for (let i = 0; i < this.n; i++) {
      this.ref[3 * i] = atoms[i].x;
      this.ref[3 * i + 1] = atoms[i].y;
      this.ref[3 * i + 2] = atoms[i].z;
    }
    this.forces = new Float64Array(this.n * 3);
    this.energy = 0;

    // Fold-metric masks for rmsd()/rmsdLig(), built once so the HUD's per-frame
    // calls allocate nothing: _maskProt → protein [0, nProt), _maskLig → ligand
    // block [ligandStart, n), the slice rmsdLig has always summed (ligandStart +
    // nLigAtoms === n by construction). Formula: physics/observables.js.
    this._maskProt = new Uint8Array(this.n);
    this._maskProt.fill(1, 0, this.nProt);
    this._maskLig = new Uint8Array(this.n);
    this._maskLig.fill(1, this.ligandStart, this.n);

    // Topology (FP5: prebuilt chunked topo skips the sync O(n²) build).
    const topo = opts?.topo ?? buildTopology(atoms);
    const coordPairs = buildMetalCoordination(atoms);
    const L = buildLists(atoms, topo, coordPairs);
    this.bonds = new Float64Array(L.bonds);
    this.angles = new Float64Array(L.angles);
    this.impropers = new Float64Array(L.impropers);
    this.propers = new Float64Array(L.propers);
    this.coord = new Float64Array(L.coord);
    this.holoSprings = L.holoSprings;
    this.nHolo = 0;
    this.holoOn = false;

    // Non-bonded exclusions
    this._excluded = new Set();
    this._scale14 = new Map();
    for (let a = 0; a < this.bonds.length; a += 3) {
      this._excluded.add(pairKey(this.bonds[a], this.bonds[a + 1]));
    }
    for (let a = 0; a < this.angles.length; a += 4) {
      this._excluded.add(pairKey(this.angles[a], this.angles[a + 2]));
    }
    for (let a = 0; a < this.propers.length; a += 5) {
      const k = pairKey(this.propers[a], this.propers[a + 3]);
      this._scale14.set(k, 0.5);
    }
    // Exclude intra-ligand non-bonded interactions so small molecules preserve their true geometry
    for (let i = this.ligandStart; i < this.n; i++) {
      for (let j = i + 1; j < this.n; j++) {
        this._excluded.add(pairKey(i, j));
      }
    }

    this.covalentBonds = new Float64Array(L.bonds);
    this.ligandAtoms = atoms.slice(this.ligandStart);
    const ligBonds = [];
    for (let a = 0; a < L.bonds.length; a += 3) {
      if (L.bonds[a] >= this.ligandStart && L.bonds[a + 1] >= this.ligandStart) {
        ligBonds.push(L.bonds[a], L.bonds[a + 1], L.bonds[a + 2]);
      }
    }
    this.ligandBonds = new Float64Array(ligBonds);
    this.springs = new Float64Array(0);
    this.springK = new Float64Array(0);
    this.springScaleActive = false;
    this.springU = 0;
    this.nativeContacts = new Float64Array(0);

    this.kBond = 200.0;
    this.kAngle = 40.0;
    this.kImproper = 20.0;
    this.kProper = 2.0;
    this.metalK = METAL_K;

    this.funnel = null;
    this.funnelOn = false;

    // Per-atom physical mass table (AMBER ff14SB masses in Da)
    this.masses = new Float64Array(this.n);
    for (let i = 0; i < this.n; i++) {
      this.masses[i] = heavyMass(atoms[i].element);
    }

    this.bindingU = 0;
    this.desolvU = 0;
    // Loop-2 S4 (R4 §5 item 1): per-term binding-accumulator opt-in.
    // trackTerms = true fills bindLJU/bindCoulU/bindHBU/desolvU + the
    // bindU component vector {lj, coul, hb, desolv, pi, cpi, xb} on each
    // compute(). DEFAULT false — legacy scalar path, zero overhead.
    this.trackTerms = false;
    this.bindLJU = 0; this.bindCoulU = 0; this.bindHBU = 0;
    this.bindU = { lj: 0, coul: 0, hb: 0, desolv: 0, pi: 0, cpi: 0, xb: 0 };
    this.repU = 0;
    this.bondU = 0;
    this.angleU = 0;
    this.improperU = 0;
    this.properU = 0;
    this.coordU = 0;
    this.elecU = 0;
    this.gbU = 0;
    this.sasaU = 0;
    this.hbondU = 0;
    this.membraneU = 0;
    this._nanStrikes = 0;

    // ---- Phase 1 opt-in adapters (defaults preserve legacy behavior) ----
    // useAmber14: per-bond/angle stiffness from amber14sb.js (else uniform
    //   kBond/kAngle). gbModel "obc2": GB-OBC2 radii + forces instead of HCT.
    //   sasaModel "lcpo": LCPO SASA instead of SasaModel. membrane: {on,
    //   thickness, width, epsWater, epsMem, zCenter} adds slab term.
    // All fall back to legacy tables when the new modules are unavailable.
    this.atoms = atoms;
    this.useAmber14 = par.useAmber14 ?? false;
    this.gbModel = par.gbModel ?? "hct";
    this.sasaModel = par.sasaModel ?? "sasa";
    this.membraneOpts = par.membrane ?? null;
    this.gbEpsIn = par.epsIn ?? 4.0;
    this.gbEpsOut = par.epsOut ?? 78.5;
    this.gbSaltM = par.saltM ?? 0.15;
    this._bondK = null;
    this._angleK = null;
    this._obc2Radii = null;
    this._lcpoElements = atoms.map((a) => a.element ?? "C");
    if (this.useAmber14) {
      try { this._buildAmberStiffness(); } catch (e) {
        console.warn(`[HeavyForceField] AMBER14 tables unavailable (${e.message}) — uniform k fallback`);
        this.useAmber14 = false;
      }
    }

    // ---- Loop-2 S3 opt-in weak interactions (R3 §6 items 1–3 + §1e) ----
    // par.weak: "off" (default, bit-identical legacy) | "on" adds π-stack,
    //   cation-π and halogen σ-hole terms after the LJ/GB grid pass.
    // par.metalAngles: true (default false) swaps metal distance springs for
    //   chem/metals.js enforceCoordination (radial k=40 + cross-angle k=20)
    //   for metals with a detected coordination geometry — R3 §1e.
    this.weakOn = par.weak === "on";
    this.metalAngles = par.metalAngles === true;
    this.weakU = 0; this.piU = 0; this.cpiU = 0; this.xbU = 0;
    this.coordAngleU = 0;
    this._weakRings = [];
    this._weakCations = [];
    this._weakHalogens = [];
    this._weakRingAtoms = new Set();
    this._metalEnforce = null; // { metals:[{index,element,donors}], elements, hasGeometry:Set }
    if (this.weakOn || this.metalAngles) this._buildWeakAndMetalLists(atoms, L.bonds);
  }

  /**
   * Build per-bond / per-angle stiffness from AMBER ff14SB tables.
   * Falls back to uniform kBond/kAngle entries when a lookup misses.
   * Non-breaking: only consumed by compute() when useAmber14 is true.
   */
  _buildAmberStiffness() {
    const nb = this.bonds.length / 3;
    this._bondK = new Float64Array(nb);
    for (let b = 0; b < nb; b++) {
      const i = this.bonds[3 * b], j = this.bonds[3 * b + 1];
      const A = this.atoms[i], B = this.atoms[j];
      const key1 = A?.atomName ?? A?.element ?? "C";
      const key2 = B?.atomName ?? B?.element ?? "C";
      let p = null;
      try { p = getBondParams(key1, key2); } catch { p = null; }
      if (!p) { try { p = getNonbondedParams(A?.element ?? "C"); } catch { p = null; } }
      this._bondK[b] = p?.k ?? this.kBond;
    }
    const na = this.angles.length / 4;
    this._angleK = new Float64Array(na);
    for (let a = 0; a < na; a++) {
      const i = this.angles[4 * a], j = this.angles[4 * a + 1], k = this.angles[4 * a + 2];
      const A = this.atoms[i], B = this.atoms[j], C = this.atoms[k];
      let p = null;
      try {
        p = getAngleParams(A?.atomName ?? A?.element ?? "C", B?.atomName ?? B?.element ?? "C", C?.atomName ?? C?.element ?? "C");
      } catch { p = null; }
      this._angleK[a] = p?.k ?? this.kAngle;
    }
  }

  /**
   * Refresh OBC-II Born radii for the current positions (Phase 1 helper).
   * @param {ArrayLike<number>} pos  flat 3n
   * @returns {Float64Array} effective radii
   */
  refreshOBC2Radii(pos) {
    const intrinsic = new Float64Array(this.n);
    for (let i = 0; i < this.n; i++) {
      const el = this.atoms[i]?.element ?? "C";
      intrinsic[i] = GB_RADII[el] ?? GB_RADII.DEFAULT ?? 1.6;
    }
    this._obc2Radii = computeOBC2Radii(pos, intrinsic, {});
    return this._obc2Radii;
  }

  setFunnel(fn) { this.funnel = fn; }

  setSpringScale(contacts, alpha = 1) {
    const S = [], K = [];
    for (const [i, j, p] of contacts) {
      if (p <= 0) continue;
      const a = this._resCa[i], b = this._resCa[j];
      if (a === undefined || b === undefined || a === b) continue;
      const i3 = 3 * a, j3 = 3 * b;
      const r0 = Math.hypot(
        this.ref[j3] - this.ref[i3],
        this.ref[j3 + 1] - this.ref[i3 + 1],
        this.ref[j3 + 2] - this.ref[i3 + 2],
      );
      S.push(a, b, r0);
      K.push(this.gamma * (1 + alpha * p));
    }
    this.springs = new Float64Array(S);
    this.springK = new Float64Array(K);
    this.springScaleActive = true;
  }

  clearSpringScale() {
    this.springs = new Float64Array(0);
    this.springK = new Float64Array(0);
    this.springScaleActive = false;
  }

  rebuildHoloSprings() {}

}

// ── Kernel attachment ───────────────────────────────────────────────────
// compute(), the grid kernels, the weak pass, the fold metrics and the physics
// descriptor live in their own modules (see the header map above). Attaching
// them here keeps HeavyForceField a single object with no subclass and no
// behavioural difference from the pre-split class.
attachMethods(HeavyForceField.prototype, physicsLevelMethods);
attachMethods(HeavyForceField.prototype, energyMethods);
attachMethods(HeavyForceField.prototype, nonBondedKernels);
attachMethods(HeavyForceField.prototype, weakKernels);
attachMethods(HeavyForceField.prototype, observableMethods);
