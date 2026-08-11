/**
 * heavy.js — all-atom "heavy mode" for the simulator (item 3).
 *
 * The default coarse-grained model represents a protein as one bead per
 * residue at its Cα (elastic-network / ENM springs) and a ligand as explicit
 * united atoms. Heavy mode is the complementary, higher-resolution view: every
 * heavy (non-hydrogen) atom of the protein side-chains *and* backbone is kept
 * explicit, with a full covalent topology (bonds, angles, impropers, proper
 * dihedrals) rebuilt from the crystallographic coordinates, plus a simple
 * AMBER-style non-bonded term (Lennard-Jones + screened electrostatics) and
 * metal-ion coordination springs.
 *
 * Nothing here alters the coarse-grained force field — heavy mode is a
 * *separate* ForceField-compatible class (HeavyForceField) that main.js can
 * build instead of (or alongside) the CG ForceField. It exposes the same
 * public surface the integrator / funnel / viewer / UI expect:
 *     n, nProt, nLigAtoms, masses, ref, forces, energy,
 *     compute(pos), rmsd(pos), kineticTemp(vel, mass),
 *     ligandAtoms, ligandBonds, holoSprings, nHolo, holoOn, springs,
 *     springK, setSpringScale, clearSpringScale, setFunnel, funnelOn, bonds, angles.
 *
 * Design notes
 * ------------
 *  - Solvent / water / hydrogen are dropped (united-atom convention, matches
 *    the CG model). Only heavy atoms survive.
 *  - Covalent bonds are detected geometrically: two heavy atoms within BOND_SLACK
 *    × (sum of covalent radii) are bonded. Angles/impropers/propers are then
 *    derived from the bond graph (1–2, 1–3, 1–4 atom triples/quartets).
 *  - Metal ions (ZN, FE, MG, CA, CU, MN, NI, CO, NA, K) do NOT form covalent
 *    bonds. Instead heavy.js detects donor atoms (N/O/S) within the metal's
 *    coordR and adds harmonic coordination springs (1–2 style) of moderate k,
 *    capped at coordN donors. This keeps the geometry stable without forcing
 *    a specific coordination chemistry.
 *  - Non-bonded pairs are excluded for 1–2 (bonded) and 1–3 (angle) pairs;
 *    1–4 (dihedral) pairs use a scale factor on the LJ/electrostatic terms.
 *  - The heavy mode is intentionally parameter-light and geometry-driven: it
 *    is a visual/structural exploration layer, not a production MD force field.
 *
 * The module is pure (imports only ff-params.js + ligand.js helpers) and has
 * no DOM access, so it is unit-testable under Node.
 */
import {
  METAL_ELEMENT, METAL_ELEMENT_DEFAULT,
  COVALENT_RADIUS, BOND_SLACK,
  LIG_ELEMENT, LIG_ELEMENT_DEFAULT,
  KB_KCAL, KCONV,
} from "./ff-params.js?v=10";
import { improperAngle } from "./ligand.js?v=10";

// Electrostatics: screened Coulomb prefactor (kcal/mol/Å per e²) and the
// Debye screening length (Å). The product q_i·q_j·K_ELEC / (ε(r)·r) with a
// distance-dependent dielectric ε(r) = r (a common protein-model choice that
// both softens short-range charge contacts and removes the 1/r singularity)
// and a moderate screening length keeps metal (2+) contacts from dominating
// while still giving a real charge-charge signal.
export const K_ELEC = 332.0;   // kcal·Å/(mol·e2)
export const SCREEN_LEN = 6.0; // Å (Debye screening length)

// Non-bonded switching cutoff (Å): interactions between R_SWITCH_ON and
// R_CUT are scaled by a smooth AMBER-style cubic switch that reaches 0 at
// R_CUT. This both keeps heavy mode stable (thousands of charged protein
// atoms would otherwise accumulate unphysical long-range energy) and bounds
// the O(n²) pair scan cost.
export const R_CUT = 8.0;
export const R_SWITCH_ON = 6.0;

/**
 * AMBER-style switching function: 1 below R_SWITCH_ON, smooth 1→0 between
 * R_SWITCH_ON and R_CUT, 0 beyond. Both S and dS/dr vanish at R_CUT, so the
 * energy and force are continuous at the cutoff.
 */
export function switchFunc(r) {
  if (r <= R_SWITCH_ON) return 1;
  if (r >= R_CUT) return 0;
  const rsq = r * r, rOn = R_SWITCH_ON, rCut = R_CUT;
  const on2 = rOn * rOn, cut2 = rCut * rCut;
  const denom = (cut2 - on2) ** 3;
  const num = (cut2 - rsq) * (cut2 - rsq) * (cut2 + 2 * rsq - 3 * on2);
  return num / denom;
}

/**
 * Derivative dS/dr of the switch (negative in the ramp region), used to keep
 * forces consistent with the switched energy.
 */
export function switchDeriv(r) {
  if (r <= R_SWITCH_ON || r >= R_CUT) return 0;
  const rsq = r * r, rOn = R_SWITCH_ON, rCut = R_CUT;
  const on2 = rOn * rOn, cut2 = rCut * rCut;
  const denom = (cut2 - on2) ** 3;
  // d/dr [ (cut2−r²)²(cut2+2r²−3on2) ]
  const a = (cut2 - rsq) * (cut2 - rsq);          // (cut2−r²)²
  const b = cut2 + 2 * rsq - 3 * on2;             // (cut2+2r²−3on2)
  const da = -4 * r * (cut2 - rsq);               // d/dr[a]
  const db = 4 * r;                                // d/dr[b]
  return (da * b + a * db) / denom;
}

// Metal coordination spring constant (kcal/mol/Å2) and its target distance:
// the metal's coordR (so the spring holds the ion at its native coordination
// distance rather than collapsing it onto the donor).
export const METAL_K = 40.0;

// LJ params for heavy atoms that are NOT in LIG_ELEMENT (e.g. protein atoms
// are all C/N/O/S; the table already covers them, but metals and any exotic
// element fall back here).
const HEAVY_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0 };

/** Solvent / crystallographic water residue names dropped by the parser. */
const SOLVENT = new Set([
  "HOH", "WAT", "H2O", "DOD", "HHO", "TIP", "TIP3", "TIP3P", "SPC", "SPCE",
]);

/**
 * Parse a PDB text into a flat heavy-atom system (protein + ligand + ions).
 *
 * Every ATOM / HETATM heavy atom becomes one particle; the crystallographic
 * connectivity is NOT used to define bonds (it is frequently incomplete) —
 * bonds are rebuilt geometrically by buildTopology().
 *
 * @param {string} pdbText  PDB-format text (ATOM/HETATM/CONECT records)
 * @returns {{atoms: Array<{x,y,z,element,resName,chain,resSeq,atomName,serial,isProtein,isMetal}>, n:number}}
 */
export function parseHeavy(pdbText) {
  const atoms = [];
  const seen = new Set();
  const bySerial = new Map();
  let n = 0;

  for (const line of pdbText.split(/\r?\n/)) {
    const rec = line.slice(0, 6);
    if (rec !== "ATOM  " && rec !== "HETATM") continue;
    const atomName = line.slice(12, 16).trim();
    const resName = line.slice(17, 20).trim();
    if (SOLVENT.has(resName)) continue;   // drop water + free ions (they are not part of the complex topology)

    const chain = (line.charAt(21) || " ").trim() || "_";
    const resSeq = parseInt(line.slice(22, 26), 10);
    const iCode = line.charAt(26).trim();
    const x = parseFloat(line.slice(30, 38));
    const y = parseFloat(line.slice(38, 46));
    const z = parseFloat(line.slice(46, 54));
    if (Number.isNaN(resSeq) || Number.isNaN(x + y + z)) continue;

    let element = line.slice(76, 78).trim().toUpperCase();
    if (!element) element = elementFromName(atomName);
    if (!element || element === "H") continue;   // united-atom: drop hydrogen

    const serial = parseInt(line.slice(6, 11), 10);
    const key = `${chain}|${resSeq}|${iCode}|${atomName}`;
    if (seen.has(key)) continue;   // alt-loc duplicates
    seen.add(key);

    const isProtein = rec === "ATOM  ";
    const isMetal = !!METAL_ELEMENT[element];
    const a = { x, y, z, element, atomName, resName, chain, resSeq, serial, isProtein, isMetal };
    atoms.push(a);
    bySerial.set(serial, a);
    n++;
  }

  if (atoms.length === 0) throw new Error("No heavy atoms found — is this a valid PDB file?");
  return { atoms, n, bySerial };
}

/** Infer an element from a PDB atom name (e.g. "CA" backbone α-carbon vs Ca ion). */
function elementFromName(name) {
  if (!name) return null;
  // Backbone + side-chain standard atom names: strip leading digits, keep the
  // first 1-2 letters. "CA" is genuinely ambiguous (α-carbon vs calcium) but
  // this path only runs when the element column is blank, which is rare.
  const m = name.match(/^[0-9]*([A-Za-z]{1,2})/);
  if (!m) return null;
  let el = m[1].toUpperCase();
  if (el.length === 2) {
    // Two-letter element symbols (CL, BR, ZN, FE, MG, CA...) are kept as-is;
    // everything else is a single-letter element followed by a digit/suffix.
    const two = new Set(["CL", "BR", "ZN", "FE", "MG", "CA", "CU", "MN", "NI", "CO", "NA", "K", "SE", "SI", "AL"]);
    if (two.has(el)) return el;
    return el[0];
  }
  return el;
}

/**
 * Build a covalent bond topology from parsed atoms.
 * Two heavy atoms are bonded when their distance ≤ BOND_SLACK × (sum of their
 * covalent radii). Metals are excluded from covalent bonding (they coordinate
 * instead, see buildMetalCoordination()).
 *
 * @param {Array} atoms   output of parseHeavy()
 * @returns {{bonds:Array<[number,number]>, angles:Array<[number,number,number]>,
 *            impropers:Array<[number,number,number,number]>,
 *            propers:Array<[number,number,number,number]>}}
 */
export function buildTopology(atoms) {
  const n = atoms.length;
  const bonds = [];
  const nbond = new Array(n).fill(0).map(() => []);   // neighbor list

  // --- 1. Geometric bond detection (excluding metals) ----------------------
  for (let i = 0; i < n; i++) {
    const A = atoms[i];
    if (A.isMetal) continue;
    const rA = COVALENT_RADIUS[A.element] ?? 0.77;
    for (let j = i + 1; j < n; j++) {
      const B = atoms[j];
      if (B.isMetal) continue;
      const rB = COVALENT_RADIUS[B.element] ?? 0.77;
      const dx = A.x - B.x, dy = A.y - B.y, dz = A.z - B.z;
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (r < BOND_SLACK * (rA + rB) && r < 2.2) {   // hard cap avoids spurious long bonds
        bonds.push([i, j]);
        nbond[i].push(j);
        nbond[j].push(i);
      }
    }
  }

  // --- 2. Angles: 1-2-3 triples (i-j-k with j bonded to both) --------------
  const angles = [];
  const angleSet = new Set();
  for (let j = 0; j < n; j++) {
    const nb = nbond[j];
    for (let a = 0; a < nb.length; a++) {
      for (let b = a + 1; b < nb.length; b++) {
        const i = nb[a], k = nb[b];
        const key = i < k ? `${i}-${j}-${k}` : `${k}-${j}-${i}`;
        if (angleSet.has(key)) continue;
        angleSet.add(key);
        angles.push([i, j, k]);
      }
    }
  }

  // --- 3. Proper dihedrals: 1-2-3-4 (i-j-k-l) ------------------------------
  const propers = [];
  const properSet = new Set();
  for (let j = 0; j < n; j++) {
    for (const i of nbond[j]) {
      for (const k of nbond[j]) {
        if (k === i) continue;
        for (const l of nbond[k]) {
          if (l === i || l === j) continue;
          const key = [i, j, k, l].join("-");
          if (properSet.has(key)) continue;
          properSet.add(key);
          propers.push([i, j, k, l]);
        }
      }
    }
  }

  // --- 4. Impropers (out-of-plane): central atom j with 3 neighbors --------
  // j bonded to i, k, l (all three) → planar-center improper.
  const impropers = [];
  for (let j = 0; j < n; j++) {
    const nb = nbond[j];
    if (nb.length < 3) continue;
    for (let a = 0; a < nb.length; a++) {
      for (let b = a + 1; b < nb.length; b++) {
        for (let c = b + 1; c < nb.length; c++) {
          impropers.push([nb[a], j, nb[b], nb[c]]);
        }
      }
    }
  }

  return { bonds, angles, propers, impropers };
}

/**
 * Detect metal–donor coordination and return a list of [metalIdx, donorIdx]
 * pairs (the metal's N highest-affinity N/O/S donors within coordR, capped at
 * coordN). Used to build coordination springs.
 *
 * @param {Array} atoms  output of parseHeavy()
 * @returns {Array<[number,number]>} metal–donor index pairs
 */
export function buildMetalCoordination(atoms) {
  const pairs = [];
  for (let i = 0; i < atoms.length; i++) {
    const M = atoms[i];
    if (!M.isMetal) continue;
    const p = METAL_ELEMENT[M.element] ?? METAL_ELEMENT_DEFAULT;
    const donors = [];
    for (let j = 0; j < atoms.length; j++) {
      if (j === i) continue;
      const D = atoms[j];
      if (D.isMetal) continue;
      if (!(D.element === "N" || D.element === "O" || D.element === "S")) continue;
      const dx = M.x - D.x, dy = M.y - D.y, dz = M.z - D.z;
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (r < p.coordR) donors.push({ j, r });
    }
    donors.sort((a, b) => a.r - b.r);
    const take = Math.min(p.coordN, donors.length);
    for (let k = 0; k < take; k++) pairs.push([i, donors[k].j]);
  }
  return pairs;
}

/**
 * Build the flat interaction lists HeavyForceField uses, mirroring the layout
 * of the CG ForceField (bonds [i,j,r0×3], angles [i,j,k,θ0×4], impropers
 * [i,j,k,l,φ0×5]) so the shared ff-harmonic kernels can be reused verbatim.
 *
 * @param {Array} atoms      output of parseHeavy()
 * @param {object} topo      output of buildTopology()
 * @param {Array} coordPairs output of buildMetalCoordination()
 * @returns {{bonds, angles, impropers, propers, coord, holoSprings, nHolo}}
 */
export function buildLists(atoms, topo, coordPairs) {
  const bonds = [];
  for (const [i, j] of topo.bonds) {
    const A = atoms[i], B = atoms[j];
    const r0 = Math.hypot(A.x - B.x, A.y - B.y, A.z - B.z);
    bonds.push(i, j, r0);
  }

  const angles = [];
  for (const [i, j, k] of topo.angles) {
    const A = atoms[i], B = atoms[j], C = atoms[k];
    const th0 = angleAt(atoms, i, j, k);
    angles.push(i, j, k, th0);
  }

  const impropers = [];
  for (const [i, j, k, l] of topo.impropers) {
    const phi0 = improperAngleFlat(atoms, i, j, k, l);
    impropers.push(i, j, k, l, phi0);
  }

  // Proper dihedrals are stored as [i,j,k,l,φ0] for the torsion kernel.
  const propers = [];
  for (const [i, j, k, l] of topo.propers) {
    const phi0 = improperAngleFlat(atoms, i, j, k, l);
    propers.push(i, j, k, l, phi0);
  }

  // Coordination springs as [i,j,r0×3] harmonic pairs (metal, donor).
  const coord = [];
  for (const [i, j] of coordPairs) {
    const A = atoms[i], B = atoms[j];
    const r0 = Math.hypot(A.x - B.x, A.y - B.y, A.z - B.z);
    coord.push(i, j, r0);
  }

  return { bonds, angles, impropers, propers, coord, holoSprings: new Float64Array(0), nHolo: 0 };
}

/** Interior angle θ at atom j between (i,j) and (j,k). */
export function angleAt(atoms, i, j, k) {
  const A = atoms[i], B = atoms[j], C = atoms[k];
  const ax = A.x - B.x, ay = A.y - B.y, az = A.z - B.z;
  const bx = C.x - B.x, by = C.y - B.y, bz = C.z - B.z;
  const la = Math.hypot(ax, ay, az) || 1e-12;
  const lb = Math.hypot(bx, by, bz) || 1e-12;
  let c = (ax * bx + ay * by + az * bz) / (la * lb);
  c = Math.min(1, Math.max(-1, c));
  return Math.acos(c);
}

/** Dihedral φ for a flat position array (self-contained, no typed pos). */
function improperAngleFlat(atoms, i, j, k, l) {
  const p = new Float64Array(atoms.length * 3);
  for (let m = 0; m < atoms.length; m++) {
    p[3 * m] = atoms[m].x; p[3 * m + 1] = atoms[m].y; p[3 * m + 2] = atoms[m].z;
  }
  return improperAngle(p, i, j, k, l);
}

/**
 * HeavyForceField — ForceField-compatible all-atom force field.
 *
 * Energy terms (units kcal/mol, Å):
 *   U = Σ_b ½ k_b (r − r0)²            covalent bonds        (k_b = 200)
 *     + Σ_θ ½ k_θ (θ − θ0)²            angles               (k_θ = 40)
 *     + Σ_φ ½ k_φ (φ − φ0)²            impropers            (k_φ = 20)
 *     + Σ_ψ ½ k_ψ (ψ − ψ0)²            proper dihedrals     (k_ψ = 2)
 *     + Σ_m ½ k_m (r − r0)²            metal coordination   (k_m = 40)
 *     + Σ_{i<j non-excl} 4ϵ[(σ/r)¹²−(σ/r)⁶]                  LJ
 *     + Σ_{i<j non-excl} q_i q_j K_ELEC exp(−r/λ)/r         screened Coulomb
 *
 * The public interface mirrors ForceField so the integrator, funnel, viewer,
 * and UI require no changes.
 */
export class HeavyForceField {
  /**
   * @param {object} system   { atoms } from parseHeavy() — or a full {atoms, n}
   * @param {object} par      { rc, gamma, ... } (unused here; kept for symmetry)
   * @param {Array}  ligands  reserved for API symmetry (heavy mode builds its
   *                          own topology; pass [] )
   */
  constructor(system, par = {}, ligands = []) {
    const atoms = system.atoms;
    this.n = atoms.length;
    this.nProt = atoms.filter((a) => a.isProtein).length;
    this.nLigAtoms = this.n - this.nProt;

    // Element → LJ/charge params (protein atoms reuse LIG_ELEMENT since the
    // elements are the same C/N/O/S).
    this._elem = new Array(this.n);
    let totalQ = 0;
    for (let i = 0; i < this.n; i++) {
      const el = atoms[i].element;
      if (METAL_ELEMENT[el]) this._elem[i] = { ...METAL_ELEMENT[el] };
      else this._elem[i] = { ...(LIG_ELEMENT[el] ?? LIG_ELEMENT_DEFAULT) };
      totalQ += this._elem[i].q;
    }
    // Net-charge neutralization: the CG ligand charges (q(O)=−0.5, q(N)=−0.3)
    // were tuned for a single ligand, but applied to ~thousands of protein
    // atoms they sum to a large net charge with no counterions, which would
    // dominate every energy term. Subtract the mean charge from every atom so
    // the system is monopole-free (relative charges — including the metal's
    // +2 — are preserved) and the electrostatics is a local dipole/quadrupole
    // flavor term, as it should be in implicit-solvent protein modeling.
    const qCorr = totalQ / this.n;
    for (let i = 0; i < this.n; i++) this._elem[i].q -= qCorr;

    // Masses (Da): heavy atoms by element (protein ~full atom masses).
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

    // Topology
    const topo = buildTopology(atoms);
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

    // Non-bonded exclusions: 1-2 (bonded) and 1-3 (angle) pairs are excluded;
    // 1-4 (dihedral) pairs get a scale factor.
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

    // Building blocks: metal coordination pairs are harmonic 1-2 springs but
    // are NOT excluded from non-bonded (so the ion also feels LJ around it).
    // Covalent bonds are exposed as `ligandBonds` so the viewer draws the
    // full all-atom bond graph and the integrator's stiff-bond stability
    // guard sees the fastest protein/ligand vibrations (k = 200 backbone
    // C–N, aromatic 200, single 300-estimates) instead of the empty CG list.
    this.ligandBonds = new Float64Array(L.bonds);
    this.ligandAtoms = null;
    this.springs = new Float64Array(0);
    this.springK = new Float64Array(0);
    this.nativeContacts = new Float64Array(0);
    this.kBond = 200.0;
    this.kAngle = 40.0;
    this.kImproper = 20.0;
    this.kProper = 2.0;
    this.metalK = METAL_K;

    // Funnel (optional, same as CG)
    this.funnel = null;
    this.funnelOn = false;

    this.bindingU = 0;
    this.desolvU = 0;
    this.repU = 0;
    this.bondU = 0;
    this.angleU = 0;
    this.improperU = 0;
    this.properU = 0;
    this.coordU = 0;
    this.elecU = 0;
    this._nanStrikes = 0;
  }

  /** Attach an optional Funnel instance (interface symmetry with ForceField). */
  setFunnel(fn) { this.funnel = fn; }

  /** CG-only spring-scaling API — no-op in heavy mode (kept for symmetry). */
  setSpringScale() {}
  clearSpringScale() {}
  rebuildHoloSprings() {}

  /** Total potential energy (kcal/mol) given positions; fills this.forces. */
  compute(pos) {
    const f = this.forces;
    f.fill(0);
    let U = 0;

    // 1. Covalent bonds (k=200) + metal coordination (k=40)
    this.bondU = harmonicFlat(pos, f, this.bonds, 3, this.kBond);
    this.coordU = harmonicFlat(pos, f, this.coord, 3, this.metalK);
    U += this.bondU + this.coordU;

    // 2. Angles (k=40)
    this.angleU = angleFlat(pos, f, this.angles, 4, this.kAngle);
    U += this.angleU;

    // 3. Impropers (k=20) + propers (k=2)
    this.improperU = improperFlat(pos, f, this.impropers, 5, this.kImproper);
    this.properU = properFlat(pos, f, this.propers, 5, this.kProper);
    U += this.improperU + this.properU;

    // 4. Non-bonded: LJ + screened Coulomb over non-excluded pairs
    const nb = this._nonBonded(pos, f);
    this.elecU = nb.elec;
    this.repU = nb.lj;
    U += nb.lj + nb.elec;

    // 5. Funnel bias (optional)
    if (this.funnel && this.funnelOn) U += this.funnel.addForces(pos, f);

    // NaN guard
    if (!Number.isFinite(U)) {
      for (let i = 0; i < f.length; i++) if (!Number.isFinite(f[i])) f[i] = 0;
      this._nanStrikes++;
      U = NaN;
    }
    this.energy = U;
    return U;
  }

  /** Spatial-hash-free non-bonded: O(n²) — fine for the sizes heavy mode targets. */
  /** Non-bonded: LJ + screened Coulomb over non-excluded pairs, with an
   *  AMBER-style switching cutoff (R_CUT) and distance-dependent dielectric
   *  ε(r) = r. O(n²) — fine for the sizes heavy mode targets. */
  _nonBonded(pos, f) {
    const n = this.n;
    const nProt = this.nProt;
    let lj = 0, elec = 0, bindLJ = 0, bindElec = 0;
    for (let i = 0; i < n - 1; i++) {
      const xi = 3 * i;
      const ei = this._elem[i];
      const qi = ei.q;
      const si = ei.sigma, epsi = ei.eps;
      for (let j = i + 1; j < n; j++) {
        const k = pairKey(i, j);
        if (this._excluded.has(k)) continue;
        const s14 = this._scale14.get(k) ?? 1.0;
        const xj = 3 * j;
        const dx = pos[xj] - pos[xi];
        const dy = pos[xj + 1] - pos[xi + 1];
        const dz = pos[xj + 2] - pos[xi + 2];
        const r2 = dx * dx + dy * dy + dz * dz;
        const r = Math.sqrt(r2) || 1e-12;
        if (r >= R_CUT) continue;   // switching cutoff

        // LJ (geometric mixing)
        const ej = this._elem[j];
        const s = 0.5 * (si + ej.sigma);
        const eps = Math.sqrt(epsi * ej.eps);
        const sr = s / r, sr6 = sr * sr * sr * sr * sr * sr;
        const ljE = 4 * eps * (sr6 * sr6 - sr6);
        const ljF = 4 * eps * (12 * sr6 * sr6 - 6 * sr6) / r;

        // Screened Coulomb with distance-dependent dielectric ε(r) = r:
        //   U = q_i q_j K_ELEC exp(−r/λ) / r2
        let ee = 0, ef = 0;
        if (qi !== 0 && ej.q !== 0) {
          const qq = qi * ej.q;
          const fac = Math.exp(-r / SCREEN_LEN) / (r * r);
          ee = K_ELEC * qq * fac;
          // dU/dr = U·(−1/λ − 2/r)  ⇒  ef = dU/dr (used with F = −dU/dr·r̂)
          ef = ee * (-1 / SCREEN_LEN - 2 / r);
        }

        // Apply the switching function to energy AND force. With U = S·U0,
        // F = −dU/dx = −(S·dU0/dx + U0·dS/dx). Here ljF = −dU_lj/dr is a
        // repulsive force magnitude while ef = +dU_elec/dr is a derivative, so
        // the combined radial term is (−S·ljF + S·ef + dS·U0) and the Cartesian
        // force is F_i = (totF)·(dx/r) with r̂ pointing from i to j.
        const S = switchFunc(r);
        const dS = switchDeriv(r);
        const totE = s14 * (S * ljE + S * ee);
        const totF = s14 * (-S * ljF + S * ef + dS * (ljE + ee));
        lj += s14 * S * ljE;
        elec += s14 * S * ee;
        // Protein–ligand cross term (i < nProt ≤ j): the bindingU the HUD's
        // U_bind readout and the NN pose scorer consume (mirrors ForceField's
        // protein–ligand nonbonded energy).
        if (i < nProt && j >= nProt) {
          bindLJ += s14 * S * ljE;
          bindElec += s14 * S * ee;
        }

        // F_i = −dU/dx_i = (totF)·(dx/r)  (dx = x_j − x_i ⇒ r̂ from i to j)
        const fx = totF * dx / r, fy = totF * dy / r, fz = totF * dz / r;
        f[xi] += fx; f[xi + 1] += fy; f[xi + 2] += fz;
        f[xj] -= fx; f[xj + 1] -= fy; f[xj + 2] -= fz;
      }
    }
    this.bindingU = bindLJ + bindElec;
    return { lj, elec };
  }

  /** Instantaneous kinetic temperature (K); matches ForceField.kineticTemp. */
  kineticTemp(vel, mass) {
    let ke = 0;
    const m = mass || this.masses;
    for (let i = 0; i < this.n * 3; i++) ke += m[(i / 3) | 0] * vel[i] * vel[i];
    ke *= 0.5 / KCONV;
    return ke / (1.5 * this.n * KB_KCAL);
  }

  /** RMSD (Å) of pos vs native ref over all heavy atoms. */
  rmsd(pos) {
    let s = 0;
    for (let i = 0; i < pos.length; i++) {
      const d = pos[i] - this.ref[i];
      s += d * d;
    }
    return Math.sqrt(s / this.n);
  }
}

/** Harmonic pair kernel (bonds + coordination) — flat [i,j,r0] list. */
function harmonicFlat(pos, f, list, stride, k) {
  let U = 0;
  for (let a = 0; a < list.length; a += stride) {
    const i = 3 * list[a], j = 3 * list[a + 1], r0 = list[a + 2];
    const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
    const dr = r - r0;
    U += 0.5 * k * dr * dr;
    const s = (k * dr) / r;
    const fx = s * dx, fy = s * dy, fz = s * dz;
    f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
    f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
  }
  return U;
}

/** Harmonic angle kernel — flat [i,j,k,θ0] list. */
function angleFlat(pos, f, list, stride, k) {
  let U = 0;
  for (let a = 0; a < list.length; a += stride) {
    const i = 3 * list[a], j = 3 * list[a + 1], kk = 3 * list[a + 2], th0 = list[a + 3];
    const ax = pos[i] - pos[j], ay = pos[i + 1] - pos[j + 1], az = pos[i + 2] - pos[j + 2];
    const bx = pos[kk] - pos[j], by = pos[kk + 1] - pos[j + 1], bz = pos[kk + 2] - pos[j + 2];
    const la = Math.hypot(ax, ay, az) || 1e-12;
    const lb = Math.hypot(bx, by, bz) || 1e-12;
    let c = (ax * bx + ay * by + az * bz) / (la * lb);
    c = Math.min(1, Math.max(-1, c));
    const th = Math.acos(c);
    const dth = th - th0;
    U += 0.5 * k * dth * dth;
    const sinTh = Math.sqrt(Math.max(1e-12, 1 - c * c));
    const pref = (k * dth) / sinTh;
    const ga = 1 / la, gb = 1 / lb;
    const axy = ax * ga, ayy = ay * ga, azy = az * ga;
    const bxy = bx * gb, byy = by * gb, bzy = bz * gb;
    let fix = pref * (bxy - c * axy) * ga;
    let fiy = pref * (byy - c * ayy) * ga;
    let fiz = pref * (bzy - c * azy) * ga;
    let fkx = pref * (axy - c * bxy) * gb;
    let fky = pref * (ayy - c * byy) * gb;
    let fkz = pref * (azy - c * bzy) * gb;
    f[i] += fix; f[i + 1] += fiy; f[i + 2] += fiz;
    f[kk] += fkx; f[kk + 1] += fky; f[kk + 2] += fkz;
    f[j] -= fix + fkx; f[j + 1] -= fiy + fky; f[j + 2] -= fiz + fkz;
  }
  return U;
}

/** Improper torsions via finite differences of improperAngle (shared kernel). */
function improperFlat(pos, f, list, stride, k) {
  const h = 1e-5;
  let U = 0;
  for (let a = 0; a < list.length; a += stride) {
    const i = list[a], j = list[a + 1], kk = list[a + 2], l = list[a + 3], phi0 = list[a + 4];
    const phi = improperAngle(pos, i, j, kk, l);
    U += 0.5 * k * (phi - phi0) * (phi - phi0);
    const pref = -k * (phi - phi0);
    for (const m of [i, j, kk, l]) {
      const c = 3 * m;
      for (let ax = 0; ax < 3; ax++) {
        const ci = c + ax, save = pos[ci];
        pos[ci] = save + h;
        const phiP = improperAngle(pos, i, j, kk, l);
        pos[ci] = save - h;
        const phiM = improperAngle(pos, i, j, kk, l);
        pos[ci] = save;
        f[ci] += pref * (phiP - phiM) / (2 * h);
      }
    }
  }
  return U;
}

/** Proper dihedrals: same functional form as impropers but with kProper. */
function properFlat(pos, f, list, stride, k) {
  return improperFlat(pos, f, list, stride, k);
}

/** Heavy-atom mass (Da) by element; falls back to a reasonable default. */
export function heavyMass(element) {
  switch (element) {
    case "C": return 12.011;
    case "N": return 14.007;
    case "O": return 15.999;
    case "S": return 32.06;
    case "P": return 30.974;
    case "F": return 18.998;
    case "CL": return 35.45;
    case "BR": return 79.904;
    case "I": return 126.904;
    case "ZN": return 65.38;
    case "FE": return 55.845;
    case "MG": return 24.305;
    case "CA": return 40.078;
    case "CU": return 63.546;
    case "MN": return 54.938;
    case "NI": return 58.693;
    case "CO": return 58.933;
    case "NA": return 22.99;
    case "K": return 39.098;
    default: return 14.0;
  }
}

/** Canonical pair key for exclusion/1-4 maps. */
function pairKey(i, j) { return i < j ? i * 1e6 + j : j * 1e6 + i; }