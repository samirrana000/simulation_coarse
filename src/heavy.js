/**
 * heavy.js — All-atom "heavy mode" force field with biophysics & high performance.
 *
 * Capabilities:
 *  - Explicit heavy atoms for protein, ligands, cofactors, and metal ions
 *  - Geometric covalent bond detection + metal coordination springs (N/O/S donors)
 *  - Fast O(N) spatial hashing grid for non-bonded interactions (LJ + GB + electrostatics)
 *  - Fast analytic proper and improper dihedral gradients (Blondel-Karplus / Bekker vector formulas)
 *  - Generalized Born + Debye-Hückel implicit solvent & AMBER partial charges
 *  - Hydrophobic SASA burial & directional H-bonding
 *  - Full compatibility with MOL2 ligand files, library placement, and PDB HETATM
 *
 * Handshake with CG mode: see docs/CG_HEAVY.md — CG (Cα ENM) and heavy
 * share the same reference coordinates (Å) and KCONV/KB_KCAL unit contract
 * (src/units.js:418.4). Switching is `selectSystem` vs `selectHeavy` +
 * `ForceField` vs `HeavyForceField`; viewer center/radius is protein-only
 * in both modes (src/viewer.js:163-182) so a distant ligand does not
 * inflate the camera (docs/CG_HEAVY.md §Handshake center/radius fix).
 */

import {
  METAL_ELEMENT, METAL_ELEMENT_DEFAULT,
  COVALENT_RADIUS, BOND_SLACK,
  LIG_ELEMENT, LIG_ELEMENT_DEFAULT,
  KB_KCAL, KCONV,
} from "./ff-params.js?v=10";
import { springForces, dihedralForcesAnalytic } from "./ff-harmonic.js?v=10";
import { assignCharges, GB_RADII } from "./physics/charges.js?v=10";
import { GeneralizedBorn, COULOMB_CONST } from "./physics/gb.js?v=10";
import { SasaModel } from "./physics/sasa.js?v=10";
import { DirectionalHBond } from "./physics/hbond.js?v=10";
import { SpatialGrid } from "./spatial-grid.js?v=10";
import { getBondParams, getAngleParams, getNonbondedParams } from "./physics/forcefield/amber14sb.js?v=10";
import { computeBornRadii as computeOBC2Radii, gbEnergyForces as gbOBC2Forces, debyeKappa } from "./physics/solvation/gb_obc2.js?v=10";
import { lcpoSasa } from "./physics/solvation/lcpo_sasa.js?v=10";
import { membraneEnergyForces, transferDgFor } from "./physics/solvation/membrane_slab.js?v=10";
import { typeMolecule, assignCharges as gaffAssignCharges } from "./chem/gaff2_mapper.js?v=10";
import {
  piStackForces, cationPiForces, halogenForces,
  buildRingFrames, buildCationList, buildHalogenList, HALOGEN_EPS,
} from "./physics/weakint.js?v=10";
import { detectCoordination, enforceCoordination } from "./chem/metals.js?v=10";

export const K_ELEC = 332.0;
export const SCREEN_LEN = 8.0; // Å
export const R_CUT = 8.5;      // Å
export const R_SWITCH_ON = 6.5;// Å

export function switchFunc(r) {
  if (r <= R_SWITCH_ON) return 1;
  if (r >= R_CUT) return 0;
  const rsq = r * r, rOn = R_SWITCH_ON, rCut = R_CUT;
  const on2 = rOn * rOn, cut2 = rCut * rCut;
  const denom = (cut2 - on2) ** 3;
  const num = (cut2 - rsq) * (cut2 - rsq) * (cut2 + 2 * rsq - 3 * on2);
  return num / denom;
}

export function switchDeriv(r) {
  if (r <= R_SWITCH_ON || r >= R_CUT) return 0;
  const rsq = r * r, rOn = R_SWITCH_ON, rCut = R_CUT;
  const on2 = rOn * rOn, cut2 = rCut * rCut;
  const denom = (cut2 - on2) ** 3;
  const a = (cut2 - rsq) * (cut2 - rsq);
  const b = cut2 + 2 * rsq - 3 * on2;
  const da = -4 * r * (cut2 - rsq);
  const db = 4 * r;
  return (da * b + a * db) / denom;
}

export const METAL_K = 40.0;
const HEAVY_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0 };

const SOLVENT = new Set([
  "HOH", "WAT", "H2O", "DOD", "HHO", "TIP", "TIP3", "TIP3P", "SPC", "SPCE", "SOL",
]);

/**
 * Count validation warnings helper (A06).
 * @param {string[]|object} warnings
 * @returns {number}
 */
export function countWarnings(warnings) {
  if (Array.isArray(warnings)) return warnings.length;
  if (warnings && Array.isArray(warnings.warnings)) return warnings.warnings.length;
  return 0;
}

/**
 * Parse a PDB text into a flat heavy-atom system with classified hetero groups.
 *
 * Input validation / warnings (A06): malformed lines are skipped with
 * `warnings` + `console.warn`; result includes `.warnings` (string[]).
 * Use `countWarnings(result.warnings)` to count.
 *
 * @param {string} pdbText
 * @returns {{atoms: Array, n: number, bySerial: Map, heteroGroups: Array, warnings: string[]}}
 */
export function parseHeavy(pdbText) {
  const atoms = [];
  const seen = new Set();
  const bySerial = new Map();
  const heteroGroupMap = new Map();
  const warnings = [];
  let n = 0;

  for (const line of pdbText.split(/\r?\n/)) {
    const rec = line.slice(0, 6);
    if (rec !== "ATOM  " && rec !== "HETATM") continue;
    const atomName = line.slice(12, 16).trim();
    const resName = line.slice(17, 20).trim();
    if (SOLVENT.has(resName)) continue;

    const chain = (line.charAt(21) || " ").trim() || "_";
    const resSeq = parseInt(line.slice(22, 26), 10);
    const iCode = line.charAt(26).trim();
    const x = parseFloat(line.slice(30, 38));
    const y = parseFloat(line.slice(38, 46));
    const z = parseFloat(line.slice(46, 54));
    if (Number.isNaN(resSeq) || Number.isNaN(x + y + z)) {
      const msg = `parseHeavy: malformed line skipped (resSeq=${resSeq} x=${x}) line="${line.slice(0, 66).trim()}"`;
      warnings.push(msg); console.warn(msg);
      continue;
    }

    let element = line.slice(76, 78).trim().toUpperCase();
    if (!element) element = elementFromName(atomName, rec);
    if (!element || element === "H") continue;

    const serial = parseInt(line.slice(6, 11), 10);
    const key = `${chain}|${resSeq}|${iCode}|${atomName}`;
    if (seen.has(key)) {
      const msg = `parseHeavy: duplicate atom ${key} skipped`;
      warnings.push(msg); console.warn(msg);
      continue;
    }
    seen.add(key);

    const isProtein = rec === "ATOM  ";
    const isHetero = !isProtein;
    const isMetal = !!METAL_ELEMENT[element];
    const heteroKey = isHetero ? `${chain}|${resSeq}|${resName}` : null;
    const a = {
      x, y, z, element, atomName, resName, chain, resSeq, serial,
      isProtein, isMetal, isHetero, isWater: false, heteroKey,
    };
    atoms.push(a);
    bySerial.set(serial, a);
    n++;

    if (isHetero) {
      let g = heteroGroupMap.get(heteroKey);
      if (!g) {
        g = { key: heteroKey, resName, chain, resSeq, atomIndices: [], isMetal: false, elements: [] };
        heteroGroupMap.set(heteroKey, g);
      }
      g.atomIndices.push(atoms.length - 1);
      g.elements.push(element);
      if (isMetal) g.isMetal = true;
    }
  }

  if (atoms.length === 0) throw new Error("No heavy atoms found in PDB file.");

  const heteroGroups = [...heteroGroupMap.values()].map((g) => ({
    key: g.key,
    resName: g.resName,
    chain: g.chain,
    resSeq: g.resSeq,
    atomIndices: g.atomIndices,
    isMetal: g.isMetal,
    element: g.elements[0] ?? null,
  }));

  if (warnings.length) console.warn(`[parseHeavy] ${warnings.length} warning(s) total`);
  return { atoms, n, bySerial, heteroGroups, warnings };
}

/**
 * Filter parseHeavy() output by chain / residue range and hetero-group selection.
 */
export function selectHeavy(parsedHeavy, { chains = null, resFrom = null, resTo = null, heteroSelection = null, includePdbLigands = true, hasExternalLigand = false } = {}) {
  const inRange = (a) =>
    (!chains || chains.includes(a.chain)) &&
    (resFrom === null || a.resSeq >= resFrom) &&
    (resTo === null || a.resSeq <= resTo);
  const proteinAtoms = parsedHeavy.atoms.filter((a) => a.isProtein && inRange(a));
  const heteroAtoms = parsedHeavy.atoms.filter((a) =>
    !a.isProtein &&
    inRange(a) &&
    (heteroSelection == null || heteroSelection[a.heteroKey] === true)
  );

  const metalsAndCofactors = [];
  const pdbLigandAtoms = [];

  for (const a of heteroAtoms) {
    if (a.isMetal || hasExternalLigand || !includePdbLigands) {
      metalsAndCofactors.push({ ...a, isLigand: false });
    } else {
      pdbLigandAtoms.push({ ...a, isLigand: true });
    }
  }

  const atoms = proteinAtoms.concat(metalsAndCofactors).concat(pdbLigandAtoms);
  const beads = atoms.map((a) => ({ ...a, x: a.x, y: a.y, z: a.z }));
  return { atoms, beads, segments: [], heavy: true, pdbLigandAtoms };
}

/**
 * Append resolved external ligand molecules to selection.
 *
 * Phase 2 opt-in: opts.gaff === true runs the GAFF2-lite typer
 * (src/chem/gaff2_mapper.js typeMolecule + Gasteiger/AM1-BCC-lite charges)
 * per molecule and stores the result on the appended atoms
 * (atom.charge, atom.gaffType). Any per-molecule failure falls back to the
 * parsed charges so a bad ligand can never break system construction.
 * Default (opts.gaff falsy) preserves legacy behavior exactly.
 */
export function appendHeavyLigands(sel, molecules, opts = {}) {
  if (!molecules || molecules.length === 0) return sel;
  const atoms = sel.atoms.slice();
  const beads = sel.beads.slice();
  let resSeq = 1;
  for (const mol of molecules) {
    const resName = (mol.resName || "LIG").slice(0, 3).toUpperCase();
    const chain = mol.chain || "L";
    if (opts.gaff) {
      try {
        typeMolecule(mol.atoms, mol.bonds ?? []);
        gaffAssignCharges(mol.atoms, mol.bonds ?? [], { writeBack: true });
      } catch (e) {
        console.warn(`[appendHeavyLigands] GAFF2 fallback failed for ${resName} (${e?.message}) — parsed charges kept`);
      }
    }
    for (const at of mol.atoms) {
      const atom = {
        x: at.x, y: at.y, z: at.z,
        element: at.element,
        resName, chain, resSeq,
        atomName: at.atomName || at.element,
        serial: at.serial ?? 0,
        // GAFF2 fields exist only on the opt-in path so legacy charge
        // assignment (element defaults via assignCharges) is untouched.
        ...(opts.gaff ? { charge: at.charge ?? 0, gaffType: at.gaffType ?? null } : {}),
        isProtein: false, isMetal: false, isLigand: true, isHetero: false, heteroKey: null,
      };
      atoms.push(atom);
      beads.push({ ...atom });
    }
    resSeq++;
  }
  return { atoms, beads, segments: sel.segments, heavy: true };
}

function elementFromName(name, rec) {
  if (!name) return null;
  const m = name.match(/^[0-9]*([A-Za-z]{1,2})/);
  if (!m) return null;
  let el = m[1].toUpperCase();
  // Guard: ATOM CA is alpha carbon, not calcium (HETATM CA is calcium)
  if (rec === "ATOM  " && name.trim() === "CA") return "C";
  if (el.length === 2) {
    const two = new Set(["CL", "BR", "ZN", "FE", "MG", "CA", "CU", "MN", "NI", "CO", "NA", "K", "SE", "SI", "AL"]);
    if (two.has(el)) return el;
    return el[0];
  }
  return el;
}

/**
 * Build covalent topology from geometry.
 *
 * Covalent radii source: CSD surveys (Allen et al., Cordero et al. Dalton
 * 2008) / Bondi 1964 — see src/ff-params.js:COVALENT_RADIUS table. Two heavy
 * atoms i,j are bonded iff  r_ij < BOND_SLACK*(r_cov(i)+r_cov(j))  AND
 * r_ij < 2.2 Å hard cap. BOND_SLACK=1.15 (15% slack) absorbs thermal spread;
 * cap 2.2 Å rejects spurious long contacts (Ca–N 2.9 Å) while correctly
 * capturing S–S disulfide 2.04 Å (S 1.02+1.02=2.04, 2.04*1.15=2.35→capped 2.20,
 * still passes). Validated on crambin 1CRN (3 SSBOND): tests/test_topology.js
 * asserts ≥3 S–S bonds recovered and no Ca–N 2.9 Å false bond.
 */
export function buildTopology(atoms) {
  const n = atoms.length;
  const bonds = [];
  const nbond = Array.from({ length: n }, () => []);

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
      if (r < BOND_SLACK * (rA + rB) && r < 2.2) {
        bonds.push([i, j]);
        nbond[i].push(j);
        nbond[j].push(i);
      }
    }
  }

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

  const impropers = [];

  // Ring detection: find 5- and 6-membered rings to add cross-ring distance braces
  const findRings = () => {
    const rings = [];
    const visited = new Set();
    for (let start = 0; start < n; start++) {
      const path = [start];
      const inPath = new Set([start]);
      const dfs = (node, parent) => {
        for (const nbr of nbond[node]) {
          if (nbr === parent) continue;
          if (nbr === start) {
            if (path.length === 5 || path.length === 6) {
              const sortedKey = [...path].sort((a,b)=>a-b).join("-");
              if (!visited.has(sortedKey)) {
                visited.add(sortedKey);
                rings.push([...path]);
              }
            }
          } else if (!inPath.has(nbr) && path.length < 6) {
            inPath.add(nbr);
            path.push(nbr);
            dfs(nbr, node);
            path.pop();
            inPath.delete(nbr);
          }
        }
      };
      dfs(start, -1);
    }
    return rings;
  };

  const detectedRings = findRings();
  for (const ring of detectedRings) {
    const len = ring.length;
    // Cross-ring distance restraints to rigidly maintain planar geometry
    for (let i = 0; i < len; i++) {
      for (let j = i + 2; j < len; j++) {
        if (i === 0 && j === len - 1) continue;
        const idxA = ring[i], idxB = ring[j];
        bonds.push([idxA, idxB]);
      }
    }
  }

  return { bonds, angles, propers, impropers };
}

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

export function buildLists(atoms, topo, coordPairs) {
  const bonds = [];
  for (const [i, j] of topo.bonds) {
    const A = atoms[i], B = atoms[j];
    const r0 = Math.hypot(A.x - B.x, A.y - B.y, A.z - B.z);
    bonds.push(i, j, r0);
  }

  const angles = [];
  for (const [i, j, k] of topo.angles) {
    const th0 = angleAt(atoms, i, j, k);
    angles.push(i, j, k, th0);
  }

  const impropers = [];
  for (const [i, j, k, l] of topo.impropers) {
    const phi0 = improperAngleFlat(atoms, i, j, k, l);
    impropers.push(i, j, k, l, phi0);
  }

  const propers = [];
  for (const [i, j, k, l] of topo.propers) {
    const phi0 = improperAngleFlat(atoms, i, j, k, l);
    propers.push(i, j, k, l, phi0);
  }

  const coord = [];
  for (const [i, j] of coordPairs) {
    const A = atoms[i], B = atoms[j];
    const r0 = Math.hypot(A.x - B.x, A.y - B.y, A.z - B.z);
    coord.push(i, j, r0);
  }

  return { bonds, angles, impropers, propers, coord, holoSprings: new Float64Array(0), nHolo: 0 };
}

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

function improperAngleFlat(atoms, i, j, k, l) {
  const A = atoms[i], B = atoms[j], C = atoms[k], D = atoms[l];
  // vectors ij, jk, kl
  const ij_x = B.x - A.x, ij_y = B.y - A.y, ij_z = B.z - A.z;
  const jk_x = C.x - B.x, jk_y = C.y - B.y, jk_z = C.z - B.z;
  const kl_x = D.x - C.x, kl_y = D.y - C.y, kl_z = D.z - C.z;

  const mx = ij_y * jk_z - ij_z * jk_y;
  const my = ij_z * jk_x - ij_x * jk_z;
  const mz = ij_x * jk_y - ij_y * jk_x;

  const nx = jk_y * kl_z - jk_z * kl_y;
  const ny = jk_z * kl_x - jk_x * kl_z;
  const nz = jk_x * kl_y - jk_y * kl_x;

  const m_len = Math.hypot(mx, my, mz) || 1e-12;
  const n_len = Math.hypot(nx, ny, nz) || 1e-12;
  const jk_len = Math.hypot(jk_x, jk_y, jk_z) || 1e-12;

  const cos_phi = (mx * nx + my * ny + mz * nz) / (m_len * n_len);
  const mxn_x = my * nz - mz * ny;
  const mxn_y = mz * nx - mx * nz;
  const mxn_z = mx * ny - my * nx;
  const sin_phi = (mxn_x * jk_x + mxn_y * jk_y + mxn_z * jk_z) / (m_len * n_len * jk_len);

  return Math.atan2(sin_phi, Math.max(-1, Math.min(1, cos_phi)));
}

/**
 * HeavyForceField — High performance all-atom heavy force field.
 */
export class HeavyForceField {
  constructor(system, par = {}, ligands = []) {
    const atoms = system.atoms;
    this.n = atoms.length;
    this.nProt = atoms.filter((a) => a.isProtein).length;
    this.nLigAtoms = atoms.filter((a) => a.isLigand).length;
    this.nHetero = this.n - this.nProt - this.nLigAtoms;
    this.ligandStart = this.nProt + this.nHetero;
    this.heteroAtoms = atoms.slice(this.nProt, this.ligandStart);
    this.gamma = par.gamma ?? 1.0;
    this.heavy = true;

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
      if (METAL_ELEMENT[el]) {
        this._elem[i] = { ...METAL_ELEMENT[el], q: charges[i] };
      } else {
        const base = LIG_ELEMENT[el] ?? LIG_ELEMENT_DEFAULT;
        this._elem[i] = { ...base, q: charges[i] };
      }
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
   * Build the once-per-topology weak-interaction lists + metal coordination
   * upgrade state (Loop-2 S3).
   * Ring frames from buildRingFrames (protein name sets + GAFF2 types +
   * geometric fallback); cation list from Lys NZ / Arg CZ / HIP / charged
   * ligand N; halogen list from Cl/Br/I with a bonded C; metal upgrade from
   * detectCoordination + classifyGeometry at the construction pose.
   * @param {Array} atoms
   * @param {Array} bonds  flat [i,j,r0,...] list from buildLists
   */
  _buildWeakAndMetalLists(atoms, bonds) {
    const bondPairs = [];
    for (let a = 0; a < bonds.length; a += 3) bondPairs.push([bonds[a], bonds[a + 1]]);
    if (this.weakOn) {
      this._weakRings = buildRingFrames(atoms, bondPairs);
      this._weakRingAtoms = new Set();
      for (const r of this._weakRings) for (const i of r.atomIdx) this._weakRingAtoms.add(i);
      this._weakCations = buildCationList(atoms);
      this._weakHalogens = buildHalogenList(atoms, bondPairs);
      // Acceptor list for halogen bonds: hbond acceptor classification
      // (backbone/sidechain O, S, aromatic N) — same set as the H-bond pass.
      const acc = this._hbClassification.isAcceptor;
      this._weakAcceptors = [];
      for (let i = 0; i < this.n; i++) if (acc[i]) this._weakAcceptors.push(i);
    }
    if (this.metalAngles) {
      // R3 §1e metal upgrade: detect coordination at the native pose, classify
      // the polyhedron, and pin ideal angles. Metals without a detected
      // geometry (coordinationNumber ≤ 1) keep the legacy k=40 springs.
      const elements = atoms.map((a) => a.element ?? "C");
      const metals = [];
      const hasGeometry = new Set();
      for (let i = 0; i < this.n; i++) {
        if (!atoms[i].isMetal) continue;
        const det = detectCoordination(this.ref, i, elements);
        metals.push({ index: i, element: atoms[i].element, donors: det.donorIndices });
        if (det.coordinationNumber >= 2) hasGeometry.add(i);
      }
      this._metalEnforce = { metals, elements, hasGeometry };
      // Remove those metals' radial springs from the legacy coord list so
      // enforceCoordination is the ONLY radial term for them (no double radial).
      if (hasGeometry.size > 0 && this.coord.length > 0) {
        const keep = [];
        for (let a = 0; a < this.coord.length; a += 3) {
          if (hasGeometry.has(this.coord[a])) continue; // metal side of [metal, donor]
          keep.push(this.coord[a], this.coord[a + 1], this.coord[a + 2]);
        }
        this.coord = new Float64Array(keep);
      }
    }
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
  }

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
  }

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
  }

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
  }

  /**
   * Weak-interaction pass (π-stack, cation-π, halogen σ-hole) — opt-in,
   * runs after the LJ/GB grid kernel. Cutoffs: ring-ring 5.5 Å (centroid),
   * cation-ring 6 Å, X···D 4 Å (all inside the kernels' Gaussian tails;
   * pre-screened by centroid distance to skip far pairs cheaply).
   * Exclusions: pairs present in _excluded (bonded/1-3/intra-ligand) and
   * same-ring pairs are skipped. Energies returned per term.
   * @param {ArrayLike<number>} pos
   * @param {Float64Array} f
   * @returns {{pi:number, cpi:number, xb:number}}
   */
  _weakInteractions(pos, f) {
    const rings = this._weakRings;
    const excl = this._excluded;
    let pi = 0, cpi = 0, xb = 0;

    // --- π-stack: ring-ring pairs within 5.5 Å centroid cutoff ---
    for (let a = 0; a < rings.length; a++) {
      const ra = rings[a];
      const a0 = 3 * ra.atomIdx[0];
      const ax = pos[a0], ay = pos[a0 + 1], az = pos[a0 + 2];
      for (let b = a + 1; b < rings.length; b++) {
        const rb = rings[b];
        // cheap prescreen on first atom (within ring diameter of centroid)
        const b0 = 3 * rb.atomIdx[0];
        const dx = pos[b0] - ax, dy = pos[b0 + 1] - ay, dz = pos[b0 + 2] - az;
        if (dx * dx + dy * dy + dz * dz > 121) continue; // 11 Å atom prescreen ≫ 5.5 + 2×2.8 ring radius
        // fused rings share a bond (e.g. Trp 5+6 rings) — those pairs sit in
        // _excluded and must not double-stack; same-molecule NON-bonded rings
        // (Phe–Phe′ stacking) keep the term (R3 §4). Intra-ligand pairs are
        // all-excluded → ligand-internal stacking off (geometry already fixed).
        let skip = false;
        for (const ia of ra.atomIdx) {
          for (const ib of rb.atomIdx) {
            if (ia === ib || excl.has(pairKey(ia, ib))) { skip = true; break; }
          }
          if (skip) break;
        }
        if (skip) continue;
        pi += piStackForces(pos, f, ra, rb);
      }
    }

    // --- cation-π: cation-ring pairs within 6 Å ---
    for (const ci of this._weakCations) {
      const cx = pos[3 * ci], cy = pos[3 * ci + 1], cz = pos[3 * ci + 2];
      for (const ring of rings) {
        // cations inside their own ring (pyridinium N) skip
        if (ring.atomIdx.includes(ci)) continue;
        const i0 = 3 * ring.atomIdx[0];
        const dx = pos[i0] - cx, dy = pos[i0 + 1] - cy, dz = pos[i0 + 2] - cz;
        if (dx * dx + dy * dy + dz * dz > 100) continue; // 10 Å prescreen
        // cation covalently tied to the ring (aniline-type N) — excluded pair
        let skip = false;
        for (const ia of ring.atomIdx) {
          if (excl.has(pairKey(ci, ia))) { skip = true; break; }
        }
        if (skip) continue;
        cpi += cationPiForces(pos, f, ci, ring);
      }
    }

    // --- halogen σ-hole: C–X···D triples, X···D within 4 Å ---
    for (const { x, c } of this._weakHalogens) {
      const xx = pos[3 * x], xy = pos[3 * x + 1], xz = pos[3 * x + 2];
      const el = String(this.atoms[x]?.element ?? "").toUpperCase();
      const epsX = HALOGEN_EPS[el] ?? 1.2;
      for (const d of this._weakAcceptors) {
        if (d === x || d === c) continue;
        const dx = pos[3 * d] - xx, dy = pos[3 * d + 1] - xy, dz = pos[3 * d + 2] - xz;
        const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > 16) continue; // 4 Å
        // bonded X···D / C···D pairs skipped (intra-ligand included)
        if (excl.has(pairKey(x, d)) || excl.has(pairKey(c, d))) continue;
        xb += halogenForces(pos, f, c, x, d, { eps: epsX });
      }
    }

    return { pi, cpi, xb };
  }

  kineticTemp(vel, mass) {
    let ke = 0;
    const m = (mass && mass.length === this.n * 3) ? mass : (this.masses || mass);
    if (m && m.length === this.n * 3) {
      for (let i = 0; i < this.n * 3; i++) ke += m[i] * vel[i] * vel[i];
    } else if (m) {
      for (let i = 0; i < this.n * 3; i++) ke += m[(i / 3) | 0] * vel[i] * vel[i];
    }
    ke *= 0.5 / KCONV;
    return ke / (1.5 * this.n * KB_KCAL);
  }

  rmsd(pos) {
    let s = 0;
    for (let i = 0; i < pos.length; i++) {
      const d = pos[i] - this.ref[i];
      s += d * d;
    }
    return Math.sqrt(s / this.n);
  }
}

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

/**
 * Per-bond harmonic energy/forces with per-entry k (AMBER14 opt-in path).
 * @param {ArrayLike<number>} pos
 * @param {Float64Array} f
 * @param {ArrayLike<number>} list  flat [i,j,r0...]
 * @param {number} stride
 * @param {ArrayLike<number>} kArr  per-bond stiffness
 * @returns {number}
 */
function harmonicFlatPerK(pos, f, list, stride, kArr) {
  let U = 0;
  for (let a = 0, b = 0; a < list.length; a += stride, b++) {
    const i = 3 * list[a], j = 3 * list[a + 1], r0 = list[a + 2];
    const k = kArr[b] ?? 200;
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

/**
 * Per-angle bending with per-entry k (AMBER14 opt-in path).
 * @param {ArrayLike<number>} pos
 * @param {Float64Array} f
 * @param {ArrayLike<number>} list  flat [i,j,k,th0...]
 * @param {number} stride
 * @param {ArrayLike<number>} kArr
 * @returns {number}
 */
function angleFlatPerK(pos, f, list, stride, kArr) {
  let U = 0;
  for (let a = 0, b = 0; a < list.length; a += stride, b++) {
    const k = kArr[b] ?? 40;
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
    const fix = pref * (bxy - c * axy) * ga;
    const fiy = pref * (byy - c * ayy) * ga;
    const fiz = pref * (bzy - c * azy) * ga;
    const fkx = pref * (axy - c * bxy) * gb;
    const fky = pref * (ayy - c * byy) * gb;
    const fkz = pref * (azy - c * bzy) * gb;
    f[i] += fix; f[i + 1] += fiy; f[i + 2] += fiz;
    f[kk] += fkx; f[kk + 1] += fky; f[kk + 2] += fkz;
    f[j] -= fix + fkx; f[j + 1] -= fiy + fky; f[j + 2] -= fiz + fkz;
  }
  return U;
}

/**
 * Physical atomic masses (Da) for heavy atoms — AMBER ff14SB / IUPAC.
 * Exhaustive over METAL_ELEMENT (ZN,FE,MG,CA,CU,MN,NI,CO,NA,K) and
 * LIG_ELEMENT (C,N,O,S,F,CL,BR,I,P) plus common extras B, SE, SI, AL.
 * Fallback for unknown elements warns via console.warn (not silent 14.0)
 * and returns 14.0 (approx N mass) so dynamics remain stable while the
 * caller is alerted — see tests for coverage.
 */
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
    case "B": return 10.81;
    case "SE": return 78.971;
    case "SI": return 28.085;
    case "AL": return 26.982;
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
    default:
      console.warn(`[heavyMass] unknown element "${element}" — fallback 14.0 Da (N mass); add to METAL_ELEMENT/LIG_ELEMENT if needed`);
      return 14.0;
  }
}

function pairKey(i, j) { return i < j ? i * 1e6 + j : j * 1e6 + i; }