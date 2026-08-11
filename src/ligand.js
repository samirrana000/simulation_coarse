/**
 * ligand.js — United-atom internal force field for small-molecule ligands.
 *
 * Model: ligand heavy atoms are explicit beads and hydrogens are implicit;
 * each heavy atom carries a united-atom mass (pdb.unitedAtomMass) that folds
 * the hydrogen mass in. All units follow the project convention: Å, ps,
 * kcal/mol, Da.
 *
 * Harmonic topology terms (energy in kcal/mol):
 *   U_bond = Σ ½ k_b (r − r0)²    k_b = 300 kcal/mol/Å² (general bond)
 *                                     200 kcal/mol/Å² (aromatic ring bond)
 *   U_ang  = Σ ½ k_θ (θ − θ0)²    k_θ = 40 kcal/mol/rad²
 *   U_imp  = Σ ½ k_φ (φ − φ0)²    k_φ = 20 kcal/mol/rad², φ0 = 0 (planarity)
 *
 * Rationale: aromatic rings (C/N 5- and 6-membered rings whose bonds all
 * measure within [1.30, 1.48] Å, the resonance window for aromatic C–C/C–N)
 * get a stiff bond spring pinned near the resonance-averaged distance
 * ([1.35, 1.44] Å), an angle spring that resists ring distortion, and an
 * improper torsion (about the ring plane) that enforces planarity. Non-
 * aromatic bonds simply freeze the measured geometry, like the Cα backbone
 * terms in forcefield.js. Proper dihedrals are intentionally NOT added here —
 * torsional flexibility is left to the ENM-style native contacts in
 * forcefield.js rather than locked by a local potential.
 *
 * Packing mirrors the forcefield.js flat-array convention: bonds as
 * [i, j, r0, ...], angles as [i, j, k, theta0, ...], impropers as
 * [i, j, k, l, phi0, ...] — all with ABSOLUTE particle indices (ligand atom k
 * → global index nProt + k, where nProt counts the preceding protein beads).
 */

import { unitedAtomMass } from "./pdb.js?v=10";

/* Force constants: kcal/mol/Å² (bonds) and kcal/mol/rad² (angles/impropers) */
const KB_AROMATIC = 200;
const KB_GENERAL = 300;
const K_THETA = 40;
const K_PHI = 20;

/* Aromatic bond-length window (Å) used for ring detection */
const AROM_MIN = 1.30;
const AROM_MAX = 1.48;
/* Clamp window for aromatic ring bond equilibrium distance (Å) */
const R0_AROM_MIN = 1.35;
const R0_AROM_MAX = 1.44;

/** Euclidean distance between two parsed atoms (Å). */
function atomDist(a, b) {
  const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Angle (a–b–c) centered at b, in radians. Same geometry as forcefield._angle:
 * vectors from the central atom b to a and c, cosθ clamped to [−1, 1].
 * Returns null for degenerate (zero-length) legs so callers can skip them.
 */
function angleAt(atoms, a, b, c) {
  const ax = atoms[a].x - atoms[b].x, ay = atoms[a].y - atoms[b].y, az = atoms[a].z - atoms[b].z;
  const bx = atoms[c].x - atoms[b].x, by = atoms[c].y - atoms[b].y, bz = atoms[c].z - atoms[b].z;
  const la = Math.hypot(ax, ay, az), lb = Math.hypot(bx, by, bz);
  if (la < 1e-6 || lb < 1e-6) return null; // degenerate (collapsed neighbours)
  let cval = (ax * bx + ay * by + az * bz) / (la * lb);
  cval = Math.min(1, Math.max(-1, cval));
  return Math.acos(cval);
}

/**
 * Locate simple 5- and 6-membered aromatic rings within one molecule.
 * A small DFS starts from each atom; every atom is visited at most once per
 * path and the parent is tracked so the walk never bounces straight back along
 * the edge it came from. A "ring" is a simple path that returns to its start.
 * A ring is marked aromatic when all members are C or N AND every ring bond
 * measures within [AROM_MIN, AROM_MAX] Å.
 * @returns {{atoms: Set<number>, bonds: Set<string>}}
 *   atoms — aromatic atom indices; bonds — sorted "i|j" keys of aromatic bonds
 */
function findAromaticRings(mol) {
  const n = mol.atoms.length;
  const adj = Array.from({ length: n }, () => []);
  for (const [a, b] of mol.bonds) { adj[a].push(b); adj[b].push(a); }

  const atoms = new Set();
  const bonds = new Set();
  const bondKey = (a, b) => (a < b ? `${a}|${b}` : `${b}|${a}`);

  for (let start = 0; start < n; start++) {
    const path = [start];
    const inPath = new Set([start]);

    const dfs = (node, parent) => {
      for (const nbr of adj[node]) {
        if (nbr === parent) continue;
        if (nbr === start) {
          const size = path.length; // ring size = atoms already on the path
          if (size === 5 || size === 6) {
            const allCN = path.every((i) => {
              const el = mol.atoms[i].element.toUpperCase();
              return el === "C" || el === "N";
            });
            let allBondsArom = true;
            for (let t = 0; t < size && allBondsArom; t++) {
              const d = atomDist(mol.atoms[path[t]], mol.atoms[path[(t + 1) % size]]);
              if (d < AROM_MIN || d > AROM_MAX) allBondsArom = false;
            }
            if (allCN && allBondsArom) {
              for (const i of path) atoms.add(i);
              for (let t = 0; t < size; t++) {
                bonds.add(bondKey(path[t], path[(t + 1) % size]));
              }
            }
          }
        } else if (!inPath.has(nbr)) {
          if (path.length < 6) { // only simple rings of size ≤ 6 are of interest
            inPath.add(nbr);
            path.push(nbr);
            dfs(nbr, node);
            path.pop();
            inPath.delete(nbr);
          }
        }
      }
    };
    dfs(start, -1);
  }
  return { atoms, bonds };
}

/**
 * Build the united-atom internal force field for a set of ligand molecules.
 * Atoms are concatenated across ligands in file order; ligand atom k gets the
 * global particle index nProt + k.
 * @param {Array} ligands  output of pdb.parseLigands()
 * @param {number} nProt   number of protein particles preceding the ligands
 * @returns {{
 *   masses: Float64Array, bonds: Float64Array, angles: Float64Array,
 *   impropers: Float64Array, atoms: Array, nLigAtoms: number
 * }}
 *   masses    — united-atom mass (Da) per ligand atom, in concatenation order
 *   bonds     — [i, j, r0, ...] absolute indices (r0 = Å)
 *   angles    — [i, j, k, theta0, ...] absolute indices (theta0 = rad)
 *   impropers — [i, j, k, l, phi0, ...] absolute indices (phi0 = 0 rad)
 *   atoms     — [{idx, element, charge, resName, ring, serial}] one per ligand
 *   nLigAtoms — total ligand atom count
 */
export function buildLigandInternalFF(ligands, nProt) {
  const masses = [];
  const atoms = [];
  const bonds = [];
  const angles = [];
  const impropers = [];

  let k = 0; // running ligand-atom index (global index = nProt + k)
  for (const mol of ligands) {
    const nMol = mol.atoms.length;
    // Cumulative atom offset BEFORE this molecule: bonds/angles/impropers use
    // molecule-LOCAL indices and MUST be shifted by this amount (single-
    // molecule systems have base=0; with several ligands the second molecule's
    // local indices would otherwise collide with the first's particles).
    const base = k;
    const { atoms: arom, bonds: aromBonds } = findAromaticRings(mol);
    const bondKey = (a, b) => (a < b ? `${a}|${b}` : `${b}|${a}`);

    const adj = Array.from({ length: nMol }, () => []);
    for (const [ia, ib] of mol.bonds) { adj[ia].push(ib); adj[ib].push(ia); }

    // ---- masses & atom table ------------------------------------------
    for (let a = 0; a < nMol; a++, k++) {
      const at = mol.atoms[a];
      masses.push(unitedAtomMass(at.element));
      atoms.push({
        idx: nProt + k,
        element: at.element,
        charge: at.charge,
        resName: mol.resName,
        ring: arom.has(a),
        serial: at.serial,
      });
    }

    // ---- bonds ---------------------------------------------------------
    for (const [ia, ib] of mol.bonds) {
      const d = atomDist(mol.atoms[ia], mol.atoms[ib]);
      const inAromaticRing = aromBonds.has(bondKey(ia, ib));
      // Aromatic bonds are pinned inside the resonance window; others freeze
      // the measured length (k differs too, used by the force-field module).
      const r0 = inAromaticRing
        ? Math.min(R0_AROM_MAX, Math.max(R0_AROM_MIN, d))
        : d;
      bonds.push(nProt + base + ia, nProt + base + ib, r0);
    }

    // ---- angles (a–b–c centered at b) ----------------------------------
    for (let b = 0; b < nMol; b++) {
      const nbrs = adj[b];
      for (let p = 0; p < nbrs.length; p++) {
        for (let q = p + 1; q < nbrs.length; q++) {
          const th0 = angleAt(mol.atoms, nbrs[p], b, nbrs[q]);
          if (th0 === null) continue; // skip degenerate geometry
          angles.push(nProt + base + nbrs[p], nProt + base + b, nProt + base + nbrs[q], th0);
        }
      }
    }

    // ---- impropers: ring planarity -------------------------------------
    // For each aromatic ring atom x pick its two ring neighbours p, q and a
    // fourth atom o (any bonded neighbour of x other than p, q). A bare ring
    // atom — e.g. benzene carbon with only 2 heavy neighbours — falls back to
    // the next atom along the ring so the planarity restraint stays defined.
    for (let x = 0; x < nMol; x++) {
      if (!arom.has(x)) continue;
      const p = adj[x].find((nb) => arom.has(nb));
      if (p === undefined) continue;
      const q = adj[x].find((nb) => nb !== p && arom.has(nb));
      if (q === undefined) continue;
      let o = adj[x].find((nb) => nb !== p && nb !== q);
      if (o === undefined) o = adj[p].find((nb) => nb !== x);
      if (o === undefined) continue;
      // Pin at the NATIVE improper angle, not 0: the |atan2| measure equals 0
      // or π for a planar 4-atom set depending on which side of the central
      // bond the fourth atom sits on (cis vs trans — both perfectly flat in
      // fused aromatics like the indole of 1BMA/0QH, where φ0 = 0 produced
      // ≈ 98 kcal/mol of false strain per restraint). Bond/angle terms above
      // already use the native reference for the same reason; only out-of-
      // plane motion is penalized.
      // Build via a scratch coordinate array (atoms are objects here).
      const scr = new Float64Array(nMol * 3);
      for (let a = 0; a < nMol; a++) {
        scr[3 * a] = mol.atoms[a].x; scr[3 * a + 1] = mol.atoms[a].y; scr[3 * a + 2] = mol.atoms[a].z;
      }
      const phi0 = improperAngle(scr, p, x, q, o);
      impropers.push(nProt + base + p, nProt + base + x, nProt + base + q, nProt + base + o, phi0);
    }
  }

  return {
    masses: new Float64Array(masses),
    bonds: new Float64Array(bonds),
    angles: new Float64Array(angles),
    impropers: new Float64Array(impropers),
    atoms,
    nLigAtoms: k,
  };
}

/**
 * Improper torsion angle for (i, j, k, l) with central atom j, in radians in
 * [0, π]. Convention: the improper measures the dihedral between the planes
 * (i–j–k) and (j–k–l) about the central j–k bond, so a planar arrangement of
 * the four atoms gives φ ≈ 0 and moving the "l" atom out of the plane
 * increases φ.
 *
 *   b1 = r_j − r_i,   b2 = r_k − r_j,   b3 = r_l − r_k
 *   n1 = b1 × b2,     n2 = b2 × b3
 *   φ  = atan2( (n1 × n2) · b̂2 ,  n1 · n2 )
 *
 * Normals are normalized and the dot-product argument is clamped to [−1, 1] so
 * a degenerate / collinear geometry can never produce NaN; the absolute value
 * maps the signed atan2 result onto [0, π].
 * @param {Float64Array} pos  flat 3n coordinates
 * @returns {number} improper angle in radians, [0, π]
 */
export function improperAngle(pos, i, j, k, l) {
  const i3 = 3 * i, j3 = 3 * j, k3 = 3 * k, l3 = 3 * l;
  const b1x = pos[j3] - pos[i3], b1y = pos[j3 + 1] - pos[i3 + 1], b1z = pos[j3 + 2] - pos[i3 + 2];
  const b2x = pos[k3] - pos[j3], b2y = pos[k3 + 1] - pos[j3 + 1], b2z = pos[k3 + 2] - pos[j3 + 2];
  const b3x = pos[l3] - pos[k3], b3y = pos[l3 + 1] - pos[k3 + 1], b3z = pos[l3 + 2] - pos[k3 + 2];

  // n1 = b1 × b2,  n2 = b2 × b3  (both ⊥ the central j–k bond)
  const n1x = b1y * b2z - b1z * b2y;
  const n1y = b1z * b2x - b1x * b2z;
  const n1z = b1x * b2y - b1y * b2x;
  const n2x = b2y * b3z - b2z * b3y;
  const n2y = b2z * b3x - b2x * b3z;
  const n2z = b2x * b3y - b2y * b3x;

  const ln1 = Math.hypot(n1x, n1y, n1z);
  const ln2 = Math.hypot(n2x, n2y, n2z);
  const lb2 = Math.hypot(b2x, b2y, b2z);
  if (ln1 < 1e-12 || ln2 < 1e-12 || lb2 < 1e-12) return 0; // degenerate

  const u1x = n1x / ln1, u1y = n1y / ln1, u1z = n1z / ln1;
  const u2x = n2x / ln2, u2y = n2y / ln2, u2z = n2z / ln2;
  const ux = b2x / lb2, uy = b2y / lb2, uz = b2z / lb2;

  // (n̂1 × n̂2)·b̂2  and  n̂1·n̂2 (clamped to avoid NaN)
  const cx = u1y * u2z - u1z * u2y;
  const cy = u1z * u2x - u1x * u2z;
  const cz = u1x * u2y - u1y * u2x;
  const y = cx * ux + cy * uy + cz * uz;
  const x = Math.min(1, Math.max(-1, u1x * u2x + u1y * u2y + u1z * u2z));

  return Math.abs(Math.atan2(y, x));
}
