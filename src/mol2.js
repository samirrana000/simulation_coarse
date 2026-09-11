/**
 * mol2.js — Tripos MOL2 ligand parser.
 *
 * Shape contract: parseMol2() returns an array of molecule objects:
 *   { resName, chain, atoms: [{x,y,z, element, charge, serial}], bonds: [[i,j]] }
 *
 * Fidelity (I82 — SDF/MOL2 full fidelity):
 *  - bond order is preserved for topology only — MOL2 bond types ("1","2","3","ar","am", etc.)
 *    are read but stored as unweighted edges [[i,j]]; topology is preserved while bond order
 *    is ignored for the CG force field (all bonds share the same harmonic force constant).
 *  - triple bond not collapsed — a bond entry with type "3" (e.g. C≡C) is kept as a single
 *    topological edge [[i,j]] without collapsing or expanding to multiple bonds; same for
 *    "2" and "ar" — no special-casing, connectivity is bond-order-agnostic.
 *  - Du/H dropped united-atom — dummy atoms (SYBYL type "Du") and hydrogens (element "H")
 *    are dropped (united-atom model); mol2Element("Du") returns null and ATOM filter
 *    `element === "H"` skips protons so heavy-atom indices stay compact.
 *
 * Input validation / warnings (A06):
 *  - Malformed ATOM/BOND lines are skipped with `warnings` and `console.warn`.
 *    The returned array has a `.warnings` property (string[]). Use
 *    `countWarnings(result.warnings)` to count skipped lines. The `continue`
 *    on malformed lines is lenient but now warns.
 */

const TWO_LETTER_ELEMENTS = new Set([
  "CL", "BR", "ZN", "FE", "MG", "CA", "CU", "MN", "NI", "CO", "NA", "SE", "SI", "AL", "LI",
]);

/**
 * Element symbol from a Tripos MOL2 atom type ("C.ar"→"C", "O.co2"→"O",
 * "Cl"→"CL", "Br"→"BR", "Du"→null for dummy atoms).
 * @param {string} atomType  raw SYBYL atom type (field 6 of an ATOM line)
 * @returns {string|null} uppercase element symbol, or null for dummies
 */
export function mol2Element(atomType) {
  const base = String(atomType || "").split(".")[0].trim();
  if (!base || base.toLowerCase() === "du") return null;   // dummy atom — skip
  const first = base.charAt(0).toUpperCase();
  const second = base.charAt(1).toLowerCase();
  const two = (first + second).toUpperCase();
  if (TWO_LETTER_ELEMENTS.has(two)) return two;
  return first;
}

/**
 * Helper to count validation warnings.
 * @param {string[]|object} warnings
 * @returns {number}
 */
export function countWarnings(warnings) {
  if (Array.isArray(warnings)) return warnings.length;
  if (warnings && Array.isArray(warnings.warnings)) return warnings.warnings.length;
  return 0;
}

/**
 * Parse a Tripos MOL2 file into ligand molecules:
 *   { resName, chain, atoms: [{x,y,z, element, charge, serial}], bonds: [[i,j]] }
 * Warnings: malformed lines are collected in `molecules.warnings` and warned.
 * @param {string} mol2Text
 * @returns {Array} molecule list with .warnings
 */
export function parseMol2(mol2Text) {
  const molecules = [];
  const warnings = [];
  let section = null;      // current @<TRIPOS>… section
  let mol = null;          // molecule under construction
  let idMap = null;        // MOL2 atom id -> heavy-atom index (per molecule)
  let droppedIds = null;   // ids of dummy/H atoms dropped (for silent BOND skip)

  const finishMol = () => {
    if (mol && mol.atoms.length >= 1) molecules.push(mol);
    mol = null; idMap = null; droppedIds = null;
  };

  for (const raw of mol2Text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line.startsWith("#")) continue;

    if (line.startsWith("@<TRIPOS>")) {
      section = line.slice(9).trim().toUpperCase();
      if (section === "MOLECULE") finishMol();
      continue;
    }
    if (!section) continue;

    const f = line.split(/\s+/);

    if (section === "MOLECULE") {
      if (!mol) {
        mol = { resName: (f[0] || "LIG").slice(0, 3).toUpperCase(), chain: "L", atoms: [], bonds: [] };
        idMap = new Map();
        droppedIds = new Set();
      }
    } else if (section === "ATOM") {
      if (!mol || f.length < 6) {
        if (mol && f.length < 6) {
          const msg = `parseMol2: malformed ATOM line skipped (fields <6) line="${line}"`;
          warnings.push(msg); console.warn(msg);
        }
        continue;
      }
      const id = parseInt(f[0], 10);
      const x = parseFloat(f[2]), y = parseFloat(f[3]), z = parseFloat(f[4]);
      if (Number.isNaN(id) || Number.isNaN(x + y + z)) {
        const msg = `parseMol2: malformed ATOM coords skipped id=${f[0]} line="${line}"`;
        warnings.push(msg); console.warn(msg);
        continue;
      }
      const element = mol2Element(f[5]);
      if (!element || element === "H") {
        if (droppedIds) droppedIds.add(id);
        continue;   // drop dummy / hydrogen
      }
      const charge = f.length >= 9 ? parseFloat(f[8]) : NaN;
      idMap.set(id, mol.atoms.length);
      mol.atoms.push({
        x, y, z,
        element,
        charge: Number.isNaN(charge) ? 0 : charge,
        serial: id,
        atomName: f[1] || element,
      });
    } else if (section === "BOND") {
      if (!mol || !idMap || f.length < 4) {
        if (mol && f.length < 4) {
          const msg = `parseMol2: malformed BOND line skipped line="${line}"`;
          warnings.push(msg); console.warn(msg);
        }
        continue;
      }
      // I82: bond order is preserved for topology only — f[3] (e.g. "1","2","3","ar","am")
      // is intentionally ignored; triple bond not collapsed, stored as [[i,j]].
      const aId = parseInt(f[1], 10), bId = parseInt(f[2], 10);
      // const bondOrder = f[3]; // kept for topology only, not used for force constant
      const i = idMap.get(aId);
      const j = idMap.get(bId);
      if (i === undefined || j === undefined || i === j) {
        // Dropped H/dummy produce missing ids — silent skip (not malformed)
        const droppedA = droppedIds && droppedIds.has(aId);
        const droppedB = droppedIds && droppedIds.has(bId);
        if ((i === undefined && !droppedA) || (j === undefined && !droppedB)) {
          // Truly unknown heavy-atom id — warn (malformed MOL2)
          const msg = `parseMol2: BOND references unknown atom id line="${line}"`;
          warnings.push(msg); console.warn(msg);
        } else if (i === j && i !== undefined) {
          const msg = `parseMol2: self-bond skipped line="${line}"`;
          warnings.push(msg); console.warn(msg);
        }
        continue;
      }
      const dup = mol.bonds.some(([p, q]) => (p === i && q === j) || (p === j && q === i));
      if (!dup) mol.bonds.push([i, j]);
    }
  }
  finishMol();
  molecules.warnings = warnings;
  if (warnings.length) console.warn(`[parseMol2] ${warnings.length} warning(s) total`);
  return molecules;
}
