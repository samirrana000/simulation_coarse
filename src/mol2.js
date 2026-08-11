/**
 * mol2.js — Tripos MOL2 ligand parser (item 5 modularization: moved verbatim
 * from pdb.js; re-exported there so existing imports are unchanged).
 *
 * Shape contract: parseMol2() returns the SAME molecule shape as
 * pdb.js#parseLigands() —
 *   { resName, chain, atoms: [{x,y,z, element, charge, serial}], bonds: [[i,j]] }
 */

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
  // two-letter halogens/metals written lowercase or mixed case in SYBYL files
  const two = (first + second).toUpperCase();
  if (two === "CL" || two === "BR") return two;
  return first;
}

/**
 * Parse a Tripos MOL2 file into ligand molecules (same shape as parseLigands):
 *   { resName, chain, atoms: [{x,y,z, element, charge, serial}], bonds: [[i,j]] }
 *
 * One file may hold several @<TRIPOS>MOLECULE blocks — each becomes a separate
 * molecule. Hydrogens are DROPPED (united-atom model: ligand.js folds H mass
 * into the heavy atom, and LIG_ELEMENT in forcefield.js has no H entry), and
 * so are "Du" dummy atoms. Atom/bond IDs are remapped per molecule to the
 * contiguous heavy-atom indices used in `bonds`. The MOL2 charge column is
 * parsed into atom.charge for information only — the force field derives its
 * charges from the element table (same as parseLigands). Bond *types* (1, 2,
 * ar, am) are ignored; connectivity alone drives bond/angle/ring detection in
 * ligand.js.
 *
 * @param {string} mol2Text
 * @returns {Array} molecule list (empty if the file has no usable molecules)
 */
export function parseMol2(mol2Text) {
  const molecules = [];
  let section = null;      // current @<TRIPOS>… section
  let mol = null;          // molecule under construction
  let idMap = null;        // MOL2 atom id -> heavy-atom index (per molecule)

  const finishMol = () => {
    if (mol && mol.atoms.length >= 2) molecules.push(mol);
    mol = null; idMap = null;
  };

  for (const raw of mol2Text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line.startsWith("#")) continue;

    if (line.startsWith("@<TRIPOS>")) {
      section = line.slice(9).trim().toUpperCase();
      if (section === "MOLECULE") finishMol();   // previous block (if any) is complete
      continue;
    }
    if (!section) continue;

    const f = line.split(/\s+/);

    if (section === "MOLECULE") {
      // first non-blank line of the block = molecule name
      if (!mol) {
        mol = { resName: (f[0] || "LIG").slice(0, 3).toUpperCase(), chain: "L", atoms: [], bonds: [] };
        idMap = new Map();
      }
      // remaining MOL2 header lines (counts, type, …) are ignored
    } else if (section === "ATOM") {
      if (!mol || f.length < 6) continue;
      // id name x y z atom_type [subst_id subst_name charge …]
      const id = parseInt(f[0], 10);
      const x = parseFloat(f[2]), y = parseFloat(f[3]), z = parseFloat(f[4]);
      if (Number.isNaN(id) || Number.isNaN(x + y + z)) continue;
      const element = mol2Element(f[5]);
      if (!element || element === "H") continue;   // dummy / hydrogen: drop
      const charge = f.length >= 9 ? parseFloat(f[8]) : NaN;
      idMap.set(id, mol.atoms.length);
      mol.atoms.push({ x, y, z, element, charge: Number.isNaN(charge) ? 0 : charge, serial: id });
    } else if (section === "BOND") {
      if (!mol || !idMap || f.length < 4) continue;
      // id a1 a2 type — connectivity only is used
      const i = idMap.get(parseInt(f[1], 10));
      const j = idMap.get(parseInt(f[2], 10));
      if (i === undefined || j === undefined || i === j) continue;
      const dup = mol.bonds.some(([p, q]) => (p === i && q === j) || (p === j && q === i));
      if (!dup) mol.bonds.push([i, j]);
    }
    // UNIFIED_ATOM / CRYSIN / FF_PBC etc. sections are ignored
  }
  finishMol();
  return molecules;
}
