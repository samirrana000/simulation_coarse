/**
 * pdb.js — Minimal PDB parser & loader with selection support.
 *
 * Responsibilities:
 *  - fetchPdb(id): download a structure from the public RCSB / PDBe APIs
 *    (tries https first candidate then falls back; plain fetch, no backend libs).
 *  - parseCa(pdbText): extract C-alpha (plus glycine CA) ATOM records from all
 *    standard models, keeping one CA per (chain, resSeq, iCode) key. Works for
 *    apo single-chain structures and multi-subunit complexes alike.
 *  - selectSystem(parsed, {chains, resFrom, resTo}): filter beads by chain IDs
 *    and inclusive residue-number range; splits each chain into contiguous
 *    segments (gaps > 4.8 Å between consecutive Cα ⇒ chain break; typical
 *    Cα–Cα distance is 3.8 Å).
 *
 * The Tripos MOL2 ligand parser now lives in mol2.js (item 5 modularization)
 * and is re-exported here so existing `import { parseMol2 } from "./pdb.js"`
 * call-sites are unchanged.
 *
 * No external units conversion happens here: coordinates stay in Ångström.
 */

import { parseMol2, mol2Element } from "./mol2.js?v=10";
export { parseMol2, mol2Element };

const RCSB_URL = (id) => `https://files.rcsb.org/download/${id}.pdb`;
const PDBE_URL = (id) => `https://www.ebi.ac.uk/pdbe/entry-files/download/pdb${id.toLowerCase()}.ent`;

/** Fetch a PDB file by 4-letter ID. Tries RCSB, then PDBe as a fallback. */
export async function fetchPdb(id) {
  const clean = id.trim();
  if (!/^[0-9A-Za-z]{4}$/.test(clean)) {
    throw new Error(`"${id}" is not a valid 4-character PDB ID.`);
  }
  for (const url of [RCSB_URL(clean.toUpperCase()), PDBE_URL(clean)]) {
    try {
      const res = await fetch(url);
      if (res.ok) {
        const text = await res.text();
        if (text.includes("ATOM")) return text;
      }
    } catch (_) {
      /* network error — try the fallback mirror */
    }
  }
  throw new Error(`Could not download PDB ${clean} from RCSB or PDBe (network/CORS?).`);
}

/**
 * Parse Cα atoms out of raw PDB text (first MODEL only, to avoid NMR bundles).
 * Handles insertion codes via a unique key per residue.
 * @returns {{beads: Array, chains: string[], nAtoms: number}}
 *   beads: [{x,y,z, chain, resSeq, resName, serial, bfac}]
 */
export function parseCa(pdbText) {
  const beads = [];
  const seen = new Set();
  const chains = new Set();
  let nAtoms = 0;
  let inModel = true; // becomes false after ENDMDL (first atomic model only)

  for (const line of pdbText.split(/\r?\n/)) {
    const rec = line.slice(0, 6);
    if (rec === "ENDMDL") {            // first model of an NMR ensemble ends here
      inModel = false;
      break;
    }
    if (!inModel) break;
    if (rec !== "ATOM  ") continue;

    // PDB fixed-column layout (keep leading spaces!)
    const atomName = line.slice(12, 16).trim();
    if (atomName !== "CA") continue;   // coarse-grain: one bead per residue at Cα

    const resName = line.slice(17, 20).trim();
    const chain = (line.charAt(21) || " ").trim() || "_"; // blank chain ID -> "_"
    const resSeq = parseInt(line.slice(22, 26), 10);
    const iCode = line.charAt(26).trim();
    const x = parseFloat(line.slice(30, 38));
    const y = parseFloat(line.slice(38, 46));
    const z = parseFloat(line.slice(46, 54));
    const bfac = parseFloat(line.slice(60, 66)); // PDB temperature factor (columns 61–66)
    if (Number.isNaN(resSeq) || Number.isNaN(x + y + z)) continue;

    const key = `${chain}|${resSeq}|${iCode}`;
    if (seen.has(key)) continue;       // alt-loc duplicates
    seen.add(key);

    beads.push({ x, y, z, chain, resSeq, resName, serial: nAtoms + 1, bfac: Number.isNaN(bfac) ? 0 : bfac });
    chains.add(chain);
    nAtoms++;
  }

  if (beads.length === 0) {
    throw new Error("No Cα ATOM records found — is this a valid (protein) PDB file?");
  }
  return { beads, chains: [...chains].sort(), nAtoms };
}

/**
 * Apply user selection criteria to a parsed structure.
 * @param {object} parsed   output of parseCa()
 * @param {object} sel      { chains: string[]|null, resFrom: number|null, resTo: number|null }
 * @returns {{beads, segments, nChains}}
 *   beads    — flat selected bead array (order preserved)
 *   segments — Array<[start, endExclusive]> contiguous stretches within one
 *              chain used to define peptide bonds (gaps are NOT bonded).
 */
export function selectSystem(parsed, { chains = null, resFrom = null, resTo = null } = {}) {
  const chainSet = chains && chains.length ? new Set(chains.map((c) => c.trim())) : null;

  const beads = parsed.beads.filter(
    (b) =>
      (!chainSet || chainSet.has(b.chain)) &&
      (resFrom === null || b.resSeq >= resFrom) &&
      (resTo === null || b.resSeq <= resTo)
  );
  if (beads.length < 3) {
    throw new Error("Selection must contain at least 3 Cα beads (need bonded + angle terms).");
  }

  // Build contiguous segments per chain: peptide backbone is continuous only
  // within a segment; a Cα–Cα gap > 4.8 Å marks a chain break / missing loop.
  const segments = [];
  let segStart = 0;
  for (let i = 1; i < beads.length; i++) {
    const prev = beads[i - 1];
    const cur = beads[i];
    const dx = cur.x - prev.x, dy = cur.y - prev.y, dz = cur.z - prev.z;
    const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
    const sameChain = cur.chain === prev.chain;
    if (!sameChain || d > 4.8) {
      segments.push([segStart, i]);
      segStart = i;
    }
  }
  segments.push([segStart, beads.length]);

  const nChains = new Set(beads.map((b) => b.chain)).size;
  return { beads, segments, nChains };
}

/**
 * Summarize a parsed structure for the UI:
 * apo vs. complex, chains and their residue ranges.
 */
export function summarizeStructure(parsed) {
  const byChain = new Map();
  for (const b of parsed.beads) {
    if (!byChain.has(b.chain)) byChain.set(b.chain, []);
    byChain.get(b.chain).push(b.resSeq);
  }
  const lines = [`${parsed.nAtoms} residues (Cα), ${parsed.chains.length} chain(s) — ` +
    (parsed.chains.length === 1 ? "apo / single-chain structure" : "multi-subunit complex")];
  for (const c of parsed.chains) {
    const r = byChain.get(c);
    lines.push(`  Chain ${c}: res ${Math.min(...r)}–${Math.max(...r)} (${r.length} Cα)`);
  }
  return lines.join("\n");
}

const SOLVENT = new Set(["HOH", "WAT", "SOL", "DOD"]);

/**
 * Derive an element symbol from a PDB atom name when the element column is blank.
 * @param {string} atomName
 * @returns {string} uppercase element symbol
 */
function elementFromName(atomName) {
  const first = atomName.trim().charAt(0).toUpperCase();
  if (first === "C") return "C";
  if (first === "N") return "N";
  if (first === "O") return "O";
  if (first === "S") return "S";
  if (first === "P") return "P";
  if (first === "F") return "F";
  if (first === "I") return "I";
  if (first === "B") return "BR";
  if (first === "H") return "H";   // protonated PDB HETATM: keep as H so the
                                   // united-atom filters can DROP it instead of
                                   // misreading it as carbon and placing a fake
                                   // C at a 1.0 Å bond length (explosion risk)
  return "C";
}

/**
 * Parse HETATM records into ligand molecules with bonds.
 * Molecules are grouped by chain|resSeq|resName in file order; bonds come from
 * CONECT records, falling back to a 1.9 Å heavy-atom distance cutoff when none.
 * @param {string} pdbText
 * @returns {Array} [{ resName, chain, atoms: [{x,y,z, element, charge, serial}], bonds: [[i,j]] }]
 */
export function parseLigands(pdbText) {
  const molecules = [];
  const byKey = new Map();

  for (const line of pdbText.split(/\r?\n/)) {
    const rec = line.slice(0, 6);
    if (rec !== "HETATM") continue;

    const atomName = line.slice(12, 16).trim();
    const resName = line.slice(17, 20).trim();
    if (SOLVENT.has(resName)) continue;
    const chain = (line.charAt(21) || " ").trim() || "_";
    const resSeq = parseInt(line.slice(22, 26), 10);
    const x = parseFloat(line.slice(30, 38));
    const y = parseFloat(line.slice(38, 46));
    const z = parseFloat(line.slice(46, 54));
    if (Number.isNaN(resSeq) || Number.isNaN(x + y + z)) continue;

    let element = line.slice(76, 78).trim().toUpperCase();
    if (!element) element = elementFromName(atomName);
    if (element === "H") continue;   // united-atom model: drop hydrogens
    const serial = parseInt(line.slice(6, 11), 10);

    const key = `${chain}|${resSeq}|${resName}`;
    let mol = byKey.get(key);
    if (!mol) {
      mol = { resName, chain, atoms: [], bonds: [], _serialToIdx: new Map() };
      byKey.set(key, mol);
      molecules.push(mol);
    }
    mol._serialToIdx.set(serial, mol.atoms.length);
    mol.atoms.push({ x, y, z, element, charge: 0, serial });
  }

  for (const line of pdbText.split(/\r?\n/)) {
    const rec = line.slice(0, 6);
    if (rec !== "CONECT") continue;

    const a = parseInt(line.slice(6, 11), 10);
    for (const [s, e] of [[11, 16], [16, 21], [21, 26], [26, 31]]) {
      const b = parseInt(line.slice(s, e), 10);
      if (Number.isNaN(a) || Number.isNaN(b) || a === b) continue;
      const mol = molecules.find((m) => m._serialToIdx.has(a) && m._serialToIdx.has(b));
      if (!mol) continue;
      const i = mol._serialToIdx.get(a);
      const j = mol._serialToIdx.get(b);
      if (i === j) continue;
      const existing = mol.bonds.some(([p, q]) =>
        (p === i && q === j) || (p === j && q === i));
      if (!existing) mol.bonds.push([i, j]);
    }
  }

  return molecules
    .filter((m) => m.atoms.length >= 2)
    .map((m) => {
      delete m._serialToIdx;
      // CONECT records in real PDBs are frequently INCOMPLETE (e.g. buffer
      // molecules like HEPES) — an omitted covalent bond leaves two heavy atoms
      // at ~1.4 Å with no repulsion exclusion, which explodes the CG force
      // field. Gap-fill: any heavy-atom pair closer than 1.8 Å that is not yet
      // directly bonded is treated as covalent (1.8 Å ≈ the longest C–N/C–O
      // bond; longer contacts are non-covalent and stay repulsive).
      const CUT = 1.8;
      for (let i = 0; i < m.atoms.length; i++) {
        for (let j = i + 1; j < m.atoms.length; j++) {
          const A = m.atoms[i], B = m.atoms[j];
          const dx = A.x - B.x, dy = A.y - B.y, dz = A.z - B.z;
          if (Math.sqrt(dx * dx + dy * dy + dz * dz) > CUT) continue;
          const existing = m.bonds.some(([p, q]) =>
            (p === i && q === j) || (p === j && q === i));
          if (!existing) m.bonds.push([i, j]);
        }
      }
      return m;
    });
}

/**
 * Map a PDB element symbol to a coarse united-atom mass in Da.
 * @param {string} element
 * @returns {number} mass in Da
 */
export function unitedAtomMass(element) {
  const e = element.toUpperCase();
  switch (e) {
    case "C": return 12.01;
    case "N": return 14.01;
    case "O": return 16.00;
    case "S": return 32.06;
    case "P": return 30.97;
    case "F": return 19.00;
    case "CL": return 35.45;
    case "BR": return 79.90;
    case "I": return 126.90;
    default: return 14.0;
  }
}
