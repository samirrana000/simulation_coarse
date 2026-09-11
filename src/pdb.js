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
 *
 * Input validation / warnings (A06):
 *  - parseCa, parseLigands, and parseHeavy (see heavy.js) are tolerant parsers:
 *    malformed ATOM/HETATM lines are skipped with a `warnings` entry and a
 *    `console.warn`. Callers may inspect `result.warnings` (array of strings)
 *    and `countWarnings(result.warnings)` for a count. The `continue` on
 *    malformed lines is intentional (lenient PDB spec), but now always warns.
 *  - Helper `countWarnings(warnings)` returns warnings.length (0 if none).
 */

import { parseMol2, mol2Element } from "./mol2.js?v=10";
export { parseMol2, mol2Element };

/**
 * Count validation warnings from a parser result.
 * @param {string[]|object} warnings  array of warning strings or result with .warnings
 * @returns {number}
 */
export function countWarnings(warnings) {
  if (Array.isArray(warnings)) return warnings.length;
  if (warnings && Array.isArray(warnings.warnings)) return warnings.warnings.length;
  return 0;
}

const RCSB_URL = (id) => `https://files.rcsb.org/download/${id}.pdb`;
const PDBE_URL = (id) => `https://www.ebi.ac.uk/pdbe/entry-files/download/pdb${id.toLowerCase()}.ent`;

/**
 * Fetch a PDB file by 4-letter ID. Tries local cache, RCSB, then PDBe.
 *
 * REST fetch cache (I88):
 *  - Remote fetches send `Cache-Control: max-age=86400` (24 h) so the browser
 *    may serve a fresh disk-cache entry without revalidation for a day.
 *  - In addition, fetch uses browser cache + ETag fallback — the browser's
 *    HTTP cache (disk/memory) + conditional ETag/Last-Modified revalidation
 *    handles 304 Not Modified from RCSB/PDBe, so repeated 4W52 fetches do not
 *    re-download when the server says the file is unchanged.
 *  - Local `./*.pdb` and `./data/*.pdb` try first and are always cache-hit;
 *    remote only on cache-miss for unknown IDs.
 * Docs: see fetch header below and note in docs/SCORER.md (I88 cross-ref).
 */
export async function fetchPdb(id) {
  const clean = id.trim();
  if (!/^[0-9A-Za-z]{4}$/.test(clean)) {
    throw new Error(`"${id}" is not a valid 4-character PDB ID.`);
  }
  const urls = [
    `./${clean.toLowerCase()}.pdb`,
    `./data/${clean.toLowerCase()}.pdb`,
    RCSB_URL(clean.toUpperCase()),
    PDBE_URL(clean),
  ];
  for (const url of urls) {
    try {
      // I88: Cache-Control: max-age=86400 via fetch headers;
      // fetch uses browser cache + ETag fallback for revalidation.
      const isRemote = url.startsWith("http");
      const res = await fetch(url, isRemote ? { headers: { "Cache-Control": "max-age=86400" }, cache: "default" } : undefined);
      if (res.ok) {
        const text = await res.text();
        if (text.includes("ATOM")) return text;
      }
    } catch (_) {
      /* try next mirror */
    }
  }
  throw new Error(`Could not download PDB ${clean} from local presets, RCSB, or PDBe.`);
}

/**
 * Parse Cα atoms out of raw PDB text (first MODEL only, to avoid NMR bundles).
 * Handles insertion codes via a unique key per residue.
 * Warnings: malformed lines are collected in `warnings` and emitted via
 * `console.warn`; see module header for warnings behavior. Caller may use
 * `countWarnings(result.warnings)` to count skipped lines.
 * Note on disulfides & PTMs (C28): this CG parser keeps only Cα beads, so
 * Cys SG atoms are intentionally dropped. Disulfide S–S <2.2 Å detection is
 * handled in the heavy-atom path (see src/heavy.js:248 buildTopology with
 * COVALENT_RADIUS S=1.02 Å, BOND_SLACK=1.15, hard cap 2.2 Å; SG–SG 2.04 Å
 * passes with 0.16 Å margin, validated in tests/test_topology.js on 1crn).
 * If SG coordinates are present in the PDB, heavy mode recovers all three
 * crambin disulfides; the CG ENM instead captures the Cα–Cα restraint
 * (Rc=10 Å) and does not need an explicit SG warning. Native PTMs (e.g.
 * phosphoserine) are likewise coalesced to Cα in CG mode and explicit in
 * heavy mode. Should a future CG variant retain SG, a check
 * `if (resName==="CYS" && sgDist<2.2) warnings.push("Cys SG–SG <2.2 Å")`
 * would be added here.
 * @returns {{beads: Array, chains: string[], nAtoms: number, warnings: string[]}}
 *   beads: [{x,y,z, chain, resSeq, resName, serial, bfac}]
 *   warnings: array of human-readable warning strings (may be empty)
 */
export function parseCa(pdbText) {
  const beads = [];
  const seen = new Set();
  const chains = new Set();
  const warnings = [];
  let nAtoms = 0;
  let inModel = true; // becomes false after ENDMDL (first atomic model only)

  for (const line of pdbText.split(/\r?\n/)) {
    const rec = line.slice(0, 6);
    if (rec === "ENDMDL") {            // first model of an NMR ensemble ends here
      inModel = false;
      break;
    }
    const isAtom = rec.startsWith("ATOM") || rec.startsWith("HETATM");
    if (!isAtom) continue;

    // PDB fixed-column layout (keep leading spaces!)
    const atomName = line.slice(12, 16).trim().toUpperCase();
    if (atomName !== "CA") continue;   // coarse-grain: one bead per residue at Cα

    const resName = line.slice(17, 20).trim();
    const chain = (line.charAt(21) || " ").trim() || "_"; // blank chain ID -> "_"
    const resSeq = parseInt(line.slice(22, 26), 10);
    const iCode = line.charAt(26).trim();
    const x = parseFloat(line.slice(30, 38));
    const y = parseFloat(line.slice(38, 46));
    const z = parseFloat(line.slice(46, 54));
    const bfac = parseFloat(line.slice(60, 66)); // PDB temperature factor (columns 61–66)
    if (Number.isNaN(resSeq) || Number.isNaN(x + y + z)) {
      const msg = `parseCa: malformed CA line skipped (resSeq=${resSeq} x=${x} y=${y} z=${z}) line="${line.slice(0, 66).trim()}"`;
      warnings.push(msg);
      console.warn(msg);
      continue;
    }

    const key = `${chain}|${resSeq}|${iCode}`;
    if (seen.has(key)) {
      const msg = `parseCa: duplicate residue ${key} skipped (altLoc)`;
      warnings.push(msg);
      console.warn(msg);
      continue;       // alt-loc duplicates
    }
    seen.add(key);

    beads.push({ x, y, z, chain, resSeq, resName, serial: nAtoms + 1, bfac: Number.isNaN(bfac) ? 0 : bfac });
    chains.add(chain);
    nAtoms++;
  }

  if (beads.length === 0) {
    throw new Error("No Cα ATOM records found — is this a valid (protein) PDB file?");
  }
  if (warnings.length) console.warn(`[parseCa] ${warnings.length} warning(s) total`);
  return { beads, chains: [...chains].sort(), nAtoms, warnings };
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
 * Warnings: malformed HETATM lines are skipped with console.warn and collected
 * in `molecules.warnings` (see module header). Use `countWarnings()` to count.
 * @param {string} pdbText
 * @returns {Array} [{ resName, chain, atoms: [{x,y,z, element, charge, serial}], bonds: [[i,j]] }] with .warnings
 */
export function parseLigands(pdbText) {
  const molecules = [];
  const byKey = new Map();
  const warnings = [];

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
    if (Number.isNaN(resSeq) || Number.isNaN(x + y + z)) {
      const msg = `parseLigands: malformed HETATM skipped (resSeq=${resSeq} x=${x}) line="${line.slice(0, 54).trim()}"`;
      warnings.push(msg);
      console.warn(msg);
      continue;
    }

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

  const filtered = molecules
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
  // Attach warnings array for callers that care (lenient parser — non-fatal)
  filtered.warnings = warnings;
  if (warnings.length) console.warn(`[parseLigands] ${warnings.length} warning(s) total`);
  return filtered;
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
