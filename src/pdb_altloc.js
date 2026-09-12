/**
 * pdb_altloc.js — Shared altloc / duplicate-atom cleaner (Stage-4).
 *
 * Single choke-point helper for all PDB parsing paths (`parseCa`,
 * `parseLigands` in src/pdb.js and `parseHeavy` in src/heavy.js).
 *
 * Rule (documented, default-ON):
 *  - Dedup key: (chain, resSeq, iCode, atom name). For Cα-only parsing the
 *    atom name is always CA so the key reduces to (chain, resSeq, iCode).
 *  - Occupancy from PDB columns 55–60 (`line.slice(54, 60)`); blank or
 *    unparseable → 1.0 (assume fully occupied; matches PDB convention that
 *    blank occupancy means 1.0).
 *  - altLoc from column 17 (`line.charAt(16)`); blank → '' (single conformer).
 *  - Winner: highest occupancy wins. Zero-occupancy copies are dropped when a
 *    non-zero copy exists (and a zero-occupancy existing entry is replaced by
 *    any non-zero incoming copy).
 *  - Tie (within 1e-9): prefer altLoc 'A' over any non-'A' (covers B-then-A
 *    file order); otherwise keep the first record encountered. Exact
 *    duplicates (same altLoc, same coordinates, repeated serial) therefore
 *    collapse to a single atom (first wins).
 *  - Bit-identical on clean files: 4W52 altloc pairs are all 0.50/0.50 ties
 *    with 'A' first, so "highest occupancy, 'A' on tie, else first" keeps the
 *    same 'A' record the old first-wins code kept (verified: 4W52 164 Cα /
 *    1308 heavy, 1CRN 46 Cα / 327 heavy unchanged).
 *
 * Zero dependencies; pure functions. Units: occupancy unitless, Å untouched.
 */

/** Tolerance for occupancy tie detection. */
export const ALTLOC_OCC_TIE = 1e-9;

/**
 * Parse the altLoc character (PDB column 17).
 * @param {string} line  raw PDB record line
 * @returns {string} altLoc code ('' for blank / single conformer)
 */
export function parseAltLoc(line) {
  const c = line.charAt(16) || " ";
  const t = c.trim();
  return t || "";
}

/**
 * Parse occupancy (PDB columns 55–60). Blank/unparseable → 1.0.
 * @param {string} line  raw PDB record line
 * @returns {number} occupancy (0 if the field explicitly reads 0)
 */
export function parseOccupancy(line) {
  const raw = line.slice(54, 60).trim();
  if (!raw) return 1.0;
  const v = parseFloat(raw);
  return Number.isFinite(v) ? v : 1.0;
}

/**
 * Decide whether an incoming duplicate record should replace the stored one.
 * @param {{occupancy:number, altLoc:string}} existing  stored winner so far
 * @param {{occupancy:number, altLoc:string}} incoming  new duplicate record
 * @returns {boolean} true → replace stored with incoming; false → drop incoming
 */
export function shouldReplaceAltloc(existing, incoming) {
  const eo = existing.occupancy, io = incoming.occupancy;
  // Zero-occupancy rule: non-zero always beats zero.
  if (io === 0 && eo > 0) return false;
  if (eo === 0 && io > 0) return true;
  // Strictly higher occupancy wins.
  if (io > eo + ALTLOC_OCC_TIE) return true;
  if (eo > io + ALTLOC_OCC_TIE) return false;
  // Tie: prefer 'A' over non-'A' (order-independent); else keep first.
  if (incoming.altLoc === existing.altLoc) return false;
  if (incoming.altLoc === "A" && existing.altLoc !== "A") return true;
  if (existing.altLoc === "A" && incoming.altLoc !== "A") return false;
  return false;
}
