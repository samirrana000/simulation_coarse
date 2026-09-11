/**
 * mmcif.js — PDBx/mmCIF stub (I81).
 *
 * Primary structure I/O in this app is PDB via src/pdb.js:parseCa.
 * mmCIF is reserved for future native support; this stub detects the
 * mmCIF `data_` block header and either delegates to the PDB parser
 * when ATOM-compatible records exist or throws a clear guidance error.
 *
 * See docs/MMCIF.md for status and migration notes.
 */

import { parseCa } from "./pdb.js?v=10";

/**
 * Parse a PDBx/mmCIF string into the same {beads, chains, nAtoms} shape
 * as parseCa, or throw a user-facing error guiding to PDB.
 *
 * Detection: mmCIF files start with a `data_<entryId>` block (case-insensitive
 * check for `data_` at top-of-file). If detected, we attempt a best-effort
 * ATOM fallback; otherwise we throw "mmCIF not yet supported, use PDB".
 *
 * @param {string} text raw file text
 * @returns {{beads: Array, chains: string[], nAtoms: number}}
 * @throws {Error} if mmCIF is requested without ATOM fallback
 */
export function parseMMCIF(text) {
  if (typeof text !== "string" || !text.trim()) {
    throw new Error("parseMMCIF: empty input");
  }
  const hasDataBlock = /^\s*data_/im.test(text) || text.includes("data_");
  if (!hasDataBlock) {
    // Not an mmCIF file — delegate to PDB parser so callers can try both
    return parseCa(text);
  }

  // I81 — mmCIF detection hit: file contains `data_` (PDBx/mmCIF block)
  // Future: parse _atom_site.* loop properly (mmCIF spec: _atom_site.group_PDB,
  // _atom_site.label_atom_id, _atom_site.Cartn_x/y/z, _atom_site.B_iso_or_equiv).
  // For now, if the mmCIF text also carries ATOM-like lines (some exporters do),
  // delegate to parseCa so coarse-grained Cα path keeps working.
  if (text.includes("ATOM") && text.includes(" CA ")) {
    try {
      return parseCa(text);
    } catch (_) {
      // fall through to guidance error
    }
  }

  // Preferred error message required by I81 spec (grep-visible)
  throw new Error("mmCIF not yet supported, use PDB — convert to PDB via pdb_extract or gemmi: gemmi convert input.cif output.pdb");
}

/**
 * Quick detector for caller UI (optional).
 * @param {string} text
 * @returns {boolean} true if text looks like mmCIF (data_ block)
 */
export function isMMCIF(text) {
  return typeof text === "string" && /^\s*data_/im.test(text);
}
