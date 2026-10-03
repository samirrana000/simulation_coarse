/**
 * mmcif.js — PDBx/mmCIF: NOT SUPPORTED (absent, not pending).
 *
 * Primary structure I/O in this app is PDB via src/pdb.js:parseCa.
 * mmCIF is absent by scope, not queued for implementation: ROADMAP.md §1
 * lists the browser-v1 hard no's, and a native `_atom_site` loop parser
 * would be a format project of its own. This module therefore does exactly
 * two things and nothing more: detect a mmCIF `data_` block so the caller
 * can say so precisely, and delegate to the PDB parser when the input
 * already carries ATOM-compatible records (some exporters emit both).
 *
 * See docs/MMCIF.md for the honest-scope entry and what to use instead.
 */

import { parseCa } from "./pdb.js";
import { ignore } from "./errors.js";

/**
 * Parse a PDBx/mmCIF string into the same {beads, chains, nAtoms} shape
 * as parseCa, or throw a user-facing error guiding to PDB.
 *
 * Detection: mmCIF files start with a `data_<entryId>` block (case-insensitive
 * check for `data_` at top-of-file). If detected, we attempt a best-effort
 * ATOM fallback; otherwise we throw an honest-scope error naming the converter.
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

  // mmCIF detection hit: file contains `data_` (PDBx/mmCIF block).
  // A native _atom_site.* loop parser is ABSENT, not pending — see ROADMAP.md §1
  // for the project's hard-no list and docs/MMCIF.md for the honest-scope entry.
  // Do not read "not supported" as a queued feature. If the mmCIF text also
  // carries ATOM-like lines (some exporters do), delegate to parseCa so the
  // coarse-grained Cα path keeps working.
  if (text.includes("ATOM") && text.includes(" CA ")) {
    try {
      return parseCa(text);
    } catch (e) {
      // The file IS mmCIF, so surfacing parseCa's "no Cα found in this PDB"
      // would be a misleading message; fall through to the honest-scope error
      // below. Counted, because a silent fallback here is how a user ends up
      // believing a structure was parsed when it was not.
      ignore(e, "parseCa fallback@parseMMCIF", "input is mmCIF (has a data_ block), so parseCa's PDB-shaped error would misattribute the cause; the honest-scope throw below replaces it");
    }
  }

  // Preferred error message (the string is grep-visible and asserted by
  // tests/test_error_surfacing.js): ABSENT capability, stated as absent.
  throw new Error("mmCIF is not supported by this simulator — PDB is the ingest format. Convert with pdb_extract or gemmi: gemmi convert input.cif output.pdb");
}

/**
 * Quick detector for caller UI (optional).
 * @param {string} text
 * @returns {boolean} true if text looks like mmCIF (data_ block)
 */
export function isMMCIF(text) {
  return typeof text === "string" && /^\s*data_/im.test(text);
}
