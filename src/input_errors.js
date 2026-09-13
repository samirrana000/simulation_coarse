/**
 * input_errors.js — Central input validator: failure class → actionable copy.
 *
 * FP2 (finish-product loop, input robustness): bad PDB text, missing/empty
 * ligand, unparsable MOL2, clash-blast placement, and oversized systems used
 * to surface raw `err.message` dumps or blank states. This module maps every
 * load/build/place failure to one stable failure class whose message states
 * WHAT happened plus the EXACT next click, keeping the raw technical detail
 * as a secondary `(detail: …)` suffix (the single-line summary/caption
 * elements have no room for a collapsed <details>, so the parenthetical is
 * the collapsed-equivalent — primary copy stays actionable).
 *
 * Failure classes (all zero-dep, pure, headless-safe):
 *   EMPTY_PDB | NO_ATOM | NO_CA | LIGAND_PARSE_FAIL | NO_POCKET |
 *   CLASH_HIGH | SYSTEM_TOO_LARGE | SELECTION_EMPTY | GENERIC
 *
 * Thresholds (documented, defaults unchanged — validators only, no physics):
 *   CLASH_RESIDUAL_THRESHOLD = 2.0 (unitless ratio 1/minRatio; a converged
 *     relaxClash pose has residual ≤ ~1.18, so > 2.0 means the worst contact
 *     sits below 50% of its vdW sum — severe steric overlap).
 *   MAX_HEAVY_ATOMS = 20000, MAX_CA_BEADS = 5000 (interactive state limit;
 *     4W52 reference is 1308 heavy / 164 Cα, so the limit is ~15×/30× headroom
 *     and never fires on the bundled systems).
 *
 * Mapped errors ARE Error instances (`.code`, `.nextStep`, `.technical`
 * attached), so existing `catch (err)` paths can `throw`/`format` them with
 * no control-flow changes. `classifyInputError` passes already-mapped errors
 * through untouched (idempotent) and never throws itself.
 */

import { parseCa, parseLigands } from "./pdb.js?v=10";
import { parseMol2 } from "./mol2.js?v=10";

/** Stable failure-class codes (frozen contract — tests assert membership). */
export const INPUT_ERROR_CODES = [
  "EMPTY_PDB",
  "NO_ATOM",
  "NO_CA",
  "LIGAND_PARSE_FAIL",
  "NO_POCKET",
  "CLASH_HIGH",
  "SYSTEM_TOO_LARGE",
  "SELECTION_EMPTY",
  "GENERIC",
];

/**
 * Clash residual over this unitless threshold ⇒ CLASH_HIGH.
 * (residual = 1/minRatio from relaxClash; converged ⇒ ≤ ~1.18.)
 */
export const CLASH_RESIDUAL_THRESHOLD = 2.0;

/** Interactive state limits (atom counts, not bytes — see header). */
export const MAX_HEAVY_ATOMS = 20000;
export const MAX_CA_BEADS = 5000;
/** Alias for the combined state limit (heavy-atom scale). */
export const MAX_SYSTEM_ATOMS = MAX_HEAVY_ATOMS;

/**
 * Actionable copy per failure class: what happened + the exact next click.
 * Every `message` contains an imperative next-step verb ("Click") so the
 * primary surface is never a dead end. `nextStep` repeats the click alone
 * for callers that render it separately.
 */
export const INPUT_ERROR_COPY = {
  EMPTY_PDB: {
    message:
      "Empty input — no PDB text to parse. Click “Load 4W52 sample (1 click)” or Fetch a PDB ID to load a structure.",
    nextStep: "Click “Load 4W52 sample (1 click)” in Structure.",
  },
  NO_ATOM: {
    message:
      "No ATOM/HETATM records found — this file is not a PDB structure. Click “Load 4W52 sample (1 click)” or drop a valid .pdb file.",
    nextStep: "Click “Load 4W52 sample (1 click)” or drop a valid .pdb file.",
  },
  NO_CA: {
    message:
      "No Cα backbone found — no protein Cα records to coarse-grain. Click “Load 4W52 sample (1 click)” for a protein structure, or switch Model & Selection → Simulation Model → All-atom heavy mode, then click Build System.",
    nextStep: "Click “Load 4W52 sample (1 click)”, then click Build System.",
  },
  LIGAND_PARSE_FAIL: {
    message:
      "Ligand file could not be parsed — no usable molecules found. Click Ligand MOL2 file and choose a valid .mol2, or click Place in Pocket (Auto) with a library ligand.",
    nextStep: "Click Ligand MOL2 file and choose a valid .mol2 file.",
  },
  NO_POCKET: {
    message:
      "No binding pocket found — the protein selection is empty. Click Build System with at least 3 Cα beads selected, then click Place in Pocket (Auto).",
    nextStep: "Click Build System, then click Place in Pocket (Auto).",
  },
  CLASH_HIGH: {
    message:
      "Placed pose still clashes (residual over threshold) — steric overlap remains. Click Place in Pocket (Auto) to retry clash-free, or click Random Surface for a fresh start.",
    nextStep: "Click Place in Pocket (Auto) to retry clash-free.",
  },
  SYSTEM_TOO_LARGE: {
    message:
      `System too large — exceeds the interactive state limit (${MAX_HEAVY_ATOMS} heavy / ${MAX_CA_BEADS} Cα atoms). Click Model & Selection → restrict Chains or Residue range, then click Build System.`,
    nextStep: "Click Model & Selection → restrict Chains or Residue range, then click Build System.",
  },
  SELECTION_EMPTY: {
    message:
      "Selection too small — need at least 3 Cα beads for bonded terms. Click Model & Selection → clear the Chains filter or widen Residue range, then click Build System.",
    nextStep: "Click Model & Selection → clear the Chains filter, then click Build System.",
  },
  GENERIC: {
    message:
      "Something failed to load (detail below). Click “Load 4W52 sample (1 click)” to restore a known-good system, then retry.",
    nextStep: "Click “Load 4W52 sample (1 click)” to restore a known-good system.",
  },
};

/**
 * Build a mapped input error (an Error with .code/.nextStep/.technical).
 * Never throws — unknown codes fall back to GENERIC.
 * @param {string} code one of INPUT_ERROR_CODES
 * @param {string|Error|null} [technical] raw underlying detail (kept secondary)
 * @returns {Error} actionable error (message = what + next click)
 */
export function createInputError(code, technical = null) {
  const c = INPUT_ERROR_CODES.includes(code) ? code : "GENERIC";
  const copy = INPUT_ERROR_COPY[c];
  const tech =
    technical instanceof Error ? technical.message : technical == null ? "" : String(technical);
  const err = new Error(copy.message);
  err.code = c;
  err.nextStep = copy.nextStep;
  err.technical = tech || copy.message;
  return err;
}

/**
 * Map any load/build/place failure to its failure class.
 * Idempotent (already-mapped errors pass through) and total: never throws,
 * never returns null — unrecognized input yields GENERIC.
 * @param {Error|string|null|undefined} err raw failure
 * @param {object} [context] optional disambiguator ({ stage })
 * @returns {Error} mapped actionable error (.code/.nextStep/.technical)
 */
export function classifyInputError(err, context = {}) {
  try {
    if (err && typeof err === "object" && INPUT_ERROR_CODES.includes(err.code)) {
      if (err instanceof Error) return err;
      return createInputError(err.code, err.technical ?? err.message);
    }
    const raw = err instanceof Error ? err.message : err == null ? "" : String(err);
    const t = raw.toLowerCase();
    const stage = (context && context.stage) || "";

    if (!raw || !raw.trim() || /empty (input|pdb)|no pdb text|no text to parse/.test(t)) {
      return createInputError("EMPTY_PDB", raw);
    }
    if (/residual|clash|steric|degenerate pose/.test(t)) {
      return createInputError("CLASH_HIGH", raw);
    }
    if (/too large|exceeds.*limit|state limit|maximum.*atoms/.test(t)) {
      return createInputError("SYSTEM_TOO_LARGE", raw);
    }
    if (/at least 3|selection (must contain|too small)|need.*beads/.test(t)) {
      return createInputError("SELECTION_EMPTY", raw);
    }
    if (
      stage === "mol2" ||
      /mol2|no usable molecules?|ligand.*(pars|fail|missing|empty)|no .*molecule/.test(t)
    ) {
      return createInputError("LIGAND_PARSE_FAIL", raw);
    }
    if (/pocket|nprot.*0|protein selection is empty|no protein/.test(t)) {
      return createInputError("NO_POCKET", raw);
    }
    if (/no cα|no ca\b|backbone.*missing|coarse-grain/.test(t)) {
      return createInputError("NO_CA", raw);
    }
    if (
      /no (usable atoms|atom|heavy atoms)|no atom\/hetatm|not a (valid |protein )?pdb|not a pdb structure|no .*records found/.test(
        t
      )
    ) {
      return createInputError("NO_ATOM", raw);
    }
    if (/valid 4-character pdb|could not download/.test(t)) {
      return createInputError("NO_ATOM", raw);
    }
    return createInputError("GENERIC", raw);
  } catch (_) {
    return createInputError("GENERIC", null);
  }
}

/**
 * Render a mapped (or raw) failure for the existing single-line
 * summary/caption surfaces: actionable copy primary, raw detail secondary
 * in a `(detail: …)` suffix (the collapsed-equivalent — these elements
 * cannot host a <details> disclosure).
 * @param {Error|string|null} err mapped or raw failure
 * @returns {string} "⚠ <what + next click> (detail: <technical>)"
 */
export function formatInputError(err) {
  let mapped;
  try {
    mapped = classifyInputError(err);
  } catch (_) {
    mapped = createInputError("GENERIC", null);
  }
  const tech = mapped.technical && mapped.technical !== mapped.message ? ` (detail: ${mapped.technical})` : "";
  return `⚠ ${mapped.message}${tech}`;
}

/**
 * Pre-validate raw PDB text before parsing (pure, never throws).
 * @param {string} text raw file/paste/drop text
 * @returns {Error|null} mapped EMPTY_PDB / NO_ATOM error, or null when usable
 */
export function validatePdbText(text) {
  try {
    if (typeof text !== "string" || !text.trim()) {
      return createInputError("EMPTY_PDB", typeof text === "string" ? "empty string" : `typeof ${typeof text}`);
    }
    if (!/^(ATOM  |HETATM)/m.test(text)) {
      return createInputError("NO_ATOM", "no ATOM/HETATM records");
    }
    return null;
  } catch (_) {
    return createInputError("GENERIC", null);
  }
}

/**
 * Total guard for parser calls: never throws on garbage — returns a
 * structured { ok, data, error } result instead.
 * @param {Function} fn parser function
 * @param {*} text parser input
 * @param {string} emptyCode failure class when fn throws on this input
 * @returns {{ok:boolean, data:*, error:Error|null}}
 */
export function guardParse(fn, text, emptyCode = "GENERIC") {
  try {
    const data = fn(text);
    return { ok: true, data, error: null };
  } catch (e) {
    return { ok: false, data: null, error: classifyInputError(e, { stage: "load" }) || createInputError(emptyCode, String(e)) };
  }
}

/**
 * parseCa that never throws uncaught on garbage inputs.
 * @param {*} text raw PDB text (any type)
 * @returns {{ok:boolean, data:object|null, error:Error|null}}
 */
export function safeParseCa(text) {
  const pre = validatePdbText(text);
  if (pre) return { ok: false, data: null, error: pre };
  return guardParse(parseCa, text, "NO_CA");
}

/**
 * parseMol2 that never throws uncaught on garbage inputs. An empty molecule
 * list is a LIGAND_PARSE_FAIL here (a MOL2 file must carry ≥1 molecule —
 * unlike PDB HETATM, where zero ligands is a valid apo protein).
 * @param {*} text raw MOL2 text (any type)
 * @returns {{ok:boolean, data:Array|null, error:Error|null}}
 */
export function safeParseMol2(text) {
  try {
    if (typeof text !== "string" || !text.trim()) {
      return { ok: false, data: null, error: createInputError("LIGAND_PARSE_FAIL", "empty MOL2 text") };
    }
    const mols = parseMol2(text);
    if (!mols || mols.length === 0) {
      return { ok: false, data: null, error: createInputError("LIGAND_PARSE_FAIL", "no usable molecules in MOL2 file") };
    }
    return { ok: true, data: mols, error: null };
  } catch (e) {
    return { ok: false, data: null, error: classifyInputError(e, { stage: "mol2" }) };
  }
}

/**
 * parseLigands that never throws uncaught on garbage inputs. Zero molecules
 * is OK:true (valid apo protein — absence of HETATM is not a ligand error;
 * the LIGAND_PARSE_FAIL class fires at MOL2-load / place time instead).
 * @param {*} text raw PDB text (any type)
 * @returns {{ok:boolean, data:Array|null, error:Error|null}}
 */
export function safeParseLigands(text) {
  try {
    if (typeof text !== "string" || !text.trim()) {
      return { ok: false, data: null, error: createInputError("EMPTY_PDB", "empty string") };
    }
    return { ok: true, data: parseLigands(text), error: null };
  } catch (e) {
    return { ok: false, data: null, error: classifyInputError(e, { stage: "load" }) };
  }
}

/**
 * System-size check against the interactive state limit (pure, never throws).
 * @param {object} counts { nHeavy, nCa } atom counts (missing ⇒ 0)
 * @returns {Error|null} mapped SYSTEM_TOO_LARGE error, or null when within limit
 */
export function checkSystemSize(counts = {}) {
  try {
    const nHeavy = Number(counts.nHeavy) || 0;
    const nCa = Number(counts.nCa) || 0;
    if (nHeavy > MAX_HEAVY_ATOMS || nCa > MAX_CA_BEADS) {
      return createInputError(
        "SYSTEM_TOO_LARGE",
        `nHeavy=${nHeavy} (limit ${MAX_HEAVY_ATOMS}), nCa=${nCa} (limit ${MAX_CA_BEADS})`
      );
    }
    return null;
  } catch (_) {
    return null;
  }
}

/**
 * Placement-outcome check: clash score over threshold or degenerate pose
 * (pure, never throws). A null/degenerate `placed` maps to CLASH_HIGH with
 * the degenerate detail preserved as technical.
 * @param {object|null} placed relaxClash/placeLigand result ({ converged, residualClash })
 * @returns {Error|null} mapped CLASH_HIGH error, or null when clash-free
 */
export function checkPlacement(placed) {
  try {
    if (!placed || !Number.isFinite(placed.residualClash)) {
      return createInputError("CLASH_HIGH", "degenerate pose (non-finite residual)");
    }
    if (placed.converged !== true || placed.residualClash > CLASH_RESIDUAL_THRESHOLD) {
      return createInputError(
        "CLASH_HIGH",
        `residual ${Number(placed.residualClash).toFixed(2)} over threshold ${CLASH_RESIDUAL_THRESHOLD} (converged=${placed.converged === true})`
      );
    }
    return null;
  } catch (_) {
    return createInputError("CLASH_HIGH", null);
  }
}
