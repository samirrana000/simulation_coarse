/**
 * controllers/structure-input.js — everything that gets a structure INTO the
 * app (PDB id fetch, 4W52 one-click sample, file drop, MOL2 ligand file).
 *
 * PUBLIC API
 *   initStructureInput({buildSystem, updateGuide})
 *                        Wire #fetchBtn (+Enter), the `[data-ex]` preset links,
 *                        #sampleBtn, #fileInput and #mol2File.
 *   loadStructure(text, sourceLabel)
 *                        Validate + parse `text` into `state.parsed` (Cα) and
 *                        `state.parsedHeavy`, refresh the summary caption and
 *                        chain into buildSystem(). Throws a classified
 *                        InputError (see ../input_errors.js) on bad input; the
 *                        callers below render it with formatInputError.
 *
 * NEEDS (imports): ui + state from ../ui.js; fetchPdb + parseCa + parseLigands
 * + parseMol2 + summarizeStructure from ../pdb.js; parseHeavy from ../heavy.js;
 * validatePdbText + classifyInputError + formatInputError from
 * ../input_errors.js; updateMol2PlaceButton from ../ligand-panel.js.
 * `buildSystem` and `updateGuide` are passed in by the composition root rather
 * than imported, so this module does not need to know the build/guide module
 * order.
 *
 * ZERO PHYSICS HERE: parsing and validation only; the force field is built by
 * ./system-build.js.
 */

/* ------------------------------------------------------------------ */
/*  Structure loading                                                  */
/* ------------------------------------------------------------------ */
import { ui, state } from "../ui.js";
import { fetchPdb, parseCa, parseMol2, summarizeStructure } from "../pdb.js";
import { parseHeavy } from "../heavy.js";
import { classifyInputError, formatInputError, validatePdbText } from "../input_errors.js";
import { updateMol2PlaceButton } from "../ligand-panel.js";
import { ignore } from "../errors.js";

/** Build-system entry point, injected by initStructureInput. */
let _buildSystem = null;
/** Guide repaint hook, injected by initStructureInput. */
let _updateGuide = null;

/**
 * Parse a PDB text into the app and rebuild.
 * @param {string} text raw PDB text
 * @param {string} sourceLabel caption prefix (e.g. "PDB 4W52")
 * @param {object} [source] provenance for the results record
 *   (`{origin, pdbId, fileName}`; see state.inputSource). Additive and
 *   optional: omitting it leaves state.inputSource null and the results record
 *   reports origin "unknown" rather than guessing.
 */
export async function loadStructure(text, sourceLabel, source = null) {
  // FP2: pre-validate raw text so empty/garbage inputs get an actionable
  // failure class (what + next click) instead of a raw parser dump.
  const preErr = validatePdbText(text);
  if (preErr) throw preErr;
  state.pdbText = text;
  let caErr = null, heavyErr = null;
  try {
    state.parsed = parseCa(text);
  } catch (e) {
    state.parsed = null;
    caErr = e;
  }
  try {
    state.parsedHeavy = parseHeavy(text);
  } catch (e) {
    state.parsedHeavy = null;
    heavyErr = e;
  }

  if (!state.parsed && !state.parsedHeavy) {
    // FP2: map the raw parser failure to an actionable failure class
    // (technical detail kept as secondary suffix by formatInputError).
    throw classifyInputError(caErr || heavyErr || new Error("No usable atoms found in PDB file."));
  }

  state.libraryLigand = null;
  state.heteroOverrides = {};
  // Where the bytes came from, for the exported results record
  // (src/results-record.js). The CONTENT HASH is the identity; this is the
  // human-readable route, and it is null when the caller did not say.
  state.inputSource = source && typeof source === "object"
    ? {
      origin: String(source.origin ?? "unknown"),
      pdbId: String(source.pdbId ?? ""),
      fileName: String(source.fileName ?? ""),
      label: String(sourceLabel ?? ""),
    }
    : null;
  if (ui.structSummary) {
    if (state.parsed) {
      ui.structSummary.textContent = `${sourceLabel}\n` + summarizeStructure(state.parsed);
    } else {
      ui.structSummary.textContent = `${sourceLabel}\n${state.parsedHeavy.atoms.length} heavy atoms found (Heavy Mode active)`;
    }
  }
  if (ui.chainsInput) ui.chainsInput.value = "";
  if (ui.resFrom) ui.resFrom.value = "";
  if (ui.resTo) ui.resTo.value = "";
  _buildSystem();
}

/** Render a classified input failure into the two reused caption surfaces. */
function reportInputError(err) {
  // FP2: actionable copy primary, raw detail secondary (no raw dump).
  if (ui.structSummary) ui.structSummary.textContent = formatInputError(err);
  if (ui.hud) ui.hud.textContent = formatInputError(err);
}

/** Wire every structure-acquisition control. Call once, after the build wiring. */
export function initStructureInput({ buildSystem, updateGuide }) {
  _buildSystem = buildSystem;
  _updateGuide = updateGuide;

  if (ui.fetchBtn && ui.pdbId) {
    ui.fetchBtn.addEventListener("click", async () => {
      const id = ui.pdbId.value.trim();
      if (!id) return;
      if (ui.structSummary) ui.structSummary.textContent = `Fetching ${id}…`;
      try {
        await loadStructure(await fetchPdb(id), `PDB ${id.toUpperCase()}`,
          { origin: "pdb-id-fetch", pdbId: id, fileName: `${id.toLowerCase()}.pdb` });
      } catch (err) {
        reportInputError(err);
      }
    });

    ui.pdbId.addEventListener("keydown", (e) => {
      if (e.key === "Enter") ui.fetchBtn.click();
    });
  }

  document.querySelectorAll("[data-ex]").forEach((a) =>
    a.addEventListener("click", (e) => {
      e.preventDefault();
      if (ui.pdbId && ui.fetchBtn) {
        ui.pdbId.value = a.dataset.ex;
        ui.fetchBtn.click();
      }
    })
  );

  // FP1 — one-click 4W52 sample (first-run UX). Reuses the existing fetchPdb
  // path (local ./4w52.pdb first, offline OK; RCSB/PDBe fallback) via the
  // same preset mechanism as the data-ex links above. User-initiated only —
  // no autoload, so the default view is unchanged for returning users.
  if (ui.sampleBtn) {
    ui.sampleBtn.addEventListener("click", () => {
      if (ui.pdbId) ui.pdbId.value = "4W52";
      if (ui.fetchBtn) ui.fetchBtn.click();
      try { _updateGuide(); } catch (e) { ignore(e, "_updateGuide@load", "checklist DOM absent headless"); }
    });
  }

  if (ui.fileInput) {
    ui.fileInput.addEventListener("change", async () => {
      const f = ui.fileInput.files[0];
      if (!f) return;
      try {
        await loadStructure(await f.text(), f.name,
          { origin: "file-upload", pdbId: "", fileName: f.name });
      } catch (err) {
        reportInputError(err);
      }
    });
  }

  /* ---------------------------------------------------------------- */
  /*  MOL2 ligand input (overrides PDB HETATM ligands)                 */
  /* ---------------------------------------------------------------- */
  if (ui.mol2File) {
    ui.mol2File.addEventListener("change", async () => {
      const f = ui.mol2File.files[0];
      if (!f) {
        state.mol2Ligands = null;
        state.mol2Fn = null;
        if (ui.mol2Info) ui.mol2Info.style.display = "none";
        updateMol2PlaceButton();
        _buildSystem();
        return;
      }
      try {
        const mols = parseMol2(await f.text());
        // FP2: empty MOL2 ⇒ actionable LIGAND_PARSE_FAIL (not a raw dump).
        if (!mols.length) throw classifyInputError(new Error(`${f.name}: no usable molecules in MOL2 file`), { stage: "mol2" });
        state.mol2Ligands = mols;
        state.mol2Fn = f.name;
        const nAtoms = mols.reduce((s, m) => s + m.atoms.length, 0);
        const nBonds = mols.reduce((s, m) => s + m.bonds.length, 0);
        if (ui.mol2Info) {
          ui.mol2Info.textContent =
            `${f.name}: ${mols.length} molecule(s), ${nAtoms} heavy atoms, ${nBonds} bonds — active (overrides HETATM)`;
          ui.mol2Info.style.display = "block";
        }
        updateMol2PlaceButton();
        _buildSystem();
      } catch (err) {
        state.mol2Ligands = null;
        state.mol2Fn = null;
        if (ui.mol2Info) {
          // FP2: actionable copy primary, raw detail secondary (no raw dump).
          ui.mol2Info.textContent = formatInputError(classifyInputError(err, { stage: "mol2" }));
          ui.mol2Info.style.display = "block";
        }
        updateMol2PlaceButton();
        _buildSystem();
      }
    });
  }
}