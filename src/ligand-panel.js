/**
 * ligand-panel.js — "Ligand Placement & Library" panel.
 *
 * Provides:
 *   1. Searchable built-in ligand library (benzene, indole, caffeine, etc.)
 *   2. Click-to-place on viewer with collision relaxation
 *   3. Auto-placement into detected binding pockets (clash-free)
 *   4. Random surface encounter placement
 *   5. Seamless placement of loaded MOL2 ligands in both Cα and Heavy mode.
 */

import { LIGAND_LIBRARY } from "./ligandLib.js?v=10";
import { parseLibraryLigand, placeLigand, findPocketCenter } from "./placement.js?v=10";
import { LIG_ELEMENT, LIG_ELEMENT_DEFAULT } from "./ff-params.js?v=10";
import { classifyInputError, formatInputError, createInputError, checkPlacement } from "./input_errors.js?v=10";
import { ui, state, viewer } from "./ui.js?v=10";

let _buildSystem = () => {};

/* ------------------------------------------------------------------ */
/*  Library selector                                                   */
/* ------------------------------------------------------------------ */

function renderLibraryList() {
  if (!ui.ligSelect || !ui.ligFilter) return;
  const q = ui.ligFilter.value.trim().toLowerCase();
  const opts = LIGAND_LIBRARY
    .filter((e) => !q || e.name.toLowerCase().includes(q) || e.id.toLowerCase().includes(q))
    .map((e) => `<option value="${e.id}">${e.name} (${e.id})</option>`)
    .join("");
  ui.ligSelect.innerHTML = opts || `<option value="">— no match —</option>`;
  if (![...ui.ligSelect.options].some((o) => o.value === ui.ligSelect.value)) {
    ui.ligSelect.value = LIGAND_LIBRARY[0] ? LIGAND_LIBRARY[0].id : "";
  }
}

function selectedEntry() {
  if (!ui.ligSelect) return LIGAND_LIBRARY[0];
  return LIGAND_LIBRARY.find((e) => e.id === ui.ligSelect.value) || LIGAND_LIBRARY[0];
}

/* ------------------------------------------------------------------ */
/*  Placement interaction                                               */
/* ------------------------------------------------------------------ */

let pickActive = false;
let pickTarget = "library";

function syncPlaceButtons() {
  if (ui.placeBtn) {
    ui.placeBtn.textContent = pickActive && pickTarget === "library" ? "Click the viewer…" : "Place on viewer";
  }
  if (ui.placeMol2Btn) {
    ui.placeMol2Btn.textContent = pickActive && pickTarget === "mol2" ? "Click the viewer…" : "Place loaded MOL2 ligand";
  }
  if (ui.cancelPlace) ui.cancelPlace.disabled = !pickActive;
  if (viewer?.canvas) viewer.canvas.style.cursor = pickActive ? "crosshair" : "";
}

function setPickMode(on, target = "library") {
  pickActive = on;
  pickTarget = target;
  syncPlaceButtons();
}

export function updateMol2PlaceButton() {
  if (!ui.placeMol2Btn) return;
  const ok = !!(state.mol2Ligands && state.mol2Ligands.length === 1);
  ui.placeMol2Btn.disabled = !ok;
  syncPlaceButtons();
}

function applyPlacedPose(mol, name, placed) {
  // FP2: degenerate pose ⇒ actionable CLASH_HIGH (not a bare warning).
  if (!placed || !Number.isFinite(placed.residualClash)) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("CLASH_HIGH", "degenerate pose (non-finite residual)"));
    setPickMode(false);
    return;
  }

  const nLib = mol.atoms.length;
  const ligOff = state.ff.n - nLib;
  for (let a = 0; a < nLib; a++) {
    state.integ.pos[3 * (ligOff + a)] = placed.pos[3 * a];
    state.integ.pos[3 * (ligOff + a) + 1] = placed.pos[3 * a + 1];
    state.integ.pos[3 * (ligOff + a) + 2] = placed.pos[3 * a + 2];
    state.ff.ref[3 * (ligOff + a)] = placed.pos[3 * a];
    state.ff.ref[3 * (ligOff + a) + 1] = placed.pos[3 * a + 1];
    state.ff.ref[3 * (ligOff + a) + 2] = placed.pos[3 * a + 2];
  }

  state.ff.holoOn = false;
  if (typeof state.ff.rebuildHoloSprings === "function") {
    state.ff.rebuildHoloSprings();
  }

  if (viewer) viewer.setSystem(state.sel, state.ff);
  state._prevPos = null;

  const note = placed.converged
    ? `Placed ${name} — clash-free (residual ${placed.residualClash.toFixed(2)})`
    : `Placed ${name} — relaxed (residual ${placed.residualClash.toFixed(2)}, ${placed.iterations} iters)`;
  if (ui.ligPlaceInfo) {
    // FP2: over-threshold residual keeps the pose (legacy behavior) but the
    // note becomes actionable (what + exact next click, detail secondary).
    const clashErr = checkPlacement(placed);
    ui.ligPlaceInfo.textContent = clashErr ? `${formatInputError(clashErr)} — pose kept (${note}).` : note;
  }
  setPickMode(false);
}

/**
 * Build the placement collision set (exported for headless regression tests).
 * See rev1-issue2 policy note inside.
 */
export function getProteinCoordsAndSigma() {
  const ff = state.ff;
  const nProt = ff.nProt;
  // rev1-issue2 collision policy: clash set = indices [0, collEnd) with
  // collEnd = ligandStart ?? nProt. In heavy mode this is protein + hetero
  // (PDB ligands/cofactors/metals demoted to hetero while an external ligand
  // is present — see heavy.js selectHeavy), EXCLUDING the incoming ligand
  // slot [ligandStart, n) whose stale coords applyPlacedPose overwrites.
  // CG ForceField has no ligandStart/hetero ⇒ collEnd = nProt (protein-only,
  // unchanged legacy behavior).
  const nTotal = ff.n ?? nProt;
  const collEnd = Math.max(0, Math.min(ff.ligandStart ?? nProt, nTotal));
  const posLen = Math.floor(state.integ.pos.length / 3);
  const nColl = Math.max(0, Math.min(collEnd, posLen));
  const sigma = new Float64Array(nColl);
  for (let i = 0; i < nColl; i++) {
    if (ff._elem?.[i] && Number.isFinite(ff._elem[i].sigma)) {
      sigma[i] = ff._elem[i].sigma; // heavy: per-atom size incl. metals
    } else if (ff._protSigma && i < ff._protSigma.length) {
      sigma[i] = ff._protSigma[i];
    } else if (i < nProt) {
      sigma[i] = 3.8; // legacy CG protein-bead diameter
    } else {
      // hetero/existing-ligand beyond _protSigma: element-based fallback.
      const el = state.sel?.atoms?.[i]?.element?.toUpperCase?.();
      sigma[i] = ((el && LIG_ELEMENT[el]) || LIG_ELEMENT_DEFAULT).sigma;
    }
  }
  return {
    pos: state.integ.pos.subarray(0, 3 * nColl),
    sigma,
  };
}

function placeAt(clientX, clientY) {
  if (!state.parsed && !state.parsedHeavy) {
    // FP2: actionable NO_POCKET (empty protein selection at place time).
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("NO_POCKET", "place on viewer with no structure loaded"));
    setPickMode(false);
    return;
  }

  const useMol2 = pickTarget === "mol2";
  let mol, name;
  if (useMol2) {
    const mols = state.mol2Ligands;
    if (!mols || mols.length !== 1) {
      // FP2: missing/empty ligand ⇒ actionable LIGAND_PARSE_FAIL.
      if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("LIGAND_PARSE_FAIL", "place MOL2 with no single-molecule ligand loaded"));
      setPickMode(false);
      return;
    }
    mol = mols[0];
    name = state.mol2Fn || "MOL2 ligand";
  } else {
    const entry = selectedEntry();
    if (!entry) return;
    try {
      mol = parseLibraryLigand(entry);
    } catch (e) {
      if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(classifyInputError(e, { stage: "place" }));
      setPickMode(false);
      return;
    }
    name = entry.name;
  }
  // FP2: empty ligand molecule can never place — actionable, not a throw.
  if (!mol || !mol.atoms || mol.atoms.length === 0) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("LIGAND_PARSE_FAIL", `${name}: ligand has no atoms`));
    setPickMode(false);
    return;
  }

  state.libraryLigand = useMol2 ? null : mol;
  _buildSystem();
  if (!state.integ || !state.ff) {
    state.libraryLigand = null;
    setPickMode(false);
    return;
  }

  const target = viewer ? viewer.screenToWorld(clientX, clientY, { snapToBead: true }) : null;
  if (!target) {
    // FP2: unresolvable pick point ⇒ GENERIC fallback (still actionable).
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(classifyInputError(new Error("Could not resolve world point at cursor."), { stage: "place" }));
    setPickMode(false);
    return;
  }

  const protein = getProteinCoordsAndSigma();
  // FP2: empty protein at place time ⇒ NO_POCKET (not a raw downstream throw).
  if (!protein.pos || protein.pos.length === 0) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("NO_POCKET", "place on viewer with empty protein selection"));
    setPickMode(false);
    return;
  }
  let placed;
  try {
    placed = placeLigand(mol, target, { protein, seed: Math.floor(Math.random() * 1000) });
  } catch (e) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(classifyInputError(e, { stage: "place" }));
    setPickMode(false);
    return;
  }
  applyPlacedPose(mol, name, placed);
}

function placeInPocket() {
  if (!state.ff) {
    // FP2: actionable NO_POCKET (not a bare warning).
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("NO_POCKET", "place in pocket with no built system"));
    return;
  }

  const entry = selectedEntry();
  const mol = state.mol2Ligands && state.mol2Ligands.length === 1 ? state.mol2Ligands[0] : parseLibraryLigand(entry);
  const name = state.mol2Ligands && state.mol2Ligands.length === 1 ? (state.mol2Fn || "MOL2") : (entry ? entry.name : "Ligand");
  // FP2: empty ligand molecule can never place — actionable, not a throw.
  if (!mol || !mol.atoms || mol.atoms.length === 0) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("LIGAND_PARSE_FAIL", `${name}: ligand has no atoms`));
    return;
  }

  if (state.mol2Ligands && state.mol2Ligands.length === 1) state.libraryLigand = null;
  else state.libraryLigand = mol;

  _buildSystem();
  const protein = getProteinCoordsAndSigma();
  // FP2: empty protein at place time ⇒ NO_POCKET (not a raw downstream throw).
  if (!protein.pos || protein.pos.length === 0) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("NO_POCKET", "place in pocket with empty protein selection"));
    return;
  }
  // Cavity search stays protein-only (nProt) while clash relaxation uses the
  // extended protein+hetero set in `protein` (rev1-issue2 explicit policy).
  const pocketCenter = findPocketCenter(protein.pos, state.ff.nProt);
  // FP2: no pocket (empty/degenerate center) ⇒ actionable, not a blind place.
  if (!pocketCenter || pocketCenter.some((v) => !Number.isFinite(v))) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("NO_POCKET", "pocket detector returned no center"));
    return;
  }
  let placed;
  try {
    placed = placeLigand(mol, pocketCenter, { protein, seed: 42 });
  } catch (e) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(classifyInputError(e, { stage: "place" }));
    return;
  }
  applyPlacedPose(mol, `${name} (Pocket)`, placed);
}

function placeRandomSurface() {
  if (!state.ff) {
    // FP2: actionable NO_POCKET (not a bare warning).
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("NO_POCKET", "random-surface place with no built system"));
    return;
  }

  const entry = selectedEntry();
  const mol = state.mol2Ligands && state.mol2Ligands.length === 1 ? state.mol2Ligands[0] : parseLibraryLigand(entry);
  const name = state.mol2Ligands && state.mol2Ligands.length === 1 ? (state.mol2Fn || "MOL2") : (entry ? entry.name : "Ligand");
  // FP2: empty ligand molecule can never place — actionable, not a throw.
  if (!mol || !mol.atoms || mol.atoms.length === 0) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(createInputError("LIGAND_PARSE_FAIL", `${name}: ligand has no atoms`));
    return;
  }

  if (state.mol2Ligands && state.mol2Ligands.length === 1) state.libraryLigand = null;
  else state.libraryLigand = mol;

  _buildSystem();
  const protein = getProteinCoordsAndSigma();
  const nProt = state.ff.nProt;
  const randBead = Math.floor(Math.random() * nProt);
  const bx = protein.pos[3 * randBead];
  const by = protein.pos[3 * randBead + 1];
  const bz = protein.pos[3 * randBead + 2];

  const target = [bx + (Math.random() - 0.5) * 8.0, by + (Math.random() - 0.5) * 8.0, bz + (Math.random() - 0.5) * 8.0];
  let placed;
  try {
    placed = placeLigand(mol, target, { protein, seed: Math.floor(Math.random() * 1000) });
  } catch (e) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = formatInputError(classifyInputError(e, { stage: "place" }));
    return;
  }
  applyPlacedPose(mol, `${name} (Surface)`, placed);
}

export function initLigandPanel(buildSystemFn) {
  _buildSystem = buildSystemFn;

  if (ui.ligFilter) {
    ui.ligFilter.addEventListener("input", renderLibraryList);
  }

  if (ui.placeBtn) {
    ui.placeBtn.addEventListener("click", () => {
      if (!state.parsed && !state.parsedHeavy) return;
      setPickMode(!pickActive || pickTarget !== "library", "library");
    });
  }

  if (ui.placeMol2Btn) {
    ui.placeMol2Btn.addEventListener("click", () => {
      if (!state.parsed && !state.parsedHeavy) return;
      setPickMode(!pickActive || pickTarget !== "mol2", "mol2");
    });
  }

  if (ui.placePocketBtn) {
    ui.placePocketBtn.addEventListener("click", placeInPocket);
  }

  if (ui.placeRandomBtn) {
    ui.placeRandomBtn.addEventListener("click", placeRandomSurface);
  }

  if (ui.cancelPlace) {
    ui.cancelPlace.addEventListener("click", () => setPickMode(false));
  }

  if (viewer?.canvas) {
    viewer.canvas.addEventListener("click", (e) => {
      if (!pickActive) return;
      placeAt(e.clientX, e.clientY);
    });
  }

  renderLibraryList();
  updateMol2PlaceButton();
}