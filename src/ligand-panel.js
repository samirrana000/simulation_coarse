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
  if (!Number.isFinite(placed.residualClash)) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = "⚠ Placement failed (degenerate pose).";
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
  if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = note;
  setPickMode(false);
}

function getProteinCoordsAndSigma() {
  const nProt = state.ff.nProt;
  const protSigma = state.ff._protSigma
    ?? (() => {
        const s = new Float64Array(nProt);
        for (let i = 0; i < nProt; i++) s[i] = state.ff._elem ? state.ff._elem[i].sigma : 3.8;
        return s;
      })();
  return {
    pos: state.integ.pos.subarray(0, 3 * nProt),
    sigma: protSigma,
  };
}

function placeAt(clientX, clientY) {
  if (!state.parsed && !state.parsedHeavy) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = "⚠ Load a structure first (Structure panel), then place.";
    setPickMode(false);
    return;
  }

  const useMol2 = pickTarget === "mol2";
  let mol, name;
  if (useMol2) {
    const mols = state.mol2Ligands;
    if (!mols || mols.length !== 1) {
      if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = "⚠ Load a single-molecule MOL2 ligand before placing.";
      setPickMode(false);
      return;
    }
    mol = mols[0];
    name = state.mol2Fn || "MOL2 ligand";
  } else {
    const entry = selectedEntry();
    if (!entry) return;
    mol = parseLibraryLigand(entry);
    name = entry.name;
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
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = "⚠ Could not resolve world point at cursor.";
    setPickMode(false);
    return;
  }

  const protein = getProteinCoordsAndSigma();
  const placed = placeLigand(mol, target, { protein, seed: Math.floor(Math.random() * 1000) });
  applyPlacedPose(mol, name, placed);
}

function placeInPocket() {
  if (!state.ff) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = "⚠ Load a structure first, then place.";
    return;
  }

  const entry = selectedEntry();
  const mol = state.mol2Ligands && state.mol2Ligands.length === 1 ? state.mol2Ligands[0] : parseLibraryLigand(entry);
  const name = state.mol2Ligands && state.mol2Ligands.length === 1 ? (state.mol2Fn || "MOL2") : (entry ? entry.name : "Ligand");

  if (state.mol2Ligands && state.mol2Ligands.length === 1) state.libraryLigand = null;
  else state.libraryLigand = mol;

  _buildSystem();
  const protein = getProteinCoordsAndSigma();
  const pocketCenter = findPocketCenter(protein.pos, state.ff.nProt);
  const placed = placeLigand(mol, pocketCenter, { protein, seed: 42 });
  applyPlacedPose(mol, `${name} (Pocket)`, placed);
}

function placeRandomSurface() {
  if (!state.ff) {
    if (ui.ligPlaceInfo) ui.ligPlaceInfo.textContent = "⚠ Load a structure first, then place.";
    return;
  }

  const entry = selectedEntry();
  const mol = state.mol2Ligands && state.mol2Ligands.length === 1 ? state.mol2Ligands[0] : parseLibraryLigand(entry);
  const name = state.mol2Ligands && state.mol2Ligands.length === 1 ? (state.mol2Fn || "MOL2") : (entry ? entry.name : "Ligand");

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
  const placed = placeLigand(mol, target, { protein, seed: Math.floor(Math.random() * 1000) });
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