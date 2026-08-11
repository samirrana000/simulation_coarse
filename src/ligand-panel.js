/**
 * ligand-panel.js — "Ligand Library & Placement" panel (item 2).
 *
 * Side-effect module (imported once by main.js) that wires the library
 * selector + filter + "Place on viewer" flow:
 *
 *   1. Populates the <select> from LIGAND_LIBRARY (filtered by the search box).
 *   2. "Place on viewer" arms a pick mode: the next click on the canvas drops
 *      the chosen ligand at the clicked world point, clash-relaxed against the
 *      current protein via placement.js.
 *   3. Placing ADD the library ligand to the system: state.libraryLigand is
 *      set (priority over any MOL2/HETATM ligand) and buildSystem() re-runs so
 *      the force field / integrator include it. The placed pose is then
 *      written back into the live position buffer AND the native reference.
 *
 *   A placed library ligand is a hypothesis, not a known crystallographic
 *   pose — so holo pose springs are DISABLED for it (holoOn=false, springs
 *   cleared), letting it diffuse and bind freely under the binding potentials
 *   instead of being pinned to a fixed bound pose.
 *
 * All DOM wiring mirrors the other panel modules (analysis-panel.js,
 * ml-tier.js): pure side-effect registration, imports `ui`/`state`/`viewer`
 * from ui.js.
 */

import { LIGAND_LIBRARY } from "./ligandLib.js?v=9";
import { parseLibraryLigand, placeLigand } from "./placement.js?v=9";
import { ui, state, viewer } from "./ui.js?v=9";
import { buildSystem } from "./main.js?v=9";

/* ------------------------------------------------------------------ */
/*  Library selector (filterable)                                      */
/* ------------------------------------------------------------------ */

function renderLibraryList() {
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

/** Currently selected library entry (by final <select> value). */
function selectedEntry() {
  return LIGAND_LIBRARY.find((e) => e.id === ui.ligSelect.value) || LIGAND_LIBRARY[0];
}

ui.ligFilter.addEventListener("input", renderLibraryList);

/* ------------------------------------------------------------------ */
/*  Placement interaction                                               */
/* ------------------------------------------------------------------ */

let pickActive = false;

function setPickMode(on) {
  pickActive = on;
  ui.cancelPlace.disabled = !on;
  ui.placeBtn.textContent = on ? "Click the viewer…" : "Place on viewer";
  viewer.canvas.style.cursor = on ? "crosshair" : "";
}

function placeAt(clientX, clientY) {
  const entry = selectedEntry();
  if (!entry || !state.parsed) {
    ui.ligPlaceInfo.textContent = "⚠ Load a structure first (panel 1), then place.";
    return;
  }
  const mol = parseLibraryLigand(entry);

  // ADOPT the library ligand as the system's hypothesis ligand and rebuild so
  // the force field / integrator include it (priority over MOL2/HETATM).
  state.libraryLigand = mol;
  buildSystem();
  if (!state.integ || !state.ff) {
    ui.ligPlaceInfo.textContent = "⚠ Could not build a system with the library ligand.";
    state.libraryLigand = null;
    return;
  }

  const nProt = state.ff.nProt;
  const nLig = state.ff.nLigAtoms;

  // World-space target from the pick (depth snaps to the bead under the
  // cursor → the click lands on the protein surface, not a random plane).
  const target = viewer.screenToWorld(clientX, clientY, { snapToBead: true });
  if (!target) {
    ui.ligPlaceInfo.textContent = "⚠ Could not resolve a world point at the cursor.";
    return;
  }

  // Protein bead positions + per-bead σ for the clash relaxation.
  const protein = {
    pos: state.integ.pos.subarray(0, 3 * nProt),
    sigma: state.ff._protSigma,
  };

  // Clash-relax a rigid pose of the picked molecule centred on the target.
  const placed = placeLigand(mol, target, { protein, seed: 42 });

  if (!Number.isFinite(placed.residualClash)) {
    ui.ligPlaceInfo.textContent = "⚠ Placement failed (degenerate pose).";
    return;
  }

  // Write the placed pose into the live position buffer (ligand atoms only).
  // Molecular concatenation order in mol2 = placement order in the FF, so
  // atom a → global index nProt + a (single molecule, base 0).
  for (let a = 0; a < nLig; a++) {
    state.integ.pos[3 * (nProt + a)] = placed.pos[3 * a];
    state.integ.pos[3 * (nProt + a) + 1] = placed.pos[3 * a + 1];
    state.integ.pos[3 * (nProt + a) + 2] = placed.pos[3 * a + 2];
    state.ff.ref[3 * (nProt + a)] = placed.pos[3 * a];
    state.ff.ref[3 * (nProt + a) + 1] = placed.pos[3 * a + 1];
    state.ff.ref[3 * (nProt + a) + 2] = placed.pos[3 * a + 2];
  }

  // A placed library ligand is a hypothesis → holo pose springs are OFF so it
  // can diffuse/bind freely under the binding potentials.
  state.ff.holoOn = false;
  state.ff.rebuildHoloSprings();

  viewer.setSystem(state.sel, state.ff); // re-colour / re-pivot for new pose
  state._prevPos = null;                 // avoid a bogus Δr spike on next frame

  const note = placed.converged
    ? `placed ${entry.name} — clash-free (residual ${placed.residualClash.toFixed(2)})`
    : `placed ${entry.name} — relaxed (residual ${placed.residualClash.toFixed(2)}, ${placed.iterations} iters)`;
  ui.ligPlaceInfo.textContent = note;
  setPickMode(false);
}

ui.placeBtn.addEventListener("click", () => {
  if (!state.parsed) {
    ui.ligPlaceInfo.textContent = "⚠ Load a structure first (panel 1), then place.";
    return;
  }
  setPickMode(!pickActive);
});

ui.cancelPlace.addEventListener("click", () => setPickMode(false));

viewer.canvas.addEventListener("click", (e) => {
  if (!pickActive) return;
  placeAt(e.clientX, e.clientY);
});

/* ------------------------------------------------------------------ */
/*  Init                                                                */
/* ------------------------------------------------------------------ */

renderLibraryList();