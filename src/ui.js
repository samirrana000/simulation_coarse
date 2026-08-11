/**
 * ui.js — DOM handles, application state, and shared singletons (item 5
 * modularization: moved verbatim from main.js).
 *
 * Every panel module (ml-tier.js, pmf-panel.js, analysis-panel.js) and the
 * orchestrator (main.js) import `ui` / `state` / `viewer` / `recorder` from
 * here. ESM guarantees this module is fully evaluated before any importer's
 * body runs, so event listeners registered by importers always see populated
 * handles.
 */

import { Viewer } from "./viewer.js?v=9";
import { Recorder } from "./recorder.js?v=9";

/* ------------------------------------------------------------------ */
/*  DOM handles                                                        */
/* ------------------------------------------------------------------ */
const $ = (id) => document.getElementById(id);
export const ui = {
  pdbId: $("pdbId"), fetchBtn: $("fetchBtn"), fileInput: $("fileInput"),
  structSummary: $("structSummary"), includeLig: $("includeLig"),
  mol2File: $("mol2File"), mol2Info: $("mol2Info"),
  chainsInput: $("chainsInput"), resFrom: $("resFrom"), resTo: $("resTo"),
  buildBtn: $("buildBtn"), selSummary: $("selSummary"),
  ligFilter: $("ligFilter"), ligSelect: $("ligSelect"),
  placeBtn: $("placeBtn"), cancelPlace: $("cancelPlace"), ligPlaceInfo: $("ligPlaceInfo"),
  rc: $("rc"), gamma: $("gamma"), temp: $("temp"), fric: $("fric"), mass: $("mass"),
  motionGain: $("motionGain"), v_motionGain: $("v_motionGain"),
  stridePs: $("stridePs"), maxFrames: $("maxFrames"), exportFmt: $("exportFmt"),
  recBtn: $("recBtn"), recStopBtn: $("recStopBtn"), dlBtn: $("dlBtn"), recStatus: $("recStatus"),
  playBtn: $("playBtn"), resetBtn: $("resetBtn"),
  showContacts: $("showContacts"), spheres: $("spheres"),
  bindPot: $("bindPot"), holoSprings: $("holoSprings"),
  funnelToggle: $("funnelToggle"), pmfReset: $("pmfReset"), pmfPlot: $("pmfPlot"),
  nnContacts: $("nnContacts"), poseScore: $("poseScore"),
  contactsFile: $("contactsFile"), scorerFile: $("scorerFile"), nnInfo: $("nnInfo"),
  anaBtn: $("anaBtn"), anaPmfBtn: $("anaPmfBtn"), analysisOut: $("analysisOut"),
  hud: $("hud"), canvas: $("canvas"),
};

/* ------------------------------------------------------------------ */
/*  Application state                                                  */
/* ------------------------------------------------------------------ */
export const state = {
  pdbText: null,       // raw structure text
  parsed: null,        // parseCa() output
  sel: null,           // selectSystem() output
  ff: null,            // ForceField
  funnel: null,        // Funnel (binding bias + WTM PMF)
  ligands: [],
  mol2Ligands: null,   // parseMol2() output — when set, replaces HETATM ligands
  mol2Fn: null,        // filename of the loaded MOL2 (for status lines)
  libraryLigand: null, // parsed library molecule (ligand-panel.js) — prioritized
                       // over mol2Ligands/HETATM when set; reset on new structure
  integ: null,         // LangevinIntegrator
  running: false,
  simSpeedPsPerFrame: 1.0, // ~1 ps/frame at 60 fps → 1 ns per ~17 s wall time
  fpsEMA: 60,
  contacts: null,      // parsed contacts.json (array of [i,j,p])
  contactsFn: null,    // filename of loaded contacts
  mlAlpha: 1,          // NN contact map stiffness scale
  scorer: null,       // PoseScorer instance
};

export const viewer = new Viewer(ui.canvas);
export const recorder = new Recorder();

/**
 * Slider live readouts. `onHotParam` is injected by main.js (the hot-reload
 * callback lives there) via initParamReadouts() — keeping it injectable is
 * the only change vs. the original inline version, needed because
 * onParamChange() cannot be imported here without a module cycle.
 */
let _onHotParam = () => {};
export function initParamReadouts(onHotParam) {
  _onHotParam = onHotParam;
  [["rc"], ["gamma"], ["temp", "v_temp"], ["fric", "v_fric"], ["mass", "v_mass"]].forEach(([k, labelId]) => {
    const el = ui[k];
    const lbl = $(labelId || `v_${k}`);
    if (el && lbl) el.addEventListener("input", () => { lbl.textContent = el.value; _onHotParam(); });
  });
}

/** Selection/force-field summary line (panel 2) — reads `state` + `ui`. */
export function updateSelSummary() {
  if (!state.ff) return;
  const nat = state.ff.springs.length / 3;
  let ligSummary = "";
  if (state.ff.nLigAtoms > 0) {
    ligSummary = ` · ${state.ff.nLigAtoms} ligand atom(s) in ${state.ligands.length} molecule(s)`;
    if (state.mol2Ligands && state.mol2Ligands.length) ligSummary += ` (MOL2: ${state.mol2Fn})`;
    if (state.ff.nHolo > 0) ligSummary += ` · ${state.ff.nHolo} holo contacts`;
  }
  let nnSummary = "";
  if (state.ff.springScaleActive) {
    nnSummary = ` · NN map active (${state.contactsFn || "contacts.json"})`;
  }
  ui.selSummary.textContent =
    `${state.sel.beads.length} Cα beads · ${state.sel.nChains} chain(s) · ` +
    `${state.ff.bonds.length / 3} bonds · ${state.ff.angles.length / 4} angles · ` +
    `${nat} elastic contacts (Rc=${Number(ui.rc.value)} Å)${ligSummary}${nnSummary}`;
}

/** Recorder status line (frames · span · ●) — shared by main loop + panel. */
export function updateRecStatus() {
  ui.recStatus.textContent =
    `${recorder.count} frames · ${recorder.spanNs.toFixed(3)} ns recorded` +
    (recorder.recording ? " ●" : "");
  ui.dlBtn.disabled = recorder.count === 0;
}
