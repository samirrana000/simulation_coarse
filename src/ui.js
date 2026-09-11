/**
 * ui.js — DOM handles, application state, and shared singletons.
 */

import { Viewer } from "./viewer.js?v=10";
import { Recorder } from "./recorder.js?v=10";

const $ = (id) => (typeof document !== "undefined" ? document.getElementById(id) : null);

export const ui = {
  pdbId: $("pdbId"), fetchBtn: $("fetchBtn"), fileInput: $("fileInput"),
  structSummary: $("structSummary"), includeLig: $("includeLig"),
  mol2File: $("mol2File"), mol2Info: $("mol2Info"),
  chainsInput: $("chainsInput"), resFrom: $("resFrom"), resTo: $("resTo"),
  modelMode: $("modelMode"), buildBtn: $("buildBtn"), selSummary: $("selSummary"),
  heteroPanel: $("heteroPanel"), heteroList: $("heteroList"),
  heteroAll: $("heteroAll"), heteroMetals: $("heteroMetals"), heteroNone: $("heteroNone"),
  ligFilter: $("ligFilter"), ligSelect: $("ligSelect"),
  placeBtn: $("placeBtn"), placePocketBtn: $("placePocketBtn"), placeRandomBtn: $("placeRandomBtn"), placeMol2Btn: $("placeMol2Btn"), cancelPlace: $("cancelPlace"), ligPlaceInfo: $("ligPlaceInfo"),
  rc: $("rc"), gamma: $("gamma"), temp: $("temp"), fric: $("fric"), mass: $("mass"),
  motionGain: $("motionGain"), v_motionGain: $("v_motionGain"),
  rcNum: $("rcNum"), gammaNum: $("gammaNum"), tempNum: $("tempNum"),
  fricNum: $("fricNum"), massNum: $("massNum"), motionGainNum: $("motionGainNum"),
  stridePs: $("stridePs"), maxFrames: $("maxFrames"), exportFmt: $("exportFmt"),
  recBtn: $("recBtn"), recStopBtn: $("recStopBtn"), dlBtn: $("dlBtn"), recStatus: $("recStatus"),
  playBtn: $("playBtn"), resetBtn: $("resetBtn"),
  showContacts: $("showContacts"), spheres: $("spheres"),
  showRibbon: $("showRibbon"), showHBonds: $("showHBonds"), showStates: $("showStates"),
  bindPot: $("bindPot"), holoSprings: $("holoSprings"),
  respaToggle: $("respaToggle"), respaOuter: $("respaOuter"),
  protAssign: $("protAssign"), gaffLig: $("gaffLig"),
  funnelToggle: $("funnelToggle"), pmfReset: $("pmfReset"), pmfPlot: $("pmfPlot"),
  pmfCaption: $("pmfCaption"),
  nnContacts: $("nnContacts"), poseScore: $("poseScore"),
  contactsFile: $("contactsFile"), scorerFile: $("scorerFile"), nnInfo: $("nnInfo"),
  anaBtn: $("anaBtn"), anaPmfBtn: $("anaPmfBtn"), analysisOut: $("analysisOut"),
  alaScanBtn: $("alaScanBtn"), alaRes: $("alaRes"),
  dccmBtn: $("dccmBtn"), dccmCanvas: $("dccmCanvas"), dccmCaption: $("dccmCaption"),
  smdBtn: $("smdBtn"), crypticBtn: $("crypticBtn"),
  hud: $("hud"), canvas: $("canvas"), canvasCaption: $("canvasCaption"),
  metricsHud: $("metricsHud"),
  sysState: $("sysState"), topPdb: $("topPdb"), topEngine: $("topEngine"), topStep: $("topStep"),
  scrub: $("scrub"), scrubLabel: $("scrubLabel"), cvStrip: $("cvStrip"), hudSpark: $("hudSpark"),
};

export const state = {
  pdbText: null,
  parsed: null,
  sel: null,
  ff: null,
  funnel: null,
  ligands: [],
  mol2Ligands: null,
  mol2Fn: null,
  heavyMode: false,
  parsedHeavy: null,
  protonation: null,
  heteroOverrides: {},
  libraryLigand: null,
  integ: null,
  running: false,
  simSpeedPsPerFrame: 1.0,
  fpsEMA: 60,
  contacts: null,
  contactsFn: null,
  mlAlpha: 1,
  scorer: null,
};

export const viewer = typeof document !== "undefined" && ui.canvas ? new Viewer(ui.canvas) : null;
export const recorder = new Recorder();

if (typeof window !== "undefined") {
  window.__state = state;
  window.__viewer = viewer;
}

let _onHotParam = () => {};
export function initParamReadouts(onHotParam) {
  _onHotParam = onHotParam;
  // Slider + adjacent numeric input pairs. Sliders write Float32 params
  // synchronously (Math.fround) so the live force field sees the exact
  // float32 value; numeric inputs stay in sync both ways.
  const pairs = [
    ["rc", "rcNum", "v_rc"], ["gamma", "gammaNum", "v_gamma"],
    ["temp", "tempNum", "v_temp"], ["fric", "fricNum", "v_fric"],
    ["mass", "massNum", "v_mass"], ["motionGain", "motionGainNum", "v_motionGain"],
  ];
  for (const [k, numId, labelId] of pairs) {
    const el = ui[k];
    const num = ui[numId] || (typeof document !== "undefined" ? document.getElementById(numId) : null);
    if (numId && !ui[numId] && num) ui[numId] = num;
    const lbl = $(labelId || `v_${k}`);
    if (el && lbl) el.addEventListener("input", () => {
      const f32 = Math.fround(Number(el.value));
      lbl.textContent = el.value;
      if (num && document.activeElement !== num) num.value = el.value;
      void f32;
      _onHotParam();
    });
    if (el && num) num.addEventListener("input", () => {
      const v = Number(num.value);
      if (!Number.isFinite(v)) return;
      const lo = Number(el.min), hi = Number(el.max);
      if (Number.isFinite(lo) && v < lo) return;
      if (Number.isFinite(hi) && v > hi) return;
      el.value = String(v);
      if (lbl) lbl.textContent = String(v);
      void Math.fround(v);
      _onHotParam();
    });
  }
}

export function updateSelSummary() {
  if (!state.ff) return;
  if (state.heavyMode) {
    const ff = state.ff;
    const nMetals = (ff.heteroAtoms ?? []).filter((a) => a.isMetal).length;
    let ligSummary = "";
    if (ff.nLigAtoms > 0) {
      ligSummary = ` · ${ff.nLigAtoms} ligand atom(s)`;
      if (state.mol2Ligands && state.mol2Ligands.length && !state.libraryLigand) {
        ligSummary += ` (MOL2: ${state.mol2Fn})`;
      }
    }
    let nnSummary = "";
    if (ff.springScaleActive) {
      nnSummary = ` · NN map active (${state.contactsFn || "contacts.json"})`;
    }
    if (ui.selSummary) {
      ui.selSummary.textContent =
        `${ff.n} heavy atoms · ${ff.nProt} protein · ${ff.nHetero} hetero (${nMetals} metal)` +
        ` · ${ff.bonds.length / 3} bonds · ${ff.angles.length / 4} angles` +
        ` · ${ff.coord.length / 3} metal coord.${ligSummary}${nnSummary}`;
    }
    return;
  }

  const nat = state.ff.springs.length / 3;
  let ligSummary = "";
  if (state.ff.nLigAtoms > 0) {
    ligSummary = ` · ${state.ff.nLigAtoms} ligand atom(s) in ${state.ligands.length} molecule(s) · ${state.ff.nHolo} holo contacts`;
    if (state.mol2Ligands && state.mol2Ligands.length && !state.libraryLigand) {
      ligSummary += ` (MOL2: ${state.mol2Fn})`;
    }
  }
  let nnSummary = "";
  if (state.ff.springScaleActive) {
    nnSummary = ` · NN map active (${state.contactsFn || "contacts.json"})`;
  }
  if (ui.selSummary) {
    ui.selSummary.textContent =
      `${state.ff.nProt} Cα beads · ${nat} elastic springs` +
      ` · ${state.ff.bonds.length / 3} peptide bonds · ${state.ff.angles.length / 4} angles${ligSummary}${nnSummary}`;
  }
}

export function updateRecStatus() {
  if (!ui.recStatus) return;
  const n = recorder.frames.length;
  const ps = recorder.recordedPs;
  ui.recStatus.textContent = `${n} frames · ${(ps / 1000).toFixed(2)} ns recorded`;
  if (ui.dlBtn) ui.dlBtn.disabled = n === 0;
}
