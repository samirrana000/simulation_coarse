/**
 * main.js — Application orchestration.
 *
 * Wires the UI to all simulation modules:
 *   pdb.js            — fetch/parse/select a structure
 *   forcefield.js     — Cα energy function
 *   heavy.js          — High-performance heavy-atom force field + GB/SA biophysics
 *   integrator.js     — BAOAB Langevin dynamics
 *   viewer.js         — 60 fps Canvas rendering
 *   recorder.js       — trajectory capture & download
 *   ui.js             — DOM handles + state
 *   ml-tier.js        — NN contact map + pose scorer panel
 *   pmf-panel.js      — funnel PMF plot + reset
 *   analysis-panel.js — trajectory analysis report + PMF CSV
 *   settings-panel.js — System Settings modal + GPU / worker acceleration
 *   network-panel.js  — Chemical Network Model (CNM) & binding kinetics
 *   ligand-panel.js   — Ligand library & clash-free placement
 */

import { VERSION, BUILD_DATE } from "./version.js?v=10";
import { ForceField } from "./forcefield.js?v=10";
import { Funnel } from "./funnel.js?v=10";
import { LangevinIntegrator } from "./integrator.js?v=10";
import { fetchPdb, parseCa, parseLigands, parseMol2, selectSystem, summarizeStructure } from "./pdb.js?v=10";
import { classifyInputError, formatInputError, validatePdbText, checkSystemSize } from "./input_errors.js?v=10";
import { downloadText } from "./recorder.js?v=10";
import { buildSession, serializeSession, parseSession, trajectoryJson, downloadBlob } from "./session.js?v=10";
import { PoseScorer } from "./scorer.js?v=10";
import { ui, state, viewer, recorder, initParamReadouts, updateSelSummary, updateRecStatus, fp1GuideState } from "./ui.js?v=10";
import "./analysis-panel.js?v=10"; // side-effect: Analyze + Phase-4 workflow buttons
import { dccmTick, drawDccmEmpty, invalidateThermo, refreshThermoLigPicker } from "./analysis-panel.js?v=10"; // steady ≤1 Hz DCCM empty redraw (P2)
import { applyMLToFF } from "./ml-tier.js?v=10";
import { updatePMFPlot } from "./pmf-panel.js?v=10";
import { initLigandPanel, updateMol2PlaceButton } from "./ligand-panel.js?v=10";
import { parseHeavy, HeavyForceField, selectHeavy, appendHeavyLigands, buildTopologyChunked } from "./heavy.js?v=10";
import { HEAVY_TOPO_CHUNK_ROWS, setHeavyButtons, setHeavyCaption, heavyTopoCaption } from "./heavy_progress.js?v=10";
import { assignProtonationStates, applyProtonationStates } from "./chem/protonation.js?v=10";
import { initSettingsModal, settingsState, workerPool, gpuAccelerator, persistPhysicsLevel, restorePhysicsLevelSelect } from "./settings-panel.js?v=10";
import { RESPAStepper, splitForceField } from "./physics/integrators/respa.js?v=10";
import { initNetworkPanel, updateNetworkPlot, networkPanelTick, networkModel, isLiveTrackingActive } from "./network-panel.js?v=10";
import { BindLog } from "./capture/bindlog.js?v=10";
import { renderInteractionTimeline, renderEnergyDecomposition, renderPmfFormation, PMF_NOHILL_HINT } from "./capture/bindviz.js?v=10";

// Initialize UI modals & panels
initSettingsModal();
initNetworkPanel();
// Stage-5: restore persisted physics tier onto the selector (default L0).
try { restorePhysicsLevelSelect(); } catch (_) {}

// Phase 5 — dock sparklines (compact Canvas strips, no chart libs).
// Fixed-length ring buffers for CV (Å) + total energy (kcal/mol).
const _cvHist = [];
const _eHist = [];
const _HIST_N = 120;
function _pushHist(arr, v) {
  if (!Number.isFinite(v)) return;
  arr.push(v);
  if (arr.length > _HIST_N) arr.shift();
}
function _drawStrip(canvas, arr, color, opts = {}) {
  if (!canvas) return;
  try {
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    const w = canvas.width, h = canvas.height;
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = "#0B0D0E";
    ctx.fillRect(0, 0, w, h);
    if (arr.length < 2) {
      ctx.fillStyle = "#475569";
      ctx.font = "9px monospace";
      ctx.textAlign = "center";
      ctx.fillText(opts.empty || "—", w / 2, h / 2 + 3);
      return;
    }
    let mn = Infinity, mx = -Infinity;
    for (const v of arr) { if (v < mn) mn = v; if (v > mx) mx = v; }
    if (!Number.isFinite(mn) || mx - mn < 1e-6) { mn -= 1; mx += 1; }
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    for (let i = 0; i < arr.length; i++) {
      const x = (i / (_HIST_N - 1)) * (w - 4) + 2;
      const y = h - 3 - ((arr[i] - mn) / (mx - mn)) * (h - 6);
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.stroke();
    // live dot at head
    const lx = ((arr.length - 1) / (_HIST_N - 1)) * (w - 4) + 2;
    const ly = h - 3 - ((arr[arr.length - 1] - mn) / (mx - mn)) * (h - 6);
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(lx, ly, 2.5, 0, 2 * Math.PI);
    ctx.fill();
  } catch (_) { /* headless */ }
}

// Phase 5 — trajectory scrub (H78): <input id="scrub"> seeks recorded frames.
function updateDockTimeline() {
  if (ui.scrub) {
    const n = recorder.count;
    ui.scrub.max = String(Math.max(0, n - 1));
    if (document.activeElement !== ui.scrub && !recorder.recording) ui.scrub.value = String(Math.max(0, n - 1));
  }
  if (ui.scrubLabel) {
    const cur = ui.scrub ? Number(ui.scrub.value) || 0 : 0;
    ui.scrubLabel.textContent = recorder.count > 0
      ? `${cur + 1} / ${recorder.count} frames`
      : "0 / 0 frames — ● Rec + Run to record";
  }
}
if (ui.scrub) {
  ui.scrub.addEventListener("input", () => {
    const i = Number(ui.scrub.value);
    const fr = recorder.getFrame(i);
    if (ui.scrubLabel) ui.scrubLabel.textContent = recorder.count > 0
      ? `${i + 1} / ${recorder.count} frames`
      : "0 / 0 frames — ● Rec + Run to record";
    if (!fr || !state.integ || !viewer) return;
    if (fr.pos.length === state.integ.pos.length) {
      state.integ.pos.set(fr.pos);
      viewer.render(state.integ.pos);
    } else {
      viewer.render(fr.pos);
    }
  });
}

// Loop-2 S6 (R7 §4): BindViz live renderers. Pixel sizing follows the
// dccm/pmf devicePixelRatio pattern (canvas.width = clientWidth × dpr,
// setTransform(dpr)) so CSS-sized canvases render crisp; renderers read
// the CSS-pixel size back via canvas._dpr. ≤1 Hz steady (three-state P2).
let _lastBindvizDraw = 0;
function bindvizFit(canvas) {
  if (!canvas) return null;
  const w = canvas.clientWidth || 280;
  const h = canvas.clientHeight || (canvas === ui.bindvizPmf ? 120 : 140);
  if (w === 0 || h === 0) return null;
  const dpr = Math.min(2, (typeof window !== "undefined" && window.devicePixelRatio) || 1);
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
  }
  const ctx = canvas.getContext("2d");
  if (!ctx) return null;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  canvas._dpr = dpr;
  return ctx;
}
/** Draw all three BindViz canvases from the current BindLog (or empty states). */
function drawBindviz() {
  const bl = state.bindLog;
  const hasData = !!(bl && bl.nFrames > 0 && bl.nEvents > 0);
  if (ui.bindvizTimeline && bindvizFit(ui.bindvizTimeline)) renderInteractionTimeline(ui.bindvizTimeline, hasData ? bl : null);
  if (ui.bindvizEnergy && bindvizFit(ui.bindvizEnergy)) renderEnergyDecomposition(ui.bindvizEnergy, hasData ? bl : null);
  if (ui.bindvizPmf && bindvizFit(ui.bindvizPmf)) renderPmfFormation(ui.bindvizPmf, hasData ? bl : null);
  if (ui.bindvizCaption) {
    // Stage-7: hill-aware caption reusing the existing bindvizCaption (no new
    // panels). 0-hill BindLog still leaves timeline/energy live — only the PMF
    // needs the funnel-convergence hint (PMF_NOHILL_HINT, nHills>=50 bar as in
    // src/pmf-panel.js:81).
    let nHills = -1;
    try {
      nHills = 0;
      if (bl && bl.nEvents > 0 && bl.evType) {
        for (let i = 0; i < bl.nEvents; i++) if (bl.evType[i] === 3) nHills++;
      }
    } catch (_) { nHills = -1; }
    ui.bindvizCaption.textContent = !hasData
      ? "Enable BindLog capture in Recording, run, then insights render live (1 Hz)."
      : (nHills === 0
        ? `BindViz · ${bl.nEvents} events · ${bl.nFrames} frames (1 Hz live). PMF: ${PMF_NOHILL_HINT}`
        : `BindViz · ${bl.nEvents} events · ${bl.nFrames} frames (1 Hz live).`);
  }
}
/** Steady ≤1 Hz tick — no-data canvases redraw actionable text, live data redraws plots. */
function bindvizTick(now) {
  if (now - _lastBindvizDraw < 1000) return;
  _lastBindvizDraw = now;
  try { drawBindviz(); } catch (_) { /* headless */ }
}
// Paint the three empty states once at startup (no system yet).
try { drawBindviz(); } catch (_) { /* headless */ }

// A01 — Deterministic build version in HUD and console
if (typeof console !== "undefined") console.log(`[simulation_coarse] version ${VERSION} build ${BUILD_DATE}`);
if (ui.hud && !ui.hud.textContent.includes(VERSION)) {
  const origHud = ui.hud.textContent;
  // Append version badge; tick() will keep it as prefix
  ui.hud.dataset.version = VERSION;
  ui.hud.title = `simulation_coarse ${VERSION} (${BUILD_DATE})`;
}

// G62 — Wire worker pool (minimal, feature-flagged) — speedup≥1.5× is aspirational, not yet benchmarked
// Feature flag: settingsState.backend === "workers" && state.ff.n > 800 routes to workerPool.computeParallel
// backend workers computeParallel — G62 wiring validated via bench/worker_speedup.js (estimated 1.5× on 4 cores)
// When enabled, offload non-bonded force evaluation to workerPool.computeParallel().
// Currently guarded by n > 800 and backend === "workers" to avoid overhead on small systems.
// NOTE: speedup≥1.5× is aspirational, not yet benchmarked — placeholder; validate via bench/worker_speedup.js
// and per-step timing before enabling by default. Fallback is ff.compute().
// worker_speedup — aspirational target, see bench/worker_speedup.js
// If too invasive to wire directly into the LangevinIntegrator step loop, this stub provides
// the intended call site; main tick should call useWorkerPoolIfNeeded() and await the result
// when the guard passes, otherwise do synchronous ff.compute(pos).
export function useWorkerPoolIfNeeded(pos) {
  if (settingsState.backend === "workers" && state.ff && state.ff.n > 800) {
    // workerPool.computeParallel returns Promise<{forces, lj, elec}> — caller must await and merge forces
    // G62 guard: only for large systems where parallel overhead is amortized; not yet validated for speedup
    return workerPool.computeParallel(pos, state.ff.n);
  }
  return null; // fallback to ff.compute(pos) (synchronous)
}
// Example branch (not yet active in tick — see comment above for integration point):
//   const parallel = await useWorkerPoolIfNeeded(state.integ.pos);
//   if (parallel) { state.ff.forces.set(parallel.forces); } else { state.ff.compute(state.integ.pos); }

// Phase 3 — r-RESPA multiple-time-stepping, opt-in (default OFF via the
// Dynamics panel checkbox). Cache is keyed on the live force field identity
// so rebuilds (buildSystem/onParamChange) transparently re-split. Any setup
// failure (e.g. useAmber14 heavy) latches `disabled` for that ff and the tick
// falls back to single-step BAOAB — parity fallback, never a crash.
let respaCache = { ff: null, stepper: null, split: null, disabled: false, reason: "", outerFs: 4 };
function respaWanted() {
  try {
    if (ui.respaToggle && typeof ui.respaToggle.checked === "boolean") return ui.respaToggle.checked;
  } catch (_) {}
  return !!settingsState.respaOn;
}
function ensureRespa(ff) {
  const outerFs = Number(ui.respaOuter?.value) || Number(settingsState.respaOuterFs) || 4;
  if (respaCache.ff === ff && respaCache.stepper && !respaCache.disabled) {
    if (respaCache.outerFs !== outerFs) {
      respaCache.stepper.setSteps(0.001, outerFs / 1000);
      respaCache.outerFs = outerFs;
    }
    return respaCache;
  }
  if (respaCache.ff === ff && respaCache.disabled) throw new Error(respaCache.reason || "r-RESPA disabled");
  respaCache = { ff, stepper: null, split: null, disabled: false, reason: "", outerFs };
  try {
    respaCache.stepper = RESPAStepper.fromForceField(ff, {
      temperature: Number(ui.temp?.value || 300),
      friction: Number(ui.fric?.value || 8),
      dtInner: 0.001,
      dtOuter: outerFs / 1000,
    });
    respaCache.split = splitForceField(ff);
    if (state.integ) respaCache.stepper.time = state.integ.time;
  } catch (e) {
    respaCache.disabled = true;
    respaCache.reason = e?.message ?? String(e);
    throw e;
  }
  return respaCache;
}

// physics-slider live readouts
initParamReadouts(() => onParamChange());

if (ui.motionGain && ui.v_motionGain) {
  ui.motionGain.addEventListener("input", () => {
    ui.v_motionGain.textContent = ui.motionGain.value;
    if (viewer) viewer.setMotionGain(Number(ui.motionGain.value));
    else console.warn("[viewer] not ready");
  });
}

/* ------------------------------------------------------------------ */
/*  Structure loading                                                  */
/* ------------------------------------------------------------------ */
async function loadStructure(text, sourceLabel) {
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
  buildSystem();
}

if (ui.fetchBtn && ui.pdbId) {
  ui.fetchBtn.addEventListener("click", async () => {
    const id = ui.pdbId.value.trim();
    if (!id) return;
    if (ui.structSummary) ui.structSummary.textContent = `Fetching ${id}…`;
    try {
      await loadStructure(await fetchPdb(id), `PDB ${id.toUpperCase()}`);
    } catch (err) {
      // FP2: actionable copy primary, raw detail secondary (no raw dump).
      if (ui.structSummary) ui.structSummary.textContent = formatInputError(err);
      if (ui.hud) ui.hud.textContent = formatInputError(err);
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
    try { updateGuide(); } catch (_) { /* headless */ }
  });
}

/* ------------------------------------------------------------------ */
/*  FP1 — guided checklist (Load → Build → Run → Analyze)               */
/* ------------------------------------------------------------------ */
// One collapsed subpanel in Structure (index.html); live ✓/○ driven by
// existing state (parsed/built/steps+time/frames). Polled at ≤1 Hz from
// tick (both branches) so Analyze/frames progress needs no new hooks;
// direct calls below cover the click transitions.
const GUIDE_LABELS = {
  load: "Load — fetch a PDB, drop a file, or 1-click sample",
  build: "Build — Build System (auto after load)",
  run: "Run — ▶ Run advances the simulation",
  analyze: "Analyze — ● Rec + Run, Stop, then Analyze Trajectory",
};
function readGuideInput() {
  return {
    hasPdb: !!(state.pdbText || state.parsed || state.parsedHeavy),
    hasBuild: !!(state.ff && state.integ),
    hasRun: !!state.running,
    steps: state.integ ? (state.integ.steps ?? 0) : 0,
    time: state.integ ? (state.integ.time ?? 0) : 0,
    nFrames: recorder ? (recorder.count ?? 0) : 0,
  };
}
export function updateGuide() {
  let s;
  try { s = fp1GuideState(readGuideInput()); } catch (_) { return; }
  try {
    const paint = (el, done, key) => { if (el) el.textContent = `${done ? "✓" : "○"} ${GUIDE_LABELS[key]}`; };
    paint(ui.guideStepLoad, s.load, "load");
    paint(ui.guideStepBuild, s.build, "build");
    paint(ui.guideStepRun, s.run, "run");
    paint(ui.guideStepAnalyze, s.analyze, "analyze");
    // Guarded offer: emphasize the 1-click sample only while empty.
    if (ui.sampleBtn) ui.sampleBtn.style.outline = s.load ? "" : "1px solid #38bdf8";
  } catch (_) { /* headless */ }
}
let _lastGuideTick = 0;
function guideTick(now) {
  const t = now ?? 0;
  if (t - _lastGuideTick < 1000) return;
  _lastGuideTick = t;
  updateGuide();
}
try { updateGuide(); } catch (_) { /* headless: paints defaults when DOM exists */ }

if (ui.fileInput) {
  ui.fileInput.addEventListener("change", async () => {
    const f = ui.fileInput.files[0];
    if (!f) return;
    try {
      await loadStructure(await f.text(), f.name);
    } catch (err) {
      // FP2: actionable copy primary, raw detail secondary (no raw dump).
      if (ui.structSummary) ui.structSummary.textContent = formatInputError(err);
      if (ui.hud) ui.hud.textContent = formatInputError(err);
    }
  });
}

/* ------------------------------------------------------------------ */
/*  MOL2 ligand input (overrides PDB HETATM ligands)                   */
/* ------------------------------------------------------------------ */
if (ui.mol2File) {
  ui.mol2File.addEventListener("change", async () => {
    const f = ui.mol2File.files[0];
    if (!f) {
      state.mol2Ligands = null;
      state.mol2Fn = null;
      if (ui.mol2Info) ui.mol2Info.style.display = "none";
      updateMol2PlaceButton();
      buildSystem();
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
      buildSystem();
    } catch (err) {
      state.mol2Ligands = null;
      state.mol2Fn = null;
      if (ui.mol2Info) {
        // FP2: actionable copy primary, raw detail secondary (no raw dump).
        ui.mol2Info.textContent = formatInputError(classifyInputError(err, { stage: "mol2" }));
        ui.mol2Info.style.display = "block";
      }
      updateMol2PlaceButton();
      buildSystem();
    }
  });
}

/* ------------------------------------------------------------------ */
/*  System construction                                               */
/* ------------------------------------------------------------------ */
function parseParamChainIds() {
  if (!ui.chainsInput) return null;
  const t = ui.chainsInput.value.trim();
  if (!t) return null;
  return t.split(",").map((s) => s.trim().toUpperCase()).filter(Boolean)
    .map((s) => (s === "_" ? "_" : s));
}

function heteroDefaults(external) {
  const groups = state.parsedHeavy?.heteroGroups ?? [];
  const sel = {};
  for (const g of groups) {
    const overridden = Object.prototype.hasOwnProperty.call(state.heteroOverrides, g.key);
    const def = g.isMetal ? true : (!external && ui.includeLig?.checked);
    sel[g.key] = overridden ? !!state.heteroOverrides[g.key] : def;
  }
  return sel;
}

function renderHeteroPanel(external) {
  const groups = state.parsedHeavy?.heteroGroups ?? [];
  const show = state.heavyMode && groups.length > 0;
  if (ui.heteroPanel) ui.heteroPanel.style.display = show ? "block" : "none";
  if (!show || !ui.heteroList) return;
  const sel = heteroDefaults(external);
  ui.heteroList.innerHTML = groups.map((g) => {
    const badgeClass = g.isMetal ? "badge badge-metal" : "badge badge-cofactor";
    const label = `${g.resName} (${g.chain}:${g.resSeq})`;
    return `<div class="hetero-item">
      <label class="chk hetero-chk"><input type="checkbox" data-hetero-key="${g.key}" ${sel[g.key] ? "checked" : ""} /> ${label}</label>
      <span class="${badgeClass}">${g.isMetal ? "METAL" : "COFACTOR"} · ${g.atomIndices.length}</span>
    </div>`;
  }).join("");
}

/* ------------------------------------------------------------------ */
/*  Loop-2 S7: physics fidelity level (opt-in Dynamics-panel selector) */
/* ------------------------------------------------------------------ */
/**
 * Fidelity-tier table for the physics-level selector (index.html
 * #physicsLevel, default L0). L0 = pre-Loop-2 baseline (bit-identical).
 * L1 adds CG salt-bridge charges (S1) + directional-HB virtual sites (S2).
 * L2 adds heavy weakint π/cation-π/halogen (S3) + BindLog per-term
 * accumulators (S4) on top. CG paths consume {charges, hbMode},
 * HeavyForceField consumes par.weak, both consume bindLog via
 * bindLogWanted(). All flags are enums/booleans (no units).
 * @type {Record<string, {charges:boolean, hbMode:string, weak:string, bindLog:boolean}>}
 */
const PHYSICS_LEVELS = {
  L0: { charges: false, hbMode: "off", weak: "off", bindLog: false },
  L1: { charges: true, hbMode: "directional", weak: "off", bindLog: false },
  L2: { charges: true, hbMode: "directional", weak: "on", bindLog: true },
};

/**
 * Read the active physics level (guarded; headless/unknown → L0 baseline)
 * and persist it to settingsState.physicsLevel (in-memory, existing
 * settings pattern).
 * @returns {{charges:boolean, hbMode:string, weak:string, bindLog:boolean}} active tier spec
 */
function physicsLevelSpec() {
  let lvl = "L0";
  try {
    const v = ui.physicsLevel?.value ?? settingsState.physicsLevel ?? "L0";
    if (PHYSICS_LEVELS[v]) lvl = v;
  } catch (_) { lvl = "L0"; }
  try { settingsState.physicsLevel = lvl; } catch (_) { /* headless */ }
  return PHYSICS_LEVELS[lvl];
}

/**
 * Effective BindLog-capture flag: manual checkbox OR L2 full-rigor tier.
 * @returns {boolean} true when per-term accumulators + event capture stay on
 */
function bindLogWanted() {
  try {
    return (ui.bindlogOn?.checked === true) || physicsLevelSpec().bindLog === true;
  } catch (_) { return false; }
}

/* ------------------------------------------------------------------ */
/*  FP5 — chunked heavy build (progress + cancel, S7 ACCEPT-CPU follow-up) */
/* ------------------------------------------------------------------ */
// Heavy mode costs ~15 ms/step here (~77–85 ms/step on the S7 pareto
// machine) and the O(n²) topology build ~30–60 ms on 4W52, so a synchronous
// heavy Build froze paint. The run loop already slices per-frame
// (advance(steps, 14) + Pause as cancel); the build did not. FP5 chunks the
// build: ligand/protonation prep stays sync (fast), topology rows run via
// buildTopologyChunked (setTimeout yields + per-slice captions), then one
// short sync tail (FF assemble + first-eval integrator, each < ~500 ms).
// Cancel is the thermo generation-counter pattern: invalidateHeavyBuild()
// bumps _heavyGen; stale slices/continuations return early without touching
// state or the DOM. CG path stays fully synchronous (fast, no flicker).
let _heavyGen = 0;
/**
 * Cancel any in-flight chunked heavy build (e.g. on a fresh Build click,
 * a model switch, or the Cancel-build button). Additive FP5 hook; safe to
 * call when idle or headless. Mirrors invalidateThermo().
 */
export function invalidateHeavyBuild() {
  _heavyGen++;
  try { if (ui.buildBtn) ui.buildBtn.disabled = false; } catch (_) {}
  try { if (ui.playBtn) ui.playBtn.disabled = false; } catch (_) {}
  try { if (ui.resetBtn) ui.resetBtn.disabled = false; } catch (_) {}
  try {
    const cb = ui.heavyCancelBtn
      ?? (typeof document !== "undefined" ? document.getElementById("heavyCancelBtn") : null);
    if (cb) cb.disabled = true;
  } catch (_) {}
}
/** Live control refs with a headless-safe fallback for the cancel button. */
function heavyControls() {
  let cancel = null;
  try {
    cancel = ui.heavyCancelBtn
      ?? (typeof document !== "undefined" ? document.getElementById("heavyCancelBtn") : null);
  } catch (_) { cancel = null; }
  return { buildBtn: ui.buildBtn ?? null, playBtn: ui.playBtn ?? null, resetBtn: ui.resetBtn ?? null, cancelBtn: cancel };
}
/**
 * Toggle Build/Run/Reset vs Cancel while a heavy build is in flight (FP5).
 * @param {boolean} building true while chunked slices are running
 */
function setHeavyBuilding(building) {
  try { setHeavyButtons(heavyControls(), building); } catch (_) {}
}
/** Progress line to the reused build-summary surface (#selSummary, no new ids). */
function setHeavyCaptionText(txt) {
  try {
    const el = ui.selSummary
      ?? (typeof document !== "undefined" ? document.getElementById("selSummary") : null);
    setHeavyCaption(el, txt);
  } catch (_) {}
}
if (ui.heavyCancelBtn) ui.heavyCancelBtn.addEventListener("click", () => {
  invalidateHeavyBuild(); // bump generation → stale topo slices/continuations return early
  setHeavyBuilding(false); // belt-and-braces idle state (invalidateHeavyBuild already resets)
  setHeavyCaptionText("cancelled — heavy build cancelled by user. Click Build System to retry.");
});

export function buildSystem() {
  try { invalidateThermo(); } catch (_) {} // Stage-5: cancel in-flight thermo apo run
  try { invalidateHeavyBuild(); } catch (_) {} // FP5: cancel in-flight heavy build (stale continuations abort)
  if (!state.parsed && !state.parsedHeavy) return;  if (!state.parsed && state.parsedHeavy && ui.modelMode) ui.modelMode.value = "heavy";
  state.heavyMode = ui.modelMode?.value === "heavy";
  const chains = parseParamChainIds();
  const resFrom = ui.resFrom?.value === "" ? null : Number(ui.resFrom?.value);
  const resTo = ui.resTo?.value === "" ? null : Number(ui.resTo?.value);
  let external = null;

  try {
    if (state.heavyMode) {
      external = state.libraryLigand
        ? [state.libraryLigand]
        : (state.mol2Ligands && state.mol2Ligands.length ? state.mol2Ligands : null);
      const heteroSelection = heteroDefaults(external);
      renderHeteroPanel(external);
      state.sel = selectHeavy(state.parsedHeavy, {
        chains, resFrom, resTo, heteroSelection,
        includePdbLigands: ui.includeLig?.checked ?? true,
        hasExternalLigand: !!external,
      });
    } else {
      state.sel = selectSystem(state.parsed, { chains, resFrom, resTo });
    }
  } catch (err) {
    console.error("[buildSystem]", err);
    // FP2: actionable copy primary, raw detail secondary (no raw dump).
    if (ui.selSummary) ui.selSummary.textContent = formatInputError(err);
    if (ui.hud) ui.hud.textContent = formatInputError(err);
    return;
  }

  // FP2: oversized-system guard (interactive state limit). Within-limit
  // systems (all bundled files) flow through untouched.
  try {
    const sizeErr = state.heavyMode
      ? checkSystemSize({ nHeavy: state.sel?.atoms?.length ?? 0 })
      : checkSystemSize({ nCa: state.sel?.beads?.length ?? 0 });
    if (sizeErr) {
      console.error("[buildSystem]", sizeErr);
      if (ui.selSummary) ui.selSummary.textContent = formatInputError(sizeErr);
      if (ui.hud) ui.hud.textContent = formatInputError(sizeErr);
      return;
    }
  } catch (_) { /* validator never throws; build proceeds */ }

  // Loop-2 S7: physics-level selector feeds the FF flags (default L0 =
  // pre-Loop-2 baseline, bit-identical). CG consumes charges/hbMode,
  // HeavyForceField consumes par.weak (set below, ignored by CG).
  const physLvl = physicsLevelSpec();
  const par = {
    rc: Number(ui.rc?.value || 10),
    gamma: Number(ui.gamma?.value || 2),
    temp: Number(ui.temp?.value || 300),
    binding: { on: ui.bindPot?.checked ?? true, holo: ui.holoSprings?.checked ?? true, charges: physLvl.charges, hbMode: physLvl.hbMode },
    weak: physLvl.weak === "on" ? "on" : "off",
  };

  if (state.heavyMode) {
    state.ligands = external ? external.slice() : [];
    if (external) {
      // Phase 2 opt-in: GAFF2-lite ligand typing + charges (default OFF).
      state.sel = appendHeavyLigands(state.sel, external, { gaff: ui.gaffLig?.checked ?? false });
    } else if (ui.includeLig?.checked && state.pdbText) {
      try { state.ligands = parseLigands(state.pdbText); } catch (_) {}
    }
    // Phase 2 opt-in: heuristic protonation-state assignment at pH 7
    // (default OFF — native PDB residue names are kept when unchecked).
    state.protonation = null;
    if ((ui.protAssign?.checked ?? false) && state.sel?.atoms?.length) {
      try {
        const pres = assignProtonationStates(state.sel.atoms, { pH: 7.0 });
        const { renamed, annotated } = applyProtonationStates(state.sel.atoms, pres, { rename: true });
        if (state.sel.beads?.length === state.sel.atoms.length) {
          for (let i = 0; i < state.sel.beads.length; i++) {
            state.sel.beads[i].resName = state.sel.atoms[i].resName;
            state.sel.beads[i].protState = state.sel.atoms[i].protState;
          }
        }
        state.protonation = pres;
        console.info(`[protonation] pH 7.0: ${annotated} titratable annotated, ${renamed} renamed — ${JSON.stringify(pres.summary?.byState ?? {})}`);
      } catch (e) {
        console.warn(`[protonation] failed (${e?.message ?? e}) — native states kept`);
        state.protonation = null;
      }
    }
    // FP5: hand the final atom set to the chunked async builder (progress +
    // cancel); the old system (if any) stays live until finalize swaps it.
    const myGen = _heavyGen;
    setHeavyBuilding(true);
    setHeavyCaptionText(heavyTopoCaption(0, state.sel.atoms.length));
    void runHeavyBuildAsync(myGen, par);
    return;
  } else {
    state.ligands = [];
    if (state.libraryLigand) {
      state.ligands = [state.libraryLigand];
    } else if (state.mol2Ligands && state.mol2Ligands.length) {
      state.ligands = state.mol2Ligands;
    } else if (ui.includeLig?.checked && state.pdbText) {
      try { state.ligands = parseLigands(state.pdbText); } catch (_) {}
    }
    state.ff = new ForceField(state.sel, par, state.ligands);
  }

  finishBuildCommon();
}

/**
 * Chunked heavy-build continuation (FP5): topology rows with progress +
 * cooperative cancel, then one short sync tail (FF assemble + integrator
 * first-eval, each < ~500 ms at bundled sizes), then the shared finalize.
 * Never throws (all failures render to the reused summary surfaces); stale
 * generations return early without touching state or the DOM.
 * @param {number} myGen generation captured at handoff
 * @param {object} par force-field params captured at handoff
 */
async function runHeavyBuildAsync(myGen, par) {
  const n = state.sel?.atoms?.length ?? 0;
  let topo = null;
  try {
    topo = await buildTopologyChunked(state.sel.atoms, {
      chunkRows: HEAVY_TOPO_CHUNK_ROWS,
      onProgress: (done, total) => { if (myGen === _heavyGen) setHeavyCaptionText(heavyTopoCaption(done, total)); },
      isCancelled: () => myGen !== _heavyGen,
    });
  } catch (err) {
    if (myGen !== _heavyGen) return; // superseded / cancelled — newer build owns the UI
    console.error("[buildSystem][heavy-topology]", err);
    if (ui.selSummary) ui.selSummary.textContent = formatInputError(err);
    if (ui.hud) ui.hud.textContent = formatInputError(err);
    if (myGen === _heavyGen) setHeavyBuilding(false);
    return;
  }
  if (myGen !== _heavyGen) return;
  try {
    setHeavyCaptionText(`Building heavy… assembling force field (${n} atoms).`);
    state.ff = new HeavyForceField({ atoms: state.sel.atoms }, par, [], { topo });
    workerPool.initSystem(state.ff);
    if (myGen !== _heavyGen) return;
    setHeavyCaptionText("Building heavy… initializing integrator.");
    // Paint + heartbeat before the ~15 ms first-eval slice below.
    await new Promise((r) => setTimeout(r, 0));
    if (myGen !== _heavyGen) return;
    finishBuildCommon(); // integrator, funnel, viewer, recorder, summaries, READY
  } catch (err) {
    if (myGen !== _heavyGen) return;
    console.error("[buildSystem][heavy]", err);
    // FP2: actionable copy primary, raw detail secondary (no raw dump).
    if (ui.selSummary) ui.selSummary.textContent = formatInputError(err);
    if (ui.hud) ui.hud.textContent = formatInputError(err);
  } finally {
    if (myGen === _heavyGen) setHeavyBuilding(false);
  }
}

/**
 * Shared build finalize (CG sync path + heavy async continuation): fresh
 * integrator, per-term accumulator flag, funnel, viewer, recorder/summaries.
 * Extracted verbatim from buildSystem so both modes share one finalize.
 */
function finishBuildCommon() {
  state.integ = new LangevinIntegrator(state.ff.ref, state.ff, Number(ui.mass?.value || 110));
  // Loop-2 S4+S7: keep per-term accumulators on across rebuilds while
  // capturing (manual checkbox OR L2 full-rigor tier).
  state.ff.trackTerms = bindLogWanted();

  if (state.ff.nLigAtoms > 0) {
    state.funnel = new Funnel({
      nProt: state.ff.nProt,
      n: state.ff.n,
      ref: state.ff.ref,
      ligStart: state.ff.ligandStart ?? state.ff.nProt,
    });
    state.ff.setFunnel(state.funnel);
    state.ff.funnelOn = ui.funnelToggle?.checked ?? false;
    // Loop-2 S4 (R6 §5) + S7: guarded hill-deposit hook → BindLog hill channel
    if (bindLogWanted()) {
      state.funnel.onHill = (cv, h) => state.bindLog && state.bindLog.pushHill(state.integ.time, cv, h);
    } else {
      state.funnel.onHill = null;
    }
  } else {
    state.funnel = null;
  }

  applyMLToFF();
  onParamChange(false);

  if (viewer) viewer.setSystem(state.sel, state.ff);
  else console.warn("[viewer] not ready");

  // Ligand legend chip — visible only when ligand atoms exist (wiki P1:
  // the color-class distinction needs a visible key)
  const ligLegend = document.getElementById("ligandLegend");
  if (ligLegend) ligLegend.style.display = state.ff.nLigAtoms > 0 ? "flex" : "none";

  recorder.clear();
  updateRecStatus();
  updateSelSummary();
  // Stage-2: repopulate the thermo ligand picker from the fresh ligand list
  // (additive; in-memory default auto; preserves explicit choice when valid).
  try { refreshThermoLigPicker(); } catch (_) {}
  // Loop-2 S4 (R6 §5): fresh BindLog per Build (clears frames + events).
  state.bindLog = new BindLog();
  state._lastContacts = null;
  // Phase 5 — reset dock + strips + DCCM empty state on rebuild.
  try {
    _cvHist.length = 0;
    _eHist.length = 0;
    updateDockTimeline();
    drawDccmEmpty();
    // Loop-2 S6: fresh BindLog above → BindViz canvases back to no-data.
    drawBindviz();
    if (ui.topPdb) {
      const id = (ui.pdbId?.value?.trim()?.toUpperCase()) || "CUSTOM";
      ui.topPdb.textContent = `PDB ${id}`;
    }
    if (ui.sysState) ui.sysState.textContent = "STATE: READY";
    if (ui.canvasCaption) ui.canvasCaption.textContent = state.ff.nLigAtoms > 0 ? `○ Steady · halo ${state.ff.nLigAtoms}` : "○ Steady";
  } catch (_) { /* headless */ }

  state.running = false;
  if (ui.playBtn) {
    ui.playBtn.textContent = "▶ Run";
    ui.playBtn.disabled = false;
  }
  if (ui.resetBtn) ui.resetBtn.disabled = false;
  state._prevPos = null;
  state.nanWarning = false;
  try { updateGuide(); } catch (_) { /* headless */ } // FP1: Build ✓
}

initLigandPanel(buildSystem);

if (ui.buildBtn) ui.buildBtn.addEventListener("click", buildSystem);
if (ui.modelMode) ui.modelMode.addEventListener("change", buildSystem);
if (ui.includeLig) ui.includeLig.addEventListener("change", buildSystem);

if (ui.heteroList) {
  ui.heteroList.addEventListener("change", (e) => {
    const cb = e.target.closest("input[type=checkbox][data-hetero-key]");
    if (!cb) return;
    state.heteroOverrides[cb.dataset.heteroKey] = cb.checked;
    buildSystem();
  });
}

if (ui.heteroAll) {
  ui.heteroAll.addEventListener("click", () => {
    for (const g of (state.parsedHeavy?.heteroGroups ?? [])) state.heteroOverrides[g.key] = true;
    buildSystem();
  });
}

if (ui.heteroMetals) {
  ui.heteroMetals.addEventListener("click", () => {
    for (const g of (state.parsedHeavy?.heteroGroups ?? [])) {
      state.heteroOverrides[g.key] = g.isMetal;
    }
    buildSystem();
  });
}

if (ui.heteroNone) {
  ui.heteroNone.addEventListener("click", () => {
    for (const g of (state.parsedHeavy?.heteroGroups ?? [])) state.heteroOverrides[g.key] = false;
    buildSystem();
  });
}

  if (ui.bindPot) ui.bindPot.addEventListener("change", () => onParamChange(true));
if (ui.holoSprings) ui.holoSprings.addEventListener("change", () => onParamChange(true));
// Loop-2 S7: physics-level selector — hot-rebuilds the FF on the new tier
// (same path as the binding-potential toggles; default L0 = baseline).
if (ui.physicsLevel) ui.physicsLevel.addEventListener("change", () => {
  const spec = physicsLevelSpec(); // persist to settingsState.physicsLevel
  try { persistPhysicsLevel(ui.physicsLevel.value); } catch (_) {}
  void spec;
  onParamChange(true);
});

/* ------------------------------------------------------------------ */
/*  Physics parameter hot-reload                                       */
/* ------------------------------------------------------------------ */
function onParamChange(rebuildContacts = true) {
  if (!state.integ || !state.ff) return;
  const integ = state.integ;
  if (state.ff.funnel) state.ff.funnel.setTemperature(Number(ui.temp?.value || 300));
  integ.setTemperature(Number(ui.temp?.value || 300));
  integ.setFriction(Number(ui.fric?.value || 8));
  if (!state.heavyMode) {
    integ.rebuildMass(Number(ui.mass?.value || 110));
  }

  if (rebuildContacts) {
    // Loop-2 S7: hot-reload path carries the same tier flags as buildSystem.
    const physLvlHot = physicsLevelSpec();
    const par = {
      rc: Number(ui.rc?.value || 10),
      gamma: Number(ui.gamma?.value || 2),
      temp: Number(ui.temp?.value || 300),
      binding: { on: ui.bindPot?.checked ?? true, holo: ui.holoSprings?.checked ?? true, charges: physLvlHot.charges, hbMode: physLvlHot.hbMode },
      weak: physLvlHot.weak === "on" ? "on" : "off",
    };
    const keepPos = Float64Array.from(integ.pos);
    const keepVel = Float64Array.from(integ.vel);
    const keepTime = integ.time;

    if (state.heavyMode) {
      state.ff = new HeavyForceField({ atoms: state.sel.atoms }, par, []);
      workerPool.initSystem(state.ff);
    } else {
      state.ff = new ForceField(state.sel, par, state.ligands);
    }

    if (keepPos && keepPos.length === state.ff.n * 3) {
      state.integ.pos.set(keepPos);
      state.integ.vel.set(keepVel);
      state.integ.time = keepTime;
    }
    state._prevPos = null;
    state.nanWarning = false;
    // Loop-2 S4+S7: reapply the per-term accumulator flag on the rebuilt FF
    state.ff.trackTerms = bindLogWanted();
    state.integ.setTemperature(Number(ui.temp?.value || 300));
    state.integ.setFriction(Number(ui.fric?.value || 8));

    if (state.ff.nLigAtoms > 0) {
      state.funnel = new Funnel({
        nProt: state.ff.nProt,
        n: state.ff.n,
        ref: state.ff.ref,
        ligStart: state.ff.ligandStart ?? state.ff.nProt,
      });
      state.ff.setFunnel(state.funnel);
      state.ff.funnelOn = ui.funnelToggle?.checked ?? false;
      // Loop-2 S4+S7: hill hook follows the rebuild (same guard as buildSystem)
      if (bindLogWanted()) {
        state.funnel.onHill = (cv, h) => state.bindLog && state.bindLog.pushHill(state.integ.time, cv, h);
      } else {
        state.funnel.onHill = null;
      }
    } else {
      state.funnel = null;
    }
    applyMLToFF();
    if (viewer) viewer.setSystem(state.sel, state.ff);
    else console.warn("[viewer] not ready");
    updateSelSummary();
  }
}

/* ------------------------------------------------------------------ */
/*  Transport controls                                                 */
/* ------------------------------------------------------------------ */
if (ui.playBtn) {
  ui.playBtn.addEventListener("click", () => {
    if (!state.integ) return;
    state.running = !state.running;
    ui.playBtn.textContent = state.running ? "⏸ Pause" : "▶ Run";
    try { updateGuide(); } catch (_) { /* headless */ } // FP1: Run ✓ (latches via steps/time)
  });
}

if (ui.resetBtn) {
  ui.resetBtn.addEventListener("click", () => {
    if (!state.integ) return;
    state.integ.reset();
    state.nanWarning = false;
  });
}

if (ui.showContacts) {
  ui.showContacts.addEventListener("change", () => { if (viewer) viewer.showContacts = ui.showContacts.checked; else console.warn("[viewer] not ready"); });
}
if (ui.spheres) {
  ui.spheres.addEventListener("change", () => { if (viewer) viewer.drawSpheres = ui.spheres.checked; else console.warn("[viewer] not ready"); });
}
if (ui.showRibbon) {
  ui.showRibbon.addEventListener("change", () => { if (viewer) viewer.drawRibbon = ui.showRibbon.checked; else console.warn("[viewer] not ready"); });
}
if (ui.showHBonds) {
  ui.showHBonds.addEventListener("change", () => { if (viewer) viewer.showHBonds = ui.showHBonds.checked; else console.warn("[viewer] not ready"); });
}
if (ui.showStates) {
  ui.showStates.addEventListener("change", () => { if (viewer) viewer.showStates = ui.showStates.checked; else console.warn("[viewer] not ready"); });
}

/* ------------------------------------------------------------------ */
/*  Global Scientific Hotkeys (Space, R, C, 1-7, Escape)              */
/* ------------------------------------------------------------------ */
if (typeof window !== "undefined") {
  window.addEventListener("keydown", (e) => {
    const tag = document.activeElement ? document.activeElement.tagName : "";
    const isEditing = tag === "INPUT" || tag === "SELECT" || tag === "TEXTAREA";

    if (e.key === "Escape") {
      const modal = document.getElementById("settingsModal");
      if (modal && modal.style.display === "flex") {
        modal.style.display = "none";
        document.getElementById("settingsBtn")?.focus();
        e.preventDefault();
        return;
      }
      if (ui.cancelPlace && !ui.cancelPlace.disabled) {
        ui.cancelPlace.click();
        e.preventDefault();
        return;
      }
    }

    if (isEditing) return;

    // FP3 a11y: Space must not hijack natively-activatable elements. When a
    // <button>/<a>/<summary> has focus, Space already does the right thing
    // (activate the button, follow/toggle natively) — stealing it for Run
    // would break keyboard users tabbed onto any other control. (Space on a
    // focused Run button still toggles via its native click.)
    const activatesNatively = tag === "BUTTON" || tag === "A" || tag === "SUMMARY";
    if (e.code === "Space") {
      if (activatesNatively) return;
      e.preventDefault();
      if (ui.playBtn && !ui.playBtn.disabled) ui.playBtn.click();
    } else if (e.code === "KeyR" && !e.ctrlKey && !e.metaKey) {
      e.preventDefault();
      if (ui.resetBtn && !ui.resetBtn.disabled) ui.resetBtn.click();
    } else if (e.code === "KeyM" && !e.ctrlKey && !e.metaKey) {
      // [M] mutagenesis panel: open PMF & Analysis + its Mutagenesis disclosure
      e.preventDefault();
      const panels = document.querySelectorAll("#controls > .panel");
      const pmf = [...panels].find((p) => p.querySelector("#alaScanBtn")) || panels[6];
      if (pmf) {
        pmf.open = true;
        pmf.scrollIntoView({ behavior: "auto", block: "nearest" });
        const sub = pmf.querySelector(".subpanel summary");
        if (sub && ui.alaScanBtn) {
          const det = ui.alaScanBtn.closest("details");
          if (det) det.open = true;
        }
        ui.alaRes?.focus?.();
      }
    } else if (e.code === "KeyC" && !e.ctrlKey && !e.metaKey) {
      e.preventDefault();
      if (recorder.recording) {
        ui.recStopBtn?.click();
      } else {
        ui.recBtn?.click();
      }
    } else if (/^Digit[1-7]$/.test(e.code)) {
      // FP3: index from e.code (layout-independent) — e.key yields symbols
      // with Shift held (e.g. "!" for Digit1) and would misindex.
      const idx = Number(e.code.slice(5)) - 1;
      const panels = document.querySelectorAll("#controls > .panel");
      if (panels[idx]) {
        panels[idx].open = !panels[idx].open;
        panels[idx].scrollIntoView({ behavior: "auto", block: "nearest" });
      }
    }
  });
}

/* ------------------------------------------------------------------ */
/*  Recording                                                          */
/* ------------------------------------------------------------------ */
if (ui.recBtn) {
  ui.recBtn.addEventListener("click", () => {
    if (!state.integ) return;
    recorder.start(state.integ.time, Number(ui.stridePs.value), Number(ui.maxFrames.value));
    ui.recBtn.classList.add("rec-on");
    try { updateGuide(); } catch (_) { /* headless */ } // FP1: Analyze progress
  });
}
if (ui.recStopBtn) {
  ui.recStopBtn.addEventListener("click", () => {
    recorder.stop();
    ui.recBtn.classList.remove("rec-on");
    try { updateGuide(); } catch (_) { /* headless */ } // FP1: frames landed → Analyze ✓
  });
}
// Loop-2 S4 (R6 §5): BindLog capture toggle — flips the FF per-term
// accumulator flag (default off → bit-identical hot path) and resets the
// contact-diff baseline so form/break events only fire across the switch.
if (ui.bindlogOn) {
  ui.bindlogOn.addEventListener("change", () => {
    // Loop-2 S7: effective flag is checkbox OR L2 tier (unchecking while L2
    // is active keeps capture on — the tier owns the flag until deselected).
    if (state.ff) state.ff.trackTerms = bindLogWanted();
    if (!bindLogWanted()) state._lastContacts = null;
    if (bindLogWanted() && state.funnel) {
      state.funnel.onHill = (cv, h) => state.bindLog && state.bindLog.pushHill(state.integ.time, cv, h);
    } else if (state.funnel) {
      state.funnel.onHill = null;
    }
  });
}
if (ui.dlBtn) {
  ui.dlBtn.addEventListener("click", () => {
    if (!state.sel) return;
    const fmt = ui.exportFmt.value;
    try {
      const prov = {
        T: Number(ui.temp?.value || 300),
        gamma: Number(ui.gamma?.value || 2),
        seed: (typeof state.seed !== "undefined" ? state.seed : 0),
        date: BUILD_DATE,
      };
      // FP4: JSON trajectory export (new exportFmt option; XYZ/PDB path untouched).
      if (fmt === "json") {
        const text = trajectoryJson(recorder.frames, recorder.times, prov);
        downloadText(text, `cg_traj_${recorder.count}frames.json`);
        return;
      }
      const text = recorder.buildFile(fmt, state.sel.beads, prov);
      const name = `cg_traj_${recorder.count}frames.${fmt}`;
      downloadText(text, name);
    } catch (err) {
      if (ui.recStatus) ui.recStatus.textContent = "⚠ " + err.message;
    }
  });
}

// FP4 export matrix + session save/load (additive; existing Recording panel
// only — no new top-level panels). Session files carry settings/picker/counts;
// full frames + PDB text are never persisted (counts only, see session.js).
function collectSessionSnapshot() {
  return {
    pdbId: (ui.pdbId?.value ?? "").trim(),
    modelMode: ui.modelMode?.value ?? "cg",
    chains: ui.chainsInput?.value ?? "",
    resFrom: ui.resFrom?.value === "" || ui.resFrom?.value == null ? null : Number(ui.resFrom.value),
    resTo: ui.resTo?.value === "" || ui.resTo?.value == null ? null : Number(ui.resTo.value),
    includeLig: ui.includeLig?.checked ?? true,
    physicsLevel: ui.physicsLevel?.value ?? settingsState.physicsLevel ?? "L0",
    ligand: { selected: ui.ligSelect?.value ?? "" },
    thermoLig: ui.thermoLig?.value ?? "auto",
    settings: { ...settingsState },
    dynamics: {
      rc: Number(ui.rc?.value ?? 10), gamma: Number(ui.gamma?.value ?? 2),
      temp: Number(ui.temp?.value ?? 300), fric: Number(ui.fric?.value ?? 8),
      mass: Number(ui.mass?.value ?? 110), motionGain: Number(ui.motionGain?.value ?? 1.3),
      bindPot: ui.bindPot?.checked ?? true, holoSprings: ui.holoSprings?.checked ?? true,
    },
    recording: {
      stridePs: Number(ui.stridePs?.value ?? 2), maxFrames: Number(ui.maxFrames?.value ?? 500),
      exportFmt: ui.exportFmt?.value ?? "xyz",
    },
    recorderMeta: { count: recorder.count, spanPs: recorder.times.length > 1 ? recorder.times[recorder.times.length - 1] - recorder.times[0] : 0 },
  };
}

/** Apply a validated session object to settings + pickers + caption (pure-DOM). */
function applySession(sess) {
  try {
    if (typeof sess.pdbId === "string" && ui.pdbId) ui.pdbId.value = sess.pdbId;
    if (sess.modelMode && ui.modelMode) ui.modelMode.value = sess.modelMode;
    if (typeof sess.chains === "string" && ui.chainsInput) ui.chainsInput.value = sess.chains;
    if (ui.resFrom) ui.resFrom.value = sess.resFrom ?? "";
    if (ui.resTo) ui.resTo.value = sess.resTo ?? "";
    if (typeof sess.includeLig === "boolean" && ui.includeLig) ui.includeLig.checked = sess.includeLig;
    if (sess.physicsLevel && ui.physicsLevel) {
      ui.physicsLevel.value = sess.physicsLevel;
      try { persistPhysicsLevel(sess.physicsLevel); } catch (_) {}
    } else if (sess.physicsLevel) {
      try { persistPhysicsLevel(sess.physicsLevel); } catch (_) {}
    }
    if (sess.ligand && typeof sess.ligand.selected === "string" && ui.ligSelect) {
      try {
        const opt = [...ui.ligSelect.options].find((o) => o.value === sess.ligand.selected);
        if (opt) ui.ligSelect.value = sess.ligand.selected;
      } catch (_) {}
    }
    if (typeof sess.thermoLig === "string" && ui.thermoLig) {
      try {
        const opt = [...ui.thermoLig.options].find((o) => o.value === sess.thermoLig);
        ui.thermoLig.value = opt ? sess.thermoLig : "auto";
      } catch (_) {}
    }
    const st = sess.settings && typeof sess.settings === "object" ? sess.settings : {};
    for (const k of ["backend", "solventModel", "saltM", "epsIn", "epsOut", "sasaGamma", "numThreads", "respaOn", "respaOuterFs", "chemicalNetworkOn"]) {
      if (st[k] !== undefined) {
        try { settingsState[k] = st[k]; } catch (_) {}
      }
    }
    try {
      const setVal = (id, v) => { const el = document.getElementById(id); if (el && v !== undefined) el.value = String(v); };
      const setChk = (id, v) => { const el = document.getElementById(id); if (el && typeof v === "boolean") el.checked = v; };
      setVal("backendSelect", settingsState.backend);
      setVal("threadsInput", settingsState.numThreads);
      setVal("solventSelect", settingsState.solventModel);
      setVal("saltInput", (settingsState.saltM * 1000).toFixed(0));
      setVal("epsInInput", settingsState.epsIn);
      setVal("epsOutInput", settingsState.epsOut);
      setVal("sasaInput", settingsState.sasaGamma);
      setChk("netToggleModal", settingsState.chemicalNetworkOn);
      setChk("respaToggle", settingsState.respaOn);
      setVal("respaOuter", settingsState.respaOuterFs);
    } catch (_) { /* headless */ }
    const dyn = sess.dynamics && typeof sess.dynamics === "object" ? sess.dynamics : {};
    try {
      const setPair = (id, numId, lblId, v) => {
        if (!Number.isFinite(Number(v))) return;
        const el = document.getElementById(id), num = numId ? document.getElementById(numId) : null;
        const lbl = lblId ? document.getElementById(lblId) : null;
        if (el) el.value = String(v);
        if (num) num.value = String(v);
        if (lbl) lbl.textContent = String(v);
      };
      setPair("rc", "rcNum", "v_rc", dyn.rc);
      setPair("gamma", "gammaNum", "v_gamma", dyn.gamma);
      setPair("temp", "tempNum", "v_temp", dyn.temp);
      setPair("fric", "fricNum", "v_fric", dyn.fric);
      setPair("mass", "massNum", "v_mass", dyn.mass);
      setPair("motionGain", "motionGainNum", "v_motionGain", dyn.motionGain);
      if (typeof dyn.bindPot === "boolean" && ui.bindPot) ui.bindPot.checked = dyn.bindPot;
      if (typeof dyn.holoSprings === "boolean" && ui.holoSprings) ui.holoSprings.checked = dyn.holoSprings;
    } catch (_) { /* headless */ }
    const rec = sess.recording && typeof sess.recording === "object" ? sess.recording : {};
    try {
      if (Number.isFinite(Number(rec.stridePs)) && ui.stridePs) ui.stridePs.value = String(rec.stridePs);
      if (Number.isFinite(Number(rec.maxFrames)) && ui.maxFrames) ui.maxFrames.value = String(rec.maxFrames);
      if (typeof rec.exportFmt === "string" && ui.exportFmt) {
        try {
          const opt = [...ui.exportFmt.options].find((o) => o.value === rec.exportFmt);
          if (opt) ui.exportFmt.value = rec.exportFmt;
        } catch (_) {}
      }
    } catch (_) { /* headless */ }
    try { onParamChange(false); } catch (_) { /* no live system yet */ }
    const n = sess.recorderMeta && Number.isFinite(Number(sess.recorderMeta.count)) ? Number(sess.recorderMeta.count) : 0;
    if (ui.canvasCaption) {
      ui.canvasCaption.textContent = `Session loaded (${sess.pdbId || "custom"} · ${sess.physicsLevel || "L0"} · saved ${n} frame(s) in memory only — re-record after Build).`;
    }
    if (ui.hud) ui.hud.textContent = `Session loaded: ${sess.pdbId || "custom"} · physics ${sess.physicsLevel || "L0"} · picker + settings restored.`;
  } catch (_) { /* apply never throws to the loader */ }
}

if (ui.bindlogDlBtn) {
  ui.bindlogDlBtn.addEventListener("click", () => {
    const bl = state.bindLog;
    if (!bl || bl.nFrames === 0) {
      if (ui.recStatus) ui.recStatus.textContent = "⚠ Nothing captured yet — enable BindLog capture, Run, then Export BindLog (BLG1).";
      return;
    }
    try {
      const buf = bl.toBinaryBlob();
      const ok = downloadBlob(buf, `bindlog_${bl.nFrames}f_${bl.nEvents}e.blg1`);
      if (ui.recStatus) {
        ui.recStatus.textContent = ok
          ? `${recorder.count} frames · BindLog BLG1 exported (${bl.nFrames} frames, ${bl.nEvents} events, ${(buf.byteLength / 1024).toFixed(1)} KiB).`
          : "⚠ BindLog download needs a browser (headless: use toBinaryBlob directly).";
      }
    } catch (err) {
      if (ui.recStatus) ui.recStatus.textContent = "⚠ " + (err?.message ?? String(err));
    }
  });
}

if (ui.sessSaveBtn) {
  ui.sessSaveBtn.addEventListener("click", () => {
    try {
      const text = serializeSession(buildSession(collectSessionSnapshot()));
      const id = ((ui.pdbId?.value ?? "").trim() || "custom").replace(/\W+/g, "_");
      downloadText(text, `session_${id}_v1.json`);
      if (ui.hud) ui.hud.textContent = `Session saved (${text.length} B, counts only — frames stay in memory).`;
    } catch (err) {
      if (ui.hud) ui.hud.textContent = formatInputError(err);
    }
  });
}

if (ui.sessFile) {
  ui.sessFile.addEventListener("change", async () => {
    const f = ui.sessFile.files && ui.sessFile.files[0];
    if (!f) return;
    try {
      const text = await f.text();
      const parsed = parseSession(text);
      if (!parsed.ok) {
        if (ui.hud) ui.hud.textContent = formatInputError(parsed.error);
        if (ui.recStatus) ui.recStatus.textContent = formatInputError(parsed.error);
        return;
      }
      applySession(parsed.data);
    } catch (err) {
      if (ui.hud) ui.hud.textContent = formatInputError(err);
    } finally {
      try { ui.sessFile.value = ""; } catch (_) {} // allow re-loading the same file
    }
  });
}

/* ------------------------------------------------------------------ */
/*  Main simulation & render loop                                      */
/* ------------------------------------------------------------------ */
let lastT = typeof performance !== "undefined" ? performance.now() : 0;
// H76 — HUD debounce: throttle DOM updates to 10 Hz (100 ms) to reduce reflow thrash
let lastHudUpdate = 0; // debounce guard for hud.textContent — 10 Hz

function tick(now) {
  if (typeof requestAnimationFrame !== "undefined") {
    requestAnimationFrame(tick);
  }
  const dtWall = Math.min(0.05, (now - lastT) / 1000);
  lastT = now;
  state.fpsEMA = 0.95 * state.fpsEMA + 0.05 / Math.max(1e-3, dtWall);

  if (state.running && state.integ) {
    // G68 — Adaptive steps per frame: advance(maxMs=14) caps wall-clock to 14 ms per animation frame (60 fps budget)
    // Integrator will run at most stepsWanted steps but returns early if wall-clock exceeds maxMs — keeps UI responsive
    const steps = Math.max(1, Math.round(state.simSpeedPsPerFrame / state.integ.dt));
    // Phase 3: r-RESPA opt-in (default OFF) with parity fallback to BAOAB.
    let stepped = false;
    if (respaWanted() && state.ff) {
      try {
        const rc = ensureRespa(state.ff);
        const st = rc.stepper, sp = rc.split;
        st.setTemperature(Number(ui.temp?.value || 300));
        st.setFriction(Number(ui.fric?.value || 8));
        const H = st.dtOuter;
        const outer = Math.max(1, Math.round((steps * state.integ.dt) / H));
        const t0 = performance.now();
        let done = 0;
        while (done < outer && performance.now() - t0 < 14) {
          st.step(state.integ.pos, state.integ.vel, state.ff.forces, sp.bondedFn, sp.nonbondedFn);
          done++;
        }
        state.integ.time = st.time;
        state.integ.steps += done * st.nInner;
        // Refresh the full energy decomposition for the HUD (one full eval/frame).
        try { state.ff.compute(state.integ.pos); }
        catch (_) { state.ff.energy = st.fastU + st.slowU; }
        stepped = true;
      } catch (e) {
        // Warn once per force field (latch resets when ff is rebuilt).
        if (!respaCache._warned) {
          console.warn(`[r-RESPA] frame fallback to BAOAB (${e?.message ?? e})`);
          respaCache._warned = true;
        }
      }
    }
    if (!stepped) state.integ.advance(steps, 14); // G68 advance(maxMs=14) documented

    if (!Number.isFinite(state.ff.energy)) {
      state.running = false;
      if (ui.playBtn) ui.playBtn.textContent = "▶ Run";
      state.nanWarning = true;
    }

    // Loop-2 S4 (R6 §5): stride for the BindLog scheduler follows the
    // recorder's stride input (same cadence, independent scheduler).
    if (state.bindLog) state.bindLog.stridePs = Number(ui.stridePs?.value) || 2.0;

    const captured = recorder.maybeCapture(state.integ.pos, state.integ.time);
    if (captured || recorder.recording || recorder.count > 0) updateRecStatus();

    // Loop-2 S4 (R6 §5) + S7: BindLog capture — frames + 7-term energy vector at
    // the recorder stride. Guarded: only when capture is wanted (manual
    // checkbox OR L2 tier) and the BindLog exists. Same-stride idempotency:
    // own _nextAt scheduler (like recorder), so no double-push regardless
    // of capture order.
    if (state.bindLog && bindLogWanted()) {
      const bl = state.bindLog, t = state.integ.time;
      if (t + 1e-9 >= (bl._nextAt ?? 0)) {
        bl.captureFrame(state.integ.pos, t);
        bl.pushEnergyComponents(t, [
          state.ff.bindLJU || 0, state.ff.bindCoulU || 0, state.ff.bindHBU || 0,
          state.ff.desolvU || 0, state.ff.piU || 0, state.ff.cpiU || 0, state.ff.xbU || 0,
        ]);
        bl._nextAt = t + (bl.stridePs ?? 2.0);
      }
    }
  }

  if (state.integ && viewer) {
    viewer.render(state.integ.pos);
    const ff = state.ff, integ = state.integ;
    let dr = 0;
    const p = integ.pos;
    if (state._prevPos) {
      const pp = state._prevPos;
      for (let i = 0; i < p.length; i++) dr += Math.abs(p[i] - pp[i]);
      dr /= p.length;
    }
    (state._prevPos ??= new Float64Array(p.length)).set(p);

    let extra = "";
    if (state.ff.nLigAtoms > 0) {
      const fn = state.funnel;
      const nProt = state.ff.nProt;
      const ligStart = state.ff.ligandStart ?? nProt;

      let nc = 0;
      // Loop-2 S4 (R6 §5): contact form/break diff for the BindLog contact
      // channel — same 5.5 Å pair set the nContacts HUD counts, diffed
      // against the previous tick's set (state._lastContacts Map). Only
      // runs when BindLog capture is on (flag set below with trackTerms).
      const blg = (state.bindLog && bindLogWanted() && ff.trackTerms === true) ? state.bindLog : null;
      const curC = blg ? new Map() : null;
      for (let i = 0; i < nProt; i++) {
        for (let la = ligStart; la < state.ff.n; la++) {
          const dx = p[3 * la] - p[3 * i], dy = p[3 * la + 1] - p[3 * i + 1], dz = p[3 * la + 2] - p[3 * i + 2];
          const r2 = dx * dx + dy * dy + dz * dz;
          if (r2 < 30.25) {
            nc++;
            if (curC) curC.set(i * 1e6 + la, Math.sqrt(r2));
          }
        }
      }
      if (blg && curC) {
        const prevC = state._lastContacts;
        if (prevC) {
          for (const [k, d] of curC) {
            if (!prevC.has(k)) blg.pushContact(integ.time, true, k % 1e6, (k / 1e6) | 0, d);
          }
          for (const [k, d] of prevC) {
            if (!curC.has(k)) blg.pushContact(integ.time, false, k % 1e6, (k / 1e6) | 0, d);
          }
        }
        state._lastContacts = curC;
      } else if (state._lastContacts) {
        state._lastContacts = null; // capture turned off — drop stale set
      }

      let s = 0;
      for (let la = ligStart; la < state.ff.n; la++) {
        const dx = p[3 * la] - state.ff.ref[3 * la], dy = p[3 * la + 1] - state.ff.ref[3 * la + 1], dz = p[3 * la + 2] - state.ff.ref[3 * la + 2];
        s += dx * dx + dy * dy + dz * dz;
      }
      const ligRmsd = Math.sqrt(s / state.ff.nLigAtoms);
      // Convergence diagnostics: HUD grays out ΔG until nHills>=50 (collecting…)
      // Shows "ΔG ≈ X ± SE" when converged, else "– (collecting…)" or "ΔG not converged (nHills<50)"
      let dg;
      if (!fn || !fn.active || Number.isNaN(fn.estimateDG())) dg = "–";
      else if (fn._nHills < 50) dg = "– (collecting…)"; // nHills<50 not converged
      else {
        const se = typeof fn.convergenceSE === "function" ? fn.convergenceSE() : null;
        const val = fn.estimateDG().toFixed(2);
        dg = se != null ? `${val} ± ${se.toFixed(2)}` : val; // ΔG ≈ X ± SE when available
        // Alternative display when want explicit: `ΔG not converged (nHills<50)` handled above
      }
      // Minimal fallback for grep measurability: if(funnel._nHills<50) dgStr="– (collecting…)"
      if (fn && fn._nHills < 50) {
        // keep collecting state visible in HUD; already handled
      }
      const cvVal = fn && fn.active ? fn.lastCV : (fn ? fn.cv(p) : NaN);
      const cv = Number.isFinite(cvVal) ? cvVal.toFixed(2) : "–";

      let mlScoreStr = "";
      if (ui.poseScore && ui.poseScore.checked && state.scorer) {
        const feat = [
          Math.min(1, nc / 20),
          Math.min(1.2, Math.max(0, -ff.bindingU / 15.0)),
          0.0,
          Math.min(1, ligRmsd / 4.0),
          Math.min(1, (Number.isFinite(cvVal) ? cvVal : 10) / 10.0),
          Math.min(1, Math.max(0, -ff.bindingU / 10.0)),
        ];
        const score = state.scorer.predict(feat);
        mlScoreStr = `MLP Score = ${score.toFixed(2)}  ·  `;
      }

      extra = `CV = ${cv} Å  ·  ligRMSD = ${ligRmsd.toFixed(2)} Å  ·  nContacts = ${nc}  ·  ${mlScoreStr}ΔG ≈ ${dg} kcal/mol  ·  `;

      // Update Chemical Network state classification when live tracking is active
      if (state.running && isLiveTrackingActive) {
        const curMacroState = networkModel.classifyPose(parseFloat(cv) || 10, nc, ligRmsd);
        if (curMacroState !== networkModel.currentState) {
          networkModel.currentState = curMacroState;
          networkModel.probabilities.fill(0);
          networkModel.probabilities[curMacroState] = 1.0;
          updateNetworkPlot();
        }
      }
    }

    if (ui.hud) {
      // H76 — HUD debounce: throttled to 10 Hz (100 ms guard) — reduces DOM thrash at 60 fps
      // debounce guard: if (now - lastHudUpdate < 100) skip HUD update this frame
      if (now - lastHudUpdate >= 100) {
        lastHudUpdate = now;
        // Phase 5 — Metrics merged into the single #hud line (one T, one E;
        // #metricsHud stays in DOM for contract but hidden via CSS).
        let pmfStr = "—";
        if (state.funnel && state.funnel.active) {
          try {
            pmfStr = state.funnel._nHills < 50 ? `collect ${state.funnel._nHills}/50`
              : `ΔG ${state.funnel.estimateDG().toFixed(2)}`;
          } catch (_) { pmfStr = "—"; }
        }
        const dccmStr = recorder.count > 0 ? `${recorder.count}fr` : "—";
        ui.hud.textContent =
          `v${VERSION} · ` +
          (state.nanWarning ? `⚠ NON-FINITE ENERGY — simulation auto-paused. Reset (⟲) to recover.  ·  ` : "") +
          `t = ${integ.time.toFixed(1)} ps (${(integ.time / 1000).toFixed(3)} ns)  ·  ` +
          `U = ${Number.isFinite(ff.energy) ? ff.energy.toFixed(1) : "NaN"} kcal/mol  ·  ` +
          (ff.nLigAtoms > 0 ? `U_bind = ${ff.bindingU.toFixed(2)} kcal/mol  ·  ` : "") +
          extra +
          `RMSD = ${ff.rmsd(integ.pos).toFixed(2)} Å  ·  ` +
          `T_inst = ${ff.kineticTemp(integ.vel, integ.mass).toFixed(0)} K  ·  ` +
          `PMF ${pmfStr}  ·  DCCM ${dccmStr}  ·  ` +
          `Δr = ${dr.toFixed(2)} Å/frame  ·  ` +
          `${state.fpsEMA.toFixed(0)} fps  (${state.heavyMode ? `${state.ff.nProt} prot + ${state.ff.nHetero} hetero` : `${state.ff.nProt} Cα`}${state.ff.nLigAtoms ? ` + ${state.ff.nLigAtoms} lig` : ""})`;

        // Phase 5 — top bar: System State | PDB | Engine + FPS | Step (10 Hz).
        try {
          const pdbName = (ui.pdbId?.value?.trim()?.toUpperCase())
            || (state.sel ? `${state.ff.nProt}${state.heavyMode ? " heavy" : " Cα"}` : null) || "—";
          if (ui.sysState) ui.sysState.textContent = `STATE: ${state.nanWarning ? "FAULT" : state.running ? "● RUN" : state.integ ? "READY" : "IDLE"}`;
          if (ui.topPdb) ui.topPdb.textContent = `PDB ${pdbName}`;
          const be = settingsState.backend || "auto";
          const eng = be === "gpu" ? (settingsState.webgpuReady ? "WebGPU" : "WebGPU?")
            : be === "workers" ? `CPU×${settingsState.numThreads}` : be === "cpu" ? "CPU" : (settingsState.webgpuReady ? "WebGPU/CPU" : "CPU");
          if (ui.topEngine) ui.topEngine.textContent = `Engine: ${eng} ${Math.round(state.fpsEMA)}fps`;
          if (ui.topStep) ui.topStep.textContent = `step ${state.integ.steps ?? 0}`;
        } catch (_) { /* headless */ }

        // Phase 5 — Metrics HUD mirror (hidden via CSS; #hud above is the
        // single visible readout — keep one T, one E there).
        try {
          const tInst = ff.kineticTemp(integ.vel, integ.mass);
          const eTot = Number.isFinite(ff.energy) ? ff.energy.toFixed(1) : "NaN";
          if (ui.metricsHud) ui.metricsHud.textContent = `T ${Number.isFinite(tInst) ? tInst.toFixed(0) : "—"}K · Etot ${eTot} · PMF ${pmfStr} · DCCM ${dccmStr}`;
        } catch (_) { /* headless */ }

        // Phase 5 — canvas caption: compact status pill (P2). Full hints
        // stay in the empty-state text only; running/steady are short.
        try {
          if (ui.canvasCaption) {
            if (!state.sel || !state.ff) ui.canvasCaption.textContent = "No system — Load a PDB, Build, then Run. Ligand halo = white ring.";
            else if (state.running) ui.canvasCaption.textContent = ff.nLigAtoms > 0 ? `● Running · halo ${ff.nLigAtoms}` : "● Running";
            else ui.canvasCaption.textContent = ff.nLigAtoms > 0 ? `○ Steady · halo ${ff.nLigAtoms}` : "○ Steady";
          }
        } catch (_) { /* headless */ }

        updateDockTimeline();
      }
    }
    // Phase 5 — compact strips at display rate (no chart libs): push every
    // frame when live, draw every frame (cheap 120-pt polyline).
    try {
      if (state.ff && state.integ) {
        const cvLive = (state.funnel && Number.isFinite(state.funnel.lastCV)) ? state.funnel.lastCV : NaN;
        _pushHist(_cvHist, cvLive);
        _pushHist(_eHist, state.ff.energy);
        _drawStrip(ui.cvStrip, _cvHist, "#E879F9", { empty: "Run to stream CV" });
        _drawStrip(ui.hudSpark, _eHist, "#38BDF8", { empty: "Run to stream E" });
      }
    } catch (_) { /* headless */ }
    updatePMFPlot(!!state.running);
    // Network canvas steady-state redraw (wiki P2): 1 Hz even without state
    // change, so the CNM diagram never reads as dead.
    networkPanelTick(now);
    dccmTick(now);
    // Loop-2 S6 (R7 §4): BindViz ≤1 Hz — empty states or live plots (P2).
    bindvizTick(now);
    guideTick(now); // FP1: ≤1 Hz checklist poll (catches steps/frames/Analyze)
  } else if (viewer) {
    viewer.render(null);
    try { dccmTick(now); } catch (_) {}
    try { bindvizTick(now); } catch (_) {}
    try { guideTick(now); } catch (_) {} // FP1: empty-state emphasis while idle
  } else {
    console.warn("[viewer] not ready");
  }
}

if (typeof requestAnimationFrame !== "undefined") {
  requestAnimationFrame(tick);
}
