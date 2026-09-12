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
import { downloadText } from "./recorder.js?v=10";
import { PoseScorer } from "./scorer.js?v=10";
import { ui, state, viewer, recorder, initParamReadouts, updateSelSummary, updateRecStatus } from "./ui.js?v=10";
import "./analysis-panel.js?v=10"; // side-effect: Analyze + Phase-4 workflow buttons
import { dccmTick, drawDccmEmpty } from "./analysis-panel.js?v=10"; // steady ≤1 Hz DCCM empty redraw (P2)
import { applyMLToFF } from "./ml-tier.js?v=10";
import { updatePMFPlot } from "./pmf-panel.js?v=10";
import { initLigandPanel, updateMol2PlaceButton } from "./ligand-panel.js?v=10";
import { parseHeavy, HeavyForceField, selectHeavy, appendHeavyLigands } from "./heavy.js?v=10";
import { assignProtonationStates, applyProtonationStates } from "./chem/protonation.js?v=10";
import { initSettingsModal, settingsState, workerPool, gpuAccelerator } from "./settings-panel.js?v=10";
import { RESPAStepper, splitForceField } from "./physics/integrators/respa.js?v=10";
import { initNetworkPanel, updateNetworkPlot, networkPanelTick, networkModel, isLiveTrackingActive } from "./network-panel.js?v=10";
import { BindLog } from "./capture/bindlog.js?v=10";
import { renderInteractionTimeline, renderEnergyDecomposition, renderPmfFormation } from "./capture/bindviz.js?v=10";

// Initialize UI modals & panels
initSettingsModal();
initNetworkPanel();

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
    ui.scrubLabel.textContent = `${recorder.count > 0 ? cur + 1 : 0} / ${recorder.count} frames`;
  }
}
if (ui.scrub) {
  ui.scrub.addEventListener("input", () => {
    const i = Number(ui.scrub.value);
    const fr = recorder.getFrame(i);
    if (ui.scrubLabel) ui.scrubLabel.textContent = `${recorder.count > 0 ? i + 1 : 0} / ${recorder.count} frames`;
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
    ui.bindvizCaption.textContent = hasData
      ? `BindViz · ${bl.nEvents} events · ${bl.nFrames} frames (1 Hz live).`
      : "Enable BindLog capture in Recording, run, then insights render live (1 Hz).";
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
    throw new Error(caErr?.message || heavyErr?.message || "No usable atoms found in PDB file.");
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
      if (ui.structSummary) ui.structSummary.textContent = "⚠ " + err.message;
      if (ui.hud) ui.hud.textContent = "⚠ " + err.message;
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

if (ui.fileInput) {
  ui.fileInput.addEventListener("change", async () => {
    const f = ui.fileInput.files[0];
    if (!f) return;
    try {
      await loadStructure(await f.text(), f.name);
    } catch (err) {
      if (ui.structSummary) ui.structSummary.textContent = "⚠ " + err.message;
      if (ui.hud) ui.hud.textContent = "⚠ " + err.message;
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
      if (!mols.length) throw new Error(`${f.name}: no usable molecules in MOL2 file`);
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
        ui.mol2Info.textContent = "⚠ " + err.message;
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

export function buildSystem() {
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
    if (ui.selSummary) ui.selSummary.textContent = "⚠ " + err.message;
    if (ui.hud) ui.hud.textContent = "⚠ " + err.message;
    return;
  }

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
    state.ff = new HeavyForceField({ atoms: state.sel.atoms }, par, []);
    workerPool.initSystem(state.ff);
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
  physicsLevelSpec(); // persist to settingsState.physicsLevel
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

    if (e.code === "Space") {
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
      const idx = parseInt(e.key, 10) - 1;
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
  });
}
if (ui.recStopBtn) {
  ui.recStopBtn.addEventListener("click", () => {
    recorder.stop();
    ui.recBtn.classList.remove("rec-on");
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
      const text = recorder.buildFile(fmt, state.sel.beads, prov);
      const name = `cg_traj_${recorder.count}frames.${fmt}`;
      downloadText(text, name);
    } catch (err) {
      if (ui.recStatus) ui.recStatus.textContent = "⚠ " + err.message;
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
        _drawStrip(ui.cvStrip, _cvHist, "#E879F9", { empty: "CV —" });
        _drawStrip(ui.hudSpark, _eHist, "#38BDF8", { empty: "E —" });
      }
    } catch (_) { /* headless */ }
    updatePMFPlot(!!state.running);
    // Network canvas steady-state redraw (wiki P2): 1 Hz even without state
    // change, so the CNM diagram never reads as dead.
    networkPanelTick(now);
    dccmTick(now);
    // Loop-2 S6 (R7 §4): BindViz ≤1 Hz — empty states or live plots (P2).
    bindvizTick(now);
  } else if (viewer) {
    viewer.render(null);
    try { dccmTick(now); } catch (_) {}
    try { bindvizTick(now); } catch (_) {}
  } else {
    console.warn("[viewer] not ready");
  }
}

if (typeof requestAnimationFrame !== "undefined") {
  requestAnimationFrame(tick);
}
