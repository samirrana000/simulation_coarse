/**
 * controllers/system-build.js — selection -> force field -> integrator.
 *
 * PUBLIC API
 *   buildSystem()          Build the system the UI currently describes:
 *                          CG (synchronous) or heavy (chunked + cancellable).
 *                          Returns void; all failures render to the reused
 *                          summary surfaces instead of throwing.
 *   initSystemBuild(seams) Wire the build surface once: the heavy Cancel-build
 *                          button, the ligand panel, #buildBtn / #modelMode /
 *                          #includeLig and the four hetero quick-actions.
 *                          `seams` is {physicsLevelSpec, ffParams} from the
 *                          composition root (src/main.js) — see params.js for
 *                          why these are injected rather than imported.
 *   invalidateHeavyBuild() Cancel any in-flight chunked heavy build.
 *
 * NEEDS (imports): ui + state + viewer + recorder + updateSelSummary +
 * updateRecStatus from ../ui.js; ForceField from ../forcefield.js;
 * HeavyForceField + selectHeavy + appendHeavyLigands +
 * buildTopologyChunked from ../heavy.js; Funnel from ../funnel.js;
 * selectSystem + parseLigands from ../pdb.js; classifyInputError +
 * formatInputError + checkSystemSize from ../input_errors.js;
 * invalidateThermo + drawDccmEmpty + refreshThermoLigPicker from
 * ../analysis-panel.js; HEAVY_TOPO_CHUNK_ROWS + setHeavyButtons +
 * setHeavyCaption + heavyTopoCaption from ../heavy_progress.js;
 * assignProtonationStates + applyProtonationStates from ../chem/protonation.js;
 * workerPool from ../settings-panel.js; initLigandPanel from
 * ../ligand-panel.js; applyMLToFF from ../ml-tier.js; BindLog from
 * ../capture/bindlog.js; bindLogWanted from ./physics-tier.js;
 * applyLiveTrackTerms from ./live-terms.js; drawBindviz from ./bindviz.js;
 * updateDockTimeline + resetDockHistory from ./dock.js; updateGuide from
 * ./guide.js; onParamChange from ./params.js.
 *
 * ZERO PHYSICS ADDED: the CG path is synchronous and unchanged; the heavy path
 * only slices the EXISTING buildTopologyChunked so paint survives. Both paths
 * end in the same finishBuildCommon().
 */

/* ------------------------------------------------------------------ */
/*  System construction                                                 */
/* ------------------------------------------------------------------ */
import { ui, state, viewer, recorder, updateSelSummary, updateRecStatus } from "../ui.js";
import { ForceField } from "../forcefield.js";
import { HeavyForceField, selectHeavy, appendHeavyLigands, buildTopologyChunked } from "../heavy.js";
import { Funnel } from "../funnel.js";
import { LangevinIntegrator } from "../integrator.js";
import { selectSystem, parseLigands } from "../pdb.js";
import { classifyInputError, formatInputError, checkSystemSize } from "../input_errors.js";
import { invalidateThermo, drawDccmEmpty, refreshThermoLigPicker } from "../analysis-panel.js";
import { HEAVY_TOPO_CHUNK_ROWS, setHeavyButtons, setHeavyCaption, heavyTopoCaption } from "../heavy_progress.js";
import { assignProtonationStates, applyProtonationStates } from "../chem/protonation.js";
import { workerPool } from "../settings-panel.js";
import { initLigandPanel } from "../ligand-panel.js";
import { applyMLToFF } from "../ml-tier.js";
import { BindLog } from "../capture/bindlog.js";
import { bindLogWanted } from "./physics-tier.js";
import { applyLiveTrackTerms } from "./live-terms.js";
import { drawBindviz } from "./binding-insights.js";
import { updateDockTimeline, resetDockHistory } from "./dock.js";
import { updateGuide } from "./guide.js";
import { onParamChange } from "./param-binding.js";

/** Cross-cutting seams injected by initSystemBuild (see module header). */
let _seams = null;

/** "#A,#b" -> ["A","B"] uppercased; null when the box is empty/absent. */
function parseParamChainIds() {
  if (!ui.chainsInput) return null;
  const t = ui.chainsInput.value.trim();
  if (!t) return null;
  return t.split(",").map((s) => s.trim().toUpperCase()).filter(Boolean)
    .map((s) => (s === "_" ? "_" : s));
}

/** Default on/off for every hetero group, honouring user overrides. */
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

/** Repaint the hetero checkbox list (heavy mode only). */
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

/** Build the system the UI currently describes (CG sync / heavy chunked). */
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
  // HeavyForceField mirrors physicsLevel/charges/hbMode queryably (Rev3/Issue1,
  // kernels unchanged) and consumes par.weak for the S3 terms (ignored by CG).
  const physLvl = _seams.physicsLevelSpec();
  const par = _seams.ffParams(physLvl);

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
  // Loop-2 S4+S7 + Rev1/Issue4: per-term accumulators stay on while a ligand
  // is present (live HUD/bindviz mirror) or while capturing (manual checkbox
  // OR L2 tier). U/forces bit-identical either way; only the event capture
  // below still requires bindLogWanted().
  applyLiveTrackTerms(state.ff);

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
  // Rev1/Issue4: fresh Build also clears the live per-term mirror (no stale stripes).
  state._liveTerms = null;
  state._liveTermsHist = [];
  // Phase 5 — reset dock + strips + DCCM empty state on rebuild.
  try {
    resetDockHistory();
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

/** Wire the build surface once (see module header for the seam). */
export function initSystemBuild(seams) {
  _seams = seams;

  if (ui.heavyCancelBtn) ui.heavyCancelBtn.addEventListener("click", () => {
    invalidateHeavyBuild(); // bump generation → stale topo slices/continuations return early
    setHeavyBuilding(false); // belt-and-braces idle state (invalidateHeavyBuild already resets)
    setHeavyCaptionText("cancelled — heavy build cancelled by user. Click Build System to retry.");
  });

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
}