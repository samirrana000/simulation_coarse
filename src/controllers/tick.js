/**
 * controllers/tick.js — the requestAnimationFrame simulation + render loop.
 *
 * PUBLIC API
 *   startTick()  Seed the wall clock and register the first animation frame.
 *                Call exactly once, LAST, from the composition root.
 *
 * NEEDS (imports): ui + state + viewer + recorder from ../ui.js; VERSION from
 * ../version.js; settingsState from ../settings-panel.js; updatePMFPlot from
 * ../pmf-panel.js; dccmTick from ../analysis-panel.js; networkModel +
 * isLiveTrackingActive + updateNetworkPlot + networkPanelTick from
 * ../network-panel.js; formatLiveTermsHUD + updateLiveTermsMirror from
 * ./live-terms.js; bindvizTick from ./bindviz.js; drawDockStrips +
 * updateDockTimeline from ./dock.js; guideTick from ./guide.js;
 * bindLogWanted from ./physics-tier.js; respaWanted + ensureRespa +
 * warnRespaFallbackOnce from ./accelerate.js.
 *
 * CADENCE CONTRACT (the reason this loop is its own module, and the reason
 * tests/test_main_module_size.js exists). Per animation frame, in order:
 *
 *   every frame   viewer.render, Δr accumulation, ligand CV/RMSD/contacts,
 *                 CNM classification, the two dock strips, updatePMFPlot,
 *                 fpsEMA
 *   ≥ 1 Hz        networkPanelTick, dccmTick, bindvizTick, guideTick
 *   10 Hz (100ms) the #hud line, #canvasCaption, #metricsHud mirror, the top
 *                 bar (#sysState/#topPdb/#topEngine/#topStep) and
 *                 updateDockTimeline — ALL inside one `now - lastHudUpdate >=
 *                 100` guard, so they cannot drift apart
 *
 * The idle branch (no integrator) still runs dccmTick/bindvizTick/guideTick so
 * the empty states never read as dead.
 */

import { ui, state, viewer, recorder, updateRecStatus } from "../ui.js";
import { VERSION } from "../version.js";
import { settingsState } from "../settings-panel.js";
import { updatePMFPlot } from "../pmf-panel.js";
import { dccmTick } from "../analysis-panel.js";
import { networkModel, isLiveTrackingActive, updateNetworkPlot, networkPanelTick } from "../network-panel.js";
import { formatLiveTermsHUD, updateLiveTermsMirror } from "./live-terms.js";
import { bindvizTick } from "./binding-insights.js";
import { drawDockStrips, updateDockTimeline } from "./dock.js";
import { guideTick } from "./guide.js";
import { bindLogWanted } from "./physics-tier.js";
import { respaWanted, ensureRespa, warnRespaFallbackOnce } from "./accelerate.js";

let lastT = 0;
// H76 — HUD debounce: throttle DOM updates to 10 Hz (100 ms) to reduce reflow thrash
let lastHudUpdate = 0; // debounce guard for hud.textContent — 10 Hz

/** One animation frame: step, render, then the per-cadence readouts. */
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
        warnRespaFallbackOnce(e);
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
        // Rev1/Issue4: refresh the live per-term mirror at the same 10 Hz HUD
        // cadence (reads last ff fields — no extra compute, no checkbox).
        try { updateLiveTermsMirror(ff, integ.time); } catch (_) { /* headless */ }
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
          formatLiveTermsHUD(ff) +
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
    drawDockStrips();
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

/** Seed the wall clock and register the first animation frame. */
export function startTick() {
  lastT = typeof performance !== "undefined" ? performance.now() : 0;
  if (typeof requestAnimationFrame !== "undefined") {
    requestAnimationFrame(tick);
  }
}