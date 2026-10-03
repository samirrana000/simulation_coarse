/**
 * controllers/binding-insights.js — Loop-2 S6: the three Binding-Insights canvases.
 *
 * PUBLIC API
 *   initBindviz()    Paint the three empty states once (no system yet). Called
 *                    once by the composition root, and again after every Build
 *                    (see finishBuildCommon in ./system-build.js).
 *   drawBindviz()    Full redraw of #bindvizTimeline / #bindvizEnergy /
 *                    #bindvizPmf plus the #bindvizCaption.
 *   bindvizTick(now) Steady ≤1 Hz wrapper around drawBindviz — no-data canvases
 *                    redraw actionable text, live data redraws plots.
 *
 * NEEDS (imports): ui + state from ../ui.js; the three renderers +
 * PMF_NOHILL_HINT from ../capture/bindviz.js; liveTermsBindLogView from
 * ./live-terms.js (the live per-term mirror, so the energy stripes are live
 * without BindLog capture).
 *
 * Pixel sizing follows the dccm/pmf devicePixelRatio pattern
 * (canvas.width = clientWidth x dpr, setTransform(dpr)) so CSS-sized canvases
 * render crisp; renderers read the CSS-pixel size back via canvas._dpr.
 * ≤1 Hz steady (three-state P2).
 */

import { ui, state } from "../ui.js";
import { renderInteractionTimeline, renderEnergyDecomposition, renderPmfFormation, PMF_NOHILL_HINT } from "../capture/bindviz.js";
import { liveTermsBindLogView } from "./live-terms.js";
import { ignore } from "../errors.js";

let _lastBindvizDraw = 0;

/** Size a canvas to its CSS box x devicePixelRatio and return its 2D context. */
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
export function drawBindviz() {
  const bl = state.bindLog;
  const hasData = !!(bl && bl.nFrames > 0 && bl.nEvents > 0);
  // Rev1/Issue4: energy canvas falls back to the live per-term mirror when the
  // BindLog carries no energy events — normal Runs render per-term stripes
  // without requiring the BindLog checkbox. Timeline/PMF still need capture
  // (contacts/hills), so they keep the hasData gate.
  let energyBl = hasData ? bl : null;
  let energyLive = false;
  if (energyBl) {
    let hasEnergy = false;
    try {
      for (let i = 0; i < energyBl.nEvents; i++) if (energyBl.evType[i] === 0) { hasEnergy = true; break; }
    } catch (_) { hasEnergy = true; }
    if (!hasEnergy) energyBl = null;
  }
  if (!energyBl) {
    try {
      const view = liveTermsBindLogView();
      if (view) { energyBl = view; energyLive = true; }
    } catch (_) { energyBl = null; }
  }
  if (ui.bindvizTimeline && bindvizFit(ui.bindvizTimeline)) renderInteractionTimeline(ui.bindvizTimeline, hasData ? bl : null);
  if (ui.bindvizEnergy && bindvizFit(ui.bindvizEnergy)) renderEnergyDecomposition(ui.bindvizEnergy, energyBl);
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
    ui.bindvizCaption.textContent = !hasData && !energyLive
      ? "Enable BindLog capture in Recording, run, then insights render live (1 Hz). Energy stripe is live without capture once running."
      : (!hasData && energyLive
        ? `Live per-term energies (no BindLog — ${state._liveTermsHist?.length ?? 0} samples, 1 Hz). Timeline/PMF need BindLog capture.`
        : (nHills === 0
          ? `BindViz · ${bl.nEvents} events · ${bl.nFrames} frames (1 Hz live). PMF: ${PMF_NOHILL_HINT}`
          : `BindViz · ${bl.nEvents} events · ${bl.nFrames} frames (1 Hz live).`));
  }
}

/** Steady ≤1 Hz tick — no-data canvases redraw actionable text, live data redraws plots. */
export function bindvizTick(now) {
  if (now - _lastBindvizDraw < 1000) return;
  _lastBindvizDraw = now;
  try { drawBindviz(); } catch (e) { ignore(e, "drawBindviz@initBindviz", "canvas is absent headless; updateSettingsUI-style UI wiring must still complete"); }
}

/** Paint the three empty states once at startup (no system yet). */
export function initBindviz() {
  try { drawBindviz(); } catch (e) { ignore(e, "drawBindviz@emptyStates", "canvas is absent headless; the caption text is set separately"); }
}