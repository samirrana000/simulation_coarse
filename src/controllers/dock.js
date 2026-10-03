/**
 * controllers/dock.js — Phase-5 dock readouts: compact strips + scrub slider.
 *
 * PUBLIC API
 *   initDock()          Wire the #scrub "seek recorded frame" listener (once).
 *   updateDockTimeline() Refresh #scrub bounds/value and the #scrubLabel
 *                       "N / M frames" caption. Callers: system build (reset),
 *                       and the RAF tick inside the 10 Hz HUD block.
 *   resetDockHistory()  Empty both ring buffers (a fresh Build must not keep
 *                       pre-Build samples on the strips).
 *   drawDockStrips()    Per-frame: push the live CV + total energy samples and
 *                       repaint #cvStrip / #hudSpark. Called from the RAF tick
 *                       EVERY frame (a 120-point polyline is cheap; this is
 *                       deliberately NOT throttled).
 *
 * NEEDS (imports): ui + state + recorder + viewer from ../ui.js.
 *
 * Phase 5 — dock sparklines (compact Canvas strips, no chart libs).
 * Fixed-length ring buffers for CV (Å) + total energy (kcal/mol).
 */

import { ui, state, viewer, recorder } from "../ui.js";
import { ignore } from "../errors.js";

// Phase 5 — dock sparklines (compact Canvas strips, no chart libs).
// Fixed-length ring buffers for CV (Å) + total energy (kcal/mol).
const _cvHist = [];
const _eHist = [];
const _HIST_N = 120;

/** Append one sample to a fixed-length ring buffer (non-finite dropped). */
function _pushHist(arr, v) {
  if (!Number.isFinite(v)) return;
  arr.push(v);
  if (arr.length > _HIST_N) arr.shift();
}

/** Paint one ring buffer as a normalised polyline with a live head dot. */
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
} catch (e) { ignore(e, "sparkline draw@drawDockStrips", "canvas absent headless; the scrub bounds are computed independently"); }
}

/** Refresh the scrub bounds + "N / M frames" caption. */
export function updateDockTimeline() {
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

/** Empty both ring buffers (a fresh Build starts with clean strips). */
export function resetDockHistory() {
  _cvHist.length = 0;
  _eHist.length = 0;
}

/**
 * Per-frame strip update: push live CV + energy, repaint both strips.
 * Runs at display rate (not throttled) — see the module header.
 */
export function drawDockStrips() {
  try {
    if (state.ff && state.integ) {
      const cvLive = (state.funnel && Number.isFinite(state.funnel.lastCV)) ? state.funnel.lastCV : NaN;
      _pushHist(_cvHist, cvLive);
      _pushHist(_eHist, state.ff.energy);
      _drawStrip(ui.cvStrip, _cvHist, "#E879F9", { empty: "Run to stream CV" });
      _drawStrip(ui.hudSpark, _eHist, "#38BDF8", { empty: "Run to stream E" });
    }
} catch (e) { ignore(e, "sparkline draw@updateDockTimeline", "canvas absent headless; the scrub label is set independently"); }
}

/** Wire the trajectory scrub slider (H78): seeks recorded frames. */
export function initDock() {
  if (!ui.scrub) return;
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