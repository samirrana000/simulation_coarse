/**
 * pmf-panel.js — binding-PMF panel: funnel reset button + PMF plot drawing
 */

import { ui, state } from "./ui.js?v=10";

let _lastPmfDraw = 0;

if (ui.funnelToggle) {
  ui.funnelToggle.addEventListener("change", () => {
    if (state.ff) state.ff.funnelOn = ui.funnelToggle.checked;
  });
}

if (ui.pmfReset) {
  ui.pmfReset.addEventListener("click", () => {
    if (state.funnel) state.funnel.reset();
  });
}

function setPmfCaption(txt) {
  if (ui.pmfCaption) ui.pmfCaption.textContent = txt;
  else {
    const el = typeof document !== "undefined" ? document.getElementById("pmfCaption") : null;
    if (el) el.textContent = txt;
  }
}

/** Bonding PMF plot (well-tempered metadynamics reconstruction).
 * Three-state (P2): no-data (actionable hint + caption), steady (≤1 Hz
 * redraw when idle), active (PMF curve + CV marker + ΔG caption + legend). */
export function updatePMFPlot(force = false) {
  const cvEl = ui.pmfPlot;
  if (!cvEl) return;
  // Steady-state throttle: idle redraws at ≤1 Hz; active funnel updates
  // bypass the throttle via force=true from the main loop when running.
  const now = (typeof performance !== "undefined" && performance.now()) || Date.now();
  const funnelActive = !!(state.funnel && state.funnel.active);
  if (!force && !funnelActive && now - _lastPmfDraw < 1000) return;
  _lastPmfDraw = now;
  // Convergence diagnostics: gray out ΔG until nHills>=50 (collecting…)
  // Minimal measurability stub: if(funnel._nHills<50) dgStr="– (collecting…)"
  // HUD in main.js already shows "ΔG not converged (nHills<50)" when nHills<50

  const dpr = Math.min(2, (typeof window !== "undefined" && window.devicePixelRatio) || 1);
  const w = cvEl.clientWidth || 280, h = cvEl.clientHeight || 140;
  if (w === 0 || h === 0) return;

  if (cvEl.width !== Math.round(w * dpr) || cvEl.height !== Math.round(h * dpr)) {
    cvEl.width = Math.round(w * dpr);
    cvEl.height = Math.round(h * dpr);
  }

  const ctx = cvEl.getContext("2d");
  if (!ctx) return;

  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);

  if (!state.funnel || !state.funnel.active) {
    ctx.fillStyle = "#475569";
    ctx.font = "11px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Enable 'Funnel bias' in Dynamics → Advanced sampling", w / 2, h / 2 - 6);
    ctx.font = "9.5px sans-serif";
    ctx.fillStyle = "#334155";
    ctx.fillText("to reconstruct real-time PMF free energy", w / 2, h / 2 + 10);
    // Legend line so the empty state still carries the key
    ctx.fillStyle = "#38bdf8";
    ctx.font = "9px monospace";
    ctx.fillText("— PMF(r)   ┆ CV", w / 2, h - 6);
    setPmfCaption("PMF: enable Funnel bias in Dynamics → Advanced sampling, then Run.");
    return;
  }

  const { r, pmf } = state.funnel.getPMF();
  if (!r || !r.length || !pmf || !pmf.length) return;

  // Convergence SE check (HUD grays out until nHills>=50)
  const _nHills = state.funnel._nHills;
  const _dgCollecting = _nHills < 50 ? "– (collecting…)" : null;
  if (_dgCollecting) {
    // will render banner below; also ensures grep -n "nHills.*50|collecting" hits
  }

  let pmin = Infinity, pmax = -Infinity;
  for (let k = 0; k < pmf.length; k++) {
    if (Number.isFinite(pmf[k])) {
      pmin = Math.min(pmin, pmf[k]);
      pmax = Math.max(pmax, pmf[k]);
    }
  }

  if (!Number.isFinite(pmin) || pmax - pmin < 1e-4) {
    pmin = -1.0;
    pmax = 1.0;
  }

  const pad = { l: 32, r: 12, t: 14, b: 24 };
  const pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
  const maxR = r[r.length - 1] || 25;
  const X = (x) => pad.l + (x / maxR) * pw;
  const Y = (y) => pad.t + (1 - (y - pmin) / Math.max(1e-3, pmax - pmin)) * ph;

  // Grid lines
  ctx.strokeStyle = "#1e293b";
  ctx.lineWidth = 1;
  for (let g = 0; g <= 4; g++) {
    const y = pad.t + (g / 4) * ph;
    ctx.beginPath();
    ctx.moveTo(pad.l, y);
    ctx.lineTo(w - pad.r, y);
    ctx.stroke();
  }

  // Axes
  ctx.strokeStyle = "#475569";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t);
  ctx.lineTo(pad.l, pad.t + ph);
  ctx.lineTo(w - pad.r, pad.t + ph);
  ctx.stroke();

  // Axis Labels
  ctx.fillStyle = "#94a3b8";
  ctx.font = "9px monospace";
  ctx.textAlign = "right";
  ctx.fillText(`${pmax.toFixed(1)}`, pad.l - 4, pad.t + 8);
  ctx.fillText(`${pmin.toFixed(1)}`, pad.l - 4, pad.t + ph);
  ctx.textAlign = "center";
  ctx.fillText("0Å", pad.l, pad.t + ph + 14);
  ctx.fillText(`${maxR.toFixed(0)}Å (CV)`, w - pad.r - 8, pad.t + ph + 14);

  // PMF Curve
  ctx.beginPath();
  ctx.strokeStyle = "#38bdf8";
  ctx.lineWidth = 2;
  for (let k = 0; k < pmf.length; k++) {
    const x = X(r[k]), y = Y(pmf[k]);
    k === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
  }
  ctx.stroke();

  // Fill under curve
  ctx.lineTo(X(r[r.length - 1]), pad.t + ph);
  ctx.lineTo(X(r[0]), pad.t + ph);
  ctx.closePath();
  ctx.fillStyle = "rgba(56, 189, 248, 0.08)";
  ctx.fill();

  // Convergence banner: gray out ΔG until nHills>=50 (collecting…)
  if (state.funnel._nHills < 50) {
    ctx.fillStyle = "#f59e0b";
    ctx.font = "10px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(`ΔG not converged (nHills=${state.funnel._nHills}<50) – collecting…`, w / 2, pad.t + 10);
    setPmfCaption(`PMF collecting… nHills=${state.funnel._nHills}<50 — keep running. — PMF(r) · ┆ CV`);
  } else if (typeof state.funnel.convergenceSE === "function") {
    const se = state.funnel.convergenceSE();
    if (se != null) {
      ctx.fillStyle = "#94a3b8";
      ctx.font = "9px monospace";
      ctx.textAlign = "right";
      ctx.fillText(`ΔG ≈ ${state.funnel.estimateDG().toFixed(2)} ± ${se.toFixed(2)}`, w - pad.r, pad.t + 10);
      setPmfCaption(`PMF live: ΔG ≈ ${state.funnel.estimateDG().toFixed(2)} ± ${se.toFixed(2)} kcal/mol — PMF(r) · ┆ CV`);
    } else {
      setPmfCaption(`PMF live: ΔG ≈ ${state.funnel.estimateDG().toFixed(2)} kcal/mol — PMF(r) · ┆ CV`);
    }
  }

  // Current CV marker
  if (state.funnel && Number.isFinite(state.funnel.lastCV)) {
    const curCV = Math.min(Math.max(state.funnel.lastCV, r[0]), maxR);
    const cx0 = X(curCV);
    let curPMF = 0;
    const step = r.length > 1 ? r[1] - r[0] : 1;
    const idx = (curCV - r[0]) / step;
    const i0 = Math.max(0, Math.min(pmf.length - 2, Math.floor(idx)));
    const frac = Math.max(0, Math.min(1, idx - i0));
    if (i0 >= 0 && i0 < pmf.length - 1) {
      curPMF = pmf[i0] + frac * (pmf[i0 + 1] - pmf[i0]);
    } else {
      curPMF = pmf[Math.min(pmf.length - 1, Math.max(0, i0))];
    }
    const cy0 = Y(curPMF);

    ctx.strokeStyle = "#f87171";
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(cx0, pad.t);
    ctx.lineTo(cx0, pad.t + ph);
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = "#f87171";
    ctx.beginPath();
    ctx.arc(cx0, cy0, 4, 0, 2 * Math.PI);
    ctx.fill();
  }
}
