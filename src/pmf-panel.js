/**
 * pmf-panel.js — binding-PMF panel: funnel reset button + PMF plot drawing
 * (item 5 modularization: moved verbatim from main.js; only the imports
 * changed — `ui` / `state` now come from ui.js).
 */

import { ui, state } from "./ui.js?v=9";

// funnel toggle + PMF reset — these listeners belong with the PMF panel
ui.funnelToggle.addEventListener("change", () => {
  if (state.ff) state.ff.funnelOn = ui.funnelToggle.checked;
});
ui.pmfReset.addEventListener("click", () => { if (state.funnel) state.funnel.reset(); });

/** Bonding PMF plot (well-tempered metadynamics reconstruction), ~once/frame. */
export function updatePMFPlot() {
  const cvEl = ui.pmfPlot;
  if (!cvEl || !state.funnel || !state.funnel.active || cvEl.width === 0) return;
  const { r, pmf } = state.funnel.getPMF();
  if (!r.length) return;
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  const w = cvEl.clientWidth, h = cvEl.clientHeight;
  if (w === 0) return;
  if (cvEl.width !== Math.round(w*dpr) || cvEl.height !== Math.round(h*dpr)) { cvEl.width = Math.round(w*dpr); cvEl.height = Math.round(h*dpr); }
  const ctx = cvEl.getContext("2d");
  ctx.setTransform(dpr,0,0,dpr,0,0);
  ctx.clearRect(0,0,w,h);
  // bounds
  let pmin = Infinity, pmax = -Infinity;
  for (let k = 0; k < pmf.length; k++) { pmin = Math.min(pmin, pmf[k]); pmax = Math.max(pmax, pmf[k]); }
  if (!Number.isFinite(pmin) || pmax - pmin < 1e-6) return;
  const pad = { l: 26, r: 8, t: 10, b: 20 };
  const pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
  const X = (x) => pad.l + (x / r[r.length-1]) * pw;
  const Y = (y) => pad.t + (1 - (y - pmin) / (pmax - pmin)) * ph;
  // grid lines + axes
  ctx.strokeStyle = "#2a3040"; ctx.lineWidth = 1;
  for (let g = 0; g <= 4; g++) {
    const y = pad.t + (g/4)*ph;
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w-pad.r, y); ctx.stroke();
  }
  ctx.strokeStyle = "#8a93a8"; ctx.fillStyle = "#8a93a8"; ctx.font = "9px sans-serif";
  ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, pad.t+ph); ctx.lineTo(w-pad.r, pad.t+ph); ctx.stroke();
  ctx.fillText("ΔG", 4, pad.t+8); ctx.fillText("0", 4, pad.t+ph);
  // curve
  ctx.beginPath(); ctx.strokeStyle = "#ffb050"; ctx.lineWidth = 1.5;
  for (let k = 0; k < pmf.length; k++) { const x = X(r[k]), y = Y(pmf[k]); k === 0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y); }
  ctx.stroke();
  // current CV marker
  if (state.funnel && state.funnel.active && Number.isFinite(state.funnel.lastCV)) {
    const cx0 = X(Math.min(Math.max(state.funnel.lastCV, r[0]), r[r.length-1]));
    ctx.fillStyle = "#6cc4ff";
    ctx.beginPath(); ctx.arc(cx0, Y(Math.max(pmin, Math.min(pmax, 0))), 3, 0, 6.2832); ctx.fill();
  }
}
