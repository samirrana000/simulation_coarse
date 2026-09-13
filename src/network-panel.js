/**
 * network-panel.js — Chemical Network Model (CNM) live visualizer & kinetics controller.
 *
 * Canvas three-state rule (wiki P2): every live canvas renders
 *   1. no-data  — actionable empty-state message
 *   2. steady   — 1 Hz redraw even without state change (not event-only)
 *   3. active   — transition glow + caption + legend so it self-explains
 */

import { ChemicalNetworkModel } from "./physics/network.js?v=10";
import { viewer } from "./ui.js?v=10";

export const networkModel = new ChemicalNetworkModel({ temperature: 300 });
export let isLiveTrackingActive = true;

let _lastSteadyDraw = 0;

export function updateNetworkPlot(force = false) {
  if (typeof document === "undefined") return; // headless/Node: no DOM to update
  const canvas = document.getElementById("networkCanvas");
  if (!canvas) return;

  const dpr0 = (typeof window !== "undefined" && window.devicePixelRatio) || 1;
  const rect = canvas.getBoundingClientRect();
  if (rect.width > 0 && canvas.width !== Math.floor(rect.width * (Math.min(2, dpr0)))) {
    const dpr = Math.min(2, dpr0);
    canvas.width = Math.floor(rect.width * dpr);
    canvas.height = Math.floor(140 * dpr);
  }

  const ctx = canvas.getContext("2d");
  if (!ctx) return;

  const dpr = Math.min(2, (typeof window !== "undefined" && window.devicePixelRatio) || 1);
  const w = canvas.width / dpr;
  const h = canvas.height / dpr;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);

  // State 1 — no-data: system not built yet, canvas must explain itself
  if (!viewer || viewer.n === 0) {
    ctx.fillStyle = "#475569";
    ctx.font = "11px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Load a structure to begin", w / 2, h / 2 - 8);
    ctx.font = "9.5px sans-serif";
    ctx.fillStyle = "#334155";
    ctx.fillText("The 4 binding macrostates appear here", w / 2, h / 2 + 10);
    const infoEl = document.getElementById("networkInfo");
    if (infoEl) infoEl.textContent = "Next: load the 4W52 sample (Structure) → Build → Run — live pose tracking drives S0–S3.";
    return;
  }

  const states = networkModel.states;
  const curState = networkModel.currentState;
  const probs = networkModel.probabilities;

  // Synchronize 3D molecular viewer with current macrostate
  if (viewer) {
    viewer.setActiveState(curState);
  }

  const paddingX = 45;
  const nodeY = h / 2 - 8;
  const spacingX = (w - 2 * paddingX) / (states.length - 1);

  // Draw transition edges with rate labels
  for (let i = 0; i < states.length - 1; i++) {
    const x1 = paddingX + i * spacingX;
    const x2 = paddingX + (i + 1) * spacingX;
    ctx.strokeStyle = "#334155";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(x1, nodeY);
    ctx.lineTo(x2, nodeY);
    ctx.stroke();

    // Arrows
    const midX = (x1 + x2) / 2;
    ctx.fillStyle = "#64748b";
    ctx.font = "8.5px sans-serif";
    ctx.textAlign = "center";
    const kFwd = networkModel.rateMatrix[i][i + 1];
    ctx.fillText(`→ ${(kFwd / 1e6).toFixed(1)} µs⁻¹`, midX, nodeY - 6);
  }

  // Draw state nodes
  for (let i = 0; i < states.length; i++) {
    const st = states[i];
    const nx = paddingX + i * spacingX;
    const isCurrent = i === curState;
    const radius = isCurrent ? 17 : 12;

    // Active glow
    if (isCurrent) {
      ctx.beginPath();
      ctx.arc(nx, nodeY, radius + 6, 0, 2 * Math.PI);
      ctx.fillStyle = st.color + "33";
      ctx.fill();
    }

    // Probability bar
    const prob = probs[i];
    const barH = Math.max(2, prob * 32);
    ctx.fillStyle = st.color;
    ctx.fillRect(nx - 7, nodeY - radius - 5 - barH, 14, barH);

    // Node Circle
    ctx.beginPath();
    ctx.arc(nx, nodeY, radius, 0, 2 * Math.PI);
    ctx.fillStyle = st.color;
    ctx.fill();
    ctx.strokeStyle = isCurrent ? "#ffffff" : "#1e293b";
    ctx.lineWidth = isCurrent ? 2.5 : 1.5;
    ctx.stroke();

    // Node Name & ID
    ctx.fillStyle = "#ffffff";
    ctx.font = "bold 9.5px monospace";
    ctx.textAlign = "center";
    ctx.fillText(`S${i}`, nx, nodeY + 3.5);

    // Name label
    ctx.fillStyle = isCurrent ? "#38bdf8" : "#94a3b8";
    ctx.font = "9px sans-serif";
    ctx.fillText(st.name.split(" ")[0], nx, nodeY + radius + 13);

    // Energy label
    ctx.fillStyle = "#64748b";
    ctx.font = "8.5px sans-serif";
    ctx.fillText(`${st.energy.toFixed(1)}k`, nx, nodeY + radius + 24);
  }

  // Text summary & metrics
  const infoEl = document.getElementById("networkInfo");
  if (infoEl) {
    const kinetics = networkModel.computeKinetics();
    const curObj = states[curState];
    const logStr = networkModel.history.length > 0
      ? `Last hop: <b>${networkModel.history[0].fromName} → ${networkModel.history[0].toName}</b> (in ${networkModel.history[0].dtNs.toFixed(2)} ns)`
      : "No transitions yet — the live 3D pose drives this diagram (place a ligand and Run)";

    const committorStr = kinetics.committors
      ? kinetics.committors.map((q, idx) => `S${idx}:${(q * 100).toFixed(0)}%`).join(" → ")
      : "";

    infoEl.innerHTML = `
      <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:4px;">
        <div><b>Active Macrostate:</b> <span style="color:${curObj.color};font-weight:600;">S${curState} (${curObj.name})</span></div>
        <div><b>Simulated Time:</b> ${(networkModel.timeNs).toFixed(2)} ns</div>
      </div>
      <div style="font-size:11px;color:#94a3b8;margin-bottom:4px;">
        &Delta;G<sub>bind</sub>: <b>${kinetics.dG.toFixed(2)} kcal/mol</b> &middot;
        K<sub>D</sub>: <b>${kinetics.KD_uM < 1000 ? `${kinetics.KD_uM.toFixed(1)} &mu;M` : `${(kinetics.KD_uM / 1000).toFixed(2)} mM`}</b> &middot;
        k<sub>on</sub>: <b>${(kinetics.k_on / 1e6).toFixed(1)} &times; 10<sup>6</sup> M<sup>-1</sup>s<sup>-1</sup></b>
      </div>
      <div style="font-size:10.5px;color:#38bdf8;margin-bottom:4px;">
        <b>TPT Committor q<sup>+</sup>:</b> ${committorStr}
      </div>
      <div style="font-size:10.5px;color:#cbd5e1;background:#0f172a;padding:4px 6px;border-radius:4px;border:1px solid #334155;">
        ${logStr}
      </div>
    `;
  }
}

/** Steady-state redraw (wiki P2): call from the render loop; redraws at ~1 Hz
 *  even without state change so the canvas never reads as dead. */
export function networkPanelTick(now) {
  if (now - _lastSteadyDraw >= 1000) {
    _lastSteadyDraw = now;
    updateNetworkPlot();
  }
}

export function initNetworkPanel() {
  if (typeof document === "undefined") return; // headless/Node: no DOM to wire
  const stepBtn = document.getElementById("netStepBtn");
  if (stepBtn) {
    stepBtn.addEventListener("click", () => {
      networkModel.stepGillespie();
      updateNetworkPlot();
    });
  }

  const multiStepBtn = document.getElementById("netMultiStepBtn");
  if (multiStepBtn) {
    multiStepBtn.addEventListener("click", () => {
      let count = 0;
      const interval = setInterval(() => {
        networkModel.stepGillespie();
        updateNetworkPlot();
        count++;
        if (count >= 10) clearInterval(interval);
      }, 100);
    });
  }

  const masterEqBtn = document.getElementById("netMasterEqBtn");
  if (masterEqBtn) {
    masterEqBtn.addEventListener("click", () => {
      networkModel.propagateMasterEquation(1e-6, 100);
      updateNetworkPlot();
    });
  }

  const resetBtn = document.getElementById("netResetBtn");
  if (resetBtn) {
    resetBtn.addEventListener("click", () => {
      networkModel.reset();
      updateNetworkPlot();
    });
  }

  const trackChk = document.getElementById("netTrackToggle");
  if (trackChk) {
    trackChk.addEventListener("change", () => {
      isLiveTrackingActive = trackChk.checked;
    });
  }

  const canvas = document.getElementById("networkCanvas");
  if (canvas) {
    const dpr = Math.min(2, (typeof window !== "undefined" && window.devicePixelRatio) || 1);
    canvas.width = Math.floor((canvas.clientWidth || 300) * dpr);
    canvas.height = Math.floor(140 * dpr);
    updateNetworkPlot();
  }
}
