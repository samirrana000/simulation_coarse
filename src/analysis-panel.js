/**
 * analysis-panel.js — panel 6: trajectory analysis report + PMF CSV download
 * (item 5 modularization: moved verbatim from main.js; only the imports
 * changed — `ui` / `state` / `recorder` now come from ui.js).
 *
 * Phase 4: opt-in translational workflows (alanine scanning, DCCM heatmap,
 * SMD unbinding ensemble, cryptic-pocket tracking) as buttons that degrade
 * gracefully when no system / trajectory / ligand is present. The original
 * Analyze flow is untouched.
 */

import { analyzeTrajectory, pmfCsv } from "./analysis.js?v=10";
import { downloadText } from "./recorder.js?v=10";
import { ui, state, viewer, recorder } from "./ui.js?v=10";
import { scanPocket, pocketResidues, formatMutationTable } from "./analysis/alanine_scanning.js?v=10";
import { computeDCCM, renderDCCMHeatmap, topCorrelations, dccmPick, highlightCorrelatedPair } from "./analysis/dccm.js?v=10";
import { trackPocketVolume, detectCryptic } from "./analysis/cryptic_pockets.js?v=10";
import { runPullingEnsemble, jarzynskiFreeEnergy, koffSurrogate } from "./analysis/unbinding_smd.js?v=10";
import { computeThermodynamics, formatThermoTable } from "./analysis/thermodynamics.js?v=10";
import { ForceField } from "./forcefield.js?v=10";
import { LangevinIntegrator } from "./integrator.js?v=10";

if (ui.anaBtn) ui.anaBtn.addEventListener("click", () => {
  if (recorder.count === 0) {
    ui.analysisOut.textContent = "⚠ Nothing recorded yet — press ● Rec, run the sim, then Stop.";
    return;
  }
  ui.analysisOut.textContent = "Analyzing… (RMSIP may take a few seconds)";
  // run off the animation frame so the status text paints first
  setTimeout(() => {
    try {
      const rep = analyzeTrajectory({
        frames: recorder.frames,
        times: recorder.times,
        ref: state.ff.ref,
        nProt: state.ff.nProt,
        n: state.ff.n,
        beads: state.sel.beads,
        ff: state.ff,
        funnel: state.funnel,
      });
      ui.analysisOut.textContent = rep.lines.join("\n");
    } catch (err) {
      ui.analysisOut.textContent = "⚠ " + err.message;
    }
  }, 20);
});

if (ui.anaPmfBtn) ui.anaPmfBtn.addEventListener("click", () => {
  if (!state.funnel) { ui.analysisOut.textContent = "⚠ No active funnel (load a ligand + build the system)."; return; }
  try {
    downloadText(pmfCsv(state.funnel), `pmf_${(state.contactsFn || "enm").replace(/\W+/g, "_")}.csv`);
  } catch (err) {
    ui.analysisOut.textContent = "⚠ " + err.message;
  }
});

/* ------------------------------------------------------------------ */
/*  Phase 4 — translational workflows (opt-in, gracefully degrading)   */
/* ------------------------------------------------------------------ */

/** Current live system as a Phase-4 ScanSystem, or null when not built. */
function currentScanSystem() {
  if (!state.ff || !state.sel) return null;
  return {
    mode: state.heavyMode ? "heavy" : "cg",
    sel: state.sel,
    ff: state.ff,
    par: undefined, // rebuilt from the live ff (parFromCgFF / gamma)
    ligands: state.ligands ?? [],
  };
}

function fmtE(v) {
  return Number.isFinite(v) ? v.toFixed(2) : "NaN";
}

function setDccmCaption(txt) {
  if (ui.dccmCaption) ui.dccmCaption.textContent = txt;
  else {
    const el = typeof document !== "undefined" ? document.getElementById("dccmCaption") : null;
    if (el) el.textContent = txt;
  }
}

let _dccmHasData = false;
let _lastDccmEmpty = 0;

/** Three-state (P2) no-data: actionable empty DCCM with caption + legend key. */
export function drawDccmEmpty(reason) {
  _dccmHasData = false;
  const cv = (ui.dccmCanvas) || (typeof document !== "undefined" ? document.getElementById("dccmCanvas") : null);
  if (!cv) return;
  try {
    const w = cv.clientWidth || 280, h = cv.clientHeight || 180;
    if (w === 0 || h === 0) return;
    const dpr = Math.min(2, (typeof window !== "undefined" && window.devicePixelRatio) || 1);
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
      cv.width = Math.round(w * dpr);
      cv.height = Math.round(h * dpr);
    }
    const ctx = cv.getContext("2d");
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = "#475569";
    ctx.font = "11px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(reason || "Record ≥ 3 frames, then DCCM Heatmap", w / 2, h / 2 - 6);
    ctx.font = "9.5px sans-serif";
    ctx.fillStyle = "#334155";
    ctx.fillText("cross-correlation −1 (blue) … +1 (amber)", w / 2, h / 2 + 10);
  } catch (_) { /* headless: caption only */ }
  setDccmCaption("DCCM: record ≥ 3 frames, then compute. Legend −1…+1.");
}

/** Steady-state (P2): ≤1 Hz empty redraw so the canvas never reads as dead. */
export function dccmTick(now) {
  if (_dccmHasData) return;
  const t = now ?? ((typeof performance !== "undefined" && performance.now()) || Date.now());
  if (t - _lastDccmEmpty < 1000) return;
  _lastDccmEmpty = t;
  drawDccmEmpty();
}

// Paint the empty state once at startup (no-data, before any trajectory).
try { drawDccmEmpty(); } catch (_) {}

// ---- Thermodynamics ΔH/ΔS (Loop-2 S5) -------------------------------
// Holo leg = recorded trajectory (needs ≥30 frames); apo leg = internal
// relaxation of a binding-off clone (same ENM, no ligand coupling).
if (ui.thermoBtn) ui.thermoBtn.addEventListener("click", () => {
  try {
    if (!state.ff) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
    if (recorder.frames.length < 30) {
      ui.analysisOut.textContent = `⚠ Need ≥30 recorded holo frames for the Schlitter covariance (have ${recorder.frames.length}). ● Rec, run ~100 ps, Stop, retry.`;
      return;
    }
    const ff = state.ff;
    const nProt = ff.nProt;
    // ligand COM in ref → pocket residues
    let lcom = [0, 0, 0];
    if (ff.nLigAtoms > 0) {
      for (let a = 0; a < ff.nLigAtoms; a++) {
        lcom[0] += ff.ref[3 * (nProt + a)] / ff.nLigAtoms;
        lcom[1] += ff.ref[3 * (nProt + a) + 1] / ff.nLigAtoms;
        lcom[2] += ff.ref[3 * (nProt + a) + 2] / ff.nLigAtoms;
      }
    } else lcom = centroidOf(ff.ref, nProt);
    const pocketIdx = [];
    for (let i = 0; i < nProt; i++) {
      if (Math.hypot(ff.ref[3 * i] - lcom[0], ff.ref[3 * i + 1] - lcom[1], ff.ref[3 * i + 2] - lcom[2]) < 8.0) pocketIdx.push(i);
    }
    // holo frames from recorder; energies from BindLog if capture was on
    const holoFrames = recorder.frames.map((f) => Float32Array.from(f));
    const holoEnergies = [];
    if (state.bindLog && state.bindLog.nEvents > 0) {
      // collect the last 7-term energy snapshot per captured frame time
      const perTime = new Map();
      for (let i = 0; i < state.bindLog.nEvents; i++) {
        if (state.bindLog.evType[i] !== 0) continue;
        const t = state.bindLog.evTime[i];
        if (!perTime.has(t)) perTime.set(t, new Array(7).fill(0));
        perTime.get(t)[state.bindLog.evA[i] % 7] = state.bindLog.evX[i];
      }
      for (const row of perTime.values()) holoEnergies.push(row);
    }
    // apo leg: internal 2000-step relaxation with binding off
    const ffApo = new ForceField(state.sel, { rc: 10, gamma: 2.0, binding: { on: false } }, state.ligands ?? []);
    const integApo = new LangevinIntegrator(ffApo.ref, ffApo, 110.0);
    integApo.setTemperature(300); integApo.setFriction(8.0);
    const apoFrames = [];
    const nApo = Math.min(2000, Math.max(600, holoFrames.length * 4));
    for (let s = 0; s < nApo; s++) {
      integApo.step();
      if (s % 2 === 0) apoFrames.push(Float32Array.from(integApo.pos));
    }
    const res = computeThermodynamics({
      holoFrames, apoFrames, holoEnergies,
      pocketIdx, nProt, mass: 110, T: 300,
    });
    ui.analysisOut.textContent = formatThermoTable(res) +
      (holoEnergies.length ? "" : "\n(note: BindLog capture was off — ΔH from recorded-frame recomputation skipped; run with BindLog on for the component split)");
    if (ui.thermoCaption) ui.thermoCaption.textContent = `ΔH/ΔS done — ${pocketIdx.length} pocket residues, ${holoFrames.length} holo frames.`;
  } catch (err) {
    ui.analysisOut.textContent = "⚠ " + (err && err.message ? err.message : String(err));
  }
});
function centroidOf(ref, n) {
  let x = 0, y = 0, z = 0;
  for (let i = 0; i < n; i++) { x += ref[3 * i]; y += ref[3 * i + 1]; z += ref[3 * i + 2]; }
  return [x / n, y / n, z / n];
}

// ---- Alanine scanning ----------------------------------------------
if (ui.alaScanBtn) ui.alaScanBtn.addEventListener("click", () => {
  const sys = currentScanSystem();
  if (!sys) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
  if ((sys.ff.n - (sys.ff.ligandStart ?? sys.ff.nProt)) <= 0) {
    ui.analysisOut.textContent = "⚠ Ala-scan needs a holo system (load/place a ligand, then Build System).";
    return;
  }
  let ids;
  const txt = ((ui.alaRes && ui.alaRes.value) || "").trim();
  try {
    if (txt) {
      ids = txt.split(/[,\s;]+/).filter(Boolean).map((s) => (/^-?\d+$/.test(s) && sys.mode === "cg" ? Number(s) : s));
    } else {
      ids = pocketResidues(sys, { rCut: 6.0, maxN: 6 }).map((p) => p.resId);
    }
  } catch (err) { ui.analysisOut.textContent = "⚠ " + err.message; return; }
  if (!ids.length) { ui.analysisOut.textContent = "⚠ No pocket residues found (ligand has no protein neighbours ≤ 6 Å)."; return; }
  ids = ids.slice(0, 8); // UI responsiveness cap; script API has no cap
  ui.analysisOut.textContent = `Alanine scanning ${ids.length} residue(s)… (short relaxations)`;
  setTimeout(() => {
    try {
      const { rows, wtHolo, wtApo } = scanPocket(sys, ids, { relaxSteps: sys.mode === "heavy" ? 25 : 80 });
      ui.analysisOut.textContent =
        `Ala-scan (${sys.mode}): WT holo ${fmtE(wtHolo)} · WT apo ${fmtE(wtApo)} kcal/mol\n` +
        formatMutationTable(rows);
    } catch (err) {
      ui.analysisOut.textContent = "⚠ " + err.message;
    }
  }, 20);
});

// ---- DCCM heatmap ----------------------------------------------------
if (ui.dccmBtn) ui.dccmBtn.addEventListener("click", () => {
  const sys = currentScanSystem();
  if (!sys) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
  if (recorder.count < 3) {
    ui.analysisOut.textContent = "⚠ DCCM needs ≥ 3 recorded frames — press ● Rec, run, then Stop.";
    return;
  }
  ui.analysisOut.textContent = "Computing DCCM…";
  setTimeout(() => {
    try {
      const dcc = computeDCCM(recorder.frames, { nProt: sys.ff.nProt, ref: sys.ff.ref });
      const cv = ui.dccmCanvas;
      if (cv) {
        cv.style.display = "block";
        renderDCCMHeatmap(cv, dcc.matrix, { n: dcc.n, size: 220 });
        _dccmHasData = true;
        setDccmCaption(`DCCM live: ${dcc.n}×${dcc.n} from ${dcc.nFrames} frames — click heatmap for pair. Legend −1…+1.`);
        if (!cv._dccmClick) {
          cv._dccmClick = true;
          cv.addEventListener("click", (e) => {
            const p = dccmPick(cv, e.clientX, e.clientY);
            if (p && viewer) {
              highlightCorrelatedPair(viewer, p.i, p.j);
              ui.analysisOut.textContent += `\n◉ pair (${p.i}, ${p.j}) C=${dcc.matrix[p.i * dcc.n + p.j].toFixed(2)} highlighted in viewer.`;
            }
          });
        }
      }
      const top = topCorrelations(dcc, 6);
      let meanAbs = 0;
      for (let i = 0; i < dcc.matrix.length; i++) meanAbs += Math.abs(dcc.matrix[i]);
      meanAbs /= dcc.matrix.length;
      ui.analysisOut.textContent =
        `DCCM: ${dcc.n}×${dcc.n} from ${dcc.nFrames} frames · mean|C| = ${meanAbs.toFixed(3)}\n` +
        `top pairs: ${top.map((t) => `(${t.i},${t.j}) ${t.c >= 0 ? "+" : ""}${t.c.toFixed(2)}`).join(" · ")}\n` +
        `click the heatmap to highlight a pair in the viewer.`;
    } catch (err) {
      ui.analysisOut.textContent = "⚠ " + err.message;
    }
  }, 20);
});

// ---- SMD unbinding ensemble ------------------------------------------
if (ui.smdBtn) ui.smdBtn.addEventListener("click", () => {
  const sys = currentScanSystem();
  if (!sys) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
  if ((sys.ff.n - (sys.ff.ligandStart ?? sys.ff.nProt)) <= 0) {
    ui.analysisOut.textContent = "⚠ SMD needs a holo system (load/place a ligand, then Build System).";
    return;
  }
  ui.analysisOut.textContent = "SMD: 4 pulls…";
  setTimeout(() => {
    try {
      const pos0 = (state.integ && state.integ.pos.length === sys.ff.n * 3)
        ? Float64Array.from(state.integ.pos) : sys.ff.ref;
      const ens = runPullingEnsemble({
        ff: sys.ff, pos0, nPulls: 4, mode: "velocity", seed: 42,
        pull: { k: 5.0, v: 4.0, nSteps: 2000, dt: 0.0015, sRupture: 6.0, T: Number(ui.temp?.value || 300) },
      });
      const je = jarzynskiFreeEnergy(ens.works, { T: Number(ui.temp?.value || 300) });
      const koff = koffSurrogate(ens.pulls, je.dF, { T: Number(ui.temp?.value || 300) });
      ui.analysisOut.textContent =
        `SMD ×4 (constant-velocity, k=2, v=2 Å/ps):\n` +
        `works = [${Array.from(ens.works).map(fmtE).join(", ")}] kcal/mol\n` +
        `Jarzynski ΔF = ${fmtE(je.dF)}${je.se !== null ? ` ± ${je.se.toFixed(2)}` : ""} · ⟨W⟩ = ${fmtE(je.meanWork)} · dissipated = ${fmtE(je.dissipated)}\n` +
        `rupture ⟨F⟩ = ${fmtE(koff.meanRuptureForce)} kcal/mol/Å · koff-score = ${fmtE(koff.meanScore)} kT (${koff.note})`;
    } catch (err) {
      ui.analysisOut.textContent = "⚠ " + err.message;
    }
  }, 20);
});

// ---- Cryptic pockets ---------------------------------------------------
if (ui.crypticBtn) ui.crypticBtn.addEventListener("click", () => {
  const sys = currentScanSystem();
  if (!sys) { ui.analysisOut.textContent = "⚠ Build a system first."; return; }
  try {
    // Pocket set: funnel pocket, else protein within 8 Å of the ligand COM.
    let pocket = (state.funnel && state.funnel.active && state.funnel.pocket)
      ? [...state.funnel.pocket] : null;
    if (!pocket) {
      pocket = pocketResidues(sys, { rCut: 8.0, maxN: 40 }).map((p) =>
        (sys.mode === "cg" ? p.resId : sys.sel.atoms.findIndex((a) =>
          a.isProtein && a.atomName === "CA" && `${a.chain}|${a.resSeq}` === p.resId)));
      pocket = pocket.filter((i) => i >= 0);
    }
    if (!pocket.length) { ui.analysisOut.textContent = "⚠ No pocket residues (needs a ligand)."; return; }
    const frames = recorder.count >= 2 ? recorder.frames
      : [(state.integ && state.integ.pos.length === sys.ff.n * 3) ? state.integ.pos : sys.ff.ref];
    const track = trackPocketVolume(frames, { pocketIndices: pocket });
    const det = detectCryptic(track, {
      kSigma: 1.5, frames: frames.length >= 2 ? frames : null,
      pocketIndices: pocket, nProt: sys.ff.nProt,
    });
    const vols = Array.from(track.volumes);
    const topRes = (det.heatmap ? [...det.heatmap.perResidue].sort((a, b) => b.score - a.score).slice(0, 5) : []);
    ui.analysisOut.textContent =
      `Cryptic pockets (${pocket.length} residues, ${track.nFrames} frame(s)):\n` +
      `V range ${fmtE(Math.min(...vols))}–${fmtE(Math.max(...vols))} Å³ · ⟨V⟩ ${fmtE(det.mean)} · open threshold ${fmtE(det.threshold)}\n` +
      `open ${(det.openFrac * 100).toFixed(1)}% · ${det.events.length} event(s)` +
      (det.events.length ? `: ${det.events.slice(0, 4).map((e) => `#${e.start}–${e.end} peak ${fmtE(e.peak)}`).join("; ")}` : "") +
      (topRes.length ? `\nprobe-accessible: ${topRes.map((r) => `#${r.index} ${(r.score * 100).toFixed(0)}%`).join(" · ")}` : "");
  } catch (err) {
    ui.analysisOut.textContent = "⚠ " + err.message;
  }
});
