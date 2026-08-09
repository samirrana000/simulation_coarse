/**
 * main.js — Application orchestration.
 *
 * Wires the UI to the modules:
 *   pdb.js        — fetch/parse/select a structure
 *   forcefield.js — build the Cα energy function
 *   integrator.js — BAOAB Langevin dynamics
 *   viewer.js     — 60 fps Canvas rendering
 *   recorder.js   — trajectory capture & download
 *
 * Main loop: requestAnimationFrame drives the render; each frame advances the
 * simulation by ~1 ps of simulation time (auto-budgeted by CPU speed so we
 * stay near 60 fps even on slower machines). All physics is numeric —
 * nothing here is keyframe animation.
 */

import { ForceField } from "./forcefield.js?v=8";
import { Funnel } from "./funnel.js?v=8";
import { LangevinIntegrator } from "./integrator.js?v=8";
import { fetchPdb, parseCa, parseLigands, parseMol2, selectSystem, summarizeStructure } from "./pdb.js?v=8";
import { Recorder, downloadText } from "./recorder.js?v=8";
import { PoseScorer, POSE_FEATURE_N } from "./scorer.js?v=8";
import { Viewer } from "./viewer.js?v=8";
import { analyzeTrajectory, pmfCsv } from "./analysis.js?v=8";

/* ------------------------------------------------------------------ */
/*  DOM handles                                                        */
/* ------------------------------------------------------------------ */
const $ = (id) => document.getElementById(id);
const ui = {
  pdbId: $("pdbId"), fetchBtn: $("fetchBtn"), fileInput: $("fileInput"),
  structSummary: $("structSummary"), includeLig: $("includeLig"),
  mol2File: $("mol2File"), mol2Info: $("mol2Info"),
  chainsInput: $("chainsInput"), resFrom: $("resFrom"), resTo: $("resTo"),
  buildBtn: $("buildBtn"), selSummary: $("selSummary"),
  rc: $("rc"), gamma: $("gamma"), temp: $("temp"), fric: $("fric"), mass: $("mass"),
  motionGain: $("motionGain"), v_motionGain: $("v_motionGain"),
  stridePs: $("stridePs"), maxFrames: $("maxFrames"), exportFmt: $("exportFmt"),
  recBtn: $("recBtn"), recStopBtn: $("recStopBtn"), dlBtn: $("dlBtn"), recStatus: $("recStatus"),
  playBtn: $("playBtn"), resetBtn: $("resetBtn"),
  showContacts: $("showContacts"), spheres: $("spheres"),
  bindPot: $("bindPot"), holoSprings: $("holoSprings"),
  funnelToggle: $("funnelToggle"), pmfReset: $("pmfReset"), pmfPlot: $("pmfPlot"),
  nnContacts: $("nnContacts"), poseScore: $("poseScore"),
  contactsFile: $("contactsFile"), scorerFile: $("scorerFile"), nnInfo: $("nnInfo"),
  anaBtn: $("anaBtn"), anaPmfBtn: $("anaPmfBtn"), analysisOut: $("analysisOut"),
  hud: $("hud"), canvas: $("canvas"),
};
// slider live readouts
[["rc"], ["gamma"], ["temp", "v_temp"], ["fric", "v_fric"], ["mass", "v_mass"]].forEach(([k, labelId]) => {
  const el = ui[k];
  const lbl = $(labelId || `v_${k}`);
  if (el && lbl) el.addEventListener("input", () => { lbl.textContent = el.value; onParamChange(); });
});

/* ------------------------------------------------------------------ */
/*  Application state                                                  */
/* ------------------------------------------------------------------ */
const state = {
  pdbText: null,       // raw structure text
  parsed: null,        // parseCa() output
  sel: null,           // selectSystem() output
  ff: null,            // ForceField
  funnel: null,        // Funnel (binding bias + WTM PMF)
  ligands: [],
  mol2Ligands: null,   // parseMol2() output — when set, replaces HETATM ligands
  mol2Fn: null,        // filename of the loaded MOL2 (for status lines)
  integ: null,         // LangevinIntegrator
  running: false,
  simSpeedPsPerFrame: 1.0, // ~1 ps/frame at 60 fps → 1 ns per ~17 s wall time
  fpsEMA: 60,
  contacts: null,      // parsed contacts.json (array of [i,j,p])
  contactsFn: null,    // filename of loaded contacts
  mlAlpha: 1,          // NN contact map stiffness scale
  scorer: null,       // PoseScorer instance
};

const viewer = new Viewer(ui.canvas);
const recorder = new Recorder();

// View-only control — must NOT trigger a force-field rebuild.
ui.motionGain.addEventListener("input", () => {
  ui.v_motionGain.textContent = ui.motionGain.value;
  viewer.setMotionGain(Number(ui.motionGain.value));
});

/* ------------------------------------------------------------------ */
/*  Structure loading                                                  */
/* ------------------------------------------------------------------ */
async function loadStructure(text, sourceLabel) {
  state.pdbText = text;
  state.parsed = parseCa(text);
  ui.structSummary.textContent = `${sourceLabel}\n` + summarizeStructure(state.parsed);
  // Default selection = everything
  ui.chainsInput.value = "";
  ui.resFrom.value = "";
  ui.resTo.value = "";
  buildSystem();
}

ui.fetchBtn.addEventListener("click", async () => {
  const id = ui.pdbId.value.trim();
  if (!id) return;
  ui.structSummary.textContent = `Fetching ${id}…`;
  try {
    await loadStructure(await fetchPdb(id), `PDB ${id.toUpperCase()}`);
  } catch (err) {
    ui.structSummary.textContent = "⚠ " + err.message;
    ui.hud.textContent = "⚠ " + err.message;
  }
});

ui.pdbId.addEventListener("keydown", (e) => { if (e.key === "Enter") ui.fetchBtn.click(); });

document.querySelectorAll("[data-ex]").forEach((a) =>
  a.addEventListener("click", (e) => {
    e.preventDefault();
    ui.pdbId.value = a.dataset.ex;
    ui.fetchBtn.click();
  })
);

ui.fileInput.addEventListener("change", async () => {
  const f = ui.fileInput.files[0];
  if (!f) return;
  try {
    await loadStructure(await f.text(), f.name);
  } catch (err) {
    ui.structSummary.textContent = "⚠ " + err.message;
    ui.hud.textContent = "⚠ " + err.message;
  }
});

/* ------------------------------------------------------------------ */
/*  Optional MOL2 ligand input (overrides PDB HETATM ligands)          */
/* ------------------------------------------------------------------ */
ui.mol2File.addEventListener("change", async () => {
  const f = ui.mol2File.files[0];
  if (!f) {
    // file picker was cleared (or cancelled after a previous load) →
    // fall back to HETATM ligands from the PDB
    state.mol2Ligands = null;
    state.mol2Fn = null;
    ui.mol2Info.style.display = "none";
    buildSystem();
    return;
  }
  try {
    const mols = parseMol2(await f.text());
    if (!mols.length) throw new Error(`${f.name}: no usable molecules (need ≥ 2 heavy atoms)`);
    state.mol2Ligands = mols;
    state.mol2Fn = f.name;
    const nAtoms = mols.reduce((s, m) => s + m.atoms.length, 0);
    const nBonds = mols.reduce((s, m) => s + m.bonds.length, 0);
    ui.mol2Info.textContent =
      `${f.name}: ${mols.length} molecule(s), ${nAtoms} heavy atoms, ${nBonds} bonds — ` +
      `active (overrides HETATM ligands)`;
    ui.mol2Info.style.display = "block";
    buildSystem();
  } catch (err) {
    state.mol2Ligands = null;
    state.mol2Fn = null;
    ui.mol2Info.textContent = "⚠ " + err.message;
    ui.mol2Info.style.display = "block";
    buildSystem();
  }
});

/* ------------------------------------------------------------------ */
/*  ML tier: NN contact map (ESM prior) + NN pose scorer               */
/* ------------------------------------------------------------------ */
// contacts.json loader (output of ml/export_esm_contacts.py)
ui.contactsFile.addEventListener("change", async () => {
  const f = ui.contactsFile.files[0];
  if (!f) return;
  try {
    const d = JSON.parse(await f.text());
    if (!Array.isArray(d.contacts)) throw new Error("contacts.json: missing contacts[]");
    state.contacts = d.contacts;       // [[i,j,p], ...]
    state.contactsFn = d.source || f.name;
    applyMLToFF();
    ui.nnInfo.textContent = `Loaded ${state.contacts.length} pairs (${state.contactsFn})`;
  } catch (err) {
    ui.nnInfo.textContent = "⚠ " + err.message;
  }
});

// scorer weights JSON loader (PoseScorer config)
ui.scorerFile.addEventListener("change", async () => {
  const f = ui.scorerFile.files[0];
  if (!f) return;
  try {
    const d = JSON.parse(await f.text());
    state.scorer = new PoseScorer(d);
    ui.nnInfo.textContent = `Scorer: ${state.scorer.layers.join("→")} (${f.name})`;
  } catch (err) {
    ui.nnInfo.textContent = "⚠ " + err.message;
  }
});

// toggles
ui.nnContacts.addEventListener("change", applyMLToFF);
ui.poseScore.addEventListener("change", () => {
  // ensure a default scorer exists; the toggle just turns the readout on
  if (!state.scorer) state.scorer = PoseScorer.docked();
});

/**
 * Apply the current ML contact prior to the ENM springs (if loaded + enabled)
 * and (re)create the default pose scorer when the toggle is on. Called after
 * every force-field rebuild so the scaling survives gamma/rc changes.
 */
function applyMLToFF() {
  if (!state.ff) return;
  if (ui.nnContacts.checked && state.contacts) {
    state.ff.setSpringScale(state.contacts, state.mlAlpha);
  } else {
    state.ff.clearSpringScale();
  }
  if (ui.poseScore.checked && !state.scorer) state.scorer = PoseScorer.docked();
  updateSelSummary();
}

/* ------------------------------------------------------------------ */
/*  System construction (selection + force field + integrator)         */
/* ------------------------------------------------------------------ */
function parseParamChainIds() {
  const t = ui.chainsInput.value.trim();
  if (!t) return null;
  return t.split(",").map((s) => s.trim().toUpperCase()).filter(Boolean)
    .map((s) => (s === "_" ? "_" : s));
}

function buildSystem() {
  if (!state.parsed) return;
  try {
    const chains = parseParamChainIds();
    const resFrom = ui.resFrom.value === "" ? null : Number(ui.resFrom.value);
    const resTo = ui.resTo.value === "" ? null : Number(ui.resTo.value);
    state.sel = selectSystem(state.parsed, { chains, resFrom, resTo });
  } catch (err) {
    ui.selSummary.textContent = "⚠ " + err.message;
    ui.hud.textContent = "⚠ " + err.message;
    return;
  }

  const par = { rc: Number(ui.rc.value), gamma: Number(ui.gamma.value), binding: { on: ui.bindPot.checked, holo: ui.holoSprings.checked } };
  state.ligands = [];
  if (state.mol2Ligands && state.mol2Ligands.length) {
    // MOL2 ligand(s) supplied as a separate file — they replace any HETATM
    // ligands from the PDB (e.g. protein-only PDB + benzene.mol2).
    state.ligands = state.mol2Ligands;
  } else if (ui.includeLig.checked && state.pdbText) {
    try { state.ligands = parseLigands(state.pdbText); }
    catch (err) { /* ligand parsing must not break system build */ }
  }
  state.ff = new ForceField(state.sel, par, state.ligands);
  state.integ = new LangevinIntegrator(state.ff.ref, state.ff, Number(ui.mass.value));
  if (state.ff.nLigAtoms > 0) {
    state.funnel = new Funnel({ nProt: state.ff.nProt, n: state.ff.n, ref: state.ff.ref });
    state.ff.setFunnel(state.funnel);
    state.ff.funnelOn = ui.funnelToggle.checked;
  } else {
    state.funnel = null;
  }
  applyMLToFF();
  onParamChange(false);

  viewer.setSystem(state.sel, state.ff);
  recorder.clear();
  updateRecStatus();
  updateSelSummary();

  state.running = false;
  ui.playBtn.textContent = "▶ Run";
  ui.playBtn.disabled = false;
  ui.resetBtn.disabled = false;
  state._prevPos = null;   // position buffer length may change on rebuild
}
function updateSelSummary() {
  if (!state.ff) return;
  const nat = state.ff.springs.length / 3;
  let ligSummary = "";
  if (state.ff.nLigAtoms > 0) {
    ligSummary = ` · ${state.ff.nLigAtoms} ligand atom(s) in ${state.ligands.length} molecule(s)`;
    if (state.mol2Ligands && state.mol2Ligands.length) ligSummary += ` (MOL2: ${state.mol2Fn})`;
    if (state.ff.nHolo > 0) ligSummary += ` · ${state.ff.nHolo} holo contacts`;
  }
  let nnSummary = "";
  if (state.ff.springScaleActive) {
    nnSummary = ` · NN map active (${state.contactsFn || "contacts.json"})`;
  }
  ui.selSummary.textContent =
    `${state.sel.beads.length} Cα beads · ${state.sel.nChains} chain(s) · ` +
    `${state.ff.bonds.length / 3} bonds · ${state.ff.angles.length / 4} angles · ` +
    `${nat} elastic contacts (Rc=${Number(ui.rc.value)} Å)${ligSummary}${nnSummary}`;
}
ui.buildBtn.addEventListener("click", buildSystem);
ui.includeLig.addEventListener("change", buildSystem);
ui.bindPot.addEventListener("change", () => onParamChange(true));
ui.holoSprings.addEventListener("change", () => onParamChange(true));
ui.funnelToggle.addEventListener("change", () => {
  if (state.ff) state.ff.funnelOn = ui.funnelToggle.checked;
});
ui.pmfReset.addEventListener("click", () => { if (state.funnel) state.funnel.reset(); });

/* ------------------------------------------------------------------ */
/*  Physics parameter hot-reload                                       */
/* ------------------------------------------------------------------ */
function onParamChange(rebuildContacts = true) {
  if (!state.integ || !state.ff) return;
  const integ = state.integ;
  if (state.ff.funnel) state.ff.funnel.setTemperature(Number(ui.temp.value));
  integ.setTemperature(Number(ui.temp.value));
  integ.setFriction(Number(ui.fric.value));
  integ.rebuildMass(Number(ui.mass.value));
  // ENM parameter change ⇒ cheap full rebuild of contact list
  if (rebuildContacts) {
  const par = { rc: Number(ui.rc.value), gamma: Number(ui.gamma.value), binding: { on: ui.bindPot.checked, holo: ui.holoSprings.checked } };
    const keepPos = Float64Array.from(integ.pos);
    const keepVel = Float64Array.from(integ.vel);
    const keepTime = integ.time;
    state.ff = new ForceField(state.sel, par, state.ligands);
    state.integ = new LangevinIntegrator(state.ff.ref, state.ff, Number(ui.mass.value));
    state.integ.pos.set(keepPos);
    state.integ.vel.set(keepVel);
    state.integ.time = keepTime;
    state._prevPos = null;
    state.integ.setTemperature(Number(ui.temp.value));
    state.integ.setFriction(Number(ui.fric.value));
    if (state.ff.nLigAtoms > 0) {
      state.funnel = new Funnel({ nProt: state.ff.nProt, n: state.ff.n, ref: state.ff.ref });
      state.ff.setFunnel(state.funnel);
      state.ff.funnelOn = ui.funnelToggle.checked;
    } else {
      state.funnel = null;
    }
    applyMLToFF();
    viewer.setSystem(state.sel, state.ff);
    updateSelSummary();
  }
}

/* ------------------------------------------------------------------ */
/*  Transport controls                                                 */
/* ------------------------------------------------------------------ */
ui.playBtn.addEventListener("click", () => {
  if (!state.integ) return;
  state.running = !state.running;
  ui.playBtn.textContent = state.running ? "⏸ Pause" : "▶ Run";
});
ui.resetBtn.addEventListener("click", () => {
  if (!state.integ) return;
  state.integ.reset();
});

ui.showContacts.addEventListener("change", () => (viewer.showContacts = ui.showContacts.checked));
ui.spheres.addEventListener("change", () => (viewer.drawSpheres = ui.spheres.checked));

/* ------------------------------------------------------------------ */
/*  Recording                                                          */
/* ------------------------------------------------------------------ */
ui.recBtn.addEventListener("click", () => {
  if (!state.integ) return;
  recorder.start(state.integ.time, Number(ui.stridePs.value), Number(ui.maxFrames.value));
  ui.recBtn.classList.add("rec-on");
});
ui.recStopBtn.addEventListener("click", () => {
  recorder.stop();
  ui.recBtn.classList.remove("rec-on");
});
ui.dlBtn.addEventListener("click", () => {
  if (!state.sel) return;
  const fmt = ui.exportFmt.value;
  try {
    const text = recorder.buildFile(fmt, state.sel.beads);
    const name = `cg_traj_${recorder.count}frames.${fmt}`;
    downloadText(text, name);
  } catch (err) {
    ui.recStatus.textContent = "⚠ " + err.message;
  }
});
function updateRecStatus() {
  ui.recStatus.textContent =
    `${recorder.count} frames · ${recorder.spanNs.toFixed(3)} ns recorded` +
    (recorder.recording ? " ●" : "");
  ui.dlBtn.disabled = recorder.count === 0;
}

/* ------------------------------------------------------------------ */
/*  Analysis (panel 6)                                                 */
/* ------------------------------------------------------------------ */
ui.anaBtn.addEventListener("click", () => {
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

ui.anaPmfBtn.addEventListener("click", () => {
  if (!state.funnel) { ui.analysisOut.textContent = "⚠ No active funnel (load a ligand + build the system)."; return; }
  try {
    downloadText(pmfCsv(state.funnel), `pmf_${(state.contactsFn || "enm").replace(/\W+/g, "_")}.csv`);
  } catch (err) {
    ui.analysisOut.textContent = "⚠ " + err.message;
  }
});

/* ------------------------------------------------------------------ */
/*  Main loop — physics & render, decoupled from display refresh       */
/* ------------------------------------------------------------------ */
let lastT = performance.now();
function tick(now) {
  requestAnimationFrame(tick);
  const dtWall = Math.min(0.05, (now - lastT) / 1000); // clamp tab-switch jumps
  lastT = now;
  state.fpsEMA = 0.95 * state.fpsEMA + 0.05 / Math.max(1e-3, dtWall);

  if (state.running && state.integ) {
    // Target simSpeed ps of dynamics per displayed frame at 60 fps.
    // dt is in ps ⇒ steps/frame = simSpeed / dt. `advance()` also
    // time-slices so slow machines degrade sim speed, not fps.
    const steps = Math.max(1, Math.round(state.simSpeedPsPerFrame / state.integ.dt));
    state.integ.advance(steps, 11); // ≤11 ms physics budget per frame

    const captured = recorder.maybeCapture(state.integ.pos, state.integ.time);
    if (captured || recorder.recording || recorder.count > 0) updateRecStatus();
  }

  if (state.integ) {
    viewer.render(state.integ.pos);
    const ff = state.ff, integ = state.integ;
    // Per-frame mean bead displacement — direct proof the dynamics are moving.
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
      // protein–ligand contacts within 5.5 Å
      let nc = 0;
      for (let i = 0; i < nProt; i++) {
        for (let la = nProt; la < state.ff.n; la++) {
          const dx = p[3*la]-p[3*i], dy = p[3*la+1]-p[3*i+1], dz = p[3*la+2]-p[3*i+2];
          if (dx*dx+dy*dy+dz*dz < 30.25) nc++;
        }
      }
      // ligand RMSD to native pose
      let s = 0;
      for (let la = nProt; la < state.ff.n; la++) {
        const dx = p[3*la]-state.ff.ref[3*la], dy = p[3*la+1]-state.ff.ref[3*la+1], dz = p[3*la+2]-state.ff.ref[3*la+2];
        s += dx*dx+dy*dy+dz*dz;
      }
      const ligRmsd = Math.sqrt(s / state.ff.nLigAtoms);
      const dg = fn && fn.active && !Number.isNaN(fn.estimateDG()) ? fn.estimateDG().toFixed(2) : "–";
      const cv = fn && fn.active ? fn.lastCV.toFixed(2) : "–";
      extra = `CV = ${cv} Å  ·  ligRMSD = ${ligRmsd.toFixed(2)} Å  ·  nContacts = ${nc}  ·  ΔG ≈ ${dg} kcal/mol  ·  `;
      // NN pose score (MLP) — feature vector per scorer.POSE_FEATURE_N
      if (ui.poseScore.checked && state.scorer) {
        let clashes = 0;
        for (let i = 0; i < nProt; i++) {
          for (let la = nProt; la < state.ff.n; la++) {
            const dx = p[3*la]-p[3*i], dy = p[3*la+1]-p[3*i+1], dz = p[3*la+2]-p[3*i+2];
            if (dx*dx+dy*dy+dz*dz < 12.25) clashes++;   // < 3.5 Å
          }
        }
        const feat = [
          Math.min(1, nc / 20),
          Math.max(0, Math.min(1.2, -ff.desolvU / 3.3)),
          Math.min(1, clashes / 5),
          Math.min(1, ligRmsd / 4),
          Math.min(1, (fn && fn.active ? fn.lastCV : 0) / 10),
          Math.max(0, Math.min(1, -ff.bindingU / 2)),
        ];
        const score = state.scorer.predict(feat);
        extra += `pose = ${score.toFixed(2)}  ·  `;
      }
    }
    ui.hud.textContent =
      `t = ${integ.time.toFixed(1)} ps (${(integ.time / 1000).toFixed(3)} ns)  ·  ` +
      `U = ${ff.energy.toFixed(1)} kcal/mol  ·  ` +
      (ff.nLigAtoms > 0 ? `U_bind = ${ff.bindingU.toFixed(2)} kcal/mol  ·  ` : '') +
      extra +
      `RMSD = ${ff.rmsd(integ.pos).toFixed(2)} Å  ·  ` +
      `T_inst = ${ff.kineticTemp(integ.vel, integ.mass).toFixed(0)} K  ·  ` +
      `Δr = ${dr.toFixed(2)} Å/frame  ·  ` +
      `${state.fpsEMA.toFixed(0)} fps  (${state.ff.nProt} Cα${state.ff.nLigAtoms ? ` + ${state.ff.nLigAtoms} lig` : ''})`;
    updatePMFPlot();
  } else {
    viewer.render(null);
  }
}
requestAnimationFrame(tick);

/* ------------------------------------------------------------------ */
/*  Binding PMF plot (well-tempered metadynamics reconstruction)        */
/* ------------------------------------------------------------------ */
function updatePMFPlot() {
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
