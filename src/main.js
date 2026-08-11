/**
 * main.js — Application orchestration.
 *
 * Wires the UI to the modules:
 *   pdb.js        — fetch/parse/select a structure
 *   forcefield.js — build the Cα energy function
 *   integrator.js — BAOAB Langevin dynamics
 *   viewer.js     — 60 fps Canvas rendering
 *   recorder.js   — trajectory capture & download
 *   ui.js         — DOM handles + shared state/viewer/recorder singletons
 *   ml-tier.js    — NN contact map + pose scorer panel
 *   pmf-panel.js  — funnel PMF plot + reset
 *   analysis-panel.js — trajectory analysis report + PMF CSV
 *
 * Main loop: requestAnimationFrame drives the render; each frame advances the
 * simulation by ~1 ps of simulation time (auto-budgeted by CPU speed so we
 * stay near 60 fps even on slower machines). All physics is numeric —
 * nothing here is keyframe animation.
 */

import { ForceField } from "./forcefield.js?v=10";
import { Funnel } from "./funnel.js?v=10";
import { LangevinIntegrator } from "./integrator.js?v=10";
import { fetchPdb, parseCa, parseLigands, parseMol2, selectSystem, summarizeStructure } from "./pdb.js?v=10";
import { downloadText } from "./recorder.js?v=10";
import { PoseScorer } from "./scorer.js?v=10";
import { ui, state, viewer, recorder, initParamReadouts, updateSelSummary, updateRecStatus } from "./ui.js?v=10";
import { applyMLToFF } from "./ml-tier.js?v=10";
import { updatePMFPlot } from "./pmf-panel.js?v=10";
import "./analysis-panel.js?v=10";  // side-effect: registers panel-6 listeners
import "./ligand-panel.js?v=10";   // side-effect: registers library+placement panel
import { parseHeavy, HeavyForceField } from "./heavy.js?v=10";

// physics-slider live readouts (rc/gamma/temp/fric/mass) → hot param reload
initParamReadouts(() => onParamChange());

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
  state.parsedHeavy = parseHeavy(text);
  state.libraryLigand = null;  // a new structure clears any placed-library ligand
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
/*  System construction (selection + force field + integrator)         */
/* ------------------------------------------------------------------ */
function parseParamChainIds() {
  const t = ui.chainsInput.value.trim();
  if (!t) return null;
  return t.split(",").map((s) => s.trim().toUpperCase()).filter(Boolean)
    .map((s) => (s === "_" ? "_" : s));
}

/**
 * Filter parseHeavy() output by chain / residue range (mirrors selectSystem).
 * Returns a { atoms, beads, segments } object shaped for buildSystem: `beads`
 * carries the heavy atoms (with element/isProtein/isMetal) so the viewer can
 * element-color them, and `segments` is empty because heavy mode has no Cα
 * backbone trace.
 */
export function selectHeavy(parsedHeavy, { chains = null, resFrom = null, resTo = null } = {}) {
  const atoms = parsedHeavy.atoms.filter((a) => {
    if (chains && !chains.includes(a.chain)) return false;
    if (resFrom !== null && a.resSeq < resFrom) return false;
    if (resTo !== null && a.resSeq > resTo) return false;
    return true;
  });
  const beads = atoms.map((a) => ({ ...a, x: a.x, y: a.y, z: a.z }));
  return { atoms, beads, segments: [], heavy: true };
}

/**
 * Append resolved EXTERNAL ligand molecules (library / MOL2 shape) onto a
 * heavy-mode selection as non-protein heavy atoms matching the parseHeavy()
 * field contract.
 *
 * Note: parseHeavy() already includes PDB HETATM cofactors/metals in
 * state.sel.atoms (it drops only water), so appending parseLigands() HETATM
 * again would double-count them. Only a separate library / MOL2 ligand needs
 * to be merged — it is added alongside the PDB cofactors, mirroring the CG
 * priority (library > MOL2 > PDB HETATM, the latter already present).
 *
 * @param {object} sel     heavy selection { atoms, beads, segments, heavy }
 * @param {Array|null} molecules  library/MOL2 molecules (or null)
 * @returns {object} a new selection with the ligand atoms/beads appended, or
 *                   the original sel unchanged when molecules is empty
 */
export function appendHeavyLigands(sel, molecules) {
  if (!molecules || molecules.length === 0) return sel;
  const atoms = sel.atoms.slice();
  const beads = sel.beads.slice();
  let resSeq = 1;
  for (const mol of molecules) {
    const resName = (mol.resName || "LIG").slice(0, 3).toUpperCase();
    const chain = mol.chain || "L";
    for (const at of mol.atoms) {
      const atom = {
        x: at.x, y: at.y, z: at.z,
        element: at.element,
        resName, chain, resSeq,
        atomName: at.element,
        serial: at.serial ?? 0,
        isProtein: false, isMetal: false, isLigand: true,
      };
      atoms.push(atom);
      beads.push({ ...atom });
    }
    resSeq++;
  }
  return { atoms, beads, segments: sel.segments, heavy: true };
}

export function buildSystem() {
  if (!state.parsed) return;
  state.heavyMode = ui.modelMode?.value === "heavy";
  try {
    const chains = parseParamChainIds();
    const resFrom = ui.resFrom.value === "" ? null : Number(ui.resFrom.value);
    const resTo = ui.resTo.value === "" ? null : Number(ui.resTo.value);
    if (state.heavyMode) {
      state.sel = selectHeavy(state.parsedHeavy, { chains, resFrom, resTo });
    } else {
      state.sel = selectSystem(state.parsed, { chains, resFrom, resTo });
    }
  } catch (err) {
    ui.selSummary.textContent = "⚠ " + err.message;
    ui.hud.textContent = "⚠ " + err.message;
    return;
  }

  const par = { rc: Number(ui.rc.value), gamma: Number(ui.gamma.value), binding: { on: ui.bindPot.checked, holo: ui.holoSprings.checked } };
  if (state.heavyMode) {
    // All-atom heavy mode: covalent topology + metal coordination, no ENM.
    // parseHeavy() already includes PDB HETATM cofactors/metals in
    // state.sel.atoms, so only an external library / MOL2 ligand is merged
    // (CG priority: library > MOL2 > PDB HETATM, the HETATM already present).
    state.ligands = [];
    let external = null;
    if (state.libraryLigand) {
      external = [state.libraryLigand];
      state.ligands = external;
    } else if (state.mol2Ligands && state.mol2Ligands.length) {
      external = state.mol2Ligands;
      state.ligands = external;
    }
    state.sel = appendHeavyLigands(state.sel, external);
    state.ff = new HeavyForceField({ atoms: state.sel.atoms }, par, []);
  } else {
    state.ligands = [];
    if (state.libraryLigand) {
      // Library ligand (placed via ligand-panel.js) takes priority over any
      // MOL2/HETATM ligand — it is the hypothesis being tested.
      state.ligands = [state.libraryLigand];
    } else if (state.mol2Ligands && state.mol2Ligands.length) {
      // MOL2 ligand(s) supplied as a separate file — they replace any HETATM
      // ligands from the PDB (e.g. protein-only PDB + benzene.mol2).
      state.ligands = state.mol2Ligands;
    } else if (ui.includeLig.checked && state.pdbText) {
      try { state.ligands = parseLigands(state.pdbText); }
      catch (err) { /* ligand parsing must not break system build */ }
    }
    state.ff = new ForceField(state.sel, par, state.ligands);
  }
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
  state.nanWarning = false;
}
ui.buildBtn.addEventListener("click", buildSystem);
ui.modelMode.addEventListener("change", buildSystem);
ui.includeLig.addEventListener("change", buildSystem);
ui.bindPot.addEventListener("change", () => onParamChange(true));
ui.holoSprings.addEventListener("change", () => onParamChange(true));

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
    if (state.heavyMode) {
      state.ff = new HeavyForceField({ atoms: state.sel.atoms }, par, []);
    } else {
      state.ff = new ForceField(state.sel, par, state.ligands);
    }
    state.integ = new LangevinIntegrator(state.ff.ref, state.ff, Number(ui.mass.value));
    state.integ.pos.set(keepPos);
    state.integ.vel.set(keepVel);
    state.integ.time = keepTime;
    state._prevPos = null;
    state.nanWarning = false;
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
  state.nanWarning = false;
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

    // Non-finite energy ⇒ the system blew past every FF guard (native clash,
    // corrupt input, …). Auto-pause instead of integrating garbage further;
    // the HUD warning below tells the user to reset or rebuild the system.
    if (!Number.isFinite(state.ff.energy)) {
      state.running = false;
      ui.playBtn.textContent = "▶ Run";
      state.nanWarning = true;
    }

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
      (state.nanWarning
        ? `⚠ NON-FINITE ENERGY — simulation auto-paused. Reset (⟲) to recover; the input structure may contain clashes.  ·  `
        : "") +
      `t = ${integ.time.toFixed(1)} ps (${(integ.time / 1000).toFixed(3)} ns)  ·  ` +
      `U = ${Number.isFinite(ff.energy) ? ff.energy.toFixed(1) : "NaN"} kcal/mol  ·  ` +
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
