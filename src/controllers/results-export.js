/**
 * controllers/results-export.js — the #resultsDlBtn export: build a
 * `results-record-v1` object from the live app state and hand it to the
 * browser's download machinery.
 *
 * PUBLIC API
 *   initResultsExport()                  Wire #resultsDlBtn.
 *   collectRunSnapshot()                 The provenance/input/health snapshot
 *                                        handed to buildResultsRecord.
 *   resultSpecs()                        The value rows, each already carrying
 *                                        the uncertainty its analysis module
 *                                        actually computed.
 *
 * NEEDS (imports): ui + state + recorder from ../ui.js; analyzeTrajectory from
 * ../analysis.js; trackPocketVolume + detectCryptic from
 * ../analysis/cryptic_pockets.js; buildResultsRecord + serializeResultsRecord +
 * resultsCsv + resultsExportGate + resultsFilename + RESULTS_VALUE_IDS from
 * ../results-record.js; errorSummary from ../errors.js; downloadText from
 * ../recorder.js; settingsState from ../settings-panel.js; ignore from
 * ../errors.js.
 *
 * WHY THIS IS A CONTROLLER AND NOT PART OF analysis-panel.js
 * ---------------------------------------------------------
 * The analysis report panel already computes most of these numbers and caches
 * the rendered text, but src/analysis-panel.js is pinned at 617 LOC by
 * tests/test_module_size.js's ratchet (it may shrink, never grow), so the export
 * cannot live there. This module re-derives the values from the same analysis
 * functions the panel uses — the numbers are the same function's output, not a
 * reimplementation — and the exported `notIncluded` list names, per family,
 * exactly what this path does NOT carry so a reader is never left guessing.
 *
 * ZERO PHYSICS HERE
 * -----------------
 * Nothing in this file touches a kernel, a force field or an integrator step. It
 * reads state, calls existing analysis functions, and writes a file.
 *
 * WHY THE GATE RUNS BEFORE THE DOWNLOAD
 * -------------------------------------
 * A record is only useful if it is true. resultsExportGate() refuses when the
 * run has recorded errors, when a value has no trust-boundary statement, or
 * when the input bytes could not be hashed. The refusal message names the
 * remedy instead of printing "⚠ error", because a user who cannot tell what to
 * do next will just reload and try again.
 */

import { ui, state, recorder } from "../ui.js";
import { analyzeTrajectory } from "../analysis.js";
import { trackPocketVolume, detectCryptic } from "../analysis/cryptic_pockets.js";
import {
  buildResultsRecord, serializeResultsRecord, resultsCsv,
  resultsExportGate, resultsFilename, RESULTS_VALUE_IDS,
} from "../results-record.js";
import { errorSummary, ignore } from "../errors.js";
import { downloadText } from "../recorder.js";
import { settingsState } from "../settings-panel.js";

/** Pocket residues for the volume track: the funnel pocket, else ≤8 Å of the ligand COM. */
const CRYPTIC_RCUT = 8.0;
const CRYPTIC_MAX_RES = 40;
const CRYPTIC_KSIGMA = 1.5;

/** Cache of the last trajectory analysis, keyed on the inputs that define it. */
let _lastReport = null;
let _lastKey = "";

/** Where the numbers below came from — reused by docs and by the record itself. */
const RESULT_SOURCES = {
  dg_bind: "src/analysis.js [5] / src/funnel.js estimateDG + convergenceSE",
  b_factor_pearson_r: "src/analysis.js [1] pearson(B_sim, B_exp)",
  b_factor_mean_sim_ang2: "src/analysis.js [1] mean B_sim = (8π²/3)⟨Δr²⟩",
  rmsf_mean_ang: "src/analysis.js [1] <RMSF>",
  var_top_k: "src/analysis.js [2] top-k variance fraction",
  rmsip_pca_enm: "src/analysis.js [2] RMSIP(PCA, ENM soft modes)",
  occupancy_bound_fraction: "src/analysis.js [3] boundFrac",
  contact_lifetime_mean_ps: "src/analysis.js [4] meanLifetimePs",
  pocket_volume_mean_ang3: "src/analysis/cryptic_pockets.js detectCryptic mean",
  pocket_volume_sd_ang3: "src/analysis/cryptic_pockets.js detectCryptic sd",
};

/** Families this export path does not carry, and exactly why. */
const NOT_INCLUDED = [
  {
    family: "ΔH/ΔS decomposition (thermo_ligand / thermo_sasa / thermodynamics)",
    why: "those values come from a separate multi-leg run (a chunked apo relaxation plus an SASA sweep) whose result is cached as rendered text in the analysis panel, not as a numeric object this controller can read without growing a module the size ratchet has pinned",
    howToGet: "run ΔH/ΔS Decomposition, then Export Thermo (TXT); the ± error bars are printed in that table",
  },
  {
    family: "alanine scanning, SMD/Jarzynski pulls, DCCM",
    why: "each is a separate button with its own output artifact and none of them writes a value this controller can read",
    howToGet: "the corresponding export button in the analysis panel (Ala-Scan / SMD Pulls / DCCM CSV)",
  },
  {
    family: "trajectory coordinates",
    why: "a results record is a summary with provenance; frames are a multi-MB artifact and would break the record's size contract",
    howToGet: "Download (XYZ/PDB/JSON) in the Recording panel",
  },
];

/** Write a message to the two caption surfaces the Recording panel already owns. */
function say(text) {
  try { if (ui.recStatus) ui.recStatus.textContent = text; } catch (e) { ignore(e, "recStatus@resultsExport", "#recStatus is null in a headless import; the export still runs"); }
  try { if (ui.hud) ui.hud.textContent = text; } catch (e) { ignore(e, "hud@resultsExport", "#hud is null in a headless import"); }
}

/**
 * Provenance + input identity + run health, read from the live app.
 *
 * Every value here comes from a module that OWNS it — the integrator for dt /
 * steps / seed, the force field for the cutoff and the system counts,
 * settingsState for the solvent and RESPA settings, src/errors.js for the
 * failure counters. Nothing is re-derived from the DOM sliders here, because the
 * force field on screen is the authority on what the run actually used.
 * @returns {object} snapshot for buildResultsRecord
 */
export function collectRunSnapshot() {
  const ff = state.ff;
  const integ = state.integ;
  const funnel = state.funnel;
  const src = state.inputSource || {};
  const es = errorSummary();
  return {
    input: {
      pdbId: src.pdbId || (ui.pdbId?.value ?? "").trim(),
      fileName: src.fileName || (ui.pdbId?.value ? `${(ui.pdbId.value || "").trim().toLowerCase()}.pdb` : ""),
      origin: src.origin || "unknown",
      pdbText: state.pdbText,
    },
    selection: {
      chains: ui.chainsInput?.value ?? "",
      resFrom: ui.resFrom?.value === "" || ui.resFrom?.value == null ? null : Number(ui.resFrom.value),
      resTo: ui.resTo?.value === "" || ui.resTo?.value == null ? null : Number(ui.resTo.value),
      includeLigand: ui.includeLig?.checked ?? true,
    },
    physics: {
      modelMode: state.heavyMode ? "heavy" : "cg",
      physicsLevel: ff?.physicsLevel ?? settingsState.physicsLevel ?? "L0",
      tierFlags: {
        charges: ff?.chargesOn === true,
        hbMode: String(ff?.hbMode ?? "off"),
        weak: String(ff?.weak ?? "off"),
        bindLog: !!state.bindLog && state.bindLog.nEvents > 0,
      },
      temperatureK: integ?.T ?? Number(ui.temp?.value ?? 300),
      frictionPerPs: integ?.zeta ?? Number(ui.fric?.value ?? 8),
      beadMassDa: Number(ui.mass?.value ?? 110),
      dtPs: integ?.dt ?? null,
      cutoffAng: ff?.rc ?? Number(ui.rc?.value ?? 10),
      springGamma: ff?.gamma ?? Number(ui.gamma?.value ?? 2),
      solventModel: settingsState.solventModel,
      saltM: settingsState.saltM,
      epsIn: settingsState.epsIn,
      epsOut: settingsState.epsOut,
      sasaGamma: settingsState.sasaGamma,
      respaOn: settingsState.respaOn === true,
      respaOuterFs: settingsState.respaOuterFs,
      nParticles: ff?.n ?? 0,
      nProtein: ff?.nProt ?? 0,
      nLigandAtoms: ff?.nLigAtoms ?? 0,
      nLigandMolecules: (state.ligands ?? []).length,
      motionGain: Number(ui.motionGain?.value ?? 1),
    },
    run: {
      // getSeed() returns null on the default Math.random path. That null is
      // carried through to the record rather than replaced with a 0, because
      // "seed 0" would read as a reproducible seed and this run is not.
      seed: typeof integ?.getSeed === "function" ? integ.getSeed() : null,
      simulatedPs: integ?.time ?? 0,
      steps: integ?.steps ?? 0,
      dtPs: integ?.dt ?? null,
    },
    trajectory: {
      nFrames: recorder.count,
      spanPs: recorder.times.length > 1 ? recorder.times[recorder.times.length - 1] - recorder.times[0] : 0,
      stridePs: recorder.stridePs,
    },
    runHealth: {
      errors: es.errors,
      ignored: es.ignored,
      lastError: es.last,
    },
    coverageNotIncluded: NOT_INCLUDED,
  };
}

/** Pocket residue indices (bead space) for the volume track, or null. */
function pocketIndices() {
  try {
    const f = state.funnel;
    if (f && f.active && Array.isArray(f.pocket) && f.pocket.length) return [...f.pocket];
    const ff = state.ff;
    if (!ff || !ff.nLigAtoms) return null;
    const nProt = ff.nProt, lig0 = ff.ligandStart ?? nProt;
    let cx = 0, cy = 0, cz = 0;
    for (let a = 0; a < ff.nLigAtoms; a++) { cx += ff.ref[3 * (lig0 + a)]; cy += ff.ref[3 * (lig0 + a) + 1]; cz += ff.ref[3 * (lig0 + a) + 2]; }
    cx /= ff.nLigAtoms; cy /= ff.nLigAtoms; cz /= ff.nLigAtoms;
    const r2 = CRYPTIC_RCUT * CRYPTIC_RCUT;
    const idx = [];
    for (let i = 0; i < nProt; i++) {
      const dx = ff.ref[3 * i] - cx, dy = ff.ref[3 * i + 1] - cy, dz = ff.ref[3 * i + 2] - cz;
      if (dx * dx + dy * dy + dz * dz <= r2) idx.push(i);
    }
    return idx.length ? idx.slice(0, CRYPTIC_MAX_RES) : null;
  } catch (e) {
    ignore(e, "pocketIndices@resultsExport", "a pocket-set failure costs the volume rows, not the rest of the record");
    return null;
  }
}

/**
 * The value rows. Each carries the uncertainty its own analysis module computed
 * — never one invented here — and an explicit reason when the module computed
 * none. `buildResultsRecord` drops a row whose value is not finite rather than
 * writing a placeholder.
 *
 * @param {object} rep analyzeTrajectory report
 * @param {object} [det] detectCryptic result (optional)
 * @returns {object[]} specs for buildResultsRecord
 */
export function resultSpecs(rep, det = null) {
  const out = [];
  const push = (id, value, unit, extra = {}) => {
    if (Number.isFinite(Number(value))) out.push({ id, value: Number(value), unit, ...extra });
  };
  const nF = rep?.nFrames ?? 0;

  // B-factors. The correlation r is a single statistic of the recorded window;
  // the analysis code computes no interval on it, so the record says so.
  if (rep?.bf) {
    push("b_factor_pearson_r", rep.bf.r, "1 (Pearson)", {
      label: `Pearson r(B_sim, B_exp) over ${rep.bf.n} residues`,
      n: rep.bf.n,
      noUncertaintyReason: "pearson() returns a point estimate with no interval; the analysis code computes no confidence band on r",
    });
    push("b_factor_mean_sim_ang2", rep.bf.meanBsim, "A^2", {
      label: "mean simulated B-factor",
      n: nF,
      noUncertaintyReason: "no error bar on the mean B is computed — the value scales directly with the recorded window length",
    });
  }
  if (rep?.lines) {
    const m = /<RMSF> = ([\d.]+) A/.exec(rep.lines.join("\n"));
    if (m) push("rmsf_mean_ang", Number(m[1]), "A", {
      label: "mean RMSF over protein beads",
      n: nF,
      noUncertaintyReason: "no error bar on the mean RMSF is computed",
    });
  }
  if (rep?.rmsip) {
    push("var_top_k", rep.rmsip.varTopK, "fraction of variance", {
      label: `variance in the top ${rep.rmsip.k} principal modes`,
      n: nF,
      noUncertaintyReason: "no error bar is computed, and an under-sampled window inflates this by construction",
    });
    // rmsip is null when the system has no elastic network (heavy mode), which
    // is absence rather than zero: no row is emitted.
    push("rmsip_pca_enm", rep.rmsip.rmsip, "1 (RMSIP)", {
      label: `RMSIP(PCA, ${rep.rmsip.k} ENM soft modes)`,
      n: nF,
      noUncertaintyReason: "RMSIP is a deterministic inner product of two mode bases; the analysis code computes no error bar for it",
    });
  }
  if (rep?.occupancy) {
    push("occupancy_bound_fraction", rep.occupancy.boundFrac, "fraction of frames", {
      label: `ligand COM inside the pocket (r ≤ ${rep.occupancy.rPocket} A)`,
      n: nF,
      noUncertaintyReason: "no binomial or bootstrapped error bar is computed on the bound-frame fraction",
    });
  }
  if (rep?.contacts) {
    push("contact_lifetime_mean_ps", rep.contacts.meanLifetimePs, "ps", {
      label: `mean longest contact streak (< 6 A, ${rep.contacts.nPairs} pairs)`,
      n: nF,
      noUncertaintyReason: "no standard error is computed over the contact set; the streak is stride-limited",
    });
  }
  if (det) {
    // The track SD is a DISPERSION, not an error bar on the mean. It travels in
    // its own field so it can never be read as one.
    push("pocket_volume_mean_ang3", det.mean, "A^3", {
      label: "mean pocket volume over recorded frames",
      n: rep?.nFrames ?? 0,
      noUncertaintyReason: "detectCryptic reports a sample SD of the per-frame track; that is a dispersion and no standard error on the mean is computed",
      dispersion: { value: det.sd, kind: "sample-sd-of-track" },
    });
    push("pocket_volume_sd_ang3", det.sd, "A^3", {
      label: "sample SD of the per-frame pocket-volume track",
      n: rep?.nFrames ?? 0,
      noUncertaintyReason: "this IS the dispersion; reporting an error bar on it would be a second-order statistic the code never computes",
    });
  }
  // ΔG from the funnel. convergenceSE() returns null below 50 hills, which
  // becomes "no error bar" rather than a 0.
  const funnel = state.funnel;
  if (funnel && funnel.active && funnel._nHills > 0) {
    const dg = funnel.estimateDG();
    const se = typeof funnel.convergenceSE === "function" ? funnel.convergenceSE() : null;
    const nHills = funnel._nHills;
    push("dg_bind", Number.isNaN(dg) ? null : dg, "kcal/mol", {
      label: `binding dG from the WTM funnel bias (${nHills} hills, ${funnel.biasFactor ?? "?"} bias factor)`,
      n: nHills,
      uncertainty: { value: se, kind: "hill-count-convergence", method: "kB*T/sqrt(nHills)" },
      noUncertaintyReason: `convergenceSE() returns null below 50 deposited hills (this run has ${nHills}), so no convergence error bar exists yet`,
    });
  }
  return out;
}

/** Run (or reuse) the trajectory analysis and the pocket-volume detector. */
function computeAnalysis() {
  const key = `${recorder.count}:${state.ff ? state.ff.n : 0}:${state.heavyMode ? "h" : "c"}`;
  if (_lastReport && _lastKey === key) return { rep: _lastReport, det: _lastReport?.cryptic ?? null };
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
  let det = null;
  const pocket = pocketIndices();
  if (pocket && pocket.length) {
    try {
      const frames = recorder.count >= 2 ? recorder.frames : [state.integ.pos];
      det = detectCryptic(trackPocketVolume(frames, { pocketIndices: pocket }), {
        kSigma: CRYPTIC_KSIGMA, frames: null, pocketIndices: pocket,
      });
    } catch (e) {
      ignore(e, "detectCryptic@resultsExport", "the volume rows are optional; the rest of the record is unaffected");
      det = null;
    }
  }
  rep.cryptic = det;
  _lastReport = rep;
  _lastKey = key;
  return { rep, det };
}

/** Drop the analysis cache (called on rebuild / new recording). */
export function invalidateResultsCache() {
  _lastReport = null;
  _lastKey = "";
}

/**
 * Build the record and hand it to the caller. Exported separately from the
 * button handler so the whole path is exercisable headlessly.
 * @returns {{ok:boolean, record?:object, text?:string, csv?:string, gate:object}}
 */
export function exportResultsRecord() {
  const snap = collectRunSnapshot();
  let computed;
  try {
    computed = computeAnalysis();
  } catch (err) {
    return { ok: false, gate: { ok: false, reason: `analysis failed: ${err && err.message ? err.message : err}`, remedy: "record a longer trajectory, or check #analysisOut", warnings: [] } };
  }
  snap.results = resultSpecs(computed.rep, computed.det);
  const record = buildResultsRecord(snap);
  const gate = resultsExportGate(record);
  if (!gate.ok) return { ok: false, record, gate };
  return { ok: true, record, text: serializeResultsRecord(record), csv: resultsCsv(record), gate };
}

/** Wire the Recording panel's results-record button. Call once at startup. */
export function initResultsExport() {
  const btn = ui.resultsDlBtn;
  if (!btn) return;
  btn.addEventListener("click", () => {
    if (!state.pdbText) { say("⚠ Load a structure first — the record hashes the input bytes, so it cannot exist without them."); return; }
    if (!state.ff || !state.sel) { say("⚠ Build a system first (Load → Build)."); return; }
    if (recorder.count < 1) { say("⚠ Nothing recorded yet — press ● Rec, run the sim, then Stop."); return; }
    say("Building the results record… (RMSIP may take a few seconds)");
    btn.disabled = true;
    setTimeout(() => {
      try {
        const out = exportResultsRecord();
        if (!out.ok) { say(`⚠ Results record NOT exported — ${out.gate.reason}. ${out.gate.remedy}`); return; }
        const name = resultsFilename(out.record);
        downloadText(out.text, name);
        const n = out.record.results.length;
        const barred = out.record.results.filter((r) => r.uncertaintyAvailable).length;
        const warn = out.gate.warnings.length ? ` ⚠ ${out.gate.warnings.join(" · ")}` : "";
        say(`${name} — ${out.text.length} B · ${n} value(s), ${barred} with an error bar · ${out.record.runHealth.errors} error(s), ${out.record.runHealth.documentedNoOps} documented no-op(s).${warn}`);
      } catch (err) {
        say("⚠ " + (err && err.message ? err.message : String(err)));
      } finally {
        try { btn.disabled = false; } catch (e) { ignore(e, "resultsDlBtn@finally", "button state only; a stuck button is recoverable by reload"); }
      }
    }, 20);
  });
}

/** Value ids this controller can emit — declared next to the emitter. */
export { RESULTS_VALUE_IDS, RESULT_SOURCES };
