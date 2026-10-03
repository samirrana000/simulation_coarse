/**
 * controllers/recording.js — the Recording panel: recorder transport,
 * trajectory + BindLog export, and session save/restore.
 *
 * PUBLIC API
 *   initRecording()   Wire #recBtn, #recStopBtn, #bindlogOn, #dlBtn,
 *                     #bindlogDlBtn, #sessSaveBtn and #sessFile.
 *   collectSessionSnapshot() -> the settings/picker/counts object handed to
 *                     session.js buildSession (never frames or PDB text).
 *   applySession(sess) Apply a validated session object to the settings,
 *                     pickers and captions (pure DOM; never throws).
 *
 * NEEDS (imports): ui + state + recorder + updateRecStatus from ../ui.js;
 * BUILD_DATE from ../version.js; downloadText from ../recorder.js;
 * buildSession + serializeSession + parseSession + trajectoryJson +
 * downloadBlob from ../session.js; formatInputError from ../input_errors.js;
 * settingsState + persistPhysicsLevel from ../settings-panel.js;
 * applyLiveTrackTerms from ./live-terms.js; bindLogWanted from
 * ./physics-tier.js; onParamChange from ./params.js; updateGuide from
 * ./guide.js.
 *
 * Session files carry settings/picker/counts; full frames + PDB text are
 * never persisted (counts only, see ../session.js).
 */

import { ui, state, recorder } from "../ui.js";
import { BUILD_DATE } from "../version.js";
import { downloadText } from "../recorder.js";
import { buildSession, serializeSession, parseSession, trajectoryJson, downloadBlob } from "../session.js";
import { formatInputError } from "../input_errors.js";
import { ignore, recordError } from "../errors.js";
import { settingsState, persistPhysicsLevel } from "../settings-panel.js";
import { applyLiveTrackTerms } from "./live-terms.js";
import { bindLogWanted } from "./physics-tier.js";
import { onParamChange } from "./param-binding.js";
import { updateGuide } from "./guide.js";

/** Settings/picker/counts snapshot handed to session.js buildSession. */
export function collectSessionSnapshot() {
  return {
    pdbId: (ui.pdbId?.value ?? "").trim(),
    modelMode: ui.modelMode?.value ?? "cg",
    chains: ui.chainsInput?.value ?? "",
    resFrom: ui.resFrom?.value === "" || ui.resFrom?.value == null ? null : Number(ui.resFrom.value),
    resTo: ui.resTo?.value === "" || ui.resTo?.value == null ? null : Number(ui.resTo.value),
    includeLig: ui.includeLig?.checked ?? true,
    physicsLevel: ui.physicsLevel?.value ?? settingsState.physicsLevel ?? "L0",
    ligand: { selected: ui.ligSelect?.value ?? "" },
    thermoLig: ui.thermoLig?.value ?? "auto",
    settings: { ...settingsState },
    dynamics: {
      rc: Number(ui.rc?.value ?? 10), gamma: Number(ui.gamma?.value ?? 2),
      temp: Number(ui.temp?.value ?? 300), fric: Number(ui.fric?.value ?? 8),
      mass: Number(ui.mass?.value ?? 110), motionGain: Number(ui.motionGain?.value ?? 1.3),
      bindPot: ui.bindPot?.checked ?? true, holoSprings: ui.holoSprings?.checked ?? true,
    },
    recording: {
      stridePs: Number(ui.stridePs?.value ?? 2), maxFrames: Number(ui.maxFrames?.value ?? 500),
      exportFmt: ui.exportFmt?.value ?? "xyz",
    },
    recorderMeta: { count: recorder.count, spanPs: recorder.times.length > 1 ? recorder.times[recorder.times.length - 1] - recorder.times[0] : 0 },
  };
}

/** Apply a validated session object to settings + pickers + caption (pure-DOM). */
export function applySession(sess) {
  try {
    if (typeof sess.pdbId === "string" && ui.pdbId) ui.pdbId.value = sess.pdbId;
    if (sess.modelMode && ui.modelMode) ui.modelMode.value = sess.modelMode;
    if (typeof sess.chains === "string" && ui.chainsInput) ui.chainsInput.value = sess.chains;
    if (ui.resFrom) ui.resFrom.value = sess.resFrom ?? "";
    if (ui.resTo) ui.resTo.value = sess.resTo ?? "";
    if (typeof sess.includeLig === "boolean" && ui.includeLig) ui.includeLig.checked = sess.includeLig;
    if (sess.physicsLevel && ui.physicsLevel) {
      ui.physicsLevel.value = sess.physicsLevel;
      try { persistPhysicsLevel(sess.physicsLevel); } catch (e) { ignore(e, "persistPhysicsLevel@applySession(sel)", "localStorage may be unavailable (private mode); the in-memory tier is already applied"); }
    } else if (sess.physicsLevel) {
      try { persistPhysicsLevel(sess.physicsLevel); } catch (e) { ignore(e, "persistPhysicsLevel@applySession(noSel)", "localStorage may be unavailable (private mode); persistence is best-effort"); }
    }
    if (sess.ligand && typeof sess.ligand.selected === "string" && ui.ligSelect) {
      try {
        const opt = [...ui.ligSelect.options].find((o) => o.value === sess.ligand.selected);
        if (opt) ui.ligSelect.value = sess.ligand.selected;
      } catch (e) { ignore(e, "ligSelect@applySession", "<select>.options is unavailable when ui.js captured no element (headless)"); }
    }
    if (typeof sess.thermoLig === "string" && ui.thermoLig) {
      try {
        const opt = [...ui.thermoLig.options].find((o) => o.value === sess.thermoLig);
        ui.thermoLig.value = opt ? sess.thermoLig : "auto";
      } catch (e) { ignore(e, "thermoLig@applySession", "picker is optional markup; thermoLigValue() falls back to auto"); }
    }
    const st = sess.settings && typeof sess.settings === "object" ? sess.settings : {};
    for (const k of ["backend", "solventModel", "saltM", "epsIn", "epsOut", "sasaGamma", "numThreads", "respaOn", "respaOuterFs", "chemicalNetworkOn"]) {
      if (st[k] !== undefined) {
        try { settingsState[k] = st[k]; } catch (e) { ignore(e, `settingsState[${k}]@applySession`, "settingsState is a plain object; only reachable via a frozen import under a bundler"); }
      }
    }
    try {
      const setVal = (id, v) => { const el = document.getElementById(id); if (el && v !== undefined) el.value = String(v); };
      const setChk = (id, v) => { const el = document.getElementById(id); if (el && typeof v === "boolean") el.checked = v; };
      setVal("backendSelect", settingsState.backend);
      setVal("threadsInput", settingsState.numThreads);
      setVal("solventSelect", settingsState.solventModel);
      setVal("saltInput", (settingsState.saltM * 1000).toFixed(0));
      setVal("epsInInput", settingsState.epsIn);
      setVal("epsOutInput", settingsState.epsOut);
      setVal("sasaInput", settingsState.sasaGamma);
      setChk("netToggleModal", settingsState.chemicalNetworkOn);
      setChk("respaToggle", settingsState.respaOn);
      setVal("respaOuter", settingsState.respaOuterFs);
} catch (e) { ignore(e, "panel restore@applySession", "individual UI controls are optional markup; each block restores independently"); }
    const dyn = sess.dynamics && typeof sess.dynamics === "object" ? sess.dynamics : {};
    try {
      const setPair = (id, numId, lblId, v) => {
        if (!Number.isFinite(Number(v))) return;
        const el = document.getElementById(id), num = numId ? document.getElementById(numId) : null;
        const lbl = lblId ? document.getElementById(lblId) : null;
        if (el) el.value = String(v);
        if (num) num.value = String(v);
        if (lbl) lbl.textContent = String(v);
      };
      setPair("rc", "rcNum", "v_rc", dyn.rc);
      setPair("gamma", "gammaNum", "v_gamma", dyn.gamma);
      setPair("temp", "tempNum", "v_temp", dyn.temp);
      setPair("fric", "fricNum", "v_fric", dyn.fric);
      setPair("mass", "massNum", "v_mass", dyn.mass);
      setPair("motionGain", "motionGainNum", "v_motionGain", dyn.motionGain);
      if (typeof dyn.bindPot === "boolean" && ui.bindPot) ui.bindPot.checked = dyn.bindPot;
      if (typeof dyn.holoSprings === "boolean" && ui.holoSprings) ui.holoSprings.checked = dyn.holoSprings;
} catch (e) { ignore(e, "recorder fields@applySession", "individual UI controls are optional markup; each block restores independently"); }
    const rec = sess.recording && typeof sess.recording === "object" ? sess.recording : {};
    try {
      if (Number.isFinite(Number(rec.stridePs)) && ui.stridePs) ui.stridePs.value = String(rec.stridePs);
      if (Number.isFinite(Number(rec.maxFrames)) && ui.maxFrames) ui.maxFrames.value = String(rec.maxFrames);
      if (typeof rec.exportFmt === "string" && ui.exportFmt) {
        try {
          const opt = [...ui.exportFmt.options].find((o) => o.value === rec.exportFmt);
          if (opt) ui.exportFmt.value = rec.exportFmt;
        } catch (e) { ignore(e, "exportFmt@applySession", "<select>.options is unavailable when ui.js captured no element (headless)"); }
      }
} catch (e) { ignore(e, "recorder fields 2@applySession", "individual UI controls are optional markup; each block restores independently"); }
    try { onParamChange(false); } catch (e) { ignore(e, "onParamChange@applySession", "a freshly loaded session has no live system to re-parameterise yet"); }
    const n = sess.recorderMeta && Number.isFinite(Number(sess.recorderMeta.count)) ? Number(sess.recorderMeta.count) : 0;
    if (ui.canvasCaption) {
      ui.canvasCaption.textContent = `Session loaded (${sess.pdbId || "custom"} · ${sess.physicsLevel || "L0"} · saved ${n} frame(s) in memory only — re-record after Build).`;
    }
    if (ui.hud) ui.hud.textContent = `Session loaded: ${sess.pdbId || "custom"} · physics ${sess.physicsLevel || "L0"} · picker + settings restored.`;
  } catch (e) {
    // RECORDED, and no longer a silent no-op. The old comment ("apply never
    // throws to the loader") was the assertion this whole exercise exists to
    // remove: applySession restores ~20 fields, and a failure part-way through
    // leaves a half-restored session that looks fully restored. The loader
    // still does not throw — the failure is now a number on the top bar
    // instead of nothing at all.
    recordError(e, "applySession@recording");
  }
}

/** Wire the whole Recording panel. Call once at startup. */
export function initRecording() {
  if (ui.recBtn) {
    ui.recBtn.addEventListener("click", () => {
      if (!state.integ) return;
      recorder.start(state.integ.time, Number(ui.stridePs.value), Number(ui.maxFrames.value));
      ui.recBtn.classList.add("rec-on");
      try { updateGuide(); } catch (e) { ignore(e, "updateGuide@Analyze click", "checklist DOM absent headless"); } // FP1: Analyze progress
    });
  }
  if (ui.recStopBtn) {
    ui.recStopBtn.addEventListener("click", () => {
      recorder.stop();
      ui.recBtn.classList.remove("rec-on");
      try { updateGuide(); } catch (e) { ignore(e, "updateGuide@frame landed", "checklist DOM absent headless"); } // FP1: frames landed → Analyze ✓
    });
  }
  // Loop-2 S4 (R6 §5) + Rev1/Issue4: BindLog capture toggle — flips the BindLog
  // event-capture path and reapplies the live trackTerms policy (accumulators
  // stay on while a ligand is present for the HUD mirror, so unchecking the box
  // never blanks the per-term HUD/energy stripe). Resets the contact-diff
  // baseline so form/break events only fire across the switch.
  if (ui.bindlogOn) {
    ui.bindlogOn.addEventListener("change", () => {
      // Loop-2 S7: effective flag is checkbox OR L2 tier (unchecking while L2
      // is active keeps capture on — the tier owns the flag until deselected).
      if (state.ff) applyLiveTrackTerms(state.ff);
      if (!bindLogWanted()) state._lastContacts = null;
      if (bindLogWanted() && state.funnel) {
        state.funnel.onHill = (cv, h) => state.bindLog && state.bindLog.pushHill(state.integ.time, cv, h);
      } else if (state.funnel) {
        state.funnel.onHill = null;
      }
    });
  }
  if (ui.dlBtn) {
    ui.dlBtn.addEventListener("click", () => {
      if (!state.sel) return;
      const fmt = ui.exportFmt.value;
      try {
        const prov = {
          T: Number(ui.temp?.value || 300),
          gamma: Number(ui.gamma?.value || 2),
          seed: (typeof state.seed !== "undefined" ? state.seed : 0),
          date: BUILD_DATE,
        };
        // FP4: JSON trajectory export (new exportFmt option; XYZ/PDB path untouched).
        if (fmt === "json") {
          const text = trajectoryJson(recorder.frames, recorder.times, prov);
          downloadText(text, `cg_traj_${recorder.count}frames.json`);
          return;
        }
        const text = recorder.buildFile(fmt, state.sel.beads, prov);
        const name = `cg_traj_${recorder.count}frames.${fmt}`;
        downloadText(text, name);
      } catch (err) {
        if (ui.recStatus) ui.recStatus.textContent = "⚠ " + err.message;
      }
    });
  }

  if (ui.bindlogDlBtn) {
    ui.bindlogDlBtn.addEventListener("click", () => {
      const bl = state.bindLog;
      if (!bl || bl.nFrames === 0) {
        if (ui.recStatus) ui.recStatus.textContent = "⚠ Nothing captured yet — enable BindLog capture, Run, then Export BindLog (BLG1).";
        return;
      }
      try {
        const buf = bl.toBinaryBlob();
        const ok = downloadBlob(buf, `bindlog_${bl.nFrames}f_${bl.nEvents}e.blg1`);
        if (ui.recStatus) {
          ui.recStatus.textContent = ok
            ? `${recorder.count} frames · BindLog BLG1 exported (${bl.nFrames} frames, ${bl.nEvents} events, ${(buf.byteLength / 1024).toFixed(1)} KiB).`
            : "⚠ BindLog download needs a browser (headless: use toBinaryBlob directly).";
        }
      } catch (err) {
        if (ui.recStatus) ui.recStatus.textContent = "⚠ " + (err?.message ?? String(err));
      }
    });
  }

  if (ui.sessSaveBtn) {
    ui.sessSaveBtn.addEventListener("click", () => {
      try {
        const text = serializeSession(buildSession(collectSessionSnapshot()));
        const id = ((ui.pdbId?.value ?? "").trim() || "custom").replace(/\W+/g, "_");
        downloadText(text, `session_${id}_v1.json`);
        if (ui.hud) ui.hud.textContent = `Session saved (${text.length} B, counts only — frames stay in memory).`;
      } catch (err) {
        if (ui.hud) ui.hud.textContent = formatInputError(err);
      }
    });
  }

  if (ui.sessFile) {
    ui.sessFile.addEventListener("change", async () => {
      const f = ui.sessFile.files && ui.sessFile.files[0];
      if (!f) return;
      try {
        const text = await f.text();
        const parsed = parseSession(text);
        if (!parsed.ok) {
          if (ui.hud) ui.hud.textContent = formatInputError(parsed.error);
          if (ui.recStatus) ui.recStatus.textContent = formatInputError(parsed.error);
          return;
        }
        applySession(parsed.data);
      } catch (err) {
        if (ui.hud) ui.hud.textContent = formatInputError(err);
      } finally {
        try { ui.sessFile.value = ""; } catch (e) { ignore(e, "sessFile@fileInput finally", "clearing the file input is cosmetic"); } // allow re-loading the same file
      }
    });
  }
}