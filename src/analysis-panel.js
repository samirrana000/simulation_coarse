/**
 * analysis-panel.js — panel 6: trajectory analysis report + PMF CSV download
 * (item 5 modularization: moved verbatim from main.js; only the imports
 * changed — `ui` / `state` / `recorder` now come from ui.js).
 */

import { analyzeTrajectory, pmfCsv } from "./analysis.js?v=8";
import { downloadText } from "./recorder.js?v=8";
import { ui, state, recorder } from "./ui.js?v=8";

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
