/**
 * ml-tier.js — ML tier panel: NN contact map (ESM prior) + NN pose scorer
 * (item 5 modularization: moved verbatim from main.js; only the imports plus
 * the local `state`/`ui` handles changed — those now come from ui.js).
 */

import { PoseScorer } from "./scorer.js?v=9";
import { ui, state, updateSelSummary } from "./ui.js?v=9";

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
export function applyMLToFF() {
  if (!state.ff) return;
  if (ui.nnContacts.checked && state.contacts) {
    state.ff.setSpringScale(state.contacts, state.mlAlpha);
  } else {
    state.ff.clearSpringScale();
  }
  if (ui.poseScore.checked && !state.scorer) state.scorer = PoseScorer.docked();
  updateSelSummary();
}
