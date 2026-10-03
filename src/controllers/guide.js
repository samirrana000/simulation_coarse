/**
 * controllers/guide.js — FP1 first-run guided checklist (Load -> Build -> Run
 * -> Analyze).
 *
 * PUBLIC API
 *   initGuide()      Paint the checklist once at startup.
 *   updateGuide()    Repaint the four ✓/○ steps from live app state. Called
 *                    directly on the click transitions (Load / Build / Run /
 *                    Rec+Stop) and from guideTick.
 *   guideTick(now)   ≤1 Hz poll wrapper, so steps/time/frame progress (which
 *                    have no hooks) still light up.
 *
 * NEEDS (imports): ui + state + recorder + fp1GuideState from ../ui.js.
 *
 * One collapsed subpanel in Structure (index.html); live ✓/○ driven by
 * existing state (parsed/built/steps+time/frames). Polled at ≤1 Hz from
 * tick (both branches) so Analyze/frames progress needs no new hooks;
 * direct calls cover the click transitions.
 */

import { ui, state, recorder, fp1GuideState } from "../ui.js";
import { ignore } from "../errors.js";

const GUIDE_LABELS = {
  load: "Load — fetch a PDB, drop a file, or 1-click sample",
  build: "Build — Build System (auto after load)",
  run: "Run — ▶ Run advances the simulation",
  analyze: "Analyze — ● Rec + Run, Stop, then Analyze Trajectory",
};

/** Snapshot the app state the checklist state machine consumes. */
function readGuideInput() {
  return {
    hasPdb: !!(state.pdbText || state.parsed || state.parsedHeavy),
    hasBuild: !!(state.ff && state.integ),
    hasRun: !!state.running,
    steps: state.integ ? (state.integ.steps ?? 0) : 0,
    time: state.integ ? (state.integ.time ?? 0) : 0,
    nFrames: recorder ? (recorder.count ?? 0) : 0,
  };
}

/** Repaint the four checklist steps from live app state. */
export function updateGuide() {
  let s;
  try { s = fp1GuideState(readGuideInput()); } catch (_) { return; }
  try {
    const paint = (el, done, key) => { if (el) el.textContent = `${done ? "✓" : "○"} ${GUIDE_LABELS[key]}`; };
    paint(ui.guideStepLoad, s.load, "load");
    paint(ui.guideStepBuild, s.build, "build");
    paint(ui.guideStepRun, s.run, "run");
    paint(ui.guideStepAnalyze, s.analyze, "analyze");
    // Guarded offer: emphasize the 1-click sample only while empty.
    if (ui.sampleBtn) ui.sampleBtn.style.outline = s.load ? "" : "1px solid #38bdf8";
} catch (e) { ignore(e, "guide paint@initGuide", "checklist DOM absent headless; state itself is computed above and stays correct"); }
}

let _lastGuideTick = 0;

/** ≤1 Hz poll of updateGuide. */
export function guideTick(now) {
  const t = now ?? 0;
  if (t - _lastGuideTick < 1000) return;
  _lastGuideTick = t;
  updateGuide();
}

/** Paint the checklist once at startup (headless-safe: paints when DOM exists). */
export function initGuide() {
  try { updateGuide(); } catch (e) { ignore(e, "updateGuide@guideTick", "checklist DOM absent headless; a ≤1 Hz poll must not kill the loop"); }
}