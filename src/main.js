/**
 * main.js — APPLICATION COMPOSITION ROOT.
 *
 * This file wires the app together and owns nothing else. It is the single
 * module index.html loads, so it is the single point whose parse failure takes
 * down the ES-module graph (AUDIT_REPORT.md §1.1 / §2.1) — which is exactly
 * why the implementation behind it was split out by responsibility. Read this
 * file as the table of contents of the application:
 *
 *   MODULE                            RESPONSIBILITY
 *   -------------------------------   ---------------------------------------
 *   ui.js                             DOM handles + shared state singletons
 *   viewer.js                         Canvas rendering (via ui.js)
 *   pdb.js                            fetch/parse/select a structure
 *   forcefield.js                     Cα energy function
 *   heavy.js                          heavy-atom force field + GB/SA
 *   integrator.js                     BAOAB Langevin dynamics
 *   recorder.js                       trajectory capture & download
 *   ml-tier.js                        NN contact map + pose scorer panel
 *   pmf-panel.js                      funnel PMF plot + reset
 *   analysis-panel.js                 trajectory analysis + PMF/DCCM/thermo
 *   settings-panel.js                 System Settings + GPU/worker accel
 *   network-panel.js                  Chemical Network Model & kinetics
 *   ligand-panel.js                   Ligand library & clash-free placement
 *   controllers/physics-tier.js       physics fidelity tier + capture flag
 *   controllers/param-binding.js             parameter hot reload + slider readouts
 *   controllers/live-terms.js         live per-term energy mirror
 *   controllers/dock.js               dock sparklines + trajectory scrub
 *   controllers/binding-insights.js            Binding-Insights canvases (≤1 Hz)
 *   controllers/accelerate.js         r-RESPA + worker-pool force backends
 *   controllers/structure-input.js    PDB/sample/file/MOL2 acquisition
 *   controllers/system-build.js       selection → force field → integrator
 *   controllers/transport.js          run/reset/display toggles + hotkeys
 *   controllers/recording.js          recorder, exports, session save/restore
 *   controllers/guide.js              FP1 onboarding checklist
 *   controllers/tick.js               the requestAnimationFrame loop
 *
 * THE ONE INJECTION SEAM
 * ----------------------
 * `physicsLevelSpec()` and the two force-field parameter literals live here
 * rather than in controllers/param-binding.js, because both build paths (the cold
 * build and the hot rebuild) must carry the active tier and the literals must
 * never drift apart. Controllers that need them cannot import them from here
 * without closing an import cycle, so the composition root passes them in —
 * which is the textbook role of a composition root anyway.
 *
 * Order below matters and is the original startup order: modals/panels first,
 * then the per-subsystem wiring, then the animation loop last.
 */

import { VERSION, BUILD_DATE } from "./version.js";
import { ui } from "./ui.js";
import "./analysis-panel.js"; // side-effect: Analyze + Phase-4 workflow buttons
import { initSettingsModal, restorePhysicsLevelSelect } from "./settings-panel.js";
import { initNetworkPanel } from "./network-panel.js";
import { PHYSICS_LEVELS, currentLevel } from "./controllers/physics-tier.js";
import { initParamBinding } from "./controllers/param-binding.js";
import { initDock } from "./controllers/dock.js";
import { initBindviz } from "./controllers/binding-insights.js";
import { updateGuide, initGuide } from "./controllers/guide.js";
import { initStructureInput } from "./controllers/structure-input.js";
import { initSystemBuild, buildSystem } from "./controllers/system-build.js";
import { initTransport } from "./controllers/transport.js";
import { initRecording } from "./controllers/recording.js";
import { startTick } from "./controllers/tick.js";

/* ------------------------------------------------------------------ */
/*  Cross-cutting seams (see "THE ONE INJECTION SEAM" above)            */
/* ------------------------------------------------------------------ */

/**
 * Read the active physics tier and persist it to settingsState.physicsLevel.
 * CG consumes charges/hbMode, HeavyForceField mirrors physicsLevel/charges/
 * hbMode queryably (Rev3/Issue1, kernels unchanged), and both build paths
 * carry `level` into `par` so describePhysics() reports the UI tier (rev2/issue1).
 * @returns {{level:string, charges:boolean, hbMode:string, weak:string, bindLog:boolean}} active tier spec
 */
export function physicsLevelSpec() {
  let lvl = "L0";
  try { lvl = currentLevel(); } catch (_) { lvl = "L0"; }
  return { level: lvl, ...PHYSICS_LEVELS[lvl] };
}

/**
 * Cold-build force-field parameter literal (buildSystem path).
 * @param {{level:string,charges:boolean,hbMode:string,weak:string}} physLvl active tier spec
 */
export function ffParams(physLvl) {
  return {
    rc: Number(ui.rc?.value || 10),
    gamma: Number(ui.gamma?.value || 2),
    temp: Number(ui.temp?.value || 300),
    physicsLevel: physLvl.level,
    binding: { on: ui.bindPot?.checked ?? true, holo: ui.holoSprings?.checked ?? true, charges: physLvl.charges, hbMode: physLvl.hbMode },
    weak: physLvl.weak === "on" ? "on" : "off",
  };
}

/**
 * Hot-rebuild force-field parameter literal (onParamChange path). Deliberately
 * a twin of ffParams with its own tier binding: the two paths feed different
 * force-field classes and a drift between them is silent (rev2/issue1 pins
 * both spellings — tests/test_rev2_issue1_physics_level.js).
 * @param {{level:string,charges:boolean,hbMode:string,weak:string}} physLvlHot active tier spec
 */
export function ffParamsHot(physLvlHot) {
  return {
    rc: Number(ui.rc?.value || 10),
    gamma: Number(ui.gamma?.value || 2),
    temp: Number(ui.temp?.value || 300),
    physicsLevel: physLvlHot.level,
    binding: { on: ui.bindPot?.checked ?? true, holo: ui.holoSprings?.checked ?? true, charges: physLvlHot.charges, hbMode: physLvlHot.hbMode },
    weak: physLvlHot.weak === "on" ? "on" : "off",
  };
}

/* ------------------------------------------------------------------ */
/*  Startup sequence (original order preserved)                          */
/* ------------------------------------------------------------------ */

// Initialize UI modals & panels
initSettingsModal();
initNetworkPanel();
// Stage-5: restore persisted physics tier onto the selector (default L0).
try { restorePhysicsLevelSelect(); } catch (_) {}

// A01 — Deterministic build version in HUD and console
if (typeof console !== "undefined") console.log(`[simulation_coarse] version ${VERSION} build ${BUILD_DATE}`);
if (ui.hud && !ui.hud.textContent.includes(VERSION)) {
  // Append version badge; tick() will keep it as prefix
  ui.hud.dataset.version = VERSION;
  ui.hud.title = `simulation_coarse ${VERSION} (${BUILD_DATE})`;
}

// Dock strips + scrub, then the BindViz empty states, then the parameter
// surface, then the subsystem wiring, and the animation loop last.
initDock();
initBindviz();
initParamBinding({ physicsLevelSpec, ffParams, ffParamsHot });
initSystemBuild({ physicsLevelSpec, ffParams });
initStructureInput({ buildSystem, updateGuide });
initGuide();
initTransport();
initRecording();

startTick();

/* ------------------------------------------------------------------ */
/*  Historical re-exports (kept so no existing import breaks)           */
/* ------------------------------------------------------------------ */
export { useWorkerPoolIfNeeded } from "./controllers/accelerate.js";
export { invalidateHeavyBuild } from "./controllers/system-build.js";
export {
  liveTermsWanted,
  applyLiveTrackTerms,
  readLiveTerms,
  formatLiveTermsHUD,
  updateLiveTermsMirror,
  liveTermsBindLogView,
} from "./controllers/live-terms.js";