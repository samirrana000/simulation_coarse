/**
 * controllers/param-binding.js — physics-parameter hot reload + slider readouts.
 *
 * PUBLIC API
 *   initParamBinding(seams) Wire the parameter surface: the slider/numeric
 *                     readout pairs (#rc/#gamma/#temp/#fric/#mass/#motionGain
 *                     via ui.initParamReadouts), the #motionGain viewer gain,
 *                     the #bindPot / #holoSprings toggles and the
 *                     #physicsLevel selector. `seams` is the composition
 *                     root's cross-cutting tier/params triple (see
 *                     src/main.js): {physicsLevelSpec, ffParams, ffParamsHot}.
 *                     Must be called exactly once, before any Build.
 *   onParamChange(rebuildContacts = true)
 *                     Push the current slider values into the live integrator,
 *                     and (by default) hot-rebuild the force field + funnel
 *                     while preserving positions/velocities/time.
 *
 * NEEDS (imports): ui + state + viewer + initParamReadouts +
 * updateSelSummary from ../ui.js; settingsState + workerPool +
 * persistPhysicsLevel from ../settings-panel.js; ForceField from
 * ../forcefield.js; HeavyForceField from ../heavy.js; Funnel from
 * ../funnel.js; applyMLToFF from ../ml-tier.js; bindLogWanted from
 * ./physics-tier.js; applyLiveTrackTerms from ./live-terms.js.
 *
 * WHY `seams` IS INJECTED rather than imported
 *   `physicsLevelSpec()` and the two force-field parameter literals are the
 *   app's two build-path contracts and live in the composition root
 *   (src/main.js). Importing them from there would close an import cycle
 *   main.js -> params.js -> main.js, so the root passes them in instead. This
 *   is the only injection seam in the app; everything else is a plain DAG.
 *
 * ZERO PHYSICS ADDED: the hot rebuild re-constructs the SAME force field from
 * the SAME `par` shape the cold build path uses — the literals are not edited
 * here, only consumed.
 */

import { ui, state, viewer, initParamReadouts, updateSelSummary } from "../ui.js";
import { workerPool, persistPhysicsLevel } from "../settings-panel.js";
import { ForceField } from "../forcefield.js";
import { HeavyForceField } from "../heavy.js";
import { Funnel } from "../funnel.js";
import { applyMLToFF } from "../ml-tier.js";
import { bindLogWanted } from "./physics-tier.js";
import { applyLiveTrackTerms } from "./live-terms.js";

/** Cross-cutting seams injected by initParamBinding (see module header). */
let _seams = null;

/**
 * Push slider values into the live integrator; hot-rebuild the force field,
 * funnel, ML tier and viewer system by default.
 * @param {boolean} [rebuildContacts=true] false for the pure readout path
 */
export function onParamChange(rebuildContacts = true) {
  // Unreachable in the app: initParamBinding is the first wiring call in
  // src/main.js, long before any control can fire. Guarded anyway so a
  // headless harness that calls this before init fails loudly instead of
  // dereferencing an absent seam.
  if (!_seams) return;
  if (!state.integ || !state.ff) return;
  const integ = state.integ;
  if (state.ff.funnel) state.ff.funnel.setTemperature(Number(ui.temp?.value || 300));
  integ.setTemperature(Number(ui.temp?.value || 300));
  integ.setFriction(Number(ui.fric?.value || 8));
  if (!state.heavyMode) {
    integ.rebuildMass(Number(ui.mass?.value || 110));
  }

  if (rebuildContacts) {
    // Loop-2 S7: hot-reload path carries the same tier flags as buildSystem.
    const physLvlHot = _seams.physicsLevelSpec();
    const par = _seams.ffParamsHot(physLvlHot);
    const keepPos = Float64Array.from(integ.pos);
    const keepVel = Float64Array.from(integ.vel);
    const keepTime = integ.time;

    if (state.heavyMode) {
      state.ff = new HeavyForceField({ atoms: state.sel.atoms }, par, []);
      workerPool.initSystem(state.ff);
    } else {
      state.ff = new ForceField(state.sel, par, state.ligands);
    }

    if (keepPos && keepPos.length === state.ff.n * 3) {
      state.integ.pos.set(keepPos);
      state.integ.vel.set(keepVel);
      state.integ.time = keepTime;
    }
    state._prevPos = null;
    state.nanWarning = false;
    // Loop-2 S4+S7 + Rev1/Issue4: reapply the live trackTerms policy on the
    // rebuilt FF (accumulators stay on while a ligand is present for the HUD
    // mirror; BindLog event capture still requires bindLogWanted()).
    applyLiveTrackTerms(state.ff);
    state.integ.setTemperature(Number(ui.temp?.value || 300));
    state.integ.setFriction(Number(ui.fric?.value || 8));

    if (state.ff.nLigAtoms > 0) {
      state.funnel = new Funnel({
        nProt: state.ff.nProt,
        n: state.ff.n,
        ref: state.ff.ref,
        ligStart: state.ff.ligandStart ?? state.ff.nProt,
      });
      state.ff.setFunnel(state.funnel);
      state.ff.funnelOn = ui.funnelToggle?.checked ?? false;
      // Loop-2 S4+S7: hill hook follows the rebuild (same guard as buildSystem)
      if (bindLogWanted()) {
        state.funnel.onHill = (cv, h) => state.bindLog && state.bindLog.pushHill(state.integ.time, cv, h);
      } else {
        state.funnel.onHill = null;
      }
    } else {
      state.funnel = null;
    }
    applyMLToFF();
    if (viewer) viewer.setSystem(state.sel, state.ff);
    else console.warn("[viewer] not ready");
    updateSelSummary();
  }
}

/** Wire every parameter control. Call once, before any Build. */
export function initParamBinding(seams) {
  _seams = seams;

  // physics-slider live readouts
  initParamReadouts(() => onParamChange());

  if (ui.motionGain && ui.v_motionGain) {
    ui.motionGain.addEventListener("input", () => {
      ui.v_motionGain.textContent = ui.motionGain.value;
      if (viewer) viewer.setMotionGain(Number(ui.motionGain.value));
      else console.warn("[viewer] not ready");
    });
  }

  if (ui.bindPot) ui.bindPot.addEventListener("change", () => onParamChange(true));
  if (ui.holoSprings) ui.holoSprings.addEventListener("change", () => onParamChange(true));
  // Loop-2 S7: physics-level selector — hot-rebuilds the FF on the new tier
  // (same path as the binding-potential toggles; default L0 = baseline).
  if (ui.physicsLevel) ui.physicsLevel.addEventListener("change", () => {
    const spec = _seams.physicsLevelSpec(); // persist to settingsState.physicsLevel
    try { persistPhysicsLevel(ui.physicsLevel.value); } catch (_) {}
    void spec;
    onParamChange(true);
  });
}