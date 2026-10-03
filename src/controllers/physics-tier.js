/**
 * controllers/physics-tier.js — Loop-2 S7 physics fidelity tier + capture flag.
 *
 * PUBLIC API
 *   PHYSICS_LEVELS  {L0,L1,L2} -> {charges, hbMode, weak, bindLog} flag table.
 *   currentLevel()  -> "L0" | "L1" | "L2". Reads the #physicsLevel selector,
 *                     falls back to settingsState, coerces unknown/headless to
 *                     L0, and persists the resolved level back into
 *                     settingsState.physicsLevel (the existing in-memory
 *                     settings pattern).
 *   bindLogWanted() -> boolean. Effective BindLog event-capture flag: the
 *                     manual #bindlogOn checkbox OR the L2 full-rigor tier.
 *
 * NEEDS (imports): ui + state from ../ui.js, settingsState from
 * ../settings-panel.js. Nothing else — this module is the leaf of the tier
 * dependency so that every other controller can read the tier without
 * importing the composition root.
 *
 * WHY THE READ IS SPLIT FROM physicsLevelSpec()
 *   The spec-returning `physicsLevelSpec()` lives in src/main.js (see the
 *   "cross-cutting seams" section there for why), so the read+persist step it
 *   needs is factored out here. Callers that only want the level call
 *   currentLevel(); callers that want the flags call physicsLevelSpec(). The
 *   persisted side effect is identical either way.
 *
 * ZERO PHYSICS HERE: this table only chooses which FLAGS a force field is
 * constructed with. The kernels that consume them (ForceField,
 * HeavyForceField) are untouched.
 */

import { ui } from "../ui.js";
import { settingsState } from "../settings-panel.js";

/* ------------------------------------------------------------------ */
/*  Loop-2 S7: physics fidelity level (opt-in Dynamics-panel selector)  */
/* ------------------------------------------------------------------ */
/**
 * Fidelity-tier table for the physics-level selector (index.html
 * #physicsLevel, default L0). L0 = pre-Loop-2 baseline (bit-identical).
 * L1 adds CG salt-bridge charges (S1) + directional-HB virtual sites (S2).
 * L2 adds heavy weakint π/cation-π/halogen (S3) + BindLog per-term
 * accumulators (S4) on top. CG paths consume {charges, hbMode},
 * HeavyForceField mirrors physicsLevel/charges/hbMode queryably (Rev3/Issue1,
 * kernels unchanged) + consumes par.weak, both consume bindLog via
 * bindLogWanted(). All flags are enums/booleans (no units).
 * @type {Record<string, {charges:boolean, hbMode:string, weak:string, bindLog:boolean}>}
 */
export const PHYSICS_LEVELS = {
  L0: { charges: false, hbMode: "off", weak: "off", bindLog: false },
  L1: { charges: true, hbMode: "directional", weak: "off", bindLog: false },
  L2: { charges: true, hbMode: "directional", weak: "on", bindLog: true },
};

/**
 * Resolve the active physics level (guarded; headless/unknown → L0 baseline)
 * and persist it to settingsState.physicsLevel (in-memory, existing settings
 * pattern).
 * @returns {"L0"|"L1"|"L2"} active tier id
 */
export function currentLevel() {
  let lvl = "L0";
  try {
    const v = ui.physicsLevel?.value ?? settingsState.physicsLevel ?? "L0";
    if (PHYSICS_LEVELS[v]) lvl = v;
  } catch (_) { lvl = "L0"; }
  try { settingsState.physicsLevel = lvl; } catch (_) { /* headless */ }
  return lvl;
}

/**
 * Effective BindLog-capture flag: manual checkbox OR L2 full-rigor tier.
 * @returns {boolean} true when per-term accumulators + event capture stay on
 */
export function bindLogWanted() {
  try {
    return (ui.bindlogOn?.checked === true) || PHYSICS_LEVELS[currentLevel()].bindLog === true;
  } catch (_) { return false; }
}