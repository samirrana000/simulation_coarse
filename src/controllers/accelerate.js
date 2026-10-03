/**
 * controllers/accelerate.js — force-evaluation backends the tick loop may use
 * instead of a plain synchronous `integ.advance()`.
 *
 * PUBLIC API
 *   respaWanted()               -> boolean. r-RESPA opt-in (default OFF) read
 *                                 from the Dynamics checkbox, else the setting.
 *   ensureRespa(ff)             -> the cached RESPA context (throws when the
 *                                 force field cannot be split; the tick then
 *                                 falls back to single-step BAOAB).
 *   warnRespaFallbackOnce(err)  -> latch the one-per-force-field "frame
 *                                 fallback to BAOAB" console warning.
 *   useWorkerPoolIfNeeded(pos)  -> Promise<{forces,...}> | null (see below).
 *
 * NEEDS (imports): ui + state from ../ui.js; settingsState + workerPool from
 * ../settings-panel.js; RESPAStepper + splitForceField from
 * ../physics/integrators/respa.js.
 *
 * ZERO PHYSICS ADDED: r-RESPA is a multiple-time-stepping scheme over the SAME
 * kernels, opt-in, with parity fallback; the worker pool path is feature-flagged
 * and off by default. Neither is on by default.
 */

import { ui, state } from "../ui.js";
import { ignore } from "../errors.js";
import { settingsState, workerPool } from "../settings-panel.js";
import { RESPAStepper, splitForceField } from "../physics/integrators/respa.js";

/**
 * G62 — Wire worker pool (minimal, feature-flagged) — speedup≥1.5× is aspirational, not yet benchmarked
 * Feature flag: settingsState.backend === "workers" && state.ff.n > 800 routes to workerPool.computeParallel
 * When enabled, offload non-bonded force evaluation to workerPool.computeParallel().
 * Currently guarded by n > 800 and backend === "workers" to avoid overhead on small systems.
 * NOTE: speedup≥1.5× is aspirational, not yet benchmarked — validate via per-step
 * timing before enabling by default. Fallback is ff.compute().
 * NOT WIRED: useWorkerPoolIfNeeded() below is exported (src/main.js re-exports
 * it) but has NO caller — src/controllers/tick.js imports only respaWanted,
 * ensureRespa and warnRespaFallbackOnce from this module, and the call site at
 * the bottom of this file is still commented out. That is honest scope, not a
 * defect, and it is why the flag defaults are what they are.
 *
 * MEASURED 2026-10-03 (removing-redundancy pass). This block previously claimed
 * "G62 wiring validated via bench/worker_speedup.js (estimated 1.5× on 4 cores)"
 * and pointed at that script three times as the thing to run. It validates
 * NOTHING: bench/worker_speedup.js opens with
 * `console.log("estimated 1.5× on 4 cores")` — an unconditional literal — and
 * then prints an Amdahl result computed from a hardcoded parallelFrac = 0.8. It
 * never constructs a worker, never times a compute, and never imports anything
 * from src/. So "validated" was false and is now deleted rather than softened.
 * The speedup figure remains unmeasured, which is the honest state:
 * worker_speedup is not a key in bench/budget.json, so
 * tests/test_budget_coverage.js does not (and must not) treat it as measured.
 */
export function useWorkerPoolIfNeeded(pos) {
  if (settingsState.backend === "workers" && state.ff && state.ff.n > 800) {
    // workerPool.computeParallel returns Promise<{forces, lj, elec}> — caller must await and merge forces
    // G62 guard: only for large systems where parallel overhead is amortized; not yet validated for speedup
    return workerPool.computeParallel(pos, state.ff.n);
  }
  return null; // fallback to ff.compute(pos) (synchronous)
}
// Example branch (not yet active in tick — see comment above for integration point):
//   const parallel = await useWorkerPoolIfNeeded(state.integ.pos);
//   if (parallel) { state.ff.forces.set(parallel.forces); } else { state.ff.compute(state.integ.pos); }

// Phase 3 — r-RESPA multiple-time-stepping, opt-in (default OFF via the
// Dynamics panel checkbox). Cache is keyed on the live force field identity
// so rebuilds (buildSystem/onParamChange) transparently re-split. Any setup
// failure (e.g. useAmber14 heavy) latches `disabled` for that ff and the tick
// falls back to single-step BAOAB — parity fallback, never a crash.
let respaCache = { ff: null, stepper: null, split: null, disabled: false, reason: "", outerFs: 4 };

/** True when r-RESPA is requested (checkbox first, then the persisted setting). */
export function respaWanted() {
  try {
    if (ui.respaToggle && typeof ui.respaToggle.checked === "boolean") return ui.respaToggle.checked;
  } catch (e) { ignore(e, "respaToggle@respaWanted", "ui ref null in a headless import; the persisted setting below is the documented fallback"); }
  return !!settingsState.respaOn;
}

/** Build (or reuse) the RESPA stepper + force-field split for `ff`. */
export function ensureRespa(ff) {
  const outerFs = Number(ui.respaOuter?.value) || Number(settingsState.respaOuterFs) || 4;
  if (respaCache.ff === ff && respaCache.stepper && !respaCache.disabled) {
    if (respaCache.outerFs !== outerFs) {
      respaCache.stepper.setSteps(0.001, outerFs / 1000);
      respaCache.outerFs = outerFs;
    }
    return respaCache;
  }
  if (respaCache.ff === ff && respaCache.disabled) throw new Error(respaCache.reason || "r-RESPA disabled");
  respaCache = { ff, stepper: null, split: null, disabled: false, reason: "", outerFs };
  try {
    respaCache.stepper = RESPAStepper.fromForceField(ff, {
      temperature: Number(ui.temp?.value || 300),
      friction: Number(ui.fric?.value || 8),
      dtInner: 0.001,
      dtOuter: outerFs / 1000,
    });
    respaCache.split = splitForceField(ff);
    if (state.integ) respaCache.stepper.time = state.integ.time;
  } catch (e) {
    respaCache.disabled = true;
    respaCache.reason = e?.message ?? String(e);
    throw e;
  }
  return respaCache;
}

/**
 * Warn once per force field when the frame falls back to BAOAB. The latch
 * lives on the cache object, so it resets automatically when the force field is
 * rebuilt (a new cache entry has no `_warned`).
 * @param {unknown} e the setup/step error
 */
export function warnRespaFallbackOnce(e) {
  if (!respaCache._warned) {
    console.warn(`[r-RESPA] frame fallback to BAOAB (${e?.message ?? e})`);
    respaCache._warned = true;
  }
}