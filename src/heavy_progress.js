/**
 * heavy_progress.js — FP5: heavy-build progress UX helpers (headless-pure).
 *
 * Context: heavy mode costs ~15 ms/step here (~77–85 ms/step on the S7 pareto
 * machine) and the O(n²) topology build ~30–60 ms on 4W52 (1308 atoms), so a
 * synchronous heavy Build freezes the tab's paint. S7 decided ACCEPT CPU (no
 * GPU port), so FP5 keeps the physics on CPU and makes the build honest:
 * chunked topology rows (see buildTopologyChunked in src/heavy.js) with
 * progress captions, Build/Run disabled during work, and cooperative cancel
 * via a generation counter (same pattern as the thermo `_thermoGen` in
 * src/analysis-panel.js).
 *
 * This module is DOM-light on purpose: helpers take explicit element refs
 * ({ disabled } / { textContent } shape) so headless proofs can pass minimal refs —
 * src/main.js wires the real elements. Zero deps. No physics, no defaults.
 */
import { ignore } from "./errors.js";

/** Bond-loop rows per chunked-topology slice (FP5: 1308 rows → 11 slices). */
export const HEAVY_TOPO_CHUNK_ROWS = 128;

/** Cancel sentinel message thrown by buildTopologyChunked (see src/heavy.js). */
export const HEAVY_BUILD_CANCELLED = "heavy build cancelled";

/**
 * Toggle the build/run/cancel control set for a heavy build (FP5).
 * Null-safe + headless-safe (minimal `{ disabled }` refs work).
 * While building: Build/Run/Reset off, Cancel armed. When idle: inverse.
 * @param {{buildBtn?:{disabled:boolean}|null, playBtn?:{disabled:boolean}|null, resetBtn?:{disabled:boolean}|null, cancelBtn?:{disabled:boolean}|null}} els control refs
 * @param {boolean} building true while a chunked heavy build is in flight
 */
export function setHeavyButtons(els = {}, building = false) {
  const on = !!building;
  try { if (els.buildBtn) els.buildBtn.disabled = on; } catch (e) { ignore(e, "buildBtn@setHeavyButtons", "headless callers pass null or minimal element refs by design"); }
  try { if (els.playBtn) els.playBtn.disabled = on; } catch (e) { ignore(e, "playBtn@setHeavyButtons", "headless callers pass null or minimal element refs by design"); }
  try { if (els.resetBtn) els.resetBtn.disabled = on; } catch (e) { ignore(e, "resetBtn@setHeavyButtons", "headless callers pass null or minimal element refs by design"); }
  try { if (els.cancelBtn) els.cancelBtn.disabled = !on; } catch (e) { ignore(e, "cancelBtn@setHeavyButtons", "headless callers pass null or minimal element refs by design"); }
}

/**
 * Write a heavy-build progress line to the build summary surface (FP5).
 * Reuses the existing `#selSummary` element (no new caption ids).
 * Null-safe + headless-safe.
 * @param {{textContent:string}|null} el summary element (or a minimal ref)
 * @param {string} txt progress line
 */
export function setHeavyCaption(el, txt) {
  try { if (el) el.textContent = txt; } catch (e) { ignore(e, "caption@setHeavyCaption", "caption is cosmetic; the build itself already completed"); }
}

/**
 * Format a topology progress line (FP5, pure).
 * @param {number} done rows finished
 * @param {number} total total rows
 * @returns {string} e.g. "Building heavy… topology 512/1308 rows (UI stays responsive)."
 */
export function heavyTopoCaption(done, total) {
  return `Building heavy… topology ${done}/${total} rows (UI stays responsive).`;
}
