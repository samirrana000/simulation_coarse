/**
 * controllers/dg-chain-ui.js — let a user expand any ΔG the app prints into the
 * arithmetic that produced it, with every term tagged measured / assumed /
 * out-of-scope.
 *
 * THE PROBLEM
 * -----------
 * The HUD shows a bare `ΔG ≈ X ± SE`. docs/VALIDATION.md shows a different ΔG
 * for the same system with a different estimator and a different pair of legs.
 * Nothing in the app said which estimator a number came from or what it
 * assumed, so the two stories could not be compared by anyone. src/dg-chain.js
 * now holds that comparison as data; this module is the only place it touches a
 * DOM.
 *
 * WHERE IT LIVES, AND WHY THERE
 * -----------------------------
 * Inside the EXISTING top-level `PMF & Analysis` panel, as a `<details
 * class="subpanel">` next to the existing Thermodynamics / Mutagenesis /
 * Correlations / Kinetics / BindViz subpanels. Two contracts forbid the
 * alternative:
 *   - scripts/wikiskill_gate.js DOM contract: every `$("id")` in src/ui.js must
 *     exist in index.html (so the new `#dgChain` handle is registered in BOTH).
 *   - Digit1-7 hotkey contract (src/controllers/transport.js binds Digit1-7 to
 *     `#controls > .panel` in DOM order): a new top-level panel would shift every
 *     digit after it. A subpanel adds none.
 * Both are asserted by scripts/wikiskill_gate.js and by
 * tests/test_dg_drift.js's panel/id check.
 *
 * ZERO PHYSICS HERE
 * -----------------
 * No kernel, no force field, no integrator step, no re-derivation of a number.
 * It reads `state.funnel` (for the live funnel chain) and whatever
 * `thermoResult()` was last produced by the ΔH/ΔS button, and renders. If the
 * estimator is absent the panel says so — absence is a first-class state here,
 * never a 0.
 *
 * DOM-FREE CORE
 * -------------
 * All the arithmetic and the chain data structure live in src/dg-chain.js,
 * which has no DOM dependency, so `node tests/test_dg_drift.js` exercises the
 * same numbers the browser renders.
 */

import { ui, state } from "../ui.js";
import { funnelChain, thermoChain, reconcileChains, chainMarkdown, DG_ESTIMATORS, DG_BRIDGE } from "../dg-chain.js";
import { ignore } from "../errors.js";

/**
 * The last ΔH/ΔS result, published by src/analysis-panel.js as
 * `state.thermoResult` at the two places it calls `computeThermodynamics`.
 * Read here rather than passed in: analysis-panel.js is pinned at its
 * module-size ratchet, and `state` is the surface it already writes to.
 */
export function thermoResult() { return state?.thermoResult ?? null; }

/** Last funnel state key rendered, so the panel redraws only when the hills move. */
let _lastKey = null;

/** The estimator the chain panel is currently expanded for. */
let _which = "funnel_dg_hud";

/**
 * Record the ΔH/ΔS result so the chain panel can expand it. Preferred path when
 * a caller has the object; the live page path is `state.thermoResult`, which
 * src/analysis-panel.js sets itself.
 * @param {object|null} r computeThermodynamics output (null clears)
 */
export function setThermoResult(r) {
  if (state) state.thermoResult = r ?? null;
  renderChain(true);
}

/**
 * One line per estimator: what it is, whether a user can see it, and whether
 * its number is expanded here. This is the enumeration the user is owed.
 * @returns {string}
 */
function estimatorListText() {
  const L = ["Every ΔG-like number this build produces (src/dg-chain.js DG_ESTIMATORS):"];
  for (const e of DG_ESTIMATORS) {
    const shown = e.onDisplayPath ? "shown in UI" : "computed, never displayed";
    L.push(`  ${e.id.padEnd(22)} ${shown.padEnd(24)} ${e.label}`);
  }
  return L.join("\n");
}

/** Why the two stories do not meet — the honest answer, not a TODO. */
function bridgeText(hud, doc) {
  const rec = reconcileChains(hud, doc);
  const L = [];
  L.push(rec.verdict);
  L.push("");
  L.push("Why they cannot be reconciled with what this tool implements (src/dg-chain.js DG_BRIDGE):");
  for (const b of rec.bridge) {
    const moved = Number.isFinite(b.movedBy) ? ` [measured move: ${b.movedBy > 0 ? "+" : ""}${b.movedBy.toFixed(2)} ${b.movedByUnit ?? "kcal/mol"}]` : "";
    L.push(`  [${b.nature}] ${b.headline}${moved}`);
    if (b.missing) L.push(`      missing: ${b.missing}`);
  }
  return L.join("\n");
}

/**
 * Render the chain for the selected estimator into `#dgChain`.
 *
 * Deliberately not throttled and not tied to the HUD cadence: this panel is
 * collapsed by default and reads one already-computed estimator, so there is
 * nothing to rate-limit and nothing to recompute.
 * @param {boolean} [force] re-render even if the funnel has not changed
 * @returns {string|null} the rendered text (null when there is no DOM)
 */
export function renderChain(force = false) {
  const el = ui.dgChain ?? (typeof document !== "undefined" ? document.getElementById("dgChain") : null);
  if (!el) return null;
  const funnel = state?.funnel ?? null;
  const thermo = thermoResult();
  const key = `${funnel ? `${funnel._nHills}|${funnel.active}` : "none"}|${_which}|${Number.isFinite(thermo?.dG_estimate) ? thermo.dG_estimate : "none"}`;
  if (!force && key === _lastKey) return el.textContent ?? null;

  const hud = funnelChain(funnel);
  const doc = thermo ? thermoChain(thermo) : null;
  const body = _which === "thermo_dg"
    ? (doc ? chainMarkdown(doc, { title: "ΔG estimate — ΔH − T·ΔS (docs/VALIDATION.md rows)" })
      : "No ΔH/ΔS decomposition yet — run ΔH/ΔS Decomposition in this panel, then this chain fills in.\nNo value is shown rather than a 0: src/scope.js states the absence explicitly.")
    : chainMarkdown(hud, { title: "HUD ΔG — well-tempered funnel bias" });
  const which = _which === "thermo_dg" ? "ΔH − T·ΔS" : "HUD funnel bias";
  const alt = _which === "thermo_dg" ? "HUD funnel bias" : "ΔH − T·ΔS";

  el.textContent = [
    body,
    "",
    `Showing: ${which}. The other story (${alt}) is a different estimator on different legs —`,
    "they are NOT two measurements of one quantity, and this chain will not average them.",
    "",
    bridgeText(hud, doc),
    "",
    estimatorListText(),
  ].join("\n");
  _lastKey = key;
  return el.textContent;
}

/**
 * Toggle which story the panel expands. Both live in the same panel so neither
 * adds a top-level entry to the Digit1-7 order.
 * @param {"funnel_dg_hud"|"thermo_dg"} which
 */
export function showChain(which) {
  _which = which === "thermo_dg" ? "thermo_dg" : "funnel_dg_hud";
  renderChain(true);
}

/** Which story the panel is currently expanded for. */
export function currentChain() { return _which; }

/**
 * Wire the two buttons. Safe headless: both handles are optional and the
 * `ui.dgChain` fallback is a guarded lookup, not a bare `document.getElementById`.
 */
export function initChainUi() {
  const bind = (btn, which) => {
    try { btn?.addEventListener("click", () => showChain(which)); } catch (e) {
      ignore(e, "dgChainUi@bind", "the handle is absent in a headless import; the chain data module must still load");
    }
  };
  bind(ui.dgChainHud, "funnel_dg_hud");
  bind(ui.dgChainThermo, "thermo_dg");
  try { renderChain(true); } catch (e) {
    ignore(e, "dgChainUi@render", "no #dgChain in this document; the chain module is still importable and testable");
  }
}
