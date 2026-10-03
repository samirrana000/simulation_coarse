/**
 * dg-chain.js — ONE auditable chain from input to every ΔG-like number the app
 * prints. Each term tagged measured | assumed | out-of-scope, with its value
 * and its uncertainty where one exists and an explicit sentence where one does
 * not.
 *
 * THE USER PROBLEM THIS SOLVES
 * ---------------------------
 * A user runs the app and sees a ΔG in the HUD. docs/VALIDATION.md reports a
 * completely different ΔG for the same system (4W52 + benzene). Neither is wrong
 * exactly — they come from different estimators — but nothing in either said
 * which estimator produced which number, what it assumed, or how far it sits from
 * the other. A scientist comparing the two had no way to know what they were
 * looking at.
 *
 * WHAT LIVES WHERE
 * ----------------
 *   src/dg-estimators.js  WHICH numbers exist: DG_ESTIMATORS (the enumeration,
 *                         with each estimator's formula, code path, inputs,
 *                         measured/assumed status, reported uncertainty and the
 *                         source sites that print it) and DG_BRIDGE (why the two
 *                         headline stories cannot be reconciled with what this
 *                         tool implements, each entry naming the missing term).
 *   this file              WHAT ONE number is made of: `funnelChain()` expands a
 *                         live Funnel, `thermoChain()` expands a
 *                         `computeThermodynamics()` result, and the renderers
 *                         turn either into text or flat rows.
 *
 * THE CHAIN, NOT PROSE
 * -------------------
 * A chain is an ordered list of `{term, value, unit, nature, note, uncertainty}`
 * rows. `nature` is one of exactly three tags:
 *
 *   measured     — produced by sampling this system's own trajectory
 *   assumed      — a scale, a reference state, or a leg-definition choice that
 *                  no sampling in this tool can put an error bar on
 *   out-of-scope — physics this tool does not have, named so a reader is not
 *                  left inferring it is missing by accident
 *
 * THERE IS NEVER A ZERO ERROR BAR
 * ------------------------------
 * `uncertaintyOf()` applies the project's rule (src/results-record.js header,
 * rule 1) verbatim: a zero, a NaN, an Infinity or an absent error bar becomes
 * `uncertainty: null` plus a reason. It reuses that rule rather than inventing a
 * second convention, and the chain says "NO ERROR BAR: …" out loud instead of
 * printing 0. The 0-sentinel it protects against is real: `computeThermodynamics`
 * leaves `dH_se` at 0 when the 20-block bootstrap has fewer than two blocks,
 * `thermo_sasa` passes `se ?? 0`, and `convergenceSE()` returns null below 50
 * hills.
 *
 * DOM-FREE AND NODE-IMPORTABLE
 * ----------------------------
 * No `document`, no `window`, no `fetch`, no `node:` imports — importing this
 * from plain Node loads nothing that needs a browser. Rendering into the page is
 * src/controllers/dg-chain-ui.js's job; that module is the only thing here that
 * would ever care about a DOM.
 */

// The unit contract, imported rather than re-typed: tests/test_constant_ledger.js
// fails any module that re-derives KB_KCAL or STANDARD_VOLUME under any name, and
// a second copy of either is exactly how a reported number stops agreeing with
// the estimator that produced it. src/units.js is a zero-import leaf, so this edge
// cannot create a cycle.
import { KB_KCAL, STANDARD_VOLUME } from "./units.js";
import { DG_ESTIMATORS, DG_BRIDGE } from "./dg-estimators.js";

export { DG_ESTIMATORS, DG_BRIDGE };

/**
 * The three tags a chain term may carry. A term with any other tag is a bug, and
 * tests/test_dg_drift.js asserts the closed set.
 * @type {readonly ["measured", "assumed", "out-of-scope"]}
 */
export const TERM_NATURE = Object.freeze(["measured", "assumed", "out-of-scope"]);

/** SASA_GAMMA, the solvent scale in src/analysis/thermodynamics.js (kcal/mol/Å²). */
export const CHAIN_SASA_GAMMA = 0.012;

/**
 * One term of a chain.
 * @typedef {{term:string, value:number|null, unit:string, nature:string,
 *            note:string, uncertainty:{value:number, kind:string, method:string}|null,
 *            noErrorBar:string|null}} ChainTerm
 */

/**
 * The project's rule 1 from src/results-record.js, reused rather than
 * re-derived: a zero, a NaN, an Infinity or an absent error bar is
 * `uncertainty: null` plus a reason. NEVER a 0, because a 0 error bar claims a
 * precision nobody measured.
 *
 * @param {number|null|undefined} value
 * @param {string} kind e.g. "bootstrap" | "hill-count" | "block-jackknife"
 * @param {string} method how it was computed
 * @param {string} whyAbsent sentence for the case where there is no error bar
 * @returns {{uncertainty:{value:number, kind:string, method:string}|null, noErrorBar:string|null}}
 */
export function uncertaintyOf(value, kind, method, whyAbsent) {
  const v = Number(value);
  if (!Number.isFinite(v) || v <= 0) {
    return { uncertainty: null, noErrorBar: whyAbsent };
  }
  return { uncertainty: { value: v, kind, method }, noErrorBar: null };
}

/**
 * Build one chain term, normalising the nature tag and the error bar.
 * @param {{term:string, value:number|null, unit?:string, nature:string, note?:string,
 *          uncertainty?:number|null, uncertaintyKind?:string, uncertaintyMethod?:string,
 *          noErrorBar?:string}} spec
 * @returns {ChainTerm}
 */
function term(spec) {
  if (!TERM_NATURE.includes(spec.nature)) {
    throw new Error(`dg-chain: term "${spec.term}" has nature "${spec.nature}"; the closed set is ${TERM_NATURE.join(" | ")}`);
  }
  const v = Number(spec.value);
  const hasValue = spec.value !== null && spec.value !== undefined && Number.isFinite(v);
  const u = uncertaintyOf(spec.uncertainty, spec.uncertaintyKind ?? "unspecified",
    spec.uncertaintyMethod ?? "", spec.noErrorBar ?? "");
  return {
    term: String(spec.term),
    value: hasValue ? v : null,
    unit: String(spec.unit ?? "kcal/mol"),
    nature: spec.nature,
    note: String(spec.note ?? ""),
    uncertainty: u.uncertainty,
    noErrorBar: u.noErrorBar,
  };
}

/**
 * Build the funnel HUD chain from a live Funnel.
 *
 * Every value is read off the estimator at call time, so the chain's arithmetic
 * identity (term0 + term1 + … === the reported ΔG) is a property of the code,
 * not of a hand-typed table.
 *
 * @param {object|null} funnel a src/funnel.js Funnel (or null when absent)
 * @returns {{estimatorId:string, value:number|null, terms:ChainTerm[], note:string}}
 */
export function funnelChain(funnel) {
  const absent = {
    estimatorId: "funnel_dg_hud",
    value: null,
    terms: [],
    note: "no funnel in this build — there is no ΔG to expand (no ligand, or no ligand atoms in the selection)",
  };
  if (!funnel || funnel.active === false || !Number.isFinite(funnel._nHills) || funnel._nHills === 0) return absent;

  const g = funnel.getPMF();
  const rStep = funnel._rStep;
  const interp = (rWant) => {
    const x = rWant / rStep - 0.5;
    if (x < 0) return g.pmf[0];
    const k0 = Math.floor(x);
    if (k0 >= funnel.bins - 1) return g.pmf[funnel.bins - 1];
    const f = x - k0;
    return g.pmf[k0] + f * (g.pmf[k0 + 1] - g.pmf[k0]);
  };
  const rBound = funnel.cv0;
  const rFar = Math.min(rBound + 8.0, funnel.rMax - rStep);
  const deltaPMF = interp(rFar) - interp(rBound);
  const kBT = funnel.T * KB_KCAL;
  const vRest = (4.0 / 3.0) * Math.PI * Math.pow(funnel.rFlat, 3);
  const dGvol = -kBT * Math.log(Math.max(1e-6, vRest / STANDARD_VOLUME));

  const se = typeof funnel.convergenceSE === "function" ? funnel.convergenceSE() : null;
  const seWhy = funnel._nHills < 50
    ? `convergenceSE() returns null below 50 hills and this run has ${funnel._nHills} — the HUD is showing a collecting value, not a converged one`
    : "convergenceSE() is kT/√nHills, a hill-count decorrelation heuristic; it is NOT an uncertainty on ΔG and it under-reports the observed drift (see DG_BRIDGE hill-convergence)";

  const terms = [
    term({
      term: `reconstructed PMF difference, pmf(rFar=${rFar.toFixed(3)} Å) − pmf(rBound=${rBound.toFixed(3)} Å)`,
      value: deltaPMF,
      nature: "measured",
      note: `well-tempered bias γ=${funnel.biasFactor}, ${funnel._nHills} hills, σ=${funnel.sigma} Å, grid Δr=${rStep.toFixed(3)} Å; a TWO-POINT read (rFar = cv0 + 8 Å), not an integral`,
    }),
    term({
      term: `Boresch standard-state correction, −dG_vol = +kT·ln(V_rest/V°), V_rest=(4/3)π·${funnel.rFlat}³=${vRest.toFixed(1)} Å³, V°=${STANDARD_VOLUME} Å³`,
      value: -dGvol,
      nature: "assumed",
      note: "the restraint volume is a PARAMETER choice (rFlat), not a measurement; changing rFlat moves this term and nothing else",
    }),
    term({
      term: "= ΔG° as displayed (two-point estimate)",
      value: funnel.estimateDG(),
      nature: "measured",
      note: "identical to estimateDG() — this row is the reported number, the two rows above are its only arithmetic",
      uncertainty: se,
      uncertaintyKind: "hill-count heuristic",
      uncertaintyMethod: "src/funnel.js convergenceSE() = kT/√nHills; a decorrelation-time heuristic, NOT an uncertainty on ΔG",
      noErrorBar: seWhy,
    }),
    term({
      term: "Jacobian 2kT·ln r (rigorous ΔG° requires it; the two-point path omits it)",
      value: null,
      nature: "out-of-scope",
      note: `integrateDGJacobian() implements the r²-weighted quadrature and is NOT on the display path (docs/FUNNEL.md §7); on this grid it reads ${safeJson(funnel.integrateDGJacobian())}. Applying it is a physics change.`,
    }),
    term({
      term: "Tiwary–Parrinello c(t) offset",
      value: g.c_t,
      nature: "assumed",
      note: "computed by getPMF() and REPORTED in the PMF CSV provenance header, but never added to estimateDG(); at this hill count it is a non-zero shift the displayed number does not carry",
    }),
    term({
      term: "ΔH / ΔS decomposition (LJ · Coul · HB · EEF1-lite desolv · Schlitter · torsion · ΔSASA)",
      value: null,
      nature: "out-of-scope",
      note: "the funnel estimator has no decomposition at all; see the thermo_dg estimator for the one that does",
    }),
    term({
      term: "replica / block uncertainty on the reconstructed PMF",
      value: null,
      nature: "out-of-scope",
      note: "one trajectory, one seed; the code computes no block average and no between-replica spread",
    }),
  ];
  return {
    estimatorId: "funnel_dg_hud",
    value: funnel.estimateDG(),
    terms,
    note: `A vertical free-energy difference along ONE radial CV, read two points at ${funnel._nHills} hills. Not ΔH − T·ΔS, not an alchemical ΔG°, not a K_D.`,
  };
}

/**
 * Build the thermo ΔG chain from a `computeThermodynamics()` result.
 * @param {object|null} r computeThermodynamics output
 * @param {object} [opts] `{dsasaOverride, sasaSe}` for the real-burial leg
 * @returns {{estimatorId:string, value:number|null, terms:ChainTerm[], note:string}}
 */
export function thermoChain(r, opts = {}) {
  if (!r || !Number.isFinite(r.dG_estimate)) {
    return { estimatorId: "thermo_dg", value: null, terms: [], note: "no ΔH/ΔS decomposition computed" };
  }
  const T = r.meta?.T ?? 300;
  const negTdsPocket = -T * r.dS.pocket;
  const negTdsLig = -T * r.dS.ligand;
  const negTdsSolv = -T * r.dS.solvent;
  const sasaSe = opts.sasaSe ?? r.meta?.sasa?.se ?? null;
  const sasaScale = CHAIN_SASA_GAMMA; // SASA_GAMMA, src/analysis/thermodynamics.js

  const terms = [
    term({
      term: "ΔH total (per-frame mean binding enthalpy, holo leg; apo binding ≡ 0 by construction)",
      value: r.dH.total,
      nature: "measured",
      note: `components LJ ${r.dH.lj.toFixed(2)} / Coul ${r.dH.coul.toFixed(2)} / HB ${r.dH.hb.toFixed(2)} / desolv ${r.dH.desolv.toFixed(2)}; the desolv component is the ASSUMED EEF1-lite per-atom table (eef1_lig_dg)`,
      uncertainty: r.dH_se,
      uncertaintyKind: "block bootstrap",
      uncertaintyMethod: "20-block bootstrap on the per-frame 7-term total, src/analysis/thermodynamics.js",
      noErrorBar: r.dH_se > 0 ? "" : "computeThermodynamics leaves dH_se at 0 when the 20-block bootstrap has fewer than two blocks — a 0 sentinel meaning NOT COMPUTED, not zero error",
    }),
    term({
      term: "−T·ΔS_pocket (Schlitter quasi-harmonic, holo − apo)",
      value: negTdsPocket,
      nature: "measured",
      note: `frames/DOF ${(r.meta?.framesPerDof ?? 0).toFixed(1)} (need ≥10); mass model ${r.meta?.massModel ?? "?"}; the sign of this term is MODEL-DEPENDENT at Cα resolution`,
      noErrorBar: "no standard error is computed on the −TΔS leg; docs/VALIDATION.md measures a replica SD of 7.1 kcal/mol on 3 reps, which is the number that actually matters here",
    }),
    term({
      term: "−T·ΔS_ligand (torsion Shannon, bound)",
      value: negTdsLig,
      nature: "measured",
      note: `${r.meta?.ligNote ?? "?"}; a rigid ligand gives exactly 0 and this is a real value, not a missing one`,
    }),
    term({
      term: "−T·ΔS_solv (solvent, from measured ΔSASA)",
      value: negTdsSolv,
      nature: "measured",
      note: `ΔSASA ${r.meta?.dsasa} Å² × ${sasaScale} kcal/mol/Å², ±50% scale band (docs/VALIDATION.md §17 physics-sign review open)`,
      uncertainty: sasaSe === null || sasaSe === 0 ? null : sasaSe * sasaScale / T,
      uncertaintyKind: "SEM over evaluated frames",
      uncertaintyMethod: "√(SE_lig² + SE_prot²) × SASA_GAMMA / T, src/analysis/thermo_sasa.js",
      noErrorBar: "the legacy contact-count proxy path passes no SE at all (p.sasa.se ?? 0 → the 0 sentinel); with it, this row has NO error bar",
    }),
    term({
      term: "= ΔG estimate (ΔH − T·ΔS)",
      value: r.dG_estimate,
      nature: "assumed",
      note: "self-consistent as arithmetic, but the apo leg is an INTERNAL RELAXATION of the holo coordinates, not an independent apo simulation, so this is not a basin-to-basin free-energy difference",
      noErrorBar: "VALIDATION measures replica SD ≈ 7.1 kcal/mol on the −TΔS leg, which dominates any bootstrap SE on ΔH — the chain reports the dominant uncertainty rather than the small one",
    }),
  ];
  return {
    estimatorId: "thermo_dg",
    value: r.dG_estimate,
    terms,
    note: `A ΔH − T·ΔS difference between a holo leg and an internal apo relaxation of the same coordinates, at T=${T} K.`,
  };
}

/**
 * The reconciliation: the arithmetic that connects the two chains' numbers and
 * names, per term, what the gap is made of.
 *
 * `funnelChain` terms and `thermoChain` terms are NOT summable against each
 * other — they measure different things — so this function does not pretend
 * they are. It reports the size of the gap, the largest single assumption known
 * to move either side, and the measured move of each.
 *
 * @param {{value:number|null, terms:ChainTerm[]}|null} hud
 * @param {{value:number|null, terms:ChainTerm[]}|null} doc
 * @returns {{gap:number|null, bridge:object[], canReconcile:boolean, verdict:string}}
 */
export function reconcileChains(hud, doc) {
  const gap = (hud && doc && Number.isFinite(hud.value) && Number.isFinite(doc.value))
    ? hud.value - doc.value : null;
  return {
    gap,
    bridge: DG_BRIDGE.map((b) => ({ ...b })),
    canReconcile: false,
    verdict: gap === null
      ? "only one of the two stories is available in this build; there is nothing to reconcile yet"
      : `the displayed ΔG and the VALIDATION.md ΔG differ by ${gap >= 0 ? "+" : ""}${gap.toFixed(2)} kcal/mol, and DG_BRIDGE lists ${DG_BRIDGE.length} named reasons, every one of which is physics this tool does not implement. Reconciling them is not a documentation fix.`,
  };
}

/**
 * JSON without throwing on a value that is not serialisable.
 * @param {*} v
 * @returns {string}
 */
function safeJson(v) {
  try { return JSON.stringify(v); } catch (_) { return "(unserialisable)"; }
}

/**
 * THE EXHAUSTIVENESS RULE, as code: "a ΔG that appears in the UI but is not in
 * this table is a bug."
 *
 * `sites` in DG_ESTIMATORS are (file, literal source substring) pairs. This
 * function re-reads those files and reports two directions of drift:
 *   - a site whose literal is GONE from its file  (the table describes code
 *     that no longer exists)
 *   - a claimed file with no ΔG token at all      (a dead row)
 * Together with the site-count check in tests/test_dg_drift.js — which asserts
 * every estimator id is claimed by at least one LIVE site and that removing a
 * row makes this function's output disagree — that is the mechanical form of
 * the rule.
 *
 * @param {{read:(rel:string)=>string}} io file reader, so the test can feed a
 *        mutated source without touching disk
 * @param {ReadonlyArray<object>} [estimators] the table to check; defaults to
 *        DG_ESTIMATORS. Taking it as an argument is what lets
 *        tests/test_dg_drift.js delete a row and watch the checker notice.
 * @returns {{missing:string[], deadRows:string[], sitesChecked:number}}
 */
export function findUnenumeratedDGBinding(io, estimators = DG_ESTIMATORS) {
  const missing = [];
  const deadRows = [];
  let sitesChecked = 0;
  const cache = new Map();
  const read = (rel) => {
    if (!cache.has(rel)) cache.set(rel, String(io.read(rel) ?? ""));
    return cache.get(rel);
  };
  for (const est of estimators) {
    let live = 0;
    for (const site of est.sites) {
      sitesChecked++;
      const src = read(site.file);
      if (src.length === 0) { missing.push(`${est.id}: cannot read ${site.file}`); continue; }
      if (!src.includes(site.pattern)) {
        missing.push(`${est.id}: site literal not found in ${site.file} — ${JSON.stringify(site.pattern)}`);
      } else {
        live++;
      }
    }
    if (live === 0) deadRows.push(est.id);
  }
  return { missing, deadRows, sitesChecked };
}

/**
 * Render an estimator as one markdown table row (docs/VALIDATION.md table).
 * @param {object} e
 * @returns {string}
 */
export function estimatorTableRow(e) {
  const u = e.uncertainty
    ? `${e.uncertainty.kind} — ${e.uncertainty.method}`
    : `NONE — ${e.noErrorBar}`;
  return `| \`${e.id}\` | ${e.label} | \`${e.formula.replace(/\|/g, "\\|")}\` | \`${e.codePath}\` | ${e.nature} | ${String(u).replace(/\|/g, "\\|")} | ${e.onDisplayPath ? "yes" : "no (computed, never displayed)"} |`;
}

/** Header for the estimator table rendered by `estimatorTableRow`. */
export const ESTIMATOR_TABLE_HEADER =
  "| id | what the user sees | formula | code path | measured? | reported uncertainty | on the display path? |\n|---|---|---|---|---|---|---|";

/**
 * Render a chain as a plain-text block for a `<pre>` / a caption / a TXT export.
 * Each line carries its nature tag, so a number is never shown without saying
 * whether it was measured, assumed, or out of scope.
 * @param {{estimatorId:string, value:number|null, terms:ChainTerm[], note:string}} chain
 * @param {{title?:string}} [opts]
 * @returns {string}
 */
export function chainMarkdown(chain, opts = {}) {
  if (!chain || !chain.terms?.length) return `${opts.title ?? chain?.estimatorId ?? "ΔG"}: no value to expand (${chain?.note ?? "nothing computed"})`;
  const L = [];
  L.push(`── ΔG provenance chain · ${opts.title ?? chain.estimatorId} ──`);
  for (const t of chain.terms) {
    const val = t.value === null ? "—" : `${t.value >= 0 ? "+" : ""}${t.value.toFixed(4)}`;
    const bar = t.uncertainty ? ` ± ${t.uncertainty.value.toFixed(4)} [${t.uncertainty.kind}]` : "";
    const why = t.uncertainty ? "" : `  (NO ERROR BAR: ${t.noErrorBar})`;
    L.push(`  [${t.nature.padEnd(11)}] ${t.term}`);
    L.push(`  ${" ".repeat(13)}${val}${bar} ${t.unit}${why}`);
    if (t.note) L.push(`  ${" ".repeat(13)}— ${t.note}`);
  }
  if (chain.note) L.push(`  verdict: ${chain.note}`);
  return L.join("\n");
}

/**
 * The same chain as a flat, CSV-friendly row list — the shape a downstream
 * script can join on without parsing prose.
 * @param {{terms:ChainTerm[]}} chain
 * @returns {{term:string, value:number|null, unit:string, nature:string,
 *            uncertainty:number|null, uncertainty_kind:string|null,
 *            no_error_bar:string}[]}
 */
export function chainRows(chain) {
  return (chain?.terms ?? []).map((t) => ({
    term: t.term,
    value: t.value,
    unit: t.unit,
    nature: t.nature,
    uncertainty: t.uncertainty ? t.uncertainty.value : null,
    uncertainty_kind: t.uncertainty ? t.uncertainty.kind : null,
    no_error_bar: t.noErrorBar ?? "",
  }));
}
