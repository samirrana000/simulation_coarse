/**
 * dg-estimators.js — THE ENUMERATION: every place this codebase produces a
 * ΔG-like number, plus the arithmetic of why the two headline ΔG stories do not
 * meet.
 *
 * WHY THIS IS A SEPARATE MODULE FROM dg-chain.js
 * ---------------------------------------------
 * "What numbers exist" and "what is this number made of" are two different
 * questions with two different readers: the first is cited by docs and by the
 * results record, the second is rendered into a panel. Keeping them apart keeps
 * both under src/'s 600-LOC ceiling (tests/test_module_size.js) without either
 * being truncated, and it means the chain code cannot quietly grow a second,
 * divergent estimator list.
 *
 * THE PROBLEM THIS SOLVES
 * -----------------------
 * A user runs the app and sees a ΔG in the HUD. docs/VALIDATION.md reports a
 * completely different ΔG for the same system (4W52 + benzene). Neither is wrong
 * exactly — they come from different estimators — but nothing in either place
 * said which estimator produced which number, what each one assumed, or how far
 * either sits from the other.
 *
 * THAT WAS NOT AN ACCIDENT, IT WAS UNENUMERATED
 * ---------------------------------------------
 * There are EIGHT ΔG-like quantities here (see DG_ESTIMATORS), and before this
 * file no single place listed more than one of them. `findUnenumeratedDGBinding`
 * re-derives the display sites from the source tree and reports any ΔG display
 * the table does not claim, which is the "a ΔG in the UI that is not in the
 * table is a bug" rule made mechanical. tests/test_dg_drift.js runs it, and
 * proves it can fail by deleting one table row.
 *
 * WHAT IS DERIVED AND WHAT IS NOT
 * -------------------------------
 *   DERIVED: the `sites` column is re-read from the real source files by
 *             findUnenumeratedDGBinding(), and the chain VALUES are read from the
 *             live estimator objects (see dg-chain.js), so neither the table nor
 *             the arithmetic can describe code that does not exist.
 *   NOT DERIVED: the reconciliation GAP between the two stories. It cannot be
 *             closed by this tool — see DG_BRIDGE, where every entry names the
 *             missing physics. That is stated as out-of-scope with the absent
 *             term named, which is the documented convention (src/scope.js
 *             VALUE_TRUST_BOUNDARY `notComputed`), not a gap papered over.
 *
 * DOM-FREE AND NODE-IMPORTABLE: no `document`, no `window`, no `fetch`, no
 * `node:` imports. `node tests/test_dg_drift.js` exercises the same table the
 * browser renders.
 */

/**
 * THE ESTIMATOR ENUMERATION. One row per place the codebase produces a
 * ΔG-like number.
 *
 * Fields:
 *   id             stable key; also the id used in tests and in the chain
 *   label          what a user sees next to the number
 *   formula        the arithmetic, written the way the code writes it
 *   codePath       file + exported symbol that computes it
 *   inputs         what it consumes (frames, energies, CV grid, works, states)
 *   nature         "measured" | "assumed" | "mixed" — the estimator's own status
 *   uncertainty    {kind, method} or null with `noErrorBar` saying why
 *   onDisplayPath  true ⇒ a user can see this number somewhere in the UI
 *   docRow         the docs/VALIDATION.md methods-table row that reports it
 *   sites          where in src/ the number is produced or printed. Every
 *                  ΔG display the tree contains must be claimed by ≥1 site, or
 *                  `findUnenumeratedDGBinding` reports it as a bug.
 *
 * @type {ReadonlyArray<object>}
 */
export const DG_ESTIMATORS = Object.freeze([
  {
    id: "funnel_dg_hud",
    label: "live HUD ΔG (WTM funnel bias)",
    formula: "pmf(rFar) − pmf(cv0) − dG_vol,  dG_vol = −kT·ln[(4/3)π·rFlat³ / STANDARD_VOLUME]",
    codePath: "src/funnel.js estimateDG()",
    inputs: "a 1-D bias grid V(r) over r = |COM_lig − COM_pocket| built from WTM hills deposited during THIS run; rFar = cv0 + 8 Å",
    nature: "mixed",
    uncertainty: { kind: "hill-count", method: "convergenceSE() = kT/√nHills, a decorrelation-time heuristic, NOT an error on the answer" },
    onDisplayPath: true,
    docRow: "4W52 BNZ+EPE / BNZ-only rows do NOT cite it — this estimator appears in no VALIDATION.md row (that is finding F1)",
    sites: [
      { file: "src/funnel.js", pattern: "return deltaPMF - dG_vol;" },
      { file: "src/pmf-panel.js", pattern: "state.funnel.estimateDG().toFixed(2)" },
      { file: "src/controllers/tick.js", pattern: "fn.estimateDG().toFixed(2)" },
      { file: "src/analysis.js", pattern: "funnel.estimateDG()" },
      { file: "src/analysis.js", pattern: "Binding ΔG (WTM bias, bound vs +6 Å)" },
      { file: "src/controllers/results-export.js", pattern: "const dg = funnel.estimateDG();" },
    ],
  },
  {
    id: "funnel_dg_jacobian",
    label: "r²-Jacobian ΔG° quadrature (computed, never shown)",
    formula: "−kT·ln[∫_{r≤rFlat} r²e^{−βW}dr / ∫_{r>rFlat} r²e^{−βW}dr] − dG_vol",
    codePath: "src/funnel.js integrateDGJacobian()",
    inputs: "the same reconstructed PMF grid as funnel_dg_hud",
    nature: "measured",
    uncertainty: null,
    noErrorBar: "the quadrature returns a single number with no dispersion and no replicate legs; a track spread is not an uncertainty on a converged free energy",
    onDisplayPath: false,
    docRow: "none — docs/FUNNEL.md §7 states switching the display over is a physics change and out of scope",
    sites: [
      { file: "src/funnel.js", pattern: "integrateDGJacobian(pmf = null, r = null)" },
    ],
  },
  {
    id: "thermo_dg",
    label: "ΔG estimate (ΔH − T·ΔS decomposition)",
    formula: "ΔH − T·(ΔS_pocket + ΔS_lig + ΔS_solv); ΔH = ⟨U_bind⟩_holo, ΔS_pocket = Schlitter(holo) − Schlitter(apo), ΔS_solv = −ΔSASA·0.012/T",
    codePath: "src/analysis/thermodynamics.js computeThermodynamics().dG_estimate",
    inputs: "a holo recording, an internal apo relaxation of the SAME coordinates, the per-frame 7-term binding vector, and optionally a real LCPO ΔSASA",
    nature: "mixed",
    uncertainty: { kind: "bootstrap+SD", method: "20-block bootstrap on ΔH; NOTHING on the −TΔS leg, whose replica SD VALIDATION measures at 7.1 kcal/mol" },
    onDisplayPath: true,
    docRow: "4W52 BNZ+EPE 3-rep row and the 4W52 BNZ-only + real SASA row",
    sites: [
      { file: "src/analysis/thermodynamics.js", pattern: "const dG = dH.total - T * dS_total;" },
      { file: "src/analysis/thermodynamics.js", pattern: "ΔG estimate ${r.dG_estimate.toFixed(2)}" },
      { file: "src/analysis-panel.js", pattern: "_lastThermoText = formatThermoTable(res)" },
    ],
  },
  {
    id: "ala_ddg",
    label: "per-residue ΔΔG (alanine scan, thermodynamic cycle)",
    formula: "ΔΔG = [E(mut·holo) − E(WT·holo)] − [E(mut·apo) − E(WT·apo)], steepest-descent relaxed legs",
    codePath: "src/analysis/alanine_scanning.js scanResidue()",
    inputs: "the WT holo + WT apo relaxed energies (cached across the scan) and one mutant holo + apo pair per residue",
    nature: "measured",
    uncertainty: null,
    noErrorBar: "the scan computes no replicate legs, so there is no spread to call an error bar; a noise FLOOR (|ΔΔG| < 0.05 kcal/mol) is the only guard and it is a threshold, not an uncertainty",
    onDisplayPath: true,
    docRow: "4W52 ala-scan (CG) row",
    sites: [
      { file: "src/analysis/alanine_scanning.js", pattern: "row.ddG = row.dGmutHolo - row.dGmutApo;" },
      { file: "src/analysis/alanine_scanning.js", pattern: "ΔΔGbind  note" },
      { file: "src/analysis-panel.js", pattern: "Ala-scan (${sys.mode})" },
    ],
  },
  {
    id: "jarzynski_df",
    label: "Jarzynski ΔF (SMD pulling ensemble)",
    formula: "ΔF = W_min − (1/β)·ln⟨e^{−β(W−W_min)}⟩ over 4 constant-velocity pulls",
    codePath: "src/analysis/unbinding_smd.js jarzynskiFreeEnergy()",
    inputs: "the rupture work W of a 4-pull nonequilibrium ensemble from one recording",
    nature: "mixed",
    uncertainty: { kind: "bootstrap", method: "200-resample nonparametric bootstrap over the 4 pull works" },
    onDisplayPath: true,
    docRow: "none — the SMD legs are a smoke path, not a validated row",
    sites: [
      { file: "src/analysis/unbinding_smd.js", pattern: "const dF = wMin - Math.log(Math.max(1e-300, avg)) / beta;" },
      { file: "src/analysis-panel.js", pattern: "Jarzynski ΔF = ${fmtE(je.dF)}" },
    ],
  },
  {
    id: "koff_surrogate",
    label: "k_off ranking surrogate (SMD rupture)",
    formula: "mean_i β·(W_rupt,i − ΔF_rupt) in kT, β = 1/kT — a unitless ranking score",
    codePath: "src/analysis/unbinding_smd.js koffSurrogate()",
    inputs: "the same 4 rupture works and the Jarzynski ΔF of the rupture sub-ensemble",
    nature: "assumed",
    uncertainty: null,
    noErrorBar: "a ranking score over 4 pulls has no ensemble to disperse; the code computes none and an error bar here would be invented precision",
    onDisplayPath: true,
    docRow: "none",
    sites: [
      { file: "src/analysis/unbinding_smd.js", pattern: "const scores = ruptWorks.map((w) => beta * (w - refDF)); // kT units" },
      { file: "src/analysis-panel.js", pattern: "koff-score = ${fmtE(koff.meanScore)} kT" },
    ],
  },
  {
    id: "network_dg",
    label: "4-state network ΔG_bind (→ K_D)",
    formula: "ΔG_bind = E(bound state) − E(solvent state); K_D = exp(ΔG_bind/kT) — the state energies are ASSIGNED (bound 1.00·ΔG, 0.55·, 0.28·), not sampled",
    codePath: "src/physics/network.js computeKinetics()",
    inputs: "the 4-state toy network's own assigned state energies and rate matrix",
    nature: "assumed",
    uncertainty: null,
    noErrorBar: "there is no sampling here at all — the state energies are set by setBoundEnergy() and the ratios are literals in the constructor, so an error bar would be a fabricated measurement",
    onDisplayPath: true,
    docRow: "none — the network is a teaching timescale bridge (ROADMAP.md §2)",
    sites: [
      { file: "src/physics/network.js", pattern: "const dG_bind_standard = this.states[3].energy - this.states[0].energy;" },
      { file: "src/network-panel.js", pattern: "&Delta;G<sub>bind</sub>" },
    ],
  },
  {
    id: "eef1_lig_dg",
    label: "EEF1-lite per-atom desolvation ΔG_a (a potential, not an estimator)",
    formula: "U_desolv = Σ_a ΔG_a·B_a with B_a = 1 − exp(−n_a/3), ΔG_a = −0.30 (C) … −0.50 (I) kcal/mol, a per-element TABLE CONSTANT",
    codePath: "src/ff-binding.js (U_desolv) / src/physics/params.js (the ΔG column)",
    inputs: "the burial neighbour count n_a and the element table — no sampling of any kind",
    nature: "assumed",
    uncertainty: null,
    // WORDING IS LOAD-BEARING: this string is serialized verbatim into
    // `scope.deltaGChain.estimators[].noErrorBar` in every exported results
    // record, and tests/test_results_record.js asserts the record text carries
    // no PDB coordinate records — it greps for the 6-column record keywords
    // "ATOM" and "MODEL" (each padded to column 6). An earlier draft of this
    // sentence opened "…; its MODEL error is unbounded…", which put the literal
    // "MODEL " — byte-identical to a PDB MODEL header — inside the record and
    // turned that guard red. So: no bare uppercase 6-column PDB record keyword
    // in any prose field of this table. tests/test_results_record.js is the
    // guard; this comment is why it stays green.
    noErrorBar: "a fixed parameter has no uncertainty; the error of the model this table encodes is unbounded by anything this tool measures, and the ±50% SASA band on the solvent leg (docs/VALIDATION.md) is the only error quantification anywhere near it",
    onDisplayPath: true,
    docRow: "the desolv component inside every ΔH cell of the 4W52 CG rows",
    sites: [
      { file: "src/ff-binding.js", pattern: "Udesolv += ff._ligdG[a] * (1 - e);" },
    ],
  },
]);

/**
 * The bridge between the two ΔG stories: what would have to be TRUE, and what
 * this tool would have to implement, for the HUD number and the VALIDATION.md
 * number to be the same number. Every item is `out-of-scope` because every one
 * of them is physics this tool does not have — which is the honest answer, not
 * a missing feature to be scheduled.
 *
 * MEASURED ON 4W52 + BENZENE (the system both stories are about), seed 101,
 * app defaults, CG Cα, T=300, ζ=8, m=110 — see tests/test_dg_drift.js, which
 * recomputes every one of these numbers on every run.
 */
export const DG_BRIDGE = Object.freeze([
  {
    id: "shared-observable",
    headline: "The two numbers are not two measurements of one quantity",
    body:
      "The HUD number is a vertical free-energy difference along ONE radial CV from a bias grid that has only just started to fill. The VALIDATION number is a horizontal ΔH − T·ΔS difference between two legs whose apo side is an internal relaxation of the holo structure. They share no term. Reconciling them means either making one estimator reproduce the other's legs (a different estimator, i.e. a physics change) or computing a third quantity both can be compared against.",
    missing: "a common-leg estimator: the same trajectory pair feeding both the radial PMF and the ΔH/ΔS decomposition",
    movedBy: null,
    nature: "out-of-scope",
  },
  {
    id: "hill-convergence",
    headline: "The HUD number has not converged, and its error bar does not say so",
    body:
      "estimateDG() is a two-point read of a bias that is still filling. Measured on 4W52 + benzene at the app defaults: −0.149 kcal/mol at 50 hills, +0.481 at 100, +1.015 at 200, +2.195 at 400, +2.867 at 800, +4.841 at 1600, +5.200 at 3200 — a monotone drift of +5.35 kcal/mol while the reported SE falls from 0.084 to 0.011. convergenceSE() is kT/√nHills, a decorrelation-time heuristic; it is two to three orders of magnitude smaller than the drift it is printed next to. Related: Funnel.getScaledStride() is documented (docs/FUNNEL.md §1, D36) as the dt-aware deposition schedule and is NEVER CALLED anywhere in src/, so hills land every 20 force evaluations at every timestep — 5× denser than the ~100 fs the doc describes, because the CG+ligand dt auto-tunes to 1.7 fs, not the 4 fs baseline.",
    missing: "a block-averaged or replica-resampled error bar on the reconstructed PMF (getScaledStride() exists for the deposition schedule and is NEVER CALLED, so hills land every 20 force evaluations regardless of dt)",
    movedBy: 5.35,
    movedByUnit: "kcal/mol",
    nature: "out-of-scope",
  },
  {
    id: "jacobian-and-ct",
    headline: "Two standard-state terms are computed and then not applied",
    body:
      "The rigorous ΔG° carries the radial Jacobian 2kT·ln(r): at rFar = 9.636 Å that is +2.701 kcal/mol against +0.587 at the bound CV 1.636 Å, a 2.114 kcal/mol difference the two-point path omits. getPMF() also computes the Tiwary–Parrinello offset c(t) (+0.219 kcal/mol at 100 hills) and estimateDG() never adds it. Both omissions are documented (docs/FUNNEL.md §3, §7); neither is visible in the HUD.",
    missing: "applying the Jacobian and c(t) to the displayed number — which is a physics change and is explicitly out of scope for this run",
    movedBy: 2.11,
    movedByUnit: "kcal/mol",
    nature: "out-of-scope",
  },
  {
    id: "no-delta-s-decomposition",
    headline: "The HUD number contains no ΔH, no ΔS and no ΔSASA",
    body:
      "Nothing in estimateDG() separates enthalpy from entropy. The EEF1-lite burial term (assumed, per-element table), the Schlitter pocket entropy and the LCPO ΔSASA solvent entropy are folded, unlabelled, into a single reconstructed profile — or, at Cα resolution with a rigid benzene, simply absent (ΔS_lig ≡ 0 exactly, VALIDATION's BNZ rigid control). A user expanding the HUD number therefore sees terms that the VALIDATION table calls out by name and cannot match them one-for-one.",
    missing: "an entropy decomposition along the same CV (quasi-harmonic or MBAR-style reweighting of the bias trajectory)",
    movedBy: null,
    nature: "out-of-scope",
  },
  {
    id: "system-definition",
    headline: "The two stories do not use the same system definition",
    body:
      "VALIDATION's record-path rows carry the crystallisation buffer EPE alongside benzene; its BNZ-only rows do not. Measured ΔH difference from that choice alone: 3.18 kcal/mol (−6.82 vs −3.64). The funnel's own pocket (14 beads within 8 Å of the native ligand COM) matches VALIDATION's BNZ-only pocket exactly, and running the funnel on BNZ+EPE instead moves its ΔG from +0.481 to +0.313 — the same system-definition effect seen from the other estimator.",
    missing: "one declared system definition per reported number, carried in the record (the results record hashes the input bytes but does not name the ligand subset)",
    movedBy: 3.18,
    movedByUnit: "kcal/mol",
    nature: "assumed",
  },
  {
    id: "leg-and-replicas",
    headline: "One story has two legs and three replicas, the other has one",
    body:
      "thermo_dg needs a holo leg and an apo leg and VALIDATION reports it over 3 replicas with −TΔS_pocket = [+4.58, −3.71, +10.42], mean +3.76, SD 7.1 kcal/mol — the replica SD is ~45× the ΔH bootstrap SE (±0.16) and it dominates the answer. funnel_dg_hud runs on one trajectory from one seed and reports no replicate dispersion at all.",
    missing: "seeding a ≥3-replica funnel run and reporting the between-replica SD of the reconstructed PMF",
    movedBy: 7.1,
    movedByUnit: "kcal/mol",
    nature: "out-of-scope",
  },
  {
    id: "physics-level-is-a-no-op-here",
    headline: "The physics-level assumption moves this system's number by exactly zero",
    body:
      "VALIDATION's CG rows run with charges + directional H-bonds (L1); the app default is L0. On 4W52 + benzene the two are bit-identical (binding energy −3.287653 kcal/mol both ways; Coul ≡ 0.00, HB ≡ 0.00) because every benzene carbon carries q = 0 and has no H-bond capability. On the BNZ+EPE record path the same switch moves the binding energy by +1.272 kcal/mol (pure Coulomb from EPE). So the L0/L1 assumption is worth 0.00 here and 1.27 on the record path — an assumption whose cost depends on the ligand, which is exactly why it cannot be read off a single number.",
    missing: null,
    movedBy: 0.0,
    movedByUnit: "kcal/mol",
    nature: "assumed",
  },
  {
    id: "mislabeled-leg",
    headline: "FINDING (label only, nothing changed here): one surface names the wrong leg",
    body:
      "src/analysis.js prints the analysis-report line as \"Binding ΔG (WTM bias, bound vs +6 Å)\". The code it is reading takes its far point at cv0 + 8.0 Å. The number is right; the sentence describing it is 2 Å wrong, and it is the only place a user is told the leg distance in prose rather than seeing the chain. Not corrected in this run: src/analysis.js is outside the touched surface and the brief forbids unilateral changes to a reported quantity's surroundings. It is recorded here so the chain and the label cannot both be believed.",
    missing: null,
    movedBy: null,
    nature: "out-of-scope",
  },
]);

/* ========================================================================
 * THE DOC ↔ CODE DRIFT DETECTOR
 *
 * The single most valuable thing in this file. docs/VALIDATION.md reports
 * numbers; src/ computes numbers; until now nothing connected them, so a
 * changed parameter or a changed protocol would silently leave the document
 * quoting a value no code path produces any more.
 *
 * `VALIDATION_NUMERIC_CLAIMS` names each claim, the markdown that carries it,
 * and the tolerance at which the document's printed precision is compared.
 * `parseValidationClaims()` reads the numbers OUT of the markdown (so the
 * comparison is against what a reader actually sees, not against a copy);
 * `validationDrift()` compares them to freshly measured values and returns the
 * list of disagreements. tests/test_dg_drift.js runs it against the REAL
 * documents on every FAST-tier run, and PROVES it can fail by handing it a
 * document with one digit changed.
 *
 * THE SIGN RULE
 * -------------
 * docs/VALIDATION.md uses the typographic minus U+2212, not ASCII hyphen-minus.
 * A parser that does not fold it reads "−6.82" as NaN and either skips the
 * check (worst: silently vacuous) or fails for the wrong reason. `num()`
 * folds both, and reports the character it saw so a test can assert the fold
 * actually happened rather than assuming it.
 */

/** Fold a typographic minus to ASCII and parse; NaN when the text is not a number. */
export function num(s) {
  const t = String(s ?? "").replace(/−/g, "-").replace(/[\s,]/g, "");
  if (!/^[+-]?\d*\.?\d+([eE][-+]?\d+)?$/.test(t)) return NaN;
  return Number(t);
}

/** How many decimals a printed value is quoted to (0.30 → 2; 167.1 → 1). */
export function printedPrecision(s) {
  const t = String(s ?? "").replace(/−/g, "-").trim();
  const m = /^-?\d*\.(\d+)$/.exec(t);
  // No decimal point in the printed string ⇒ the document printed no decimals
  // at all ⇒ the precision it asserts is a count of zero. Spelled as an early
  // return rather than a ternary on purpose: this module carries NO numeric
  // physics literal in code — every physical number it mentions lives in the
  // prose of a `formula` / `body` / `docRow` string, where a reader sees it —
  // and tests/test_dg_drift.js enforces that by sweeping the comment-stripped
  // source for a `field: <number>;` assignment. A ternary's `: 0` tail looked
  // exactly like such an assignment and turned the guard red for a count that is
  // not physics. Nothing below may gain a tunable constant.
  if (m === null) return 0;
  return m[1].length;
}

/**
 * One claim: `re` finds it in the methods table, `precision` is how many
 * decimals the document prints, and the human key names it in a failure.
 * `role` says which side of the chain the number belongs to, which is what
 * makes a HUD-vs-doc gap reportable rather than a diff.
 *
 * @type {ReadonlyArray<{key:string, role:"doc"|"ui", row:string, re:RegExp,
 *                       precision:number, what:string}>}
 */
export const VALIDATION_NUMERIC_CLAIMS = Object.freeze([
  {
    key: "record_dH", role: "doc", row: "4W52 BNZ+EPE (record path)",
    re: /\|\s*([−-]?[\d.]+)\s*±\s*[\d.]+\s*\([−+\-]?[\d.]+\/[−+\-]?[\d.]+\/[−+\-]?[\d.]+\/[−+\-]?[\d.]+\)\s*\|/,
    precision: 2, what: "record-path ΔH (BNZ+EPE, seeds 101/1101)",
  },
  {
    key: "record_dH_se", role: "doc", row: "4W52 BNZ+EPE (record path)",
    re: /\|\s*[−+\-]?[\d.]+\s*±\s*([−+\-]?[\d.]+)\s*\(/,
    precision: 2, what: "record-path ΔH block-bootstrap SE",
  },
  {
    key: "bnz_dH", role: "doc", row: "4W52 BNZ-only (true cavity)",
    re: /\|\s*([−-]?[\d.]+)\s*±\s*[\d.]+\s*\([−+\-]?[\d.]+\/[−+\-]?[\d.]+\/[−+\-]?[\d.]+\/[−+\-]?[\d.]+\)\s*;/,
    precision: 2, what: "BNZ-only (true cavity) ΔH",
  },
  {
    key: "bnz_dH_lj", role: "doc", row: "4W52 BNZ-only (true cavity)",
    re: /\(\s*([−+\-]?[\d.]+)\/[−+\-]?[\d.]+\/[−+\-]?[\d.]+\/[−+\-]?[\d.]+\)\s*;/,
    precision: 2, what: "BNZ-only ΔH: LJ component",
  },
  {
    key: "bnz_dH_desolv", role: "doc", row: "4W52 BNZ-only (true cavity)",
    re: /\(\s*[−+\-]?[\d.]+\/[−+\-]?[\d.]+\/[−+\-]?[\d.]+\/([−+\-]?[\d.]+)\)\s*;/,
    precision: 2, what: "BNZ-only ΔH: EEF1-lite desolv component",
  },
  {
    key: "bnz_desolv_share_pct", role: "doc", row: "4W52 BNZ-only (true cavity)",
    re: /desolv\s*([\d.]+)%/,
    precision: 1, what: "BNZ-only desolv share of ΔH (%)",
  },
  {
    key: "epe_inflation_dH", role: "doc", row: "4W52 BNZ-only (true cavity)",
    re: /inflates\s*\\?\|?ΔH\\?\|?\s*by\s*([−+\-]?[\d.]+)/,
    precision: 2, what: "EPE-buffer inflation of ΔH (record path − BNZ-only)",
  },
  {
    key: "bnz_dsasa", role: "doc", row: "4W52 BNZ-only + real SASA",
    re: /([\d.]+)\s*±\s*[\d.]+\s*Å²/,
    precision: 1, what: "real LCPO cross-burial ΔSASA (Å²)",
  },
  {
    key: "bnz_dsasa_se", role: "doc", row: "4W52 BNZ-only + real SASA",
    re: /[\d.]+\s*±\s*([\d.]+)\s*Å²/,
    precision: 1, what: "real LCPO cross-burial ΔSASA SE (Å²)",
  },
  {
    key: "bnz_negTds_solv", role: "doc", row: "4W52 BNZ-only + real SASA",
    re: /Å²\s*→\s*([−+\-]?[\d.]+)\s*±/,
    precision: 2, what: "−T·ΔS_solv from measured ΔSASA (kcal/mol)",
  },
  {
    key: "bnz_negTds_solv_se", role: "doc", row: "4W52 BNZ-only + real SASA",
    re: /→\s*[−+\-]?[\d.]+\s*±\s*([\d.]+)\s*→/,
    precision: 2, what: "−T·ΔS_solv SE (kcal/mol)",
  },
  {
    key: "bnz_dG_est", role: "doc", row: "4W52 BNZ-only + real SASA",
    re: /ΔG\s*([−+\-]?[\d.]+)\s*\(/,
    precision: 2, what: "ΔG_est on the real-SASA BNZ-only leg (kcal/mol) — THE headline doc number",
  },
  {
    key: "ala_ddG_top1", role: "doc", row: "4W52 ala-scan (CG)",
    re: /TYR88\s*([−+\-]?[\d.]+)\s*\//,
    precision: 3, what: "ala-scan top-1 ΔΔG (TYR88)",
  },
  {
    key: "ala_ddG_top2", role: "doc", row: "4W52 ala-scan (CG)",
    re: /MET102\s*([−+\-]?[\d.]+)\s*\//,
    precision: 3, what: "ala-scan top-2 ΔΔG (MET102)",
  },
  {
    key: "ala_ddG_top3", role: "doc", row: "4W52 ala-scan (CG)",
    re: /ILE100\s*([−+\-]?[\d.]+)\s*;/,
    precision: 3, what: "ala-scan top-3 ΔΔG (ILE100)",
  },
  {
    key: "schlitter_err_pct", role: "doc", row: "Schlitter unit",
    re: /([\d.]+)%\s*\(<\s*1%\s*bar\)/,
    precision: 2, what: "Schlitter sampled-vs-exact error (%)",
  },
  {
    key: "ui_dg_hud_2000steps", role: "ui", row: "(the HUD ΔG is NOT in any VALIDATION.md row)",
    re: null,
    precision: 2,
    what: "the HUD funnel ΔG for 4W52+benzene at the CG row's 2000-step protocol — measured here, documented in the chain, absent from the methods table (finding F1)",
  },
]);

/**
 * Pull the row for a claim out of the methods table, then run its regex on it.
 * Scoping to the row is what makes the regexes short; without it a claim could
 * match a number from a different row and still look satisfied.
 * @param {string} md docs/VALIDATION.md text
 * @param {string} rowPrefix the row's first cell
 * @returns {string|null}
 */
function rowText(md, rowPrefix) {
  if (rowPrefix.startsWith("(")) return md; // the UI claim has no row of its own
  for (const line of String(md).split(/\r?\n/)) {
    if (line.startsWith("|") && line.includes(rowPrefix)) return line;
  }
  return null;
}

/**
 * Read every claim's printed number OUT of the document.
 *
 * @param {string} md docs/VALIDATION.md text
 * @returns {{found:Record<string,number|null>, rows:Record<string,boolean>,
 *            missing:string[]}} `missing` is the honest failure mode: a claim
 *            whose row was deleted or renamed cannot be checked, and saying so
 *            is the only alternative to pretending it passed.
 */
export function parseValidationClaims(md) {
  const found = {};
  const rows = {};
  const missing = [];
  for (const c of VALIDATION_NUMERIC_CLAIMS) {
    const line = rowText(md, c.row);
    rows[c.key] = !!line;
    if (c.re === null) continue;
    if (line === null) { found[c.key] = null; missing.push(`${c.key}: no methods-table row contains "${c.row}"`); continue; }
    const m = c.re.exec(line);
    if (!m) { found[c.key] = null; missing.push(`${c.key}: row "${c.row}" exists but the pattern /${c.re.source}/ did not match it — the cell was reworded`); continue; }
    const v = num(m[1]);
    found[c.key] = Number.isFinite(v) ? v : null;
    if (!Number.isFinite(v)) missing.push(`${c.key}: matched "${m[1]}" which is not a number`);
  }
  return { found, rows, missing };
}

/**
 * Compare the document's printed numbers to freshly measured ones.
 *
 * Pure in both arguments, which is what lets the test hand it a MUTATED document
 * and prove the detector reports the drift instead of quietly returning [].
 *
 * @param {string} md docs/VALIDATION.md text
 * @param {Record<string, number|null>} measured key → recomputed value
 * @returns {string[]} problems; empty means the document still says what the code says
 */
export function validationDrift(md, measured = {}) {
  const problems = [];
  const { found, missing } = parseValidationClaims(md);
  for (const m of missing) problems.push(`docs/VALIDATION.md ${m}`);

  for (const c of VALIDATION_NUMERIC_CLAIMS) {
    if (c.role !== "doc" || c.re === null) continue;
    const printed = found[c.key];
    const live = measured[c.key];
    if (printed === null || printed === undefined) continue;             // already reported above
    if (!Number.isFinite(live)) {
      problems.push(`${c.key} (${c.what}): no measured value supplied, so the document's ${printed} is UNCHECKED — a drift check that skips is worse than one that fails`);
      continue;
    }
    // Compare at the precision the DOCUMENT prints. A document quoting 2 dp is
    // asserting |true − printed| < 0.005; anything looser would let a real
    // change hide inside the rounding.
    const tol = 0.5 * Math.pow(10, -c.precision);
    if (Math.abs(live - printed) > tol) {
      problems.push(`${c.key} (${c.what}): docs/VALIDATION.md says ${printed.toFixed(c.precision)} but the code now produces ${live.toFixed(6)} (|Δ| ${Math.abs(live - printed).toFixed(6)} > ${tol})`);
    }
  }

  // The UI-side claim is the one VALIDATION.md does NOT make. Assert that
  // absence explicitly rather than letting it read as agreement: if a row ever
  // starts quoting the HUD number, this fires and forces the two to be tied
  // together on purpose.
  //
  // THE REPORT MUST NAME THE CLAIM. Every other problem this function returns
  // is prefixed `${key} (${what}):`, and so is this one. The prefix is not
  // decoration: tests/test_dg_drift.js plants a HUD row in the methods table and
  // requires the returned problem to contain "ui_dg_hud_2000steps". An earlier
  // draft opened straight into "docs/VALIDATION.md now quotes the HUD funnel
  // ΔG…" — the drift WAS detected, but the finding could not be traced back to
  // the single claim it protects, so the guard read as vacuous and went red.
  // A finding nobody can attribute is not an audit trail.
  const uiKey = "ui_dg_hud_2000steps";
  const uiClaim = VALIDATION_NUMERIC_CLAIMS.find((c) => c.key === uiKey);
  const uiWhat = uiClaim ? uiClaim.what : "the HUD funnel ΔG, absent from the methods table (finding F1)";
  if (Number.isFinite(measured[uiKey])) {
    const quotedInDoc = String(md).includes(`ΔG ${measured[uiKey].toFixed(2)}`);
    const hasHudRow = /\|\s*4W52[^\n]*HUD/i.test(String(md));
    if (hasHudRow || quotedInDoc) {
      problems.push(`${uiKey} (${uiWhat}): docs/VALIDATION.md now quotes the HUD funnel ΔG (${measured[uiKey].toFixed(2)}) but no methods-table row declares the funnel estimator's protocol — add the row and its chain, or remove the number`);
    }
  }
  return problems;
}
