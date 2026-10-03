/**
 * tests/suites.js — the declarative suite registry for tests/test_all.js.
 *
 * WHY THIS FILE EXISTS
 * --------------------
 * The registry lists *which* scripts run in *which* tier. It deliberately
 * does NOT list how many assertions each script makes. The old tables in
 * test_all.js hardcoded `expect: 20` per suite, which meant that adding one
 * assertion to a suite produced a confusing "expected 20/0" failure with no
 * hint that the *harness* was stale, not the test. Counts are now derived
 * from each script's own output at run time (see test_all.js::runSuiteFile):
 *
 *   - `N PASSED, M FAILED`  -> parsed (the modern convention)
 *   - else `✓` lines        -> counted
 *   - else one `PASS` gate  -> counted as 1
 *
 * The gate is the child's EXIT CODE plus any parsed `M FAILED > 0`. Nothing
 * here can make a test fail by being out of date.
 *
 * TIER CONTRACT
 * -------------
 *   FAST   (always)             unit-scale, no external services.
 *                               Budget: must stay under ~60 s — `npm test`
 *                               and every CI gate depend on it. Enforced by
 *                               tests/test_suite_registry.js (FAST_BUDGET_S).
 *   MEDIUM  (--medium, MEDIUM=1, or --slow/--slow-full)
 *                               10-60 s of seeded/ensemble work. Opt-in so
 *                               the FAST gate stays cheap.
 *   SLOW    (--slow, SLOW=1, or --slow-full)
 *                               multi-minute validation legs. --slow runs
 *                               the SMOKE protocol; --slow-full runs FULL.
 *                               MEDIUM is implied by SLOW (it is cheaper).
 *
 * tools/manual: anything that cannot be part of a zero-dependency Node gate
 * (e.g. Playwright) lives in tests/manual/ and is covered by the self-check
 * in tests/test_suite_registry.js, so it cannot rot unnoticed either.
 *
 * NOT tests: tests/generate_golden.js (golden generator), scripts/pareto_bench.mjs
 * (benchmark), scripts/fetch_coreset.mjs (data fetch), scripts/smoke_browsers.mjs
 * (Playwright smoke), scripts/check.sh + scripts/wikiskill_gate.js (gates).
 * They are inventoried in evolve/reports/test-inventory.md as tools.
 */

export const FAST_BUDGET_S = 60;

/** Tiers, cheapest first. */
export const TIERS = ["FAST", "MEDIUM", "SLOW"];

/**
 * @typedef {{file: string, tier: "FAST"|"MEDIUM"|"SLOW", timeout?: number,
 *            args?: string[], note?: string}} Suite
 * @type {Suite[]}
 */
export const SUITES = [
  // ---------------------------------------------------------------- FAST
  { file: "tests/test_charges.js", tier: "FAST" },
  { file: "tests/test_virtual_sites.js", tier: "FAST" },
  { file: "tests/test_seeded_integrator.js", tier: "FAST" },
  { file: "tests/test_altloc_cleaner.js", tier: "FAST" },
  { file: "tests/test_rotbonds.js", tier: "FAST" },
  { file: "tests/test_ala_noise_floor.js", tier: "FAST" },
  { file: "tests/test_input_errors.js", tier: "FAST" },
  { file: "tests/test_session_roundtrip.js", tier: "FAST" },
  { file: "tests/test_results_record.js", tier: "FAST", note: "results-record-v1 round-trip (JSON + CSV), the uncertainty-is-never-zero rule, the run-health gate, and the ROADMAP.md §1 scope-drift detector with fault injection" },
  { file: "scripts/test_bindviz.mjs", tier: "FAST" },
  { file: "scripts/test_bindlog.mjs", tier: "FAST" },
  { file: "scripts/test_bindlog_integration.mjs", tier: "FAST" },
  // honest-scope doc/code citation validator. Registered so it can actually
  // gate: it was previously unwired AND failing (exit 1) on a premise that
  // phase R3 had obsoleted, which is why nobody noticed for so long. 0.02 s,
  // reads 12 files, so FAST is the correct tier.
  { file: "scripts/validate_binding_physics_r1.mjs", tier: "FAST", note: "R1 §0 code citations + R3 weak-interaction presence + live ROADMAP/LIMITATIONS scope guard" },
  // numerically-hard physics checks, each ~0.05-1 s
  { file: "tests/test_gb_fd.js", tier: "FAST", note: "GB analytic force vs central difference" },
  { file: "tests/test_forces_fd.js", tier: "FAST", timeout: 120000, note: "HeavyForceField force-field FD" },
  { file: "tests/test_weakint.js", tier: "FAST", timeout: 180000, note: "~11 s, the single slowest FAST entry" },
  // structural / forcefield invariants
  { file: "tests/test_topology.js", tier: "FAST", note: "disulfide detection, Ca-N rejection" },
  { file: "tests/test_exclusions.js", tier: "FAST", note: "1-4 scaling + intra-ligand exclusion" },
  { file: "tests/test_hbond.js", tier: "FAST", note: "directional H-bond term" },
  { file: "tests/test_bond_dist.js", tier: "FAST", note: "backbone 3.81 A held under Langevin" },
  { file: "tests/test_cutoff.js", tier: "FAST", note: "ENM cutoff monotonicity" },
  { file: "tests/test_dt.js", tier: "FAST", note: "auto-tuned dt per system size" },
  { file: "tests/test_enm_seq.js", tier: "FAST", note: "sequence-weighted ENM springs" },
  { file: "tests/test_nve.js", tier: "FAST", note: "NVE energy conservation" },
  { file: "tests/test_golden.js", tier: "FAST", note: "10-step golden trajectory regression (CG)" },
  { file: "tests/test_heavy_golden.js", tier: "FAST", note: "bit-exact heavy-mode regression: tests/golden/4w52_10steps.json is CG-only and never imports heavy.js, so before this the heavy engine had FD tests at a 1e-3 RELATIVE tolerance and no bit-level net at all. 6 constructor configs x 12 poses x 24 accumulators + a 10-step heavy Langevin run, compared as 4-lane FNV-1a over raw IEEE-754 bytes. COMPARE-only; tests/generate_heavy_golden.js regenerates and is not a test." },
  { file: "tests/test_mol2_fidelity.js", tier: "FAST", note: "MOL2 united-atom handling" },
  { file: "tests/test_b_factors.js", tier: "FAST", note: "B-factor <-> contact Pearson R" },
  { file: "tests/test_parity.js", tier: "FAST", note: "determinism / headless parity" },
  { file: "tests/test_negative.js", tier: "FAST", note: "zero-charge control must differ" },
  { file: "tests/test_pmf_export.js", tier: "FAST", note: "PMF CSV header provenance" },
  { file: "tests/test_funnel.js", tier: "FAST", note: "well-tempered funnel sanity" },
  { file: "tests/test_funnel_grid.js", tier: "FAST", note: "deposition grid quadrature" },
  { file: "tests/test_gpu_clamp.js", tier: "FAST", note: "GPU sigma-over-r clamp" },
  { file: "tests/test_recorder_cap.js", tier: "FAST", note: "maxFrames memory guard" },
  { file: "tests/test_provenance.js", tier: "FAST", note: "file provenance headers" },
  { file: "tests/test_ck.js", tier: "FAST", note: "detailed balance + CK" },
  { file: "tests/test_gillespie.js", tier: "FAST", note: "Gillespie dwell distribution" },
  { file: "tests/test_master_stability.js", tier: "FAST", note: "master-equation normalisation" },
  { file: "tests/test_classify.js", tier: "FAST", note: "pose classifier thresholds" },
  { file: "tests/test_live_tracking.js", tier: "FAST", note: "live COM/contact tracking" },
  { file: "tests/test_picking.js", tier: "FAST", note: "screen->world pick accuracy" },
  { file: "tests/test_ligand_colors.js", tier: "FAST", note: "class-before-element colouring" },
  { file: "tests/test_viewer_view_state.js", tier: "FAST", note: "view transform / system parity" },
  { file: "tests/test_placement_hetero.js", tier: "FAST", note: "hetero-atom sigma + clash escape" },
  { file: "tests/test_l0_default_exposure.js", tier: "FAST", note: "L0 default exposure surface" },
  { file: "tests/test_rev1_issue4_live_terms.js", tier: "FAST" },
  { file: "tests/test_rev1_issue5_rmsd_split.js", tier: "FAST" },
  { file: "tests/test_rev2_issue1_physics_level.js", tier: "FAST" },
  { file: "tests/test_rev2_issue5_placement_escape.js", tier: "FAST" },
  { file: "tests/test_heavy_obc2_trackterms.js", tier: "FAST", note: "obc2 + trackTerms must NOT silently substitute the HCT force field. Pins the 2026-10-03 TDZ fix in heavy/nonbonded.js: HEADLINE evidence is the ON-vs-OFF property (per-term tracking is an accounting switch, so energy and forces must agree BIT-for-bit with it on or off — they did not, because the catch returned pure-HCT energy over OBC2+HCT forces), plus fault injection proving a throwing kernel's partial force contribution is rolled back and the substitution is loud (ff.physicsFallbacks / ff.lastPhysicsFallback), plus a structural pin that every substitutable catch in src/heavy/ is snapshot/restore-guarded." },
  { file: "tests/test_rev3_issue1_heavy_physics.js", tier: "FAST" },
  { file: "tests/test_observables_parity.js", tier: "FAST", note: "CG/heavy kineticTemp + rmsd parity" },
  { file: "tests/test_unit_contract.js", tier: "FAST", note: "one unit contract, one value (KB_KCAL was defined twice)" },
  { file: "tests/test_constant_ledger.js", tier: "FAST", note: "every SCREAMING_CASE constant has one home; catches a contract constant re-derived under ANY name, including inside a WGSL template string; rule 5 extends the same idea to object-literal PARAMETER TABLES, which rules 1-4 cannot see" },
  { file: "tests/test_element_params.js", tier: "FAST", note: "one element, one sigma/eps/q: CG (forcefield.js) and heavy (heavy.js) resolve byte-identically; re-planting the deleted HEAVY_ELEMENT_DEFAULT fails; pins + prints the live CG-vs-heavy metal coverage gap" },
  { file: "tests/test_cache_contract.js", tier: "FAST", note: "no version query literal on any module edge (121 were hand-typed); sw.js exists, is wired up, and re-fetches with cache:no-store — executed in a node:vm sandbox, not grepped" },
  { file: "tests/test_main_module_size.js", tier: "FAST", note: "main.js stays a composition root: <= 400 LOC (was 1686), <= 20 module specifiers (was 23), no physics/rendering imports, and every controllers/ module its header advertises must exist. A refactor without this guard regrows." },
  { file: "tests/test_module_size.js", tier: "FAST", note: "class-level god-module guard, because a per-file limit only guards the file it names: main.js was split in 2026-09 and heavy.js (1674 LOC) silently inherited the crown in 2026-10 without tripping that guard. Every src/ module <= 600 LOC unless it is in a ratchet budget pinned to its exact measured size; the 8 pre-existing over-size modules may shrink but never grow; src/main.js and src/heavy.js additionally <= 400/200 as facades; largest/median module ratio <= 6x (the ratio evolve/evolve.mjs penalises); and the src/ import graph is acyclic (Tarjan SCC over comment-stripped specifiers). A 700-LOC probe and a 2-module cycle both turn it red." },
  // anti-rot self-checks: registry wiring + docs/ index no-orphans
  { file: "tests/test_suite_registry.js", tier: "FAST" },
  { file: "tests/test_docs_index.js", tier: "FAST", note: "every file under docs/ linked from docs/README.md" },
  // honest-scope gates: the two categories of claim that cannot be falsified
  { file: "tests/test_doc_citations.js", tier: "FAST", note: "every path:line claim in every tracked .md resolves to a real, non-blank line (13 allowlisted, all in dated artifacts)" },
  { file: "tests/test_error_surfacing.js", tier: "FAST", note: "no bare catch in src/ (was 56), the recorder counts and surfaces, and a deliberately-failing operation reaches the top status bar" },
  { file: "tests/test_budget_coverage.js", tier: "FAST", note: "every bench/budget.json key is measured by a named producer or explicitly unmeasured with a reason; no published surface claims an fps number" },
  { file: "tests/test_evolve_planner.js", tier: "FAST", note: "plan emits 15 schema-valid goals whose read_only paths all exist" },
  { file: "tests/test_evolve_gate.js", tier: "FAST", note: "the gate can fail: structure not volume, sees uncommitted work, exits non-zero" },

  // -------------------------------------------------------------- MEDIUM
  { file: "scripts/test_pocket_entropy.mjs", tier: "MEDIUM", timeout: 300000, note: "~19 s seeded pocket-entropy pilot" },

  // ---------------------------------------------------------------- SLOW
  { file: "scripts/test_thermo.mjs", tier: "SLOW", timeout: 600000 },
  { file: "scripts/calibration_4w52.mjs", tier: "SLOW", timeout: 600000 },
  { file: "scripts/validate_flexlig.mjs", tier: "SLOW", timeout: 600000 },
  { file: "scripts/validate_1crn_null.mjs", tier: "SLOW", timeout: 600000 },
  { file: "scripts/test_thermo_heavy.mjs", tier: "SLOW", timeout: 900000, note: "--smoke unless --slow-full" },
];

/**
 * Files that are allowed to exist without being in a tier, and are neither in
 * tests/manual/. Each entry MUST carry a reason — an unjustified exemption is
 * exactly the rot this registry exists to prevent.
 * @type {{match: RegExp, why: string}[]}
 */
export const EXEMPT = [];

/** Suites of one tier, in declaration order. */
export function suitesIn(tier) {
  return SUITES.filter((s) => s.tier === tier);
}
