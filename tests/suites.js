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
  { file: "scripts/test_bindviz.mjs", tier: "FAST" },
  { file: "scripts/test_bindlog.mjs", tier: "FAST" },
  { file: "scripts/test_bindlog_integration.mjs", tier: "FAST" },
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
  { file: "tests/test_golden.js", tier: "FAST", note: "10-step golden trajectory regression" },
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
  { file: "tests/test_rev3_issue1_heavy_physics.js", tier: "FAST" },
  // the registry's own dogfood: proves every tests/test_*.js is accounted for
  { file: "tests/test_suite_registry.js", tier: "FAST" },

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
