# Test inventory — what actually gates, and what pretends to

Measured 2026-09-30 on linux, node v22.23.0, repo root as cwd. Every runtime
below is a real `node <file>` wall clock (process spawn included), taken with a
90 s timeout; `TIMEOUT` means the run was killed at the cap.

**How to reproduce:** `node /tmp/…/measure` equivalent — or simply
`for f in tests/test_*.js; do /usr/bin/time -f "%e" node "$f"; done`.

## The defect this closes

Before this change, `tests/test_all.js` ran **17** of the **55** test scripts
in the repository. The other 38 either never ran anywhere, or ran only via a
handful of one-off CI steps. A test file that exists but never runs is a false
claim of coverage: it looks rigorous, it gets cited in the wiki, and it gates
nothing. Two of the unwired files had in fact **already been deleted** in the
working tree by a concurrent change while remaining listed as "29/29 passing"
in `WHAT_CHANGED.md` — a concrete demonstration of the rot.

The second half of the defect was the harness itself: `FAST_SUITES` and
`SLOW_SUITES` hardcoded `expect: <n>` per suite, so any suite that gained an
assertion produced a bare `expected 20/0` failure that looked like a physics
regression. Both halves are now fixed (see "Wiring" and "Expect fragility").

## Inventory

`wired` = runs in a tier declared in `tests/suites.js`. `other gate` = runs
somewhere else (CI yaml, a wikiskill gate, an evolve gate).

### tests/test_*.js — the suite directory

| file | asserts | runtime | wired | other gate | tier |
|---|---:|---:|:---:|---|---|
| test_all.js | 32 inline + children | 16.6 s | harness | npm test, CI, wikiskill_gate, evolve gate | — |
| tests/test_suite_registry.js | 12 | 0.05 s | **y** | — | FAST (new) |
| tests/suites.js | — (registry) | — | n/a | — | — |
| test_ala_noise_floor.js | 16 | 0.6 s | y | — | FAST |
| test_altloc_cleaner.js | 19 | 0.04 s | y | — | FAST |
| test_b_factors.js | 1 (PASS gate) | 0.3 s | **y** (was unwired) | — | FAST |
| test_bond_dist.js | 1 (PASS gate) | 0.1 s | **y** (was unwired) | — | FAST |
| test_charges.js | 20 | 0.1 s | y | — | FAST |
| test_ck.js | 8 | 0.03 s | **y** (was unwired) | — | FAST |
| test_classify.js | 8 | 0.1 s | **y** (was unwired) | — | FAST |
| test_cutoff.js | 1 (PASS gate) | 0.03 s | **y** (was unwired) | — | FAST |
| test_dt.js | 1 (PASS gate) | 0.03 s | **y** (was unwired) | — | FAST |
| test_enm_seq.js | 1 (PASS gate) | 0.2 s | **y** (was unwired) | — | FAST |
| test_exclusions.js | 11 | 0.03 s | **y** (was unwired) | — | FAST |
| test_forces_fd.js | 1 (PASS gate) | 1.0 s | **y** (was CI-only) | CI `check.yml` | FAST |
| test_funnel.js | 10 | 0.03 s | **y** (was unwired) | — | FAST |
| test_funnel_grid.js | 3 | 0.03 s | **y** (was unwired) | — | FAST |
| test_gb_fd.js | 1 (PASS gate) | 0.03 s | **y** (was CI-only) | CI `check.yml` | FAST |
| test_gillespie.js | 1 (PASS gate) | 0.03 s | **y** (was unwired) | — | FAST |
| test_golden.js | 1 (PASS gate) | 0.1 s | **y** (was unwired) | — | FAST |
| test_gpu_clamp.js | 1 (PASS gate) | 0.02 s | **y** (was unwired) | — | FAST |
| test_hbond.js | 8 | 0.03 s | **y** (was unwired) | — | FAST |
| test_input_errors.js | 75 | 0.04 s | y | — | FAST |
| test_l0_default_exposure.js | 40 | 0.1 s | **y** (was unwired) | — | FAST |
| test_ligand_colors.js | 30 | 0.1 s | **y** (was unwired) | — | FAST |
| test_live_tracking.js | 2 | 0.1 s | **y** (was unwired) | — | FAST |
| test_master_stability.js | 3 | 0.03 s | **y** (was unwired) | — | FAST |
| test_mol2_fidelity.js | 15 | 0.03 s | **y** (was unwired) | — | FAST |
| test_negative.js | 3 | 0.2 s | **y** (was unwired) | — | FAST |
| test_nve.js | 1 (PASS gate) | 0.1 s | **y** (was unwired) | — | FAST |
| test_parity.js | 1 (PASS gate) | 0.1 s | **y** (was unwired) | — | FAST |
| test_picking.js | 31 | 0.03 s | **y** (was unwired) | — | FAST |
| test_placement_hetero.js | 16 | 0.03 s | **y** (was unwired) | — | FAST |
| test_pmf_export.js | 12 | 0.03 s | **y** (was unwired) | — | FAST |
| test_provenance.js | 8 | 0.03 s | **y** (was unwired) | — | FAST |
| test_recorder_cap.js | 1 (PASS gate) | 0.02 s | **y** (was unwired) | — | FAST |
| test_rev1_issue4_live_terms.js | 59 | 0.4 s | **y** (was unwired) | — | FAST |
| test_rev1_issue5_rmsd_split.js | 14 | 0.03 s | **y** (was unwired) | — | FAST |
| test_rev2_issue1_physics_level.js | 16 | 0.1 s | **y** (was unwired) | — | FAST |
| test_rev2_issue5_placement_escape.js | 14 | 0.03 s | **y** (was unwired) | — | FAST |
| test_rev3_issue1_heavy_physics.js | 22 | 0.2 s | **y** (was unwired) | — | FAST |
| test_rotbonds.js | 17 | 0.04 s | y | — | FAST |
| test_seeded_integrator.js | 16 | 0.2 s | y | — | FAST |
| test_session_roundtrip.js | 46 | 0.04 s | y | — | FAST |
| test_topology.js | 4 | 0.03 s | **y** (was unwired) | — | FAST |
| test_viewer_view_state.js | 24 | 0.05 s | **y** (new file from concurrent work) | — | FAST |
| test_virtual_sites.js | 8 | 0.1 s | y | — | FAST |
| test_weakint.js | 49 | 11.0 s | y | — | FAST (slowest entry) |

`test_rev2_issue3_viewergl_sync.js` and `test_rev3_issue3_viewergl_view_sync.js`
existed on disk and in `git ls-files` at the start of this task and were
deleted mid-task by a concurrent change (they imported the now-deleted
`src/viewer-gl.js`). They had **never** been wired into any tier, so nothing
noticed the deletion — the exact failure mode this inventory exists to stop.
They are not counted in the totals above.

### tests/ — non-`test_*` files

| file | kind | runtime | note |
|---|---|---:|---|
| tests/generate_golden.js | **tool**, not a test | 0.05 s | writes `tests/golden/4w52_10steps.json`; the reader (`test_golden.js`) is in FAST |
| tests/golden/4w52_10steps.json | fixture | — | consumed by `test_golden.js` |
| tests/manual/browser_test.js | **test**, manual | 14.5 s | Playwright; moved out of the zero-dep gate, see `tests/manual/README.md` |

### scripts/*.mjs

| file | asserts | runtime | wired | other gate | tier |
|---|---:|---:|:---:|---|---|
| test_bindlog.mjs | 18 | 0.1 s | y | — | FAST |
| test_bindlog_integration.mjs | 14 | 0.1 s | y | — | FAST |
| test_bindviz.mjs | 22 | 0.04 s | y | — | FAST |
| test_pocket_entropy.mjs | 16 | 19.3 s | **y** (was unwired) | — | **MEDIUM** (new tier) |
| test_thermo.mjs | 7 | 2.5 s | y | — | SLOW |
| test_thermo_heavy.mjs | 13 | 10.2 s `--smoke` / ~200 s FULL | y | — | SLOW (TIMEOUT at 90 s without `--smoke`) |
| calibration_4w52.mjs | 14 | 1.9 s | y | — | SLOW |
| validate_flexlig.mjs | 10 | 0.5 s | y | — | SLOW |
| validate_1crn_null.mjs | 8 | 0.1 s | y | — | SLOW |
| **validate_binding_physics_r1.mjs** | — | **FAILS, exit 1** | **n** | **none** | see Finding 1 |
| pareto_bench.mjs | tool | 23–30 s | n/a | — | benchmark, writes `docs/pareto_frontier.csv` |
| smoke_browsers.mjs | manual test | 9.8 s | n/a | — | Playwright, not in gate |
| fetch_coreset.mjs | tool | 0.3 s | n/a | — | offline data fetch/restore |
| wikiskill_gate.js | gate | ~17 s | n/a | CI `check.yml` | DOM + regression gate |
| check.sh | gate | ~3 s | n/a | CI `check.yml` (`npm run check`) | `node --check` over all of `src/` |

## Tier contract

| tier | when | contents | budget |
|---|---|---|---|
| **FAST** | always (`npm test`) | 49 unit-scale suites, no external services | **measured 16.2 s**, contract ≤ 60 s, asserted by the harness at run time and by `test_suite_registry.js` |
| **MEDIUM** | `--medium` / `MEDIUM=1`, implied by `--slow` | 1 suite: `test_pocket_entropy.mjs` (19.3 s seeded pilot) | 10–60 s |
| **SLOW** | `--slow` (SMOKE) / `--slow-full`\|`--long`\|`--full` (FULL) | 5 validation legs: thermo, thermo_heavy, calibration, flexlig, 1crn-null | minutes |
| **manual** | by hand | `tests/manual/browser_test.js` | needs Playwright + browser binaries |

`MEDIUM` exists because `test_pocket_entropy.mjs` is real coverage (16 asserts
on pocket-entropy decomposition) but costs 19.3 s — folding it into FAST would
have pushed the gate from 16 s to ~36 s, still inside budget but spending a
quarter of it on one suite. Anything slower than the budget gets promoted to
MEDIUM rather than silently left unwired.

## Expect-count fragility — how it was fixed

**Chosen: auto-discovery from a declarative registry, counts derived at run
time.** The registry (`tests/suites.js`) declares `{file, tier, timeout, args}`
and nothing else. `test_all.js` no longer holds any expected count.

Why this over "advisory warn instead of fail": an advisory count is a warning
nobody reads — it fires on every legitimate test addition, so it becomes noise
and then gets ignored, which is the failure mode we are trying to remove. A
derived count is always correct by construction, and it can never be stale
because it is not stored.

The gate is now: **child exit code is 0, and the child did not self-report
`N FAILED > 0`.** Counts are reporting only, obtained by
`countAsserts()` in this order:

1. `N PASSED, M FAILED` summary line (the modern convention),
2. otherwise the number of `✓` lines,
3. otherwise 1, if the script printed a bare `PASS`,
4. otherwise 0, printed as `⚠ UNCACHED` so a suite that has silently stopped
   asserting is visible instead of contributing a silent zero.

`test_suite_registry.js` asserts the registry contains **zero** `expect`/
`asserts` keys, so the fragile pattern cannot be reintroduced.

Grand totals therefore move freely: adding an assertion can no longer break
the harness, and `wikiskill_gate.js` / `evolve.mjs`'s `>= 352` floor is now a
genuine regression floor rather than a number that must be hand-edited in
three places every time a test grows.

## Self-check (the durable part)

`tests/test_suite_registry.js`, registered in FAST, so it runs in `npm test`
and in CI. It fails if:

- any `tests/test_*.js` or `scripts/test_*.mjs` file exists that is neither
  registered in a tier, present in `tests/manual/`, or matched by an `EXEMPT`
  entry with a written reason;
- a file in `tests/manual/` is not named in `tests/manual/README.md`;
- the registry declares a hardcoded assertion count, duplicates an entry,
  points at a file that does not exist, or uses an unknown tier name;
- `FAST` or `SLOW` is empty, or the FAST budget constant exceeds 60 s;
- an `EXEMPT` entry has a trivial reason or matches nothing;
- any of the five physics-critical suites (`test_gb_fd`, `test_forces_fd`,
  `test_golden`, `test_nve`, `test_exclusions`) is not in FAST — so nobody can
  quietly demote real coverage to `tests/manual/`.

## Findings

**Finding 1 — pre-existing failure, not caused by this change.**
`scripts/validate_binding_physics_r1.mjs` exits 1:

```
AssertionError [ERR_ASSERTION]: heavy lacks halogen
    at scripts/validate_binding_physics_r1.mjs:38
```

It greps `src/heavy.js` for the *absence* of `halogen`/`cation-pi`/etc. —
i.e. it enforces the R1 honest-scope claim that the heavy kernel has no
halogen term. `src/heavy.js` now mentions halogens, so the claim the script
encodes is out of date (either the physics grew, or the doc-comment changed).
It is not wired into any tier and is not in `check.yml`, which is exactly how
it went unnoticed. Left in place as a standing finding, not deleted and not
"fixed" by weakening the assertion — resolving it is a documentation/physics
decision, not a test-harness decision.

**Finding 2 — the 352 baseline is now a floor, not a target.** The grand total
moved 352 → 737 purely from wiring up files that already existed and already
passed. No assertion was removed, added, or weakened; nothing in `src/` was
touched.
