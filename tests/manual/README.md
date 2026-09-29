# tests/manual/ — real tests that cannot be part of the zero-dependency gate

Everything in this directory is a **real, working test** — none of it is dead
code and nothing here is "obsolete coverage we forgot to delete". These files
are here for one reason only: they need something the project deliberately does
not depend on, so they cannot run inside `npm test` / CI without turning a
zero-dependency repository into a Playwright-installation-dependent one.

`tests/test_suite_registry.js` (registered in the FAST tier) fails if any file
lands in `tests/` that is neither wired into a tier in `tests/suites.js` nor
present here. That is what keeps this directory from becoming a graveyard:
**moving a test here is a decision with a reason attached, and the reason is
written below and asserted by the self-check.**

## browser_test.js

Playwright-driven end-to-end run of the real page: loads `index.html` over a
local `node:http` server, drives the UI (load 4W52, start simulation, toggle
heavy mode, pick a ligand, open the settings modal, exercise the chemical
network panel), captures screenshots into `screenshots/`, and fails on any
browser console error.

**Why it is not in the gate**

- `playwright` is **not** in `package.json` — the project's contract is
  "zero runtime dependencies by design". Adding it to `npm test` would make
  the gate un-runnable on a bare checkout.
- It needs browser binaries (`npx playwright install firefox`) and a
  ~15 s wall clock for one run.
- It is stateful against the live DOM, so a headless run is exactly the class
  of test the wiki's pattern P6 says Node gates cannot substitute for.

**How to run it manually**

```bash
npx playwright install firefox          # one-time, if not already installed
node tests/manual/browser_test.js       # writes screenshots/*.png
```

It is still part of the project's verification story: the wiki skill-impact log
records it running 8/8 stages with zero console errors, and it is the *only*
thing that catches DOM-contract breakage that `scripts/wikiskill_gate.js`
cannot see from the id-subset check alone. If you change `index.html`
structure or panel nesting, run it.

## Not tests — where else non-gating scripts live

- `tests/generate_golden.js` — the golden-file **generator** (not a test; it
  writes `tests/golden/4w52_10steps.json`). The test that *reads* that file,
  `test_golden.js`, is in the FAST tier.
- `scripts/pareto_bench.mjs` — benchmark, writes `docs/pareto_frontier.csv`.
- `scripts/smoke_browsers.mjs` — cross-browser Playwright smoke; same
  dependency reason as `browser_test.js` (run
  `node scripts/smoke_browsers.mjs`).
- `scripts/fetch_coreset.mjs` — offline data fetch/restore, not a test.
- `scripts/check.sh`, `scripts/wikiskill_gate.js` — gates, not tests.
- `scripts/validate_binding_physics_r1.mjs` — doc-citation checker for
  `docs/BINDING_PHYSICS_R1.md` §0. **It currently fails** (see
  `evolve/reports/test-inventory.md`); it is left in place as a standing
  finding rather than deleted, because the doc it validates still exists.
