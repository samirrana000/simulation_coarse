# Evolution 4 — S3 (short)

**Goal:** Close the test gap: 49 test files, only ~13 actually gated
**Verdict:** ACCEPTED
**Note:** 36 dead tests wired, expect-counts replaced by registry + orphan self-check, 352->738, FAST 16.1s/60s

**Acceptance:** An inventory lists every tests/*.js as WIRED or RETIRED with a reason; every WIRED file appears in test_all.js; RETIRED files are deleted or moved to tests/manual/ with a README explaining how to run them; test_all.js auto-discovers suites so a new test file cannot be silently unwired (a test that asserts every test_*.js is either wired or in tests/manual/).

**Touched:** tests/test_all.js, tests/, evolve/reports/test-inventory.md
