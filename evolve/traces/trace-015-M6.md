# Evolution 15 — M6 (mid)

**Goal:** Break up the main.js god module (1685 LOC, 30 fns, 23 imports)
**Verdict:** ACCEPTED
**Note:** main.js 1686->160 LOC, 12 controllers, tick fingerprint identical, size+shape guard registered

**Acceptance:** main.js under 400 LOC and only wiring; at least 4 extracted modules (load, build, hud, session); no behaviour change (test_all.js >= baseline, DOM contract green); each extracted module is independently node --check-able and documented in README project layout.

**Touched:** src/main.js, src/controllers/
