# Evolution 8 — P01 (short)

**Goal:** scripts/validate_binding_physics_r1.mjs exits 1 and is wired into NO gate that can fail
**Verdict:** ACCEPTED
**Note:** R1 validator was unwired+red; docs/BINDING_PHYSICS_R1.md affirmatively false about implemented weak terms; now structurally asserted and gating

**Acceptance:** `node scripts/validate_binding_physics_r1.mjs` exits 0, OR the file is moved to a tier that `npm test`/CI actually executes (justified by its measured runtime). The failing assertion is fixed at its source — not weakened, not deleted. `node evolve/evolve.mjs gate` goes red if it regresses again.

**Touched:** scripts/, tests/suites.js, src/, docs/
