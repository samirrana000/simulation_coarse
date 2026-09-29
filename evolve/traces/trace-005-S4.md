# Evolution 5 — S4 (short)

**Goal:** De-duplicate the physics helpers duplicated across forcefield.js and heavy.js
**Verdict:** ACCEPTED
**Note:** found real KB_KCAL 5.03e-8 CG/heavy divergence; observables.js single kernel + parity test

**Acceptance:** src/physics/observables.js exports kineticTemp/rmsdTo/rmsdAll; forcefield.js and heavy.js delegate with no duplicated body; a new test asserts CG and heavy kineticTemp agree on identical input to 1e-12; gate stays >= baseline.

**Touched:** src/physics/observables.js, src/forcefield.js, src/heavy.js, tests/
