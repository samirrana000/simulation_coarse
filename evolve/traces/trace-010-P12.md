# Evolution 10 — P12 (mid)

**Goal:** 8 declared constants (SCREAMING_CASE = number) are defined in more than one src/ file
**Verdict:** ACCEPTED
**Note:** found live k_B 5e-8 split in thermodynamics and 1.9e-4 Coulomb split in worker+GPU; single-sourced at CODATA 332.06371

**Acceptance:** Every duplicated declared constant is defined in exactly one module and imported by the others, and a parity test asserts CG and heavy agree to < 1e-12 for each one. The duplication scan itself becomes a test so the class cannot regrow.

**Touched:** src/, tests/
