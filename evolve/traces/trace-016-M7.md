# Evolution 16 — M7 (mid)

**Goal:** Single source of truth for LJ/charge parameter tables
**Verdict:** ACCEPTED
**Note:** canonical physics/params.js; 0 ulp divergence for 9 elements; dead duplicate deleted; LIVE CG-metal gap measured+pinned

**Acceptance:** One element parameter table in src/physics/params.js; ff-params.js and heavy.js both import it; a test asserts element resolution is byte-identical across CG and heavy paths for every element in the table; no duplicated literal remains.

**Touched:** src/physics/params.js, src/ff-params.js, src/heavy.js, tests/
