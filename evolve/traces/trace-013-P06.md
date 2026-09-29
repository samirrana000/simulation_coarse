# Evolution 13 — P06 (short)

**Goal:** science=1.0000 is saturated, but the overclaim check reads only 4000 of README.md's 23443 chars and passes on a keyword regex
**Verdict:** ACCEPTED
**Note:** heavy_compute_ms 2.0 was fiction at birth - one-commit file history, 15.17ms measured at that same commit; rebaselined 16.0 with provenance

**Acceptance:** The overclaim check examines 100% of README.md (and the docs it summarises) and asserts the ABSENCE of specific overclaim patterns rather than the presence of a disclaimer phrase. A test plants a known overclaim sentence and asserts science < 1.0.

**Touched:** evolve/evolve.mjs, tests/, README.md
