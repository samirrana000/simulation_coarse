# Evolution 6 — M10 (mid)

**Goal:** Docs index + scientist entry path: 37 md files, no map
**Verdict:** ACCEPTED
**Note:** docs routed by question, 38 files linked exactly once, no-orphan test proven to fail on planted orphan

**Acceptance:** docs/README.md exists, routes by question, every file in docs/ is linked from it exactly once, and the trust-boundary docs (ROADMAP, LIMITATIONS, APPLICABILITY, VALIDATION) are the first thing listed. A test asserts no orphan doc.

**Touched:** docs/README.md, README.md, tests/
