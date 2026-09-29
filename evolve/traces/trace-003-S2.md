# Evolution 3 — S2 (short)

**Goal:** Delete the two dead stub modules that ship warnings on load
**Verdict:** ACCEPTED
**Note:** stubs had 0 app importers; real defect was docs claiming capabilities that don't exist + H71 checked off on a stub

**Acceptance:** scorer-onnx.js and viewer-gl.js deleted; zero dangling imports (grep + node --check all); zero 'not yet implemented' console.warn in src/; README/CHANGELOG updated to say the capability is absent, not pending.

**Touched:** src/, docs/SCORER.md, docs/VIEWER.md, README.md, CHANGELOG.md
