# Trace 2026-09-02 — UI truth audit (user-reported symptoms → causes)

Raw execution record. Write-only. All numbers verbatim from session output.

## Attempted (prior sessions, summarized with their measured outcomes)

1. Blank-render audit (2026-08-29): user saw nothing on load.
   - Cause: missing `}` at src/settings-panel.js:78 (syntax → whole module graph dead).
   - Cause: viewer canvas 0×0 under flex (fixed via getBoundingClientRect +
     ResizeObserver + DPR query, src/viewer.js:266-296).
   - Cause: `#viewerWrap { min-height:240px }` absent (css/style.css:120).
   - Verified: tests/test_all.js 32 PASSED after fixes.
2. 100-item transformation plan executed in subagent waves (38 + 62 items).
   - All items marked [x]; dashboard: CG 0.20 ms, Heavy 14.7 ms/step,
     scale exponent 2.06, B-factor median R 0.27 (< 0.45 target), NVE drift 0.05%,
     forces FD maxRel 7.4e-6, 15/15 checks.
   - Honest stubs present: GPU R>0.999 aspirational; barriers hard-coded heuristic.

## Today's user-reported symptoms (verbatim intent)

- "cant differentiate between ligand and heavy atom molecules in heavy atom
  representation as both coloured in same colouring schemes" → confirmed:
  src/viewer.js:135-139 assigns ELEMENT_COLOR to ALL atoms i < n when heavy,
  so ligand C == protein C visually. Ligand atoms are i >= nProt (known).
- "too many features but what actually improved in 100 steps cant seems to
  understand" → no user-facing changelog exists; value buried in
  TRANSFORMATION_PLAN_100.md / dashboard JSON.
- "what chemical network model doing, no idea about it and visually its not
  showing anything at all" → two causes found:
  (a) updateNetworkPlot() only called on state CHANGE (src/main.js:652-659)
      and once at init when canvas may be display:none (details collapsed at
      load? No — panel 5 has `open`; but initial draw happens before system
      exists, and if ligand never placed, classifyPose never fires → stale).
  (b) No explanation surface in default view; "How the Chemical Network
      Works" is a nested details (index.html:190-202), invisible unless
      expanded; canvas has no axis/legend/empty-state message.
- "too many features ... overcrowding features that the user may not want to
  use" → 7 numbered `<details open>` panels force-feeding every feature;
  toolbar shows 6 checkboxes; equations panel always visible.

## Root-cause synthesis (for Maintainer compilation)

- R1: Rendering identity problem — color semantics conflate molecule class
  (protein vs ligand) with element. Data needed to separate them (nProt) is
  already in scope at color-assignment time.
- R2: Perceived-value problem — changes were real (dashboard numbers) but
  never surfaced as a user-scannable diff. Docs exist; nobody reads docs
  to answer "what changed for me".
- R3: Dead-visualization problem — event-driven-only redraw + no empty-state
  + no legend ⇒ canvas reads as "showing nothing" even when wired correctly.
- R4: Feature-load problem — every feature gets equal visual weight;
  progressive disclosure missing; defaults open everything.

## Validation status at trace close

- tests/test_all.js: 32 PASSED, 0 FAILED (baseline R_best for gating).
- Files changed this trace: .wikiskill/* (bootstrap only, no src changes yet).
