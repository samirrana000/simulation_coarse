# Skill Impact Log — proposals, verdicts, and the numbers that decided them

Per WikiSkill: skills may be rolled back; this log never is. Every proposal
cites the wiki pattern that motivated it and the gating outcome.

---

## Proposal 2026-09-02 #1 — CREATE skill `ui-ux-decluttering`

- **Motivated by:** P1 (identity-blind coloring), P2 (dead canvas), P3
  (unsurfaced value), P4 (feature hierarchy), P5 (UI contract fragility).
- **Proposal:** New executable skill codifying the ranked-flow UI rules:
  progressive disclosure, three-canvas-states rule, value-diff shipping,
  DOM-contract test before refactor, color-class-before-element.
- **Validation plan:** baseline R_best = tests/test_all.js 32 PASSED (measured
  2026-09-02). Gate = full test suite + node --check on all touched files.
- **Verdict:** ACCEPTED
- **Gating numbers (2026-09-02):**
  - node --check: 35 src files clean
  - tests/test_all.js: 32 PASSED, 0 FAILED (baseline held)
  - DOM contract: 62/62 ui ids present in index.html
  - Hotkey contract: Digit1-7 = 7 top-level panels ✓
  - Full suite: 29/29 tests/test_*.js PASS
  - New regression test: tests/test_ligand_colors.js 8/8 PASS
  - Browser E2E (Playwright/Firefox): 8/8 stages, console errors NONE
  - Pixel verification (headless render): ligand vivid 3293 px vs protein
    grey 20442 px (distinct), legend visible, network canvas 5277 colored
    px (drawing), toolbar permanent checkboxes 0, default-open panels 2
  - Headless physics: node cli.js 4W52 10 steps, finite energies (11.907
    kcal/mol final), RMSD 0.043 Å
- **Notes:** one REJECTED-then-fixed sub-attempt during execution: initial
  panel restructure left orphaned `</div></details>` (index.html) that
  broke the sidebar layout — caught by the browser test, not by node-only
  tests. Lesson recorded as new pattern P6.

---

## Proposal 2026-09-02 #2 — CREATE pattern P6 (wiki addition, not a skill)

- **Motivated by:** the browser-test catch above.
- **Verdict:** ACCEPTED (wiki additions are ungated; they only require
  file:line evidence + measured numbers, both present in P6).

