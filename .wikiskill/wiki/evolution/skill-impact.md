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

---

## Proposal 2026-09-12 #3 — Stage-7 Loop-2 closeout (perf honesty + GPU decision + docs)

- **Motivated by:** Pareto heap ±50% GC noise with no allocation tracking;
  undocumented GPU-port decision (bindU/weak not in WGSL); README Loop-2 gap.
- **Proposal:** Additive-only closeout, zero deps, no retuning, no UI changes:
  allocation tracking in `scripts/pareto_bench.mjs` (heapUsed delta per tier +
  optional `--expose-gc` gc bracketing, history-preserving CSV writer);
  accept-CPU GPU decision note; README Loop-2 + Stages 1–7 section.
- **Validation plan:** baseline R_best = tests/test_all.js 215 PASSED (measured
  2026-09-12). Gate = `node scripts/wikiskill_gate.js` (OPEN required).
- **Verdict:** ACCEPTED
- **Gating numbers (2026-09-12):**
  - node --check: touched files clean (`scripts/pareto_bench.mjs` et al.)
  - tests/test_all.js: 215 PASSED, 0 FAILED fast (~13 s); --slow 245/245
  - DOM contract: 98/98 ui ids present in index.html (untouched, no new ids)
  - Hotkey contract: Digit1-7 vs top-level panels ✓
  - Bench FULL x2 + gc run (~25 s each): ms stable ±5% (L0 0.22/L1 0.20-0.21/
    L2 78.8-84.8/L3 76.1-83.1); gc-bracketed heap deltas ≤0.3 MB all tiers
    vs ±15 MB un-bracketed (GC-noise debt closed with numbers)
  - Pose recovery unchanged: 10/10 CG, 9/10 heavy
- **Notes:** full record in `docs/BINDING_LOOP2_DONE.md` §14 (allocation
  method + numbers, GPU accept-CPU verdict, README/wiki refs); no GPU code
  changes; Stages 1–6 debts (heavy verdict NOT-restored, seeding, 215/245
  tiers, altloc+rotbonds, async thermo+persist, 4W52 anchor) all carried.

