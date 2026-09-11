---
name: ui-ux-decluttering
description: Rules for decluttering, restyling, or restructuring this simulator's browser UI — ranked progressive disclosure, live-canvas three-state rule, color-class-before-element, DOM-contract safety. Grounded in wiki patterns P1–P5.
---

# UI/UX Decluttering Protocol (anti-AI-slop)

Grounded in `.wikiskill/wiki/patterns/P1..P5`. Apply when touching
index.html, css/style.css, src/ui.js, any *-panel.js, or src/viewer.js
visual paths.

## 1. Rank, never remove (P4)
- Exactly ONE primary flow visible by default: **Load → Build → Run → Record**.
- Everything else: collapsed `<details>` (NO `open` attribute) or moved
  behind an explicit "Advanced" disclosure. Analysis/kinetics/ML/equations
  are collapsed by default.
- No numbered panels ("1 —", "2 —") — numbering implies mandatory sequence
  and reads as generated boilerplate. Order + summary text carries meaning.
- Toolbar: primary actions as buttons; view toggles collapsed into a single
  "View" popover/dropdown, not 6 permanent checkboxes.
- Anti-slop test: state the app's primary action within 5 seconds of first
  load; max 2 accent colors; no decorative gradients/borders on dark theme;
  real content density, not filler cards.

## 2. Three-state canvas rule (P2)
Every live canvas must explicitly render:
1. **no-data** — actionable empty-state text ("Load a structure to begin"),
2. **steady** — redrawn on a low cadence (≤1 Hz) even without state change,
3. **active** — visible transition indicator + one-line caption + legend.
A canvas that only redraws on events is indistinguishable from a dead one.

## 3. Color class before element (P1)
- Semantic classes (protein / ligand / metal / water) get DISTINCT channels
  before element shading. Ligand = distinct saturated palette + white halo
  ring (`i >= nProt`), never the same CPK map as protein.
- Colorblind-safe pairs; test distinctness headlessly (assert colors[i]
  differs from same-element protein color).

## 4. Ship the value diff (P3)
- Any batch of improvements MUST produce a user-scannable changelog
  (WHAT_CHANGED.md), linked from the app header. Numbers-first; separate
  [measured] / [docs & infra] / [honest misses]. Never claim aspirational
  stubs as wins.

## 5. DOM contract first (P5)
- BEFORE restructuring index.html: run `node scripts/wikiskill_gate.js`
  (DOM + hotkey contracts). After restructuring: run it again.
- Every id in src/ui.js must exist in index.html; hotkey Digit1-N range
  must match the count of `#controls > .panel` direct children; update the
  range in src/main.js when panel count changes.
- Keep `if (ui.x)` null-guards — they are load-bearing.

## 6. Verification (always)
1. `node scripts/wikiskill_gate.js` — gate must be OPEN.
2. Extended: full `tests/test_*.js` suite + `node cli.js --pdb 4w52 --steps 10`.
3. Visual: serve `python3 -m http.server 8123`, load 4W52 heavy + benzene,
   confirm ligand visually distinct, network canvas shows three states,
   default view shows only the core flow.
