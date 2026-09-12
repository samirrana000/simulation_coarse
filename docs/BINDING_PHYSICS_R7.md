# Phase R7 — Binding-Physics Visualization (BindViz)

Loop-1, step 7 (final). Implementation: `src/capture/bindviz.js`. Tests: `scripts/test_bindviz.mjs` (22/22, incl. D1 null-bindlog case).
Additive only; `tests/test_all.js` 32/32; gate OPEN. DOM wiring DONE in Loop 2 (see §4).

## 1. Renderers (pure Canvas2D, no chart libs, anti-slop flat)

All follow the three-state rule: no-data actionable text → steady → active, and are
null-safe (null canvas / empty BindLog → draw empty state or return, never throw — tested).

### (a) Interaction timeline — `renderInteractionTimeline(canvas, bindlog, {tFrom, tTo})`
- Rows = (ligand-atom, residue) contact pairs from contact+/− events, sorted by total
  lifetime; colored runs from form→break; class by mean distance: <3.5 Å H-bond cyan
  `#38BDF8`, <4.6 Å hydrophobic amber `#FBBF24`, else slate. Right column: top-5 contacts
  by lifetime + legend. Time axis in ps.
### (b) Energy decomposition — `renderEnergyDecomposition(canvas, bindlog)`
- Per-term polylines over time (7 R4 terms), legend, enthalpy table top-right with
  per-term means and total ΔH. Verified on synthetic data: ΔH renders −8.22 vs analytic −8.225.
### (c) PMF formation — `renderPmfFormation(canvas, bindlog, {binCount, cvRange, cursorT, hillWidth})`
- Rebuilds the deposition-bias curve from hill events up to `cursorT` (Gaussian sum,
  width 0.5 CV units default), magenta fill @18%, hill count + max-bias annotation,
  supports time animation of PMF formation.
### (d) Pareto frontier — `renderParetoFrontier(canvas, tiers)`
- Log-x ms/step vs top-1 accuracy scatter; frontier line through Pareto-optimal tiers
  (magenta), est-tiers slate. Consumes R5 `docs/pareto_frontier.csv` data.

## 2. Test verification (21 assertions)

- Timeline: 16 labels, ≥4 contact bars, residue-99 row present, ps axis, empty state, null-safe.
- Energy: 28 lineTo/4 paths, ΔH −8.22 ✓ (analytic −8.225), term legends.
- PMF: hills-count annotation, filled curve, max bias −0.32 kcal/mol @CV 3.3, cursorT
  slicing (5→2 hills), empty state.
- Pareto: 4 points + frontier, labels, empty state, null-safe.
- Integration: blob-round-tripped BindLog (R6) renders identical timeline.

## 3. BindLog consumption mapping

| Renderer | Events consumed | Source (Loop-2 wiring) |
|---|---|---|
| timeline | contact+ (1), contact− (2) | ff-binding contact tracker / analysis lifetimes |
| energy | energy (0, a=term id) | R2/R3 per-term accumulators |
| PMF | hill (3) | funnel.js deposit callback |
| (future states) | state (4) | physics/network.js hop callback |
| (future pockets) | pocketVol (5) | cryptic_pockets.trackPocketVolume |

## 4. Wiring — DONE (Loop-2 S6)

One collapsed `<details class="subpanel">` "Binding Insights (BindViz)" inside the
existing PMF & Analysis panel (no new top-level panels; Digit1-7 contract intact).
Canvases `bindvizTimeline`/`bindvizEnergy`/`bindvizPmf` + caption
`bindvizCaption`, registered in `src/ui.js` (gate DOM contract 97 ids), rendered
from `src/main.js` `drawBindviz()`/`bindvizTick()` at ≤1 Hz steady, empty states
painted at startup and on `buildSystem()` reset (three-state P2). Pixel sizing
follows the dccm/pmf `devicePixelRatio` pattern (`bindvizFit`), with the renderers
reading the CSS-pixel size back via `canvas._dpr`. D1 fixed (Loop-1 review):
`renderPmfFormation` null-bindlog guard → actionable empty state; test added
(22/22).

Verified (headless Chromium, 1CRN + 4W52): empty states painted at load; 1CRN
BindLog run → 1239 events / 177 frames, timeline + energy canvases live
(pixel-diversity 74/233 distinct colors), caption switches to
"BindViz · N events · N frames (1 Hz live)."; 4W52 + funnel → 1371 hills +
4521 contact events rendered; Build → all three canvases back to no-data text,
caption restored. Zero page errors. `test_bindviz` 22/22; `test_all` 32/32;
gate OPEN before and after.

## 5. Honest limitations

- Class-by-distance is a proxy; true classes (salt bridge vs H-bond) need the R2/R3
  term ids — renderer already accepts arbitrary term ids for energy; contact class will
  switch to term-tagged events once Loop-2 accumulators exist.
- Timeline rows capped by canvas height; scrolling/virtualization is a Loop-2+ concern.
- PMF cursor animation assumes uniform hill width (documented default).
