# What changed — and what it's worth

A plain-language answer to "what actually improved?" Full technical log:
`TRANSFORMATION_PLAN_100.md`; verified metrics: `TRANSFORMATION_DASHBOARD.json`.

## Fixed (measured)

| What | Before | After |
|---|---|---|
| App rendered blank on load | syntax error killed all modules; canvas 0×0 | fixed; 32/32 regression tests pass |
| Picking misplaced atoms (motion amplification) | clicks landed Å away from atoms | gain-inverted unproject; <0.5 Å error at 5× amplification |
| Heavy-mode forces wrong | GB & SASA derivative bugs | finite-difference verified to maxRel 7.4e-6 |
| Energy conservation (CG) | uncontrolled drift | NVE drift 0.05% |
| Speed | unmeasured | CG 0.20 ms/step; heavy-atom 14.7 ms/step (≈60 fps for ≤3k atoms) |
| Crash behavior | silent NaN freeze | auto-pause + "NON-FINITE ENERGY" HUD warning |
| Ligand vs protein look (heavy mode) | identical CPK greys | ligand = vivid palette + white halo + on-screen legend |
| Chemical network panel | drew once, looked dead | explains itself; redraws live (1 Hz steady + on every state hop) |

## Added

- Deterministic runs: seeded RNG, provenance headers (T, γ, seed, build) in
  every exported trajectory; SHA-256 manifest for bundled data.
- Testing: 32-assertion regression suite, 31 additional targeted tests,
  benches (PDBbind GB RMSE 1.53 kcal/mol; pose recovery to 0.50 Å).
- Headless runs: `node cli.js --pdb 4w52 --steps 10`, Docker, CI (GitHub Actions).
- Docs: `docs/` (tutorial, viewer, placement, network, performance,
  limitations, applicability), tools (AMBER export), CITATION, MIT license.

## Honest misses (do not oversell)

- B-factor correlation: median R 0.27 vs 0.45 target — flexible-loop regions
  still under-predicted by ENM.
- Cost scaling exponent 2.06 — not the O(N) goal for the heavy-atom mode.
- GPU path is a validated stub (clamped, correct) — not yet a speedup.
- CNM barriers (2.5/3.8/4.5 kcal/mol) are heuristic defaults until you
  enable Funnel bias and let the PMF derive them.

## UI redesign (this release)

- One visible flow: **Load → Build → Run → Record**. Advanced panels
  (dynamics tuning, kinetics, recording, analysis, ML) collapsed by default.
- Toolbar: Run/Reset buttons + one "View" menu instead of 6 checkboxes.
- Equations reference moved out of the default view.
- Ligand is now visually distinct in heavy-atom mode (palette + halo + legend).

## Phase 5 — Anti-slop interface deployment

### [measured]

| Check | Result |
|---|---|
| `node scripts/wikiskill_gate.js` | OPEN — 41 files `--check` clean; 32 PASSED / 0 FAILED (≥ 32 baseline); 90 ui ids ⊆ index.html; Digit1-7 vs 7 `#controls > .panel` (+ equations) |
| Full suite `tests/test_*.js` | 29/29 files exit 0, incl. `test_all.js` 32/32, `test_ligand_colors.js` 8/8 |
| `node cli.js --pdb 4w52 --steps 10` | 10 steps, t = 0.017 ps, final E = 12.13 kcal/mol, T_inst = 265 K, RMSD = 0.044 Å |
| Ligand-distinct (headless assert) | protein C `[180,180,180]` vs ligand C `[255,121,98]` — differ on same element; CG chain palette ≠ ligand palette; unknown ligand element → `[255,159,128]` |
| Default view (static) | 2 `<details open>` (Structure + Model = Load → Build; Run lives in toolbar), 5 collapsed, 0 numbered panels, 1 View popover |
| Canvas three-state (static) | viewer empty ≤ 1 Hz + actionable text + halo legend; network 1 Hz steady + glow + caption; PMF idle ≤ 1 Hz + ΔG/CV caption + legend; DCCM empty + heatmap + caption |
| Slider→param sync | 6 slider + numeric pairs (rc/gamma/temp/fric/mass/motionGain) bidirectional, `Math.fround` Float32, hot-rebuild on input |

### [docs & infra]

- Header keeps the `What changed` link → this file.
- New chrome (all flat `#0B0D0E/#121517/#1B1F22` + 1 px `#2A3036`, no shadows/blur):
  top status bar (`#sysState #topPdb #topEngine #topStep`),
  Metrics HUD (`#metricsHud`: T / Etot / PMF / DCCM),
  canvas caption (`#canvasCaption`), PMF/DCCM captions, bottom dock
  (`#scrub #scrubLabel #cvStrip #hudSpark` + `[Space]/[R]/[M]/[1-7]` hints).
- Phase 1–4 controls preserved and collapsed, not removed: protonation/GAFF2 in
  Structure → Advanced chemistry; r-RESPA/funnel in Dynamics → Advanced sampling;
  ML prior in ML Tier; Ala-scan in Mutagenesis; DCCM in Correlations; SMD/cryptic in Kinetics.
- Hotkeys: `[Space]` run/pause, `[R]` reset, `[C]` record toggle, `[M]` mutagenesis panel,
  `[1-7]` instant (`behavior: auto`) panel toggles; `[Esc]` closes modal / cancels placement.
- Numbered-panel references removed (`Panel 4` → `Dynamics → Advanced sampling`;
  `(panel 1)` → `(Structure panel)`).

### [honest misses]

- Visual check was static/headless only (DOM + canvas-code paths + color assert);
  no live-browser screenshot of 4W52-heavy + benzene this run — load
  `python3 -m http.server 8123`, build heavy 4W52, place benzene, confirm halo + legend manually.
- Metrics-HUD `DCCM` shows recorded-frame count, not live mean|C| (full DCCM per frame is too heavy).
- Top-bar `Engine` reflects the selected backend + fps, not measured FLOPS or WebGPU speedup.
- CV/E strips draw every animation frame but text HUD/metrics stay at 10 Hz by design (DOM thrash guard).
