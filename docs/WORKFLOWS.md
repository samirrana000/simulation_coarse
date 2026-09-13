# Phase 4 — Translational Workflows

Exposed in PMF & Analysis panel (Mutagenesis / Correlations / Kinetics subpanels); original Analyze/PMF flow intact.
Also fixed latent bug: `analysis-panel.js` was never imported, so Analyze/PMF buttons were dead — now imported in `main.js`.

- `src/analysis/alanine_scanning.js` — `scanResidue/scanPocket/formatMutationTable/buildAlanineMutant`;
  ΔΔG=[G_holo(mut)−G_holo(WT)]−[G_apo(mut)−G_apo(WT)] from steepest-descent minima. CG: springs k×0.7, r0−0.4 Å;
  heavy: delete non-{N,CA,C,O,OXT,CB}. Click-residue helper `resIdFromBeadIndex` ready (unwired to viewer pick).
  4W52: TYR88 +0.00 CG / +0.59 heavy, LEU84 −0.01, ALA99 0.00 (apolar cavity → CG≈0 correct).
- `src/analysis/dccm.js` — `computeDCCM` (Kabsch-aligned, Cii≡1), `renderDCCMHeatmap` (null-safe, `canvas._dccm` picking),
  `highlightCorrelatedPair`. 4W52 164×164: asym 0, mean|C|=0.24.
- `src/analysis/cryptic_pockets.js` — `trackPocketVolume/detectCryptic (mean+1.5σ)/probeAccessibility/convexHullVolume/PocketMetaD`.
  4W52: ⟨V⟩≈619, open 5–7.5%. `addForces` live-ready, not yet in tick loop.
- `src/analysis/unbinding_smd.js` — `runConstantVelocityPull/runConstantForcePull/runPullingEnsemble/jarzynskiFreeEnergy (shift+bootstrap)/koffSurrogate`.
  4 pulls: works [173.8,160.7,157.2,158.4], ΔF=157.95±2.35 (v=4 Å/ps fast protocol, dissipative by design).

Runs synchronously (CG: seconds) — move to workerPool/GPU backend in follow-up.

## Physics fidelity level (Loop-2 S7)

Dynamics & Force Field → "Physics fidelity level (L0 / L1 / L2)" subpanel
(`#physicsLevel`, default L0 = baseline, opt-in hot-rebuild):
L0 fast (CG isotropic, charges off) · L1 balanced (CG charges + directional
H-bonds, ≈1.09×) · L2 full-rigor (L1 + heavy weakint + BindLog capture,
≈1.0× over tier base). L2 keeps BindLog capture on even with the Recording
checkbox off. Full bench: `docs/pareto_frontier.csv`, `docs/BINDING_LOOP2_DONE.md`.

## CG-interactive / heavy-offline (FP5)

S7 decided ACCEPT CPU (no GPU port: Loop-2 binding terms stay CPU-only), so
heavy mode stays ~400× slower per step than CG. The honest workflow is tiered:
CG runs interactively in the tab; heavy runs short/offline with progress, then
exports for analysis. Heavy builds are chunked (topology rows with progress
captions + Cancel; Build/Run disabled during work) — the tab never freezes.

### Expected cost (4W52, `docs/pareto_frontier.csv` 2026-09-12 FULL rows)

| Tier | ms/step | dt (ps) | ≈ steps/s in tab¹ | ≈ ps/s in tab |
|---|---|---|---|---|
| L0 CG | 0.237 | 0.004 (no lig) / 0.0017 (lig) | ~3500 / ~3500 | ~14 / ~6 |
| L1 CG+ | 0.235 | same as L0 | ~3500 | ~14 / ~6 |
| L1 CG+ +BindLog | 0.216 | same as L0 | ~3500 | ~14 / ~6 |
| L2 heavy | 84.8 | 0.001 | ~10 | ~0.01 |
| L3 heavy+R3 | 83.1 | 0.001 | ~10 | ~0.01 |
| L2 full-rigor | 85.5 | 0.001 | ~10 | ~0.01 |

¹ Tick budget `advance(steps, 14)` caps compute at 14 ms/frame; CG fits
~59 steps/frame, heavy fits ~1 step per several frames. Timings are
machine-dependent (S7 §14: ±5% across runs on one machine); ratio ≈ 400:1
holds across machines. L4 heavy+OBC2+RESPA (est 46.6) is unchanged by Loop 2.

### Recommended budgets per tier

| Tier | Use for | Recommended run | Recording |
|---|---|---|---|
| L0/L1 CG | interactive explore, placement, funnel/PMF collecting, thermo holo legs | minutes (10³–10⁴ steps) | stride 2 ps, 200–500 frames (0.4–1 ns) |
| L2 heavy | short relaxations, pocket checks, single-point rigor | ≤ 1000 steps in-tab (seconds–minutes) | stride ≥ 2 ps, ≤ 100 frames |
| heavy long sampling | thermo/SMD/pocket-entropy depth | headless scripts only (`test_thermo_heavy` pattern), never the tab | export, don't record live |

A 1 ns heavy recording (10⁶ steps) is ~1 day at pareto speed — that is the
offline case by design, not a bug.

### Export-then-analyze path (FP4 matrix, one click each)

Trajectory XYZ/PDB/JSON (`#dlBtn` + `#exportFmt`) · PMF CSV (`#anaPmfBtn`) ·
Thermo TXT (`#thermoDlBtn`) · BindLog BLG1 (`#bindlogDlBtn`) · DCCM CSV
(`#dccmDlBtn`) · Settings JSON (`#settingsDlBtn`) · Session JSON v1
(`#sessSaveBtn` + `#sessFile`, counts only — frames stay in memory).
Analyze offline in `notebooks/`; reload exact setups via session files.
