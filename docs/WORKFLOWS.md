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
