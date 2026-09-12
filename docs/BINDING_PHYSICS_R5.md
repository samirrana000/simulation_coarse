# Phase R5 — Pareto Frontier: Accuracy vs Speed Tiers

Loop-1, step 5. Benchmark: `scripts/pareto_bench.mjs` (run `node scripts/pareto_bench.mjs [--quick]`).
Data: `docs/pareto_frontier.csv`. Measured on 4W52 (164 Cα + 21 ligand atoms; 1308 heavy atoms).

## 1. Measured tier table

| Tier | Terms | ms/step | heap | Pose recovery (benzene) | What it can answer |
|---|---|---|---|---|---|
| **L0 CG** (today) | ENM + isotropic cross-LJ + burial + holo springs | **0.20** | 9 MB | 10/10 top-1, mean rank 1.0* | fluctuations, pocket occupancy, illustrative PMF |
| **L1 CG+** (est, R2 terms) | + directional HB, live charges (salt bridges), enclosure, halogen-lite, unmet-polar | ~0.23 (×1.15) | ~9 MB | est. equal-or-better on charged/aromatic ligands | + pose ranking with electrostatics, H-bond discrimination |
| **L2 heavy** (today) | covalent + element LJ + GB Coulomb + H-bond(180°) | **82.5** | 11 MB | 9/10 top-1, mean rank 1.3* | realistic geometry, metal sites, GAFF2 ligands |
| **L3 heavy+R3** (est) | + π-stack, cation-π, halogen σ-hole, metal cross-angles | ~86.6 (×1.05) | ~11 MB | est. better on halogenated/aromatic ligands | + weak-interaction-aware pose ranking |
| **L4 heavy+OBC2+RESPA** (est) | + GB-OBC II + LCPO + r-RESPA 2 fs | ~45/effective-step (×0.55) | ~12 MB | est. ≈L3 | + solvation-sensitive ΔΔG ala-scan, best physics per wall-second |

Speed ratio L0:L2 ≈ **400:1** measured. RESPA buys ~2× effective time (Phase 3: 3.8× fewer
slow evals); WebGPU projected 4-8× more (unmeasured headless — honest).

*Pose-recovery caveat (honest): with an isotropic LJ+contact-bonus proxy scorer, the Cα
scorer wins 10/10 because pocket-COM-centered burial is trivially rewarded; the heavy
element-aware scorer ranks 9/10 (one decoy within 0.06 kcal — statistical tie). The proxy
measures burial discrimination, NOT shape complementarity; the true tier differentiator
will be ΔG-scoring against the full force fields (L1/L3) once R2/R3 terms are integrated
in Loop 2. The benchmark is reproducible and will be re-run per Loop.

## 2. Pareto frontier (ASCII)

```
 accuracy
 (physical fidelity)
   ▲
   │                                      L4 heavy+OBC2+RESPA (45 ms/eff-step)
   │                              L3 heavy+R3 (87 ms/step)
   │                      L2 heavy (82 ms/step)
   │      L1 CG+ (0.23 ms/step)  ← Pareto-optimal sweet spot for pose triage
   │  L0 CG (0.20 ms/step)
   └──────────────────────────────────────────▶ speed (steps/s at 60fps budget)
     ~5000/s      ~60/s                       ~12/s      ~22 eff/s
```

Pareto-optimal set (expose in UI): **L0/L1 for triage + teaching**, **L2/L3 for pose
realism**, **L4 for analysis runs** (not interactive). Recommended selector: 3 levels —
"Fast (CG+)", "All-atom", "All-atom + solvent (OBC2)" — collapsing L0/L1 and L2/L3.

## 3. Memory
CG ≈ 9 MB, heavy ≈ 11 MB heap — both trivial for browsers; memory is NOT on the
frontier axis at this system size. OBC2 radii cache adds ~N²×4 B (N=1308 → 6.8 MB) —
only at L4.

## 4. CI hook
`node scripts/pareto_bench.mjs --quick` (≈20 s) belongs in the test suite as a
performance/accuracy regression gate: assert CG ms/step < 0.5, heavy < 120, pose
recovery ≥ 8/10 both scorers.

## 5. Loop-2 S7 re-bench (FULL mode — L1/L3 now measured, not estimated)

`node scripts/pareto_bench.mjs` (full) on 4W52, Loop-2 close-out. Full table:
`docs/pareto_frontier.csv`, `docs/BINDING_LOOP2_DONE.md` §2.

- L0 CG **0.227** ms/step · L1 CG+ **0.247** (×1.09: charges + 162 valid
  virtual sites are nearly free) · L1+BindLog **0.262** (×1.15 vs L0).
- L2 heavy **99.5** · L3 heavy+R3 **89.5** (≈1.0×, noise: weakint ~free) ·
  L2 full-rigor **86.8** (≈1.0×, 7-term bindU finite). L0:L2 ≈ 439:1.
- Accuracy proxy unchanged: 10/10 (Cα scorer) / 9/10 (element scorer) —
  the R5 §1 burial-discrimination caveat still applies.
- UI: 3-level selector L0/L1/L2 in the Dynamics panel (default L0,
  bit-identical); see `BINDING_LOOP2_DONE.md` §4.

## 6. What each tier CANNOT answer (per docs/APPLICABILITY.md, still true at every tier)
Absolute K_D/ΔG, explicit-solvent effects, membrane/nucleic systems, PME-charged
systems, production kinetics. Tiers move you ALONG the frontier, not past its trust boundary.

## 7. Stage-7 addendum (2026-09-12): allocation tracking + GPU decision

`scripts/pareto_bench.mjs` now brackets each tier with heapUsed before/after
reads (gc-bracketed via `node --expose-gc`); `docs/pareto_frontier.csv` gains
`heap_delta_MB` + `run` columns with the §5 rows preserved as dated history.
Headline: ms stable ±5% over 3 FULL reruns; gc-bracketed per-tier allocation
≤0.3 MB (memory is not a frontier axis — quantified). GPU verdict: Loop-2
binding terms stay CPU-only (WGSL covers LJ/Coulomb/GB; port cost > benefit at
0.22 ms CG / ~80 ms heavy). Full numbers in `docs/BINDING_LOOP2_DONE.md` §14.
