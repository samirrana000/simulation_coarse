# Performance — HPC Wiring (G62–G70)

This document tracks HPC-related performance wiring. All targets marked aspirational are not yet benchmarked and must be validated before enabling by default.

## Worker Pool (G62)

- Flag: `settingsState.backend==="workers" && n>800` routes to `workerPool.computeParallel` (`src/main.js:39`, `src/worker-pool.js:91`).
- Guard: `_busy` prevents concurrent `computeParallel` calls (`src/worker-pool.js:35`).
- Aspirational speedup ≥1.5× on 4 cores — not yet benchmarked. Validate via `bench/worker_speedup.js` (prints "estimated 1.5× on 4 cores").
- See `bench/worker_speedup.js:8` for simulated estimate.

## GPU WebGPU Compute (G64)

- Clamp: `sr<5` present via `Math.min(sigma/r,5)` / `min(sigma/r,5.0)` in WGSL (`src/gpu.js:51` WGSL, `src/gpu.js:214` JS helper `gpuClampSr`).
- Test: `tests/test_gpu_clamp.js` verifies `sr=200→5`.
- Aspirational correlation: GPU vs CPU R>0.999 — not yet benchmarked, placeholder target (`src/gpu.js:6`).

## Neighbor List Skin (G66)

- **Verlet skin 2Å, rebuild every 10 steps, 20% cut** — aspirational target, not yet implemented.
- Current: `SpatialGrid` rebuilds every `compute()` call with `R_CUT=8.5Å` and no skin (`src/heavy.js:473`, `src/spatial-grid.js:8`).
- Desired: introduce 2 Å skin (effective cutoff 10.5 Å), rebuild every 10 steps, yielding ~20% reduction in pair-list rebuild cost. Until implemented, the O(N) grid still beats O(N²) but incurs per-step hash cost.

## dt Auto-Tuning Honest (G67)

- Ligand-aware dt in `src/integrator.js:143` (`_pickDt`): CG alone 4fs vs CG+ligand 1.7fs vs heavy 1fs.
- `maxDt` = 0.004 (CG alone) | 0.0017 (CG+ligand) | 0.001 (heavy).
- Test: `tests/test_dt.js` prints those values.

## Adaptive Steps per Frame (G68)

- `advance(maxMs=14)` caps wall-clock per animation frame to 14 ms (`src/main.js:525` `state.integ.advance(steps,14)`, `src/integrator.js:262` `advance(stepsWanted, maxMs=12)` default, caller passes 14).
- Keeps 60 fps UI responsive under heavy load.

## Memory Leak Guards (G69)

- `src/recorder.js:31` `maxFrames` cap (default 500) documented; `maybeCapture` auto-stops when `frames.length >= maxFrames` (`src/recorder.js:80`).
- Test: `tests/test_recorder_cap.js` verifies recorder caps at `maxFrames`.

## Honest vs GROMACS (G70)

- See `bench/vs_gromacs.md` — honest ~50× slower per force evaluation than GROMACS, but ~100× faster to first visualization. No further wiring needed.

## Performance Budget (J99)

Budgets are enforced in `bench/budget.json:1` and checked in CI (see `.github/workflows/check.yml:19` `CI budget warn` comment — warns, does not fail).

```json
{"heavy_compute_ms": 2.0, "fps": 30, "cg_compute_ms": 0.5}
```

- `bench/budget.json:1` `heavy_compute_ms: 2.0` — p95 `HeavyForceField.compute(pos)` on 4W52 heavy (1308 atoms) must be ≤2.0 ms. Measured via `bench/perf.js:108` `Heavy: n=1308 ... ms/compute`. Exceed ⇒ `::warning` in CI, not error.
- `bench/budget.json:1` `cg_compute_ms: 0.5` — p95 `ForceField.compute(pos)` on CG 164 beads must be ≤0.5 ms (`bench/perf.js:107` `CG: n=164`).
- `bench/budget.json:1` `fps: 30` — Interactive threshold. `Viewer.render` + `advance(maxMs=14)` must sustain ≥30 fps on 4W52 heavy integrated GPU. See `docs/VIEWER.md:27` and `bench/perf.js` scale note. Below 30 fps ⇒ CI warns.
- **Budget source:** `bench/budget.json:1` is the single source of truth. `docs/PERFORMANCE.md:47` budget section and `manuscript/reproduce.sh:49` `budget check` both import from it. CI step is `CI budget warn` (soft gate) — see `.github/workflows/check.yml:19`.

To reproduce locally:

```bash
node bench/perf.js | tee /tmp/perf.log
python3 -c "import json; b=json.load(open('bench/budget.json')); print(b)"
# CI warn logic is also in manuscript/reproduce.sh:50 budget check (warn not fail)
```

If a change pushes `heavy_compute_ms` above 2.0 or drops fps below 30, the PR CI will emit `::warning :: Heavy compute X ms exceeds budget 2.0 ms` — fix by optimizing `src/heavy.js` spatial grid or reducing `R_CUT` overhead before merge.

---
*Teams: all aspirational figures must be validated via `bench/perf.js` and `bench/alloc.js` before claim.*
