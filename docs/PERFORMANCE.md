# Performance — HPC Wiring (G62–G70)

This document tracks HPC-related performance wiring. All targets marked aspirational are not yet benchmarked and must be validated before enabling by default.

## Worker Pool (G62)

- Flag: `settingsState.backend==="workers" && n>800` routes to `workerPool.computeParallel` (`src/controllers/accelerate.js:37`, `src/worker-pool.js:91`).
- Guard: `_busy` prevents concurrent `computeParallel` calls (`src/worker-pool.js:35`).
- Aspirational speedup ≥1.5× on 4 cores — not yet benchmarked. Validate via `bench/worker_speedup.js` (prints "estimated 1.5× on 4 cores").
- See `bench/worker_speedup.js:8` for simulated estimate.

## GPU WebGPU Compute (G64)

- Clamp: `sr<5` present via `Math.min(sigma/r,5)` / `min(sigma/r,5.0)` in WGSL (`src/gpu.js:68` WGSL, `src/gpu.js:229` JS helper `gpuClampSr`).
- Test: `tests/test_gpu_clamp.js` verifies `sr=200→5`.
- Aspirational correlation: GPU vs CPU R>0.999 — not yet benchmarked, placeholder target (`src/gpu.js:6`).

## Neighbor List Skin (G66)

- **Verlet skin 2Å, rebuild every 10 steps, 20% cut** — aspirational target, not yet implemented.
- Current: `SpatialGrid` rebuilds every `compute()` call with `R_CUT=8.5Å` and no skin (`src/heavy.js:473`, `src/spatial-grid.js:9`).
- Desired: introduce 2 Å skin (effective cutoff 10.5 Å), rebuild every 10 steps, yielding ~20% reduction in pair-list rebuild cost. Until implemented, the O(N) grid still beats O(N²) but incurs per-step hash cost.

## dt Auto-Tuning Honest (G67)

- Ligand-aware dt in `src/integrator.js:143` (`_pickDt`): CG alone 4fs vs CG+ligand 1.7fs vs heavy 1fs.
- `maxDt` = 0.004 (CG alone) | 0.0017 (CG+ligand) | 0.001 (heavy).
- Test: `tests/test_dt.js` prints those values.

## Adaptive Steps per Frame (G68)

- `advance(maxMs=14)` caps wall-clock per animation frame to 14 ms (`src/controllers/tick.js:91` `state.integ.advance(steps,14)`, `src/integrator.js:262` `advance(stepsWanted, maxMs=12)` default, caller passes 14).
- Keeps the UI responsive under heavy load. The 14 ms cap is the mechanism; the resulting frame rate is **not measured here** (see the budget section).

## Memory Leak Guards (G69)

- `src/recorder.js:31` `maxFrames` cap (default 500) documented; `maybeCapture` auto-stops when `frames.length >= maxFrames` (`src/recorder.js:80`).
- Test: `tests/test_recorder_cap.js` verifies recorder caps at `maxFrames`.

## Honest vs GROMACS (G70)

- See `bench/vs_gromacs.md` — honest ~50× slower per force evaluation than GROMACS, but ~100× faster to first visualization. No further wiring needed.

## Performance Budget (J99)

`bench/budget.json` is the single source of truth, evaluated by exactly one
implementation — `bench/budget_check.js` — which both
`.github/workflows/check.yml` and `manuscript/reproduce.sh` call. (They used to
carry private copies of the same python; they had already drifted, and both
printed an fps budget for a quantity neither measured.)

Every key is either **measured** — with the bench file and metric that produce
it — or **explicitly not measured**, with the reason. `tests/test_budget_coverage.js`
fails if a new key arrives with neither.

```json
{
  "heavy_compute_ms": 16.0,
  "cg_compute_ms": 0.5,
  "fps": null
}
```

- **`heavy_compute_ms: 16.0`** — measured. `bench/perf.js` `heavy.meanMs`:
  `HeavyForceField.compute(ref)` on 4W52 heavy (1308 atoms), 20 warmup + 30 timed.
  **Measured 14.11 ms** (range 13.86–14.23 over 8 independent runs; worst
  single-run max 15.39). The 16.0 is that value plus ~14 % headroom, so it is a
  **regression tripwire, not a target**.

- **`cg_compute_ms: 0.5`** — measured. `bench/perf.js` `cg.meanMs`:
  `ForceField.compute(ref)` on the 164 Cα beads. **Measured 0.174 ms**
  (range 0.164–0.195). Unchanged from birth; it was reachable then and is now.

- **`fps: null`** — **NOT MEASURED. No frame rate is claimed anywhere.** A frame
  rate needs a frame: a canvas, a compositor, a display refresh, a GPU. This repo
  has zero dependencies by design, there is no canvas polyfill and no jsdom, and
  `new Viewer({})` in bare Node throws `TypeError: canvas.getContext is not a
  function` — `src/viewer.js` binds `window.addEventListener("resize")` and reads
  `devicePixelRatio` at construction. `requestAnimationFrame` does not exist in
  Node either. Until a real-browser leg in `tests/manual/` writes a
  frames/second figure to a file CI reads, this stays `null` and no script prints
  a number for it. That is an honest absence, and it is better than the fake it
  replaces.

  What a headless bench *can* say about interactivity, and does: heavy compute is
  **14.1 ms**, and `advance(steps, 14)` caps physics at 14 ms of wall clock per
  animation frame (`src/controllers/tick.js:91`). **14.1 > 14**, so heavy mode fits at most
  one force evaluation per frame at the reference position. That is a real,
  measured, falsifiable statement about frame pacing. It is not a frame rate.

  **Known remaining frame-rate strings in the tick loop (`src/controllers/tick.js`)** — the module map, the
  G68 comment above the `advance` call, the H76 HUD-debounce comment, and the
  phase note on the top bar all still name a nominal refresh rate. They are
  design notes to the next maintainer, not published claims, and `src/` was
  outside the write scope of the change that produced this section, so they were
  left alone rather than silently deleted or silently ignored. They are recorded
  here so the next person to touch `src/controllers/tick.js` can settle them.
  `tests/test_budget_coverage.js` scans the publishing surfaces (docs, bench,
  .github, manuscript, scripts, root markdown) and fails if any of those grows a
  new frame-rate number; widening it to `src/` is a one-line change to
  `SCAN_DIRS`.

### Where 2.0 ms came from (it was never real)

The budget file's entire history is **one commit**, `3ee4914` (2026-09-11), and it
*creates* the file already containing `heavy_compute_ms: 2.0`. There is no earlier
value for the code to have regressed from. Running that same commit's own
`bench/perf.js` today gives **15.17 ms**; running its `src/heavy.js` directly
gives **51.2 ms** at `de10a73` (2026-08-11), falling to 14.1 ms now. So 2.0 was
fiction at birth: the O(N) `SpatialGrid` work in `3ee4914` took heavy from ~51 ms
to ~15 ms, and the budget was written against neither number. It has never once
been met, and CI emitted `::warning Heavy compute ... exceeds budget 2.0` on
**every run since** — a permanently lit warning that gated nothing.

Reproduce any of it:

```bash
node bench/perf.js | tee /tmp/perf.log    # the measurement
node bench/budget_check.js                # the contract check (exit 1 if over)
git log -p --follow bench/budget.json     # the whole history: one commit
```

If a change pushes `heavy_compute_ms` past 16.0, `bench/budget_check.js` reports
`OVER` and CI emits a `::warning`. If `fps` is ever added back as a number, a
browser leg must produce that number first — `tests/test_budget_coverage.js` will
reject the key otherwise.

---
*Teams: all aspirational figures must be validated via `bench/perf.js` and `bench/alloc.js` before claim.*
