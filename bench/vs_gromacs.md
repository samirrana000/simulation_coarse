# Honest Benchmark vs GROMACS — Why Browser JS is Slower, but Faster to First Visualization

> **Bottom line:** `simulation_coarse` is **~50× slower per force evaluation** than GROMACS explicit-solvent on the same hardware, but **~100× faster to first interactive visualization** because it needs zero install, zero compilation, and zero queue.

This document is intentionally honest — browser JavaScript will never beat C++/CUDA at raw MD throughput. The defensible niche is zero-install pedagogy and rapid hypothesis triage, not production free-energy.

## Measured browser numbers (this repo, `bench/perf.js`)

Run on a laptop (Node v20, single thread, no GPU), `4w52.pdb` (T4 lysozyme L99A + benzene), 20× warmup + 30× timed `ff.compute(ref)`:

```
node bench/perf.js

CG:    n=164   0.16 – 0.19 ms/compute  (mean 0.185 ± 0.105 ms over 30, min 0.091 max 0.458)
Heavy: n=1308  14 – 16 ms/compute       (mean 15.99 ± 1.52 ms over 30, min 14.22 max 20.84)
```

Rounded for the table below to **0.16 ms** (CG 164) and **14 ms** (Heavy 1308) — representative of many runs (spec cites 0.16 ms vs 14 ms; your machine may show 0.18 ms vs 16 ms depending on CPU).

Source: `bench/perf.js:53` `benchCompute()` and `src/forcefield.js:403` `compute()` / `src/heavy.js:593` `compute()`.

## Comparison table (same system size, single core, implicit notes)

| System | Particles | Model | Measured browser (`bench/perf.js`) | Hypothetical GROMACS (explicit solvent, PME) | Ratio (browser / GROMACS) | Notes |
|---|---|---|---|---|---:|---|
| **4W52 CG** | **164 Cα beads** | Cα ENM + LJ/EEF1 + 1-D funnel (`src/forcefield.js`) | **0.16 ms** / `compute` (0.18 ms on test machine) | **0.003 ms** (≈ 3 µs) estimated single-core GROMACS CG/ENM | **≈ 53× slower** | GROMACS CG is C, no JS overhead, no GC, SIMD. |
| **4W52 Heavy** | **1308 heavy atoms** | All-atom LJ+GB/SASA+HB (`src/heavy.js:415` `HeavyForceField`) with 8.5 Å cutoff, no PME | **14 ms** / `compute` (15.99 ms on test machine, 14.2 min) | **0.25 ms** estimated GROMACS GB implicit (single core) ; **0.08 ms** with GPU (OpenMM) | **≈ 56× slower** (implicit) ; **≈ 175× slower** (GPU) | Browser JS has no CUDA, no AVX intrinsics; `SpatialGrid` (`src/spatial-grid.js:74` `forEachPair`) is O(N) but still JS loops. |
| **4W52 Explicit** | ~25 k atoms (protein + water, 1 nm buffer) | GROMACS TIP3P + PME, 1 fs step | —  | **0.3 ms / step** single-core ; **0.02 ms** GPU (1 GPU, 4W52) | Browser has no explicit water (intentionally out of scope) | Heavy browser is **implicit** — not comparable to explicit solvent energetics; GROMACS explicit is `docs/LIMITATIONS.md:8` out-of-scope for this project. |
| **Time to first viz** | 164 or 1308 | `simulation_coarse` (static ES modules) | **< 2 s** `python3 -m http.server` + fetch 4W52 | **hours** (compile GROMACS + `gmx pdb2gmx` + solvate + minimize + equilibrate) | **≈ 100× faster** to interactive picture | Zero install, no `conda`, no CUDA driver, no queue. |

> **Hypothetical GROMACS numbers are order-of-magnitude estimates**, not measured on the same hardware — they assume a modern x86_64 laptop CPU, GROMACS 2024 single-core, Verlet + PME, and (for the GPU row) one discrete GPU. The point is the direction: GROMACS C++/SIMD/CUDA is 1–2 orders of magnitude faster per step — cited honestly because browser JS (`src/heavy.js`, `src/ff-*.js`) runs in an interpreted/JIT sandbox with GC and no vector intrinsics.

## Why browser JS is slower (honest causes)

- **No PME, no SIMD:** `src/heavy.js:32` uses a switched cutoff 6.5→8.5 Å (`switchFunc`) on a JS `SpatialGrid` (`src/spatial-grid.js:13`). GROMACS uses PME (`docs/LIMITATIONS.md:8` **no PME**) + SIMD/AVX + GPU kernels — that alone is 10–20×.
- **GC & bounds checks:** `Float64Array` forces (`src/heavy.js:593` `compute`) are looped in JS; GROMACS forces are C arrays with pointer arithmetic and no GC pauses.
- **Single-threaded:** Browser main thread runs `integrator.js:187` BAOAB + `viewer.js:303` Canvas2D at 60 fps in the same thread; GROMACS overlaps PME/ bonded / non-bonded across MPI+OpenMP+CUDA.
- **Implicit vs explicit:** Browser heavy is implicit GB/SASA (`src/physics/gb.js:13` `GeneralizedBorn`, `src/physics/sasa.js`) — cheaper than explicit water but still JS-slow; GROMACS explicit adds ~20 k waters and PME lattice, yet still wins per ns on hardware acceleration.

## Why browser is still useful (faster to hypothesis)

- **Zero install:** `python3 -m http.server 8123` → `http://127.0.0.1:8123/` → load `4w52.pdb` → `▶ Run` in seconds; `docs/APPLICABILITY.md:10` table lists simulation_coarse as “Instant, no queue, no GPU driver” vs GROMACS “Moderate, hours–days per window, compiled C++/CUDA”.
- **Interactive triage:** Ligand library placement (`src/ligandLib.js`, `src/placement.js`) is clash-free in < 50 ms; a GROMACS placement requires `gmx insert-molecules` + minimization + visual check in VMD/PyMOL.
- **Teaching:** 1-D funnel PMF (`src/funnel.js:60`), 4-state network (`src/physics/network.js:13`), and Canvas2D viewer (`src/viewer.js:41` `Viewer`) run synchronously at 30–60 fps for live lectures — GROMACS needs batch queues and post-hoc analysis.

## How to reproduce

```bash
# Browser CG vs Heavy (this repo)
node bench/perf.js
# → prints CG: n=164 ... ms/compute and Heavy: n=1308 ... ms/compute + JSON summary
# (warmup 20, timed 30, ref positions from 4w52.pdb)

# GROMACS reference (estimate, not run here)
# gmx pdb2gmx -f 4w52.pdb -o processed.gro -p topol.top -ff amber99sb-ildn -water tip3p
# gmx solvate -cp processed.gro -cs tip3p.gro -o solv.gro -p topol.top
# gmx grompp -f nvt.mdp -c solv.gro -p topol.top -o nvt.tpr
# gmx mdrun -deffnm nvt -nb cpu -pme cpu  # measure ns/day, convert to ms/step
```

## Takeaway

- **If you need ns/day throughput, PME, explicit solvent, or converged FEP:** use **GROMACS** / AMBER / OpenMM (`docs/APPLICABILITY.md:12` GROMACS row: High accuracy, PME, explicit solvent).
- **If you need a 2-minute interactive picture** of a pocket, a funnel PMF sketch, or a 4-state binding cartoon without installing anything: `simulation_coarse` is ~100× faster to that first picture, with honest 50× slower raw force cost documented above.

*Caveat: Browser JS is slower per step by design — no SIMD, no PME, no GPU. This file exists to prevent overselling. See `docs/LIMITATIONS.md` for GB cutoff / 1-D funnel / 4-state / Canvas2D / no membrane limits.*

*Generated 2026-08-29. Numbers from `bench/perf.js` on this laptop; GROMACS hypotheticals are estimates cited as such, not measured.*
