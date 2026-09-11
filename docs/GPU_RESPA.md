# Phase 3 — WebGPU WGSL Kernels & r-RESPA

- `src/compute/wgsl/cell_list.wgsl` — two-pass counting-sort cell list (`bin_count`→CPU prefix-sum→`bin_fill`), cell edge = RCUT 12 Å, workgroup 64.
- `src/compute/wgsl/nonbonded_forces.wgsl` — one-thread-per-atom tiled kernel, 27-cell walk, LJ (Lorentz–Berthelot) + screened
  Coulomb (332.0637) + Still-GB reaction field w/ κ; heavy.js-exact switch; sorted exclusions, 1-4 ×0.5.
- `src/compute/webgpu_backend.js` — `initWebGPU/dispatchNonbonded/isSupported`; unified Float32 4N layout, persistent buffers,
  mapped zero-copy readback; fallback GPU→worker-pool→single-core, never throws.
- `src/physics/integrators/respa.js` — `RESPAStepper.step()`; fast=bonded/angles/torsions/ENM, slow=grid LJ/Coulomb/GB+SASA/LCPO+membrane+funnel;
  B_slow/2 → n×inner-VV (1 fs) → outer exact-OU → refresh → B_slow/2. UI: `respaToggle` + `respaOuter` (2/4 fs, default 2 fs; 4 fs has ~2% shadow bias).

Validation: `node --check` clean; CG+heavy splits bit-exact (max|dF|≤9e-13); 1CRN RESPA vs BAOAB <U> diff 0.130% (<1%);
4HHB N=4558 single-core ~11 steps/s, RESPA ~3.8× fewer slow evals. Headless has no GPU (`isSupported()=false`) — 4–8× GPU is a projection, re-measure in Chrome+GPU.
