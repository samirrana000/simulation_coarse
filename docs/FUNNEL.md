# Funnel & Well-Tempered Metadynamics — simulation_coarse

> Collective variable: **r = |COM_lig − COM_pocket|** (Å). The funnel is flat inside `r ≤ rFlat` so the bound-state equilibrium is undisturbed; beyond `rFlat` a harmonic wall `U_wall = ½ kf (r−rFlat)²` pulls the ligand toward the pocket.

## 1. Well-Tempered Metadynamics

- Gaussian hills: `g(r') = w·exp(−(r'−r)²/2σ²)` with `σ = sigma` (default 0.3 Å), deposited on a 1-D grid `r ∈ [0, rMax]` (`bins=96`, `Δr = rMax/bins = 0.25 Å`).
- Well-tempered height: `w = w0·exp(−V(r)/(kB·T·(γ−1)))` where `γ = biasFactor` (default 6). `γ=1` recovers plain metadynamics (`w=w0`).
- Deposition schedule: every `hillStride` force calls. **D36** — `hillStride` auto-scales with integrator `dt` to keep ~100 fs physical spacing:
  ```
  scaledStride = hillStride * (0.004/dt)   // dt in ps, 0.004 ps = 4 fs CG baseline
  ```
  At `dt=0.004 ps`, `hillStride=20` → 80 fs; at `dt=0.001 ps` → `scaledStride≈80` → ~80 fs. Call `funnel.getScaledStride(dt)` to compute.

## 2. Bias Grid Resolution (D35)

- Default `Δr = 0.25 Å` just below `σ=0.3 Å`, barely ok. Grid recovers Gaussian integral within 1 %:
  ```
  Σ_k g(r_k) Δr ≈ w √(2π) σ
  ```
  Tested in `tests/test_funnel_grid.js` with `w=0.4, σ=0.3, Δr=0.25`.

## 3. Reweighting: Tiwary-Parrinello c(t) (D32)

- Unbiased free energy: `F(r) = −γ/(γ−1)·V(r) + c(t)`
- `c(t)` from **Tiwary & Parrinello, J. Phys. Chem. B 2015, 119, 736–742, Eq. 8**:
  ```
  c(t) = (1/β) ln [ ∫ ds exp(β·γ/(γ−1)·V(s,t)) / ∫ ds exp(β/(γ−1)·V(s,t)) ]
  ```
  Discretized over grid bins: `sumNum = Σ exp(β·γ/(γ−1)·V_k)`, `sumDen = Σ exp(β/(γ−1)·V_k)`, `c_t = (1/β) ln(sumNum/sumDen)`.
- `Funnel.getPMF()` returns `{r, pmf, dG_vol, c_t}` where `pmf = −γ/(γ−1)·V` shifted so minimum in `[0, 1.5·rFlat]` is 0; `c_t` reported separately. `estimateDG()` uses **pmf** (not raw bias) so the shift and `c(t)` are consistently applied: `ΔG ≈ pmf(rFar)−pmf(rBound) − dG_vol`.

## 4. Volume Correction (D31/D39)

- Restraint volume `V_rest = 4/3 π rFlat³`. Standard-state correction (Boresch et al., J. Phys. Chem. B 2003, Eq. 6):
  ```
  dG_vol = −kT ln(V_rest / V°),  V° = STANDARD_VOLUME = 1660.54 Å³ (1 M)
  ```
  `STANDARD_VOLUME` imported from `src/units.js` in both `src/funnel.js` and `src/physics/network.js`.

## 5. Multiple Walkers (D34)

- Constructor param `nWalkers` (default 1):
  ```js
  const funnel = new Funnel({ nProt, n, ref, nWalkers: 2 });
  // or
  const f1 = new Funnel({ nProt, n, ref, nWalkers: 2 });
  const f2 = new Funnel({ nProt, n, ref, nWalkers: 2 });
  f1.mergeBias(f2); // sums _bias and _biasForce, increments _nHills
  ```
- Each walker deposits hills onto its own grid; periodically call `mergeBias(otherFunnel)` to share bias (summing grids). This is a stub for true shared-memory walkers — callers must synchronize deposition schedules. For production multi-walker, share a single underlying bias array or broadcast hills via messaging.

## 6. PMF Export Provenance (D40)

- `Funnel.exportPMF()` / `Funnel.getPMFcsv()` returns CSV with header:
  ```
  # T=300, gamma=6, hills=123, V0=1660.54
  # dG_vol=..., c_t=...
  r_Ang,pmf_kcal_per_mol
  ...
  ```
- `analysis.js:pmfCsv(funnel)` also emits `# T=..., gamma=..., hills=..., V0=1660.54` as first line.

## 7. Convergence

- HUD grays out `ΔG` until `nHills ≥ 50`; `convergenceSE()` ≈ `kT/√nHills`.
- Rigorous `ΔG°` integrates `r² exp(−βW)` over bound vs unbound (see `integrateDGJacobian()`).

*Last updated: 2026-09-01 — D32/D34/D35/D36/D39/D40.*
