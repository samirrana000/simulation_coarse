# Chemical Network Model & Markov Kinetics — simulation_coarse

> **4-state toy, not full MSM; use PyEMMA for production.** This module is a pedagogical 4-state toy (Bulk / Encounter / Intermediate / Native Bound) with Kramers rates and detailed balance, intended to illustrate timescale separation and reactive flux—not to replace a converged MSM built with PyEMMA/MSMBuilder on explicit-solvent trajectories.

## 1. States (src/physics/network.js:13)

| ID | Name | ΔG (kcal/mol) | Description |
|----|------|---------------|-------------|
| 0 | Bulk Solvated | 0.0 | Free 3D diffusion in bulk aqueous solvent |
| 1 | Encounter Complex | −1.8 | Non-specific electrostatic surface association |
| 2 | Binding Intermediate | −3.4 | Vestibule / pocket entrance channel |
| 3 | Native Bound Pose | −6.2 | Specific pocket lock & hydrogen bonded pose |

The energies scale with `setBoundEnergy(dG)` and obey `G1≈0.28·G3, G2≈0.55·G3` so intermediate wells track the bound well after `rebuildRateMatrix()`.

## 2. Detailed Balance & Rate Matrix (src/physics/network.js:104)

Construction enforces microscopic reversibility:

```
pi_i · K_ij = pi_j · K_ji
pi_i ∝ exp(−G_i / kBT)   with bulk chemical potential G0 + kBT·ln([L]/1M), V°=STANDARD_VOLUME=1660.54 Å³
K_ij = ν0·exp(−(dG_barrier − G_i)/kBT)   forward; reverse from pi ratio
K_ii = −Σ_{j≠i} K_ij
```

- Attempt frequency `ν0 = 1e10 s⁻¹`.
- Barriers `2.5 / 3.8 / 4.5 kcal/mol` (Bulk→Encounter diffusion, Encounter→Intermediate desolvation, Intermediate→Native steric gating) are heuristic; real barriers should be derived from PMF via `setBarriersFromPMF(pmf,r)` as `peak − well`.
- Bulk association is pseudo-first-order: `k01 = k_on_base·[L]` with `[L]=concentrationM` (default 1 mM); `k10` follows from detailed balance `pi0/pi1`.
- Checked in `tests/test_ck.js` and `tests/test_all.js:148` (`|pi0·K01 − pi1·K10|<1e-8`).

## 3. Master Equation Propagation (src/physics/network.js:222)

```
dp/dt = p·K
```

`propagateMasterEquation(totalTimeSec, maxSubSteps)` uses **small dt** via CFL condition `dt ≤ 0.4 / max|K_ii|` and adaptive substepping `numSubSteps = ceil(totalTime/dtMax)` (min 50, capped by `maxSubSteps`). Explicit Euler with renormalization:

```
dp[j] += Σ_i p[i]·K[i][j] · dt
p ← max(0,p) / Σ p
```

Chapman-Kolmogorov is validated in `tests/test_ck.js`: `‖exp(K·2t) − exp(K·t)²‖ < 1e-6` (matrix exponent) and `propagate(2t) ≈ propagate(t)∘propagate(t)` within `1e-6`. Stability is validated in `tests/test_master_stability.js` (500 steps, `Σ p =1±1e-10`, no negative).

## 4. Gillespie SSA (src/physics/network.js:171)

`stepGillespie()` draws:

```
τ = −ln(r1)/k_tot ,  k_tot = −K_ii
next = sample j ∝ K_ij
```

Dwell-time distribution is exponential with mean `⟨τ⟩ = 1/k_tot`. Validated in `tests/test_gillespie.js` (1000 steps, per-state `⟨τ⟩` within 10% of `1/k_tot`, histogram printed).

## 5. Pose Classification (src/physics/network.js:262)

```js
classifyPose(comDist, nContacts, ligRMSD)
 3 if comDist≤5.5 && nContacts≥12 && ligRMSD≤2.5  → Native Bound
 2 else if comDist≤9.0 && nContacts≥4              → Intermediate
 1 else if comDist≤15.0 || nContacts≥1            → Encounter
 0 else                                           → Bulk
```

Validated in `tests/test_classify.js` (200-frame CG Langevin traj via `ForceField`+`Langevin`, histogram, F1 for Native or at least no-crash + histogram). Live-tracking linkage `isLiveTrackingActive` is exercised in `tests/test_live_tracking.js` and in `src/main.js:648` where `networkModel.classifyPose(cv, nc, ligRMSD)` drives `networkModel.currentState` when `isLiveTrackingActive` is true.

## 6. Transition Path Theory (src/physics/network.js:278)

Solves ` (K·q)_i =0` for intermediates with boundaries `q0=0, q3=1`, iterative solver (100 sweeps). Reactive flux:

```
J_ij = max(0, pi_i·K_ij·(q_j − q_i))
totalReactiveFlux = Σ_j J_0j
```

Returned as `{committors, flux, totalReactiveFlux}`.

## 7. Kinetics Observables (src/physics/network.js:318) — Units

`computeKinetics()` returns:

| Field | Symbol | Unit | Formula |
|-------|--------|------|---------|
| `k_on` | k_on | **M⁻¹ s⁻¹** | `totalReactiveFlux / (pi0·[L])` (apparent second-order); fallback `|K01|/[L]` |
| `k_off` | k_off | **s⁻¹** | `k_on · K_D(M)` with `K_D(M)=exp(dG/kBT)` |
| `KD_uM` | K_D | **µM** | `exp(dG/kBT)·1e6`; `K_D(M)=exp(dG/kBT)` |
| `dG` | ΔG_bind | **kcal/mol** | `G3 − G0` (standard state) |
| `mfpt_ns` | MFPT | **ns** | `Σ T_{i→i+1}` with `T_{i→i+1}=(1+ k_{i,i-1}T_{i-1→i})/k_{i,i+1}` |
| `pi` | π | — | Boltzmann populations |
| `committors` | q⁺ | — | from TPT |

`STANDARD_VOLUME = 1660.54 Å³` provides `V°` provenance (imported from `src/units.js:32`, checked in `tests/test_provenance.js` and `tests/test_pmf_export.js`). Values are **order-of-magnitude teaching numbers**, not experimental rates—barriers and `ν0` are toy parameters.

## 8. Scope & Limitations

- **4-state toy, not full MSM; use PyEMMA for production.** For publishable kinetics build a full MSM on explicit-solvent trajectories with PyEMMA (https://github.com/markovmodel/PyEMMA) or MSMBuilder, with proper featurization, TICA, and validation (Chapman-Kolmogorov, implied timescales).
- Barriers are heuristic unless set via PMF.
- Rates scale with `ν0` and `concentrationM`; changing either rescales `k_on/k_off` linearly.
- No membrane / nucleic-acid / QM terms; not for production FEP.

## 9. Tests

```
node tests/test_ck.js                 # Chapman-Kolmogorov CK: propagate(2t)≈propagate(t)² <1e-6
node tests/test_classify.js           # 200-frame CG traj + classify histogram, F1>0.6 or no-crash
node tests/test_gillespie.js          # 1000 Gillespie steps, dwell ⟨τ⟩≈1/k_tot ±10%
node tests/test_master_stability.js   # 500 propagations, Σ p=1±1e-10, no negative
node tests/test_live_tracking.js      # isLiveTrackingActive mock classification vs true COM/contacts
```

*Last updated: 2026-09-01 — E42/E44/E46/E47/E48/E49/E50.*
