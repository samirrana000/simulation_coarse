# Phase 1 — Physics Rigor (opt-in)

All off by default; legacy paths unchanged unless flags set.

- `src/physics/forcefield/amber14sb.js` — ff14SB-style bond/angle tables, proper Fourier
  `V=Σ Vn/2·[1+cos(nφ−γ)]` multi-term, harmonic impropers. Getters: `getBondParams/getAngleParams/getDihedralParams/getImproperParams`.
- `src/physics/forcefield/tirion_anm.js` — `γij=γ0·(R0/rij0)^6`, R0=3.81 Å; SS classify via r(i,i+3)+pseudo-φ/θ;
  enable via `par.enmModel="tirion"` / `useTirionNetwork()`. Note: rescale γ0 ~20× vs uniform γ (Kmean≈0.09 vs 2.0).
- `src/physics/solvation/gb_obc2.js` — GB-OBC II (α=1.0,β=0.8,γ=4.85), HCT Born radii, Still `f_GB`, Debye-Hückel
  `κ=0.329·√(I·298.15/T·78.5/εout)`; `computeBornRadii/gbEnergyForces`. Enable via `par.gbModel="obc2"`.
- `src/physics/solvation/lcpo_sasa.js` — LCPO P1–P4 areas + `F=−γ·dA/dr`; `lcpoSasa()`. Enable via `par.sasaModel="lcpo"`.
  Offset differs from legacy SASA — use ΔG only.
- `src/physics/solvation/membrane_slab.js` — IMM1/HDGB slab |z|≤15 Å, `ε(z)`, transfer profile; `par.membrane={on}`.
  Keep zCenter away from soluble proteins.

Validation (1CRN, 10k steps): NVE drift 0.061% (gate 0.05%; 50-step gate passes), NVT RMSD 1.199 Å PASS, 25–33k steps/s.
Finite-difference: GB-OBC2 1.5e-9, LCPO ≤7.5e-12. `test_all` 32/32 PASS.
