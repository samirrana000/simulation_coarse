# Validation — methods & numbers (FP6)

Trust boundary, with numbers: **ranking-only use; NOT FEP, NOT absolute Kd.**
CG ΔH replays ±0.3 seeded but single-rep ΔG_est underbinds by 2–3 kcal/mol
against experiment with bootstrap ±0.16 ⊕ replica-SD ±7.1 (entropy-driven) —
absolute ΔG is sampling-noise-limited. Cα pocket-entropy sign is
model-dependent; CG ala-scan has no hotspot resolution (|ΔΔG| < 0.05 ≈ noise);
the solvent term's sign convention is code-exact UNfavorable (+γ·ΔSASA) with a
physics-sign review open (§17). Absolute-ΔG claims additionally require:
BNZ-only (buffer-free) legs, ≥10 frames/DOF Schlitter, nonzero measured ΔSASA.

## Methods table

| System | Protocol (seeds/steps) | Observable | Computed ± error | Experimental / literature | Verdict |
|---|---|---|---|---|---|
| 4W52 BNZ+EPE (record path) | CG 2000 steps stride 2, T300 ζ8.0, charges+dir-HB, seeds 101/1101 | ΔH (LJ/Coul/HB/desolv) | −6.82 ± 0.16 (−1.88/+1.04/−0.13/−5.85) | — (record replay ±0.30) | REPLAYS record |
| 4W52 BNZ+EPE | same, 3 reps 101/202/303 | ΔH / −TΔS_pocket / ΔG_est | ΔH [−6.82,−7.50,−6.00] mean −6.77 SD 0.75; −TΔS [+4.58,−3.71,+10.42] mean +3.76 SD 7.1 | ITC −5.20±0.20 / NMR −4.20±0.10 (Mondal 2018 Tab.1) | RANKING-ONLY (underbinds 2–3, SD dominates) |
| 4W52 BNZ-only (true cavity) | CG same, seeds 101/1101, 14-res pocket | ΔH | −3.64 ± 0.08 (−0.98/0/0/−2.66); desolv 73.1% | same as above | TRUE-CAVITY anchor; EPE inflates \|ΔH\| by −3.18 (reported, not corrected) |
| 4W52 BNZ-only + real SASA | + LCPO cross-burial stride 10, 100+100 evals | ΔSASA → −TΔS_solv → ΔG_est | 167.1 ± 3.5 Å² → +2.01 ± 0.04 → ΔG −0.16 (±50% band [−1.17,+0.84]) | benzene-burial scale | NONZERO solvent term; sign review open |
| 4W52 EPE (flexible lig) | CG 500 steps stride 2, seeds 501/1501 | ΔS_lig / ΔH / pocket / ΔSASA | ΔS_lig 0.00305 NONZERO (4 rotatable); ΔH −2.58±0.13; pocket 7; ΔSASA 429.6±11.2 | BNZ rigid control ΔS_lig exactly 0 | rotbond path LIVE end-to-end |
| 1CRN apo-vs-apo (null) | CG 500 steps stride 2, seeds 701/1701, COM-8Å pocket (16 res) | ΔH / −TΔS_pocket | 0.00 ± 0.00 (\|ΔH\|<0.5 vs 4W52 −3.6/−6.8); −TΔS +5.67 bounded <20 | null: no ligand ⇒ no signal | NULL PASSES (noise bounded, never forced zero) |
| 4W52 heavy FULL | 3 reps × 300+1500 steps stride 2, seeds 1001/2002/3003 | ΔH / −TΔS full | ΔH −19.73 (LJ −8.2/desolv −11.4); −TΔS [22.63,−24.27,13.96] mean +4.11 | same ITC/NMR | ΔH trustworthy; sign NOT-RESTORED (variance-dominated) |
| 4W52 heavy SMOKE | 1 rep × 50+300 steps stride 2 | same (pipeline smoke) | 13/13 in 12.4 s; rep0 −TΔS −20.97, ΔH −17.02 | — | SMOKE parity (same 13 asserts) |
| 4W52 heavy LONG | 10 ps × 3 (5000 f/leg, f/DOF 15.7) + 26 χ torsions | −TΔS full/sc/χ | full +6.76 SD 14.06; sc +1.85; χ +0.24 ± 0.33 (≈0, ≤0.5 bound) | — | NO converged restriction signal separable from rattle |
| 4W52 ala-scan (CG) | BNZ-only, rCut 8Å/maxN 12, relax 80 | ΔΔG top-3 | TYR88 +0.014 / MET102 +0.012 / ILE100 +0.009; all `~noise` (< 0.05 floor) | T4L liners (Merski 2015): 1/3 nominal (MET102) | NO CG hotspot resolution (encoded in code) |
| Schlitter unit | 2-DOF Gaussian, 50k samples | sampled-vs-exact err | 0.20% (< 1% bar) | analytic | UNIT ANCHOR |
| Torsion unit | locked vs 12-bin uniform rotor | S | locked ≈0; uniform = kB·ln12 ± 0.001 | kB·ln12 exact | UNIT ANCHOR |
| Bench tiers (4W52 FULL) | pareto_bench, gc-bracketed | ms/step | L0 0.237 / L1 0.235 / L1+BindLog 0.216 / L2 84.8 / L3 83.1 / L2-full 85.5; L0:L2 ≈ 400:1 | — | Loop-2 ≈1.0× on heavy; GPU port declined (CPU ACCEPT) |

## ΔG audit trail — every ΔG-like number, and what it assumed

The methods table above reports one observable per row and nothing else. This
section is the audit trail behind it: **every** ΔG-like number this build
produces, so a reader can expand a ΔG into the assumptions behind it instead of
taking it on trust. Each entry carries one of three nature tags — the same three
words the in-app "ΔG audit chain" subpanel prints, so a doc reader and a UI
reader see one vocabulary:

- `measured` — produced by sampling this system's own trajectory
- `assumed` — a scale, a reference state, or a leg-definition choice that no
  sampling in this tool can put an error bar on
- `out-of-scope` — physics this tool does not have, named so a reader is not
  left inferring it is missing by accident

The table below is generated from `src/dg-estimators.js` (DG_ESTIMATORS) and
re-checked against the live code paths on every FAST-tier run by
`tests/test_dg_drift.js`, so it is never typed twice. `mixed` in the "measured?"
column means the number is arithmetic over both measured and assumed terms:
read the terms, not the headline.

**No number in this section revises a number above.** This adds provenance to
results that are already reported; it does not revise them.

### Every ΔG-like estimator in this build

| id | what the user sees | formula | code path | measured? | reported uncertainty | on the display path? |
|---|---|---|---|---|---|---|
| `funnel_dg_hud` | live HUD ΔG (WTM funnel bias) | `pmf(rFar) − pmf(cv0) − dG_vol,  dG_vol = −kT·ln[(4/3)π·rFlat³ / STANDARD_VOLUME]` | `src/funnel.js estimateDG()` | mixed | hill-count — convergenceSE() = kT/√nHills, a decorrelation-time heuristic, NOT an error on the answer | yes |
| `funnel_dg_jacobian` | r²-Jacobian ΔG° quadrature (computed, never shown) | `−kT·ln[∫_{r≤rFlat} r²e^{−βW}dr / ∫_{r>rFlat} r²e^{−βW}dr] − dG_vol` | `src/funnel.js integrateDGJacobian()` | measured | NONE — the quadrature returns a single number with no dispersion and no replicate legs; a track spread is not an uncertainty on a converged free energy | no (computed, never displayed) |
| `thermo_dg` | ΔG estimate (ΔH − T·ΔS decomposition) | `ΔG = ΔH − T·(ΔS_pocket + ΔS_lig + ΔS_solv); ΔH = ⟨U_bind⟩_holo, ΔS_pocket = Schlitter(holo) − Schlitter(apo), ΔS_solv = −ΔSASA·0.012/T` | `src/analysis/thermodynamics.js computeThermodynamics().dG_estimate` | mixed | bootstrap+SD — 20-block bootstrap on ΔH; NOTHING on the −TΔS leg, whose replica SD this document measures at 7.1 kcal/mol | yes |
| `ala_ddg` | per-residue ΔΔG (alanine scan, thermodynamic cycle) | `ΔΔG = [E(mut·holo) − E(WT·holo)] − [E(mut·apo) − E(WT·apo)], steepest-descent relaxed legs` | `src/analysis/alanine_scanning.js scanResidue()` | measured | NONE — the scan computes no replicate legs, so there is no spread to call an error bar; a noise FLOOR (\|ΔΔG\| < 0.05 kcal/mol) is the only guard and it is a threshold, not an uncertainty | yes |
| `jarzynski_df` | Jarzynski ΔF (SMD pulling ensemble) | `ΔF = W_min − (1/β)·ln⟨e^{−β(W−W_min)}⟩ over 4 constant-velocity pulls` | `src/analysis/unbinding_smd.js jarzynskiFreeEnergy()` | mixed | bootstrap — 200-resample nonparametric bootstrap over the 4 pull works | yes |
| `koff_surrogate` | k_off ranking surrogate (SMD rupture) | `mean_i β·(W_rupt,i − ΔF_rupt) in kT, β = 1/kT — a unitless ranking score` | `src/analysis/unbinding_smd.js koffSurrogate()` | assumed | NONE — a ranking score over 4 pulls has no ensemble to disperse; the code computes none and an error bar here would be invented precision | yes |
| `network_dg` | 4-state network ΔG_bind (→ K_D) | `ΔG_bind = E(bound state) − E(solvent state); K_D = exp(ΔG_bind/kT) — the state energies are ASSIGNED (bound 1.00·ΔG, 0.55·, 0.28·), not sampled` | `src/physics/network.js computeKinetics()` | assumed | NONE — there is no sampling here at all — the state energies are set by setBoundEnergy() and the ratios are literals in the constructor, so an error bar would be a fabricated measurement | yes |
| `eef1_lig_dg` | EEF1-lite per-atom desolvation ΔG_a (a potential, not an estimator) | `U_desolv = Σ_a ΔG_a·B_a with B_a = 1 − exp(−n_a/3), ΔG_a = −0.30 (C) … −0.50 (I) kcal/mol, a per-element TABLE CONSTANT` | `src/ff-binding.js (U_desolv) / src/physics/params.js (the ΔG column)` | assumed | NONE — a fixed parameter has no uncertainty; the error of the model this table encodes is unbounded by anything this tool measures, and the ±50% SASA band on the solvent leg (this document) is the only error quantification anywhere near it | yes |

Five of the eight carry no error bar at all, and each one says why in its own
cell. That is this project's never-write-zero rule: an estimator with no
dispersion reports `NONE` plus the reason, never `0`. A downstream consumer that
sees `NONE` must not substitute a zero.

### Finding F1 — the live HUD ΔG is in no row of the methods table

`funnel_dg_hud` is on the display path: a user sees it in the HUD after any run
long enough to deposit hills. It is **absent from every methods-table row above**,
and that absence is asserted rather than left to inference —
`tests/test_dg_drift.js` recomputes the HUD number from `src/funnel.js` on every
run and turns red if this document ever starts quoting it without also
declaring the funnel protocol (seeds, step count, hill budget, leg definition)
in its own Protocol cell. Until such a row exists, this audit trail is where
that number lives, and the methods table above stays silent about it.

### The two ΔG stories are irreconcilable with what this tool implements

Stated plainly, because a reader should not have to work it out: **the ΔG a user
sees in the HUD and the ΔG reported in the methods table above are
irreconcilable with what this tool implements.** They are not two measurements
of one quantity — they share no term. The HUD number is a vertical free-energy
difference along ONE radial CV read two points from a bias grid that is still
filling; the table's number is a horizontal ΔH − T·ΔS difference between a holo
leg and an *internal relaxation of the same coordinates* as its apo leg.

Closing the gap is **not a documentation fix**. It would require physics this
tool does not have, named term by term below, and none of those terms is
scheduled work. The gap is therefore reported, never averaged away.

| bridge id | nature | reason | measured move | what is missing |
|---|---|---|---|---|
| `shared-observable` | out-of-scope | the two numbers are not two measurements of one quantity — they share no term | — (a different observable, so no single-pair number describes it) | a common-leg estimator: the same trajectory pair feeding both the radial PMF and the ΔH/ΔS decomposition |
| `hill-convergence` | out-of-scope | the HUD number has not converged, and its printed error bar does not say so | measured: +5.35 kcal/mol of monotone drift while the reported SE falls by ~100× | a block-averaged or replica-resampled error bar on the reconstructed PMF |
| `jacobian-and-ct` | out-of-scope | two standard-state terms are computed and then not applied (the radial Jacobian, and the Tiwary–Parrinello c(t) offset) | measured: +2.11 kcal/mol | applying either term to the displayed number — a physics change, explicitly out of scope |
| `no-delta-s-decomposition` | out-of-scope | the HUD number contains no ΔH, no ΔS and no ΔSASA of its own | — | an entropy decomposition along the same CV |
| `system-definition` | assumed | the two stories do not use the same system definition (one carries the crystallisation buffer EPE, the other does not) | measured: 3.18 kcal/mol of ΔH on the same seed pair | one declared system definition per reported number, carried in the record |
| `leg-and-replicas` | out-of-scope | one story has two legs and three replicas, the other has one trajectory from one seed | measured: 7.1 kcal/mol replica SD on −TΔS_pocket, against a ±0.16 kcal/mol ΔH bootstrap SE | a ≥3-replica funnel run reporting the between-replica SD of the reconstructed PMF |
| `physics-level-is-a-no-op-here` | assumed | the L1-vs-L0 assumption moves this system's number by exactly zero, because benzene carries no charge and no H-bond capability | measured: 0.00 kcal/mol on BNZ-only; +1.27 kcal/mol on the BNZ+EPE record path | — (nothing missing; an assumption whose cost depends on the ligand, which is why it cannot be read off one number) |
| `mislabeled-leg` | out-of-scope | FINDING (label only, nothing changed): one surface names the wrong leg distance | — (label only; no number changed) | — (a reporting-surface correction, deliberately not made inside a results run) |

The measured moves dominate the attribution: `leg-and-replicas` (7.1 kcal/mol)
and `hill-convergence` (5.35 kcal/mol) are each an order of magnitude larger
than the disagreement between the two headline numbers above, which is why
reconciling them is a sampling-and-protocol question rather than an arithmetic
one. `movedBy` is measured, not asserted — each entry is a pair of runs on
4W52 + benzene, seed 101, app defaults, CG Cα, T=300, ζ=8, m=110, whose only
difference is the choice named in that row. The live gap between whichever two
estimators a user is actually looking at is recomputed by the in-app ΔG audit
chain and shown beside this list; it is reported, never averaged into a single
number.

## How to run

- Fast (gate/CI): `node tests/test_all.js` → 352/352 (~16 s).
- Slow-smoke: `node tests/test_all.js --slow` → 404/404 (~40 s; heavy SMOKE).
- Slow-full: `node tests/test_all.js --slow-full` → 404/404 (~230 s; heavy FULL).
- Alone: `node scripts/validate_1crn_null.mjs` (8/8, <1 s) ·
  `node scripts/calibration_4w52.mjs` (14) · `node scripts/validate_flexlig.mjs` (10).

## References

Protocols: `scripts/test_thermo.mjs` (CG) · `scripts/test_thermo_heavy.mjs`
(heavy SMOKE/FULL) · `scripts/calibration_4w52.mjs` (anchor) ·
`scripts/validate_flexlig.mjs` (EPE) · `scripts/validate_1crn_null.mjs` (null).
Record: `docs/BINDING_LOOP2_DONE.md` §§8–21, §27 (FP6).
