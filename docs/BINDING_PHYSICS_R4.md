# Phase R4 — Entropy & Enthalpy Decomposition (ΔG = ΔH − TΔS)

Loop-1, step 4. Prototype: `/tmp/opencode/r4_entropoproto.mjs` (zero-dep, verified).

## 1. Method survey (browser-engine-applicable)

| Method | Formula | Frames needed | Cost | Accuracy expectation | Reference |
|---|---|---|---|---|---|
| **Schlitter QH** | S ≤ (kB/2)·ln det(I + kBT·e²·σ̄/kB²) | ≥10× DOF | O(f·d² + d³), d=3N | upper bound; good for ΔS(holo−apo) differences | Schlitter 1993 |
| Andricioaei–Karplus QH | exact 1D harmonic sum from σ̄ eigenvalues | ≥10× DOF | same + eigensolve | tighter, mode-count sensitive | Andricioaei & Karplus 2001 |
| **Torsion Shannon (ligand)** | S = −kB·Σ p_i ln p_i over discretized rotamers (30° bins) | ~100 frames/rotor | O(f·n_tors) | converges FAST; ±0.1-0.3 kcal/mol | standard rotamer analysis |
| Rigid-body / standard-state | ΔG°_trans+rot ≈ 2.4 kcal/mol @1M (fixed analytic) | — | — | constant, fold into PMF offset | Gilson et al. 1997 |
| Solvent entropy proxy | −TΔS_solv ≈ −ΔSASA·(0.005–0.02 kcal/mol/Å²) | 1 frame pair | O(LCPO) | ±30% (scale uncertainty) | hydrophobic-release scales (Dill/Schneider lineage) |
| Dihedral PCA / MIST | mutual-information over torsions | 1000+ | O(f·n_tors²) | better than QH for flexible ligands | Numata et al. 2012 |

MM/GBSA practice (Genheden & Ryde 2015): QH/NMA on ~1000 snapshots, errors ±1-3 kcal/mol
absolute, better for relative; enthalpy from pairwise-decomposed ⟨U⟩.

## 2. Engine decomposition pipeline (design)

**Inputs** (all already recordable): holo + apo trajectories (positions per frame);
per-frame pairwise energy components (U_LJ, U_Coul, U_HB, U_desolv, U_new-terms);
pocket residue set (within 8 Å of ligand COM in native frame).

**ΔH** = ⟨U_bind⟩_holo − ⟨U_bind⟩_apo. Apo U_bind ≡ 0 by construction in the engine
(binding terms defined protein–ligand). Pairwise components give the interaction
breakdown (vdW / electrostatic / H-bond / desolvation / π-X / cation-π / halogen —
each term is separately accumulated in the R2/R3 kernels already).

**ΔS_lig** = torsion-Shannon over ligand rotatable bonds (bound trajectory) − same
over the bulk/funnel-bulk region of the same trajectory (ligand in solvent far from
protein; the funnel-PMF runs already visit bulk). Benzene (rigid): ~0.

**ΔS_pocket** = Schlitter(holo pocket, 3N DOF) − Schlitter(apo pocket). Kabsch-align
frames to native before covariance (removes drift). Requires ≥10 frames/DOF: for a
14-residue pocket (42 DOF) ≥420 frames — our 2 ps stride recorder default (500
frames) just meets it.

**ΔS_solv** = −ΔSASA(holo−apo buried) × scale. LCPO SASA already exists (Phase 1);
use 0.012 kcal/mol/Å² (mid-scale) with stated ±50% band.

**Consistency check**: ΔG_PMF (funnel metadynamics, panel 5) vs ΔH − TΔS_tot. The
metadynamics ΔG implicitly contains ALL entropy; agreement within combined error
bars validates both. Persistent disagreement = sampling error, report both.

## 3. Prototype verification

### (a) Schlitter vs exact covariance (2 Gaussian DOFs, 50k samples)
- exact 0.0680, sampled 0.0681 kcal/mol/K → **err 0.20%** ✓
### (b) Torsion Shannon vs analytic (3 rotors: uniform/biased/locked)
- analytic 0.00378, sampled 0.00377 kcal/mol/K → **err 0.30%** ✓
### (c) 4W52 CG pipeline (holo 6000 steps vs apo 6000 steps, ENM + native-contact binding)
- Unit-system audit en route caught TWO prototype bugs (velocity conversion 1 kcal/mol
  = 418.98 Da·Å²/ps² backward in two places) — after fix: σ_v = 1.498 Å/ps (theory 1.51 ✓),
  per-axis pocket σ 0.23-0.25 Å (ENM expectation ✓), apo pocket more mobile than holo ✓.
- **ΔH = ⟨U_bind⟩ = 3.26 kcal/mol** (equipartition cross-check 5.96 — same order;
  difference = ligand partially escaping the harmonic well + protein side absorbing
  momentum, honest ±).
- **ΔS_pocket = S_holo − S_apo = −0.011 kcal/mol/K → penalty +TΔS = 3.21 kcal/mol**
  (literature pocket-restriction range 2–12 ✓).
- Benzene rigid → ΔS_lig ≈ 0 (would be computed from funnel-bulk torsions for flexible
  ligands).
- Consistency: ΔG ≈ 3.26 − 3.21 ≈ 0.05 (native-contact model — nearly athermal by
  construction; the real engine's LJ/desolvation terms shift ΔH substantially; that
  is the Loop-2 integration job).

## 4. Per-frame data requirements (for R6 capture-port design)

| Quantity | Bytes/frame (164 Cα + 21 lig) | Notes |
|---|---|---|
| positions (Cα+ligand) | 3×4×N ≈ 2.2 KB Float32 | already recorded (recorder.js) |
| pairwise U components (7 terms) | 7×4 = 28 B | needs per-term accumulators in ff-binding/R2/R3 kernels (Loop 2) |
| ligand torsion bins | 1-2 B/torsion | derived from positions, no extra storage |
| pocket SASA | 4 B | from LCPO each frame (or stride-10) |
| Born-radii snapshot | optional | heavy mode only, stride-100 |

## 5. Loop-2 integration status (S5 DONE)

`src/analysis/thermodynamics.js` implemented (Schlitter via Cholesky ln-det, signed
torsion Shannon, Horn-quaternion Kabsch alignment, block bootstrap ΔH SE, ΔSASA proxy)
+ `Thermodynamics (ΔH/ΔS)` panel section (`thermoBtn`, analysis-panel.js handler:
holo = recorded frames, apo = internal 2000-step binding-off relaxation).

**Measured (4W52 CG, charges+directional-HB on, 2000 steps × 3 replicas):**
- ΔH = −6.9 ± 0.16 kcal/mol: LJ −1.9 / Coul +0.5 / HB −0.1 / desolv −5.7
- −TΔS_pocket per-run [−3.42, −1.98, −1.48] — **negative and systematic**.
- **Honest physics finding:** at Cα-only resolution the pocket-entropy sign is
  model-dependent: the ligand is a 63-DOF thermal bath that *injects* pocket-Cα
  mobility through the binding terms, while true pocket restriction lives in
  sidechain rotors (invisible at Cα). The trustworthy observables at this tier are
  ΔH + component split, ligand torsion ΔS, and bounded converged |ΔS_pocket|.
  Heavy mode (sidechain atoms present) should restore the physical sign — S7 bench.
- test_thermo: 7/7 (Schlitter 0.20%, torsion kB·ln12 exact, ΔH finite, bounded ΔS).

### Stage-1 verdict (2026-09-12): heavy mode does NOT restore the physical sign — NOT-RESTORED

Protocol (`scripts/test_thermo_heavy.mjs`, 13/13 asserts on finiteness/bounds — the
sign itself is never asserted): 4W52 heavy, BNZ-only ligand (EPE buffer excluded via
`heteroSelection {BNZ:true, EPE:false}` — same selection as `tests/test_negative.js`;
EPE would corrupt the ligand-COM pocket with a second diffusing body). Holo =
`HeavyForceField` weak:on + trackTerms (charges inherent); apo = ligand-free clone
(same protein atom indices; `bindingU ≡ 0` verified — the heavy analogue of CG
`binding:{on:false}`, cf. `apoSystemOf()` heavy branch). `LangevinIntegrator` reads
per-atom `ff.masses` (dt auto-tunes to 1 fs); 300-step equil discard + 1500 production
steps, stride 2 → 750 frames/leg × 3 replicas. Pocket = 106 heavy protein atoms within
8 Å of BNZ COM (56 backbone N/CA/C/O + 50 sidechain, incl. 14 Cα); per-atom masses in
Schlitter via the new additive `masses` option in `computeThermodynamics` (uniform
fallback untouched; `meta.massModel` records which). ΔH from the 4 exact
protein↔ligand cross terms (lj/coul/hb/desolv); S3 weak totals (pi/cpi/xb) are
whole-system and reported as context only, not folded into ΔH.

Per-replica −TΔS_pocket (kcal/mol), three independent runs of the committed protocol:
- run A: [8.23, 2.84, 3.63] → mean +4.90
- run B: [10.12, −3.62, 24.26] → mean +10.26
- run C: [−11.67, −13.67, 21.24] → mean −1.37
- run D (post-verdict confirmation): [−2.78, 59.93, −7.73] → mean +16.47

Pooled 12 replicas span −13.7…+59.9 (one apo leg wandered far in 1.5 ps): at 750 frames / 318 DOF (frames/DOF 2.4) the
Schlitter ΔS is variance-dominated — the sign is undetermined, and the early positive
means (runs A/B) were sampling noise, not restoration.

Deep-sampling pilot (single pair, 500 equil + 6000 production steps stride 1 → 6000
frames/leg, /tmp pilot): block curve full-pocket −TΔS at 750/1500/3000/6000 frames =
−36.20/−0.95/+9.64/−4.78; at 6000 frames (frames/DOF 18.9): full −4.78, backbone
−5.82, sidechain +1.04, Cα −1.26. The Cα value reproduces the CG scale (−2.5); the net
sign at converged sampling is NEGATIVE.

Diagnosis: (1) the ligand-bath effect persists in heavy dynamics — backbone −5.82 /
Cα −1.26 at depth, i.e. binding terms pump ligand thermal motion into pocket backbone
exactly as in CG (no ENM involved: heavy has none; the covalent/angle network is still
rattled); (2) the hypothesized sidechain-restriction mechanism is REAL but small
(sc +1.04 at depth, positive-leaning in 2/3 committed runs) — overwhelmed by the
backbone bath term at an 8 Å pocket / ps timescale; (3) sampling: 750-frame Schlitter
over 318 DOF is noise (±15 spread); frames/DOF ≥ 10 is necessary but not sufficient —
stride-2 1 fs frames decorrelate on ~100 fs vibration timescales (≈15 independent
samples per leg), so multi-ps production is required.

What WOULD restore an honest positive sign: (a) 10+ ps production × ≥3 replicas to
converge the sidechain-restriction term; (b) an explicit pocket-χ torsion-Shannon term
(heavy analogue of ΔS_lig — converges in ~100 frames/rotor per the §1 table);
(c) sidechain-only pocket entropy (excluding the rattled backbone shell). NOT sidechain
restraints (would fake the sign).

Trustworthy observables unchanged, plus one addition: heavy ΔH = −19.7 ± 0.3 kcal/mol
(LJ −8.2 / desolv −11.4, Coul/HB ≈ 0 — physical for nonpolar benzene), tight across all
12 replicas. UI: no thermo-handler change — step-4 wiring was YES-conditional, and a heavy
UI leg would block the browser for minutes; a caption-only note was considered and
rejected as premature without a converged number.

## 5 (original). Loop-2 integration order

1. **DONE (Loop-2 S4)** — Per-term energy accumulators in binding kernels:
   `src/ff-binding.js` (CG: LJ/Coul/HB + desolv, flat fields `bindLJU/bindCoulU/
   bindHBU` + `bindU` vector, guarded by `ff.trackTerms`, default OFF =
   bit-identical) and `src/heavy.js` (heavy: grid kernels fill per-pair
   `bindTerms {lj, coul, hb}` for protein↔ligand pairs; SASA cross-burial
   `bindSasaE` → desolv; S3 weak terms pi/cpi/xb join the vector; split sums
   exactly to `bindingU` on HCT, OBC2 and LCPO paths). BindLog wiring in
   `src/main.js` tick (frames at recorder stride + 7-term energy events +
   contact form/break diff at 5.5 Å + funnel hills via `Funnel.onHill`),
   UI checkbox `bindlogOn` (default OFF). Tests: `scripts/test_bindlog_
   integration.mjs` 14/14 — split sums exact, OFF bit-identical, 300-step run
   → 30 frames / 210 energy events / contact flicker / blob 0-drift /
   ⟨U_LJ⟩ = −1.6 kcal/mol at native; CG ms/step ratio OFF = 1.03× (gate 1.1×).
2. Recorder: optional per-frame component vector + pocket residue list (R6
   format) — covered by the BindLog energy channel for S4; pocket list S5.
3. Analysis panel "Thermodynamics" section (S5): ΔH table + Schlitter
   ΔS_pocket (apo vs holo recorded pair) + torsion ΔS_lig + SASA entropy;
   explicit error bars (frames/DOF ratio, bootstrap over frames).
4. Funnel-PMF cross-check row (S5): ΔG_PMF vs ΔH−TΔS_tot, discrepancy flagged.
5. Standard-state offset 2.4 kcal/mol constant applied to PMF ΔG (1 M reference).

Non-goals: explicit-water entropy, full 3N QH for the whole protein (browser cost;
pocket-only is the useful signal), mode-coupling corrections (MIST) unless Loop 3.

## References

- Schlitter J. 1993, Mol. Simul. (QH upper-bound form).
- Andricioaei I., Karplus M. 2001, J. Chem. Phys. (AK QH).
- Genheden S., Ryde U. 2015, Expert Opin. Drug Discov. (MM/GBSA practice/errors).
- Numata J., Wan M., Knapp E.-W. 2012, Genome Inform. (MIST torsion entropy).
- Gilson M. et al. 1997, Biophys. J. (standard-state binding).
- Chang C.-E., Chen W., Gilson M. 2007 (SASA-scaled solvent entropy calibration).
