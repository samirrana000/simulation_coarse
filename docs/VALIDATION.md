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
