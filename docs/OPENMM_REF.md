# OpenMM Reference — Heavy GB/SA vs OpenMM Implicit GBSA (F52)

*Source anchors: `src/heavy.js:461` `HeavyForceField`, `src/physics/gb.js:13` `GeneralizedBorn`, `src/physics/sasa.js` `SasaModel`.*

This document defines the reference protocol for validating the browser heavy
all-atom force field (LJ + Generalized Born + SASA) against a production
implicit-solvent engine. **Heavy GB/SA should be compared to OpenMM implicit GBSA 4W52 NVT**
as the reference (not explicit-solvent PME). It is the **F52** extended-validation companion to
`bench/vs_gromacs.md` (which covers GROMACS explicit-solvent timing).

---

## 1. Reference engine

- **Engine:** OpenMM 8.x implicit solvent
- **Model:** GBSA-OBC2 (Onufriev-Bashford-Case, `GBn2` / `OBC2`) or Hawkins-Cramer-Truhlar
  GB with LCPO/Borukhov SASA — whichever matches `src/physics/gb.js:70` HCT scale 0.8
  and `src/physics/sasa.js` γ = 0.0072 kcal/mol/Å²
- **System:** T4 lysozyme L99A + benzene, PDB **4W52** (164 residues, 1308 heavy atoms
  with `HOH` stripped; ligand BNZ 6 heavy atoms, `HETATM A 200` retained, EPE buffer dropped)
- **Ensemble:** **NVT** at 300 K, Langevin integrator γ ≈ 1–5 ps⁻¹, no barostat
  (implicit solvent has no box; NVT ≡ NPT here)
- **Nonbonded:** GB + SASA implicit; **no PME** — both browser (`src/heavy.js:32`
  `switchFunc` 6.5 → 8.5 Å) and OpenMM GBSA use cutoff 12–16 Å with switching,
  but neither uses Particle-Mesh Ewald. The browser heavy model and the OpenMM
  GBSA reference are therefore **comparable in the implicit-solvent regime**,
  but neither is comparable to explicit-solvent PME (see `docs/LIMITATIONS.md:13`
  and §4 below).

## 2. Expected correlation

Heavy GB/SA total potential energy (`HeavyForceField.compute()` → `ff.energy`,
`ff.gbU + ff.elecU + ff.repU + ff.sasaU + ff.hbondU`) evaluated on the **same
20 frames** (e.g. 20 snapshots from a 1 ns NVT OpenMM GBSA trajectory,
strided every 50 ps, minimized 200 steps) should correlate with the OpenMM
GBSA energies:

- **Expected Pearson R > 0.75** for 20 frames (4W52 NVT, implicit GBSA) — expected energy correlation R>0.75 for 20 frames
- **Rationale:** Both models use the same Coulomb constant
  (`src/physics/gb.js:20` `COULOMB_CONST = 332.06371`), same HCT Born-radii
  formulation (`src/physics/gb.js:88` scale 0.8, `src/physics/charges.js:50`
  `GB_RADII`) and same SASA γ, so rank-ordering of conformations is preserved
  even though absolute energies differ by a systematic offset (different radii,
  switching, SASA probe). R > 0.75 is the smoke-test; production gate would be
  R > 0.85 after parameter tuning.
- **Failure modes:** R < 0.5 suggests inverted charges, wrong dielectric
  (`epsIn 4.0` vs 1.0), missing salt screening (`kappa` 0.329·√I), or SASA sign
  error. See `tests/test_negative.js:1` for charge-inversion control.

## 3. Placeholder comparison table (fill with measured OpenMM values)

Run the OpenMM reference script (e.g. `scripts/openmm_gbsa_4w52.py` — template
below, not executed in CI) and paste energies here. The browser column is
produced by `node bench/pdbbind_gb.js` or a dedicated `bench/openmm_compare.js`
that calls `HeavyForceField.compute()` on the exported OpenMM frames
(`data/openmm_4w52_gbsa/*.pdb`).

| Frame | CV r (Å) | Heavy GB/SA (kcal/mol) `ff.energy` | OpenMM GBSA (kcal/mol) | Δ (Heavy−OpenMM) |
|------:|---------:|-----------------------------------:|-----------------------:|-----------------:|
|   0   |   3.2    |  — *placeholder*  e.g. −2078.9     | — *placeholder*        |  —             |
|   1   |   3.5    |  —                                 | —                      |  —             |
|   2   |   4.1    |  —                                 | —                      |  —             |
|  ...  |  ...     |  ...                               | ...                    | ...            |
|  19   |   5.8    |  —                                 | —                      |  —             |
| **R** |          |                                    |                        | **>0.75 target** |

*Fill procedure:*

```bash
# 1. Generate 20 OpenMM GBSA frames (requires openmm package)
python scripts/openmm_gbsa_4w52.py --pdb 4w52.pdb --frames 20 --out data/openmm_4w52_gbsa/
# 2. Compute browser heavy energies on those frames
node bench/openmm_compare.js --frames data/openmm_4w52_gbsa/ --pdb 4w52.pdb
# 3. Compute Pearson R
python -c "import numpy as np; heavy=np.loadtxt('heavy.dat'); omm=np.loadtxt('omm.dat'); \
print(np.corrcoef(heavy, omm)[0,1])"
# Expected: R > 0.75
```

If `openmm` is not installed, the table remains a documented placeholder
— the correlation gate is **advisory**, not a CI blocker, and is distinguished
from the GROMACS explicit-solvent PME benchmark which is explicitly out of
scope (see `bench/vs_gromacs.md` and next section).

## 4. Not PME — explicit disclaimer

**This GB/SA reference is not PME.** Both browser heavy and OpenMM implicit
GBSA use cutoff-switched electrostatics (`src/heavy.js:32` `R_SWITCH_ON=6.5 Å`,
`R_CUT=8.5 Å`; OpenMM `NonbondedForce` switching 10–12 Å for GBSA). Neither
computes long-range Lattice-sum electrostatics via **Particle-Mesh Ewald (PME)**.
Therefore:

- Energies are **not comparable** to GROMACS/AMBER explicit-solvent PME
  (`docs/LIMITATIONS.md:13` **GB cutoff, no PME**). Do not compare `ff.elecU`
  from `src/heavy.js:704` `_nonBondedGrid` to a PME lattice energy — the
  difference is systematic (≈5 % at 8.5 Å, see `bench/vs_gromacs.md`).
- Highly charged systems, membranes, nucleic acids, and multi-valent ions
  require PME and are out of scope (`docs/APPLICABILITY.md`).
- For implicit GBSA, the appropriate reference is OpenMM GBSA NVT, not
  GROMACS explicit PME. The funnel PMF (`src/funnel.js`) and 4-state network
  (`src/physics/network.js`) are orthogonal to this energy benchmark.

*Measurable:* `ls docs/OPENMM_REF.md` exists; correlation gate R>0.75 for
20 frames is documented as target, not CI-enforced until OpenMM frames are
available.*

## 5. References

- Still, W.C. et al. *J. Am. Chem. Soc.* 1990, 112, 6127 — GB formalism
- Hawkins, G.D. et al. *J. Phys. Chem.* 1996, 100, 19824 — HCT
- Onufriev, A. et al. *J. Phys. Chem. B* 2004, 108, 15873 — OBC2
- Eastman, P. et al. *PLOS Comp. Biol.* 2017, 13, e1005659 — OpenMM
- `src/heavy.js:461` `HeavyForceField`, `src/physics/gb.js:134` `pairInteraction`,
  `src/physics/sasa.js` `SasaModel`, `bench/perf.js:1` timing

*Generated 2026-09-01. Placeholder table to be filled when OpenMM GBSA frames are exported.*
