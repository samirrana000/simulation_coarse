# Limitations

This document enumerates known model limitations that affect interpretation of results. For an honest assessment of when the model is appropriate, see `docs/APPLICABILITY.md` and the `README.md` Limitations section which both reference this file.

## Explicit honest limitations (bullet list)

- **GB cutoff, no PME:** Electrostatics uses a screened Generalized Born / Coulomb kernel with a smooth switching cutoff 6.5 → 8.5 Å (`src/heavy/params.js:36` `switchFunc`) and **no PME** — long-range interactions beyond 8.5 Å are neglected, so highly charged systems accrue 5 % truncation error and must use GROMACS/AMBER PME instead.
- **1D funnel CV only:** Binding PMF is reconstructed only along the single radial collective variable `r = |COM_lig − COM_pocket|` (`src/funnel.js:165` `cv()`) — orthogonal barriers (rotation, pocket dehydration) are invisible, so the 1-D PMF can hide hidden slow degrees of freedom.
- **4-state kinetics toy:** The Chemical Network Model discretizes binding into exactly 4 macrostates (Bulk, Encounter, Intermediate, Bound) with Kramers rates (`src/physics/network.js:13` `NETWORK_STATES`) — this **4-state** toy illustrates timescale separation but is not a converged MSM and cannot replace PyEMMA/MSMBuilder.
- **Canvas2D vs WebGL:** Rendering is 2-D Canvas2D with painter's depth sorting (`src/viewer.js:481` `order.sort` by `pz`) and no depth buffer — overlapping spheres at high depth complexity produce z-fighting that WebGL depth buffering would eliminate. Whether that costs frame rate on a given integrated GPU is **not measured in this repo** (no canvas in the headless gate), so no fps figure is quoted; the tradeoff is qualitative and the cost/benefit reasoning is in `docs/VIEWER.md`.
- **Opt-in aromatic/halogen weak terms are approximate:** heavy mode *can* add π-stack, cation-π and halogen σ-hole terms (`src/physics/weakint.js`, wired in `src/heavy.js` `_weakInteractions`), but they are **off by default** and must be requested with `par.weak: "on"`. They are ring-centroid / angle-gated Gaussians with fixed per-element ε (`HALOGEN_EPS`), not fitted multipole models: they capture the *direction* of these interactions, not their magnitude in a specific protein. Do not read absolute ΔG contributions from them. CG mode has **no** equivalent (no ring planes at Cα). Chalcogen bonds remain unimplemented.
- **No membrane / nucleic acids / QM:** The force field has no lipid bilayer, no DNA/RNA backbone terms, and no quantum Hamiltonian (`src/ff-params.js`, `src/heavy.js`) — membrane proteins, nucleic-acid folding, or bond-breaking chemistry are out of scope and require GROMACS/CHARMM/NAMD/QM-MM. (`src/physics/solvation/membrane_slab.js` is an IMM1-style *implicit* flat slab — a z-only desolvation/dielectric profile with no lipid atoms — not a bilayer.)

## Long-range electrostatics: cutoff 8.5 Å, no PME

- **Model:** Electrostatics uses a screened Coulomb / Generalized Born kernel with a smooth switching cutoff from **6.5 Å → 8.5 Å** (AMBER-style switch `switchFunc` in `src/heavy/params.js:36`) and **no Particle-Mesh Ewald (PME)**. There is **no PME** — long-range interactions beyond 8.5 Å are neglected.
- **Implication:** Not for highly charged systems (nucleic acids, membranes, poly-electrolytes, high ionic strengths that alter screening). Errors grow with net charge and system size; use GROMACS/AMBER/OpenMM with PME for those.
- **Convergence check:** Increasing the cutoff from 8.5 Å to 12 Å changes the total non-bonded + GB energy by **≈5 %** for the 4W52 test system (164 residues + benzene), consistent with an approximate truncation error. The 8.5 Å choice trades accuracy for browser O(N) grid speed; energy drift is monitored in tests but not converged to production-MD tolerances.
- **Recommendation:** If your system is highly charged, has multivalent ions, or requires <1 % electrostatic convergence, do not use this browser model — use explicit-solvent PME in a production MD engine.

## 1D funnel metadynamics

- **Model:** Well-tempered metadynamics bias is deposited only along the 1-D distance CV described above (`src/funnel.js:200` `addForces`); `bins=96`, `rMax=24 Å`, `sigma=0.3 Å`.
- **Impact:** Conformational states that are orthogonal to this CV appear degenerate (e.g., flipped benzene orientation at same `r`), so free-energy barriers in those dimensions are not sampled.
- **Mitigation:** For production binding pathways use PLUMED 2-D CVs or funnel metadynamics in GROMACS + PLUMED; this browser PMF is illustrative.

## 4-state kinetics model

- **Model:** `ChemicalNetworkModel` (`src/physics/network.js:66` `rebuildRateMatrix`) enforces detailed balance on a 4×4 rate matrix with attempt frequency `ν₀=1e10 s⁻¹` and toy barriers 1–5 kcal/mol.
- **Impact:** The 4-state discretization is a pedagogical illustration — kinetics scale with the chosen barrier heights and `ν₀`, so reported `k_on/k_off/K_D` are order-of-magnitude teaching numbers, not experimental rates.
- **Recommendation:** For publishable kinetics use a full MSM (PyEMMA, MSMBuilder) on explicit-solvent trajectories.

## Rendering: Canvas2D only (no WebGL)

- **Current:** `Viewer` (`src/viewer.js:303` `render`) uses Canvas2D `arc` + manual painter's sort (`src/viewer.js:481`); no GPU depth buffer, no instancing.
- **Impact:** Canvas2D overdraw scales linearly with atom count and fails z-correctness for dense heavy-mode scenes, so overlapping spheres can z-fight. This is a standing trade-off of the shipped renderer, not a gap awaiting a second renderer.
- **Status:** Canvas2D is the only renderer in this project. There is no WebGL renderer — no `WebGLRenderingContext` is created and no GPU library is a dependency. A WebGL path would trade away the zero-dependency / anywhere-runs property for a frame-rate win this demo does not need; see `docs/VIEWER.md` for the full reasoning and the cost/benefit note.

## Out of scope: membranes, nucleic acids, QM

- **Model scope:** Force fields are protein Cα ENM (`src/forcefield.js:91` `ForceField`) or all-atom heavy with LJ/GB/SASA/h-bonds (`src/heavy/forcefield.js:82` `HeavyForceField`) — no lipid force fields, no nucleic-acid sugar-pucker/dihedral corrections, no QM/MM.
- **Impact:** Simulating a membrane protein, DNA folding, or catalytic bond breaking with this model produces qualitatively wrong energetics and should not be attempted.
- **Use instead:** GROMACS/CHARMM-GUI for membranes, AMBER OL3/bcs1 for nucleic acids, CP2K/ORCA + NAMD/OpenMM for QM/MM.

## Other limitations (summary)

- Coarse-grained / heavy-atom approximate force field (LJ + GB + SASA + directional H-bonds) — heuristic parameters, not calibrated for free-energy.
- Implicit solvent only (GB + SASA, Debye-Hückel κ ≈0.329·√I); no explicit water entropy.
- Langevin BAOAB dynamics at coarse time scales — illustrative sampling, not microsecond-converged.

See also `docs/APPLICABILITY.md` for the full trust-boundary comparison table and `docs/VIEWER.md` for Canvas2D vs WebGL rendering details.
