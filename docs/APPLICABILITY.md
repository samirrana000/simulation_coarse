# Applicability & Scope — When to Use simulation_coarse

> **Trust boundary:** `simulation_coarse` is a browser-based, coarse-grained / heavy-atom educational and hypothesis-generation tool. It is **NOT** a replacement for atomistic molecular dynamics (MD) or free-energy methods. Do **NOT** use it for production binding-affinity prediction, free-energy perturbation (FEP), or any regulatory / clinical decision.

---

## 1. At a Glance — Comparative Table

| Tool | Accuracy | Speed (typical system ~160 aa + ligand) | Primary use case | Install / Distribution |
|---|---|---|---|---|
| **simulation_coarse** (this repo) | **Low–Medium** — Cα ENM + implicit LJ/GB/SASA + heuristic H-bonds; no PME, no explicit solvent, no nucleic-acid / membrane physics. Qualitative fluctuations, pocket occupancy, and illustrative PMF only. | **Instant** — 30–60 fps in browser, ~1 ps/frame, whole 4W52 demo runs in seconds on a laptop. No queue, no GPU driver install. | Interactive teaching, rapid pocket exploration, hypothesis generation, lightweight PMF illustration, browser demos. | **Zero install** — static ES modules, `python3 -m http.server`, runs offline after load. |
| **GROMACS** | **High** — all-atom, PME electrostatics, validated force fields (CHARMM, AMBER, OPLS), explicit solvent, replica exchange, FEP/TI. | Moderate — ns/day on CPU, µs/day on GPU cluster; hours–days per free-energy window. | Production MD, rigorous FEP / TI, membrane/nucleic-acid simulations, large-scale sampling. | Compiled C++/CUDA, MPI, manual install or conda/spack; non-trivial build matrix. |
| **AMBER** | **High** — gold-standard all-atom MD, PME, GPU-accelerated (pmemd.cuda), well-tested GAFF/GAFF2 + TIP3P/TIP4P, rigorous alchemical FEP/TI. | Moderate–Fast on GPU — similar to GROMACS for production FEP campaigns. | Production MD and alchemical binding free energies in academic/industrial pipelines. | Licensed (AMBERTools free, pmemd licensed), compiled CUDA build. |
| **OpenMM** | **High** — all-atom, highly extensible Python API, PME, custom forces, supports AMBER/CHARMM force fields, explicit solvent, alchemical. | **Fast (GPU)** — best-in-class single-GPU throughput for all-atom MD; FEP via Python plugins. | Custom MD engines, method development, high-throughput GPU MD, prototyping new FEP protocols. | `conda install openmm`, Python + CUDA/OpenCL, relatively easy. |
| **NAMD** | **High** — all-atom, PME, Charm++ scalable to 10k+ cores, explicit solvent, FEP / TI, QM/MM via ORCA/MOPAC interface. | Scales to supercomputers; single-workstation slower than OpenMM/GROMACS. | Very large systems (ribosomes, membranes, viruses), supercomputer campaigns, QM/MM. | Prebuilt binaries, Charm++ dependency for scaling; config-file driven. |
| **Rosetta** | **Medium** — knowledge-based + physics score functions, excellent for design/docking/folding, not MD; no PME, implicit solvent approximations. | Fast for docking/design; Monte Carlo rather than MD time series. | Protein design, loop modeling, docking (RosettaLigand, GALigandDock), structure prediction. | Large C++ build, PyRosetta pip, steep learning curve. |
| **Mol*** | **N/A (viewer)** — no physics, experimental CIF/PDB rendering, high-quality visualization, density maps. | Instant — pure WebGL viewer. | Structural visualization, density inspection, publication figures, not simulation. | Web component / npm, zero physics. |
| **NGL Viewer** | **N/A (viewer)** — no physics, lightweight WebGL viewer, excellent for embedding in web apps. | Instant — pure WebGL viewer. | Lightweight web visualization, trajectory playback, not simulation. | JS library, CDN or npm, zero physics. |

---

## 2. When to Use simulation_coarse

- You need a **two-minute interactive demo** (e.g., T4 lysozyme L99A + benzene, PDB 4W52) to illustrate ENM fluctuations, ligand pocket occupancy, or a funnel-metadynamics PMF sketch.
- You are **teaching** Langevin BAOAB integration, ENM vs. all-atom, or Kramers binding kinetics and want students to run something without installing GROMACS.
- You want to **triage ligand placement hypotheses** (library placement, clash relaxation) before committing a supercomputer allocation to OpenMM/GROMACS FEP.
- You need an **offline, zero-backend web component** that works on a conference laptop with no internet.

## 3. When NOT to Use simulation_coarse — Trust Boundaries

> **This is NOT for FEP.** Any number labeled `ΔG`, `K_D`, `k_on/k_off`, or PMF produced by simulation_coarse is a qualitative, implicit-solvent, coarse-grained estimate. It has **not** been calibrated against experimental binding data, does not sample explicit water or entropy correctly, and does not converge like a rigorous alchemical free-energy calculation. Do not report it as a binding affinity. For quantitative affinity, use GROMACS/AMBER/OpenMM FEP, TI, or MBAR with explicit solvent and proper error analysis.

Additional boundaries:

- **No quantitative kinetics:** The Chemical Network Model (CNM) rates use Kramers theory on toy barriers (1–5 kcal/mol) and an attempt frequency `ν₀ = 1e10 s⁻¹`. They illustrate timescale separation, not experimental `k_on/k_off`.
- **No force-field transferability claims:** LJ/GB parameters are heuristic, charges are mean-neutralized, and the SASA term is EEF1-lite. Results are system-dependent.
- **No PME / long-range electrostatics:** Electrostatics is screened Coulomb + distance-dependent dielectric, cut at 9 Å. Do not use for highly charged systems.
- **Limited sampling:** A browser `advance(steps, 14 ms)` loop cannot reach microseconds. Metadynamics hills are illustrative, not converged.

## 4. What We Will NOT Do — Explicit Out-of-Scope List

The project will **NOT** add the following. If you need them, use the tools in the table above.

| Scope item we will NOT do | Why not / Use instead |
|---|---|
| **Membrane simulations** (lipid bilayers, membrane proteins, insertion energetics) | Requires explicit lipids, anisotropic pressure coupling, and validated lipid force fields. Use GROMACS / NAMD / CHARMM-GUI + OpenMM. |
| **Nucleic acids** (DNA/RNA folding, base-pairing, sugar puckers) | No nucleic-acid bonded terms, stacking, or backbone dihedral corrections. Use AMBER (OL3/bsc1) or GROMACS. |
| **QM/MM or reactive chemistry** (bond breaking/formation, catalysis, metal redox) | No quantum Hamiltonian, no SCF. Use NAMD/QM/MM, CP2K, or QM package + OpenMM. |
| **PME (Particle-Mesh Ewald) long-range electrostatics** | Explicit PME is O(N log N) and CPU/GPU heavy; browser target cannot host it. Use GROMACS/AMBER/OpenMM with PME. |
| **Rigorous FEP / TI / MBAR binding free energies** | Requires alchemical intermediates, soft-cores, explicit solvent, and converged sampling. Use GROMACS/AMBER/OpenMM + MBAR. |
| **Glycosylation, post-translational mods, cofactors beyond simple LJ** | No specialized PTM libraries. Use Rosetta or all-atom MD with appropriate patches. |
| **Large assemblies (> ~2000 residues) at high accuracy** | O(N) grid still degrades and Cα ENM loses accuracy; no domain decomposition. Use NAMD/GROMACS on HPC. |
| **NMR/X-ray refinement, cryo-EM fitting** | No experimental restraint potentials. Use Rosetta, Phenix, or ISOLDE. |

## 5. How to Cite vs. Production Tools

If you use simulation_coarse for a figure or a hypothesis, cite it as an **interactive coarse-grained model** and clearly state its limitations (list above). Any quantitative claim about a ligand's affinity or residence time must be followed up with at least one of:

- GROMACS / AMBER / OpenMM **explicit-solvent FEP** (with cycle closure and error bars), or
- Experimental measurement (ITC, SPR, MST).

Do not cite simulation_coarse alone as evidence that a compound binds.

---

## 6. FAQ

**Q: Can I get a publishable ΔG from the PMF panel?**
A: No. The PMF is a *well-tempered metadynamics on a 1-D COM distance* with heuristic ligand–protein LJ and no explicit water entropy. It is useful for comparing *relative* pocket accessibility between ligands in the same browser session, not for publishing an absolute ΔG.

**Q: The Chemical Network says k_on = 10⁸ M⁻¹s⁻¹ — is that real?**
A: No. It is a Kramers rate with a fixed attempt frequency and toy barriers, plus a 1 mM standard-state assumption. It exists to connect the 3D trajectory (ps) to network kinetics (ns–s) for teaching.

**Q: What *is* trustworthy here?**
A: Backbone fluctuation patterns (B-factors), qualitative pocket shape, clash-free placement, and the pedagogical link between ENM, Langevin dynamics, and binding-state discretization. All else is exploratory.

---

*Last updated: 2026-08-29. This document is intentionally critical — the goal is to prevent over-interpretation of a lightweight browser model.*
