# Documentation index — route by question, not by filename

This repo ships **37 markdown files (5 102 lines) and 2 non-document assets
(326 lines)** under `docs/`. Before this file existed, there was no map: a
reader arriving from the top-level `README.md` could not tell which of those
files answered *“can I trust this ΔG number?”*.

So this index is organised by the **question you actually have**, not by
filename. Pick a row, open one file, done. Every file under `docs/` appears
exactly once below; that property is enforced by
`tests/test_docs_index.js` (FAST tier), so this table cannot silently rot.

If you only read one thing: read **Q0**. The honest out-of-scope contract is
the most load-bearing document in the repository and it now has a link here.

---

## Q0 — “Can I trust this number?” — the trust boundary · read first

Nothing else in this repo is more important than these four. The shortest
complete answer to *"is this ΔG/K_D/k_on a real result?"* is: **no, not as an
absolute number** — it is a ranking-only, implicit-solvent, coarse-grained
estimate.

| Document | The question it actually answers | Read |
|---|---|---|
| [VALIDATION.md](VALIDATION.md) | *What did you actually measure, with what seeds, and what was the error bar?* The methods table with per-leg ΔH / −TΔS / ΔG numbers, replica SDs, and ITC/NMR comparison. Ends with the verdict: **ranking-only, NOT FEP, NOT absolute Kd**. | 43 lines · 3 min |
| [LIMITATIONS.md](LIMITATIONS.md) | *What breaks, and how far?* The 8.5 Å cutoff with no PME (≈5 % residual), the 1-D funnel CV, the 4-state kinetics toy, Canvas2D z-fighting, and the no-membrane/no-nucleic-acid/no-QM wall. Each with the tool to use instead. | 50 lines · 3 min |
| [APPLICABILITY.md](APPLICABILITY.md) | *Should I be using this tool at all?* A deliberately critical comparison table against **GROMACS, AMBER, OpenMM, NAMD, Rosetta, Mol\*, NGL** (accuracy / speed / use case / install), plus an explicit “What we will NOT do” table and a FAQ on whether the PMF ΔG is publishable. | 79 lines · 5 min |
| [../ROADMAP.md](../ROADMAP.md) — *repo root* | *What is permanently out of scope?* **“Browser v1 will NOT do: No QM/MM, no explicit membrane, no PME.”** §1 is the hard contract every PR is checked against, and §5 is the rule for proposing a scope change. | 93 lines · 5 min |

> The out-of-scope contract in [ROADMAP.md](../ROADMAP.md) §1 is deliberately
> duplicated (with the same wording) as the trust-boundary banner in
> APPLICABILITY and as the Limitations summary, so that no single unlinked file
> can hide it.

---

## Q1 — “What does the physics actually model?”

| Document | The question it actually answers | Read |
|---|---|---|
| [CG_HEAVY.md](CG_HEAVY.md) | *What are the two force-field resolutions, and how do they hand off without double-counting?* The Cα Tirion ENM and the all-atom heavy model side by side, plus the `nativeContacts` double-coverage guard. | 163 lines · 8 min |
| [PHYSICS_RIGOR.md](PHYSICS_RIGOR.md) | *What does “opt-in rigor” actually switch on?* ff14SB bond/angle/Fourier-dihedral tables, Tirion distance-weighted ANM, GB-OBC II + HCT Born radii, LCPO SASA, implicit membrane slab. All default-OFF. | 17 lines · 2 min |
| [CHEMISTRY.md](CHEMISTRY.md) | *What chemistry/topology/protonation is available?* PROPKA-style pKa→tautomer assignment, GAFF2-lite typing + Gasteiger/BCC-lite charges, chirality/planarity checks, metal-coordination polyhedra. All default-OFF. | 15 lines · 2 min |
| [CHARGES.md](CHARGES.md) | *Where do the partial charges come from, and how accurate are they?* An **approximate united-atom** ff14SB mapping with a per-residue net-charge audit table and an explicit “approximate, not full ff14SB” statement. | 47 lines · 3 min |
| [EXCLUSIONS.md](EXCLUSIONS.md) | *Which non-bonded pairs are excluded or scaled, and why?* The 1-2 / 1-3 / 1-4 policy (1-4 scaled 0.5) plus full intra-ligand exclusion, audited against the `pairKey` code. | 42 lines · 3 min |
| [UNITS.md](UNITS.md) | *What are the units and where are they defined?* Å / ps / kcal/mol / Da / K as the single source of truth in `src/units.js`, plus derived conversion constants. | 88 lines · 4 min |
| [FUNNEL.md](FUNNEL.md) | *How is the 1-D PMF actually reconstructed?* Well-tempered metadynamics hills on `r = |COM_lig − COM_pocket|`, grid quadrature, stride auto-scaling with `dt`, and the Jacobian correction. | 69 lines · 4 min |
| [NETWORK.md](NETWORK.md) | *Can I believe the k_on/k_off?* (No.) The 4-state Bulk→Encounter→Intermediate→Bound model, detailed-balance rate-matrix construction, and the explicit “toy, not a converged MSM — use PyEMMA” warning. | 114 lines · 6 min |

---

## Q2 — “How do I run it, and what does the UI do?”

| Document | The question it actually answers | Read |
|---|---|---|
| [TUTORIAL.md](TUTORIAL.md) | *Walk me through the 2-minute 4W52 demo — and show me where it fails.* Half the document is a worked success, half is the honest failure mode (charged-ligand ΔΔG ≈ ±2–4 kcal/mol). | 142 lines · 7 min |
| [VIEWER.md](VIEWER.md) | *What renders the atoms, and is there a WebGL path?* **WebGL is NOT implemented and is out of scope** — Canvas2D painter's sort, the z-fighting trade-off, and the cost/benefit note for *not* adding a second renderer. | 106 lines · 6 min |
| [PLACEMENT.md](PLACEMENT.md) | *How do I drop a ligand on the protein clash-free?* The four placement entry points, rigid clash relaxation, the hetero-inclusive collision set, and the pocket highlight. | 54 lines · 3 min |
| [EXPORT.md](EXPORT.md) | *How do I get my trajectory out?* XYZ / PDB multi-MODEL file layout, the `maxFrames` memory guard, provenance `REMARK` headers, and the recommended DCD conversion. | 93 lines · 5 min |
| [ACCESSIBILITY.md](ACCESSIBILITY.md) | *Is it usable without a mouse?* The measurable keyboard contract (ESC closes modal + restores focus, Space toggles run), logical tab order, and the axe-core audit with 0 violations. | 83 lines · 4 min |
| [UI_COCKPIT.md](UI_COCKPIT.md) | *What is the “anti-slop cockpit” layout?* One visible Load→Build→Run→Record flow, 5 collapsed panels, the top status bar, the bottom dock, the three-state canvas rule, and the keyboard map. | 16 lines · 2 min |

---

## Q3 — “What can I actually do with it?” (workflows)

| Document | The question it actually answers | Read |
|---|---|---|
| [WORKFLOWS.md](WORKFLOWS.md) | *What analysis modules are wired into panel 6?* Alanine-scanning ΔΔG, DCCM + heatmap, cryptic-pocket volume tracking, and SMD pulling with a Jarzynski ΔF (with its dissipative-protocol caveat). | 69 lines · 4 min |
| [SCORER.md](SCORER.md) | *Can it learn a pose score?* The dependency-free JSON MLP the browser actually ships (6 features → `pose` in the HUD), weights format, and the **ONNX is NOT implemented** note. | 109 lines · 6 min |

---

## Q4 — “How fast is it, and will it run on my machine?”

| Document | The question it actually answers | Read |
|---|---|---|
| [PERFORMANCE.md](PERFORMANCE.md) | *What is the speed story and what is merely aspirational?* Worker-pool wiring, WebGPU clamp + placeholder targets, neighbour-list skin — with every un-benchmarked target labelled as such, and the CI budget that gates regressions. | 68 lines · 4 min |
| [GPU_RESPA.md](GPU_RESPA.md) | *Is there a GPU path and a multi-timestep integrator?* The WGSL cell-list + tiled nonbonded kernels, the GPU→worker→CPU fallback chain, and r-RESPA (1 fs inner / 2–4 fs outer). | 12 lines · 2 min |
| [DOCKER.md](DOCKER.md) | *Can I reproduce the environment in a container?* A `python:3.11-slim` image that mirrors the zero-dependency quick-start, with any port mapping. | 32 lines · 2 min |

---

## Q5 — “How do I cite this, and what is the trust boundary I must state?”

| Document | The question it actually answers | Read |
|---|---|---|
| [CITATION.md](CITATION.md) | *What exactly do I put in the methods section?* The BibTeX for tag `v1.0-jpcb`, the Zenodo DOI placeholder, and the required disclosure that any quantitative affinity claim needs FEP or experiment behind it. | 54 lines · 3 min |

The trust-boundary table itself is in **Q0** (APPLICABILITY) — start there, not here.

---

## Q6 — “How do I build on this, or hand work to a real engine?”

| Document | The question it actually answers | Read |
|---|---|---|
| [PLUMED.md](PLUMED.md) | *Can I use my own collective variables?* How the browser CV `cv = |COM_lig − COM_pocket|` maps onto PLUMED `DISTANCE`, the COLVAR export, and a minimal umbrella-input snippet for GROMACS/OpenMM follow-up. | 137 lines · 7 min |
| [BRIDGE.md](BRIDGE.md) | *How do I load an exported trajectory into a real analysis stack?* The MDAnalysis/MDTraj path to DCD, and how to reproduce the browser's RMSF→B-factor analysis externally. | 135 lines · 7 min |
| [OPENMM_REF.md](OPENMM_REF.md) | *What is the reference protocol for validating the heavy GB/SA kernel?* The OpenMM 8.x implicit-GBSA (OBC2) comparison — the correct reference is implicit solvent, **not** explicit-solvent PME. | 116 lines · 6 min |
| [MMCIF.md](MMCIF.md) | *Can I load mmCIF?* (**Not really.**) PDB is primary; `src/mmcif.js` is a detection-only stub kept for forward compatibility. | 67 lines · 3 min |

---

## Q7 — “What is the history of the binding-physics decisions?”

The nine documents below are **historical design and review records, not
current instructions.** They record how the binding physics was designed,
independently audited, and then actually landed. Read the **final** state of
the series, not the intermediate proposals.

**If you read nothing else here, read R7 and then the LOOP2 close-out.**

| Document | The question it actually answers | Read |
|---|---|---|
| [BINDING_PHYSICS_REVIEW_LOOP1.md](BINDING_PHYSICS_REVIEW_LOOP1.md) | *Did an independent reviewer find anything wrong with R1–R7?* The adversarial Loop-1 audit of the whole program, with its findings and dispositions. | 243 lines · 12 min |
| [BINDING_PHYSICS_R1.md](BINDING_PHYSICS_R1.md) | *What should we model, at what energy scale, in which resolution?* The literature survey that grounded the whole program to the code, with a verified “what the code does today” baseline. | 342 lines · 17 min |
| [BINDING_PHYSICS_R2.md](BINDING_PHYSICS_R2.md) | *The CG upgrade proposal* — virtual sites + the term set — and its explicit CG non-goals (no ring planes at Cα). | 419 lines · 20 min |
| [BINDING_PHYSICS_R3.md](BINDING_PHYSICS_R3.md) | *The heavy-mode weak-interaction proposal* — π-stacking, cation-π, halogen σ-hole. All terms heavy-mode only. | 156 lines · 8 min |
| [BINDING_PHYSICS_R4.md](BINDING_PHYSICS_R4.md) | *How do you get ΔG = ΔH − TΔS out of a browser trajectory?* Entropy/enthalpy decomposition, Schlitter and torsion methods, and the noise floor that limits the result. | 183 lines · 9 min |
| [BINDING_PHYSICS_R5.md](BINDING_PHYSICS_R5.md) | *Where is the accuracy/speed Pareto frontier?* The tier table that motivated the L0/L1/L2 physics-level selector. | 81 lines · 4 min |
| [BINDING_PHYSICS_R6.md](BINDING_PHYSICS_R6.md) | *How do you capture per-term energies without killing the frame rate?* The BindLog design. | 90 lines · 5 min |
| [BINDING_PHYSICS_R7.md](BINDING_PHYSICS_R7.md) | **The final state of the design series.** Read this one for the current design rationale. | 73 lines · 4 min |
| [BINDING_LOOP2_DONE.md](BINDING_LOOP2_DONE.md) | **What actually shipped, and what did it measure?** The S1–S7 close-out record: the re-bench table, the per-step validation, and the calibration boundary. The largest document in the repo; skim its §1–§3. | 1 325 lines · 20 min skim |

---

## Q8 — “How does it compare to a real engine on speed?”

This one lives outside `docs/` because it is a benchmark against an external
tool rather than a description of this one.

| Document | The question it actually answers | Read |
|---|---|---|
| [../bench/vs_gromacs.md](../bench/vs_gromacs.md) | *Are you faster or slower than GROMACS?* The honest answer: **~50× slower per force evaluation, ~100× faster to first interactive visualization.** The defensible niche is zero-install pedagogy, not throughput. | 68 lines · 4 min |
| [pareto_frontier.csv](pareto_frontier.csv) | *The raw accuracy/speed tier data* behind the Q7 Pareto claims — ms/step, heap MB, pose-recovery top-1, and mean rank per physics tier. Machine-readable; not a document. | 22 lines data |

---

## Non-document assets under `docs/`

| Asset | What it is | Read |
|---|---|---|
| [blueprint/PmfTrajectoryScrubber.tsx](blueprint/PmfTrajectoryScrubber.tsx) | **A blueprint, not shipped code.** A React/TS design for a PMF-linked trajectory scrubber, unwired to `index.html` and depending on React, which this zero-dependency project does not use. Read for design intent only; do not import it. | 304 lines · 10 min |

---

## How this index is kept honest

`tests/test_docs_index.js` (FAST tier, registered in `tests/suites.js`)
re-derives the file inventory from disk on every run and asserts:

- **No orphans** — every file under `docs/` is linked from this file, and in
  particular every `*.md` under `docs/`. A new doc that nobody links is a
  failing test, not a silent omission.
- **No duplicates** — no file under `docs/` is linked more than once, so this
  index stays a map rather than becoming a second copy of the docs.
- **No broken links** — every relative link target here resolves to a file
  that exists.
- **The trust boundary stays routed** — the repo-root
  [ROADMAP.md](../ROADMAP.md) and the external benchmark
  [vs_gromacs.md](../bench/vs_gromacs.md) are both reachable from here, even
  though neither lives under `docs/`.

`docs/README.md` itself is exempt from the no-orphan rule (it is the index;
it is not expected to link to itself).

To add a doc: write it, add exactly one row above, run `npm test`.
