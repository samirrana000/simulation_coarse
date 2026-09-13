# CG Protein Simulator

A fully client-side coarse-grained (CG) molecular-dynamics app that simulates
protein structure and **protein–ligand binding** in the browser. Everything runs
in your browser — physics, rendering, and analysis are plain ES modules with no
build step and no backend. Built around the classic T4 lysozyme L99A + benzene
system (PDB **4W52**) as a working binding demo.

- One bead per residue at Cα, plus explicit heavy atoms for ligands
- **All-atom heavy mode** (panel 2): full PDB structures with covalent
  topology, metal coordination, and element-wise LJ + screened Coulomb
- Langevin dynamics integrated with the BAOAB scheme
- Elastic network model (ENM) + backbone bonds/angles (Cα mode)
- Protein–ligand binding potentials: cross 12-6 LJ, screened electrostatics,
  H-bonds, EEF1-style desolvation, and native "holo" pose springs
- **Ligand library + viewer placement** (panel 3): pick a molecule, click the
  viewer, and it is dropped rigidly clash-free at the protein surface
- Funnel bias + well-tempered metadynamics to reconstruct the **binding PMF**
- Optional ML tier: an ESM contact-prior (from Python) and an MLP pose scorer
- Post-run **analysis**: B-factor correlation, essential dynamics + RMSIP,
  ligand occupancy, contact lifetimes, PMF export
- **Phase 1 — Physics rigor (opt-in):** AMBER ff14SB-style bonds/angles/Fourier dihedrals (`src/physics/forcefield/amber14sb.js`), Tirion distance-weighted ANM + SS dihedrals (`tirion_anm.js`), GB-OBC II + HCT Born radii (`src/physics/solvation/gb_obc2.js`), analytical LCPO SASA (`lcpo_sasa.js`), implicit membrane slab (`membrane_slab.js`)
- **Phase 2 — Chemistry (opt-in):** heuristic pKa/tautomer assigner HIE/HID/HIP, ASH/GLH, CYX (`src/chem/protonation.js`), GAFF2-lite typer + Gasteiger/BCC-lite charges (`src/chem/gaff2_mapper.js`), chirality/planarity checks (`src/chem/stereo.js`), metal coordination polyhedra (`src/chem/metals.js`)
- **Phase 3 — Compute:** WGSL cell-list + tiled nonbonded kernels (`src/compute/wgsl/`), `webgpu_backend.js` orchestrator (GPU→worker→CPU fallback), r-RESPA multi-timestep integrator (`src/physics/integrators/respa.js`, 1 fs inner / 2–4 fs outer)
- **Phase 4 — Workflows:** alanine scanning ΔΔG (`src/analysis/alanine_scanning.js`), DCCM + heatmap (`src/analysis/dccm.js`), cryptic-pocket MetaD + volume tracking (`src/analysis/cryptic_pockets.js`), SMD pulling + Jarzynski ΔF (`src/analysis/unbinding_smd.js`)
- **Phase 5 — Anti-slop cockpit:** Load → Build → Run → Record visible by default, everything else collapsed; top status bar, metrics HUD, bottom dock with timeline/CV strips, keyboard `[Space]/[R]/[M]/[1-7]`; see `WHAT_CHANGED.md`

## When to use / when not to use

> **Trust boundary:** `simulation_coarse` is a browser-based educational / hypothesis tool, **NOT** a replacement for rigorous MD or free-energy perturbation (FEP). Numbers labeled `ΔG`, `K_D`, or `k_on/k_off` are qualitative — do not report them as binding affinities.

- **Use it for:** interactive teaching, rapid pocket exploration and clash-free ligand placement, illustrative PMF sketches, and offline browser demos.
- **Do NOT use it for:** production FEP/TI, quantitative `K_D`/`ΔG`, membrane or nucleic-acid systems, QM/MM, or PME electrostatics.
- **Long-range electrostatics:** cutoff 8.5 Å, no PME — not for highly charged systems (5% convergence test, see docs/LIMITATIONS.md).

See **[docs/APPLICABILITY.md](docs/APPLICABILITY.md)** for the full critical comparison vs **GROMACS, AMBER, OpenMM, NAMD, Rosetta, Mol*, NGL** (Accuracy / Speed / Use case / Install) and the explicit **“What we will NOT do”** scope. See also **[docs/LIMITATIONS.md](docs/LIMITATIONS.md)** and **[docs/CHARGES.md](docs/CHARGES.md)** (AMBER ff14SB approximate).

---

## Quick start

The app uses ES modules, so it must be served over HTTP (not opened as `file://`).

```bash
cd simulation_coarse
python3 -m http.server 8123
# open http://127.0.0.1:8123/
```

### The 4W52 binding demo (2 minutes)

0. **First run?** Open **Structure → Getting started — 4-step checklist** and
   click **Load 4W52 sample (1 click)** (works offline; the ✓/○ steps track
   Load → Build → Run → Analyze live, and every empty panel tells you the
   next step until data arrives).
1. Leave **"Include ligands (HETATM + CONECT)"** checked.
2. Type `4W52` in the **PDB ID** box and click **Fetch (RCSB / PDBe)**.
   (Alternatively drop a local `.pdb` file into **PDB file**.)
   The structure summary should read `164 Cα beads · … · 21 ligand atom(s) in 2 molecule(s) · 63 holo contacts`
   (benzene BNZ + a HEPES buffer molecule).
3. **Alternative: separate ligand file.** Use a protein-only PDB (HETATM
   ligands unchecked) and supply the ligand as a Tripos MOL2 file via
   **Ligand MOL2 file** — it *overrides* any PDB HETATM ligands. A ready-made
   `benzene.mol2` for the 4W52 pocket ships in the project root. MOL2
   hydrogens and `Du` dummy atoms are dropped (united-atom model); element is
   taken from the SYBYL atom type (`C.ar` → C, `Cl` → CL) and bond types are
   ignored (connectivity only). One file may hold several `@<TRIPOS>MOLECULE`
   blocks. Clearing the file picker falls back to HETATM ligands.
4. In panel **3 — Force Field & Integrator**, tick **Funnel bias
   (well-tempered metadynamics → PMF)**.
5. Click **▶ Run**. Watch the HUD: `CV`, `ligRMSD`, `U_bind`, `nContacts` and the
   growing `ΔG ≈` value; the PMF curve accumulates in panel **5 — Binding PMF**.
6. Click **● Rec**, keep running, then **Stop** (panel **4 — Recording**).
7. In panel **6 — Analysis** click **Analyze recorded trajectory** for the
   full report (B-factors, RMSIP, occupancy, contact lifetimes, ΔG).

Other built-in examples: `1UBQ` (apo monomer), `1CRN`, `4HHB` (multi-subunit),
`1AKE` (hinge motion).

---

## What the UI does

| Panel | What it controls |
|-------|------------------|
| **1 — Structure** | Fetch a PDB by ID or load a file; toggle HETATM/CONECT ligand parsing; optional separate **MOL2 ligand file** (overrides HETATM) |
| **2 — Selection** | Restrict to chains and residue ranges; builds the system. **Model** toggle chooses Cα elastic network vs **all-atom heavy mode** |
| **3 — Force Field & Integrator** | ENM cutoff `Rc`, spring strength `γ`, temperature, friction `ζ`, bead mass; toggles for binding potentials, holo contacts, funnel bias; the collapsed **ML tier** sub-panel; motion-gain is view-only |
| **Ligand Library & Placement** | Pick a built-in molecule (benzene, phenol, indole, caffeine, …), filter, and **Place on viewer** to drop it clash-free at the protein surface |
| **4 — Recording** | Capture a trajectory (frame stride in ps, max frames); export as multi-frame XYZ or PDB |
| **5 — Binding PMF** | Live well-tempered metadynamics PMF reconstruction; **Reset PMF** |
| **6 — Analysis** | Run the post-hoc analysis on the recorded frames; export the PMF as CSV |

The HUD (right, above the viewer) shows, live: simulation time, total potential
`U`, `U_bind` (protein–ligand energy), collective variable `CV`, ligand RMSD,
`nContacts` (< 5.5 Å), `ΔG` estimate, protein `RMSD`, instantaneous temperature,
mean displacement/frame, and FPS.

---

## Units & key parameters

| Quantity | Unit |
|----------|------|
| Length | Å |
| Time | ps |
| Energy | kcal/mol |
| Mass | Da |

| Parameter | Default | Meaning |
|-----------|---------|---------|
| ENM cutoff `Rc` | 10 Å | residue pairs within `Rc` get an elastic spring |
| Spring strength `γ` | 2.0 kcal/mol/Å² | ENM stiffness |
| Temperature | 300 K | Langevin heat bath (well-tempered metaD uses it too) |
| Friction `ζ` | 8.0 ps⁻¹ | Langevin friction (sets the timestep) |
| Bead mass | 110 Da | uniform Cα mass (ligand atoms use element masses) |

The timestep is auto-set from the fastest bond and `ζ`; the sim targets
~1 ps of dynamics per displayed frame.

---

## Energy model

```
U = Σ ½k_b(r−r₀)²                          peptide bonds
  + Σ ½k_θ(θ−θ₀)²                           backbone + ligand angles
  + Σ ½γ·H(Rc−r₀)(r−r₀)²                   elastic-network contacts
  + Σ ε[(r_e/r)¹²−2(r_e/r)⁶+1], r<r_e       excluded volume (repulsion)
  + Σ 4ε[(σ/r)¹²−(σ/r)⁶]                    protein–ligand LJ
  + Σ 332·q₁q₂/(ε_r·r), ε_r = 4+76·tanh(r/8)  screened electrostatics
  + Σ −ε_HB·exp(−(r−3.2)²/0.72)             H-bonds (polar donor/acceptor)
  + Σ ΔG_t·(1−e^(−n/3))                     EEF1-lite desolvation (burial)
  + Σ ½γ_lig(r−r₀)², r₀ ≤ 6 Å               holo pose springs (bound pose)
  + funnel wall ½k_f(r−rFlat)² (r > rFlat)  funnel bias (flat in the pocket)
  + well-tempered metaD hills               PMF reconstruction
```

Ligand internal geometry (bonds, angles, aromatic-ring planarity) is built from
CONECT records (or from the MOL2 bond list when the ligand comes from a
separate `.mol2` file) with automatic gap-filling for heavy atoms closer than
1.8 Å.
The funnel bias is exactly flat inside the bound state, so the bound equilibrium
is undisturbed while the ligand is pulled into the pocket from far away.

---

## Heavy mode (all-atom)

Set panel 2's **Model** to *All-atom heavy mode* to switch from the Cα elastic
network to a full heavy-atom model (`src/heavy.js`). It parses **every** heavy
atom (protein + HETATM ligands + **metal ions**; water and free ions are
dropped), builds covalent topology by geometry, adds metal-coordination
springs, and runs the same BAOAB integrator via a `HeavyForceField` that
shares `ForceField`'s public interface.

```
U = Σ 1⁄2·200(r−r0)2            covalent bonds (geometric detection)
  + Σ 1⁄2·40(θ−θ0)2              angles
  + Σ 1⁄2·20(φ−φ0)2              improper dihedrals (planarity)
  + Σ 1⁄2·2(φ−φ0)2               proper dihedrals
  + Σ 1⁄2·40(r−r0)2              metal coordination (N/O/S donors → M)
  + Σ 4ε[(σ/r)12−(σ/r)6]         LJ (geometric mixing)
  + Σ 332·qi qj e^(−r/6)/r2      screened Coulomb (ε(r)=r, neutralized)
```

Bonds are detected when two heavy atoms sit within `1.15 × (covalent radii)`
(with a 2.2 Å hard cap); angles/impropers/propers follow from the bond graph.
Non-bonded 1-2/1-3 pairs are excluded and 1-4 pairs scaled 0.5. A smooth
AMBER-style switching cutoff (6→8 Å) bounds the O(n2) pair scan. Metal ions
(the `METAL_ELEMENT` table: Zn, Fe, Mg, Ca, Cu, Mn, Ni, Co, Na, K) never form
covalent bonds — instead each metal collects up to `coordN` nearby N/O/S donors
within `coordR` and holds them with k=40 springs, so the ion keeps its
coordination geometry while still feeling LJ. Charges are mean-neutralized so
the net charge is zero and electrostatics is a local dipole/quadrupole term.

The viewer element-colors heavy atoms (metals get CPK-inspired colors) and
draws the full covalent bond graph. Heavy mode is intentionally
parameter-light and geometry-driven — a structural/visual exploration layer,
not production MD.

---

## Ligand library & placement

Panel 3 ships a built-in **ligand library** (`src/ligandLib.js`): benzene,
phenol, toluene, chlorobenzene, indole, imidazole, acetate, ethanolamine,
DMSO, and caffeine as inline MOL2. Filter the list, pick one, click **Place on
viewer**, then click anywhere on the protein — the ligand is added to the
system and lowered onto the protein surface with a rigid clash-relaxation so
the final pose is clash-free (holo pose springs are OFF for placed poses, since
they are hypotheses rather than the crystallographic pose).

---

## ML tier

Two optional neural ingredients live in the collapsed **ML tier** sub-panel:

### 1. NN contact map (ESM prior)

A JSON file with a contact prior used to **scale the ENM spring constants**:

```json
{ "source": "heuristic", "chain": "A", "nRes": 164, "contacts": [[i, j, p_ij], ...] }
```

Spring `(i, j)` gets stiffness `γ·(1 + α·p_ij)`. Load it via **contacts.json**,
then tick **NN contact map (ESM prior)**; `… · NN map active` appears in the
selection summary. A ready-made prior for 4W52 ships in
`data/4W52_contacts.json`.

Generate your own prior for a structure with the Python exporter:

```bash
python3 ml/export_esm_contacts.py --fetch 4W52 --out data/4W52_contacts.json
# heuristic mode (no ESM weights needed on this machine)
python3 ml/export_esm_contacts.py --pdb my.pdb --out contacts.json --min-p 0.05 --max-pairs 0
# ESM mode (requires fair-esm / torch):  --model esm2_t6_8M_UR50D
```

### 2. NN pose score (MLP)

A small feed-forward MLP that scores the instantaneous pose from 6 features:
`nContacts`, burial (`−desolvU`), clash count, ligand RMSD, `CV`, and `−U_bind`.
Load a weights file via **scorer weights (JSON)**, tick **NN pose score (MLP)**,
and the live `pose =` value appears in the HUD. Weights format:

```json
{
  "layers":   [6, 4, 1],
  "W":        [[row-major 4×6 weights], [1×4 weights]],
  "b":        [[4 biases], [0]],
  "act":      "silu",              // "silu" | "relu" | "tanh" | "none"
  "featMean": [0,0,0,0,0,0],       // optional
  "featStd":  [1,1,1,1,1,1]        // optional
}
```

With no file loaded, ticking the box uses a built-in linear "docked-pose" scorer.

---

## Analysis (panel 6)

Run on the recorded trajectory (need ≥ 2 frames — press **● Rec**, run, **Stop**):

- **B-factors** — per-residue mean-square fluctuation → `B = 8π²/3 ⟨Δr²⟩`,
  Pearson-correlated against the experimental temperature factors from the PDB.
- **Essential dynamics + RMSIP** — frames are Kabsch-superposed onto the native
  structure; the top-k principal components of the covariance are compared with
  the k softest ENM normal modes (rigid-body modes removed):

  ```
  RMSIP = sqrt( (1/k) Σ_{i,j≤k} (v_i·w_j)² )
  ```

- **Ligand occupancy** — a 1 Å atom-density grid over the ligand volume: the
  peak occupancy cell vs. the crystal ligand COM (pose recovery), plus the
  fraction of frames the ligand COM stays inside the pocket (bound fraction).
- **Contact lifetimes** — protein–ligand pairs within 6 Å tracked frame to
  frame; mean/max uninterrupted residence streaks in ps and the most persistent
  contacts.
- **Binding free energy** — the well-tempered metadynamics ΔG and a
  downloadable **PMF CSV** (`r, pmf` + ΔG footer).

---

## Binding-physics Loop-2 + Stages 1–7

Loop-2 upgraded the binding physics in 7 additive, opt-in stages (default L0
path bit-identical; full record in `docs/BINDING_LOOP2_DONE.md`):

- **Physics-level selector** (panel 3, persisted to `localStorage`): **L0**
  fast CG (default) · **L1** balanced (CG formal charges + directional-HB
  virtual sites, ≈1.09×) · **L2** full-rigor (heavy π-stack / cation-π /
  halogen weak terms + per-term accumulators, ≈1.0× over tier base).
- **BindLog** (panel 5/6 wiring): per-term `bindU` vector (CG 4-term, heavy
  7-term), tick frames, energy events, contact diffs, funnel hills — off
  unless enabled or L2 selected.
- **BindViz** (collapsed subpanels): ≤1 Hz timeline/energy/PMF strips with a
  null-bindlog guard.
- **Thermo panel**: ΔH ± SE with component split + Schlitter/torsion −TΔS
  (async chunked apo leg, progress in the caption).
- **Seeded determinism**: opt-in `seed` on the Langevin integrator
  (mulberry32); default unseeded path unchanged.
- **Calibration boundary** (`scripts/calibration_4w52.mjs`): 4W52 anchor
  ΔG_est ≈ −2.2 vs experimental −5.2 (ITC) / −4.2 (NMR) kcal/mol — replica SD
  (±7, entropy-driven) dominates, so ΔG is **ranking-only, NOT FEP, NOT
  absolute Kd**.
- **Perf honesty**: `node [--expose-gc] scripts/pareto_bench.mjs` →
  `docs/pareto_frontier.csv` (ms/step + gc-bracketed heap deltas ≤0.3 MB;
  history rows preserved). Loop-2 binding terms are CPU-only by decision
  (see `docs/BINDING_LOOP2_DONE.md` §14) — GPU covers LJ/Coulomb/GB.
- **Tests**: `node tests/test_all.js` → **215/215 fast** (~13 s, gate);
  `--slow` → **245/245** (+ thermo 7 + heavy 13 + calibration 10).

---

## Limitations

> **Honest scope:** This browser model is not production MD — see the explicit list in **[docs/LIMITATIONS.md](docs/LIMITATIONS.md)** for full detail.

- **GB cutoff, no PME** — 6.5→8.5 Å switching, no PME (`src/heavy.js:32`); not for highly charged systems.
- **1D funnel CV only** — PMF along `r = |COM_lig − COM_pocket|` (`src/funnel.js:165`); orthogonal barriers invisible.
- **4-state kinetics toy** — Bulk/Encounter/Intermediate/Bound Kramers network (`src/physics/network.js:13`); illustrative, not a converged MSM. **4-state toy, not full MSM; use PyEMMA for production** — see `docs/NETWORK.md` (`src/physics/network.js:13`).
- **Canvas2D vs WebGL** — 2-D Canvas painter's sort (`src/viewer.js:481`), no depth buffer; see `docs/VIEWER.md` and `src/viewer-gl.js`.
- **No membrane / nucleic acids / QM** — protein-only force fields; use GROMACS/CHARMM/NAMD/QM-MM for those.

For the comparative table vs GROMACS/AMBER/OpenMM see **[docs/APPLICABILITY.md](docs/APPLICABILITY.md)**. For performance vs GROMACS see **[bench/vs_gromacs.md](bench/vs_gromacs.md)**. Full limitations also in **[docs/LIMITATIONS.md](docs/LIMITATIONS.md)**.

---

## Project layout

```
index.html          UI (panels 1–6)
css/style.css       styling
4w52.pdb            demo structure (also the t22 fixture)
benzene.mol2        demo ligand (Tripos MOL2) — separate-file MOL2 input
src/
  main.js           orchestration: wiring, main loop, HUD, ML tier
  pdb.js            PDB parsing (Cα beads, ligands/HETATM, B-factors) + MOL2 ligands
  forcefield.js     the CG energy function (ENM, binding, desolvation, springs)
  heavy.js          all-atom heavy mode (parseHeavy + HeavyForceField + metals)
  ligand.js         ligand internal geometry (bonds/angles/ring planarity)
  ligandLib.js      built-in ligand library (10 molecules as inline MOL2)
  ligand-panel.js   library selector + viewer placement flow
  placement.js      rigid clash-free placement (placement.js)
  integrator.js     BAOAB Langevin integrator
  funnel.js         funnel bias + well-tempered metadynamics → PMF/ΔG
  scorer.js         MLP pose scorer
  analysis.js       post-run analysis (B-factors, RMSIP, occupancy, lifetimes, PMF)
  recorder.js       trajectory capture + XYZ/PDB export
  viewer.js         Canvas rendering (Cα beads + heavy atoms, picking)
  ff-params.js      element/LJ/charge tables + metal & covalent-radius tables
  ff-harmonic.js    harmonic kernels (bonds/angles/springs)
  ff-repulsion.js   excluded-volume repulsion
  ff-binding.js     protein–ligand binding potentials
  ml-tier.js        NN contact prior + MLP pose scorer wiring
  pmf-panel.js / analysis-panel.js   panel 5 / 6 side-effect modules
  physics/forcefield/amber14sb.js + tirion_anm.js   Phase-1 tables & CG ANM
  physics/solvation/gb_obc2.js + lcpo_sasa.js + membrane_slab.js
  physics/integrators/respa.js   r-RESPA stepper
  chem/protonation.js + gaff2_mapper.js + stereo.js + metals.js
  compute/wgsl/cell_list.wgsl + nonbonded_forces.wgsl + webgpu_backend.js
  analysis/alanine_scanning.js + dccm.js + cryptic_pockets.js + unbinding_smd.js
ml/export_esm_contacts.py   Python contact-prior exporter (ESM or heuristic)
data/4W52_contacts.json     ready-made contact prior for the demo
docs/PHYSICS_RIGOR.md + CHEMISTRY.md + GPU_RESPA.md + WORKFLOWS.md + UI_COCKPIT.md
  Phase 1–5 subsystem docs (see docs/)
```

---

## Development notes

- ES modules are cache-busted with a query suffix (`?v=10`); bump it after
  editing any module so browsers fetch fresh code.
- Modules import each other with that suffix; Node strips the query string on
  `file://` URLs, so the modules also run under plain Node for unit tests.
- Unit/browser harnesses live under `/tmp/opencode/` (`t2`–`t22`, `probe4w52`).
  The browser tests drive the app over `http://127.0.0.1:8123/` with
  puppeteer-core and the bundled Chrome.
- The ML exporter runs in **heuristic mode** when the ESM weights are absent
  (`fair-esm` is not installed); it still produces a valid contact prior.

---

## Version, license & citation

- **Version:** `1.1.0-fp7` (2026-09-13) — canonical constant in
  `src/version.js` (this repo has no `package.json`; the constant feeds the
  HUD prefix, trajectory `REMARK` provenance, and the dock footer). Release
  history in `CHANGELOG.md`.
- **License:** MIT — see [LICENSE](LICENSE) (© 2026 Samir Rana).
- **Citation:** see [CITATION.cff](CITATION.cff) (`v1.0-jpcb`, Zenodo DOI
  placeholder `10.5281/zenodo.XXXXXXX` until deposition) and
  [docs/CITATION.md](docs/CITATION.md) for BibTeX; validation numbers in
  [docs/VALIDATION.md](docs/VALIDATION.md). Trust boundary: ranking-only use —
  NOT FEP, NOT absolute Kd.
