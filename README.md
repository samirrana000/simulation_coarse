# CG Protein Simulator

A fully client-side coarse-grained (CG) molecular-dynamics app that simulates
protein structure and **protein–ligand binding** in the browser. Everything runs
in your browser — physics, rendering, and analysis are plain ES modules with no
build step and no backend. Built around the classic T4 lysozyme L99A + benzene
system (PDB **4W52**) as a working binding demo.

- One bead per residue at Cα, plus explicit heavy atoms for ligands
- Langevin dynamics integrated with the BAOAB scheme
- Elastic network model (ENM) + backbone bonds/angles
- Protein–ligand binding potentials: cross 12-6 LJ, screened electrostatics,
  H-bonds, EEF1-style desolvation, and native "holo" pose springs
- Funnel bias + well-tempered metadynamics to reconstruct the **binding PMF**
- Optional ML tier: an ESM contact-prior (from Python) and an MLP pose scorer
- Post-run **analysis**: B-factor correlation, essential dynamics + RMSIP,
  ligand occupancy, contact lifetimes, PMF export

---

## Quick start

The app uses ES modules, so it must be served over HTTP (not opened as `file://`).

```bash
cd simulation_coarse
python3 -m http.server 8123
# open http://127.0.0.1:8123/
```

### The 4W52 binding demo (2 minutes)

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
| **2 — Selection** | Restrict to chains and residue ranges; builds the Cα system |
| **3 — Force Field & Integrator** | ENM cutoff `Rc`, spring strength `γ`, temperature, friction `ζ`, bead mass; toggles for binding potentials, holo contacts, funnel bias; the collapsed **ML tier** sub-panel; motion-gain is view-only |
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
  ligand.js         ligand internal geometry (bonds/angles/ring planarity)
  integrator.js     BAOAB Langevin integrator
  funnel.js         funnel bias + well-tempered metadynamics → PMF/ΔG
  scorer.js         MLP pose scorer
  analysis.js       post-run analysis (B-factors, RMSIP, occupancy, lifetimes, PMF)
  recorder.js       trajectory capture + XYZ/PDB export
  viewer.js         Canvas rendering
ml/export_esm_contacts.py   Python contact-prior exporter (ESM or heuristic)
data/4W52_contacts.json     ready-made contact prior for the demo
```

---

## Development notes

- ES modules are cache-busted with a query suffix (`?v=8`); bump it after
  editing any module so browsers fetch fresh code.
- Modules import each other with that suffix; Node strips the query string on
  `file://` URLs, so the modules also run under plain Node for unit tests.
- Unit/browser harnesses live under `/tmp/opencode/` (`t2`–`t22`, `probe4w52`).
  The browser tests drive the app over `http://127.0.0.1:8123/` with
  puppeteer-core and the bundled Chrome.
- The ML exporter runs in **heuristic mode** when the ESM weights are absent
  (`fair-esm` is not installed); it still produces a valid contact prior.
