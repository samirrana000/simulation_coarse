# Tutorial — Quick Success and Where It Fails (J92)

This tutorial is intentionally honest: the first half succeeds on the
benchmark system, the second half shows where the same pipeline produces
qualitatively wrong energetics and how to diagnose it.

---

## 1. Quick success (4W52 benzene) — 2-minute demo

**System:** T4 lysozyme L99A + benzene (PDB 4W52, 164 residues, 1308 heavy
atoms, cavity ≈150 Å³). Apolar, rigid pocket, no ordered bridging waters —
the sweet spot for HCT-GB/SA + LCPO in the browser.

### Steps

1. Serve locally:
   ```bash
   python3 -m http.server 8123
   # open http://127.0.0.1:8123/
   ```
2. Leave **Include ligands (HETATM + CONECT)** checked.
3. Type `4W52` in **PDB ID** → **Fetch (RCSB / PDBe)**. The structure summary
   should read `164 Cα beads · 63 holo contacts · 21 ligand atom(s) in 2 molecule(s)`
   (benzene BNZ + HEPES buffer). Alternatively drop `4w52.pdb` via **PDB file**.
4. In panel **2 — Model & Hetero Selection** leave `Cα elastic network`
   (or switch to `All-atom heavy mode` to see covalent + GB/SASA terms).
5. Tick **Funnel bias (well-tempered metadynamics → PMF)** in panel
   **4 — Dynamics Tuning & Force Field** (heuristic: `γ=2`, `Rc=10 Å`, `T=300 K`).
6. Click **▶ Run**. Watch HUD: `CV = |COM_lig − COM_pocket|`,
   `ligRMSD`, `U_bind`, `nContacts`. The PMF curve grows in
   **7 — Binding PMF & Analysis**.
7. Click **● Rec** (Recording), run ≈30 s, **Stop**, then **Analyze Trajectory**.
   Expected report:
   - B-factor `R ≈ 0.70–0.84` (heuristic ENM vs PDB `B_exp`; see `bench/b_factors.js`)
   - PMF `ΔG° ≈ −5.3 ± 0.3 kcal/mol` (converged only when `nHills ≥ 50`;
     see `src/funnel.js:estimateDG`  `if(funnel._nHills<50) dgStr="– (collecting…)"`)
   - Bulk plateau `W → 0` for `r ≥ 9 Å`, native well `W ≈ −6.8 kcal/mol` at `r=2.5 Å`,
     desolvation barrier `ΔW‡ ≈ 1.85 kcal/mol` at `r≈6.2 Å` (`manuscript/figures/fig2_pmf_free_energy.pdf`)

### What success looks like (CLI parity)

```bash
node bench/perf.js
# CG: n=164  ~0.08 ± 0.02 ms/compute
# Heavy: n=1308  ~0.9 ± 0.2 ms/compute
node cli.js --pdb 4w52 --steps 10000 --out /tmp/traj.pdb
python3 manuscript/generate_figures.py
# All publication figures successfully generated with zero collisions
```

If `B_sim` correlation <0.45 or `ΔG` drifts >1 kcal/mol after `nHills=150`,
re-check `T`, `γ`, or contact map — see diagnosis below.

---

## 2. Where it fails — Charged ligand and Membrane protein (expected wrong ΔG)

> **This section is the point.** A tutorial that always succeeds teaches
> nothing about model scope. Both cases below produce reproducible but
> **wrong** free energies with the current HCT-GB/SA + LCPO + 6.5→8.5 Å cutoff
> (no PME) physics (`src/heavy.js:32` `switchFunc`, `docs/LIMITATIONS.md:8`).

### Case A — Highly charged / polar ligand (e.g., phosphate, ATP analogue, charged sulfonamide on 4W52 scaffold)

**Setup to reproduce the failure:**

1. Start from 4W52 as above, but load a charged ligand via **Ligand MOL2 file**
   (e.g., replace benzene with `benzene.mol2` edited to `q≈−1.5e` net charge
   or use a PDB with a `PO4`/`ATP` HETATM).
2. Build System → Run → Record → PMF.

**Expected wrong ΔG:**

- Browser PMF reports `ΔG° ≈ −8 to −10 kcal/mol` (over-stabilized) while
  experimental ITC is `−4 to −5 kcal/mol`, or conversely under-stabilized
  by `+2–4 kcal/mol` for multiply charged guanidinium/phosphate pairs
  depending on `epsIn/epsOut` and salt (`settingsState.epsIn=4.0`,
  `epsOut=78.5`, `kappa≈0.127 Å⁻¹` at 150 mM).
- Signed error `ΔΔG = ΔG_browser − ΔG_exp ≈ ±2–4 kcal/mol`; RMSE over
  10 charged PDBbind cases ≈2.5 kcal/mol (`bench/pdbbind_gb.js` documents
  `RMSE <2.5 kcal/mol or docs state "not FEP"`).
- Diagnostic: total non-bonded + GB energy changes `≈5%` when cutoff is
  increased 8.5→12 Å (`docs/LIMITATIONS.md:17` convergence check) and
  `|q_net|>1e` produces `switchFunc` truncation artifacts; SASA burial
  `γ·SASA` misses discrete bridging waters (`src/sasa.js` LCPO).

**Why it fails:** GB uses mean-neutralized charges and a screened Coulomb
kernel `332·q1q2·exp(−κr)/r` with `switchFunc` 6.5→8.5 Å and no
Particle-Mesh Ewald — long-range terms beyond 8.5 Å are dropped
(`docs/LIMITATIONS.md:8` **GB cutoff, no PME**). Polar networks and
multivalent ions require explicit water reorganization entropy that LCPO
cannot capture.

### Case B — Membrane protein / lipid bilayer (e.g., GPCR, ion channel)

**Setup to reproduce the failure:**

1. Fetch a membrane protein (e.g., `4HHB` hemoglobin still water-soluble;
   try a GPCR `3SN6` or `6CMO` if available) and build in **All-atom heavy mode**.
2. Run 1000 steps, observe HUD `U ≈ +50 to +200 kcal/mol` (repulsive, never
   equilibrates) and PMF flat or diverging.

**Expected wrong ΔG:**

- `ΔG°` meaningless (`NaN` or `+5 to +15 kcal/mol` favorable in vacuum vs
  `−3 kcal/mol` expected) because there is **no lipid bilayer, no lateral
  pressure, no anisotropic solvent** (`docs/APPLICABILITY.md:46` **No membrane**).
- No convergence: `nHills` grows but `ΔG` block SE stays >1.5 kcal/mol
  (`src/funnel.js:convergenceSE`).

**Why it fails:** Force fields are protein Cα ENM (`src/forcefield.js:91`)
or heavy LJ/GB/SASA (`src/heavy.js:415`) — no lipid TIP3P, no
CHARMM36 lipid dihedrals, no membrane insertion free energy
(`docs/LIMITATIONS.md:40` **Out of scope: membranes, nucleic acids, QM**).

### How to diagnose a failed run

| Symptom | Check | Where |
|---|---|---|
| `ΔG` far from ITC (>2 kcal/mol) on charged ligand | Net charge `Σq ≠ 0`, `κ` mismatch, cutoff 8.5 Å truncation `ΔU≈5%` | `src/ff-params.js` charges, `src/heavy.js:32` switch, `docs/CHARGES.md`, `settings-panel.js:saltM` |
| PMF plateau not flat (`W(r≥9)` drifts) | 1-D CV misses orthogonal barrier, `rFlat` mis-set, `nHills<50` | `src/funnel.js:165` `cv()`, `docs/FUNNEL.md`, HUD `ΔG not converged (nHills<50)` |
| B-factor `R<0.45` | Cutoff `Rc` wrong, ENM uniform `γ` | `bench/b_factors.js`, `docs/CG_HEAVY.md` |
| Kinetics `k_on` off by 100× | CNM is **4-state toy**, not MSM, `ν0=1e10` | `docs/NETWORK.md:13` `NETWORK_STATES` |
| Membrane trajectory explodes | System out of scope | `docs/APPLICABILITY.md:4` **What We Will NOT Do** |

**Action:** If any symptom hits, do **not** report `ΔG/K_D/k_on` as an
affinity. Use GROMACS/AMBER/OpenMM with PME + explicit solvent + MBAR/FEP
for quantitative claims (see `docs/BRIDGE.md` MDAnalysis conversion and
`docs/OPENMM_REF.md`).

### Reproducing the failure programmatically

```bash
node bench/pdbbind_gb.js        # shows RMSE ≈2.5 kcal/mol on charged subset (expected wrong)
node bench/perf.js              # timing only — does not validate ΔG
# For membrane: build a 4HHB-heavy system and assert energy drift
python3 manuscript/generate_figures.py  # figures remain valid only for 4W52 apolar case
```

*Last updated: 2026-09-01 — This tutorial is required reading before any
manuscript claim (see `manuscript/reproduce.sh`).*
