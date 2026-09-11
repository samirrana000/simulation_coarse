# Roadmap — Honest Scope for simulation_coarse (J100)

> **Browser v1 will NOT do:** No QM/MM, no explicit membrane, no PME in browser v1.
> This is the single-sentence scope guard that every future PR is checked against.
> For full trust boundaries see `docs/APPLICABILITY.md` and `docs/LIMITATIONS.md`.

## 1. What we will NOT do in browser v1 (explicit out-of-scope)

These are **hard no's** for `v1` (tag `v1.0-jpcb`, see `CITATION.cff:13`).
If you need them, use GROMACS / AMBER / OpenMM / NAMD / CHARMM as noted in
`docs/APPLICABILITY.md:4`.

- **No QM/MM** — No quantum Hamiltonian, no SCF, no bond breaking/formation, no
  catalysis or metal redox. The heavy force field is LJ + screened Coulomb +
  GB/SA + harmonic bonds/angles/dihedrals (`src/heavy.js:415`), not a
  semi-empirical or DFT engine. *Use instead:* CP2K, ORCA, NAMD/QM-MM, OpenMM
  + psi4 for reactive chemistry.
- **No explicit membrane** — No lipid bilayer, no lateral pressure coupling,
  no anisotropic barostat, no CHARMM36 lipid parameters. The model is
  protein-only Cα ENM (`src/forcefield.js:91`) or heavy-atom GB/SA
  (`src/heavy.js:415`) in implicit solvent. Simulating a GPCR, channel, or
  transporter with this model produces qualitatively wrong energetics.
  *Use instead:* GROMACS/CHARMM-GUI + explicit lipids, NAMD, OpenMM membrane
  builder.
- **No PME in browser v1** — No Particle-Mesh Ewald. Long-range electrostatics
  beyond 8.5 Å are truncated via `switchFunc` 6.5→8.5 Å (`src/heavy.js:32`)
  with no PME and no lattice sum. The 8.5→12 Å energy convergence test shows
  ≈5% residual error (`docs/LIMITATIONS.md:17`) — unacceptable for highly
  charged nucleic acids or membranes. *Use instead:* GROMACS/AMBER/OpenMM
  with PME for any system where `|q_net|>1e` or long-range matters.
- **No rigorous FEP/TI/MBAR** — No alchemical intermediates, no soft-cores,
  no replica exchange. The PMF is a 1-D funnel `r=|COM_lig−COM_pocket|`
  (`src/funnel.js:165`) and the kinetics are a **4-state toy**
  (`src/physics/network.js:13`), pedagogical, not converged.
- **No nucleic acids / glycosylation / PTM libraries beyond simple LJ.**

Grep check for this section:

```bash
grep -n "QM/MM" ROADMAP.md
grep -n "membrane.*no.*PME" ROADMAP.md   # "membrane, no PME" also matches
grep -n "No QM/MM, no explicit membrane, no PME" ROADMAP.md
```

## 2. What we WILL do in v1 (honest niche)

*Zero-install, interactive, multi-scale, visually steered* biophysics for
education, rapid hypothesis generation, and methods prototyping — the thesis
from `TRANSFORMATION_PLAN_100.md:1`.

- **Cα ENM + heavy-atom GB/SA that stays fast:** `bench/budget.json:1`
  `{"heavy_compute_ms":2.0,"fps":30,"cg_compute_ms":0.5}` with CI warn
  (`docs/PERFORMANCE.md:47` budget section, `.github/workflows/check.yml:19`).
- **Well-tempered funnel metadynamics on 1-D `r`** with Jacobian `2kT ln r`
  and `c(t)` convergence diagnostics (HUD grays `ΔG` until `nHills≥50`).
- **4-state Chemical Network + TPT** strictly for teaching timescale bridging,
  not publishable `k_on/k_off`.
- **Canvas2D viewer with WebGL successor stub** (`src/viewer-gl.js:1`
  `WebGL not yet, fallback to Canvas2D`, `docs/VIEWER.md`).
- **Honest tutorials that fail** (`docs/TUTORIAL.md:31` Where it fails —
  charged ligand `ΔΔG≈±2–4 kcal/mol`, membrane `ΔG` meaningless).

## 3. Timeline

| Milestone | Date | Artifact | Status |
|---|---|---|---|
| `v1.0-jpcb` | 2026-09-01 | Tag + `CITATION.cff:13` + Zenodo `10.5281/zenodo.XXXXXXX` placeholder (`docs/CITATION.md`) | Pending deposition — placeholder in `manuscript/manuscript_jpcb.tex:359` |
| Reproducibility | 2026-09-01 | `manuscript/reproduce.sh` (bench/perf + figs), `docs/TUTORIAL.md` Where it fails | Done (J92, J93) |
| Accessibility | 2026-09-01 | `docs/ACCESSIBILITY.md` axe 0 violations, ESC + Space in `index.html:280`/`src/main.js:477` | Done (J98) |
| Performance budget | 2026-09-01 | `bench/budget.json:1`, `docs/PERFORMANCE.md:47`, `check.yml:19` CI warn | Done (J99) |
| Browser v1 scope lock | 2026-09-01 | This file `ROADMAP.md:3` **No QM/MM, no explicit membrane, no PME** | Done (J100) |

## 4. What might come after v1 (not promised)

- WebGL instanced spheres with depth buffer (`docs/VIEWER.md:16`) if ≥1.5× fps at 1308 atoms.
- Neighbor-list skin 2 Å, rebuild every 10 steps, 20% cut (`docs/PERFORMANCE.md:18`).
- PLUMED CV compatibility note (`docs/BRIDGE.md:??`) — 2-D CVs remain out of scope.
- None of the above overrides **No QM/MM, no explicit membrane, no PME in browser v1**.

## 5. How to propose a scope change

A PR that adds QM/MM, membrane, or PME must also:

1. Update this file, `docs/APPLICABILITY.md`, `docs/LIMITATIONS.md`, and
   `README.md:29` Limitations, and
2. Add a benchmark vs GROMACS/AMBER PME with RMSE and performance numbers in
   `bench/vs_gromacs.md`, and
3. Pass the same `bench/budget.json:1` budget without regressing fps <30.

Absent that, the PR will be closed as out-of-scope — intentionally.

*Last updated: 2026-09-01 — J100. Measurable: `grep -n "QM/MM\|membrane.*no.*PME" ROADMAP.md` must hit.*
