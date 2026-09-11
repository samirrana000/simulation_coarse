# MDAnalysis Bridge — DCD/PDB Export & RMSF Parity (I83)

This note documents how to load trajectories exported from the browser
(`src/recorder.js:buildFile` → **PDB multi-MODEL** or **XYZ**, plus the
recommended **DCD** via CPPTRAJ/MDAnalysis conversion) in **MDAnalysis**
and reproduce the same **RMSF → B-factor** analysis that the in-browser
`src/analysis.js:analyzeTrajectory` computes.

## Export formats from the recorder

| Format | Produced by | Typical file |
|--------|-------------|--------------|
| PDB multi-MODEL | `recorder.buildFile("pdb", beads)` | `cg_traj_500frames.pdb` |
| XYZ multi-frame | `recorder.buildFile("xyz", beads)` | `cg_traj_500frames.xyz` |
| DCD (after conversion) | `MDAnalysis` / `MDTraj` write | `trajectory.dcd` |

The recorder stores each frame as `Float32Array(3n)` + `times[]` (ps).
PDB export writes `MODEL/ENDMDL` with `REMARK time = <ps> ps` and one
`CA` ATOM per Cα bead using `beads[i].resName/chain/resSeq` so re-import
into PyMOL/VMD/MDAnalysis preserves residue identity.

## Loading in MDAnalysis

### Option A — PDB multi-MODEL directly (simplest, no conversion)

```python
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF
import numpy as np

# 1) Export from browser: Recording panel → Format=PDB → Download
#    Save as cg_traj_500frames.pdb (contains MODEL/ENDMDL blocks)
#    Reference (native) structure: original PDB used to Build System

u = mda.Universe("cg_traj_500frames.pdb")  # topology+trajectory in one PDB
# or, if you exported topology separately:
# u = mda.Universe("topology.pdb", "cg_traj_500frames.pdb")

ca = u.select_atoms("name CA")             # Cα selection matches src/analysis.js nProt
print(f"{len(ca)} Cα atoms, {len(u.trajectory)} frames")

# Align trajectory to first frame (Kabsch in MDAnalysis == src/analysis.js:kabsch)
from MDAnalysis.analysis import align
align.AlignTraj(u, u, select="name CA", in_memory=True).run()

# RMSF per residue (Å)
rmsf = RMSF(ca).run()
rmsf_vals = rmsf.results.rmsf  # shape (nProt,), Å

# B-factors matching src/analysis.js: B = (8π²/3) * <Δr²>
# where <Δr²> == RMSF²  (mean-square fluctuation about mean)
B_sim_mda = (8 * np.pi**2 / 3.0) * (rmsf_vals ** 2)
print(B_sim_mda[:5])

# Pearson vs experimental B-factors (if topology has tempfactors)
B_exp = ca.tempfactors  # Å², from PDB column 61-66 (beads[i].bfac)
from scipy.stats import pearsonr
mask = B_exp > 0
r, pval = pearsonr(B_sim_mda[mask], B_exp[mask])
print(f"Pearson r(B_sim, B_exp) = {r:.3f} over {mask.sum()} residues")
# Should match report line [1] in src/analysis.js:analyzeTrajectory
```

### Option B — PDB topology + DCD (recommended for large trajectories)

If you need DCD (more compact, faster), convert once:

```python
# One-time conversion: PDB multi-MODEL → DCD + topology PDB
import MDAnalysis as mda
u = mda.Universe("cg_traj_500frames.pdb")
# write first frame as topology
u.trajectory[0]
ca = u.select_atoms("all")
ca.write("topology.pdb")
# write full trajectory as DCD
with mda.Writer("trajectory.dcd", n_atoms=len(ca)) as W:
    for ts in u.trajectory:
        W.write(ca)
```

Then analysis is identical with `Universe("topology.pdb", "trajectory.dcd")`.

### Minimal XYZ path

```python
u = mda.Universe("topology.pdb", "cg_traj_500frames.xyz")  # XYZ reader is built-in
# then same AlignTraj + RMSF as above
```

## RMSF → B-factor parity check

`src/analysis.js:371` defines:

```js
const RMSF_TO_B = (8 * Math.PI * Math.PI) / 3; // B = 8π²/3 ⟨Δr²⟩ ≈ 26.32
Bsim[i] = RMSF_TO_B * (msf[i] / nF);          // msf = Σ Δr² over superposed frames
```

`analysis.js` superposes every frame onto `ref` via
`kabsch(fr, ref, nProt)` before accumulating `msf`. The MDAnalysis
snippet above replicates this via `AlignTraj(..., select="name CA")`
(Kabsch least-squares). If you skip alignment, RMSF will be inflated by
rigid-body diffusion — always align first.

Expected equivalence (within Float32 round-off):

```
B_sim (analysis.js)  ≈  (8π²/3) * RMSF(MDA)²
Pearson r(B_sim, B_exp) should agree to ±0.01
```

## Troubleshooting & validation

- **Frame count mismatch:** `len(u.trajectory)` should equal `recorder.count`
  and `recorder.times.length`. Check `REMARK time =` lines in the PDB.
- **Units:** CG model units are Å / ps / kcal/mol — MDAnalysis reads Å
  natively, no conversion needed.
- **Chains:** blank chain `_` is exported as `" "` (PDB convention) and
  parsed back as chain `_` — selection `name CA` is chain-agnostic.
- **Ligand:** for Cα+ligand systems (`nProt + nLig`), use
  `u.select_atoms("name CA or resname BEN")` analogously to
  `src/analysis.js` `nProt`/`ligStart` split.

## References

- MDAnalysis RMSF: https://docs.mdanalysis.org/stable/documentation_pages/analysis/rms.html
- In-browser report: `src/analysis.js:analyzeTrajectory` — sections [1] B-factors,
  [2] essential dynamics/RMSIP, [3] ligand occupancy, [4] contacts, [5] ΔG
- Recorder I/O: `src/recorder.js:Recorder`, `src/recorder.js:getFrame`

## Grep hits

- `docs/BRIDGE.md` — this file (MDAnalysis bridge)
- `notebooks/mdanalysis_bridge.ipynb` — optional notebook mirror (if present)
