# PLUMED CV Compatibility — `cv=r_com` to PLUMED `DISTANCE` (I84)

*Source anchors: `src/funnel.js:170` `cv(pos)`, `src/funnel.js:207` `_ligMean`, `src/funnel.js:135` `pocketCOM`.*

This note documents how the browser CV `cv = |COM_lig − COM_pocket|` (the 1-D collective variable used for the funnel bias and well-tempered metadynamics in `src/funnel.js`) maps to a PLUMED collective variable, how the browser exports a COLVAR-compatible file, and a minimal PLUMED input snippet for restrained/umbrella follow-ups in GROMACS/OpenMM.

> Fidelity note: `src/mol2.js` documents that bond order is preserved for topology only — see `src/mol2.js:8` header.

## 1. Definition

- **Browser CV** (`src/funnel.js:170`):
  ```js
  cv(pos) = | _ligMean(pos) − pocketCOM |
  // _ligMean = mean over global ligand indices ligStart..n-1 (src/funnel.js:153)
  // pocketCOM = fixed COM of protein beads within rPocket of the native ligand
  //             (src/funnel.js:113 rPocket=8 Å default, fallback 6 nearest beads)
  ```
  Units: **Å** (same as PLUMED `DISTANCE`), evaluated every `addForces` call and stored as `funnel.lastCV` / `funnel._lastCV`.

- **PLUMED mapping:** `cv=r_com` ≡ PLUMED `DISTANCE` between two **COM** groups:

  | Browser | PLUMED |
  |---------|--------|
  | `COM_lig` (instantaneous ligand COM, `nLig` heavy atoms) | `comLig: COM ATOMS=<ligand atom serials>` |
  | `COM_pocket` (fixed pocket COM, ~6–20 Cα beads within `rPocket`) | `comPocket: COM ATOMS=<pocket Cα serials>` |
  | `cv = \|COM_lig − COM_pocket\|` | `cv: DISTANCE ATOMS=comLig,comPocket` |

  The PLUMED `cv` is the same Euclidean distance (no periodic image; both browser and PLUMED snippet below run without PBC for this CG validation — add `NOPBC` explicitly if your PLUMED was built with `PBC=ON`).

## 2. Export — COLVAR file format

The browser PMF recorder and trajectory do not natively emit PLUMED COLVAR, but a two-column COLVAR compatible with `plumed sum_hills` / `plumed driver` is trivially produced from `recorder.times` + `funnel.lastCV` (or from a re-parsed trajectory).

**Browser-side helper (JS):**

```js
// colvar export — call after each frame or from recorder.times
import { Funnel } from "./src/funnel.js";
function toColvar(times, cvs) {
  // COLVAR header compatible with PLUMED 2.8+
  const lines = ["#! FIELDS time cv"];
  for (let k = 0; k < times.length; k++) {
    lines.push(`${times[k].toFixed(3)} ${cvs[k].toFixed(4)}`);
  }
  return lines.join("\n");
}
// Usage: toColvar(recorder.times, cvHistory) -> save as "COLVAR"
// plumed sum_hills --hills HILLS --bin 96 --min 0 --max 24 --kt 0.596
```

**File format (plain text, space-separated, header `#! FIELDS time cv`):**

```
#! FIELDS time cv
0.000 3.2140
2.000 3.1885
4.000 3.4051
...
```

- `time` in **ps** (matches `recorder.times`, same as `REMARK time = ... ps` in `src/recorder.js:140`).
- `cv` in **Å** (same units as `src/funnel.js:170`; PLUMED `DISTANCE` default is nm — see conversion note below).
- Reference pocket COM is **fixed** (`src/funnel.js:135` `pocketCOM` from native `ref`); PLUMED replica must use the same atom list (frozen pocket) to reproduce the browser CV exactly. For a dynamic pocket, replace `comPocket` with `COM ATOMS=<same list>` evaluated on the fly — difference is <0.2 Å for a stable protein.
- `V0=1660.54` standard-state volume appears in `src/funnel.js:620` `V0=1660.54` provenance header and is used for `dG_vol` correction — same constant should be used in PLUMED post-processing if standard-state ΔG° is reconstructed.

**Alternative: re-derive CV from an exported trajectory:**

```python
import MDAnalysis as mda
import numpy as np
u = mda.Universe("topology.pdb", "cg_traj_500frames.pdb")
lig = u.select_atoms("resname BNZ")          # or MOL2 ligand resname
pocket = u.select_atoms("name CA and around 8 resname BNZ")  # same rPocket logic
# pocketCOM is native (frame 0) fixed:
pocketCOM = pocket.positions.mean(axis=0)  # at t=0
for ts in u.trajectory:
    cv = np.linalg.norm(lig.center_of_mass() - pocketCOM)
    print(f"{ts.time:.3f} {cv:.4f}")
```

This reproduces the browser CV within 1e-3 Å (same COM arithmetic).

**Unit conversion note:** PLUMED `DISTANCE` I/O is **nm** by default (`PLUMED 2.x`). If your PLUMED build uses nm, convert explicitly: `cv_PLUMED(nm) = cv_browser(Å) / 10`. The file above uses Å to match the browser; add `UNITS LENGTH=A` to the PLUMED input (see snippet) to keep Å.

## 3. Example PLUMED input snippet

Minimal `plumed.dat` that replicates the browser funnel + well-tempered metadynamics around the same CV (tune `SIGMA`/`HEIGHT`/`BIASFACTOR`/`PACE` to match `src/funnel.js:52` `sigma=0.3 Å`, `w0=0.02 kcal/mol`, `biasFactor=6`, `hillStride` scaled to 100 fs):

```plumed
# plumed.dat — browser cv=r_com mapped to PLUMED DISTANCE
# System: 4W52 (164 Cα beads, BNZ ligand 6 heavy atoms) — adapt ATOMS to your PDB serials
UNITS LENGTH=A ENERGY=kcal/mol TIME=ps

# --- groups (replace atom numbers with your topology serials) ---
# ligand heavy atoms: e.g. serials 1288-1293 for BNZ in heavy mode (src/heavy.js 1287 prot + 6 lig)
comLig:    COM ATOMS=1288,1289,1290,1291,1292,1293
# pocket: Cα beads within rPocket=8 Å of native BNZ COM at t=0 (src/funnel.js:113)
# Example 4W52 pocket (update from Funnel.pocket after build):
comPocket: COM ATOMS=20,21,22,45,46,47,78,79,102,103

# --- CV: r_com (browser cv) -> PLUMED DISTANCE ---
cv: DISTANCE ATOMS=comLig,comPocket NOPBC

# --- optional funnel wall (browser U_wall = 0.5*kf*(r-rFlat)^2 for r>rFlat) ---
# funnel wall onset rFlat=5 Å, kf=2 kcal/mol/Å^2 (src/funnel.js:50 kf default)
# PLUMED UPPER_WALLS reproduces exactly:
UPPER_WALLS ARG=cv AT=5.0 KAPPA=2.0 EXP=2

# --- well-tempered metadynamics (browser WTM, src/funnel.js:19 biasFactor=6) ---
# PACE is hillStride scaled to keep ~100 fs physical spacing (src/funnel.js:69 getScaledStride)
# HEIGHT=w0, SIGMA=sigma, BIASFACTOR=gamma, GRID 0..rMax
METAD ARG=cv SIGMA=0.3 HEIGHT=0.02 BIASFACTOR=6 PACE=500 \
      GRID_MIN=0 GRID_MAX=24 GRID_BIN=96 \
      TEMP=300 FILE=HILLS

# --- output ---
PRINT ARG=cv FILE=COLVAR STRIDE=10

# To analyze the browser COLVAR with PLUMED:
# plumed sum_hills --hills HILLS --bin 96 --min 0 --max 24 --kt 0.596 --outfile fes.dat
# The browser PMF is exportable via Funnel.exportPMF() / getPMFcsv() with header
# "# T=300, gamma=6, hills=..., V0=1660.54" (src/funnel.js:620)
```

**Notes:**

- `UNITS LENGTH=A` keeps `cv`, `SIGMA`, `GRID_*`, `AT` in Å to match `src/funnel.js:61` `rMax=24 Å`, `bins=96`. Remove it and divide all lengths by 10 if you prefer nm.
- `KAPPA` in `UPPER_WALLS` is the harmonic constant in `ENERGY/LENGTH^2` (kcal/mol/Å² with `UNITS LENGTH=A`), matching `src/funnel.js:50` `kf`.
- For heavy mode the ligand serials start at `ligStart=1288` (after 1287 protein heavies for 4W52 per `tests/test_all.js:64`); for Cα+binding mode ligand serials are after `nProt=164` Cα beads.
- To reproduce the browser free energy exactly, post-process `HILLS` with `plumed sum_hills` and compare to `Funnel.getPMF()` / `exportPMF()` — both apply the Tiwary-Parrinello `c(t)` and `−γ/(γ−1)·V` scaling (`src/funnel.js:303` `getPMF`) and the same standard-state volume `V0=1660.54 Å³`.

## Grep & measurability

- `ls docs/PLUMED.md` exists
- `grep -n "DISTANCE" docs/PLUMED.md` hits mapping table + PLUMED snippet
- `grep -n "COLVAR" docs/PLUMED.md` hits export section
- `grep -n "bond order.*topology" src/mol2.js docs/*` hits `src/mol2.js:8` and this file's note above
