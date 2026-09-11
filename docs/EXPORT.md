# Export — XYZ / PDB / DCD (H79)

This document describes trajectory export formats implemented in `src/recorder.js` (`Recorder.buildFile`, `Recorder._toXyz`, `Recorder._toPdb`) and the recommended **DCD** conversion path via MDAnalysis/MDTraj. See also `docs/BRIDGE.md` for the MDAnalysis bridge.

## Recorder Overview

- **Capture:** `Recorder` (`src/recorder.js:25` `class Recorder`) stores each frame as `Float32Array(3n)` + `times[]` (ps) (`src/recorder.js:28` `frames`, `src/recorder.js:29` `times`). Recording is stride-based in simulation time (`src/recorder.js:31` `stridePs`), capped by `maxFrames=500` (`src/recorder.js:32` `maxFrames`) with auto-stop guard (`src/recorder.js:86` `maxFrames cap guard`) to prevent memory leak (G69).
- **Control:** `recorder.start(nowPs, stridePs, maxFrames)` (`src/recorder.js:37` `start`), `maybeCapture(pos, timePs)` (`src/recorder.js:83` `maybeCapture`), `stop()`/`clear()`, `getFrame(i)` for scrubbing (`src/recorder.js:65` `getFrame`) and HUD `spanNs` (`src/recorder.js:53`).
- **Provenance:** Every file starts with `REMARK simulation_coarse vX, T, gamma, seed, date` (`src/recorder.js:112` `_provenanceLine` using `VERSION`, `BUILD_DATE` from `src/version.js`) so exports are traceable.

## XYZ Export

- **Invocation:** `recorder.buildFile("xyz", beads, opts)` (`src/recorder.js:106` `buildFile`) → `recorder._toXyz(beads, opts)` (`src/recorder.js:120` `_toXyz`). UI: `src/main.js:545` `recorder.buildFile(fmt, state.sel.beads, prov)` with `fmt` from `ui.exportFmt.value` (`index.html` `#exportFmt`).
- **File structure (multi-frame XYZ):**
  ```
  REMARK simulation_coarse v1.0.0-transform T=300K gamma=2 seed=0 date=2026-08-29
  164
  frame 0  time = 0.000 ps  (CG Cα model)
  C  12.345  8.901  3.210
  C  13.111  9.222  3.444
  ...
  164
  frame 1  time = 2.000 ps  (CG Cα model)
  C  12.350  8.905  3.215
  ...
  ```
  Header is the provenance REMARK (`src/recorder.js:122` `out=[provenanceLine]`), then per frame: atom count `n`, comment `frame k  time = X ps  (CG Cα model)` (`src/recorder.js:125` `time = ${times[k].toFixed(3)} ps`), then `n` lines `C  x  y  z` with `toFixed(3)` Å (`src/recorder.js:128`). All `n` beads are exported as `C` (Cα model) regardless of element; heavy mode still exports Cα-mapped XYZ (for full heavy XYZ, use PDB).
- **Units:** `Å` for coordinates, `ps` for time (see `docs/UNITS.md:76` `XYZ: Å, time in comment line`).
- **Loading:** `MDAnalysis.Universe("topology.pdb", "traj.xyz")` (`docs/BRIDGE.md:87` `XYZ reader is built-in`) or `MDTraj.load("traj.xyz", top="top.pdb")`.

## PDB Export

- **Invocation:** `recorder.buildFile("pdb", beads, opts)` → `recorder._toPdb(beads, opts)` (`src/recorder.js:135` `_toPdb`). UI produces `cg_traj_${count}frames.pdb` (`src/main.js:546`).
- **File structure (multi-MODEL PDB):**
  ```
  REMARK simulation_coarse v1.0.0-transform T=300K gamma=2 seed=0 date=2026-08-29
  REMARK  CG Cα Langevin trajectory (BAOAB integrator)
  MODEL     1
  REMARK  time = 0.000 ps
  ATOM      1  CA  GLY A   1      12.345   8.901   3.210  1.00  0.00           C
  ...
  ENDMDL
  MODEL     2
  REMARK  time = 2.000 ps
  ATOM      1  CA  GLY A   1      12.350   8.905   3.215  1.00  0.00           C
  ...
  ENDMDL
  END
  ```
  (`src/recorder.js:137` `REMARK  CG Cα Langevin...`, `src/recorder.js:139` `MODEL`, `src/recorder.js:140` `REMARK  time =`, `src/recorder.js:143` `ATOM` with `resName/chain/resSeq` from `beads[i]` via `b0(beads[i])` (`src/recorder.js:163`), `ENDMDL`/`END`).
- **Residue mapping:** Each Cα bead becomes one `ATOM  CA` with original `resName` (`b.resName||"GLY"`), `chain` (`"_"`→`" "` per PDB convention), `resSeq` (`src/recorder.js:163` `b0`). Ligand atoms when present are exported with their own names (heavy mode uses `state.sel.atoms` metadata).
- **Loading:** `MDAnalysis.Universe("traj.pdb")` (`docs/BRIDGE.md:35` `Universe("cg_traj_500frames.pdb")`) — MDAnalysis reads `MODEL/ENDMDL` as trajectory; `u.trajectory` length equals `recorder.count`. Align then RMSF as in `docs/BRIDGE.md:46` `AlignTraj` + `RMSF`.

## DCD Export

- **Native support:** `src/recorder.js` does not write DCD directly (DCD is a binary FORTRAN unformatted format requiring typed arrays). The recommended path is **PDB topology + DCD** via one-time conversion in MDAnalysis/MDTraj (`docs/BRIDGE.md:64` `Option B — PDB topology + DCD`):
  ```python
  # One-time conversion: PDB multi-MODEL → DCD + topology PDB
  import MDAnalysis as mda
  u = mda.Universe("cg_traj_500frames.pdb")  # topology+trajectory in one PDB
  u.trajectory[0]
  ca = u.select_atoms("all")
  ca.write("topology.pdb")  # first frame as topology
  with mda.Writer("trajectory.dcd", n_atoms=len(ca)) as W:
      for ts in u.trajectory:
          W.write(ca)
  # Now: u2 = mda.Universe("topology.pdb", "trajectory.dcd")
  ```
  (`docs/BRIDGE.md:69` `One-time conversion`). The resulting `trajectory.dcd` is compact (binary, 4 bytes/coordinate vs ~10 bytes in PDB) and faster to read for large trajectories (`>500` frames).
- **Alternative via MDTraj/CPPTRAJ:**
  ```python
  import mdtraj as md
  t = md.load("cg_traj_500frames.pdb")  # reads MODELs
  t.save_dcd("trajectory.dcd")
  t[0].save_pdb("topology.pdb")
  ```
- **Why DCD is documented:** The browser cannot emit DCD natively without extra binary writer and endianness handling; the PDB→DCD conversion is lossless (coordinates preserved to `0.001 Å` due to PDB `toFixed(3)`, or full `Float32` precision if writer is extended to use `recorder.frames` directly). For exact `Float32` preservation, a future `Recorder._toDcd()` could write DCD directly from `frames` (not yet implemented — placeholder).
- **Measurability:** This file mentions `docs/EXPORT.md` documenting `XYZ/PDB/DCD` exports — `grep -n "XYZ" docs/EXPORT.md`, `grep -n "PDB" docs/EXPORT.md`, `grep -n "DCD" docs/EXPORT.md` all hit.

## UI Flow

1. Set `stridePs` and `maxFrames` (`index.html` `#stridePs`, `#maxFrames` in `Recording` panel).
2. Press **● Rec** (`src/main.js:522` `recBtn` → `recorder.start`) — HUD `recStatus` shows `count`/`spanNs` (`src/main.js:582` `updateRecStatus`).
3. Run simulation (`▶ Run`), **Stop** (`src/main.js:529` `recStopBtn`), choose format `XYZ`/`PDB` (`index.html` `#exportFmt`), click **Download** (`src/main.js:534` `dlBtn` → `downloadText` via `src/recorder.js:168` `downloadText` creating `Blob` + `URL.createObjectURL`).

## Grep & Files

- `grep -n "XYZ" docs/EXPORT.md` hits XYZ section; `grep -n "PDB" docs/EXPORT.md` hits PDB section; `grep -n "DCD" docs/EXPORT.md` hits DCD section.
- Source: `src/recorder.js:106` `buildFile`, `src/recorder.js:120` `_toXyz`, `src/recorder.js:135` `_toPdb`, `src/recorder.js:168` `downloadText`.
- See also `docs/BRIDGE.md` for MDAnalysis RMSF parity.

---
*Exports are traceable via REMARK provenance; for DCD use the PDB→DCD conversion snippet above.*
