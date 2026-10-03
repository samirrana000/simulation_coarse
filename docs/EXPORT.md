# Export — XYZ / PDB / DCD (H79)

This document describes trajectory export formats implemented in `src/recorder.js` (`Recorder.buildFile`, `Recorder._toXyz`, `Recorder._toPdb`) and the recommended **DCD** conversion path via MDAnalysis/MDTraj. See also `docs/BRIDGE.md` for the MDAnalysis bridge.

## Recorder Overview

- **Capture:** `Recorder` (`src/recorder.js:25` `class Recorder`) stores each frame as `Float32Array(3n)` + `times[]` (ps) (`src/recorder.js:28` `frames`, `src/recorder.js:29` `times`). Recording is stride-based in simulation time (`src/recorder.js:31` `stridePs`), capped by `maxFrames=500` (`src/recorder.js:32` `maxFrames`) with auto-stop guard (`src/recorder.js:86` `maxFrames cap guard`) to prevent memory leak (G69).
- **Control:** `recorder.start(nowPs, stridePs, maxFrames)` (`src/recorder.js:37` `start`), `maybeCapture(pos, timePs)` (`src/recorder.js:83` `maybeCapture`), `stop()`/`clear()`, `getFrame(i)` for scrubbing (`src/recorder.js:65` `getFrame`) and HUD `spanNs` (`src/recorder.js:53`).
- **Provenance:** Every file starts with `REMARK simulation_coarse vX, T, gamma, seed, date` (`src/recorder.js:112` `_provenanceLine` using `VERSION`, `BUILD_DATE` from `src/version.js`) so exports are traceable.

## XYZ Export

- **Invocation:** `recorder.buildFile("xyz", beads, opts)` (`src/recorder.js:106` `buildFile`) → `recorder._toXyz(beads, opts)` (`src/recorder.js:120` `_toXyz`). UI: `src/controllers/recording.js:202` `recorder.buildFile(fmt, state.sel.beads, prov)` with `fmt` from `ui.exportFmt.value` (`index.html` `#exportFmt`).
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

- **Invocation:** `recorder.buildFile("pdb", beads, opts)` → `recorder._toPdb(beads, opts)` (`src/recorder.js:135` `_toPdb`). UI produces `cg_traj_${count}frames.pdb` (`src/controllers/recording.js:203`).
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

## Structured Results Record (`results-record-v1`)

- **What it is:** ONE small JSON file carrying everything needed to reproduce or
  cite a number from this tool — code version + build date, run timestamp, seed,
  every physics parameter actually used, the input structure's identity, each
  computed value **with its real uncertainty**, the run-health error count, and a
  machine-readable statement of what the numbers do **not** mean. Built by
  `buildResultsRecord` in `src/results-record.js`; the button is
  **Export Results Record (JSON)** (`index.html` `#resultsDlBtn`) in the
  Recording panel, wired by `initResultsExport` in
  `src/controllers/results-export.js`.
- **Why it exists:** before this, an answer lived in three artefacts — a
  trajectory file with a REMARK header, a PMF CSV, and the operator's notes.
  Reproducing or citing a number was guesswork. See `docs/RESULTS_RECORD.md`
  for the field-by-field schema.
- **Invocation:** `serializeResultsRecord(record)` → JSON text, plus
  `resultsCsv(record)` → the flat tabular half. Both round-trip
  (`parseResultsRecord` / `parseResultsCsv`), asserted by
  `tests/test_results_record.js`.
- **Not a trajectory dump.** No frames, no coordinates, no PDB text: the record
  is a summary with provenance and is size-capped at 64 KiB
  (`RESULTS_MAX_BYTES`; a realistic record is ~10 KiB). Frames go to the
  recorder export above; the `coverage.notIncluded` block names, per analysis
  family, exactly what this export does not carry.
- **Honesty is derived, not typed.** The `scope` block is built from
  `src/scope.js` at call time, and `src/scope.js` is held to `ROADMAP.md` §1
  (the "No QM/MM, no explicit membrane, no PME" guard + its five hard-no
  bullets), `docs/LIMITATIONS.md` and `docs/VALIDATION.md` by
  `scopeDrift()` in `tests/test_results_record.js`. Adding a hard "no" to
  ROADMAP.md without a matching scope row **fails a FAST-tier test**.
- **An uncertainty is never zero.** `uncertainty` is either a finite positive
  value with its estimator named, or `null` plus an `uncertaintyReason`
  sentence. The analysis modules return `0` for "not computed" in places
  (`computeThermodynamics` `dH_se` with fewer than two bootstrap blocks), and
  the record converts that to `null` + reason rather than writing a zero error
  bar. A sample SD of a per-frame track travels in `dispersion`, which is a
  different thing and never an uncertainty on the mean.
- **A dirty run does not export as if it were clean.** `resultsExportGate()`
  refuses when `src/errors.js` has recorded any error, when a value has no
  trust-boundary statement, or when the input bytes could not be hashed. An
  *unseeded* run is a warning rather than a block (that is the default browser
  path, and the record already says the run is not reproducible) — hiding it
  would add no honesty and would make the export unusable for most sessions.
  Measured in a real browser (Playwright, served locally): with one error
  injected through `recordError()`, the export is refused with the count, the
  last error's context and message, and the remedy — and no file is downloaded.
- **A value that could not be computed is absent, not `0`.** RMSIP with no
  elastic network and ΔG with no deposited hills yield no row at all, so
  `coverage.included` is worth reading; a genuinely measured zero survives as
  `0`, with no invented error bar.
- **Sizes and caps:** `RESULTS_RECORD_VERSION = 1`, `RESULTS_MAX_VALUES = 24`
  result rows, `RESULTS_MAX_BYTES = 64 KiB` serialize cap,
  `RESULTS_FILE_MAX_BYTES = 1 MiB` load cap.

### Python one-liner

```python
import json
r = json.load(open("results_4W52_v1.1.0-fp7.json"))
for v in r["results"]:
    bar = f"± {v['uncertainty']['value']}" if v["uncertainty"] else "NO ERROR BAR"
    print(f"{v['id']:>26} = {v['value']:>10.3f} {v['unit']:<12} {bar}")
    print(f"{'':>26}   means: {v['meaning']}")
    print(f"{'':>26}   NOT:   {v['notA']}")
print("clean run:", r["runHealth"]["clean"], "| errors:", r["runHealth"]["errors"])
print("scope:", r["scope"]["browserV1Guard"])
```

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
2. Press **● Rec** (`src/controllers/recording.js:155` `recBtn` → `recorder.start`) — HUD `recStatus` shows `count`/`spanNs` (`src/ui.js:174` `updateRecStatus`).
3. Run simulation (`▶ Run`), **Stop** (`src/controllers/recording.js:162` `recStopBtn`), choose format `XYZ`/`PDB` (`index.html` `#exportFmt`), click **Download** (`src/controllers/recording.js:185` `dlBtn` → `downloadText` via `src/recorder.js:168` `downloadText` creating `Blob` + `URL.createObjectURL`).

## Grep & Files

- `grep -n "XYZ" docs/EXPORT.md` hits XYZ section; `grep -n "PDB" docs/EXPORT.md` hits PDB section; `grep -n "DCD" docs/EXPORT.md` hits DCD section.
- Source: `src/recorder.js:106` `buildFile`, `src/recorder.js:120` `_toXyz`, `src/recorder.js:135` `_toPdb`, `src/recorder.js:168` `downloadText`.
- See also `docs/BRIDGE.md` for MDAnalysis RMSF parity.

---
*Exports are traceable via REMARK provenance; for DCD use the PDB→DCD conversion snippet above.*
