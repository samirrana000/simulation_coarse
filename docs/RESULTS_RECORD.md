# Results Record — one file that makes a number citable

**The problem this closes.** A scientist runs a simulation in this browser app
and gets their answer spread across three artefacts: a trajectory file whose
header carries a `REMARK` provenance line, a PMF CSV, and their own notes about
what parameters they used. The parameters, the seed, the code version, the input
structure's identity, the computed values **and their uncertainty** are never in
one place. Reproducing or citing a number from this tool is guesswork. This
project's entire thesis is honest, reproducible, interactive biophysics, so an
un-reproducible result undermines it directly.

**What it is.** One small JSON file — `results-record-v1`, schema version 1 —
carrying provenance, input identity, every computed value with its real
uncertainty (or an explicit statement that it has none), the run's error count,
and a machine-readable statement of what those numbers do *not* mean. Built by
`src/results-record.js`; the button is **Export Results Record (JSON)**
(`#resultsDlBtn`) in the Recording panel, wired by `src/controllers/results-export.js`.

**It is not a trajectory dump.** No frames, no coordinates, no PDB text. Frames
go to the recorder export ([EXPORT.md](EXPORT.md)); this is the summary that says
what was computed, with what, and how much to trust it. Hard-capped at 64 KiB
(`RESULTS_MAX_BYTES`); a realistic record is ~10 KiB.

**It is not a substitute for rigorous FEP.** The record says so itself, in
machine-readable form, in every file it writes. See the honesty section below.

---

## Field provenance

| Field | Supplied by | Notes |
|---|---|---|
| `provenance.code.version` | `src/version.js` `VERSION` | equals `package.json` `version`; `scripts/check.sh` fails on drift |
| `provenance.code.buildDate` | `src/version.js` `BUILD_DATE` | ISO date of this transform |
| `generatedAt` | `new Date().toISOString()` at build time | the one wall-clock field; not asserted in tests |
| `provenance.run.seed` | `LangevinIntegrator.getSeed()` (`src/integrator.js`) | **`null` on the unseeded `Math.random` path** — never coerced to 0 |
| `provenance.run.seedSource` | derived from whether `seed` is null | names `SeededRNG` (mulberry32) or says the run is not reproducible |
| `provenance.run.simulatedPs` / `.steps` / `.dtPs` | `integ.time`, `integ.steps`, `integ.dt` | |
| `provenance.physics.modelMode` | `state.heavyMode` | `cg` or `heavy` |
| `provenance.physics.physicsLevel` | `ff.physicsLevel` (`src/forcefield.js`) | `L0` / `L1` / `L2` |
| `provenance.physics.tierFlags` | `ff.chargesOn`, `ff.hbMode`, `ff.weak`, live `state.bindLog` | the flags actually in force, not the selector value |
| `provenance.physics.temperatureK` / `.frictionPerPs` / `.beadMassDa` | `integ.T`, `integ.zeta`, the mass control | read from the integrator, not re-derived from sliders |
| `provenance.physics.dtPs` | `integ.dt` (auto-tuned per system size) | |
| `provenance.physics.cutoffAng` / `.springGamma` | `ff.rc`, `ff.gamma` | ENM cutoff and spring force constant |
| `provenance.physics.solventModel` / `.saltM` / `.epsIn` / `.epsOut` / `.sasaGamma` | `settingsState` (`src/settings-panel.js`) | |
| `provenance.physics.respa` | `settingsState.respaOn` / `respaOuterFs` | |
| `provenance.physics.system` | `ff.n`, `ff.nProt`, `ff.nLigAtoms`, `state.ligands.length` | system composition |
| `provenance.physics.electrostatics` | constant naming the 6.5→8.5 Å truncation and **NO PME** | so a reader cannot assume a lattice sum |
| `provenance.notAPhysicsParameter.motionGain` | the `#motionGain` slider | filed **separately on purpose**: it is a Canvas2D display amplification (`viewer.setMotionGain`), not a force. Listing it as a physics parameter would be a lie about the run. |
| `input.pdbId` / `.fileName` / `.origin` | `state.inputSource`, set in `loadStructure` (`src/controllers/structure-input.js`) | `pdb-id-fetch` or `file-upload` |
| `input.contentHash` | `sha256Text(state.pdbText)` (`src/sha256.js`) | `sha256:<64 hex>` — **the same convention and the same digest as `data/manifest.json`**, so a record and the manifest are comparable with `===`. Verified against the manifest's own value for the bundled `4w52.pdb`. |
| `input.selection` | the chain / residue-range / ligand controls | the window actually simulated |
| `trajectory.*` | `recorder.count`, `recorder.times`, `recorder.stridePs` | |
| `results[]` | `resultSpecs()` in the controller, from `analyzeTrajectory`, `Funnel.estimateDG`/`convergenceSE`, `computeThermodynamics`, `detectCryptic` | see below |
| `runHealth.errors` / `.documentedNoOps` / `.lastError` | `errorSummary()` from `src/errors.js` | `clean` is **derived** (`errors === 0`); a caller cannot pass `clean: true` for a dirty run |
| `scope` | `scopeBlock()` from `src/scope.js`, at call time | see the honesty section |

---

## Results: value + uncertainty

Each `results[]` row:

```json
{
  "id": "dg_bind",
  "label": "binding dG from the WTM funnel bias (128 hills, 6 bias factor)",
  "value": -3.41,
  "unit": "kcal/mol",
  "uncertainty": {
    "value": 0.42,
    "kind": "hill-count-convergence",
    "method": "kB*T/sqrt(nHills)",
    "n": 128
  },
  "uncertaintyAvailable": true,
  "uncertaintyReason": null,
  "n": 128,
  "meaning": "a ranking indicator between systems simulated with identical settings, seed protocol and hill budget",
  "notA": "an absolute binding free energy, a ΔG° prediction, or a K_D",
  "notComputed": ["alchemical intermediates or soft cores", "replica exchange", "…"],
  "substitute": "alchemical FEP/TI or umbrella sampling in GROMACS/AMBER/OpenMM + PLUMED",
  "documentedIn": "src/scope.js VALUE_TRUST_BOUNDARY",
  "honestyMissing": false
}
```

### The uncertainty rule

`uncertainty` is **either a finite positive value with its estimator named, or
`null`**. There is no third option, and in particular there is no `0`.

This matters because the analysis modules already return `0` for *not computed*:

- `computeThermodynamics` leaves `dH_se` at `0` when its 20-block bootstrap has
  fewer than two populated blocks (`src/analysis/thermodynamics.js`);
- `p.sasa.se ?? 0` does the same for the real-burial term;
- `Funnel.convergenceSE()` returns `null` below 50 deposited hills.

`normalizeUncertainty()` in `src/results-record.js` maps any of those to
`uncertainty: null` plus an `uncertaintyReason` sentence that says *why* it is
unknown. A zero error bar is a lie: it claims a precision nobody measured.

**What a value with no uncertainty looks like:**

```json
{
  "id": "occupancy_bound_fraction",
  "value": 0.86,
  "unit": "fraction of frames",
  "uncertainty": null,
  "uncertaintyAvailable": false,
  "uncertaintyReason": "no binomial or bootstrapped error bar is computed on the bound-frame fraction",
  "n": 51,
  "honestyMissing": false
}
```

The field is present, explicitly `null`, with a reason — not omitted, and not
zero. The export gate reports which ids have no error bar in its warnings, so
the omission is visible at export time too.

### Dispersions are not uncertainties

`detectCryptic` returns a sample standard deviation of the per-frame
pocket-volume track. That is a **dispersion of the sample**, not an error bar on
the mean, and conflating the two is the classic false precision. It travels in
its own field:

```json
"dispersion": {
  "value": 88.1,
  "kind": "sample-sd-of-track",
  "note": "a dispersion of the underlying sample, NOT an uncertainty on the mean above"
}
```

with `uncertainty: null` on the same row.

### A value that could not be computed is absent, not zero

The same rule applies to the value itself. `resultRow()` drops a row whose value
is not finite, so RMSIP with no elastic network (heavy mode) or ΔG with no
deposited hills (`estimateDG()` → NaN) never appear as `0 kcal/mol`. A zero
would be a measurement; absence is the truth. Read `coverage.included` to see
which rows the record actually holds — and a *genuine* measured zero (a pocket
volume of 0 Å³, say) does survive as `0`, with no invented error bar.

### Error bars actually used (all from existing analysis code, none invented)

| Quantity | Estimator | Source |
|---|---|---|
| `dg_bind` | `kB*T/sqrt(nHills)` hill-count convergence bound | `Funnel.convergenceSE` (`src/funnel.js`) |
| `dH_total` | 20-block standard error | `computeThermodynamics` (`src/analysis/thermodynamics.js`) |
| `dS_total` | B-block jackknife (when χ torsions supplied) | `chiEntropyDelta` |
| `dG_jarzynski` | seeded nonparametric bootstrap | `jarzynskiFreeEnergy` (`src/analysis/unbinding_smd.js`) |

---

## Honesty: derived, not typed

The `scope` block is **built from `src/scope.js` at call time** — a record can
never disagree with the table it ships with. `src/scope.js` in turn is held to
the documents by `scopeDrift()`, which `tests/test_results_record.js` runs
against the real `ROADMAP.md`, `docs/LIMITATIONS.md` and `docs/VALIDATION.md`:

- the `> **Browser v1 will NOT do:**` guard sentence is **parsed verbatim** out
  of ROADMAP.md and compared to `BROWSER_V1_GUARD` (the line the project's own
  documented grep, `ROADMAP.md:42`, hits);
- ROADMAP §1 is parsed for its `- **No …**` bullets (scoped to §1, so §2's "what
  we WILL do" list cannot leak in), and **every** bullet must have a matching
  `OUT_OF_SCOPE` row — and every row must still match a bullet;
- each `limitsTerms` phrase must still appear in `docs/LIMITATIONS.md`;
- `docs/VALIDATION.md` must still contain `ranking-only` and `NOT FEP`.

**The drift test proves it can fail.** It feeds `scopeDrift()` mutated documents
and asserts a problem is reported for each of: a new §1 bullet with no table row,
a rewritten guard sentence, a deleted §1 bullet, a `LIMITATIONS.md` edit that
drops a cited statement, a `VALIDATION.md` that loses the ranking-only verdict,
and an unparsable ROADMAP (which must not silently read as "no drift").

**What is *not* derived at runtime, stated plainly.** The browser cannot read
`ROADMAP.md` at export time: there is no build step to inline it, and a `file://`
origin cannot read it at all. So the table ships as code and a failing test is
what catches an edit to the document — the cheap direction to fail in. The record
says this itself, in `scope.derivation.runtimeDocRead: false` plus the reason, so
a downstream reader is never misled about how the statement was obtained.

### Per-value trust boundary

Every exported value carries `meaning`, `notA`, `notComputed`, and `substitute`.
`dg_bind`, for example, states that it is *a ranking indicator between systems
run with identical settings*, that it is **not** an absolute binding free energy
or a `K_D`, that no alchemical intermediates / soft cores / replica exchange
exist behind it, and that rigorous FEP/TI in GROMACS/AMBER/OpenMM + PLUMED is
the substitute.

**No value exports without one.** `resultsExportGate()` fails closed when a
result id has no `VALUE_TRUST_BOUNDARY` entry (`honestyMissing: true`). "We
forgot to say what ΔG means" cannot ship as a clean-looking JSON file.

---

## Run health

`runHealth` carries `errors`, `documentedNoOps`, `lastError` and a **derived**
`clean` (`errors === 0`), sourced from the `src/errors.js` recorder — the same
counter shown on the top status bar and `#topErrors`.

A run with recorded errors is **not exportable as if it were clean**:
`resultsExportGate()` refuses, naming the count, the last error and the remedy.
The counter is cumulative from page load (it is not per-run), so the remedy says
reload rather than pretending a fresh run would clear it.

Documented no-ops (`errors.js` `ignore()` — the defensive catches) are counted
and carried too, but do not block: they are not failures, and hiding them would
be as dishonest as hiding real errors. They appear as an export **warning** and
as `runHealth.documentedNoOps`.

An **unseeded** run is a warning, not a block. Blocking it would make the export
unusable for every default browser session (the default path is
`Math.random`) while adding no honesty, because `provenance.run.seedSource`
already states that the run is not reproducible.

---

## Formats and round-trip

**JSON** — `serializeResultsRecord(record)` → pretty JSON, newline-terminated,
size-capped. `parseResultsRecord(text)` → `{ok, data, error}`, rejecting a wrong
version, a missing `results` array, malformed JSON, an empty file, or a file
over `RESULTS_FILE_MAX_BYTES` (1 MiB, `SYSTEM_TOO_LARGE`).

**CSV** — `resultsCsv(record)` → the flat tabular half, the columns a
downstream script joins on:

```
id,label,value,unit,uncertainty,uncertainty_kind,n,uncertainty_reason,meaning,not_a
```

`parseResultsCsv(text)` reads it back. A cell with no error bar is an **empty
field**, which parses back to `null` — so "no error bar" survives the round trip
as "no error bar". A consumer that writes `Number(cell) || 0` gets a zero only by
its own choice, which is the correct place for that decision to live.

Both round-trips are asserted in `tests/test_results_record.js` (FAST tier),
following the conventions of `tests/test_session_roundtrip.js`: a real
`serialize → parse → assert every value` cycle for provenance, input identity,
every result row (value, unit, uncertainty object, null-reason, meaning, notA)
and the scope block.

---

## Grep & files

- `grep -n "results-record-v1" docs/RESULTS_RECORD.md src/results-record.js` — schema id.
- Source: `src/results-record.js` (`buildResultsRecord`, `serializeResultsRecord`,
  `parseResultsRecord`, `resultsCsv`, `parseResultsCsv`, `resultsExportGate`,
  `resultsFilename`), `src/scope.js` (the table + the ROADMAP parser +
  `scopeDrift`), `src/sha256.js` (the content hash),
  `src/controllers/results-export.js` (the button wiring).
- Tests: `tests/test_results_record.js`, registered in `tests/suites.js` FAST tier.
- See also [EXPORT.md](EXPORT.md) for trajectory formats,
  [LIMITATIONS.md](LIMITATIONS.md) for the per-value trust boundary this record
  quotes, and [VALIDATION.md](VALIDATION.md) for the measured numbers and the
  ranking-only verdict. The failure counters come from `src/errors.js` and are
  shown on the top status bar (`#topErrors`).

---

*Every record states what it does not compute. A number from this tool is a
ranking indicator until an alchemical calculation or an experiment says
otherwise.*
