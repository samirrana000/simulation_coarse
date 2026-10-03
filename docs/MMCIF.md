# PDBx/mmCIF Support Status (I81)

**Status: PDB is the ingest format; mmCIF is absent (out of scope).**

This browser simulator loads structures as **PDB** via `src/pdb.js:parseCa` and
`src/pdb.js:fetchPdb`. The PDB pathway is the only fully supported ingest
today (ATOM/HETATM fixed-column parsing, chain/residue filtering in
`src/main.js:buildSystem`). For reproducibility and widest tool
compatibility the exported trajectories are also PDB multi-MODEL or XYZ
(`src/recorder.js:buildFile`).

## mmCIF — absent (out of scope)

There is **no mmCIF support in this project**, and no queued work to add some. A native
`_atom_site` loop parser is absent *by scope* — the same convention as ONNX in
[SCORER.md](SCORER.md) and WebGL in [VIEWER.md](VIEWER.md), per ROADMAP.md §1.
`src/mmcif.js` exists only to detect a mmCIF block and refuse:

- **Detection:** checks for the mmCIF `data_` block header (`/^\s*data_/im` or
  `text.includes("data_")`) — this is the canonical PDBx/mmCIF entry prefix.
- **Delegation:** if the file also contains ATOM ` CA ` records, delegates to
  `parseCa(text)` so coarse-grained Cα ingest keeps working during migration.
- **Refusal:** otherwise throws `mmCIF is not supported by this simulator — PDB is
  the ingest format. Convert with pdb_extract or gemmi: …`. The wording is
  deliberate: *not supported*, not *not yet supported*.

```js
import { parseMMCIF } from "./src/mmcif.js";
try {
  const parsed = parseMMCIF(await file.text());
} catch (e) {
  alert(e.message); // "mmCIF is not supported by this simulator — PDB is the ingest format. …"
}
```

## Why PDB first

- PDB fixed-column ATOM/HETATM parsing is trivial in the browser, zero deps,
  and matches the `1ake.pdb .. 4w52.pdb` presets shipped in the repo root.
- mmCIF requires a full `_atom_site` loop parser (categories, quoted fields,
  multi-line values) and handling of `label_atom_id` vs `auth_atom_id` — a
  format project of its own, and out of scope for a coarse-grained browser tool.

## Workaround

Download the **PDB** format from RCSB (`https://files.rcsb.org/download/<ID>.pdb`)
or convert locally:

```bash
gemmi convert input.cif output.pdb
# or
pdb_extract -e input.cif -o output.pdb
```

Then load via **PDB file** input in the UI. PDB is the ingest format, full stop.

## Grep hits

- `src/mmcif.js:parseMMCIF` — `data_` detection + delegation/refusal
- `docs/MMCIF.md` — this file (PDB is primary, mmCIF is absent)
