# PDBx/mmCIF Support Status (I81)

**Status: PDB is primary; mmCIF is future.**

This browser simulator loads structures as **PDB** via `src/pdb.js:parseCa` and
`src/pdb.js:fetchPdb`. The PDB pathway is the only fully supported ingest
today (ATOM/HETATM fixed-column parsing, chain/residue filtering in
`src/main.js:buildSystem`). For reproducibility and widest tool
compatibility the exported trajectories are also PDB multi-MODEL or XYZ
(`src/recorder.js:buildFile`).

## mmCIF stub

`src/mmcif.js:parseMMCIF` is a **stub** added for forward compatibility:

- **Detection:** checks for the mmCIF `data_` block header (`/^\s*data_/im` or
  `text.includes("data_")`) — this is the canonical PDBx/mmCIF entry prefix.
- **Delegation:** if the file also contains ATOM ` CA ` records, delegates to
  `parseCa(text)` so coarse-grained Cα ingest keeps working during migration.
- **Guidance error:** otherwise throws `mmCIF not yet supported, use PDB` with
  conversion advice (`gemmi convert` / RCSB PDB download).

```js
import { parseMMCIF } from "./src/mmcif.js";
try {
  const parsed = parseMMCIF(await file.text());
} catch (e) {
  alert(e.message); // "mmCIF not yet supported, use PDB"
}
```

## Why PDB first

- PDB fixed-column ATOM/HETATM parsing is trivial in the browser, zero deps,
  and matches the `1ake.pdb .. 4w52.pdb` presets shipped in the repo root.
- mmCIF requires a full `_atom_site` loop parser (categories, quoted fields,
  multi-line values) and handling of `label_atom_id` vs `auth_atom_id` —
  not yet worth the bundle size for a coarse-grained demo.

## Migration plan (future)

1. Native `_atom_site` loop parser mapping
   `_atom_site.group_PDB == "ATOM"` + `_atom_site.label_atom_id == "CA"`
   → Cα beads (same shape as `parseCa` output).
2. Wire `fetchPdb` → RCSB `.../download/<id>.cif` fallback when `.pdb`
   returns 404 (many new entries are CIF-only).
3. Expand `tests/` to round-trip a CIF via `parseMMCIF` vs `parseCa` on
   PDB-converted fixtures and assert identical bead counts / chains.

## Workaround today

Download the **PDB** format from RCSB (`https://files.rcsb.org/download/<ID>.pdb`)
or convert locally:

```bash
gemmi convert input.cif output.pdb
# or
pdb_extract -e input.cif -o output.pdb
```

Then load via **PDB file** input in the UI (`src/main.js:loadStructure`) —
PDB remains the guaranteed path until native mmCIF lands.

## Grep hits

- `src/mmcif.js:parseMMCIF` — stub + `data_` detection + delegation/throw
- `docs/MMCIF.md` — this file (PDB is primary, mmCIF is future)
