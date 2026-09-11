# Data — Contact Prior Versioning (I86)

## 4W52_contacts.json header

`data/4W52_contacts.json` ships a **versioned contact prior** for the Elastic Network Model scaling `γ·(1+α·p_ij)` (see `src/ml-tier.js`, `ml/export_esm_contacts.py`).

Header fields:

```json
{
  "source": "heuristic",
  "chain": "A",
  "sequence": "MNIFEM...KNL",
  "nRes": 164,
  "meta": {
    "created": "2026-08-05T15:22:25.947030+00:00",
    "min_p": 0.05,
    "model_version": "heuristic-v1",
    "version": "1.0.0",
    "model": "heuristic",
    "contact_prior_version": "1.0.0-heuristic-v1"
  },
  "contacts": [[i, j, p_ij], ...]
}
```

- `source` — ML source (`heuristic` or ESM model e.g. `esm2_t33_650M_UR50D` when `ml/export_esm_contacts.py --model` is used)
- `model_version` — **contact prior version** (`heuristic-v1` for the checked-in file; ESM mode writes `esm2_t33_650M_UR50D-v1`)
- `version` — semantic file format version (`1.0.0`)
- `contact_prior_version` — combined `version-model_version` for provenance (`1.0.0-heuristic-v1`)
- `min_p` — threshold that filtered `contacts` (`0.05`)
- `created` — ISO8601 creation timestamp

## Versioning policy

- **heuristic-v1** — deterministic Kyte-Doolittle + Gaussian pseudo-prior (no ESM weights). Ships with repo; `ml/export_esm_contacts.py` produces it when `fair-esm` is absent (`[esm unavailable - using heuristic prior]`).
- **esm2_t6_8M_UR50D-v1**, **esm2_t33_650M_UR50D-v1**, etc. — real ESM-derived attention prior from `ml/export_esm_contacts.py --model <name>`.

Bumping `version` (e.g. `1.0.1`, `2.0.0`) signals a breaking change in `contacts` format; bumping `model_version` signals new weights/thresholds.

## Generating a new versioned prior

```bash
# heuristic-v1 (no ESM)
python3 ml/export_esm_contacts.py --pdb 4w52.pdb --out data/4W52_contacts.json --min-p 0.05

# ESM mode (requires fair-esm + torch)
pip install fair-esm torch
python3 ml/export_esm_contacts.py --pdb 4w52.pdb --out contacts.json --model esm2_t6_8M_UR50D --min-p 0.05
# -> source=esm2_t6_8M_UR50D, model_version=esm2_t6_8M_UR50D-v1
```

The exporter writes `meta.model_version` and `meta.version` automatically; older files without a `model_version` should be re-exported or manually annotated.

## Grep hits (measurable)

```bash
grep -n "version" data/4W52_contacts.json   # hits model_version, version, contact_prior_version
grep -n "version" data/README.md            # hits this doc
```
