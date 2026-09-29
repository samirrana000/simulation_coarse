# Evolution 2 — S5 (short)

**Goal:** Cut repository bloat: 404MB data/coreset and 131MB .git
**Verdict:** ACCEPTED
**Note:** tracked 427MB -> 15.6MB (96.3% drop), 1712 files hash-verified restorable, all benches pass with data absent

**Acceptance:** A data/README.md documents exactly which files the app needs; a scripts/fetch_coreset.mjs regenerates the rest; runtime-needed files remain tracked and small; total tracked bytes drop by >80%; bench scripts that need the full set fail LOUDLY with a one-line instruction rather than silently mis-measuring; data/manifest.json hashes stay valid.

**Touched:** data/, scripts/, .gitignore, bench/
