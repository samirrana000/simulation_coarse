# CHANGELOG — simulation_coarse

Release history with commit hashes, gate status, and test tiers.
Trust boundary (all releases): ranking-only use — NOT FEP, NOT absolute Kd.
Test tiers: `node tests/test_all.js` (FAST, gate) · `--slow` (SLOW-SMOKE) ·
`--slow-full` (SLOW-FULL). Version flow: `src/version.js` (no `package.json`).

## 1.1.0-fp7 (2026-09-13) — release closeout [FP7]

- Version: `1.0.0-transform` → `1.1.0-fp7` (`src/version.js`; minor bump —
  Loop-2 + FP1–FP7 were backward-compatible additive-only, no breaking
  changes; patch would understate seven finish-product increments).
  Footer version line in `index.html` dock (text only, no new DOM ids);
  HUD prefix + trajectory `REMARK` provenance follow automatically.
- License/citation: `LICENSE` (MIT, © 2026 Samir Rana) present; `CITATION.cff`
  (`v1.0-jpcb`, Zenodo DOI placeholder) present; README gains
  `## Version, license & citation`; secret grep over
  `src scripts tests ml index.html CITATION.cff` → zero hits.
- Gate BEFORE (working tree, pre-FP7): OPEN — 45 files clean, 352 PASSED /
  0 FAILED (≥ 352 baseline), 112 ui ids, Digit1-7 vs 8 panels.
- Gate AFTER: **OPEN — 45 files clean, 352 PASSED / 0 FAILED, 112 ids**
  (measured 2026-09-13; footer adds zero DOM ids).
- FAST `node tests/test_all.js` → **352/352 in ~13.6 s** (before) /
  **352/352 in ~13.4 s** (after, bit-identical tier).
- SLOW-SMOKE `node tests/test_all.js --slow` → **404/404 in ~29.5 s**
  (FAST 352 + SLOW 52: thermo 7 + heavy SMOKE 13 + calibration 14 +
  flexlig 10 + 1crn-null 8).
- `node --check` clean on all touched files (version.js display-only value
  bump; index.html parses); serve smoke: `/`, `src/main.js`,
  `src/version.js`, `4w52.pdb` all 200, footer version string present.
- Full record: `docs/BINDING_LOOP2_DONE.md` §28; wiki
  `.wikiskill/wiki/evolution/skill-impact.md` Proposal FP7 (ACCEPTED).

## FP6 — methods/validation table + 1CRN negative control (2026-09-13, working tree)

- NEW `docs/VALIDATION.md`: trust boundary with numbers (ΔG_est underbinds
  2–3 kcal/mol, bootstrap ±0.16 ⊕ replica-SD ±7.1) + methods table
  (4W52-BNZ −6.82 ± 0.16, BNZ-only −3.64 ± 0.08, ΔSASA 167.1 ± 3.5 →
  solvent +2.01 ± 0.04 → ΔG −0.16 vs ITC −5.20 ± 0.20 / NMR −4.20 ± 0.10,
  EPE ΔS_lig 0.00305 vs BNZ exactly 0, L0:L2 ≈ 400:1).
- NEW `scripts/validate_1crn_null.mjs` (8 asserts, SLOW, <1 s): 1CRN null
  ΔH 0.00 ± 0.00 (|ΔH| < 0.5 vs 4W52 −3.6/−6.8), −TΔS +5.67 bounded < 20.
- Tests: FAST stays 352; SLOW 44 → 52, `--slow`/`--slow-full` 396 → 404.

## FP5 — heavy-build progress UX + CG/heavy workflow (2026-09-13, working tree)

- Chunked heavy topology (`buildTopologyChunked`, 128 rows/slice, 11 slices
  on 4W52) with captions + `heavyCancelBtn` (+1 id → 112); CG path untouched.
- `docs/WORKFLOWS.md` CG-interactive/heavy-offline section (1 ns heavy ≈
  1 day → headless only). Tests: FAST stays 352; FP5 proof 18/18.

## FP4 — export matrix + session save/load (2026-09-13, working tree)

- One-click exports (traj XYZ/PDB/JSON, PMF CSV, thermo TXT, BindLog BLG1,
  DCCM CSV, settings JSON, session JSON v1 ≤ 256 KiB, frames never embedded).
- NEW `src/session.js` + `tests/test_session_roundtrip.js` (46 asserts).
- Tests: FAST 306 → **352** (320 FAST + Tier-0 32); gate baseline → 352;
  111 ids.

## FP3 — flake root-cause + browsers + a11y (2026-09-13, working tree)

- Seeded bindlog-integration (seed 101, 17/30 negatives, 5× 14/14) + seeded
  temperature mean (303.7 K); Chromium + Firefox smoke PASS, Safari
  static-only; `:focus-visible` rings, canvas `role="img"`, Space guard.
- Tests: FAST stays 306; 105 ids.

## FP2 — input robustness (2026-09-13, working tree)

- NEW `src/input_errors.js`: 9 frozen failure classes, every message carries
  the next click; `tests/test_input_errors.js` (75 asserts).
- Tests: FAST 231 → **306**; gate baseline → 306; 105 ids.

## FP1 — first-run UX (2026-09-13, working tree)

- 1-click 4W52 sample, 4-step checklist (Load → Build → Run → Analyze),
  empty-state audit; +5 ids → 105. Tests: FAST stays 231.

## followup1-7 (2026-09-13) — `173ccad`

- Deep entropy LONG (10 ps × 3: full +6.76 SD 14.06, χ +0.24 ± 0.33 —
  variance-dominated, NOT-RESTORED stands); thermo ligand picker (EPE
  inflates |ΔH| by −3.18, reported not corrected); real LCPO burial
  (ΔSASA 167.1 ± 3.5 → −TΔS_solv +2.01 ± 0.04, sign review open);
  ala-scan noise floor 0.05 (all `~noise`, NO CG resolution); EPE flexlig
  end-to-end (ΔS_lig 0.00305 vs BNZ 0); heavy SMOKE/FULL split
  (`--slow` ~36 s vs `--slow-full` ~230 s); thermo Cancel + PMF 0-hill hint.
- Tests: FAST 215 → 231 (ala floor +16); SLOW 34 → 44; slow 275; 99 → 100 ids.

## stages1-7 (2026-09-12) — `b3cf17d`

- Heavy pocket-ΔS verdict NOT-RESTORED; seeded Langevin (mulberry32,
  default path bit-identical); tiered `test_all` (FAST 179); altloc cleaner +
  rotbonds (FAST → 215); async thermo + physics-level persist; 4W52
  calibration anchor (ΔG_est −2.23 vs ITC −5.20/NMR −4.20, ranking-only);
  alloc-tracked bench (gc deltas ≤ 0.3 MB, GPU-port declined, CPU ACCEPT).
- Tests: FAST 32 → 179 → **215**; SLOW 199 → 245 → 249.

## binding-loop2 (2026-09-12) — `2c29d25`

- S1–S7: CG charges, directional-HB virtual sites, heavy weakint, BindLog
  (CG 4-term/heavy 7-term), thermo ΔH/ΔS panel, BindViz, L0/L1/L2 selector
  (default L0 bit-identical). ΔH −6.9 ± 0.16; CG cost ×1.09; Loop-2 ≈1.0×
  on heavy. Tests: `test_all` **32/32** + thermo 7 + bindviz 22 + bindlog 14.

## phases1-5 (2026-09-12) — `3ee4914`

- Physics rigor (AMBER-lite/Tirion-GB/LCPO/membrane), chemistry
  (protonation/GAFF2-lite/stereo/metals), GPU/RESPA, workflows
  (ala-scan/DCCM/cryptic/SMD), cockpit UI. Docs + README.

## items 2–4 (2026-09-11/12) — `4eeaffc` `de10a73` `57eb2e8` `d21becf` `b76c019`

- Ligand library + clash-free placement; all-atom heavy mode (covalent
  topology + metal coordination); MOL2/heavy parity; sidebar restyle.
  `7de5505` `ae1961d` bootstrap the CG app + force-field modularization.

## Test-count lineage (FAST gate baseline)

`32` (Loop-2 S7) → `179` (tiered) → `215` (altloc+rotbonds) →
`231` (ala floor) → `306` (FP2 input_errors +75) → `352` (FP4 session +46).
SLOW: `199` → `245` → `249` → `275` (FAST 231 + SLOW 44) →
`404` (FAST 352 + SLOW 52: thermo 7 + heavy 13 + calibration 14 +
flexlig 10 + 1crn-null 8).
