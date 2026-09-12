# Loop-2 DONE — Binding-Physics Integration (S1–S7)

All 7 steps landed, additive-only, zero new deps, vanilla ES modules, JSDoc,
units Å/ps/kcal/mol/Da. Gate OPEN before and after; full suite green.

## 1. S1–S7 table

| # | Step | What shipped | Files | Validation |
|---|---|---|---|---|
| S1 | CG salt-bridge charges (R2 term b) | `CG_FORMAL_CHARGES` (ASP/GLU −1, LYS/ARG +1, HIS 0); opt-in `_protQ` override, default OFF | `src/ff-params.js`, `src/forcefield.js` | charged-ligand ΔU matches R2 screened profile |
| S2 | CG virtual sites + directional HB + 52%-flag fix (R2 term a) | Cα-triplet virtual O-sites for ALL residues; angular gates; flag = site-existence | `src/forcefield.js`, `src/ff-binding.js` | natives ≥ −1.6, decoys 0.000; CG cost ×1.09 (S7) |
| S3 | Heavy π-stack + cation-π + halogen σ-hole (R3 a–c) | NEW `src/physics/weakint.js`; `par.weak` opt-in, default OFF; aromatic LJ-ε ×0.70 + cation charge gate | `src/physics/weakint.js`, `src/heavy.js` | R3 §5 prototype numbers; FD < 1e-6; cost ≈1.0× (S7) |
| S4 | Per-term accumulators → BindLog wiring (R4 §5.1, R6 §5) | `bindU` vector (CG 4-term, heavy 7-term); tick frames + energy events + contact diff + funnel hills; `bindlogOn` default OFF | `src/ff-binding.js`, `src/heavy.js`, `src/main.js` | `test_bindlog_integration` 14/14; OFF bit-identical; CG ratio 1.03× |
| S5 | Thermodynamics ΔH/ΔS panel (R4 pipeline) | `src/analysis/thermodynamics.js` (Schlitter Cholesky, torsion Shannon, Kabsch, bootstrap, SASA proxy) + panel section | `src/analysis/thermodynamics.js`, `src/analysis-panel.js` | `test_thermo` 7/7; ΔH −6.9±0.16 (see §3) |
| S6 | BindViz Binding Insights subpanel + D1 fix (R7 §4) | Collapsed `<details>` in PMF & Analysis; `bindvizTimeline/Energy/Pmf` ≤1 Hz; null-bindlog guard | `index.html`, `src/main.js`, `src/capture/bindviz.js`, `src/ui.js` | `test_bindviz` 22/22; headless pixel-diversity verified |
| S7 | Pareto re-bench (FULL) + physics-level selector + docs (this step) | `pareto_bench.mjs` measures L0/L1/L2-full; `#physicsLevel` selector; cation-π NaN guard; this doc | `scripts/pareto_bench.mjs`, `index.html`, `src/main.js`, `src/ui.js`, `src/settings-panel.js`, `src/physics/weakint.js` | gate OPEN (98 ids); 32/32, 7/7, 22/22, 14/14; selector check 24/24 |

## 2. S7 re-bench (FULL mode, 4W52: 164 Cα + 21 lig atoms; 1308 heavy atoms)

`node scripts/pareto_bench.mjs` → `docs/pareto_frontier.csv` (canonical rows
kept; L1/L3 promoted from Loop-1 estimates to measurements; L1+BindLog and
L2-full-rigor rows added; prior estimates in git history).

| Tier | ms/step | heap | top-1 | mean rank | Terms |
|---|---|---|---|---|---|
| L0 CG | 0.227 | 7.4 MB | 10/10 | 1.0 | ENM + isotropic-LJ + burial |
| L1 CG+ | 0.247 (×1.09) | 9.2 MB | 10/10 | 1.0 | + charges + directional-HB v-sites (162 valid, S1+S2 measured) |
| L1 CG+ +BindLog | 0.262 (×1.15 vs L0) | 10.8 MB | — | — | + trackTerms accumulators, bindU finite (S4 measured) |
| L2 heavy | 99.521 | 12.1 MB | 9/10 | 1.3 | covalent + LJ + GB |
| L3 heavy+R3 | 89.490 (≈1.0×, noise) | 24.2 MB | 9/10 | 1.3 | + weakint, 19 rings, weakU finite (S3 measured) |
| L2 full-rigor | 86.829 (≈1.0×, noise) | 17.6 MB | 9/10 | 1.3 | heavy + weakint + 7-term bindU (S7 measured) |
| L4 heavy+OBC2+RESPA (est) | 54.736 | — | est | est | OBC2 + LCPO + RESPA 2 fs (unchanged by Loop 2) |

Speed ratio L0:L2 ≈ **439:1** (≈400:1 in Loop 1 — reproduced). Heap column is
GC-snapshot ±50% (review D4); memory is NOT a frontier axis at this size.
Pose-recovery caveat stands (R5 §1 footnote): the proxy scorer measures burial
discrimination, not shape complementarity; L1/L3 differentiation needs full-FF
ΔG scoring (Loop 3).

S7 bench hardening (additive, 2 lines): `cationPiEnergy` divided by an
unguarded `dist` — coincident cation/centroid (overlapping duplicate-altloc
atoms in the raw crystal; LJ already reports 1e155 there) yielded NaN that
poisoned total U. Added the `dist < 1e-8 → zero-cut` guard in the existing
halogen-guard style (`src/physics/weakint.js:185`). S3 prototype tests
(dist ≥ 2 Å) and `test_all` (weak default OFF) unaffected; heavy L0/L2
energies now finite headless.

## 3. Carried S5 numbers (unchanged, for the Loop-2 record)

4W52 CG, charges + directional-HB on, 2000 steps × 3 replicas:
**ΔH = −6.9 ± 0.16 kcal/mol** (LJ −1.9 / Coul +0.5 / HB −0.1 / desolv −5.7);
**−TΔS_pocket per-run [−3.42, −1.98, −1.48]** — negative and systematic.
Honest finding: at Cα-only resolution the pocket-entropy sign is
model-dependent (ligand = 63-DOF bath injecting pocket-Cα mobility; true
restriction lives in sidechain rotors). Trustworthy at this tier: ΔH +
component split, ligand torsion ΔS, bounded |ΔS_pocket|. Heavy mode
(sidechains present) should restore the physical sign — open bench item.
Unit anchors: Schlitter 0.20%, torsion kB·ln12 exact, BindLog 1.95×.

## 4. Physics-level selector (S7)

One opt-in `<select id="physicsLevel">` in a `details.subpanel` inside the
Dynamics & Force Field panel (no new top-level panel; Digit1-7 contract
intact; 98 ui ids ⊆ index.html). Default **L0** = current baseline,
bit-identical (charges off, isotropic HB, weak off, BindLog off).

| Level | CG flags | Heavy flag | BindLog | Cost |
|---|---|---|---|---|
| L0 fast (default) | charges off, `hbMode: "off"` | `weak: "off"` | off unless checkbox | 1.0× |
| L1 balanced | charges on, `hbMode: "directional"` | off | off unless checkbox | ≈1.09× (CG) |
| L2 full-rigor | charges on, directional | `weak: "on"` | ON (effective OR with checkbox) | ≈1.0× over tier base |

Wiring (`src/main.js`): `PHYSICS_LEVELS` table + `physicsLevelSpec()`
(guarded, persists to `settingsState.physicsLevel` in-memory) +
`bindLogWanted()` (checkbox OR L2). Both build paths (`buildSystem`,
`onParamChange` hot-rebuild) carry `{charges, hbMode, weak}`; trackTerms,
funnel-hill hook, tick frame/energy/contact capture all route via
`bindLogWanted()`; guarded `change` handler hot-rebuilds. Unchecking
`bindlogOn` at L2 keeps capture on (tier owns the flag until deselected).
Headless selector contract check: 24/24 (`/tmp/opencode/s7_selector_check.mjs`:
tier table from source, real-FF flag flips, finite energies, DOM nesting).

## 5. Physics trust boundary (unchanged by Loop 2)

Qualitative ΔG/Kd ranking only — NOT FEP, NOT absolute Kd: implicit solvent,
Cα-tier entropy sign model-dependent (§3), proxy pose scorer, ±50% heap,
session-dependent timings. Tiers move along the frontier, never past
`docs/APPLICABILITY.md`.

## 6. Heavy-mode future (sidechain restores pocket sign)

Re-run the S5 apo/holo paired protocol in heavy mode with L2 full-rigor and
check −TΔS_pocket sign flip positive (2–12 kcal/mol lit band); needs
≥10 frames/DOF over sidechain-inclusive pocket + funnel-bulk torsions for
ΔS_lig on a flexible ligand (benzene ΔS_lig ≈ 0 by construction).

## 7. Verification (S7 close-out)

- Gate before: OPEN (41 files clean, 97 ids, 32/32). After: OPEN (98 ids).
- `tests/test_all.js` 32/32 · `test_thermo` 7/7 · `test_bindviz` 22/22 ·
  `test_bindlog_integration` 14/14 · selector contract 24/24.
- `node --check` on all changed files; page + changed modules served 200.
- No new deps, no new top-level panels, no sidebar/toolbar touch.

## 8. Stage-1 verdict (2026-09-12): heavy pocket-ΔS sign — NOT-RESTORED

S7 hypothesized heavy all-atom mode restores the physical pocket-entropy sign (§3,
§6 open bench item). Tested by `scripts/test_thermo_heavy.mjs` (13/13; sign never
asserted) + a 6000-frame deep-sampling pilot. Full record in
`docs/BINDING_PHYSICS_R4.md` §5 (Stage-1 verdict).

- Committed protocol (750 frames/leg, 318 DOF, f/DOF 2.4), 4 runs × 3 replicas:
  run means +4.90 / +10.26 / −1.37 / +16.47; pooled span −13.7…+59.9 →
  variance-dominated, sign undetermined (early positive means were sampling noise).
- Pilot block curve (750→6000 frames): full −36.20/−0.95/+9.64/−4.78; at 6000 frames
  (f/DOF 18.9): full −4.78, backbone −5.82, sidechain +1.04, Cα −1.26 (≈ CG scale).
- Verdict: NOT-RESTORED. Bath effect persists in heavy backbone (no ENM involved);
  sidechain restriction is real but small (+1.04, overwhelmed at 8 Å / ps scale).
- What would restore it honestly: 10+ ps × ≥3 replicas; pocket-χ torsion-Shannon term;
  sidechain-only pocket entropy. Not sidechain restraints.
- Bonus trustworthy observable: heavy ΔH = −19.7 ± 0.3 (LJ −8.2 / desolv −11.4,
  Coul/HB ≈ 0 for benzene), tight over all 12 replicas.
- Code: additive `masses` option + `S_pocket` absolutes in
  `src/analysis/thermodynamics.js` (CG path bit-identical, `test_thermo` still 7/7);
  no UI change (step-4 was YES-conditional).
- Verification at close-out: gate OPEN; `test_all` 32/32; `test_thermo` 7/7;
  `test_thermo_heavy` 13/13; `node --check` clean on all touched files.

## 9. Stage-2 note (2026-09-12): deterministic Langevin seeding — additive, opt-in

Pain: `src/integrator.js` `_fillGaussian` used unseeded `Math.random()` → every
CG/heavy replica drew a fresh noise stream, so per-run −TΔS signs flipped run
to run and `scripts/test_thermo.mjs` could only assert a 3-replica mean bound.

What changed (additive only, default unseeded path bit-identical, zero deps):
- `src/integrator.js`: `LangevinIntegrator` accepts opt-in seeding —
  `new LangevinIntegrator(ref, ff, mass, { seed })` / `{ rng }` plus
  `setSeed(n)` / `setRng(rng)` / `getSeed()` / `_randUniform()`; `_fillGaussian`
  draws from the injected uniform source when present, else `Math.random()`
  exactly as before. Backed by `src/seeded-rng.js` mulberry32; JSDoc; units
  unchanged. The `(ref, ff, { seed })` 3-arg overload is also accepted.
- `scripts/test_thermo.mjs`: fixed seeds per replica `SEEDS = [101, 202, 303]`
  (holo `SEEDS[rep]`, apo `SEEDS[rep]+1000`; `THERMO_SEED_BASE` env override);
  replica-mean assertion approach unchanged (no new single-run sign assert).
- `scripts/test_thermo_heavy.mjs`: same pattern with `SEEDS = [1001, 2002, 3003]`
  (`THERMO_HEAVY_SEED_BASE` override); `runLeg(ff, collectE, seed)`; sign still
  reported, never asserted.
- `tests/test_seeded_integrator.js` (new, 16 asserts): same seed → bit-identical
  positions/energy; different seeds → diverge; frozen golden (4W52 CG, gamma 1.0,
  T 300, zeta 5.0, seed 12345, 200 steps → energy 136.69571481112268) matches at
  1e-9; API forms (`setSeed`/`setRng`/overloads/`getSeed() null` default) work.

Deterministic numbers at close-out (seeded, replayed bit-identically twice):
- `test_thermo` 7/7: ΔH −6.82 ± 0.16; 3-rep −TΔS_pocket [4.58, −3.71, 10.42] →
  mean 3.77 (bounded < 20; sign still model-dependent at Cα res).
- `test_thermo_heavy` 13/13: ΔH mean −19.73; full −TΔS [22.63, −24.27, 13.96] →
  mean 4.11; backbone 8.96 / sidechain −4.85 / Cα 1.25 (UNDER-SAMPLED f/DOF 2.4;
  Stage-1 NOT-RESTORED verdict stands — seeding changes reproducibility only).
- `test_seeded_integrator` 16/16 × 3 consecutive runs; `test_all` 32/32; gate
  OPEN; `node --check` clean. No UI changes (headless-only; no new UI ids).

## 10. Stage-3 note (2026-09-12): tiered integration — Loop-2 guarded by CI

Pain: `tests/test_all.js` stayed 32/32 — none of the Loop-2 suites (charges,
virtual-sites, weakint, thermo, thermo-heavy, bindviz, bindlog, seeded
integrator) ran in the gate, so CI guarded none of Loop-2.

What changed (additive only, zero deps, no UI changes):
- `tests/test_all.js`: Tier-0 32 asserts untouched; new `runTier` helper spawns
  each standalone suite via `child_process.execFileSync` with a timeout and
  folds its `N PASSED, M FAILED` counts into the grand total. Child stdout is
  captured — only a lowercase one-liner per suite is printed, so the single
  uppercase results line stays the grand total. `runSuiteFile` surfaces child
  FAIL lines on error.
- `scripts/wikiskill_gate.js`: checks the LAST results line (grand total) with
  baseline 179 (was: first line, baseline 32); gate timeout 120s → 180s.

| Tier | How to run | Suites (asserts) | Runtime measured |
|---|---|---|---|
| FAST (default, gate/CI) | `node tests/test_all.js` | Tier-0 (32) + charges (20) + virtual-sites (8) + weakint (49) + seeded-integrator (16) + bindviz (22) + bindlog (18) + bindlog-integration (14) = **179** | **~16s** (weakint ~15s dominates) |
| SLOW (opt-in) | `node tests/test_all.js --slow` or `SLOW=1` | FAST 179 + test_thermo (7, seeded CG) + test_thermo_heavy (13, seeded heavy) = **199** | **~218s** (thermo 3.6s + heavy 198.7s) |

Verification at close-out: gate OPEN (before: OPEN at 32/32; after: OPEN at
179/179); `node tests/test_all.js` → 179/179 in 16.0s; `--slow` → 199/199 in
218.4s; `node --check` clean on both touched files. Default fast path stays
well under the ~60s budget, so gate/CI stays fast.

## 11. Stage-4 note (2026-09-12): altloc/dedup cleaner + rotatable-bond detection

Two known input problems closed (additive only, zero deps, default-ON,
bit-identical on clean inputs; no sidebar/toolbar/top-level-panel changes,
no new UI ids):

(a) Altloc/dedup cleaner — root-causes the S7 coincident-cation/centroid NaN
(dirty input: overlapping duplicate-altloc atoms in the raw crystal PDB) at
the input instead of guarding the symptom. Single shared rule in NEW
`src/pdb_altloc.js` (`parseAltLoc` col 17 / `parseOccupancy` cols 55–60,
blank occupancy → 1.0; `shouldReplaceAltloc`), wired into all three PDB
parsing paths at their existing `seen`-duplicate choke points:
`parseCa` + `parseLigands` (`src/pdb.js`) and `parseHeavy` (`src/heavy.js`).
Rule: dedup key (chain, resSeq, iCode, atom name) keeps the
highest-occupancy copy; tie (≤1e-9) prefers altLoc 'A' order-independently,
else first; zero-occupancy copies dropped when a non-zero copy exists;
exact duplicates (same altLoc/coords, repeated serial) collapse to one atom;
every skip/replace logged to `warnings` + `console.warn` (count preserved in
the existing `[parse*] N warning(s)` summary). Default-ON is safe because
4W52's 160 altloc lines are all 0.50/0.50 ties with 'A' first — the new rule
keeps the same record the old first-wins code kept:
4W52 164 Cα / 1308 heavy / 2 hetero + 6-atom BNZ / 15-atom EPE ligands and
1CRN 46 Cα / 327 heavy are byte-identical before/after (asserted in the new
test). `parseLigands` output shape unchanged (internal `_atomKeys` map only,
deleted before return) so downstream CG/heavy paths are untouched.
`tests/test_altloc_cleaner.js` (new, 19 asserts): synthetic occupancy-winner,
A-on-tie (both orders), zero-occupancy drop, heavy sidechain winner,
exact-duplicate HETATM collapse (heavy + ligand, incl. no-coincident-pair
check d = 1.000 Å), warning counts, and the 4W52/1CRN regression counts.

(b) Rotatable-bond detection + flexible-ligand ΔS — NEW
`src/analysis/rotbonds.js` (`findRotatableBonds` / `autoTorsions`, zero-dep,
JSDoc): single non-aromatic bond (order≈1; assumed single when untyped —
rings/terminals cover the risk), non-ring (per-bond BFS cycle test, any ring
size), non-amide (gaffType 'nh' → explicit C=O orders → <1.35 Å geometry
fallback), both sides heavy-degree > 1 (excludes methyl/hydroxyl trivials),
plus symmetric-top fan exclusion (≥3 leaf neighbors ⇒ tert-butyl-like /
S(=O)₃ head excluded — documented). One deterministic torsion quadruplet per
rotatable bond (most-substituted arms, smallest index on ties). Measured:
benzene (MOL2 and 4W52 BNZ) → 0; 4W52 HEPES/EPE (15 bonds) → 4
([[0,9],[3,6],[6,7],[9,10]]; ring + S(=O)₃ fan excluded). Wired into
`computeThermodynamics` as an autoTorsions fallback: explicit `p.torsions`
always wins; else a supplied ligand graph (`p.ligand` or
`p.ligandAtoms/Bonds` + offset) auto-detects; no graph ⇒ legacy rigid path
(ΔS_lig = 0), bit-identical. `meta.rotatableBonds` added and the
`formatThermoTable` ligand line now reports the count
(e.g. `(1 auto rotatable bonds (1 torsions) [1 rotatable])`).
`tests/test_rotbonds.js` (new, 17 asserts): rigid 0s, EPE 4, chain/amide/
t-butyl rules, locked S = 0.00000 vs flexible S = 0.00494 (600-frame
12-bin ensemble, kB·ln12 scale), auto offset, ΔS_lig > 0 via fallback
(0.00494, `meta.rotatableBonds` 1), rigid/legacy exactly 0, table-line count.

Verification at close-out: gate OPEN (before: OPEN at 179; after: OPEN at
215); `node tests/test_all.js` → 215/215 (Tier-0 32 + FAST 183, ~13s);
`scripts/test_thermo.mjs` still 7/7 deterministic (ΔH −6.82 ± 0.16, 3-rep
−TΔS_pocket [4.58, −3.71, 10.42] — identical to the Stage-2 record; ligand
line now reads `rigid … ⇒ 0 [0 rotatable]`); `node --check` clean on all
9 touched/new files. Gate baseline bumped 179 → 215
(`scripts/wikiskill_gate.js`, `tests/test_all.js` comment).

## 12. Stage-5 note (2026-09-12): async thermo apo leg + physics-level persistence

Two known UI debts closed (additive only, zero deps, defaults unchanged:
default L0 + default thermo numbers; no sidebar/toolbar/top-level-panel
changes; no new DOM ids — all 98 ui ids ⊆ index.html, `src/ui.js` untouched):

(a) Async/chunked thermo (`src/analysis-panel.js:129-236`, `src/main.js:451-453`):
pain was the ΔH/ΔS click handler running the full ~2000-step apo relaxation
synchronously → blocked the UI. The apo leg now runs in `setTimeout(0)`
slices of `THERMO_CHUNK_STEPS = 150` steps (exported,
`src/analysis-panel.js:136`; 600-step min → 4 slices, 2000-step max → 14).
Each slice posts progress to the existing `#thermoCaption`
(`relaxing apo… N/M steps`), `thermoBtn` is disabled during the run and
re-enabled in the `finally` (or on setup error), errors render `⚠ …` to
`#analysisOut` without stranding the button, and a module generation counter
(`_thermoGen`, `invalidateThermo()`) cancels a stale run on rebuild
(`buildSystem` calls it first) or on a fresh click. Holo leg still uses the
recorded frames; the final render is the identical `formatThermoTable(res)` +
BindLog-off note; integrator params, `nApo` formula, pocket rule unchanged;
null-safe + guarded throughout (headless import clean).

(b) Selector persistence (`src/settings-panel.js:49-90`,
`src/main.js:35,45-47,645-650`): the L0/L1/L2 selector was in-memory only
(`settingsState.physicsLevel`) → reset on reload. No other setting persisted
(all in-memory), so this follows the `specs/PLAN.md` "persist to
`localStorage`" convention: NEW key `sim.physicsLevel`
(`PHYSICS_LEVEL_KEY`), guarded `readStoredPhysicsLevel()` /
`persistPhysicsLevel(lvl)` (validates L0/L1/L2, never throws — private-mode /
headless safe) / `restorePhysicsLevelSelect()` (syncs state + select element).
Restored at module load and in `initSettingsModal()` + once more from
`main.js` after panel init; the existing `change` handler persists on every
switch. Default L0 when absent/invalid; bogus values ignored.

Verification at close-out: gate OPEN before (215) and after (215);
`node tests/test_all.js` → 215/215 fast (~13s); `node --check` clean on all
3 touched files; headless imports clean (`analysis-panel` exports
`THERMO_CHUNK_STEPS=150` + `invalidateThermo=function`, `settings-panel`
exports `KEY=sim.physicsLevel`, fresh default L0); page + all touched modules
serve 200; DOM contract still 98 ids. Responsiveness proof (real 4W52 CG apo
`ForceField` + `LangevinIntegrator`, 600 steps, chunk 150): 4 slices,
heartbeat timer fired 3× mid-run (event loop free vs 0 for a sync loop),
captions `relaxing apo… 150/600` → `600/600`, button disabled→re-enabled.
Persistence proof (localStorage stub): fresh → L0; persist L2 → stored
`sim.physicsLevel=L2`, reload read → L2, select restores L2; bogus/invalid
(`L9`) → ignored/null → fall back L0.

## 13. Stage-6 note (2026-09-12): 4W52 calibration anchor — measured, not fitted

Pain: the physics trust boundary (§5) said "qualitative ΔG/Kd, not FEP" with
no number, ΔH was desolv-dominated with no share quoted, the SASA ±50% band
had no kcal value, and alanine-scan was never checked against literature.
`scripts/calibration_4w52.mjs` (new, 10 asserts, SLOW tier, seeded 101/1101,
deterministic bit-identical twice, ~2s CG-thermo class) measures all four.
Zero deps, no parameter retuning, no UI changes (headless-only, no new ids).

(a) Identity correction (measured from `4w52.pdb`): 4W52 is **T4 lysozyme L99A
+ benzene** (TITLE/COMPND; ligands BNZ 6 atoms + EPE buffer 15 atoms; 164
residues, no Zn, single HIS31) — NOT carbonic anhydrase II. The CA-II hotspot
premise (His94/96/119, Thr199, Leu198) is absent by construction (here
VAL94/ARG96/ARG119; resSeq 198/199 do not exist). Judged against the true T4L
cavity liners (Merski et al. 2015 PNAS 112:5039, the 4W52 primary citation).

(b) Anchor table (kcal/mol; A/B measured live here, heavy cited from §8–9):

| Row | ΔH ± SE (LJ/Coul/HB/desolv) | −TΔS_pocket | ΔG_est | Experimental ΔG |
|---|---|---|---|---|
| A record path (BNZ+EPE, seeded rep0) | −6.82 ± 0.16 (−1.88/+1.04/−0.13/−5.85) | +4.58 | −2.23 | ITC −5.20 ± 0.20 / NMR −4.20 ± 0.10 |
| B BNZ-only control (true cavity) | −3.64 ± 0.08 (−0.98/0.00/0.00/−2.66) | +1.47 | −2.17 | (same) |
| Heavy record (Stage-1/2, cited) | −19.73 mean (−8.2/≈0/≈0/−11.4) | mean +4.11 | — | (same) |

Experimental sources (stated, not re-measured; PDBbind/MOAD carry no direct
4W52 Kd — T4L is not a drug target): Mondal et al., PLoS Comput Biol
14:e1006180 (2018), Table 1, quoting T4L-L99A/benzene ITC (−5.2 ± 0.2) and NMR
(−4.2 ± 0.1, kon≈1e6 M−1s−1/koff≈950 s−1). Cross-check at 300 K
(RT = 0.596): Kd 0.15 mM → −5.25; Kd 0.8 mM → −4.25. Consistent.
3-replica context (`test_thermo.mjs`, measured 2026-09-12, seeds 101/202/303):
ΔH [−6.82, −7.50, −6.00] (mean −6.77, SD 0.75); −TΔS [+4.58, −3.71, +10.42]
(mean +3.76, SD 7.1); ΔG_est [−2.23, −11.20, +4.42] (mean −3.00).
Single-rep error vs experiment: rep0 ΔG_est −2.23 underbinds by +2.97 (ITC) /
+1.97 (NMR) with bootstrap ±0.16 — but the replica SD (±7.1, entropy-driven)
dominates: absolute ΔG is sampling-noise-limited, ranking-only use confirmed.
Two honest calibration findings: (i) the record-path pocket (11 residues
ALA74/MET102…GLY110, all-mol COM pulled by surface EPE) misses the benzene
cavity (B pocket: 14 residues ILE78/LEU84/…/LEU118 around the BNZ COM);
(ii) EPE inflates |ΔH| by −3.18 (A−B) — reported, not corrected.

(c) Desolv share (measured): A 85.8% (−5.85/−6.82), B 73.1% (−2.66/−3.64),
heavy 57.9% (−11.4/−19.7, §8). ΔH is desolv-dominated at both tiers; the CG
Coul term (+1.04, repulsive, EPE-driven — B Coul is 0.00) is the other
buffer-pollution signature.

(d) SASA sensitivity (measured): as-run ΔSASA = 0 Å² (no contact-count series)
→ −TΔS_solv 0.00, ±50% swing **±0.00**. Illustrative (ASSUMED ΔSASA 180 Å²
benzene-burial upper bound): −2.16, band [−1.08, −3.24], swing ±1.08 —
i.e. even full burial moves ΔG by ≈2 kcal/mol, below the replica SD (7.1).

(e) Ala-scan spot-check (CG BNZ-only, rCut 8 Å/maxN 12, relax 80,
deterministic): top-3 TYR88:A +0.014 / MET102:A +0.012 / ILE100:A +0.009;
max|ΔΔG| 0.015; ALA99/ALA98 defined 0. Verdict: **PARTIAL nominal (1/3,
MET102 is a genuine cavity liner) but effectively NO hotspot resolution** —
all |ΔΔG| < 0.02, ranking is noise. Why (CG limits): the ENM perturbation
(0.7× + 0.4 Å r0 shift) relaxes back and cancels in the holo−apo cycle, and
sidechain contacts are invisible at Cα — magnitudes, not just order, fail.

(f) Hardened trust boundary (numbers): CG ΔH replays the Loop-2 record to
±0.3 (seeded replay lock) and lands within 1.4–2.5× of the experimental
ΔG scale, but ΔG_est carries ±0.16 bootstrap ⊕ ±7 replica-SD error and
underbinds by 2–3 kcal/mol single-rep — **ranking-only use; NOT FEP, NOT
absolute Kd** (prior §5 stands, now quantified). Absolute-ΔG claims additionally
require: BNZ-only (buffer-free) legs, ≥10 frames/DOF Schlitter, and a nonzero
measured ΔSASA series.

Verification at close-out: gate OPEN before (215) and after (215);
`node tests/test_all.js` → 215/215 fast (~13s); `--slow` → 245/245
(Tier-0 32 + FAST 183 + SLOW 30: thermo 7 in 2.7s + heavy 13 in 169.8s +
calibration 10 in 1.9s); calibration log diff-clean twice;
`node --check` clean on all 4 touched/new files. SLOW_SUITES + gate comment
updated; FAST tier untouched (still 215).

## 14. Stage-7 note (2026-09-12): perf honesty (allocation tracking) + GPU-port decision + docs closeout

(a) Allocation tracking (`scripts/pareto_bench.mjs:13-24,102-170,300-360`,
additive, zero deps, timing logic untouched): each tier now brackets its timed
loop with `heapUsed` before/after reads (`heapDelta`); with `node --expose-gc`
`gc()` is forced at each bracket (near-true allocation), without it the
GC-noise caveat is printed WITH the measured numbers. FULL reruns (25 s each,
so L2/L3 re-measured, none estimated) on 4W52:

| Tier | run1 ms | run2 ms | gc-run ms | gc heap Δ (MB) | no-gc heap Δ run1/run2 (MB) |
|---|---|---|---|---|---|
| L0 CG | 0.225 | 0.218 | 0.237 | +0.29 | −0.25 / −0.64 |
| L1 CG+ | 0.198 | 0.208 | 0.235 | +0.09 | +3.22 / −1.81 |
| L1 CG+ +BindLog | 0.201 | 0.215 | 0.216 | +0.06 | +1.27 / −0.02 |
| L2 heavy | 78.82 | 79.65 | 84.80 | +0.14 | −7.27 / +1.56 |
| L3 heavy+R3 | 76.13 | 77.51 | 83.06 | +0.18 | −1.41 / +15.57 |
| L2 full-rigor | 77.07 | 80.14 | 85.49 | −0.25 | −4.93 / −3.98 |

ms stable ±5% across runs; pose recovery unchanged (10/10 CG, 9/10 heavy,
meanRank 1.0/1.3). Without gc, deltas swing ±15 MB (floating garbage —
this WAS the ±50% snapshot debt: snapshots 7–27 MB). With gc, every tier
allocates ≤0.3 MB per timed block (≈0–2 KB/step; CG force buffers are
reused `Float64Array`s). Memory is NOT a frontier axis at this size — now
quantified, not asserted. `docs/pareto_frontier.csv` carries the new
`heap_delta_MB` + `run` columns: 7 Loop-2 S7 history rows preserved verbatim
(empty delta, run `Loop2-S7-FULL history pre-alloc-tracking`) + 7 fresh dated
rows (`2026-09-12 FULL`, gc-bracketed); the writer merges by run tag, so
reruns replace only their own tag and never delete history.

(b) GPU-port decision — ACCEPT CPU, no port, no GPU code changes. Grep over
`src/compute/wgsl/`: zero Loop-2 binding terms in WGSL (no bindU, weakint,
desolv, virtual-sites; the single "directional" hit is the header comment at
`src/compute/wgsl/nonbonded_forces.wgsl:16-17` stating directional H-bonds,
SASA, membrane and bonded terms stay CPU-side by design). WGSL covers LJ +
screened Coulomb + Still-GB only. Port cost: virtual-site angular gates and
weakint ring/cation detection are branchy gather-heavy kernels — the same
class already deliberately kept on CPU (Born-radii HCT precedent,
`nonbonded_forces.wgsl:8-12`); a port buys a second physics implementation to
hold in parity. Benefit ≈ 0: CG tiers cost 0.22–0.24 ms/step on CPU (nothing
to win) and Loop-2 adds ~0% on heavy (≈1.0× L3/L2-full vs L2), so binding
terms are not the bottleneck the GPU would relieve. Verdict: all Loop-2
binding terms stay CPU-only; GPU stays LJ/Coulomb/GB.

(c) Docs closeout: README gains `## Binding-physics Loop-2 + Stages 1–7`
(L0/L1/L2 selector, BindLog/Viz, thermo, seeding, calibration ranking-only
NOT-FEP boundary, test tiers fast 215 / slow 245); verdict appended to
`.wikiskill/wiki/evolution/skill-impact.md` (gate OPEN, 215 fast / 245 slow,
Stages 1–7); R5 addendum in `docs/BINDING_PHYSICS_R5.md` §7 points here.

Verification at close-out: gate OPEN before (215) and after (215);
`node tests/test_all.js` → 215/215 fast; `node --check` clean on all touched
files; no UI changes (no new ids, DOM contract untouched, no curl target).
