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

## 15. Stage-1 deep-sampling follow-up (2026-09-12): 10 ps × 3 + χ-term + sidechain-only

Prescription from §8 ((a) 10+ ps × ≥3 replicas, (b) pocket-χ torsion-Shannon
term, (c) sidechain-only pocket entropy). New `scripts/test_pocket_entropy.mjs`
(`test_thermo_heavy.mjs` untouched); additive thermo helpers in
`src/analysis/thermodynamics.js` (`BACKBONE_ATOM_NAMES` /
`isBackboneAtomName` / `splitPocketByBackbone`, `CHI_DEFS` /
`pocketChiTorsions`, `chiEntropyDelta` with jackknife SE, `formatChiLine` +
opt-in `p.chiTorsions → r.chi` passthrough — CG path bit-identical, no-χ
callers get `r.chi === null` and byte-identical tables).

Protocol ps: heavy dt = 0.001 ps = 1 fs (`src/integrator.js:_pickDt`,
measured on this system), so LONG (`--long`/`--slow`, opt-in, never default)
= 500 equil + 10000 production steps, stride 2 → 5000 frames/leg =
**10.0 ps production + 0.5 ps equil** per leg × 3 seeded replicas
(holo SEEDS 1001/2002/3003, apo +1000; `THERMO_CHI_SEED_BASE` override).
PILOT (default, fast) = 50 + 600 → 300 frames/leg (0.60 ps), 1 replica.
Pocket 106 heavy atoms (56 backbone N/CA/C/O + 50 sidechain CB-outward, 14 Cα;
318 DOF) → LONG frames/DOF **15.7** (≥10 target).

χ coverage (documented, ALA-free): 25 pocket residues → 16 with sidechain in
pocket → **14 χ residues** (ALA98/99 skipped — CB but no rotor; 0 unresolved)
→ **26 χ torsions** (14 χ1 + 11 χ2 + MET χ3: ILE78/100 ×2, LEU84/91/118/121/
133 ×2, VAL87/103/111 ×1, TYR88 ×2, MET102 ×3, PHE114/153 ×2) → LONG
frames/rotor **192** (target ~100).

Per-replica LONG table (5000 f/leg, −TΔS kcal/mol):

| rep | full | bb | sc (sidechain-only) | Cα | χ −TΔS ± jackSE | ΔH |
|---|---|---|---|---|---|---|
| 0 | −8.47 | −4.94 | −3.53 | −2.68 | +0.27 ± 0.30 (S_holo 0.0456 vs S_apo 0.0465) | −20.32 |
| 1 | +9.53 | +5.88 | +3.65 | +1.57 | −0.18 ± 0.26 | −20.40 |
| 2 | +19.23 | +13.79 | +5.44 | +4.76 | +0.62 ± 0.44 | −19.57 |
| mean | **+6.76 (SD 14.06)** | +4.91 | **+1.85 (repSD 4.75)** | +1.22 | **+0.24 (±0.33)** | **−20.10** (LJ −7.60 / Coul 0.00 / HB 0.00 / desolv −12.50) |

Block-convergence curve (replica-mean −TΔS, first-N frames):

| N | f/DOF | fr/rot | full | bb | sc | ca | χ |
|---|---|---|---|---|---|---|---|
| 750 | 2.4 | 29 | −2.12 | −0.54 | −1.58 | −1.06 | +0.28 |
| 1500 | 4.7 | 58 | −11.41 | −3.99 | −7.42 | −1.85 | +0.06 |
| 3000 | 9.4 | 115 | +5.14 | +7.10 | −1.96 | +0.85 | +0.05 |
| 5000 | 15.7 | 192 | +6.76 | +4.91 | +1.85 | +1.22 | +0.24 |

Verdict: **NO — no converged positive-restriction signal separable from
backbone rattle at 10 ps × 3** (no forced sign; either outcome publishable).
Sampling floor dominates: replica SD 14.06 (full) / 4.75 (sc) ≫ means, and
block signs flip through 5000 frames (sc −1.58→−7.42→−1.96→+1.85). χ is
stable-positive across blocks but +0.24 ± 0.33 — consistent with zero (upper
bound ~0.5 kcal/mol); backbone +4.91 exceeds sidechain +1.85 (sc−bb −3.06),
so nothing is separable. Consistent with the 6000f pilot (full −4.78): both
are variance-dominated with undetermined sign — now quantified at adequate
conditioning. What remains precisely: (i) sampling floor — need ~10× longer
or orthogonal enhanced sampling to shrink replica SD below ~2 kcal/mol;
(ii) possible force-field ceiling — pocket-χ restriction genuinely ≤0.5
kcal/mol at the 8 Å / ps scale while the backbone ligand-bath drives the net.
ΔH stays the trustworthy observable (−20.10, tight over replicas, reproduces
the §8 heavy −19.7 ± 0.3).

Verification at close-out: gate OPEN before and after; `node
tests/test_all.js` → 215/215 fast (~19 s, unchanged — new script in no tier,
LONG opt-in only); PILOT 16/16 twice bit-identical (modulo the wall line);
LONG 16/16 once (wall 1410 s); `node --check` clean on both touched files.
No UI changes, no parameter retuning, zero deps.

## 16. Stage-2 note (2026-09-12): thermo ligand picker — BNZ cavity default

Pain (Stage-6 finding, §13): the thermo record path used ALL HETATM groups
(BNZ+EPE) for the pocket COM + binding-energy attribution. Surface EPE pulls
the all-mol COM ~10.8 Å off the benzene cavity (record pocket 11 residues
ALA74/MET102…GLY110 surface loop vs BNZ cavity 14 residues
ILE78/LEU84/…/LEU118) and inflates |ΔH| by −3.18 (row A −6.82 vs row B
BNZ-only −3.64). History is NOT rewritten — the record stands; new runs
default to the buffer-free cavity.

What changed (additive only, zero deps, pocket 8 Å rule unchanged, no new
top-level panels, Digit1-7 contract intact):
- NEW `src/analysis/thermo_ligand.js` (helpers, JSDoc): `resolveThermoLigand`
  (rule, stated: "auto" = BNZ-first fallback — all BNZ molecules when present,
  else first group, else none; "all" = every group = record path; explicit
  index or resName = that subset), `ligandAtomOffsets` /
  `selectedLigandAtomIndices` (FF ref-suffix concatenation order),
  `selectedLigandCom` / `pocketFromCom` (8 Å) / `sliceSelectedPositions`
  (full-frame → selected-only slicing for energy recomputation),
  `refreshThermoLigOptions` (rebuilds options, preserves choice, headless-safe).
- `index.html` thermo subpanel: NEW `<select id="thermoLig">` (auto/cavity
  BNZ-only default + per-group options + all) following the existing
  subpanel/select pattern (+ one-line summary); registered in `src/ui.js`
  (`thermoLig`, in-memory default auto, persist NOT required).
- `src/analysis-panel.js`: handler reads the picker (headless-safe, default
  auto), pocket COM from the SELECTED ligand only, ΔH attribution from the
  SELECTED ligand only ("all"/full-coverage reuses the BindLog channel
  bit-identically; a proper subset is recomputed per recorded frame with a
  selected-only FF copying the live binding flags; apo leg built with the
  subset), result caption names the ligand + pocket count. `all` reproduces
  the Loop-2 record bit-for-bit.
- `src/main.js`: `buildSystem` repopulates the picker via exported
  `refreshThermoLigPicker()` (guarded, additive).

Numbers (fast headless check, calibration protocol 2000 steps/stride 2,
seeds 101/1101, charges + directional-HB — `/tmp/opencode/stage2_verify.mjs`):
- Picker: auto → `auto-bnz`, BNZ #1 (6 atoms) PASS; "0" → BNZ, "1" → EPE #2
  (15 atoms); all → both groups (21 atoms).
- Pocket: auto n=14 (ILE78 LEU84 LYS85 VAL87 TYR88 ARG96 ALA98 ALA99 ILE100
  MET102 VAL103 VAL111 ALA112 LEU118) vs all n=11 (ALA74 MET102 VAL103 PHE104
  GLN105 MET106 GLY107 GLU108 THR109 GLY110 VAL111) — matches §13.
- ΔH ALL (record) −6.82 ± 0.16 (LJ −1.88/Coul +1.04/HB −0.13/desolv −5.85)
  replays row A ±0.30 PASS (Δ=0.00); AUTO/BNZ-only −3.64 ± 0.08 (LJ
  −0.98/Coul 0.00/HB 0.00/desolv −2.66) replays row B PASS (Δ=0.00);
  EPE-only −3.76 ± 0.17 (Coul +0.77 — the buffer-pollution signature);
  EPE inflation ALL−BNZ = −3.18 confirmed.
- Slice/recompute sanity: 200 ALL-trajectory frames sliced 21→6 atoms,
  selected-only mean U_bind finite PASS.

Verification at close-out: gate OPEN before (215, 98 ids) and after (215,
99 ids — only +thermoLig); `node tests/test_all.js` → 215/215 fast (one
flaky 201/1 run mid-session, then 215/215 twice + gate 215 — temperature-trial
flake, not a code regression); `node --check` clean on all touched files;
page + touched modules serve 200 (index.html, analysis-panel.js,
thermo_ligand.js); no new top-level panels (8 unchanged).

## 17. Stage-3 note (2026-09-13): real solvent burial — LCPO wired into holo−apo

Pain (open debt from §13): the solvent term was dead — as-run ΔSASA = 0 (no
contact-count series ever passed) → −TΔS_solv 0.00 with a ±50% swing of
±0.00. The LCPO module (`src/physics/solvation/lcpo_sasa.js`: `lcpoSasa`,
`capAreaDeriv`, `vdwRadiusFor`, `PROBE_RADIUS`) existed but was never wired
into the holo−apo thermo legs.

What changed (additive only, zero deps, JSDoc; existing ΔSASA = 0 callers
unchanged when SASA absent; no new top-level panels; no new DOM ids):
- NEW `src/analysis/thermo_sasa.js` (exports `THERMO_SASA_STRIDE = 10`
  `:55`, `THERMO_SASA_CHUNK = 25` `:58`, `cgBeadExtendedRadius` `:65`,
  `selectedLigandElements` `:75`, `ligandExtendedRadii`/`isolatedAreas`,
  `frameLigandBurial` `:125`, `frameProteinSasa` `:171`,
  `sasaBurial` `:219` + chunked `sasaBurialChunked` `:277`).
  Method: ΔSASA = Δ_lig + Δ_prot with ΔSASA = ⟨apo⟩ − ⟨holo⟩ (positive =
  burial). Δ_lig (dominant): per holo frame, free-ligand LCPO areas
  (ligand-only subsystem, in-regime for the small molecule) × cross-burial
  survival from protein neighbours
  (`buried_a = A_free,a · (1 − Π_p max(0, 1 − C_pa/S_a))`, caps from
  `capAreaDeriv`); Δ_prot: ⟨S_prot⟩ stripped-protein full-LCPO apo minus
  holo (same regime both legs, self-burial cancels). SEs are SEMs over
  evaluated frame means; dsasa SE = √(SE_lig² + SE_prot²), legs treated as
  independent (documented). Why NOT full-complex LCPO totals (measured 4W52
  CG BNZ seeds 101/1101): bead-scale P3/P4 fits are out of regime
  (bound-ligand LCPO 419.8 > free 340.8 by 79 Å² — anti-burial inversion;
  full SEP −49.0 ± 3.9 Å² unphysical; literal apo−holo totals −15.0 Å² —
  apo keeps the pinned ligand per protocol). The cross-burial route is
  monotonic by construction. Cost (documented `:40-44`): per holo frame one
  ligand-only LCPO O(nL²) + O(nProt·nSel) cross distances; per
  stripped-protein frame one full LCPO O(nProt²). CG 4W52 (164+6) ≈
  0.3 ms/frame → 100 frames/leg sub-second; stride 10 keeps 10k-frame
  recordings tractable; the chunked variant keeps the UI responsive.
- `src/analysis/thermodynamics.js`: `SASA_GAMMA = 0.012` surfaced as a named
  const `:24` (value unchanged); opt-in `p.sasa` override `:402-412` +
  solvent path `:498-518` (`r.meta.sasa` carried; absent ⇒ legacy proxy
  path bit-identical, ΔSASA = 0 default, `r.meta.sasa` null); real-burial
  solvent table line `:562-569` (ΔSASA ± SE, method, stride, frame counts;
  legacy line byte-identical otherwise). −TΔS_solv = GAMMA × ΔSASA with SE
  (`solvSE = γ·se`; `dS_solv = −dsasa·γ/T`).
- `src/analysis-panel.js` thermo handler `:295-365`: after the chunked apo
  leg, computes SASA for recorded holo + relaxed apo via
  `sasaBurialChunked` (chunked/async like the apo leg, progress to
  `#thermoCaption` — `solvent SASA… d/t frames` — stale generations
  discarded), protein blocks from both legs + SELECTED-ligand subset
  (Stage-2 picker via `selectedLigandAtomIndices`/`selectedLigandElements`;
  whole-tail-block fallback + honest note on layout mismatch; heavy mode
  uses element-derived radii, CG the bead radius), renders the solvent line
  with real numbers via `formatThermoTable` + caption
  (`ΔSASA x.x ± y.y Å²`); SASA failure falls back to the legacy ΔSASA = 0
  path with an honest note (never strands the report).
- `scripts/calibration_4w52.mjs` re-run (`:191-225`, asserts `:254-257`):
  real burial on the BNZ-only legs (seeds 101/1101, `SASA_STRIDE = 10`
  `:70-71`, `sasaBurial` + `computeThermodynamics` sasa override, CG bead
  radius 3.4 Å, BNZ 6-C selected elements).

Numbers (BNZ-only, 1000 frames/leg → stride 10 → 100+100 evals, live):
- ΔSASA **167.1 ± 3.5 Å²** (dLig 158.9 ± 2.6 + dProt 8.2 ± 2.4;
  free-ligand ⟨S⟩ 340.8 Å²) — positive BNZ-scale burial (assert 100–250).
  Stride-5 check: 167.3 ± 2.5 (Δ 0.2 Å² — subsampling-insensitive).
- −TΔS_solv **+2.01 ± 0.04** (= γ·ΔSASA, γ = 0.012) — nonzero with SE
  (assert > 0.5).
- Revised anchor-B ΔG_est **−0.16** (legacy −2.17 + 2.01; assert exact to
  1e-9) with the ±50% band on the REAL base: solvent [+1.00, +3.01],
  ΔG [−1.17, +0.84], swing ±1.00.
- Sign-convention flag (honest, open for physics review): the code gives
  UNfavorable +2.01 (`dS_solv = −ΔSASA·γ/T`); the §13 illustrative row used
  favorable −2.16 (−γ·ΔSASA). Hydrophobic release argues favorable, so a
  sign review is owed — this stage reports the code-exact numbers and does
  NOT silently flip physics (legacy ΔSASA = 0 callers are exactly 0 either
  way). Absolute-ΔG reading additionally stays sampling-limited (replica SD
  ±7.1, §13) — ranking-only use stands.

Verification at close-out: gate OPEN before (215 PASSED, 99 ids) and after
(215 PASSED, 99 ids); `node tests/test_all.js` → 215/215 fast before and
after (~13 s); `node scripts/calibration_4w52.mjs` → 14/14 (10 legacy +
4 new) with the burial block above; calibration log diff-clean twice
(bit-identical, no RNG); `node --check` clean on all 8 touched files
(calibration, test_all, thermodynamics, thermo_sasa, thermo_ligand,
analysis-panel, ui, main). SLOW-tier total moves 245 → 249
(Tier-0 32 + FAST 183 + SLOW 34: thermo 7 + heavy 13 + calibration 14);
FAST/gate baseline untouched (still 215). No new top-level panels (8
unchanged); no new DOM ids (all 99 ⊆ index.html, `src/ui.js` untouched
this stage).

## 18. Binding-Loop2 Stage-4 note (2026-09-13): ala-scan honest display — Option B

Pain (open debt from §13e): the CG alanine-scan top-3 (TYR88 +0.014,
MET102 +0.012, ILE100 +0.009, max|ΔΔG| 0.015) sit an order of magnitude
below any plausible relaxation/convergence precision — the Cα-ENM
perturbation (0.7× + 0.4 Å r0 shift) relaxes back and cancels in the
holo−apo cycle with sidechains invisible. Verdict was prose-only
("PARTIAL nominal, effectively NO"); nothing in code stopped a reader
from quoting +0.014 as a hotspot ΔΔG.

Option tried and rejected (measured, cheap): sidechain-heavy-atom-count
scaling (Option A) — ddG × count (TYR×8, MET×4, ILE×4…) gives TYR88
0.112 / MET102 0.048 / ILE100 0.037 with the SAME top-3 in the SAME
order and unchanged 1/3 liner overlap. Arbitrary inflation with no new
discriminating power — rejected as dishonest; documented here instead
of shipped. (Buried-SASA / contact-persistence variants not pursued:
same verdict class — multipliers on a cancelled number.)

What shipped (Option B, formal ranking-only demotion; additive only,
zero deps, JSDoc; no force-field retuning — display only; defaults
unchanged except honest display):
- `src/analysis/alanine_scanning.js`: NEW `ALA_DDG_NOISE_FLOOR = 0.05`
  kcal/mol `:73` + pure helpers `isNoiseDdG` (non-finite → noise, never
  claims signal on NaN) / `annotateScanNoise` (tags `{noise, noiseFloor}`
  in place, deterministic, no RNG); `scanPocket` auto-tags rows
  (`opts.noiseFloor` override, default 0.05); `formatMutationTable`
  appends ` ~noise` per below-floor row plus a footer
  `[ala-scan: |ΔΔG| < 0.05 kcal/mol ≈ noise — ranking only, no CG
  hotspot resolution]`. Energies untouched — magnitudes identical.
- `src/analysis-panel.js:412-415`: ala-scan header gains `(ranking only
  — |ΔΔG| < 0.05 ≈ noise, no CG hotspot claim)` (hotspot-recall
  disclaimer in caption; no new DOM ids).
- NEW `tests/test_ala_noise_floor.js` (16 asserts, ~1 s, deterministic):
  floor value, boundary/custom-floor/NaN flagging, synthetic tagging +
  table flags, live 4W52 BNZ cavity (rCut 8 Å/maxN 12/relax 80):
  finite ×12, max|ΔΔG| < floor, every row noise, ordering bit-identical
  across two runs, top-3 still contains MET102, live table flagged.
- `tests/test_all.js` FAST gains the suite (183 → 199 FAST, total
  215 → 231); `scripts/wikiskill_gate.js` baseline 215 → 231.

Top-3 before/after (4W52 BNZ-only, relax 80, live, deterministic twice):
before TYR88:A +0.014 / MET102:A +0.012 / ILE100:A +0.009
(max|ΔΔG| 0.0148); after IDENTICAL numbers, every row `~noise` +
ranking-only footer. Hotspot verdict: 1/3 nominal (MET102 genuine T4L
liner; TYR88/ILE100 not liners) but formally **NO resolution** — now
encoded in code (`noise === true` ×12, table + caption disclaimers),
not just prose. Honest resolution path stays open: heavy-mode explicit
sidechain surgery (bench item, not this stage).

Verification at close-out: gate OPEN before (215 PASSED, 99 ids) and
after (231 PASSED, 99 ids); `node tests/test_all.js` → 231/231 fast
(~16 s); new suite 16/16 twice bit-identical; live top-3 replayed
twice bit-identical (TYR88 +0.014 / MET102 +0.012 / ILE100 +0.009);
`node --check` clean on all 5 touched/new files (alanine_scanning,
analysis-panel, test_ala_noise_floor, test_all, wikiskill_gate). No new
top-level panels (8 unchanged); no new DOM ids (`src/ui.js` untouched).

## 19. Stage-5 note (2026-09-13): second validation system — 4W52 EPE flexible ligand

Pain (open debt): everything validated on 4W52/T4L only (rigid BNZ);
the rotbonds flexible-ligand path (EPE 4 rotatable, chain 1) was
unit-tested (`test_rotbonds` 17/17) but never end-to-end in-app.

System choice: **4W52 EPE pocket — same protein, different ligand
(cheapest)**. EPE (HEPES buffer, 15 atoms / 15 bonds) → **4 rotatable
[[0,9],[3,6],[6,7],[9,10]]** (ring + S(=O)₃ fan excluded) vs BNZ rigid
0. Same protein isolates the flexibility variable; crystal pose needs no
placement; pocket differs honestly (EPE surface site 7 residues vs BNZ
cavity 14), so this is a genuine second system. Rejected: a MOL2 library
ligand with ≥2 rotatable bonds placed headlessly — measured on
`src/ligandLib.js`, max is ethanolamine 1 (benzene/phenol/toluene/
chlorobenzene/indole/imidazole/acetate/DMSO/caffeine all 0; methyl rotors
excluded by the degree>1 rule), so no library ligand reaches ≥2.

What shipped (additive only, zero deps, JSDoc; fast tier runtime
unchanged; no new top-level panels; no new DOM ids — `src/ui.js`
untouched):
- NEW `scripts/validate_flexlig.mjs` (10 asserts, SLOW tier, seeded
  501/1501, deterministic bit-identical twice, ~0.6 s < 60 s): EPE-only +
  BNZ-only CG systems (charges + directional-HB), 500 steps / stride 2 →
  250 frames/leg. Reports rotatable count, ligand torsion ΔS (locked 0 vs
  sampled >0), ΔH component split, pocket residues, ΔSASA via the
  `thermo_sasa.sasaBurial` helper (stride 10, CG bead radius).
- `tests/test_all.js` SLOW_SUITES gains the suite (SLOW 34 → 44;
  full `--slow` 231 + 44 = 275); FAST/gate baseline untouched (still 231;
  `scripts/wikiskill_gate.js` unchanged — FAST-only gate).
- `src/analysis-panel.js` thermo handler surfaces the rotatable count in
  the ligand-note line (`(ligand … [N rotatable]; pocket …)` via NEW
  `selectedRotatableCount` + `findRotatableBonds` import) — additive
  display only (EPE → 4, BNZ → 0, all → 4); no physics change, no new ids.

Numbers (live, seeded 501/1501, twice identical):
- EPE **ΔS_lig 0.00305 NONZERO** end-to-end (`4 auto rotatable bonds
  (4 torsions) [4 rotatable]`) vs locked 0.00000 vs **BNZ 0.00000 exactly**
  (`rigid ligand (0 rotatable bonds) ⇒ 0 [0 rotatable]` — zero control).
- EPE **ΔH −2.58 ± 0.13 (LJ −0.82 / Coul +1.68 / HB −0.36 / desolv −3.08)**
  (BNZ control −3.38 ± 0.02: LJ −0.83 / Coul 0.00 / HB 0.00 / desolv −2.55;
  EPE Coul +1.68 is the surface-site signature).
- EPE pocket 7 @8Å: HIS31 LEU32 VAL103 PHE104 GLN105 MET106 GLY107
  (BNZ cavity 14: ILE78 LEU84 … LEU118).
- EPE **ΔSASA 429.6 ± 11.2 Å²** (dLig 430.0 ± 10.9 + dProt −0.3 ± 2.8;
  free-ligand ⟨S⟩ 1372.9 Å²; 25+25 frames, stride 10, ~30 ms).
- What it proves: the rotbond autoTorsions path is live end-to-end —
  sampled flexible-ligand ΔS > 0 while the rigid control stays exactly 0.

Verification at close-out: gate OPEN before (231 PASSED, 99 ids) and
after (231 PASSED, 99 ids); `node tests/test_all.js` → 231/231 fast
(~17 s, unchanged); `validate_flexlig` 10/10 twice bit-identical (~0.6 s);
`node --check` clean on all 3 touched/new files (validate_flexlig,
test_all, analysis-panel). No new top-level panels (8 unchanged); no new
DOM ids (all 99 ⊆ index.html, `src/ui.js` untouched).

## 20. Stage-6 note (2026-09-13): heavy-thermo smoke/full split — unblock --slow CI

Pain (open debt from §§13/19): the SLOW tier carried the FULL heavy
protocol (3 reps × 300 equil + 1500 production steps stride 2 → 750
frames/leg, 6 legs × 1800 steps) at ~170–200 s, so every `--slow` CI run
paid the manual-depth cost. FAST stayed 231 but SLOW was CI-hostile.

What shipped (additive only, zero deps, no UI changes, FAST runtime
unchanged):
- `scripts/test_thermo_heavy.mjs`: NEW SMOKE mode via `--smoke` flag or
  `HEAVY_SMOKE=1` env (either triggers; default with no flag/env stays
  FULL). SMOKE = **1 replica × (50 equil + 300 production steps, stride
  2 → 150 frames/leg = 0.30 ps production + 0.05 ps equil at heavy dt
  0.001 ps = 1 fs)**. Same 13 asserts (setup 5 + pipeline 8:
  finiteness/boundedness only — sign reported, NEVER asserted); frame-count
  assert tracks `STEPS/STRIDE` so SMOKE expects 150/leg, FULL 750/leg.
  Header/setup logs print `[SMOKE]`/`[FULL]` + protocol + ps; final line
  reads `test_thermo_heavy (SMOKE|FULL)`. Cα message now reports
  frames/DOF neutrally (no stale ≥10 claim). FULL path byte-logic
  untouched (same seeds 1001/2002/3003, same EQUIL/STEPS/REPS when
  unflagged).
- `tests/test_all.js`: NEW `SLOW_FULL` (`--slow-full` | `--long` |
  `--full` | `SLOW_FULL=1`; implies SLOW) + `runSuiteFile` forwards
  `suite.args`/`suite.env`. Heavy entry carries `args: SLOW_FULL ? [] :
  ["--smoke"]`, so `--slow` = SMOKE heavy, `--slow-full`/`--long` = FULL
  heavy. Tier label prints `SLOW-SMOKE` vs `SLOW-FULL`. Skip message +
  summary comment updated. Assert counts identical either way (heavy 13,
  SLOW +44, totals 231/275).
- `scripts/wikiskill_gate.js`: comment refreshed (FAST 231 baseline,
  SLOW 44, `--slow` total 275); baseline check unchanged
  (FAST-only gate stays 231).

Numbers (live, this stage):
- SMOKE `--smoke`: **13/13 in 12.4 s** (wall 12 s): rep0 −TΔS full
  **−20.97** | bb −5.49 | sc −15.48 | Cα −3.14 | ΔH **−17.02** (LJ −7.73 /
  Coul 0.00 / HB 0.00 / desolv −9.29); frames/DOF 0.5 UNDER-SAMPLED by
  design (pipeline smoke only — verdict still NOT-RESTORED, §15).
  Deterministic twice via `--smoke` + once via `HEAVY_SMOKE=1` —
  bit-identical all three runs.
- FULL (default, rerun this stage): **13/13 in 205 s** (3m24s):
  [22.63, −24.27, 13.96] → mean 4.11; bb 8.96 / sc −4.85 / Cα 1.25;
  ΔH mean −19.73 — replays the Stage-2 record (§9) exactly.
- Tiers:

| Tier | How to run | Suites (asserts) | Runtime measured |
|---|---|---|---|
| FAST (default, gate/CI) | `node tests/test_all.js` | Tier-0 (32) + FAST (199) = **231** | **~16 s** (weakint ~15 s dominates, unchanged) |
| SLOW-SMOKE (fast CI) | `node tests/test_all.js --slow` or `SLOW=1` | FAST 231 + thermo (7, 3.5 s) + heavy SMOKE (13, 12.2 s) + calibration (14, 2.4 s) + flexlig (10, 0.6 s) = **275** | **~36 s** total (FAST 16.9 + SLOW 18.7) |
| SLOW-FULL (manual) | `node tests/test_all.js --slow-full` or `--long` | same 275 asserts, heavy FULL | **~230 s** (FAST ~17 + FULL heavy ~205 + rest ~7) |

How to run each: fast `node tests/test_all.js`; slow-smoke `node
tests/test_all.js --slow`; slow-full `node tests/test_all.js
--slow-full` (or `--long`); heavy alone SMOKE `node
scripts/test_thermo_heavy.mjs --smoke` (or `HEAVY_SMOKE=1 node
scripts/test_thermo_heavy.mjs`), FULL `node
scripts/test_thermo_heavy.mjs`.

Verification at close-out: gate OPEN before (231 PASSED, 99 ids) and
after (231 PASSED, 99 ids); `node tests/test_all.js` → 231/231 fast
(~16 s, unchanged); `--slow` → 275/275 in 36.1 s (SLOW-SMOKE tier 18.7
s, heavy-smoke leg 12.2 s); FULL heavy rerun 13/13 in 205 s (untouched
path re-verified live, not cited); `node --check` clean on all 3
touched files (test_thermo_heavy, test_all, wikiskill_gate). No new
top-level panels (8 unchanged); no new DOM ids (`src/ui.js` untouched).

## 21. Stage-7 note (2026-09-13): thermo cancel UX + PMF 0-hill hint + wiki-path fix

Three open UX debts from §§15–20 closed (additive only, zero deps, no new
top-level panels, Digit1-7 contract intact, every new id in `src/ui.js`):

(a) Thermo Cancel button — the Stage-5 generation counter
(`invalidateThermo`, `src/analysis-panel.js:185`) could cancel a stale
chunked run, but the UI exposed no Cancel control (button stayed disabled
for the whole apo+SASA run with no way to abort). NEW
`<button id="thermoCancelBtn">` in the thermo subpanel row
(`index.html:307`, same `.row.btns` pattern, ships `disabled`) +
registration (`src/ui.js:40`) + `setThermoRunning(running)` helper
(`src/analysis-panel.js:197`, headless-safe with `getElementById`
fallback) + click handler (`:210-214`: `invalidateThermo()` so stale
`stepChunk`/SASA continuations return early, idle button state, caption
`cancelled — thermo run cancelled by user.`). Cancel is enabled only
while a run is in flight (`:333` apo leg start, `:382` SASA leg) and
disabled otherwise (load idle, all three completion paths `:416,426,433`,
`invalidateThermo` itself `:188` so rebuild also disarms it). Physics and
defaults unchanged — display/control only.

(b) PMF 0-hill hint — the BindViz PMF canvas drew a dead-end
(`no metadynamics hills deposited yet`) with no convergence guidance when
the BindLog carried events but zero funnel hills. NEW exported
`PMF_NOHILL_HINT` (`src/capture/bindviz.js:229`):
`no metadynamics hills — enable Funnel bias in Dynamics → Advanced
sampling, run ≥ 50 hills for PMF convergence; timeline/energy live
regardless`. Canvas draws the two-line form (`:259-262`, follows the
`pmf-panel.js` empty pattern); `drawBindviz` reuses the existing
`bindvizCaption` for the single-line form when `hasData && nHills === 0`
(`src/main.js:150-167`, hill count null-safe, no new panels/ids). The
≥50 bar matches the PMF panel convergence gate (`src/pmf-panel.js:81`,
`nHills<50` collecting) — a display threshold, not a measured convergence
claim. Null bindlog still takes the guarded D1 empty state.

(c) Wiki-path fix — the gate verdict named `wiki/evolution/skill-impact.md`
(relative, no such path from root) instead of the real
`.wikiskill/wiki/evolution/skill-impact.md`. Fixed in the gate message
(`scripts/wikiskill_gate.js:92`) and both stale refs in
`.wikiskill/README.md:24,28`. `docs/BINDING_LOOP2_DONE.md:413` already
named the real path.

Verification at close-out: gate OPEN before (231 PASSED, 42 files clean,
99 ids, stale path in message) and after (231 PASSED, 42 files clean,
**100 ids**, correct path in output). `node tests/test_all.js` → 231/231
fast; headless Stage-7 check (`/tmp/opencode/stage7_verify.mjs`) 15/15 —
cancel pair idle/reset, start→cancel caption + button states, 0-hill
canvas hint text, null-bindlog guard. `node --check` clean on all 5
touched JS files; page + all 6 touched servables 200. One mid-session
gate CLOSED flake (bindlog-integration HB temperature trial 13/1 —
`0/30 negative HB samples at the exact native reference gate`;
untouched physics, same flake class as §16) — failing suite alone 14/14,
 `test_all` 231/231, gate rerun OPEN. No new top-level panels (8
unchanged); only +1 DOM id (`thermoCancelBtn`, ⊆ index.html).

## 22. FP1 note (2026-09-13): first-run UX — sample autoload, guided checklist, empty-state audit

Pain: science-complete app, but a new user landed on an empty viewer and
blank summaries with no path forward (Load → Build → Run → Analyze
undiscoverable; `structSummary`/`selSummary`/`networkInfo` rendered blank).

(a) Sample autoload — one-click `Load 4W52 sample` button (`index.html`
Structure panel, id `sampleBtn` ⊆ `src/ui.js`). Reuses the existing
`fetchPdb` path (local `./4w52.pdb` first, offline OK; RCSB/PDBe fallback)
via the same preset mechanism as the `data-ex` links — no new blobs, no
autoload on boot, so the default view is unchanged for returning users;
emphasized (outline) only while empty (`updateGuide`, `src/main.js`).

(b) Guided checklist — ONE collapsed subpanel in Structure (no new
top-level panels; 8 unchanged, Digit1-7 intact) with 4 live steps
(`guideStepLoad/Build/Run/Analyze`, all ⊆ `src/ui.js`). Pure state machine
`fp1GuideState` (`src/ui.js`) maps existing state
(parsed/built/steps+time/frames) to ✓/○; `updateGuide` (`src/main.js`)
paints on build/run/stop clicks + ≤1 Hz `guideTick` in both tick branches
(catches Analyze/frames with no new hooks).

(c) Empty-state audit — every key canvas/caption now actionable when empty
(viewer overlay + `#hud` point at the 1-click sample; timeline/scrub,
`recStatus`, CV/E strips gained next-step hints; PMF/thermo/DCCM/BindViz/
analysis captions already actionable, kept); fixed 3 dead/blank states
(`structSummary`, `selSummary`, `networkInfo` was `""`).

Verification at close-out: gate OPEN before (231 PASSED, 42 files clean,
100 ids) and after (231 PASSED, 42 files clean, **105 ids**).
`node tests/test_all.js` → 231/231 fast; headless FP1 check
(`/tmp/opencode/fp1_guide_verify.mjs`) 16/16 — full Load→Build→Run→Analyze
plus rebuild-reset, run twice. `node --check` clean on all 4 touched JS
files; serve smoke: page + 4 modules + `4w52.pdb` 200, 5 FP1 ids present,
8 top-level panels intact. Additive only, zero deps, defaults unchanged.

## 23. FP2 note (2026-09-13): input robustness — failure classes with next-step copy

Pain (open gap from FP1): input robustness — bad PDB text, missing/empty
ligand, unparsable MOL2, clash-blast placement, and oversized systems threw
raw `err.message` dumps (or blank states) instead of actionable messages.

What shipped (additive only, zero deps, JSDoc; defaults unchanged; no new
top-level panels — 8 unchanged, Digit1-7 intact; zero new DOM ids, so
`src/ui.js` untouched and the DOM contract stays 105):

- NEW `src/input_errors.js` (central validator, dependency-light: imports
  only the existing `pdb.js`/`mol2.js` parsers): 9 frozen failure classes,
  each message = what happened + the exact next click (every message carries
  an imperative "Click" verb); mapped errors ARE `Error` instances
  (`.code`/`.nextStep`/`.technical`), so existing `catch (err)` paths needed
  no control-flow changes. `classifyInputError` is idempotent + total (never
  throws, GENERIC fallback); `formatInputError` renders
  `⚠ <actionable> (detail: <technical>)` — the single-line
  summary/caption surfaces cannot host a `<details>` disclosure, so the
  parenthetical is the collapsed-equivalent (primary copy stays actionable).
  Pure guards: `validatePdbText` (EMPTY_PDB/NO_ATOM pre-check),
  `safeParseCa/safeParseMol2/safeParseLigands` (never throw uncaught on
  garbage — structured `{ ok, data, error }`), `checkSystemSize`,
  `checkPlacement`. Thresholds (validators only, no physics):
  `CLASH_RESIDUAL_THRESHOLD = 2.0` (unitless 1/minRatio; converged ⇒ ≤ ~1.18,
  so > 2.0 = worst contact below 50% vdW sum) and the interactive state
  limit `MAX_HEAVY_ATOMS = 20000` / `MAX_CA_BEADS = 5000` (~15×/30× headroom
  over 4W52's 1308/164 — never fires on bundled systems).
- Wiring (`src/main.js`: load pre-check + fetch/file/MOL2/build catch paths
  now `formatInputError(classifyInputError(err))`, oversized-system guard
  after selection; `src/ligand-panel.js`: all place guards + try/catch round
  `placeLigand`, degenerate/over-threshold poses map to CLASH_HIGH while the
  pose-kept behavior is unchanged) — reuses existing
  structSummary/selSummary/mol2Info/ligPlaceInfo/hud surfaces only.
- NEW `tests/test_input_errors.js` (75 asserts, ~0.1 s, deterministic):
  frozen classes + thresholds, every class message/nextStep contains Click,
  raw-throw → class mapping, classify totality (null/undefined/number/
  object/novel), format primary/secondary contract, `validatePdbText`,
  garbage matrix (safe parsers never throw), size/placement boundaries,
  live 4W52 (164 Cα) + benzene.mol2 (6 atoms/6 bonds) + 4W52 ligands (2).
  `tests/test_all.js` FAST gains the suite (274 FAST, total 306);
  `scripts/wikiskill_gate.js` baseline 231 → 306.

Failure-class table (input → message → next step):

| Input | Code | Message (what) | Next step (exact click) |
|---|---|---|---|
| empty/null/non-string PDB text | EMPTY_PDB | Empty input — no PDB text to parse | Click “Load 4W52 sample (1 click)” or Fetch a PDB ID |
| header-only/prose, no ATOM/HETATM | NO_ATOM | Not a PDB structure | Click “Load 4W52 sample (1 click)” or drop a valid .pdb |
| HETATM-only (no Cα) | NO_CA | No Cα backbone to coarse-grain | Click “Load 4W52 sample”, or heavy mode → Build System |
| empty/unparsable MOL2; place with no ligand | LIGAND_PARSE_FAIL | No usable molecules | Click Ligand MOL2 file → valid .mol2, or Place in Pocket (Auto) |
| place with no structure/selection; bad pocket center | NO_POCKET | Protein selection empty | Click Build System (≥3 Cα), then Place in Pocket (Auto) |
| residual > 2.0 / unconverged / degenerate pose | CLASH_HIGH | Steric overlap remains | Click Place in Pocket (Auto), or Random Surface |
| nHeavy > 20000 or nCa > 5000 | SYSTEM_TOO_LARGE | Exceeds interactive state limit | Click Model & Selection → restrict Chains/range → Build System |
| selection < 3 Cα beads (bonus class) | SELECTION_EMPTY | Need ≥3 Cα for bonded terms | Click Model & Selection → clear filter → Build System |
| anything unrecognized | GENERIC | Failed, detail below | Click “Load 4W52 sample (1 click)” to restore known-good |

Verification at close-out: gate OPEN before (231 PASSED, 42 files clean,
105 ids) and after (**306 PASSED**, 43 files clean, 105 ids);
`node tests/test_all.js` → 306/306 fast (~14 s, new suite 75/75 in 0.0 s);
garbage matrix proof (live `node --input-type=module` check): empty →
EMPTY_PDB, garbage prose → NO_ATOM, header-only → NO_ATOM, HETATM-only →
NO_CA, bad MOL2 → LIGAND_PARSE_FAIL, clash 5.0 → CLASH_HIGH — all mapped
(⚠ + Click), zero raw throws; `node --check` clean on all 6 touched/new
files; serve smoke: page + `input_errors.js`/`main.js`/`ligand-panel.js`/
`ui.js` + `4w52.pdb` + `benzene.mol2` all 200, 5 FP1 ids present, 8
 top-level panels intact. No new top-level panels; no new DOM ids.

## 24. FP3 note (2026-09-13): flake root-cause, browser matrix, a11y/keyboard audit

Three QA debts from §§16/21 closed (additive only, zero deps, defaults
unchanged; no top-level-panel changes — 8 unchanged; no new DOM ids, so
`src/ui.js` untouched and the DOM contract stays 105):

(a) Flake fix — root cause, not retry-hiding. The §§16/21 flake
(bindlog-integration HB trial 13/1: `0/30 negative HB samples`, then 14/14
alone) is unseeded-Langevin tail risk: at the exact native reference the
directional-HB gate reads 0 by construction (physical), so negative HB
samples appear only on thermal excursions during the 300-step run — usually
~20/30, once 0/30. Same stochastic class as the `test_all.js` §8
temperature block, which hid it with an up-to-10 unseeded-trial
keep-first-in-band retry loop. Both now seeded (thermo-SEEDS family,
`src/integrator.js` opt-in `{ seed }`, default path bit-identical):
- `scripts/test_bindlog_integration.mjs:32-45`: single trajectory seeded
  (`BINDLOG_SEED`, default 101, `BINDLOG_SEED_BASE` env override mirroring
  `THERMO_SEED_BASE`). Seed survey (10 seeds) all gave ≥5/30 negatives;
  seed 101 locks a representative 17/30. The `hbSeen > 0` bar is UNCHANGED
  (no lowered standard); assert count stays 14/14 so the FAST gate
  expectation is untouched.
- `tests/test_all.js:300-324`: retry loop replaced by 3 seeded replicas
  (`TEMP_SEEDS = [101, 202, 303]`) with a replica-mean assertion in the same
  260–340 K band. Physics reason: the OU thermostat is exact in
  distribution, but one 200-step kinetic-T sample carries O(1/√steps) noise
  while the ensemble mean sits at the bath. Measured: reps
  [299.9, 280.8, 330.4] → mean 303.7 K. Assert count unchanged (2 in §8).
- 5-run record (final tree, consecutive): 14/14, 14/14, 14/14, 14/14, 14/14
  (bit-identical 17/30 negatives each run); override check
  `BINDLOG_SEED_BASE=202` → 14/14 (18/30, still passing with margin).

(b) Cross-browser matrix — NEW `scripts/smoke_browsers.mjs` (node:http
static server + dev-playwright, headless): per browser load → open the
"Getting started" subpanel (real-user path; Playwright treats closed-details
content as hidden for clickability) → Load 4W52 sample (local `./4w52.pdb`
first, offline OK) → build (STATE: READY) → Run ~3.5 s → assert step > 0,
HUD live `t =`, no ⚠ banner, zero uncaught page errors.

| Browser | Version | Steps | HUD | Page errors | Verdict |
|---|---|---|---|---|---|
| Chromium (headless) | 151.0.7922.34 (playwright chromium-1234) | step 12115, t = 21.2 ps | live, clean | 0 | PASS |
| Firefox (headless) | 153.0 (playwright firefox-1538) | step 11400, t = 19.9 ps | live, clean | 0 | PASS |
| Safari | n/a (no WebKit engine on linux; no Safari binary, no playwright webkit build) | — | — | — | STATIC-ONLY: 200 on `/index.html`, `src/main.js`, `src/ui.js`, `src/integrator.js`, `src/seeded-rng.js`, `4w52.pdb`, `benzene.mol2`, `css/style.css`; real run owed on macOS |

(c) A11y/keyboard audit — all 30 buttons already named (text or
`aria-label="Close settings"`); Tab order is natural DOM (What-changed →
Settings → panel summaries → inputs → viewer → canvas); the Digit1-7
`isEditing` (INPUT/SELECT/TEXTAREA) guard already existed. Gaps fixed:
- Visible focus: NEW `:focus-visible` accent outline for
  btn/a/summary/input/select/#canvas (`css/style.css:149-157`; mouse clicks
  stay ring-free, Tab always shows).
- Canvas semantics (`index.html`): `#canvas:376` gains `tabindex="0"`
  `role="img"` + viewer `aria-label`; 6 data canvases
  (`:239,303,345,361-363`) + 2 dock strips (`:410,412`) gain `role="img"` +
  `aria-label` (display-only, no new tab stops).
- Tiny JS guards (`src/main.js:893-895,926` + same-shape fallback
  `index.html:517`): Space no longer hijacks focused BUTTON/A/SUMMARY
  (native activation wins; Space on a focused Run button still toggles via
  its native click); Digit index derives from `e.code` (layout-independent;
  `e.key` yields symbols with Shift held).
- Keyboard spot-check (temp harness, both browsers, 16/16): Tab order sane
  incl. canvas reachability; 2px focus ring on Tab focus; Digit1 while typing
  leaves panels untouched and types "1"; Space on focused Fetch keeps native
  behavior; Digit1 outside inputs toggles panel; zero page errors.

Verification at close-out: gate OPEN before (306 PASSED, 43 files clean,
105 ids) and after (**306 PASSED**, 43 files clean, 105 ids);
`node tests/test_all.js` → 306/306 fast (seeded §8 mean 303.7 K);
integration 5× 14/14 above; smoke 2/2 PASS above; `node --check` clean on
all 4 touched JS files (test_bindlog_integration, test_all, main,
smoke_browsers) + `index.html` parses. No new top-level panels
(8 unchanged); no new DOM ids (`src/ui.js` untouched).

## 25. FP4 note (2026-09-13): export/persistence — matrix audit, session save/load, roundtrip

Pain (open gap from FP1–FP3): export/persistence — the BindLog BLG1 blob
existed (`src/capture/bindlog.js:166` `toBinaryBlob`), PMF CSV + thermo
table existed ad hoc (`src/analysis.js:619` `pmfCsv`,
`src/analysis/thermodynamics.js:548` `formatThermoTable`), but no session
save/load roundtrip and no documented export matrix (thermo TXT, BLG1,
DCCM CSV, settings JSON had no one-click path).

What shipped (additive only, zero deps, JSDoc; defaults unchanged; no new
top-level panels — 8 unchanged, Digit1-7 intact; +6 DOM ids, all in
`src/ui.js` ⊆ `index.html`):

(a) Export-matrix audit — every exportable artifact now has a one-click
button inside its existing subpanel, all downloads via the existing
`downloadText`/blob pattern (`src/recorder.js:168`):

| Artifact | Button | Format | Wiring |
|---|---|---|---|
| Trajectory XYZ / PDB | Download (`#dlBtn`, Recording) | `.xyz` / `.pdb` | pre-existing (`src/main.js` `recorder.buildFile`) |
| Trajectory JSON | Download (`#dlBtn` + new `json` `#exportFmt` option) | `.json` (`trajectoryJson`) | `src/main.js` json branch → `src/session.js` `trajectoryJson` |
| PMF | Export PMF (CSV) (`#anaPmfBtn`, PMF & Analysis) | `.csv` (`pmfCsv`) | pre-existing (`src/analysis-panel.js`) |
| Thermo table | Export Thermo (TXT) (`#thermoDlBtn`, thermo subpanel) | `.txt` (cached `formatThermoTable`) | `src/analysis-panel.js` (`_lastThermoText` set on both success paths) |
| BindLog | Export BindLog (BLG1) (`#bindlogDlBtn`, Recording) | `.blg1` (`toBinaryBlob`) | `src/main.js` → `src/session.js` `downloadBlob` |
| DCCM | Export DCCM (CSV) (`#dccmDlBtn`, DCCM subpanel) | `.csv` (`dccmCsv` i,j,C) | `src/analysis-panel.js` (`_lastDccm` cache) |
| Settings | Export Settings (JSON) (`#settingsDlBtn`, settings modal footer) | `.json` (`settingsJson`, sim.* keys) | `src/settings-panel.js` |
| Session | Save Session (JSON) (`#sessSaveBtn`, Recording) + Load via `#sessFile` | `.json` (schema v1, §b) | `src/main.js` + NEW `src/session.js` |

Empty-guarded throughout: thermo/DCCM/BindLog buttons render an actionable
`⚠ … run first, then Export` hint to the existing `analysisOut`/`recStatus`
surfaces when there is nothing to export (no silent empty files;
`trajectoryJson([])` throws `No frames recorded yet`).

(b) Session save/load — NEW `src/session.js` (headless-pure, imports only
`src/input_errors.js`): `buildSession` / `serializeSession` /
`validateSession` / `parseSession` + `dccmCsv` / `parseDccmCsv` /
`trajectoryJson` / `settingsJson` / `downloadBlob` /
`estimateFramesBytes`; `SESSION_VERSION = 1`,
`SESSION_MAX_BYTES = 256 KiB`, `SESSION_FILE_MAX_BYTES = 1 MiB`.
Schema `{version: 1, app, savedAt, pdbId, modelMode, chains, resFrom,
resTo, includeLig, physicsLevel, ligand:{selected}, thermoLig,
settings:{backend, numThreads, solventModel, saltM, epsIn, epsOut,
sasaGamma, respaOn, respaOuterFs, chemicalNetworkOn},
dynamics:{rc, gamma, temp, fric, mass, motionGain, bindPot, holoSprings},
recording:{stridePs, maxFrames, exportFmt},
recorderMeta:{count, spanPs, includeFrames:false, framesOmitted, note}}`.
Load restores settings + pickers + caption (`applySession` in
`src/main.js`: pdbId/model/chains/range/includeLig/physicsLevel (+persist),
ligand + thermo pickers, settingsState + modal inputs, dynamics sliders,
recording inputs, `canvasCaption` session note + `hud` confirm; file input
reset so the same file re-loads). Load validation reuses input_errors
codes: malformed JSON / bad version / bad enum → GENERIC, oversized file →
SYSTEM_TOO_LARGE, empty → GENERIC.
Limits (documented cutoff): full PDB text and full trajectory frames are
NEVER persisted — `recorderMeta` carries counts only (`includeFrames`
always false, no `frames` key); sessions serialize to ~1 KiB against the
256 KiB cap; nothing session-shaped goes to localStorage (files only —
`sim.physicsLevel` remains the single localStorage key).

(c) Roundtrip test — NEW `tests/test_session_roundtrip.js` (46 asserts,
~0.0 s, deterministic, no `Math.random`): save→load restores
physicsLevel/picker/settings (L2/all/cpu/temp-300/exportFmt-json/count-42,
no `frames` key, `< 256 KiB`); guards (empty/malformed→GENERIC, version
99→GENERIC, `L9`→invalid + build-fallback L0, 1 MiB+1 B file→
SYSTEM_TOO_LARGE, `estimateFramesBytes(500,164)` exact); BLG1 blob
parse-back header ok (magic `BLG1`, version 1, `nF=3 n=4`, `nE=11`,
roundtrip 3 frames/11 events); PMF CSV parse-back row count (header
`# T=300, gamma=6, hills=5, V0=1660.54`, 96 data rows = `bins`, footer
`hills,5`); DCCM CSV parse-back (4×4, diagonal 1, off-diag 1e-6);
trajectory JSON + settings JSON parse-back (`sim.physicsLevel` L2).
`tests/test_all.js` FAST gains the suite (320 FAST, total 352);
`scripts/wikiskill_gate.js` baseline 306 → 352.

Verification at close-out: gate OPEN before (306 PASSED, 43 files clean,
105 ids) and after (**352 PASSED**, 44 files clean, **111 ids**);
`node tests/test_all.js` → 352/352 fast (~13 s, new suite 46/46 in 0.0 s);
`node --check` clean on all 8 touched/new files (session, ui, main,
analysis-panel, settings-panel, test_session_roundtrip, test_all,
wikiskill_gate); serve smoke: page + 6 JS modules + `4w52.pdb` +
`benzene.mol2` + `css/style.css` all 200, 6 new ids present, 8 top-level
panels intact. No new top-level panels; defaults unchanged; large-frame
guard holds (sessions ~1 KiB, frames never embedded).

## 26. FP5 note (2026-09-13): heavy-build progress UX + CG-interactive/heavy-offline workflow

Pain (open gap from FP1–FP4): heavy mode costs ~77–85 ms/step (S7 pareto) and
the O(n²) topology build ~30–60 ms on 4W52 (1308 atoms), so a synchronous
heavy Build froze paint with no progress and no cancel. S7 decided ACCEPT CPU
(no GPU port — binding terms stay CPU-only), so the fix is honest progress +
a documented tiered workflow, never a frozen tab.

Blocking points read before the change (`src/main.js`): `buildSystem` built
the heavy FF synchronously (`new HeavyForceField`) and the hot-rebuild path
did the same; the run loop already slices per-frame (`integ.advance(steps,
14)`) with Pause as cancel, so only the build needed chunking. Existing UX
reused: thermo chunked captions + `_thermoGen` cancel
(`src/analysis-panel.js:183-218`), BindViz 1 Hz (`src/main.js` bindvizTick),
`guideStep` checklist (`src/main.js` guideTick).

What shipped (additive only, zero deps, JSDoc; defaults unchanged; no new
top-level panels — 8 unchanged, Digit1-7 intact; +1 DOM id, in `src/ui.js`):

- Chunked topology (`src/heavy.js`): `buildTopology` refactored into shared
  `topologyBondRows` (row-range bond loop) + `findTopologyRings` +
  `finishTopology` (angles/propers/ring-brace tail) with a sync wrapper that
  is bit-identical; NEW `buildTopologyChunked` (`:475`, `chunkRows` default
  128) runs the same rows in slices with `setTimeout(0)` yields, per-slice
  `onProgress`, and cooperative `isCancelled()` polling (throws `heavy build
  cancelled`). `HeavyForceField` takes an optional 4th `opts.topo` (prebuilt
  topo skips the sync O(n²); absent → legacy path bit-identical, all 3-arg
  callers untouched).
- Build wiring (`src/main.js:540-600,680-778`): CG path stays fully
  synchronous (fast, no flicker). Heavy prep (ligands/protonation) stays
  sync, then `runHeavyBuildAsync` chunks topology (11 slices on 4W52) with
  `Building heavy… topology d/n rows` captions to the reused `#selSummary`,
  one short sync tail (FF assemble + first-eval integrator, each < ~500 ms),
  then the shared `finishBuildCommon()` finalize (extracted verbatim — CG
  and heavy share one finalize). Cancel is the thermo gen-counter pattern:
  `invalidateHeavyBuild()` bumps `_heavyGen`; stale slices/continuations
  return early without touching state/DOM. `setHeavyBuilding` disables
  Build/Run/Reset during work and re-enables after (gen-guarded); NEW
  `Cancel build` button (`index.html:112`, id `heavyCancelBtn` ⊆
  `src/ui.js:19`, ships disabled, armed only mid-build).
- NEW `src/heavy_progress.js` (headless-pure, zero deps):
  `HEAVY_TOPO_CHUNK_ROWS = 128` (`:19`), `HEAVY_BUILD_CANCELLED` (`:22`),
  `setHeavyButtons` (`:31`), `setHeavyCaption` + `heavyTopoCaption` (`:56`).
- Choice documented: run slicing stays `advance(14)` (a single heavy
  `ff.compute` is atomic — finer slicing infeasible without touching the
  zero-alloc grid kernels; each step < 500 ms) with HUD 10 Hz + Pause as the
  run progress/cancel; slider hot-rebuild stays sync (single slice, benign
  ~100 ms race: build-start params win, retune after READY). No modal — the
  chunked captions + disabled buttons are the progress surface (P1–P5).
- Workflow doc (`docs/WORKFLOWS.md` §"CG-interactive / heavy-offline"):
  L0/L1-interactive vs L2-heavy-offline table (ms from
  `docs/pareto_frontier.csv` 2026-09-12 FULL: L0 0.237 / L1 0.235 /
  L1+BindLog 0.216 / L2 84.8 / L3 83.1 / L2-full 85.5; ratio ≈ 400:1),
  recommended budgets (CG minutes + stride 2 ps/200–500 frames; heavy in-tab
  ≤ 1000 steps/≤ 100 frames; 1 ns heavy ≈ 1 day → headless scripts only),
  export-then-analyze path (FP4 one-click matrix → `notebooks/`).

Verification at close-out: gate OPEN before (352 PASSED, 44 files clean,
111 ids) and after (**352 PASSED**, 45 files clean, **112 ids** — only
+heavyCancelBtn); `node tests/test_all.js` → 352/352 fast (untouched, FAST
baseline stays 352); headless FP5 proof (`/tmp/opencode/fp5_heavy_proof.mjs`)
18/18 — sync-vs-chunked topo bit-identical (1507 bonds / 1789 angles / 4160
propers), 11 progress slices, heartbeat 5× mid-build (vs 0 sync), max slice
6.4 ms (< 500), cancel throws + gen-counter newest-only finalize,
button-pair transitions on stubs, `{topo}` FF energy bit-identical
(−2106.29), single heavy step 14.9 ms (< 500); `node --check` clean on all
touched files; serve smoke: page + all modules + `4w52.pdb` 200, new id
present, 8 top-level panels intact.

## 27. FP6 note (2026-09-13): methods/validation table + 1CRN negative control

Pain (open gap from FP1–FP5): numbers scattered across §§8–21, R-docs,
calibration, and the pareto CSV with no single methods/validation table;
and every system was 4W52/T4L (BNZ rigid + EPE flexible) — no
third/negative control, so the pipeline had never been shown to return
null on a ligand-free system.

What shipped (additive only, zero deps, analysis + docs, NO retuning;
no UI changes — no new ids, `src/ui.js` untouched):

(a) `docs/VALIDATION.md` (new, concise): trust boundary WITH numbers up
front (ranking-only, NOT FEP: single-rep ΔG_est underbinds 2–3 kcal/mol
with bootstrap ±0.16 ⊕ replica-SD ±7.1) + one methods table (system,
protocol seeds/steps, observable, computed ± error,
experimental/literature, verdict) covering 4W52-BNZ record
(−6.82 ± 0.16) + BNZ-only true cavity (−3.64 ± 0.08, EPE inflation −3.18
reported not corrected) + real SASA (ΔSASA 167.1 ± 3.5 → solvent +2.01 ±
0.04 → revised ΔG −0.16 vs ITC −5.20 ± 0.20 / NMR −4.20 ± 0.10) + EPE
flexible (ΔS_lig 0.00305 NONZERO vs BNZ exactly 0) + 1CRN null + heavy
SMOKE/FULL/LONG + ala-scan (floor 0.05, all `~noise`, NO CG resolution)
+ Schlitter 0.20% / torsion kB·ln12 unit anchors + bench tiers (L0:L2 ≈
400:1, Loop-2 ≈1.0× on heavy). Run-matrix + script refs included.

(b) Negative control `scripts/validate_1crn_null.mjs` (new, 8 asserts,
SLOW tier, seeded 701/1701 with `NULL_SEED_BASE` override, <1 s):
1CRN crambin (46 Cα, 0 HETATM ligands — ligand-free by parse) apo-vs-apo
CG (both legs binding-off, record protocol T300/ζ8.0/mass-110, 500 steps
stride 2 → 250 frames/leg; zeros binding vectors by construction; pocket
= COM-8Å 16 residues, no ligand COM exists). Asserts are
finiteness/boundedness only — sign never asserted, zero never forced.
`tests/test_all.js` SLOW_SUITES gains the suite (SLOW 44 → 52, full
`--slow`/`--slow-full` 396 → 404); FAST/gate baseline untouched (352;
`scripts/wikiskill_gate.js` comment refreshed, check unchanged).

Numbers (live, seeded 701/1701, three runs bit-identical):
- 1CRN **ΔH 0.00 ± 0.00** (|ΔH| < 0.5 null-scale bound vs 4W52 signal
  −3.6/−6.8) — null returns null on the energy channel.
- 1CRN **−TΔS_pocket +5.67** (bounded < 20, test_thermo bar) — honestly
  nonzero: apo-vs-apo Schlitter noise floor (S_holo 1.429 vs S_apo 1.448),
  same variance class as the §15 LONG finding, never a forced zero.
- f/DOF 5.2 UNDER-SAMPLED by design (pipeline null-smoke only).

Verification at close-out: gate OPEN before (352 PASSED, 45 files clean,
112 ids) and after (352 PASSED, 45 files clean, 112 ids);
`node tests/test_all.js` → 352/352 fast (~13 s, unchanged);
`validate_1crn_null` 8/8 three runs bit-identical (<1 s);
`node --check` clean on all 4 touched/new files (validate_1crn_null,
test_all, wikiskill_gate, VALIDATION.md is docs-only). No new top-level
panels (8 unchanged); no new DOM ids.

## 28. FP7 note (2026-09-13): release closeout — CHANGELOG, version, license/citation, final gate

Gap closed: Loop-2 → followup1-7 → FP1–FP6 shipped science + UX with no
CHANGELOG/version, an unchecked license/citation surface, and no final
gate record. Docs/metadata only — zero code-behavior changes (one
display-only version-string bump; defaults, physics, panels untouched).

(a) BEFORE baseline (working tree, pre-FP7): gate OPEN (45 files clean,
352 PASSED / 0 FAILED ≥ 352 baseline, 112 ui ids, Digit1-7 vs 8 panels);
FAST `node tests/test_all.js` → 352/352 in ~13.6 s. No `package.json` in
this repo (confirmed by glob) — the version flow IS `src/version.js`
(HUD prefix `src/main.js:1491`, `REMARK` provenance
`src/recorder.js:117`, console badge); README had usage but no
version/license/citation section; `LICENSE` (MIT, © 2026 Samir Rana) and
`CITATION.cff` (`v1.0-jpcb`, Zenodo placeholder) present.

(b) What shipped (additive, docs/metadata only, zero deps):
- NEW `CHANGELOG.md`: releases phases1-5 → Loop-2 → stages1-7 →
  followup1-7 → FP1–FP6 → 1.1.0-fp7 with commit hashes
  (`7de5505`…`173ccad`; FP1–FP6 working tree), FAST lineage
  32→179→215→231→306→352, SLOW lineage 199→245→249→275→404, one-line
  key numbers per release, and the FP7 triple gate record.
- Version `1.0.0-transform` → **`1.1.0-fp7`** (`src/version.js:13`,
  BUILD_DATE → 2026-09-13): minor, not patch/major — Loop-2 + FP1–FP7
  were backward-compatible additive-only (no breaking changes → no
  major; seven finish-product increments → more than a patch). Display
  only: HUD/REMARK follow automatically; static dock-footer line in
  `index.html` (text only, zero new DOM ids, no panels touched).
  `CITATION.cff` `v1.0-jpcb` archival tag intentionally untouched.
- License/citation: README gains `## Version, license & citation`
  (version pointer + MIT link + CFF/BibTeX/VALIDATION pointers +
  ranking-only boundary); secret grep over `src scripts tests ml
  index.html CITATION.cff` (api-key/secret/password/private-key/token
  patterns) → zero hits.

(c) VERIFY (all 2026-09-13, final tree): gate AFTER **OPEN** (45 files
clean, 352 PASSED / 0 FAILED, 112 ids — footer adds zero ids); FAST
**352/352 in ~13.4 s**; SLOW-SMOKE `--slow` **404/404 in ~29.5 s**
(FAST 352 + SLOW 52: thermo 7 + heavy SMOKE 13 + calibration 14 +
flexlig 10 + 1crn-null 8); `node --check` clean (`version.js`,
`index.html` parses); serve smoke `/`, `src/main.js`,
`src/version.js`, `4w52.pdb`, `CHANGELOG.md` all 200, footer
`v1.1.0-fp7` string present in served HTML. No new top-level panels
(8 unchanged); no new DOM ids (`src/ui.js` untouched).

## 29. Revolution-2 Issue-5 note (2026-09-18): deterministic coincident escape + placement-policy docs

Pain: the coincident degenerate branch in `clashGrad` (`src/placement.js:89`)
pushed a fixed `+x+y+z` diagonal (`Fx += 1, Fy += 1, Fz += 1`, `minRatio = 0`)
— brittle when `+x+y+z` is walled in a dense hetero shell (net push walks
into the wall while the step halves to `<1e-10`), and `minRatio = 0`
reported `residualClash = Infinity` (non-finite degenerate pose) for a
finite clash energy. Docs were stale on the
`ligandStart`/`excludeFrom`/`bindingTermsActive`/live-terms/`rmsdLig`
policies and carried pre-rev1-issue2 line numbers.

What shipped (placement escape branch + docs + focused test only; no
caller, kernel, UI, or gate changes):
- `src/placement.js:89-104` (escape branch only): deterministic
  index-hashed escape direction per `(ligand a, collider i)` with its torque
  arm (rotation joins the escape), replacing the fixed diagonal. Fixed given
  indices, so placement stays deterministic given `seed` (seed still owns
  the initial SO(3) rotation, `src/placement.js:246`); `minRatio` floor
  (`1e-3/rE`) keeps `residualClash` (`src/placement.js:220`) finite.
- Policies (re-stated with correct line refs; see `docs/PLACEMENT.md`):
  placeInPocket clashes hetero-inclusive (`getProteinCoordsAndSigma`,
  `src/ligand-panel.js:116`, `collEnd = ligandStart ?? nProt`) while the
  cavity search stays protein-only (`findPocketCenter(protein.pos,
  state.ff.nProt)`, `src/ligand-panel.js:251`); slot exclusion via
  `opts.excludeFrom` / `protein.excludeFrom` (`src/placement.js:157`);
  hetero-excluded viewer pocket (`src/viewer.js:211`) and H-bond overlay
  (`src/viewer.js:505`); RMSD split (`rmsd` protein-only `src/heavy.js:1395`,
  `rmsdLig` `[ligandStart, n)` `src/heavy.js:1413`, legacy `rmsdAll`);
  live per-term mirror without BindLog (`liveTermsWanted` et al.,
  `src/main.js:571`) with the wiring descriptor `bindingTermsActive`
  (`src/ff-binding.js:51`).
- NEW `tests/test_rev2_issue5_placement_escape.js` (14 asserts, <1 s,
  deterministic; NOT wired into `tests/test_all.js`, gate untouched):
  single-atom coincident escape leaves off-diagonal (cos with `+x+y+z` −0.93
  vs old exactly 1.0), still-coincident pose reports finite residual 3400
  (old: `Infinity`) with honest `converged:false`, dense hetero cage
  (octahedral ±2.5 A + coincident center) converges finite + bit-identical
  twice, `placeLigand` same-seed bit-identical. Fails 2/14 on the pre-fix
  tree (diagonal cos = 1.0, residual `Infinity`), passes 14/14 after.

Verification at close-out: `node --check src/placement.js` clean; new
suite 14/14 twice bit-identical; `node tests/test_placement_hetero.js`
16/16 (hetero-aware collision set unaffected — no coincidence there);
`git diff --check` clean. No new top-level panels; no new DOM ids.
