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
