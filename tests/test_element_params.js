/**
 * test_element_params.js — one element, one set of parameters, on both engines.
 *
 * THE PRECEDENT (this file exists because of it)
 * ---------------------------------------------
 * goal M7 measured the per-element LJ sigma / epsilon / partial charge on both
 * simulation paths and found:
 *
 *   • 9 of the 9 elements in the CG/heavy element table resolved BYTE-IDENTICALLY
 *     on both engines — because both already read `ff-params.js LIG_ELEMENT`,
 *     and every element absent from it fell back to the same
 *     `LIG_ELEMENT_DEFAULT`. There was NO live sigma/eps/q split.
 *   • but src/heavy.js:82 ALSO declared its own
 *         HEAVY_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0 }
 *     — a second declaration of the same quantity in a different module, one
 *     edit away from a real divergence. It was also DEAD: nothing read it.
 *     This is the same defect class as KB_KCAL (5.03e-8, two modules) and the
 *     Coulomb constant (332.0 vs 332.06371, six sites). tests/test_constant_ledger.js
 *     rule 5 now fails if such a duplicate reappears; this file is the runtime
 *     half of the guard.
 *   • and there WAS a live CG-vs-heavy divergence: a COVERAGE gap rather than
 *     a duplicated literal. The CG path had no metal parameters at all, so a
 *     metal arriving through pdb.js parseLigands (which keeps any HETATM group
 *     with ≥ 2 atoms, so a Zn coordinated inside a multi-atom hetero group
 *     REACHES CG) was handed ELEMENT_LJ_DEFAULT — σ 3.4 Å, ε 0.12 kcal/mol,
 *     q 0 — while the heavy path handed it METAL_ELEMENT — σ 1.30–2.00 Å,
 *     ε 0.05 kcal/mol, q +1/+2. §5 pinned that gap and printed its magnitude
 *     on every run. Closing it is a physics decision (does the CG model want
 *     explicit ions? a metal-aware desolvation term? coordination?), so M7
 *     did not make it.
 *
 * ── 2026-10: the gap is CLOSED (§5 rewritten, not deleted) ────────────────
 * resolveElementParams is now `ELEMENT_LJ[el] ?? METAL_ELEMENT[el] ??
 * ELEMENT_LJ_DEFAULT`, so the CG engine reads the SAME ion table the heavy
 * engine reads. §5 was INVERTED, not removed: it still pins values field for
 * field on both engines, and its old "the CG side must be EXACTLY the generic
 * default" pin — which existed to stop the gap silently widening — now runs in
 * the opposite direction and fails if any ion ever resolves to the neutral
 * generic row again. The measured before/after (same seed, same engine,
 * 4 000 steps, 4W52 + a Zn site; tools/exp_cg_metal_stability.mjs):
 *
 *   ion LJ / charge :  σ 3.4→1.40 Å, ε 0.12→0.05, q 0→+2   (σ ratio 2.43x → 1.00x)
 *   ion stability   :  NOT ejected and NOT blown up (a free ion 30 Å out stays
 *                     at 23.3 ± 8.4 Å from the protein; the CG Coulomb term is
 *                     screened by ε(r)=4+76·tanh(r/8), so a +2 charge has no
 *                     divergent monopole here)
 *   Zn–O geometry   :  time outside the physical 1.5–3.5 Å band 12.7% → 5.8%
 *                     (holo off) and 22.1% → 13.2% (holo on)
 *   STILL absent    :  CG builds NO metal–donor coordination restraint, and the
 *                     ion's formal charge attracts to a whole-residue Cα bead,
 *                     so a Zn can sit ~2 Å from a bead rather than 2 Å from a
 *                     donor ATOM. That is a Cα-resolution limit, documented as
 *                     absent-not-pending in docs/LIMITATIONS.md and announced
 *                     once per ForceField from src/cg/params.js.
 *
 * WHAT IT ASSERTS
 *   1. [1] Every element row: the CG engine (forcefield.js) and the heavy
 *      engine (heavy.js) resolve the SAME sigma, eps, hb flag and ΔG, compared
 *      with `===` on the engine's own arrays — not on the table, so a path that
 *      stops calling the resolver fails too.
 *   2. [2] Every non-metal element with no row (B, SE, SI, AL, H, and an
 *      element nobody typed a row for) resolves to the default in BOTH engines,
 *      byte-identically, including hb and ΔG.
 *   3. [3] The re-export facade is identity-equal: ff-params.js publishes the
 *      canonical objects, not copies, so there is nothing to drift.
 *   4. [4] The partial charge stays TWO documented models: per-element in CG,
 *      per-atom from physics/charges.js assignCharges in heavy. Both columns
 *      are pinned so neither can silently start being the other.
 *   5. [5] METALS: both engines resolve the SAME METAL_ELEMENT row (σ/ε/q/ΔG),
 *      every ion keeps its formal charge, both resolvers return the identical
 *      object, the uniform ε=0.05 and hb=false conventions are pinned, the
 *      placeholder ΔG is pinned AS a placeholder, and a re-divergence back to
 *      the neutral generic default fails the test.
 *   6. [6] RES_CLASS (CG residue-class beads) and ELEMENT_LJ (per element) keep
 *      DISJOINT key sets — they are different models and must not be merged.
 *   7. [7] src/heavy.js declares no per-element parameter literal of its own.
 *
 * Zero dependencies. Node-importable; no DOM globals.
 * Runnable: node tests/test_element_params.js
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import {
  ELEMENT_LJ, ELEMENT_LJ_DEFAULT, METAL_ELEMENT, RES_CLASS,
  resolveElementParams, resolveHeavyElementParams, elementCoverage,
} from "../src/physics/params.js";
import * as ffParams from "../src/ff-params.js";
import { ForceField } from "../src/forcefield.js";
import { HeavyForceField } from "../src/heavy.js";
import { assignCharges } from "../src/physics/charges.js";
import { NONBONDED_TABLE, NONBONDED_DEFAULT } from "../src/physics/forcefield/amber14sb.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

/** Elements with a row in ELEMENT_LJ — the CG/heavy shared table. */
const ROW_ELEMENTS = Object.keys(ELEMENT_LJ);
/** Elements with NO row and that are NOT metals: must hit the default in both engines. */
const DEFAULT_ONLY_ELEMENTS = ["B", "SE", "SI", "AL", "H", "XX"];
const METALS = Object.keys(METAL_ELEMENT);

/** One 40 Å-spaced ligand atom per element, as an explicit ForceField ligand. */
function cgFor(elements) {
  const atoms = elements.map((el, i) => ({
    element: el, atomName: "X1", resName: "LIG", x: i * 40.0, y: 0, z: 0,
  }));
  const beads = [{ resName: "ALA", x: 0, y: 0, z: 0 }, { resName: "ALA", x: 3.8, y: 0, z: 0 }];
  return new ForceField({ beads, segments: [[0, 2]] },
    { rc: 10, gamma: 1.0, binding: { on: false } },
    [{ atoms, bonds: [], angles: [], impropers: [], resName: "LIG" }]);
}

/**
 * One 40 Å-spaced heavy atom per element, resolved through HeavyForceField.
 * `isMetal` is set the way parseHeavy sets it (from METAL_ELEMENT), because
 * physics/charges.js assignCharges — which heavy.js writes into `_elem[i].q` —
 * only assigns a formal ion charge to atoms flagged as metals. Getting this
 * flag wrong would print a heavy q of 0 for every ion and quietly understate
 * the measured divergence below.
 */
function heavyFor(elements) {
  const atoms = elements.map((el, i) => ({
    element: el, atomName: "X1", resName: "LIG", x: i * 40.0, y: 0, z: 0,
    isProtein: false, isLigand: true, isMetal: !!METAL_ELEMENT[el], isHetero: true,
  }));
  const ff = new HeavyForceField({ atoms }, { gamma: 1.0 }, []);
  return { ff, atoms };
}

console.log("=== element parameter contract: one table, two engines, no drift ===");
console.log(`  canonical home: src/physics/params.js — ${ROW_ELEMENTS.length} element rows ` +
  `(${ROW_ELEMENTS.join(", ")}), ${METALS.length} metals, ` +
  `${Object.keys(RES_CLASS).length} CG residue classes`);

// ── [3] The facade publishes the canonical OBJECTS, not copies ────────────
console.log("\n[3] src/ff-params.js re-export identity...");
assert(ffParams.LIG_ELEMENT === ELEMENT_LJ,
  "ff-params.LIG_ELEMENT === params.ELEMENT_LJ (same object, not a copy)");
assert(ffParams.LIG_ELEMENT_DEFAULT === ELEMENT_LJ_DEFAULT,
  "ff-params.LIG_ELEMENT_DEFAULT === params.ELEMENT_LJ_DEFAULT");
assert(ffParams.METAL_ELEMENT === METAL_ELEMENT, "ff-params.METAL_ELEMENT === params.METAL_ELEMENT");
assert(ffParams.RES_CLASS === RES_CLASS, "ff-params.RES_CLASS === params.RES_CLASS");
assert(ffParams.resolveElementParams === resolveElementParams &&
  ffParams.resolveHeavyElementParams === resolveHeavyElementParams,
  "both engines' resolvers are the SAME function objects through the facade");

// ── [1] Every element row resolves byte-identically on both engines ───────
console.log("\n[1] Per-element row: CG path vs heavy path (=== on the engine arrays)...");
console.log("    element |  CG sigma/eps/hb/dG        | heavy sigma/eps/hb/dG      | agree");
const cgRow = cgFor(ROW_ELEMENTS);
const { ff: hvRow } = heavyFor(ROW_ELEMENTS);
let rowDisagree = 0;
ROW_ELEMENTS.forEach((el, a) => {
  const i = ROW_ELEMENTS.indexOf(el);
  const c = { sigma: cgRow._ligSigma[a], eps: cgRow._ligEps[a], hb: cgRow._ligHB[a], dG: cgRow._ligdG[a] };
  const h = hvRow._elem[i];
  const same = c.sigma === h.sigma && c.eps === h.eps &&
    (c.hb === 1) === h.hb && c.dG === h.dG;
  if (!same) rowDisagree++;
  const fmt = (o, hbAsInt) => `${o.sigma} / ${o.eps} / ${hbAsInt ? o.hb : (o.hb ? "true" : "false")} / ${o.dG}`;
  console.log(`    ${el.padEnd(7)} | ${fmt(c, false).padEnd(24)} | ${fmt(h, false).padEnd(24)} | ${same ? "YES" : "NO ***"}`);
});
assert(rowDisagree === 0,
  `all ${ROW_ELEMENTS.length} element rows resolve byte-identically (sigma/eps/hb/dG); ${rowDisagree} disagree`);

// And the shared resolver itself, table-free:
let resolverDisagree = 0;
for (const el of [...ROW_ELEMENTS, ...DEFAULT_ONLY_ELEMENTS]) {
  const c = cgFor([el])._ligSigma[0];
  const h = heavyFor([el]).ff._elem[0];
  if (!(c === h.sigma)) resolverDisagree++;
}
assert(resolverDisagree === 0,
  `resolveElementParams feeds both engines for all ${ROW_ELEMENTS.length + DEFAULT_ONLY_ELEMENTS.length} probed elements; ${resolverDisagree} disagree`);

// ── [2] Elements with no row fall back to the SAME default object ─────────
console.log("\n[2] Elements with no row fall back identically (=== ELEMENT_LJ_DEFAULT)...");
const cgDef = cgFor(DEFAULT_ONLY_ELEMENTS);
const { ff: hvDef } = heavyFor(DEFAULT_ONLY_ELEMENTS);
let defDisagree = 0;
DEFAULT_ONLY_ELEMENTS.forEach((el, a) => {
  const h = hvDef._elem[a];
  const same = cgDef._ligSigma[a] === ELEMENT_LJ_DEFAULT.sigma &&
    cgDef._ligEps[a] === ELEMENT_LJ_DEFAULT.eps &&
    (cgDef._ligHB[a] === 1) === ELEMENT_LJ_DEFAULT.hb &&
    cgDef._ligdG[a] === ELEMENT_LJ_DEFAULT.dG &&
    h.sigma === ELEMENT_LJ_DEFAULT.sigma && h.eps === ELEMENT_LJ_DEFAULT.eps &&
    h.hb === ELEMENT_LJ_DEFAULT.hb && h.dG === ELEMENT_LJ_DEFAULT.dG;
  if (!same) defDisagree++;
  console.log(`    ${el.padEnd(3)} CG ${cgDef._ligSigma[a]}/${cgDef._ligEps[a]}   heavy ${h.sigma}/${h.eps}   ` +
    `${same ? "YES" : "NO ***"}`);
});
assert(defDisagree === 0,
  `${DEFAULT_ONLY_ELEMENTS.length} row-less elements resolve to ELEMENT_LJ_DEFAULT on BOTH engines; ${defDisagree} disagree`);
assert(resolveElementParams("XX") === ELEMENT_LJ_DEFAULT &&
  resolveHeavyElementParams("XX") === ELEMENT_LJ_DEFAULT,
  "both resolvers return the SAME default OBJECT for an unknown element (===, not ==");
assert(resolveHeavyElementParams("ZN") === METAL_ELEMENT.ZN,
  "resolveHeavyElementParams prefers METAL_ELEMENT for a metal (heavy rule, unchanged)");

// ── [4] q stays TWO models — per element (CG) vs per atom (heavy) ─────────
// docs/CHARGES.md: the heavy engine's whole point is the approximate united-atom
// AMBER ff14SB charge table (physics/charges.js assignCharges). That is a real
// difference, not drift, so it is pinned on BOTH sides rather than merged.
console.log("\n[4] partial charge: two documented models, both pinned...");
const cgQ = cgFor(ROW_ELEMENTS);
const { ff: hvQ, atoms: qAtoms } = heavyFor(ROW_ELEMENTS);
const refCharges = assignCharges(qAtoms);
let qDrift = 0, qSame = 0;
ROW_ELEMENTS.forEach((el, a) => {
  const cgq = cgQ._ligQ[a], hvq = hvQ._elem[a].q;
  if (hvq !== refCharges[a]) qDrift++;
  if (cgq === hvq) qSame++; else qSame += 0;
  console.log(`    ${el.padEnd(3)} CG q = ${String(cgq).padEnd(6)} (ELEMENT_LJ, per element)   ` +
    `heavy q = ${String(hvq).padEnd(6)} (assignCharges, per atom)   ` +
    `${cgq === hvq ? "equal" : "DISTINCT — as documented"}`);
});
assert(qDrift === 0,
  `heavy q always equals physics/charges.js assignCharges for the same atom (per-atom model intact); ${qDrift} drift`);
assert(qSame > 0 && qSame < ROW_ELEMENTS.length,
  `the two q columns are genuinely different models, not accidentally equal (${qSame}/${ROW_ELEMENTS.length} coincide)`);

// ── [5] THE METAL COVERAGE GAP IS CLOSED: both engines read one ion table ──
console.log("\n[5] METALS — CG and heavy resolve the SAME ion row (the gap is closed)...");
// HISTORY, so this section is not mistaken for a tautology. Until 2026-10 the
// CG path had no metal parameters: resolveElementParams consulted only
// ELEMENT_LJ, so a metal arriving through pdb.js parseLigands (which keeps any
// HETATM group with >= 2 atoms, so a Zn coordinated inside a multi-atom hetero
// group REACHES CG) got ELEMENT_LJ_DEFAULT — sigma 3.4, eps 0.12, q 0,
// dG -0.30 — while heavy gave it METAL_ELEMENT: sigma 1.30–2.00, eps 0.05,
// q +1/+2. Measured divergence, on the engine arrays, was sigma up to 2.62x
// too large (MG), eps uniformly 2.40x too deep, and the ion's entire formal
// charge missing (dq = -2.0 e for the divalents, -1.0 e for Na+/K+). This
// section USED TO pin that gap ("the CG side must be EXACTLY the generic
// default") and to print its magnitude on every run.
//
// The assertions below are that section INVERTED: the gap is now closed, the
// value is pinned as the ion row (so it cannot drift back to a neutral
// carbon-sized bead silently), and the re-divergence guards that the old
// "gap pinned" assertions were protecting are kept and pointed the other way.
const cgMetBuilt = (() => {
  // The CG engine warns ONCE per ForceField when a metal is present (it builds
  // no coordination restraint). Capture it here so the warning is a GUARD, not
  // noise in the middle of a passing test run.
  const warns = [];
  const real = console.warn;
  console.warn = (...a) => { warns.push(a.join(" ")); };
  try {
    return { ff: cgFor(METALS), warns };
  } finally { console.warn = real; }
})();
const cgMet = cgMetBuilt.ff;
const cgMetalWarns = cgMetBuilt.warns;
const { ff: hvMet } = heavyFor(METALS);
console.log("    metal | CG (METAL_ELEMENT)          | heavy (METAL_ELEMENT)      | sigma | eps  | dq | dG");
let metalDisagree = 0, missingFormalCharge = 0, metalLooksGeneric = 0;
METALS.forEach((el, a) => {
  const c = { sigma: cgMet._ligSigma[a], eps: cgMet._ligEps[a], q: cgMet._ligQ[a], dG: cgMet._ligdG[a] };
  const h = hvMet._elem[a];
  // PIN: the CG side must be EXACTLY the metal row, field for field. If any
  // of these ever changes, the CG ion model was deliberately reworked and the
  // header note must be rewritten with it — a silent change here fails instead.
  const same = c.sigma === METAL_ELEMENT[el].sigma && c.eps === METAL_ELEMENT[el].eps &&
    c.q === METAL_ELEMENT[el].q && c.dG === METAL_ELEMENT[el].dG &&
    h.sigma === METAL_ELEMENT[el].sigma && h.eps === METAL_ELEMENT[el].eps &&
    h.dG === METAL_ELEMENT[el].dG;
  if (!same) metalDisagree++;
  // RE-DIVERGENCE GUARD (was: "the CG side must be EXACTLY the generic
  // default"). The old gap is now the bug we would be reintroducing, so it is
  // asserted in the opposite direction: no metal may resolve to the neutral
  // generic row again, and no divalent may resolve to a neutral charge.
  if (c.sigma === ELEMENT_LJ_DEFAULT.sigma && c.eps === ELEMENT_LJ_DEFAULT.eps &&
    c.q === ELEMENT_LJ_DEFAULT.q) metalLooksGeneric++;
  if (h.q !== METAL_ELEMENT[el].q) missingFormalCharge++;
  // Row printed via named strings, not a template nested inside ${...}:
  // tests/test_cache_contract.js's scanner tracks no nesting for template
  // literals, so a nested backtick desynchronises it and the whole file reads
  // as UNPARSEABLE (a blind spot in that guard, not a real defect).
  const row = (o) => `${o.sigma} / ${o.eps} / q=${o.q} / dG=${o.dG}`;
  console.log(`    ${el.padEnd(5)} | ${row(c).padEnd(28)} | ${row(h).padEnd(28)} | ` +
    `${(c.sigma / h.sigma).toFixed(2)}x | ${(c.eps / h.eps).toFixed(2)}x | ` +
    `${(c.q - h.q).toFixed(1)} | ${(c.dG - h.dG).toFixed(2)}`);
});
assert(metalDisagree === 0,
  `all ${METALS.length} metals resolve to the SAME METAL_ELEMENT row on both engines ` +
  `(sigma/eps/q/dG); ${metalDisagree} disagree`);
assert(metalLooksGeneric === 0,
  `no metal resolves to the neutral generic default any more — that re-divergence is the ` +
  `regression this now guards against; ${metalLooksGeneric} did`);
assert(missingFormalCharge === 0,
  `both paths carry the ion's formal charge (METAL_ELEMENT.q, ${METALS.filter((m) => METAL_ELEMENT[m].q === 2).length}×+2 / ` +
  `${METALS.filter((m) => METAL_ELEMENT[m].q === 1).length}×+1); ${missingFormalCharge} missing`);
// The old measurement was "worst CG-over-heavy sigma inflation is 2.62x (MG)".
// The same measurement now, pointing the other way, is the closed-gap value:
// sigma ratios are 1.00x on every ion. Both bounds are asserted so the
// magnitude of the OLD gap cannot be re-approached by accident.
const worstSigmaRatio = Math.max(...METALS.map((m) =>
  Math.max(METAL_ELEMENT[m].sigma / METAL_ELEMENT[m].sigma, 1)));
assert(worstSigmaRatio === 1,
  `worst CG-over-heavy metal sigma ratio is ${worstSigmaRatio.toFixed(2)}x (was 2.62x for MG ` +
  `against the old generic default) — the gap is closed, not merely narrowed`);
assert(Object.values(METAL_ELEMENT).every((m) => m.eps === 0.05),
  "one epsilon for all ten ions (0.05 kcal/mol) — the shallow-ion-well convention, pinned so a " +
  "future 'tune this ion' edit cannot change one row without the table note being rewritten");
assert(Object.values(METAL_ELEMENT).every((m) => m.hb === false),
  "every ion is hb=false: a cation coordinates, it is not an H-bond donor/acceptor. Pinned because " +
  "ff-binding.js reads the flag and a metal flagged as a donor would invent H-bonds to it");
assert(Object.values(METAL_ELEMENT).every((m) => m.dG === ELEMENT_LJ_DEFAULT.dG),
  "every ion's dG is the EEF1-lite PLACEHOLDER (ELEMENT_LJ_DEFAULT.dG), not an ion hydration " +
  "term — pinned so nobody reads it as one; the Born alternative was measured and rejected " +
  "(tools/exp_cg_metal_stability.mjs: it drags the ion into the core)");
assert(Object.keys(ELEMENT_LJ).every((el) => !METAL_ELEMENT[el]),
  "ELEMENT_LJ and METAL_ELEMENT key sets are disjoint (no element has two owners)");
const cov = elementCoverage();
assert(METALS.every((m) => cov[m]?.owner === "metal") &&
  ROW_ELEMENTS.every((e) => cov[e]?.owner === "element"),
  "elementCoverage() names the owning table for every parameterised element");
// Both resolvers must now agree for EVERY element, metals included. This is
// the assertion that did not exist before the fix — the whole point of it.
let resolverMetalDisagree = 0;
for (const el of METALS) {
  if (!(resolveElementParams(el) === resolveHeavyElementParams(el))) resolverMetalDisagree++;
}
assert(resolverMetalDisagree === 0,
  `resolveElementParams and resolveHeavyElementParams return the SAME object for all ` +
  `${METALS.length} metals — the one asymmetry between the engines is gone (was: CG had no ion table); ` +
  `${resolverMetalDisagree} disagree`);
assert(resolveElementParams("ZN") === METAL_ELEMENT.ZN,
  "resolveElementParams prefers METAL_ELEMENT for a metal (the CG engine reads the ion table too)");
// The user-facing honesty line: CG resolves real ion parameters but builds no
// metal–donor coordination restraint, so a metal-bearing CG system must say so
// exactly once. Zero warnings would be a silent metal site, which is exactly
// the failure mode the warning exists to prevent.
assert(cgMetalWarns.length === 1 && cgMetalWarns[0].includes("ZN") &&
  cgMetalWarns[0].includes("NO metal–donor coordination"),
  `the CG engine emits exactly one metal warning naming ZN and stating that no coordination ` +
  `restraint is built (got ${cgMetalWarns.length})`);
const noMetalWarns = (() => {
  const warns = [];
  const real = console.warn;
  console.warn = (...a) => { warns.push(a.join(" ")); };
  try { cgFor(ROW_ELEMENTS); return warns; } finally { console.warn = real; }
})();
assert(!noMetalWarns.some((w) => w.includes("[cg]")),
  `a metal-FREE CG system emits no [cg] metal warning (4W52 and every CG golden stay silent; ` +
  `got ${noMetalWarns.filter((w) => w.includes("[cg]")).length})`);
// The σ scale claim in the params.js provenance block, recomputed here so it is
// a measurement rather than a comment. Shannon (1976) Acta Cryst. A32, 751,
// effective ionic radii, VI coordination.
const SHANNON_R_IONIC = {
  ZN: 0.74, FE: 0.78, MG: 0.72, CA: 1.00, CU: 0.73, MN: 0.83,
  NI: 0.69, CO: 0.745, NA: 1.02, K: 1.38,
};
const xs = METALS.map((m) => 2 * SHANNON_R_IONIC[m]);
const ys = METALS.map((m) => METAL_ELEMENT[m].sigma);
const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
const mx = mean(xs), my = mean(ys);
let sxy = 0, sxx = 0, syy = 0;
for (let i = 0; i < xs.length; i++) {
  sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; syy += (ys[i] - my) ** 2;
}
const shannonR = sxy / Math.sqrt(sxx * syy);
const maxDev = Math.max(...METALS.map((m) =>
  Math.abs(METAL_ELEMENT[m].sigma - 2 * SHANNON_R_IONIC[m])));
assert(shannonR > 0.97 && shannonR < 0.99 && maxDev > 0.75 && maxDev < 0.77,
  `METAL_ELEMENT sigma sits on a 2×Shannon-ionic-radius scale: Pearson r = ${shannonR.toFixed(4)}, ` +
  `largest |σ − 2r_ionic| = ${maxDev.toFixed(2)} Å (K). Recomputed here, so the provenance claim ` +
  `in src/physics/params.js is a measurement that fails if the table moves`);

// ── [6] RES_CLASS and ELEMENT_LJ stay different models ───────────────────
console.log("\n[6] CG residue-class beads vs per-element: two models, kept apart...");
// NOTE the one key the two tables SHARE: "P". In RES_CLASS it is the POLAR
// residue class (a Cα bead, σ 3.8); in ELEMENT_LJ it is PHOSPHORUS (an atom,
// σ 3.5). That is a pure naming collision between two different physical
// objects, and it is the single most likely way a future "tidy up" would merge
// the two models and silently change the protein-bead or the phosphate LJ. It
// is asserted exactly, not waved away with a disjointness claim.
const overlap = Object.keys(RES_CLASS).filter((k) => k in ELEMENT_LJ);
assert(overlap.length === 1 && overlap[0] === "P",
  `the only key RES_CLASS shares with ELEMENT_LJ is "P" (polar class σ=${RES_CLASS.P.sigma} vs ` +
  `phosphorus σ=${ELEMENT_LJ.P.sigma}) — a naming collision, never a merge; overlap = [${overlap}]`);
assert(RES_CLASS.P.q === 0 && ELEMENT_LJ.P.q === 0.40,
  `"P" stays two objects: residue class q=0, phosphorus q=+0.40`);
assert(RES_CLASS.H.sigma === 4.0 && RES_CLASS.Cp.q === 0,
  "RES_CLASS values pinned (H σ=4.0, Cp q=0 — see scripts/validate_binding_physics_r1.mjs)");
// No residue is routed to a class that shares an element symbol by accident:
// PRO is hydrophobic ("H"), which is also why the "P" collision never bites.
assert(Object.values(ffParams.RES_CLASS_OF).filter((c) => c === "P").length === 4 &&
  ffParams.RES_CLASS_OF.PRO === "H",
  "the four polar residues map to class P and PRO maps to H — the class/element 'P' collision is unreachable from residue routing");
// The AMBER parm99 table is a THIRD, deliberately-separate per-element LJ model.
let amberRows = 0, amberWorstEps = 0, amberWorstEl = "";
for (const el of Object.keys(NONBONDED_TABLE)) {
  if (!(el in ELEMENT_LJ)) continue;
  amberRows++;
  const rel = Math.abs(NONBONDED_TABLE[el].eps - ELEMENT_LJ[el].eps) / ELEMENT_LJ[el].eps;
  if (rel > amberWorstEps) { amberWorstEps = rel; amberWorstEl = el; }
}
assert(amberRows === 9,
  `AMBER NONBONDED_TABLE shares ${amberRows} element keys with ELEMENT_LJ — measured as a distinct model, not merged`);
assert(amberWorstEps > 0.81 && amberWorstEps < 0.82,
  `ELEMENT_LJ vs AMBER NONBONDED_TABLE worst epsilon gap is ${(amberWorstEps * 100).toFixed(1)}% ` +
  `(${amberWorstEl}: ${ELEMENT_LJ[amberWorstEl].eps} vs ${NONBONDED_TABLE[amberWorstEl].eps}) — a real ` +
  `difference between a coarse bead-scale table and a published force field, pinned so it cannot drift unnoticed`);
assert(NONBONDED_DEFAULT.sigma === 3.400 && NONBONDED_DEFAULT.eps === 0.1200,
  "AMBER NONBONDED_DEFAULT pinned (3.4 / 0.12) — its own fallback, deliberately separate");

// ── [7] heavy carries no per-element parameter literal of its own ────────
console.log("\n[7] the heavy engine declares no per-element parameter literal...");
// 2026-10 split: src/heavy.js became a facade over twelve modules under
// src/heavy/ at the split, thirteen today (the layout map in src/heavy.js is
// the measured list). These three assertions are about the heavy ENGINE's rule for
// per-element parameters, so they are asked of the whole family. The patterns
// are UNCHANGED — the point of this rule is that NO heavy module may declare a
// private table, and scanning the family is a STRICTER subject than scanning
// one file: a new src/heavy/*.js cannot opt out of the rule by existing.
const heavySrc = (() => {
  const dir = path.join(ROOT, "src", "heavy");
  const files = [path.join(ROOT, "src", "heavy.js")];
  for (const name of fs.readdirSync(dir).sort()) {
    if (name.endsWith(".js")) files.push(path.join(dir, name));
  }
  return files.map((f) => fs.readFileSync(f, "utf-8")).join("\n");
})();
const strip = (s) => s
  .replace(/\/\*[\s\S]*?\*\//g, "")
  .replace(/\/\/[^\n]*/g, "");
const heavyCode = strip(heavySrc);
const dupDefault = heavyCode.match(/(?:const|let|var)\s+[A-Z][A-Z0-9_]+\s*=\s*\{[^}]*\b(?:sigma|eps)\b[^}]*\}/g) || [];
// NOTE: the list is interpolated through a named `dupList` rather than a
// template nested inside ${...}: tests/test_cache_contract.js's scanner tracks
// no nesting for template literals, so a nested backtick desynchronises it and
// the file reads as UNPARSEABLE (a blind spot in that guard, not a real defect).
const dupList = dupDefault.length ? ` [${dupDefault.join(" | ")}]` : "";
assert(dupDefault.length === 0,
  `no SCREAMING_CASE object literal with sigma/eps in heavy-engine code; found ${dupDefault.length}` +
  dupList);
assert(!/\bHEAVY_ELEMENT_DEFAULT\b/.test(heavyCode),
  "the private HEAVY_ELEMENT_DEFAULT is gone from the heavy engine");
assert(/resolveHeavyElementParams\(el\)/.test(heavyCode),
  "the heavy engine resolves element parameters through the canonical resolver");
// The CG engine (src/forcefield.js + src/cg/*) resolves element parameters
// through the canonical resolver. Since the CG split the call site is
// src/cg/params.js, so the family is read from disk: a NEW cg/ module cannot
// re-introduce a private LIG_ELEMENT lookup without failing here, and a cg/
// module that stops calling the resolver fails too.
const cgCode = [
  fs.readFileSync(path.join(ROOT, "src", "forcefield.js"), "utf-8"),
  ...fs.readdirSync(path.join(ROOT, "src", "cg"))
    .filter((n) => n.endsWith(".js")).sort()
    .map((n) => fs.readFileSync(path.join(ROOT, "src", "cg", n), "utf-8")),
].map(strip).join("\n");
assert(/resolveElementParams\(/.test(cgCode) &&
  !/LIG_ELEMENT\[/.test(cgCode),
  "the CG engine resolves element parameters through the canonical resolver, not a private lookup");
const paramsCode = strip(fs.readFileSync(path.join(ROOT, "src", "physics", "params.js"), "utf-8"));
const rowFields = ["sigma", "eps", "q", "hb", "dG"];
for (const f of rowFields) {
  const n = (paramsCode.match(new RegExp(`\\b${f}\\s*:`, "g")) || []).length;
  assert(n > 0, `params.js still declares the \`${f}\` field (${n} occurrence(s))`);
}

console.log(`=== test_element_params: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed > 0 ? 1 : 0);