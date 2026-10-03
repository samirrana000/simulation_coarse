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
 *   • and there IS a live CG-vs-heavy divergence, but it is a COVERAGE gap
 *     rather than a duplicated literal: the CG path has no metal parameters at
 *     all, so a metal arriving through pdb.js parseLigands is handed
 *     ELEMENT_LJ_DEFAULT — σ 3.4 Å, ε 0.12 kcal/mol, q 0 — while the heavy
 *     path hands it METAL_ELEMENT — σ 1.30–2.00 Å, ε 0.05 kcal/mol, q +1/+2.
 *     §5 pins and prints that gap. It is NOT fixed here: choosing a CG ion
 *     model is a physics decision, not a refactor.
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
 *   5. [5] The metal gap, measured and printed on every run, with the current
 *      values pinned so it cannot widen unnoticed.
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

// ── [5] THE LIVE DIVERGENCE: the CG path has no metal parameters ─────────
console.log("\n[5] LIVE CG-vs-heavy DIVERGENCE — CG has no metal table (measured, not fixed)...");
const cgMet = cgFor(METALS);
const { ff: hvMet } = heavyFor(METALS);
console.log("    metal | CG (fallback to ELEMENT_LJ_DEFAULT) | heavy (METAL_ELEMENT)      | sigma | eps   | dq(e)");
let metalBad = 0, missingFormalCharge = 0;
METALS.forEach((el, a) => {
  const c = { sigma: cgMet._ligSigma[a], eps: cgMet._ligEps[a], q: cgMet._ligQ[a], dG: cgMet._ligdG[a] };
  const h = hvMet._elem[a];
  // PIN: the CG side must be EXACTLY the generic default. If this ever changes,
  // the CG metal model was deliberately reworked and the header note must be
  // rewritten with it — a silent change here fails instead.
  if (!(c.sigma === ELEMENT_LJ_DEFAULT.sigma && c.eps === ELEMENT_LJ_DEFAULT.eps &&
    c.q === ELEMENT_LJ_DEFAULT.q && c.dG === ELEMENT_LJ_DEFAULT.dG)) metalBad++;
  if (h.q !== METAL_ELEMENT[el].q) missingFormalCharge++;
  console.log(`    ${el.padEnd(5)} | ${`${c.sigma} / ${c.eps} / q=${c.q} / dG=${c.dG}`.padEnd(35)} | ` +
    `${`${h.sigma} / ${h.eps} / q=${h.q}`.padEnd(24)} | ${(c.sigma / h.sigma).toFixed(2)}x | ` +
    `${(c.eps / h.eps).toFixed(2)}x | ${(c.q - h.q).toFixed(1)}`);
});
assert(metalBad === 0,
  `all ${METALS.length} metals still resolve to the generic default in CG (gap pinned, not silently narrowed); ${metalBad} changed`);
assert(missingFormalCharge === 0,
  `the heavy path does carry the ion's formal charge (METAL_ELEMENT.q, ${METALS.filter((m) => METAL_ELEMENT[m].q === 2).length}×+2 / ` +
  `${METALS.filter((m) => METAL_ELEMENT[m].q === 1).length}×+1); ${missingFormalCharge} missing`);
const worstSigma = Math.max(...METALS.map((m) =>
  (ELEMENT_LJ_DEFAULT.sigma / METAL_ELEMENT[m].sigma)));
assert(worstSigma > 2.6 && worstSigma < 2.63,
  `worst CG-over-heavy metal sigma inflation is ${worstSigma.toFixed(2)}x (MG) — the documented magnitude`);
assert(Object.keys(ELEMENT_LJ).every((el) => !METAL_ELEMENT[el]),
  "ELEMENT_LJ and METAL_ELEMENT key sets are disjoint (no element has two owners)");
const cov = elementCoverage();
assert(METALS.every((m) => cov[m]?.owner === "metal") &&
  ROW_ELEMENTS.every((e) => cov[e]?.owner === "element"),
  "elementCoverage() names the owning table for every parameterised element");

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
// 2026-10 split: src/heavy.js became a facade over eleven modules under
// src/heavy/. These three assertions are about the heavy ENGINE's rule for
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
assert(dupDefault.length === 0,
  `no SCREAMING_CASE object literal with sigma/eps in heavy-engine code; found ${dupDefault.length}` +
  `${dupDefault.length ? ` (${dupDefault.join(" | ")})` : ""}`);
assert(!/\bHEAVY_ELEMENT_DEFAULT\b/.test(heavyCode),
  "the private HEAVY_ELEMENT_DEFAULT is gone from the heavy engine");
assert(/resolveHeavyElementParams\(el\)/.test(heavyCode),
  "the heavy engine resolves element parameters through the canonical resolver");
const ffCode = strip(fs.readFileSync(path.join(ROOT, "src", "forcefield.js"), "utf-8"));
assert(/resolveElementParams\(/.test(ffCode) &&
  !/LIG_ELEMENT\[/.test(ffCode),
  "forcefield.js resolves element parameters through the canonical resolver, not a private lookup");
const paramsCode = strip(fs.readFileSync(path.join(ROOT, "src", "physics", "params.js"), "utf-8"));
const rowFields = ["sigma", "eps", "q", "hb", "dG"];
for (const f of rowFields) {
  const n = (paramsCode.match(new RegExp(`\\b${f}\\s*:`, "g")) || []).length;
  assert(n > 0, `params.js still declares the \`${f}\` field (${n} occurrence(s))`);
}

console.log(`=== test_element_params: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed > 0 ? 1 : 0);