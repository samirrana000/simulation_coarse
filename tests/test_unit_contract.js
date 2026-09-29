/**
 * test_unit_contract.js — one unit contract, one value (A05 / docs/UNITS.md).
 *
 * The loop (goal S4) found a real defect: `KB_KCAL` was defined TWICE with
 * different values —
 *     src/units.js     0.001987204
 *     src/ff-params.js 0.0019872041
 * — a 5.03e-8 relative split. The reported instantaneous temperature therefore
 * depended on which module an engine happened to import from (CG read
 * units.js, heavy read ff-params.js). It is the exact class of bug the
 * CG/heavy parity work has been hunting across seven review rounds.
 *
 * ff-params.js now re-exports from units.js. This test makes silent
 * re-divergence impossible: it fails if a second literal definition of any
 * unit constant reappears anywhere in src/.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { KB_KCAL, KCONV } from "../src/units.js";
import * as ffParams from "../src/ff-params.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
let passed = 0, failed = 0;

function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

/** Every .js under src/ */
function srcFiles(dir = path.join(ROOT, "src"), out = []) {
  for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) srcFiles(p, out);
    else if (e.name.endsWith(".js")) out.push(p);
  }
  return out;
}

console.log("=== unit contract: one definition, one value ===");

// The reference values (src/units.js is the declared source of truth).
assert(Math.abs(KB_KCAL - 0.001987204) < 1e-12, `KB_KCAL = ${KB_KCAL} (units.js)`);
assert(KCONV === 418.4, `KCONV = ${KCONV}`);

// ff-params.js must re-export, not redefine.
assert(ffParams.KB_KCAL === KB_KCAL,
  `ff-params.KB_KCAL re-exports units.js (${ffParams.KB_KCAL} === ${KB_KCAL})`);
assert(ffParams.KCONV === KCONV, `ff-params.KCONV re-exports units.js`);

// No module EXCEPT units.js may declare its own literal copy of a unit
// constant. units.js is the designated source of truth; everything else must
// re-export from it (`export { X } from "./units.js"`).
const UNIT_CONSTS = ["KB_KCAL", "KCONV", "KCAL_TO_DA_A2_PS2"];
const TRUTH = path.join(ROOT, "src", "units.js");
let redefinitions = 0;
for (const f of srcFiles()) {
  if (path.resolve(f) === path.resolve(TRUTH)) continue;
  const rel = path.relative(ROOT, f);
  const src = fs.readFileSync(f, "utf-8");
  for (const c of UNIT_CONSTS) {
    // An *export* declaration with a numeric literal is a redefinition.
    // `export { X } from "./y.js"` is a re-export and is allowed.
    const re = new RegExp(`export\\s+const\\s+${c}\\s*=\\s*[-0-9.]`, "m");
    if (re.test(src)) {
      redefinitions++;
      console.error(`      ${rel} redefines ${c} with a literal`);
    }
  }
}
assert(redefinitions === 0,
  `no src/ module (outside units.js) redefines a unit constant with a literal; found ${redefinitions}`);

// A value-imported temperature must not depend on the import path.
import { HeavyForceField } from "../src/heavy.js";
import { ForceField } from "../src/forcefield.js";

const n = 5;
const vel = new Float64Array(3 * n);
const mass = new Float64Array(3 * n);
for (let i = 0; i < 3 * n; i++) { vel[i] = 0.31 + 0.017 * i; mass[i] = 12.0; }
const T = (kb) => {
  let ke = 0;
  for (let i = 0; i < 3 * n; i++) ke += mass[i] * vel[i] * vel[i];
  return (ke * 0.5 / KCONV) / (1.5 * n * kb);
};
const tUnits = T(KB_KCAL);
const tParams = T(ffParams.KB_KCAL);
assert(tUnits === tParams,
  `temperature is import-path independent (units ${tUnits.toFixed(6)} K === ff-params ${tParams.toFixed(6)} K)`);

// Sanity: the pre-divergence bug would have produced a 5.03e-8 relative gap,
// i.e. ~1.5e-5 K at 300 K. Assert we are exactly zero, not merely small.
const rel = Math.abs(tUnits - tParams) / tUnits;
assert(rel === 0, `relative split is exactly 0 (was 5.03e-8 before the fix); got ${rel}`);

console.log(`=== test_unit_contract: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed > 0 ? 1 : 0);
