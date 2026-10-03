/**
 * test_cg_class_shape.js — prove the forcefield.js split changed NO observable
 * shape of the ForceField class.
 *
 * WHY THIS TEST EXISTS
 * --------------------
 * src/forcefield.js was one 947-line file holding one class. It is now a
 * facade over eight modules under src/cg/ at the split (nine today), and
 * src/cg/forcefield.js declares
 * the class. The obvious risk of that move is the class METHOD becoming
 * something other than a class method: if the methods were attached with
 * Object.assign, every one of them becomes an ENUMERABLE OWN property, and
 * `for (const k in ff)` starts yielding "compute", "rmsd", "_binding" and every
 * other method name. Anything that serialises an instance with for..in, walks
 * it with for..in, or does `{...ff}` would silently change behaviour — and the
 * 1304-assertion suite would not necessarily notice.
 *
 * So this asserts the shape directly, and asserts it against the FACADE
 * (src/forcefield.js) as well as the class module, because the facade must
 * hand back the identical class object — not a copy, not a subclass.
 *
 * WHAT IS PINNED
 *   1. Every class method is a non-enumerable, writable, configurable
 *      prototype property — i.e. a real class method, not an assigned own prop.
 *   2. `for..in` over an instance yields EXACTLY the own enumerable data
 *      fields and nothing from the prototype. This is the assertion that fails
 *      if anyone ever reintroduces prototype composition with Object.assign.
 *   3. The enumerated key SET matches the pre-split list, and the ORDER does
 *      too (for..in order is insertion order for string keys), so a consumer
 *      that JSON-stringifies via for..in gets a byte-identical result.
 *   4. The prototype chain is exactly [ForceField.prototype, Object.prototype].
 *   5. The facade re-exports the SAME class object as cg/forcefield.js.
 *   6. The instance is not sealed/frozen and its fields are writable, i.e. the
 *      move did not accidentally harden the object.
 *
 * Runnable: node tests/test_cg_class_shape.js
 */

import { ForceField as FromFacade, KB_KCAL, KCONV, resolvePhysicsLevel, DEFAULT_PHYSICS_LEVEL, CG_PHYSICS_LEVELS } from "../src/forcefield.js";
import { ForceField as FromClass } from "../src/cg/forcefield.js";
import { parseCa, selectSystem } from "../src/pdb.js";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

/**
 * The complete class-method name set of the PRE-SPLIT src/forcefield.js at git
 * 61374b2. Recorded by reading the original class body, NOT by importing the
 * current one — if the set silently changed, comparing it to itself would
 * agree with the regression.
 */
const EXPECTED_METHODS = [
  "constructor",
  "_pairKey", "_dist", "_angle",
  "compute", "setFunnel", "describePhysics", "rebuildHoloSprings",
  "_harmonicPairs", "_springForces", "setSpringScale", "clearSpringScale",
  "applySeqWeights", "useTirionNetwork", "_angleForces", "_ligandInternal",
  "_ligandBondForces", "_improperForces", "_binding", "_repulsion",
  "_encodeCell", "_cellKey", "_decodeX", "_decodeY", "_decodeZ",
  "kineticTemp", "rmsd",
];

/** Own enumerable field set + order the pre-split constructor produced. */
const EXPECTED_FIELDS = [
  "nProt", "nLigAtoms", "rc", "gamma", "kBond", "kAngle", "epsRep", "sigmaRep",
  "bindOn", "physicsLevel", "hbMode", "chargesOn", "bindRcut", "holoGamma",
  "holoOn", "nHolo", "funnel", "funnelOn", "ref", "isLigand", "ligandAtoms",
  "ligandMasses", "ligandBonds", "ligandAngles", "ligandImpropers", "n",
  "masses", "_maskProt", "resClass", "_protSigma", "_protEps", "_protQ",
  "_protHB", "_vSites", "_ligSigma", "_ligEps", "_ligQ", "_ligHB", "_ligdG",
  "bonds", "angles", "springs", "springK", "springScaleActive", "enmModel",
  "tirionDihedrals", "tirionDihedralK", "tirionSS", "holoSprings", "_excluded",
  "_repSigma", "nativeContacts", "forces", "rcRep", "_cell", "_grid", "_gridB",
  "energy", "bindingU", "desolvU", "trackTerms", "bindLJU", "bindCoulU",
  "bindHBU", "bindU", "_dens", "_dBdn", "_bpA", "_bpJ", "_bpR",
];

function buildFf() {
  const pdbText = fs.readFileSync(path.join(ROOT, "4w52.pdb"), "utf-8");
  const sel = selectSystem(parseCa(pdbText));
  return new FromFacade(sel, { rc: 10, gamma: 1.0 });
}

function main() {
  console.log("=== CG ForceField class shape (for..in / descriptor contract) ===");

  assert(FromFacade === FromClass,
    "the facade re-exports the SAME class object as src/cg/forcefield.js (identity, not a copy)");

  // ---- 1. every method is a real class method -------------------------
  const proto = FromFacade.prototype;
  const ownNames = Object.getOwnPropertyNames(proto).sort();
  const missing = EXPECTED_METHODS.filter((m) => !ownNames.includes(m));
  const extra = ownNames.filter((m) => !EXPECTED_METHODS.includes(m));
  assert(missing.length === 0, `all ${EXPECTED_METHODS.length} pre-split methods are still on the prototype` +
    (missing.length ? ` (missing: ${missing.join(", ")})` : ""));
  assert(extra.length === 0, `the split added no new prototype members` +
    (extra.length ? ` (added: ${extra.join(", ")})` : ""));

  let nonEnum = 0, badDesc = 0;
  for (const m of EXPECTED_METHODS) {
    const d = Object.getOwnPropertyDescriptor(proto, m);
    if (!d || d.enumerable !== false) nonEnum++;
    if (!d || d.writable !== true || d.configurable !== true) badDesc++;
  }
  assert(nonEnum === 0,
    `all ${EXPECTED_METHODS.length} methods are NON-enumerable prototype properties (a real class method, ` +
    `not Object.assign'd own props — the heavy.js failure mode)`);
  assert(badDesc === 0,
    `all ${EXPECTED_METHODS.length} methods are writable + configurable (exact class-method descriptors)`);

  // ---- 2. prototype chain --------------------------------------------
  const chain = [];
  for (let p = proto; p; p = Object.getPrototypeOf(p)) {
    chain.push(p === proto ? "ForceField.prototype" : (p === Object.prototype ? "Object.prototype" : p.constructor.name));
  }
  assert(chain.length === 2 && chain[0] === "ForceField.prototype" && chain[1] === "Object.prototype",
    `prototype chain is exactly [${chain.join(" -> ")}]`);
  assert(FromFacade.prototype.constructor === FromFacade,
    "prototype.constructor points back at the class");

  // ---- 3/4. for..in over an INSTANCE: the behavioural proof -----------
  const ff = buildFf();
  const forIn = [];
  for (const k in ff) forIn.push(k);

  assert(!forIn.some((k) => EXPECTED_METHODS.includes(k)),
    `for..in over an instance yields NO method names (the Object.assign regression this test exists to catch)` +
    ` — leaked: ${forIn.filter((k) => EXPECTED_METHODS.includes(k)).join(", ") || "none"}`);

  const ownKeys = Object.keys(ff);
  assert(JSON.stringify(forIn) === JSON.stringify(ownKeys),
    `for..in enumerates exactly Object.keys(), in the same order (${ownKeys.length} keys)`);

  const missingFields = EXPECTED_FIELDS.filter((k) => !ownKeys.includes(k));
  const newFields = ownKeys.filter((k) => !EXPECTED_FIELDS.includes(k));
  assert(missingFields.length === 0,
    `all ${EXPECTED_FIELDS.length} pre-split own fields are present` +
    (missingFields.length ? ` (missing: ${missingFields.join(", ")})` : ""));
  assert(newFields.length === 0,
    `the split added no new own fields` + (newFields.length ? ` (added: ${newFields.join(", ")})` : ""));

  assert(JSON.stringify(forIn) === JSON.stringify(EXPECTED_FIELDS),
    "for..in ORDER is identical to the pre-split constructor order (so a for..in serializer is byte-identical)");

  // ---- 5. the object is not hardened ----------------------------------
  assert(!Object.isFrozen(ff) && !Object.isSealed(ff) && Object.isExtensible(ff),
    "the instance is still extensible (not frozen/sealed by the move)");
  const before = ff.gamma;
  ff.gamma = 2.5;
  assert(ff.gamma === 2.5, "own fields are still writable");
  ff.gamma = before;

  // ---- 6. the non-class exports still resolve -------------------------
  assert(typeof KB_KCAL === "number" && typeof KCONV === "number",
    `the facade still re-exports KB_KCAL=${KB_KCAL}, KCONV=${KCONV}`);
  assert(DEFAULT_PHYSICS_LEVEL === "L0" && CG_PHYSICS_LEVELS.L0.hbMode === "off",
    "the tier constants still come through the facade");
  assert(resolvePhysicsLevel({}).level === "L0" && resolvePhysicsLevel({ physicsLevel: "L2" }).charges === true,
    "resolvePhysicsLevel still works through the facade");
  assert(typeof ff.compute === "function" && typeof ff.rmsd === "function" &&
         typeof ff.kineticTemp === "function" && typeof ff.describePhysics === "function",
    "the four consumer-facing methods are callable");

  console.log(`\n=== test_cg_class_shape: ${passed} PASSED, ${failed} FAILED ===`);
  process.exit(failed > 0 ? 1 : 0);
}

main();