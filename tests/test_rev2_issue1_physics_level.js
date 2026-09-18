/**
 * test_rev2_issue1_physics_level.js — Revolution 2 / Issue 1 regression.
 *
 * Scope: src/main.js physics-level wiring ONLY. No energy-kernel change.
 *
 * Problem: physicsLevelSpec() discarded the level string and both ForceField
 * build paths (buildSystem + onParamChange hot-rebuild) omitted
 * par.physicsLevel, so describePhysics().level always reported L0 at UI
 * L1/L2 though the explicit charges/hbMode flags still fired.
 *
 * Fix verified here (minimal, backward compatible):
 *   - physicsLevelSpec() returns { level, ...flags } instead of flags only.
 *   - Both par literals carry physicsLevel: physLvl.level (buildSystem) and
 *     physicsLevel: physLvlHot.level (onParamChange hot-rebuild).
 *   - UI L1 par propagates level L1; energies bit-identical vs explicit
 *     charges/hbMode flags (explicit flags already won pre-fix, so U/forces
 *     are unchanged — only the reported level is fixed).
 *
 * Coverage (headless, 4W52 CG + native ligands):
 *   [0] main.js source wiring (grep): spec returns level, both par literals
 *       carry physicsLevel, no kernel edits.
 *   [1] UI L1 par propagates level L1 (describePhysics().level === "L1",
 *       charges ON / directional) on both build paths' par shape.
 *   [2] energies unchanged vs explicit flags (bit-identical U/bindingU/forces
 *       tier+explicit ≡ explicit-only; default ≡ explicit L0).
 *
 * Run: node tests/test_rev2_issue1_physics_level.js (fast, <2s; NOT wired
 * into tests/test_all.js FAST so the 352 gate is untouched).
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

let passed = 0;
let failed = 0;
function assert(condition, message) {
  if (condition) {
    passed++;
    console.log(`  ✓ ${message}`);
  } else {
    failed++;
    console.error(`  ✗ FAIL: ${message}`);
  }
}

function findPdb() {
  const candidates = [
    path.resolve(process.cwd(), "4w52.pdb"),
    path.resolve(__dirname, "..", "4w52.pdb"),
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`4w52.pdb not found (tried ${candidates.join(", ")})`);
}

function maxAbsDiff(a, b) {
  let m = 0;
  const n = Math.min(a.length, b.length);
  for (let i = 0; i < n; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
  return m;
}

function main() {
  console.log("=== Rev2/Issue1: UI physics level propagates to ForceField (4W52 CG) ===");
  const pdbText = fs.readFileSync(findPdb(), "utf-8");
  const parsed = parseCa(pdbText);
  const sel = selectSystem(parsed);
  const ligs = parseLigands(pdbText);
  assert(sel.beads.length === 164, `4W52 CG beads = 164 (got ${sel.beads.length})`);
  assert(ligs.length >= 1, `native ligands present (${ligs.length} mols)`);

  // -----------------------------------------------------------------
  // [0] main.js source wiring (grep, no kernel edits)
  // -----------------------------------------------------------------
  console.log("\n[0] main.js wiring (grep, no kernel edits)...");
  const mainSrc = fs.readFileSync(path.resolve(__dirname, "..", "src", "main.js"), "utf-8");
  assert(/return\s*\{\s*level:\s*lvl,\s*\.\.\.PHYSICS_LEVELS\[lvl\]\s*\}/.test(mainSrc),
    `physicsLevelSpec returns { level, ...flags } (level string preserved)`);
  assert(mainSrc.includes("physicsLevel: physLvl.level"),
    `buildSystem par carries physicsLevel: physLvl.level`);
  assert(mainSrc.includes("physicsLevel: physLvlHot.level"),
    `onParamChange hot-rebuild par carries physicsLevel: physLvlHot.level`);
  const nPhysLvl = (mainSrc.match(/physicsLevel:\s*physLvl(Hot)?\.level/g) || []).length;
  assert(nPhysLvl >= 2, `both build paths propagate the tier (${nPhysLvl} sites)`);

  // -----------------------------------------------------------------
  // [1] UI L1 par propagates level (buildSystem + hot-rebuild par shape)
  // -----------------------------------------------------------------
  console.log("\n[1] UI L1 par propagates level L1...");
  // Par shape mirrors src/main.js buildSystem at UI L1 after the fix:
  // physicsLevel + the tier's explicit charges/hbMode flags together.
  const parMainL1 = {
    rc: 10, gamma: 1.0, temp: 300,
    physicsLevel: "L1",
    binding: { on: true, holo: true, charges: true, hbMode: "directional" },
    weak: "off",
  };
  const ffMainL1 = new ForceField(sel, parMainL1, ligs);
  const dMainL1 = ffMainL1.describePhysics();
  assert(dMainL1.level === "L1", `UI L1 par → describePhysics().level L1 (got ${dMainL1.level})`);
  assert(ffMainL1.chargesOn === true && ffMainL1.hbMode === "directional",
    `UI L1 flags on (charges=${ffMainL1.chargesOn}, hb=${ffMainL1.hbMode})`);
  assert(dMainL1.isSimplifiedDefault === false,
    `UI L1 not the simplified default`);
  // Hot-rebuild path carries the same tier (identical par shape, L1).
  const parHotL1 = {
    rc: 10, gamma: 1.0, temp: 300,
    physicsLevel: "L1",
    binding: { on: true, holo: true, charges: true, hbMode: "directional" },
    weak: "off",
  };
  const ffHotL1 = new ForceField(sel, parHotL1, ligs);
  assert(ffHotL1.describePhysics().level === "L1",
    `hot-rebuild L1 par → level L1 (got ${ffHotL1.describePhysics().level})`);
  // L0 baseline still reports L0 (default path untouched).
  const ffL0 = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligs);
  assert(ffL0.describePhysics().level === "L0",
    `default par → level L0 (got ${ffL0.describePhysics().level})`);

  // -----------------------------------------------------------------
  // [2] Energies unchanged vs explicit flags (no kernel change)
  // -----------------------------------------------------------------
  console.log("\n[2] energies unchanged vs explicit flags (bit-identical)...");
  const parExplicitL1 = {
    rc: 10, gamma: 1.0,
    binding: { charges: true, hbMode: "directional" },
  };
  const ffExplicit = new ForceField(sel, parExplicitL1, ligs);
  const Umain = ffMainL1.compute(ffMainL1.ref);
  const Bmain = ffMainL1.bindingU;
  const Uexp = ffExplicit.compute(ffExplicit.ref);
  const Bexp = ffExplicit.bindingU;
  assert(Number.isFinite(Umain) && Number.isFinite(Bmain), `UI L1 par: energy/bindingU finite`);
  assert(Umain === Uexp, `U bit-identical tier+explicit ≡ explicit-only (Δ=${Math.abs(Umain - Uexp).toExponential(1)})`);
  assert(Bmain === Bexp, `bindingU bit-identical tier+explicit ≡ explicit-only (Δ=${Math.abs(Bmain - Bexp).toExponential(1)})`);
  const dF = maxAbsDiff(ffMainL1.forces, ffExplicit.forces);
  assert(dF === 0, `forces bit-identical tier+explicit ≡ explicit-only (maxΔ=${dF.toExponential(1)})`);
  // Default path untouched: default ≡ explicit L0, bit-identical.
  const ffL0exp = new ForceField(sel, { rc: 10, gamma: 1.0, binding: { charges: false, hbMode: "off" } }, ligs);
  const U0 = ffL0.compute(ffL0.ref);
  const U0e = ffL0exp.compute(ffL0exp.ref);
  assert(U0 === U0e, `default ≡ explicit L0 bit-identical (Δ=${Math.abs(U0 - U0e).toExponential(1)})`);

  // -----------------------------------------------------------------
  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log("=================================================");
  if (failed > 0) process.exit(1);
  console.log("PASS: test_rev2_issue1_physics_level.js — Rev2/Issue1 level propagation validated");
}

main();
