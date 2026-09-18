/**
 * test_rev3_issue1_heavy_physics.js — Revolution 3 / Issue 1 regression.
 *
 * Scope: src/heavy.js + src/main.js comment ONLY. No energy-kernel change.
 *
 * Problem: HeavyForceField constructor ignored par.physicsLevel /
 * par.binding.charges / par.binding.hbMode, so UI L1/L2 silently no-opped in
 * heavy except the weak flag.
 *
 * Fix verified here (minimal, backward compatible):
 *   - HeavyForceField stores physicsLevel + mirrors chargesOn/hbMode
 *     queryably via describePhysics() (tier defaults + explicit-binding-wins,
 *     unknown → L0 — same rule as ForceField/resolvePhysicsLevel).
 *   - Energy kernels untouched: U/forces bit-identical across tiers.
 *   - Misleading main.js comment updated (heavy mirrors the tier queryably).
 *
 * Coverage (headless, 4W52 heavy):
 *   [0] source wiring (grep): heavy resolver + stored fields + describePhysics,
 *       main.js comment no longer claims heavy consumes only par.weak.
 *   [1] UI L1 par → physicsLevel L1, chargesOn true, hbMode directional,
 *       describePhysics().level L1 (required asserts).
 *   [2] default → L0/off + explicit-flags-win + unknown→L0.
 *   [3] kernels unchanged: U/bindingU/forces bit-identical tier vs default.
 *
 * Run: node tests/test_rev3_issue1_heavy_physics.js (fast, <10s; NOT wired
 * into tests/test_all.js FAST so the 352 gate is untouched).
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { parseHeavy, selectHeavy, HeavyForceField } from "../src/heavy.js";

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
  console.log("=== Rev3/Issue1: HeavyForceField mirrors physicsLevel (4W52 heavy) ===");
  const pdbText = fs.readFileSync(findPdb(), "utf-8");
  const parsed = parseHeavy(pdbText);
  const sel = selectHeavy(parsed);
  assert(sel.atoms.length > 1000, `4W52 heavy atoms present (${sel.atoms.length})`);

  // -----------------------------------------------------------------
  // [0] source wiring (grep, no kernel edits)
  // -----------------------------------------------------------------
  console.log("\n[0] source wiring (grep, no kernel edits)...");
  const heavySrc = fs.readFileSync(path.resolve(__dirname, "..", "src", "heavy.js"), "utf-8");
  assert(/resolveHeavyPhysicsLevel/.test(heavySrc),
    `heavy.js defines resolveHeavyPhysicsLevel`);
  assert(/this\.physicsLevel\s*=\s*_phys\.level/.test(heavySrc),
    `constructor stores this.physicsLevel`);
  assert(/this\.chargesOn\s*=\s*_phys\.charges/.test(heavySrc),
    `constructor mirrors this.chargesOn`);
  assert(/this\.hbMode\s*=\s*_phys\.hbMode/.test(heavySrc),
    `constructor mirrors this.hbMode`);
  assert(/describePhysics\(\)/.test(heavySrc),
    `HeavyForceField exposes describePhysics()`);
  const mainSrc = fs.readFileSync(path.resolve(__dirname, "..", "src", "main.js"), "utf-8");
  assert(mainSrc.includes("mirrors physicsLevel"),
    `main.js comment updated (heavy mirrors the tier queryably)`);

  // -----------------------------------------------------------------
  // [1] UI L1 par propagates level (required asserts)
  // -----------------------------------------------------------------
  console.log("\n[1] UI L1 par → queryable mirror...");
  // Par shape mirrors src/main.js buildSystem at UI L1 after the fix.
  const parMainL1 = {
    gamma: 2.0, temp: 300,
    physicsLevel: "L1",
    binding: { on: true, holo: true, charges: true, hbMode: "directional" },
    weak: "off",
  };
  const ffL1 = new HeavyForceField({ atoms: sel.atoms }, parMainL1, []);
  assert(ffL1.physicsLevel === "L1",
    `HeavyForceField physicsLevel L1 (got ${ffL1.physicsLevel})`);
  assert(ffL1.chargesOn === true,
    `HeavyForceField chargesOn true (got ${ffL1.chargesOn})`);
  assert(ffL1.hbMode === "directional",
    `HeavyForceField hbMode directional (got ${ffL1.hbMode})`);
  assert(ffL1.describePhysics().level === "L1",
    `describePhysics level L1 (got ${ffL1.describePhysics().level})`);
  assert(ffL1.describePhysics().isSimplifiedDefault === false,
    `L1 not the simplified default`);
  // Tier-only par (no explicit binding flags) resolves identically.
  const ffTierL1 = new HeavyForceField({ atoms: sel.atoms }, { gamma: 2.0, physicsLevel: "L1" }, []);
  assert(ffTierL1.physicsLevel === "L1" && ffTierL1.chargesOn === true && ffTierL1.hbMode === "directional",
    `tier-only L1 par resolves (level=${ffTierL1.physicsLevel}, charges=${ffTierL1.chargesOn}, hb=${ffTierL1.hbMode})`);

  // -----------------------------------------------------------------
  // [2] default L0 + explicit-win + unknown fallback
  // -----------------------------------------------------------------
  console.log("\n[2] defaults + precedence...");
  const ffDef = new HeavyForceField({ atoms: sel.atoms }, { gamma: 2.0 }, []);
  assert(ffDef.physicsLevel === "L0" && ffDef.chargesOn === false && ffDef.hbMode === "off",
    `default → L0/off (got ${ffDef.physicsLevel}/${ffDef.hbMode})`);
  assert(ffDef.describePhysics().level === "L0" && ffDef.describePhysics().isSimplifiedDefault === true,
    `default describePhysics L0 simplified default`);
  const ffWin = new HeavyForceField(
    { atoms: sel.atoms },
    { gamma: 2.0, physicsLevel: "L1", binding: { charges: false, hbMode: "off" } }, []);
  assert(ffWin.chargesOn === false && ffWin.hbMode === "off",
    `explicit binding flags win over the tier (L1 + explicit off → off)`);
  const ffBad = new HeavyForceField({ atoms: sel.atoms }, { gamma: 2.0, physicsLevel: "L9" }, []);
  assert(ffBad.physicsLevel === "L0",
    `unknown level falls back to L0 (got ${ffBad.physicsLevel})`);

  // -----------------------------------------------------------------
  // [3] kernels unchanged (bit-identical U/forces across tiers)
  // -----------------------------------------------------------------
  console.log("\n[3] kernels unchanged (bit-identical energy)...");
  const U1 = ffL1.compute(ffL1.ref);
  const B1 = ffL1.bindingU;
  const U0 = ffDef.compute(ffDef.ref);
  const B0 = ffDef.bindingU;
  assert(Number.isFinite(U1) && Number.isFinite(B1), `L1 energy/bindingU finite`);
  assert(Number.isFinite(U0) && Number.isFinite(B0), `L0 energy/bindingU finite`);
  assert(U1 === U0, `U bit-identical L1 ≡ L0 default (Δ=${Math.abs(U1 - U0).toExponential(1)})`);
  assert(B1 === B0, `bindingU bit-identical L1 ≡ L0 default (Δ=${Math.abs(B1 - B0).toExponential(1)})`);
  const dF = maxAbsDiff(ffL1.forces, ffDef.forces);
  assert(dF === 0, `forces bit-identical L1 ≡ L0 default (maxΔ=${dF.toExponential(1)})`);

  // -----------------------------------------------------------------
  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log("=================================================");
  if (failed > 0) process.exit(1);
  console.log("PASS: test_rev3_issue1_heavy_physics.js — Rev3/Issue1 heavy mirror validated");
}

main();
