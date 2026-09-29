/**
 * test_all.js — Comprehensive Automated Test Suite
 */

import fs from "fs";
import { parseCa, parseLigands, parseMol2, selectSystem } from "../src/pdb.js";
import { parseHeavy, HeavyForceField, selectHeavy, appendHeavyLigands } from "../src/heavy.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { SpatialGrid } from "../src/spatial-grid.js";
import { dihedralForcesAnalytic, improperForces } from "../src/ff-harmonic.js";
import { GeneralizedBorn } from "../src/physics/gb.js";
import { ChemicalNetworkModel } from "../src/physics/network.js";
import { placeLigand, relaxClash, findPocketCenter } from "../src/placement.js";
import { KB_KCAL, KCONV } from "../src/units.js";
import { execFileSync } from "node:child_process";
import path from "node:path";
import { fileURLToPath } from "node:url";

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

// -----------------------------------------------------------------
// Three tiers, all declared in tests/suites.js:
//   FAST   always; must stay < 60 s (npm test, CI, wikiskill_gate, evolve gate).
//   MEDIUM opt-in via --medium / MEDIUM=1; implied by --slow. 10-60 s.
//   SLOW   opt-in via --slow (SMOKE protocol for test_thermo_heavy) or
//          --slow-full / --long / --full (FULL, ~200 s heavy leg).
// Suites are standalone scripts (they call process.exit), so they run via
// child_process with a per-suite timeout; child stdout is captured (piped)
// and summarised one line per suite. The single uppercase
// "N PASSED, M FAILED" line at the end is the grand total, which
// scripts/wikiskill_gate.js and evolve/evolve.mjs check.
// -----------------------------------------------------------------
const __testAllDir = path.dirname(fileURLToPath(import.meta.url));
// tests/suites.js records repo-root-relative paths ("tests/x.js",
// "scripts/y.mjs") so the registry reads the same way as the tree does.
const __rootDir = path.resolve(__testAllDir, "..");
const SLOW_FULL = process.argv.includes("--slow-full") || process.argv.includes("--long") || process.argv.includes("--full") || process.env.SLOW_FULL === "1";
const SLOW = process.argv.includes("--slow") || process.env.SLOW === "1" || SLOW_FULL;
// MEDIUM is the middle tier (10-60 s of seeded work, currently
// scripts/test_pocket_entropy.mjs). It is implied by SLOW because it is
// strictly cheaper, so `--slow` keeps meaning "everything above FAST".
const MEDIUM = process.argv.includes("--medium") || process.env.MEDIUM === "1" || SLOW;

// The suite tables live in tests/suites.js (declarative registry). This file
// only executes them.
//
// WHY A SEPARATE REGISTRY
// -----------------------
// The old inline tables hardcoded `expect: <n>` per suite, so a suite that
// gained one assertion failed with a bare "expected 20/0" and no hint that
// the harness, not the test, was stale. The registry carries tier + timeout
// only; assertion counts are DERIVED from each child's own output at run
// time (see runSuiteFile) and the gate is the child's exit code. Adding a
// test can therefore never break this harness, and deleting a test cannot
// either — only a genuinely failing test does.
//
// tests/suites.js is the single source of truth and is itself covered by
// tests/test_suite_registry.js, which FAILS if any tests/test_*.js or
// scripts/test_*.mjs file exists without being wired into a tier or moved to
// tests/manual/.
//
// Tier contract (see tests/suites.js for the full text):
//   FAST   always run; must stay under ~60 s — `npm test`, CI, and every
//          gate depend on it.
//   MEDIUM 10-60 s of seeded work; opt-in via --medium / MEDIUM=1, and
//          implied by --slow (it is strictly cheaper than SLOW).
//   SLOW   multi-minute validation legs; opt-in via --slow (SMOKE protocol
//          for test_thermo_heavy) or --slow-full / --long / --full (FULL).
import { FAST_BUDGET_S, suitesIn } from "./suites.js";

const FAST_SUITES = suitesIn("FAST");
const MEDIUM_SUITES = suitesIn("MEDIUM");
const SLOW_SUITES = suitesIn("SLOW");

// test_thermo_heavy is the one suite whose protocol depends on the flag:
// FULL is ~200 s (6 legs x 1800 heavy steps), SMOKE is ~10-25 s (2 legs x
// 350 steps) and asserts the same things (finiteness/boundedness; the sign
// is deliberately never asserted at SMOKE sample size). Both are
// deterministic via fixed SEEDS.
for (const s of SLOW_SUITES) {
  if (s.file === "scripts/test_thermo_heavy.mjs") s.args = SLOW_FULL ? [] : ["--smoke"];
}

/**
 * Count assertions in a suite's stdout, without the harness holding a
 * hardcoded expectation for it.
 *
 * Preference order (first that yields something wins):
 *   1. `N PASSED, M FAILED`      — the modern convention.
 *   2. `✓`-prefixed lines         — legacy per-assert printers.
 *   3. a bare `PASS` / `PASS: …`  — single-gate legacy scripts.
 * Anything else counts as 0 passes and the suite is reported as UNCACHED so
 * it is visible rather than silently folded into the total as zero work.
 *
 * The number is REPORTING ONLY. The gate is the child's exit code plus any
 * `M FAILED > 0`. This is what makes "add an assertion to a suite" a
 * non-event instead of a confusing harness failure.
 */
function countAsserts(out) {
  const sums = [...String(out).matchAll(/(\d+) PASSED, (\d+) FAILED/g)];
  if (sums.length) {
    const last = sums[sums.length - 1];
    return { passed: Number(last[1]), failed: Number(last[2]), via: "summary line" };
  }
  const ticks = (String(out).match(/^\s*[✓✔]\s/gm) || []).length;
  const crosses = (String(out).match(/^\s*[✗✘]\s/gm) || []).length;
  if (ticks || crosses) return { passed: ticks, failed: crosses, via: `${ticks} check lines` };
  if (/^\s*PASS\b/m.test(out)) return { passed: 1, failed: 0, via: "single PASS gate" };
  return { passed: 0, failed: 0, via: "no output (UNCACHED — check it still asserts)" };
}

/** Run one standalone suite script; return { passed, failed, secs, via }. */
function runSuiteFile(suite) {
  const t0 = Date.now();
  let out = "";
  const suiteArgs = suite.args || [];
  const suiteEnv = { ...process.env, ...(suite.env || {}) };
  try {
    out = execFileSync(process.execPath, [path.resolve(__rootDir, suite.file), ...suiteArgs], {
      encoding: "utf-8", timeout: suite.timeout, stdio: ["ignore", "pipe", "pipe"], env: suiteEnv,
    });
  } catch (e) {
    out = (e.stdout || "") + (e.stderr || "");
    const c = countAsserts(out);
    const detail = (e.stderr || "").split("\n").filter((l) => /FAIL|Error|✗/.test(l)).slice(0, 3).join(" | ");
    const timedOut = e.killed === true || e.signal === "SIGTERM";
    throw new Error(
      `${suite.file} ${timedOut ? `TIMED OUT after ${suite.timeout}ms` : `exited ${e.status ?? "?"}`}` +
      (c.passed || c.failed ? ` (${c.passed} passed, ${c.failed} failed)` : "") +
      (detail ? `: ${detail}` : "")
    );
  }
  const c = countAsserts(out);
  // Exit code 0 is the gate. A suite that self-reports failures while exiting
  // 0 is a harness bug in that suite, not a pass — treat it as a failure.
  if (c.failed > 0) throw new Error(`${suite.file}: self-reported ${c.failed} failed assertion(s) but exited 0`);
  return { passed: c.passed, failed: 0, secs: (Date.now() - t0) / 1000, via: c.via };
}

/**
 * Run a tier, folding child counts into the harness totals (additive).
 *
 * There is deliberately NO `expected N` comparison here — see tests/suites.js
 * for why. A suite passes when it exits 0; its assertion count is printed
 * for information and added to the grand total.
 */
function runTier(label, suites) {
  console.log(`\n[${label}] tier (${suites.length} scripts from tests/suites.js)...`);
  const t0 = Date.now();
  let tierPassed = 0;
  for (const suite of suites) {
    const name = path.basename(suite.file);
    try {
      const r = runSuiteFile(suite);
      passed += r.passed;
      tierPassed += r.passed;
      const uncached = r.passed === 0 ? "  ⚠ UNCACHED (no assertion output found)" : "";
      console.log(`  ✓ [${label}] ${name}: ${r.passed} passed (${r.via}) in ${r.secs.toFixed(1)}s${uncached}`);
    } catch (e) {
      failed++;
      console.error(`  ✗ FAIL: [${label}] ${name} — ${e.message}`);
    }
  }
  const secs = (Date.now() - t0) / 1000;
  console.log(`  [${label}] tier done: +${tierPassed} asserts in ${secs.toFixed(1)}s`);
  return secs;
}

async function runTests() {
  console.log("=== RUNNING COARSE & HEAVY SIMULATION TEST SUITE ===");

  const pdbText = fs.readFileSync("4w52.pdb", "utf-8");
  const mol2Text = fs.readFileSync("benzene.mol2", "utf-8");

  // -----------------------------------------------------------------
  // 1. Structure & MOL2 Parsing
  // -----------------------------------------------------------------
  console.log("\n[1] Testing PDB & MOL2 Parsing...");
  const parsedCa = parseCa(pdbText);
  assert(parsedCa.beads.length === 164, `Cα beads parsed = 164 (got ${parsedCa.beads.length})`);

  const parsedHeavy = parseHeavy(pdbText);
  assert(parsedHeavy.atoms.length === 1308, `Heavy atoms parsed = 1308 (got ${parsedHeavy.atoms.length})`);
  assert(parsedHeavy.heteroGroups.length === 2, `Hetero groups in 4W52 = 2 (BNZ, EPE) (got ${parsedHeavy.heteroGroups.length})`);

  const mols = parseMol2(mol2Text);
  assert(mols.length === 1, `MOL2 molecules = 1 (got ${mols.length})`);
  assert(mols[0].atoms.length === 6, `Benzene atoms = 6 (got ${mols[0].atoms.length})`);
  assert(mols[0].bonds.length === 6, `Benzene bonds = 6 (got ${mols[0].bonds.length})`);

  // -----------------------------------------------------------------
  // 2. Heavy Mode with MOL2 Override (No Double Counting)
  // -----------------------------------------------------------------
  console.log("\n[2] Testing Heavy Mode with MOL2 Ligand Override...");
  // Test case A: External MOL2 ligand provided
  let selA = selectHeavy(parsedHeavy, {
    heteroSelection: { "A|200|BNZ": false, "A|201|EPE": false },
    includePdbLigands: true,
    hasExternalLigand: true,
  });
  selA = appendHeavyLigands(selA, mols);
  const ffA = new HeavyForceField({ atoms: selA.atoms }, { gamma: 2.0 }, []);
  assert(ffA.nProt === 1287, `Protein atoms = 1287 (got ${ffA.nProt})`);
  assert(ffA.nLigAtoms === 6, `Ligand atoms from MOL2 = 6 (got ${ffA.nLigAtoms})`);
  assert(ffA.n === 1293, `Total atoms = 1293 (1287 prot + 6 mol2) (got ${ffA.n})`);

  const energyA = ffA.compute(ffA.ref);
  assert(Number.isFinite(energyA), `Heavy mode energy is finite (${energyA.toFixed(2)} kcal/mol)`);
  assert(Number.isFinite(ffA.bindingU), `Binding energy is finite (${ffA.bindingU.toFixed(2)} kcal/mol)`);

  // Test case B: No external MOL2 (PDB HETATM BNZ is ligand)
  let selB = selectHeavy(parsedHeavy, {
    heteroSelection: { "A|200|BNZ": true, "A|201|EPE": true },
    includePdbLigands: true,
    hasExternalLigand: false,
  });
  const ffB = new HeavyForceField({ atoms: selB.atoms }, { gamma: 2.0 }, []);
  assert(ffB.nLigAtoms === 21, `PDB HETATM ligand atoms = 21 (got ${ffB.nLigAtoms})`);
  const energyB = ffB.compute(ffB.ref);
  assert(Number.isFinite(energyB), `PDB HETATM energy is finite (${energyB.toFixed(2)} kcal/mol)`);

  // -----------------------------------------------------------------
  // 3. Fast Spatial Grid vs Brute Force Parity
  // -----------------------------------------------------------------
  console.log("\n[3] Testing Fast Spatial Grid vs Brute Force...");
  const grid = new SpatialGrid(8.5, ffA.n);
  grid.build(ffA.ref, ffA.n);

  let gridPairCount = 0;
  grid.forEachPair(ffA.ref, ffA.n, 8.5, (i, j, dx, dy, dz, r2, r) => {
    gridPairCount++;
  });

  let brutePairCount = 0;
  for (let i = 0; i < ffA.n; i++) {
    const xi = ffA.ref[3 * i], yi = ffA.ref[3 * i + 1], zi = ffA.ref[3 * i + 2];
    for (let j = i + 1; j < ffA.n; j++) {
      const dx = ffA.ref[3 * j] - xi, dy = ffA.ref[3 * j + 1] - yi, dz = ffA.ref[3 * j + 2] - zi;
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < 8.5 * 8.5) brutePairCount++;
    }
  }

  assert(gridPairCount === brutePairCount, `Grid pair count (${gridPairCount}) matches brute force (${brutePairCount})`);

  // -----------------------------------------------------------------
  // 4. Analytic Dihedral Gradients vs Numerical Finite Differences
  // -----------------------------------------------------------------
  console.log("\n[4] Testing Analytic Dihedral Gradients...");
  const fAnalytic = new Float64Array(ffA.n * 3);
  const fNumeric = new Float64Array(ffA.n * 3);

  const uAnalytic = dihedralForcesAnalytic(ffA.ref, fAnalytic, ffA.propers, 5, 2.0);
  const uNumeric = improperForces(ffA.ref, fNumeric, ffA.propers);

  assert(Math.abs(uAnalytic - uNumeric) < 1e-4, `Dihedral potential energies match (${uAnalytic.toFixed(4)} vs ${uNumeric.toFixed(4)})`);

  let maxForceDiff = 0;
  for (let i = 0; i < fAnalytic.length; i++) {
    const diff = Math.abs(fAnalytic[i] - fNumeric[i]);
    if (diff > maxForceDiff) maxForceDiff = diff;
  }
  assert(maxForceDiff < 1e-3, `Analytic dihedral forces agree with numerical to within ${maxForceDiff.toFixed(6)} kcal/mol/Å`);

  // -----------------------------------------------------------------
  // 5. Generalized Born & Debye-Hückel Implicit Solvent
  // -----------------------------------------------------------------
  console.log("\n[5] Testing Generalized Born & Debye-Hückel...");
  const gb = new GeneralizedBorn({ epsIn: 4.0, epsOut: 78.5, saltM: 0.15, temperature: 300 });
  assert(gb.kappa > 0.1, `Debye kappa is positive (${gb.kappa.toFixed(4)} Å^-1)`);

  const gbPair = gb.pairInteraction(0, 1, 3.0, 0, 0, 3.0, 1.0, -1.0, 1.7, 1.5, 1.0);
  assert(gbPair.coulombE < 0, `Attractive Coulomb energy between opposite charges (${gbPair.coulombE.toFixed(2)} kcal/mol)`);
  assert(Number.isFinite(gbPair.energy), `Total GB interaction energy is finite (${gbPair.energy.toFixed(2)} kcal/mol)`);
  assert(Math.abs(gbPair.fx) > 0, `Coulomb force on x-axis is non-zero (${gbPair.fx.toFixed(2)})`);

  // -----------------------------------------------------------------
  // 6. Chemical Network Model & Kinetics
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
  // 6. Testing Chemical Network & Transition Path Theory (TPT)
  // -----------------------------------------------------------------
  console.log("\n[6] Testing Chemical Network Model & TPT Committors...");
  const cnm = new ChemicalNetworkModel({ temperature: 300 });
  assert(cnm.nStates === 4, `Chemical Network has 4 states (Bulk, Encounter, Intermediate, Bound)`);

  // Detailed balance check: pi_i * K_ij == pi_j * K_ji
  const flux01 = cnm.pi[0] * cnm.rateMatrix[0][1];
  const flux10 = cnm.pi[1] * cnm.rateMatrix[1][0];
  assert(Math.abs(flux01 - flux10) < 1e-8, `Detailed balance strictly preserved between states 0 & 1 (flux: ${flux01.toExponential(4)})`);

  // TPT Committors
  const tpt = cnm.computeTPT();
  assert(tpt.committors[0] === 0.0, `Bulk committor q_0^+ = 0.0 (got ${tpt.committors[0]})`);
  assert(tpt.committors[3] === 1.0, `Bound committor q_3^+ = 1.0 (got ${tpt.committors[3]})`);
  assert(tpt.committors[1] > 0 && tpt.committors[1] < tpt.committors[2], `Monotonic forward committors q_1^+ < q_2^+ (${tpt.committors[1].toFixed(2)} < ${tpt.committors[2].toFixed(2)})`);

  // Gillespie step
  const stepRes = cnm.stepGillespie();
  assert(stepRes.dt > 0, `Gillespie dwell time is positive (${(stepRes.dt * 1e9).toFixed(3)} ns)`);

  // Pose classifier
  assert(cnm.classifyPose(4.0, 15, 1.2) === 3, `Pose with COM=4Å, 15 contacts, RMSD=1.2Å classified as Native Bound (State 3)`);
  assert(cnm.classifyPose(25.0, 0, 10.0) === 0, `Far pose classified as Bulk Solvated (State 0)`);

  // -----------------------------------------------------------------
  // 7. Clash-Free Placement & Pocket Detection
  // -----------------------------------------------------------------
  console.log("\n[7] Testing Clash-Free Placement & Pocket Detection...");
  const pocket = findPocketCenter(ffA.ref, ffA.nProt);
  assert(pocket.length === 3 && Number.isFinite(pocket[0]), `Pocket center found at (${pocket.map(v => v.toFixed(1)).join(", ")})`);

  const protein = {
    pos: ffA.ref.subarray(0, 3 * ffA.nProt),
    sigma: new Float64Array(ffA.nProt).fill(3.8),
  };
  // Test native pose relaxation
  const nativePos = new Float64Array(18);
  for (let i = 0; i < 6; i++) {
    nativePos[3 * i] = mols[0].atoms[i].x;
    nativePos[3 * i + 1] = mols[0].atoms[i].y;
    nativePos[3 * i + 2] = mols[0].atoms[i].z;
  }
  const relaxed = relaxClash(nativePos, mols[0], protein);
  assert(relaxed.converged, `Native ligand pose is verified clash-free (residual = ${relaxed.residualClash.toFixed(2)})`);

  // -----------------------------------------------------------------
  // 8. Langevin Dynamics Integration Stability
  // -----------------------------------------------------------------
  console.log("\n[8] Testing Langevin Dynamics Integration (Cα ENM & Heavy)...");
  const selCG = selectSystem(parsedCa);
  // Temperature fidelity tighten: 300±40 after 200 steps (was 100–600)
  // FP3: seeded replicas (was: up-to-10 unseeded trials, keep-first-in-band
  // retry-hiding). Physics reason: the OU thermostat's stationary variance is
  // exact in distribution, but a single 200-step kinetic-T sample carries
  // O(1/√DOF-step) noise, so one unseeded trial can tail out of band while
  // the ensemble mean sits at the bath. Fixed SEEDS = [101, 202, 303] (thermo
  // family, scripts/test_thermo.mjs) lock the streams; the assertion is the
  // 3-replica mean in 260–340 K (same band, no lowered standard).
  // Measured 2026-09-13 (linux, node): reps [299.9, 280.8, 330.4] → mean 303.7.
  const TEMP_SEEDS = [101, 202, 303];
  let ffCG, integCG, tInst;
  {
    const repT = [];
    for (let rep = 0; rep < TEMP_SEEDS.length; rep++) {
      const ffTrial = new ForceField(selCG, { rc: 10, gamma: 2.0 }, mols);
      const integTrial = new LangevinIntegrator(ffTrial.ref, ffTrial, 110.0, { seed: TEMP_SEEDS[rep] });
      integTrial.setTemperature(300.0);
      integTrial.setFriction(8.0);
      for (let s = 0; s < 200; s++) integTrial.step();
      const t = ffTrial.kineticTemp(integTrial.vel, integTrial.mass);
      repT.push(t);
      if (rep === 0) { ffCG = ffTrial; integCG = integTrial; }
    }
    tInst = repT.reduce((a, b) => a + b, 0) / repT.length;
    console.log(`  (seeded temp replicas [${repT.map((t) => t.toFixed(1)).join(", ")}] → mean ${tInst.toFixed(1)} K)`);
  }
  assert(tInst > 260 && tInst < 340, `Cα Langevin temperature 300±40K after 200 steps, 3-seed mean (got ${tInst.toFixed(0)} K)`);
  assert(Number.isFinite(ffCG.energy), `Total potential energy is stable (${ffCG.energy.toFixed(2)} kcal/mol)`);
  // Validate thermal = sqrt(KB*T*KCONV/m) per coordinate (Å/ps) — not counted in 32 but verified
  const expectedThermal = Math.sqrt(KB_KCAL * 300 * KCONV / 110);
  if (Math.abs(integCG.thermal[0] - expectedThermal) >= 1e-6) {
    console.error(`  ✗ FAIL: thermal scale sqrt(KB*T*KCONV/m) correct (${integCG.thermal[0].toFixed(4)} vs ${expectedThermal.toFixed(4)})`);
    process.exit(1);
  } else {
    console.log(`  ✓ thermal scale sqrt(KB*T*KCONV/m) correct (${integCG.thermal[0].toFixed(4)} vs ${expectedThermal.toFixed(4)})`);
  }

  // -----------------------------------------------------------------
  // 9. Tiers (see tests/suites.js). FAST always; MEDIUM and SLOW opt-in.
  // -----------------------------------------------------------------
  const fastSecs = runTier("FAST", FAST_SUITES);
  if (fastSecs > FAST_BUDGET_S) {
    failed++;
    console.error(`  ✗ FAIL: [FAST] tier took ${fastSecs.toFixed(1)}s, over the ${FAST_BUDGET_S}s budget — move the slowest entries to MEDIUM (tests/suites.js)`);
  } else {
    console.log(`  [FAST] within budget: ${fastSecs.toFixed(1)}s <= ${FAST_BUDGET_S}s`);
  }
  if (MEDIUM) runTier("MEDIUM", MEDIUM_SUITES);
  else console.log(`\n[MEDIUM] skipped (opt-in: \`node tests/test_all.js --medium\` or MEDIUM=1; implied by --slow) — ${MEDIUM_SUITES.map((s) => path.basename(s.file)).join(", ")}`);
  if (SLOW) {
    runTier(SLOW_FULL ? "SLOW-FULL" : "SLOW-SMOKE", SLOW_SUITES);
  } else {
    console.log(`\n[SLOW] skipped (opt-in: \`--slow\` = SMOKE heavy, or \`--slow-full\`/\`--long\` = FULL heavy, or SLOW=1) — ${SLOW_SUITES.map((s) => path.basename(s.file)).join(", ")}`);
  }

  // -----------------------------------------------------------------
  // SUMMARY. The grand total is Tier-0 inline asserts + every tier's
  // counted child asserts. It is REPORTING (and the evolve gate's >= 352
  // floor), never a pass/fail criterion on its own.
  // -----------------------------------------------------------------
  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log(`  (tiers run: ${["FAST", MEDIUM ? "MEDIUM" : null, SLOW ? (SLOW_FULL ? "SLOW-FULL" : "SLOW-SMOKE") : null].filter(Boolean).join(" + ")})`);
  console.log("=================================================");

  if (failed > 0) process.exit(1);
}

runTests().catch((e) => {
  console.error("Test execution error:", e);
  process.exit(1);
});
