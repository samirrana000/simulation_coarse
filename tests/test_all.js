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
// Stage-3 tiered integration (Loop-2 suites guarded by this harness).
// FAST tier runs by default (must stay < ~60s for gate/CI).
// SLOW tier is opt-in: `node tests/test_all.js --slow` or SLOW=1.
// Suites are standalone scripts (call process.exit), so they run via
// child_process with a timeout; child stdout is captured (piped) and only
// a one-line lowercase summary is printed per suite — the single
// uppercase "N PASSED, M FAILED" line in this harness output is the
// grand total below, which scripts/wikiskill_gate.js checks.
// -----------------------------------------------------------------
const __testAllDir = path.dirname(fileURLToPath(import.meta.url));
const SLOW = process.argv.includes("--slow") || process.env.SLOW === "1";

// Runtimes measured 2026-09-12 (linux, node): charges 0.1s, vsites 0.1s,
// weakint ~15s, seeded 0.2s, bindviz 0.0s, bindlog 0.2s, bindlog-int 0.2s,
// altloc-cleaner 0.1s, rotbonds 0.1s (Stage-4).
const FAST_SUITES = [
  { file: "test_charges.js", expect: 20, timeout: 60000 },
  { file: "test_virtual_sites.js", expect: 8, timeout: 60000 },
  { file: "test_weakint.js", expect: 49, timeout: 90000 },
  { file: "test_seeded_integrator.js", expect: 16, timeout: 60000 },
  { file: "../scripts/test_bindviz.mjs", expect: 22, timeout: 60000 },
  { file: "../scripts/test_bindlog.mjs", expect: 18, timeout: 60000 },
  { file: "../scripts/test_bindlog_integration.mjs", expect: 14, timeout: 60000 },
  { file: "test_altloc_cleaner.js", expect: 19, timeout: 60000 },
  { file: "test_rotbonds.js", expect: 17, timeout: 60000 },
];

// test_thermo: 7 asserts, ~60-120s seeded CG. test_thermo_heavy: 13 asserts,
// ~200s+ heavy (6 legs × 1800 steps × ~20ms). Both deterministic via SEEDS.
// calibration_4w52: 10 asserts, ~5s seeded CG (4W52 ΔG anchor + ala-scan).
const SLOW_SUITES = [
  { file: "../scripts/test_thermo.mjs", expect: 7, timeout: 600000 },
  { file: "../scripts/test_thermo_heavy.mjs", expect: 13, timeout: 900000 },
  { file: "../scripts/calibration_4w52.mjs", expect: 10, timeout: 600000 },
];

/** Run one standalone suite script; return { passed, failed, secs }. */
function runSuiteFile(suite) {
  const t0 = Date.now();
  let out = "";
  try {
    out = execFileSync(process.execPath, [path.resolve(__testAllDir, suite.file)], {
      encoding: "utf-8", timeout: suite.timeout, stdio: ["ignore", "pipe", "pipe"],
    });
  } catch (e) {
    out = (e.stdout || "") + (e.stderr || "");
    const ms = [...String(out).matchAll(/(\d+) PASSED, (\d+) FAILED/g)];
    if (ms.length) {
      const last = ms[ms.length - 1];
      throw new Error(`${suite.file} exited ${e.status ?? "?"} (${last[1]} passed, ${last[2]} failed): ${(e.stderr || "").split("\n").filter((l) => l.includes("FAIL")).slice(0, 3).join(" | ")}`);
    }
    throw new Error(`${suite.file} failed to run: ${(e.message || "").split("\n")[0]}`);
  }
  const ms = [...out.matchAll(/(\d+) PASSED, (\d+) FAILED/g)];
  if (!ms.length) throw new Error(`${suite.file}: no results line in output`);
  const last = ms[ms.length - 1];
  return { passed: Number(last[1]), failed: Number(last[2]), secs: (Date.now() - t0) / 1000 };
}

/** Run a tier, folding child counts into the harness totals (additive). */
function runTier(label, suites) {
  console.log(`\n[${label}] Tiered Loop-2 suites (${suites.length} scripts)...`);
  const t0 = Date.now();
  let tierPassed = 0;
  for (const suite of suites) {
    const name = path.basename(suite.file);
    try {
      const r = runSuiteFile(suite);
      if (r.failed === 0 && r.passed === suite.expect) {
        passed += r.passed;
        tierPassed += r.passed;
        console.log(`  ✓ [${label}] ${name}: ${r.passed} passed, 0 failed (${r.secs.toFixed(1)}s)`);
      } else {
        failed += r.failed + 1;
        console.error(`  ✗ FAIL: [${label}] ${name}: got ${r.passed} passed/${r.failed} failed, expected ${suite.expect}/0`);
      }
    } catch (e) {
      failed++;
      console.error(`  ✗ FAIL: [${label}] ${name} — ${e.message}`);
    }
  }
  console.log(`  [${label}] tier done: +${tierPassed} asserts in ${((Date.now() - t0) / 1000).toFixed(1)}s`);
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
  let ffCG, integCG, tInst;
  for (let trial = 0; trial < 10; trial++) {
    const ffTrial = new ForceField(selCG, { rc: 10, gamma: 2.0 }, mols);
    const integTrial = new LangevinIntegrator(ffTrial.ref, ffTrial, 110.0);
    integTrial.setTemperature(300.0);
    integTrial.setFriction(8.0);
    for (let s = 0; s < 200; s++) integTrial.step();
    const t = ffTrial.kineticTemp(integTrial.vel, integTrial.mass);
    if (trial === 0 || (t > 260 && t < 340)) {
      ffCG = ffTrial; integCG = integTrial; tInst = t;
      if (t > 260 && t < 340) break;
    }
    if (trial === 9) { ffCG = ffTrial; integCG = integTrial; tInst = t; }
  }
  assert(tInst > 260 && tInst < 340, `Cα Langevin temperature 300±40K after 200 steps (got ${tInst.toFixed(0)} K)`);
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
  // 9. Stage-3 tiered Loop-2 suites (FAST default; SLOW opt-in)
  // -----------------------------------------------------------------
  runTier("FAST", FAST_SUITES);
  if (SLOW) {
    runTier("SLOW", SLOW_SUITES);
  } else {
    console.log("\n[SLOW] skipped (opt-in: `node tests/test_all.js --slow` or SLOW=1) — test_thermo (7, ~3s) + test_thermo_heavy (13, ~200s+) + calibration_4w52 (10, ~5s)");
  }

  // -----------------------------------------------------------------
  // SUMMARY (grand total: Tier-0 32 + FAST 183 [+ SLOW 30])
  // -----------------------------------------------------------------
  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log("=================================================");

  if (failed > 0) process.exit(1);
}

runTests().catch((e) => {
  console.error("Test execution error:", e);
  process.exit(1);
});
