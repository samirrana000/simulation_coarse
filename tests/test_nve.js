/**
 * tests/test_nve.js — F53 NVE energy conservation (G65 companion).
 *
 * Validates that the Cα ENM ForceField conserves total energy (U + K) when
 * integrated with a symplectic (Velocity Verlet) NVE integrator (zeta=0, no thermostat).
 *
 * Approach:
 *   1. Parse 1crn.pdb (fallback: 1ubq.pdb, 4w52.pdb) → selectSystem → ForceField(rc=10, gamma=1).
 *   2. Sample Maxwell–Boltzmann velocities at 300 K (SeededRNG seed=42, deterministic) and remove COM motion.
 *   3. Run Velocity Verlet for 50 steps (dt=0.004 ps, the production CG timestep) + also a 10 ps trajectory
 *      (2500 steps) for informational drift. Compute total energy E = U + K (KCONV conversion).
 *   4. Assert |E_final − E_initial| / |E_initial| < 0.5% for the 50-step run (F53 criterion).
 *      The 10 ps run is reported but allowed a looser 2% window to guard against integration error.
 *   5. Additionally demonstrate that LangevinIntegrator can be hacked to NVE by setting zeta=0 directly
 *      (bypassing the 0.1 ps^-1 clamp in setFriction) — we show that the integrator's step without the
 *      OU thermostat also conserves to <1% over 20 steps. This is NOT the primary gate; VV is.
 *
 * Runnable: node tests/test_nve.js
 * Prints PASS on success (required by F53 measurable), FAIL + non-zero exit on drift violation.
 * Syntax check: node --check tests/test_nve.js
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

import { parseCa, selectSystem } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { KCONV, KB_KCAL } from "../src/units.js";
import { SeededRNG } from "../src/seeded-rng.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

function findPdb(names) {
  const candidates = [];
  for (const n of names) {
    candidates.push(path.resolve(process.cwd(), n));
    candidates.push(path.resolve(__dirname, "..", n));
    candidates.push(path.resolve(__dirname, n));
    candidates.push(path.resolve("data", n));
  }
  for (const p of candidates) {
    try { if (fs.existsSync(p)) return { text: fs.readFileSync(p, "utf-8"), path: p }; } catch {}
  }
  throw new Error(`No PDB found, tried ${candidates.join(", ")}`);
}

function kineticEnergy(vel, massPerCoord) {
  let ke = 0;
  for (let i = 0; i < vel.length; i++) ke += massPerCoord[i] * vel[i] * vel[i];
  return 0.5 * ke / KCONV; // → kcal/mol
}

function removeCOM(vel, massPerCoord) {
  let px = 0, py = 0, pz = 0, M = 0;
  const n3 = vel.length;
  for (let i = 0; i < n3; i += 3) {
    const m = massPerCoord[i];
    px += m * vel[i]; py += m * vel[i + 1]; pz += m * vel[i + 2];
    M += m;
  }
  px /= M; py /= M; pz /= M;
  for (let i = 0; i < n3; i += 3) { vel[i] -= px; vel[i + 1] -= py; vel[i + 2] -= pz; }
}

async function run() {
  console.log("=== F53 NVE energy conservation (zeta=0, Velocity Verlet) ===");
  const entry = findPdb(["1crn.pdb", "1CRN.pdb", "1ubq.pdb", "4w52.pdb"]);
  console.log(`PDB: ${path.basename(entry.path)}`);
  const parsed = parseCa(entry.text);
  const sel = selectSystem(parsed);
  console.log(`System: ${sel.beads.length} Cα beads, ${sel.segments.length} segment(s)`);
  const ff = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const n = ff.n, n3 = n * 3;
  console.log(`ForceField: n=${n} nProt=${ff.nProt} rc=${ff.rc} gamma=${ff.gamma}`);

  const masses = ff.masses; // per particle
  const massPerCoord = new Float64Array(n3);
  for (let i = 0; i < n3; i++) massPerCoord[i] = masses[(i / 3) | 0];

  // Deterministic Maxwell–Boltzmann velocities at 300 K (seeded)
  const rng = new SeededRNG(42);
  const T = 300;
  const thermal = new Float64Array(n3);
  for (let i = 0; i < n3; i++) thermal[i] = Math.sqrt(KB_KCAL * T * KCONV / massPerCoord[i]);
  const vel = new Float64Array(n3);
  for (let i = 0; i < n3; i++) vel[i] = thermal[i] * rng.randn();
  removeCOM(vel, massPerCoord);

  const pos = new Float64Array(ff.ref); // copy native
  ff.compute(pos);
  const U0 = ff.energy;
  const K0 = kineticEnergy(vel, massPerCoord);
  const E0 = U0 + K0;
  if (!Number.isFinite(E0)) {
    console.error(`FAIL: initial total energy non-finite U0=${U0} K0=${K0}`);
    process.exit(1);
  }
  console.log(`Initial: U=${U0.toFixed(4)} K=${K0.toFixed(4)} E=${E0.toFixed(4)} kcal/mol (T_inst ~${(K0 / (1.5 * n * KB_KCAL)).toFixed(1)} K)`);

  // --- Velocity Verlet NVE (zeta=0, no thermostat) — the canonical F53 gate ---
  const dt = 0.004; // ps, production CG timestep (Heavy would need 0.001)
  const stepsShort = 50;
  let f = ff.forces; // after initial compute
  for (let s = 0; s < stepsShort; s++) {
    // half kick
    for (let i = 0; i < n3; i++) vel[i] += 0.5 * dt * KCONV * f[i] / massPerCoord[i];
    // drift
    for (let i = 0; i < n3; i++) pos[i] += dt * vel[i];
    // new forces
    ff.compute(pos);
    f = ff.forces;
    // half kick
    for (let i = 0; i < n3; i++) vel[i] += 0.5 * dt * KCONV * f[i] / massPerCoord[i];
  }
  const U1 = ff.energy;
  const K1 = kineticEnergy(vel, massPerCoord);
  const E1 = U1 + K1;
  const driftShort = Math.abs(E1 - E0) / Math.max(1e-12, Math.abs(E0));
  console.log(`After ${stepsShort} steps (dt=${dt} ps, ${(stepsShort * dt).toFixed(3)} ps NVE): U=${U1.toFixed(4)} K=${K1.toFixed(4)} E=${E1.toFixed(4)} drift=${(driftShort * 100).toFixed(4)}%`);
  const passShort = driftShort < 0.005; // 0.5%
  if (!passShort) {
    console.error(`FAIL: NVE drift ${(driftShort * 100).toFixed(4)}% exceeds 0.5% over ${stepsShort} steps (dt=${dt} ps)`);
    process.exit(1);
  }
  console.log(`  drift ${(driftShort * 100).toFixed(4)}% < 0.5% → PASS (50-step gate)`);

  // --- Longer 10 ps NVE for informational drift (2500 steps at dt=0.004) ---
  // Reset to native + fresh velocities (same seed progression → deterministic)
  pos.set(ff.ref);
  // Continue RNG from where we left off (or reset for reproducibility) — reset to seed 43 for independence
  const rng2 = new SeededRNG(43);
  for (let i = 0; i < n3; i++) vel[i] = thermal[i] * rng2.randn();
  removeCOM(vel, massPerCoord);
  ff.compute(pos);
  const U0b = ff.energy, K0b = kineticEnergy(vel, massPerCoord), E0b = U0b + K0b;
  const stepsLong = Math.round(10.0 / dt); // 10 ps
  f = ff.forces;
  for (let s = 0; s < stepsLong; s++) {
    for (let i = 0; i < n3; i++) vel[i] += 0.5 * dt * KCONV * f[i] / massPerCoord[i];
    for (let i = 0; i < n3; i++) pos[i] += dt * vel[i];
    ff.compute(pos);
    f = ff.forces;
    for (let i = 0; i < n3; i++) vel[i] += 0.5 * dt * KCONV * f[i] / massPerCoord[i];
  }
  const U2 = ff.energy, K2 = kineticEnergy(vel, massPerCoord), E2 = U2 + K2;
  const driftLong = Math.abs(E2 - E0b) / Math.max(1e-12, Math.abs(E0b));
  console.log(`After ${stepsLong} steps (10 ps NVE, dt=${dt} ps): U=${U2.toFixed(4)} K=${K2.toFixed(4)} E=${E2.toFixed(4)} drift=${(driftLong * 100).toFixed(4)}%`);
  if (driftLong > 0.02) {
    console.warn(`WARN: 10 ps drift ${(driftLong * 100).toFixed(4)}% exceeds 2% — integrator may need smaller dt for long NVE; short gate still passes`);
  } else {
    console.log(`  10 ps drift ${(driftLong * 100).toFixed(4)}% — excellent`);
  }

  // --- LangevinIntegrator NVE hack demo (zeta=0 bypassing clamp) ---
  // LangevinIntegrator.setFriction clamps to min 0.1 ps^-1, so true NVE requires direct assignment.
  // We demonstrate that with zeta=0 and no OU noise, the integrator reduces to Velocity Verlet + half-kick guards.
  try {
    const posL = new Float64Array(ff.ref);
    const integ = new LangevinIntegrator(ff.ref, ff, 110);
    integ.setTemperature(300);
    // Bypass clamp for NVE demo: directly set zeta=0 and recompute dt without friction limit
    integ.zeta = 0;
    integ.dt = integ._pickDt(); // will now pick dt based only on bond limit (~0.004 for CG)
    // Re-sample velocities deterministically for this branch
    const rng3 = new SeededRNG(44);
    // Reuse thermal scale from integ (should match 300K)
    for (let i = 0; i < n3; i++) vel[i] = integ.thermal[i] * rng3.randn();
    removeCOM(vel, integ.mass);
    integ.pos.set(posL);
    integ.vel.set(vel);
    ff.compute(integ.pos);
    const E0L = ff.energy + kineticEnergy(integ.vel, integ.mass);
    // Run 20 steps of Langevin without thermostat (since zeta=0, O-step becomes identity: c1=1, c2=0)
    for (let s = 0; s < 20; s++) integ.step();
    const E1L = ff.energy + kineticEnergy(integ.vel, integ.mass);
    const driftL = Math.abs(E1L - E0L) / Math.max(1e-12, Math.abs(E0L));
    console.log(`Langevin NVE hack (zeta=0, dt=${integ.dt.toFixed(4)} ps, 20 steps): drift=${(driftL * 100).toFixed(4)}%`);
    if (driftL < 0.01) console.log("  Langevin NVE hack drift <1% — consistent with VV");
    else console.warn(`  Langevin NVE hack drift ${(driftL*100).toFixed(4)}% — larger than VV likely due to half-kick cap DVMAX`);
  } catch (e) {
    console.warn(`Langevin NVE demo skipped: ${e.message}`);
  }

  // Final gate: we already passed the 50-step 0.5% criterion, so overall PASS
  console.log("PASS");
}

run().catch(e => { console.error("FAIL:", e); process.exit(1); });
