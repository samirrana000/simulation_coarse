/**
 * tests/generate_golden.js — Deterministic golden-file generator for Cα Langevin regression.
 *
 * Approach (documented for reproducibility):
 *   1. Read 4w52.pdb (164 Cα beads) → parseCa → selectSystem → ForceField({rc:10, gamma:1})
 *   2. Seed a mulberry32 PRNG with fixed seed=42 via src/seeded-rng.js SeededRNG.
 *      The generator's Math.random is *monkey-patched* for the duration of the
 *      run so that *all* Langevin stochasticity (initial Maxwell–Boltzmann sampling
 *      and every BAOAB O-step Gaussian draw) is deterministic. This is done via
 *      rng.install() / rng.restore() which saves and replaces global Math.random.
 *      The integrator itself is unmodified — it continues to call Math.random()
 *      internally via _fillGaussian(); the patch makes those calls replay identically.
 *   3. Construct LangevinIntegrator(ref, ff) — initial velocities sampled deterministically,
 *      compute initial energy, then advance exactly 10 BAOAB steps at T=300K, zeta=5 ps^-1.
 *   4. Capture final pos (Float64Array 3n), vel (3n), energy (kcal/mol), dt, time.
 *   5. Compute lightweight hashes posHash / velHash (32-bit FNV-inspired over rounded
 *      coordinates) so the JSON can be diffed without storing multi-KB arrays, while
 *      still optionally embedding trimmed arrays for manual inspection.
 *   6. Write tests/golden/4w52_10steps.json with {n, energy, posHash, velHash, seed, steps, dt, ...}.
 *
 * Determinism guarantee: re-running `node tests/generate_golden.js` on the same
 * Node version and source must produce byte-identical posHash/velHash because
 * only integer mulberry32 → Box-Muller → integrator math (all deterministic
 * IEEE-754) is involved. If the force field or integrator changes, the golden
 * file will mismatch and the regression suite should fail.
 *
 * Usage: node tests/generate_golden.js
 * Output: tests/golden/4w52_10steps.json (created, valid JSON)
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

import { SeededRNG } from "../src/seeded-rng.js";
import { parseCa, selectSystem } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

function findPdb() {
  const candidates = [
    path.resolve(process.cwd(), "4w52.pdb"),
    path.resolve(__dirname, "..", "4w52.pdb"),
    path.resolve(__dirname, "4w52.pdb"),
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`4w52.pdb not found (tried ${candidates.join(", ")})`);
}

function findOutput() {
  const candidates = [
    path.resolve(process.cwd(), "tests/golden/4w52_10steps.json"),
    path.resolve(__dirname, "golden/4w52_10steps.json"),
    path.resolve(__dirname, "..", "tests/golden/4w52_10steps.json"),
  ];
  // prefer first that is inside an existing tests/golden dir
  for (const p of candidates) {
    const dir = path.dirname(p);
    if (fs.existsSync(dir)) return p;
  }
  // fallback to cwd one and ensure dir exists
  const fallback = path.resolve(process.cwd(), "tests/golden/4w52_10steps.json");
  fs.mkdirSync(path.dirname(fallback), { recursive: true });
  return fallback;
}

// Simple 32-bit hash over float array (FNV-1a inspired on quantized ints)
function arrayHash(arr, scale = 1000) {
  let h = 0x811c9dc5; // FNV offset
  for (let i = 0; i < arr.length; i++) {
    const q = Math.round(arr[i] * scale); // quantize to 0.001 Å / Å/ps
    h ^= (q & 0xff);
    h = Math.imul(h, 0x01000193);
    h ^= ((q >> 8) & 0xff);
    h = Math.imul(h, 0x01000193);
    h ^= ((q >> 16) & 0xff);
    h = Math.imul(h, 0x01000193);
  }
  return (h >>> 0).toString(16).padStart(8, "0");
}

function roundArray(arr, decimals = 6) {
  const out = new Array(arr.length);
  for (let i = 0; i < arr.length; i++) out[i] = Number(arr[i].toFixed(decimals));
  return out;
}

async function main() {
  const pdbPath = findPdb();
  const outPath = findOutput();
  const pdbText = fs.readFileSync(pdbPath, "utf-8");

  const seed = 42;
  const rng = new SeededRNG(seed);
  rng.install();
  try {
    const parsed = parseCa(pdbText);
    const sel = selectSystem(parsed);
    const ff = new ForceField(sel, { rc: 10, gamma: 1.0 });
    // Integrator will sample initial velocities deterministically via patched Math.random
    const integrator = new LangevinIntegrator(ff.ref, ff, 110);
    integrator.setTemperature(300);
    integrator.setFriction(5.0);

    const STEPS = 10;
    for (let s = 0; s < STEPS; s++) integrator.step();

    const pos = integrator.pos;
    const vel = integrator.vel;
    const energy = ff.energy;
    const n = ff.n;
    const dt = integrator.dt;
    const timePs = integrator.time;

    const posHash = arrayHash(pos, 1000);
    const velHash = arrayHash(vel, 1000);
    const energyHash = arrayHash(new Float64Array([energy]), 1000);

    // Trimmed arrays for inspection: first 9 coords + full length hashes
    const posTrim = roundArray(pos.slice(0, Math.min(pos.length, 30)), 6);
    const velTrim = roundArray(vel.slice(0, Math.min(vel.length, 30)), 6);
    const posFull = roundArray(pos, 4); // 4 decimals to keep file readable but deterministic

    const payload = {
      _comment: "Golden regression for 4W52 Cα Langevin (seed=42 deterministic via mulberry32 patched Math.random). Regenerate with: node tests/generate_golden.js",
      n,
      nProt: ff.nProt,
      seed,
      steps: STEPS,
      dtPs: dt,
      timePs,
      energy,
      energyHash,
      posHash,
      velHash,
      // Primary regression keys (required by spec)
      // plus trimmed arrays for manual inspection
      pos: posFull, // full pos array rounded to 4 decimals (use posHash for fast diff)
      posTrim,
      velTrim,
      pdb: path.basename(pdbPath),
      generator: "tests/generate_golden.js via SeededRNG(mulberry32) seed=42 patching Math.random",
      timestamp: new Date().toISOString(),
      nodeVersion: process.version,
    };

    fs.mkdirSync(path.dirname(outPath), { recursive: true });
    fs.writeFileSync(outPath, JSON.stringify(payload, null, 2) + "\n", "utf-8");
    console.log(`Golden written to ${outPath}`);
    console.log(`n=${n} energy=${energy.toFixed(4)} posHash=${posHash} velHash=${velHash} seed=${seed} steps=${STEPS}`);
    console.log(JSON.stringify({ n, energy, posHash, velHash, seed, steps: STEPS, dt }, null, 2));
  } finally {
    rng.restore();
  }
}

main().catch((e) => { console.error(e); process.exit(1); });
