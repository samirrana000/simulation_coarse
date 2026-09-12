/**
 * tests/test_seeded_integrator.js — Stage-2 deterministic seeding regression.
 *
 * Covers the additive opt-in RNG injection on LangevinIntegrator:
 *   1. same seed → bit-identical trajectory (positions after N steps exact);
 *   2. different seeds → diverge (max |Δpos| >> 0);
 *   3. golden check: seeded 4W52 CG 200-step final energy/positions match
 *      stored constants (tolerance 1e-9; exact on the same Node build).
 *   4. API: default path has getSeed() null (Math.random, L0 bit-identical);
 *      setSeed()/setRng()/opts.seed + 3-arg overload work.
 *
 * Protocol (frozen — do not retune without regenerating the golden):
 *   4W52 Cα, ForceField({ rc: 10, gamma: 1.0 }), mass 110 Da, T=300 K,
 *   zeta=5.0 ps⁻¹, seed 12345, 200 BAOAB steps. Golden generated 2026-09-12
 *   via the same seeded constructor path (no Math.random patching).
 *
 * Runnable: node tests/test_seeded_integrator.js
 * Prints PASS/FAIL per assertion; exits 1 on any failure.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { parseCa, selectSystem } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { SeededRNG } from "../src/seeded-rng.js";

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

// Frozen golden (Stage-2, seed 12345, 200 steps, gamma 1.0, T 300, zeta 5.0).
const GOLDEN_SEED = 12345;
const GOLDEN_STEPS = 200;
const GOLDEN_ENERGY = 136.69571481112268;
const GOLDEN_HEAD = [
  -18.118594657563364,
  -2.0864837348196557,
  8.781627804927995,
  -20.643975439874588,
  -0.6213036438065017,
  11.072134592825348,
  -22.821380716267218,
  2.362892104104878,
  11.411841154109993,
];
const GOLDEN_TAIL = [
  -11.818294812747084,
  5.655695884898453,
  3.186587061064632,
];
const GOLDEN_SUM = -785.1490827088148;
const GOLDEN_SUMSQ = 188095.1882679323;

let passed = 0, failed = 0;
function assert(c, m) {
  if (c) { passed++; console.log(`  ✓ ${m}`); }
  else { failed++; console.error(`  ✗ FAIL: ${m}`); }
}

function buildSystem() {
  const pdbText = fs.readFileSync(findPdb(), "utf-8");
  return selectSystem(parseCa(pdbText));
}

/**
 * Run the frozen protocol with a given seed; return final state snapshot.
 * @param {object} sel selectSystem output
 * @param {number} seed mulberry32 seed
 * @param {number} [steps=GOLDEN_STEPS] BAOAB steps
 */
function runSeeded(sel, seed, steps = GOLDEN_STEPS) {
  const ff = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const integ = new LangevinIntegrator(ff.ref, ff, 110, { seed });
  integ.setTemperature(300);
  integ.setFriction(5.0);
  for (let s = 0; s < steps; s++) integ.step();
  return { ff, integ, pos: Float64Array.from(integ.pos), energy: ff.energy };
}

async function main() {
  console.log("=== test_seeded_integrator (Stage-2 determinism) ===");
  const sel = buildSystem();

  // 1. same seed → bit-identical trajectory
  const a = runSeeded(sel, GOLDEN_SEED);
  const b = runSeeded(sel, GOLDEN_SEED);
  let maxDiffSame = 0;
  for (let i = 0; i < a.pos.length; i++) {
    const d = Math.abs(a.pos[i] - b.pos[i]);
    if (d > maxDiffSame) maxDiffSame = d;
  }
  assert(maxDiffSame === 0, `same seed ${GOLDEN_SEED} → bit-identical positions (maxDiff ${maxDiffSame})`);
  assert(a.energy === b.energy, `same seed → bit-identical energy (${a.energy})`);

  // 2. different seeds → diverge
  const c = runSeeded(sel, 54321);
  let maxDiff = 0;
  for (let i = 0; i < a.pos.length; i++) {
    const d = Math.abs(a.pos[i] - c.pos[i]);
    if (d > maxDiff) maxDiff = d;
  }
  assert(maxDiff > 1e-3, `different seeds diverge (max |Δpos| ${maxDiff.toFixed(4)} > 1e-3)`);
  assert(a.energy !== c.energy, `different seeds → different energy (${a.energy.toFixed(4)} vs ${c.energy.toFixed(4)})`);

  // 3. golden check (tolerance 1e-9; exact on same Node build)
  const eDiff = Math.abs(a.energy - GOLDEN_ENERGY);
  assert(eDiff < 1e-9, `golden energy ${a.energy} vs ${GOLDEN_ENERGY} (diff ${eDiff}) < 1e-9`);
  let maxHeadDiff = 0;
  for (let i = 0; i < GOLDEN_HEAD.length; i++) {
    const d = Math.abs(a.pos[i] - GOLDEN_HEAD[i]);
    if (d > maxHeadDiff) maxHeadDiff = d;
  }
  assert(maxHeadDiff < 1e-9, `golden head[0..8] maxDiff ${maxHeadDiff} < 1e-9`);
  let maxTailDiff = 0;
  for (let i = 0; i < GOLDEN_TAIL.length; i++) {
    const d = Math.abs(a.pos[a.pos.length - GOLDEN_TAIL.length + i] - GOLDEN_TAIL[i]);
    if (d > maxTailDiff) maxTailDiff = d;
  }
  assert(maxTailDiff < 1e-9, `golden tail maxDiff ${maxTailDiff} < 1e-9`);
  let sum = 0, sumsq = 0;
  for (const v of a.pos) { sum += v; sumsq += v * v; }
  assert(Math.abs(sum - GOLDEN_SUM) < 1e-9, `golden pos sum ${sum} vs ${GOLDEN_SUM} < 1e-9`);
  assert(Math.abs(sumsq - GOLDEN_SUMSQ) < 1e-6, `golden pos sumsq ${sumsq} vs ${GOLDEN_SUMSQ} < 1e-6`);

  // 4. API: default unseeded path unchanged; injection forms work
  const ffD = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const integD = new LangevinIntegrator(ffD.ref, ffD, 110);
  assert(integD.getSeed() === null, "default constructor → getSeed() null (Math.random path)");
  assert(typeof integD._randUniform() === "number", "default _randUniform() returns number");

  const ffS = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const integS = new LangevinIntegrator(ffS.ref, ffS, 110);
  integS.setSeed(7);
  assert(integS.getSeed() === 7, "setSeed(7) → getSeed() 7");
  const ffO = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const integO = new LangevinIntegrator(ffO.ref, ffO, { seed: 9 });
  assert(integO.getSeed() === 9, "3-arg overload (ref, ff, { seed: 9 }) works");
  const ffR = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const integR = new LangevinIntegrator(ffR.ref, ffR, 110, { rng: new SeededRNG(11) });
  assert(typeof integR._randUniform() === "number", "opts.rng (SeededRNG instance) injects uniform source");
  const ffF = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const rngFn = new SeededRNG(13);
  const integF = new LangevinIntegrator(ffF.ref, ffF, 110, { rng: () => rngFn.rand() });
  assert(typeof integF._randUniform() === "number", "opts.rng (function) injects uniform source");
  // SeededRNG-object injection replays identically to opts.seed
  const r1 = runSeeded(sel, 77, 50);
  const ff77 = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const integ77 = new LangevinIntegrator(ff77.ref, ff77, 110, { rng: new SeededRNG(77) });
  integ77.setTemperature(300); integ77.setFriction(5.0);
  for (let s = 0; s < 50; s++) integ77.step();
  let maxRngDiff = 0;
  for (let i = 0; i < r1.pos.length; i++) {
    const d = Math.abs(r1.pos[i] - integ77.pos[i]);
    if (d > maxRngDiff) maxRngDiff = d;
  }
  assert(maxRngDiff === 0, `opts.seed ≡ opts.rng(SeededRNG same seed) after 50 steps (maxDiff ${maxRngDiff})`);

  console.log(`\n=== test_seeded_integrator: ${passed} PASSED, ${failed} FAILED ===`);
  process.exit(failed ? 1 : 0);
}

main().catch((e) => { console.error(e); process.exit(1); });
