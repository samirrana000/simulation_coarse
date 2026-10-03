/**
 * tests/generate_cg_ligand_golden.js — records the CG LIGAND-path golden.
 *
 * NOT PART OF npm test. tests/suites.js lists test_cg_ligand_golden.js (the
 * CONSUMER); this is its producer and is inventoried as a tool, alongside
 * generate_golden.js.
 *
 * The recorded numbers are the contract tests/test_cg_ligand_golden.js checks.
 * They were recorded from the PRE-SPLIT tree, so if a later refactor moves the
 * ligand kernels and changes a single bit, that test goes red — which is the
 * whole point of having a ligand golden at all.
 *
 * The --write flag is mandatory and deliberate: regenerating a golden to make a
 * red test green is the failure mode this file exists to make hard to do by
 * accident. Without it the script only prints what it WOULD record and exits
 * non-zero if that differs from what is on disk.
 *
 * Usage:
 *   node tests/generate_cg_ligand_golden.js --write    # record
 *   node tests/generate_cg_ligand_golden.js            # dry-run / verify
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

import { SeededRNG } from "../src/seeded-rng.js";
import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const OUT = path.resolve(__dirname, "golden/4w52_benzene_ligand.json");
const write = process.argv.includes("--write");

function arrayHash(arr, scale = 1000) {
  let h = 0x811c9dc5;
  for (let i = 0; i < arr.length; i++) {
    const q = Math.round(arr[i] * scale);
    h ^= (q & 0xff); h = Math.imul(h, 0x01000193);
    h ^= ((q >> 8) & 0xff); h = Math.imul(h, 0x01000193);
    h ^= ((q >> 16) & 0xff); h = Math.imul(h, 0x01000193);
  }
  return (h >>> 0).toString(16).padStart(8, "0");
}
function numHash(x) {
  const buf = new DataView(new ArrayBuffer(8));
  buf.setFloat64(0, x);
  let h = 0x811c9dc5;
  for (let i = 0; i < 8; i++) { h ^= buf.getUint8(i); h = Math.imul(h, 0x01000193); }
  return (h >>> 0).toString(16).padStart(8, "0");
}

const SEED = 7, STEPS = 10;
const pdbPath = [
  path.resolve(__dirname, "..", "4w52.pdb"),
  path.resolve(process.cwd(), "4w52.pdb"),
].find((p) => fs.existsSync(p));
if (!pdbPath) throw new Error("4w52.pdb not found");
const pdbText = fs.readFileSync(pdbPath, "utf-8");

const rng = new SeededRNG(SEED);
rng.install();
let out;
try {
  const parsed = parseCa(pdbText);
  const sel = selectSystem(parsed);
  const ligands = parseLigands(pdbText);
  const ff = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligands);
  ff.trackTerms = true;
  const integrator = new LangevinIntegrator(ff.ref, ff, 110);
  integrator.setTemperature(300);
  integrator.setFriction(5.0);
  for (let s = 0; s < STEPS; s++) integrator.step();

  // Capture the trajectory's own binding/desolvation split BEFORE the probe
  // below: the probe re-runs compute() at a displaced geometry and would
  // otherwise overwrite ff.bindingU / ff.desolvU, silently redefining what
  // these two numbers mean.
  const lastEnergy = ff.energy;
  const lastBindingU = ff.bindingU;
  const lastDesolvU = ff.desolvU;

  // Probe the individual kernels on a DISPLACED geometry, not the native
  // reference. At native coordinates the ENM springs, the angles and the
  // excluded volume all sit at their zero-energy points, so their hashes would
  // be the hash of 0.0 — a vacuous pin that passes no matter what those three
  // kernels compute. A fixed SeededRNG kick (never Math.random) puts every
  // term off its minimum at once, so all five hashes are live constraints.
  const probeRng = new SeededRNG(SEED + 1);
  const probe = Float64Array.from(ff.ref);
  for (let i = 0; i < probe.length; i++) probe[i] += probeRng.rand() * 1.6 - 0.8;

  const saveF = Float64Array.from(ff.forces);
  const terms = {
    harmonic: ff._harmonicPairs(probe, ff.forces, ff.bonds, 3, ff.kBond)
      + ff._springForces(probe, ff.forces)
      + (ff.holoSprings.length ? ff._harmonicPairs(probe, ff.forces, ff.holoSprings, 3, ff.holoGamma) : 0)
      + (ff.nativeContacts.length ? ff._harmonicPairs(probe, ff.forces, ff.nativeContacts, 3, 1.0) : 0),
    angle: ff._angleForces(probe, ff.forces),
    ligandInternal: ff._ligandInternal(probe, ff.forces),
    repulsion: ff._repulsion(probe, ff.forces),
    binding: ff._binding(probe, ff.forces),
  };
  ff.forces.set(saveF);
  for (const [k, v] of Object.entries(terms)) {
    if (!Number.isFinite(v) || v === 0) throw new Error(`probe term ${k} is vacuous (${v}) — the golden would pin nothing`);
  }

  out = {
    provenance: "Recorded from the PRE-SPLIT tree at git 61374b2 (src/forcefield.js, 947 LOC), before the src/cg/ extraction. Consumed by tests/test_cg_ligand_golden.js. tests/golden/4w52_10steps.json is protein-only and cannot see the ligand kernels; this one can. Regenerate with: node tests/generate_cg_ligand_golden.js --write",
    system: "4W52 — T4 lysozyme L99A + benzene (BNZ), Cα CG + united-atom ligand",
    seed: SEED, steps: STEPS, dtPs: 0.004,
    n: ff.n, nProt: ff.nProt, nLigAtoms: ff.nLigAtoms,
    nHolo: ff.nHolo, nNativeContacts: ff.nativeContacts.length / 3,
    energy: lastEnergy, energyHash: numHash(lastEnergy),
    posHash: arrayHash(integrator.pos, 1000),
    velHash: arrayHash(integrator.vel, 1000),
    forceHash: arrayHash(ff.forces, 1000),
    bindingU: lastBindingU, desolvU: lastDesolvU,
    terms: Object.fromEntries(Object.entries(terms).map(([k, v]) => [k, numHash(v)])),
  };
} finally {
  rng.restore();
}

if (write) {
  fs.writeFileSync(OUT, JSON.stringify(out, null, 2) + "\n");
  console.log(`WROTE ${path.relative(path.resolve(__dirname, ".."), OUT)}`);
} else {
  const prev = fs.existsSync(OUT) ? JSON.parse(fs.readFileSync(OUT, "utf-8")) : null;
  const a = JSON.stringify(prev), b = JSON.stringify(out);
  console.log(a === b ? "DRY-RUN: identical to the recorded golden" : "DRY-RUN: DIFFERS from the recorded golden");
  console.log(JSON.stringify(out, null, 2));
  if (a !== b) process.exit(1);
}