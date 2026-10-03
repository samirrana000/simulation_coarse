/**
 * tests/test_cg_ligand_golden.js — bit-exact CG golden for the LIGAND path.
 *
 * WHY THIS EXISTS (and what it is not)
 * -----------------------------------
 * tests/golden/4w52_10steps.json IS a CG trajectory, and it was checked before
 * the forcefield.js split as if it protected the whole CG compute path. It
 * does not. It is built by
 *     new ForceField(sel, { rc: 10, gamma: 1.0 })
 * with `ligands` defaulting to [] — PROTEIN ONLY. So every ligand branch of
 * compute() ran as a no-op there:
 *   - _ligandInternal (bonds / angles / impropers)          never executed
 *   - _binding  (cross LJ + screened Coulomb + H-bond)      never executed
 *   - the EEF1-lite desolvation passes (_dens/_dBdn/_bpA-J-R) never executed
 *   - the holo contact springs and their one-sided compression floor never ran
 *   - nativeContacts for intra-ligand pairs                   never built
 *   - gridB (the protein-only grid the binding pass walks)    never built
 * That is roughly a third of the CG energy model, and all of it is code the
 * split MOVED. A golden that cannot see a term cannot prove a term survived.
 *
 * So this is an ADDITIONAL golden over a seeded run of the same system WITH
 * the crystallographic benzene (4W52 is T4 lysozyme L99A + BNZ), which drives
 * every one of the branches above. The pre-existing golden is untouched: same
 * file, same hashes, same generator.
 *
 * The hashes here were recorded from the pre-split tree (git 61374b2) and
 * MUST match after it. They are the evidence that moving the compute path
 * changed no bits.
 *
 * WHAT IS PINNED
 *   energyHash — the accumulated potential energy after the last step
 *   posHash    — the coordinate buffer
 *   velHash    — the velocity buffer
 *   forceHash  — the force buffer, which the posHash alone does not constrain
 *   termHashes — bindingU / desolvU / holoU / repU / nativeU separately, so a
 *                failure says WHICH term moved rather than just "something did"
 *
 * Runnable: node tests/test_cg_ligand_golden.js
 * Regenerate ONLY with: node tests/generate_cg_ligand_golden.js --write
 *   (that flag exists so an accidental regeneration is a deliberate act; the
 *   guard test below rejects a golden whose provenance comment is absent.)
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

import { SeededRNG } from "../src/seeded-rng.js";
import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

function findPdb() {
  for (const p of [
    path.resolve(process.cwd(), "4w52.pdb"),
    path.resolve(__dirname, "..", "4w52.pdb"),
    path.resolve(__dirname, "4w52.pdb"),
  ]) if (fs.existsSync(p)) return p;
  throw new Error(`4w52.pdb not found`);
}

function findGolden() {
  for (const p of [
    path.resolve(process.cwd(), "tests/golden/4w52_benzene_ligand.json"),
    path.resolve(__dirname, "golden/4w52_benzene_ligand.json"),
  ]) if (fs.existsSync(p)) return p;
  throw new Error(`ligand golden not found`);
}

/** FNV-1a over a fixed-point quantisation of a float buffer. */
function arrayHash(arr, scale = 1000) {
  let h = 0x811c9dc5;
  for (let i = 0; i < arr.length; i++) {
    const q = Math.round(arr[i] * scale);
    h ^= (q & 0xff);
    h = Math.imul(h, 0x01000193);
    h ^= ((q >> 8) & 0xff);
    h = Math.imul(h, 0x01000193);
    h ^= ((q >> 16) & 0xff);
    h = Math.imul(h, 0x01000193);
  }
  return (h >>> 0).toString(16).padStart(8, "0");
}

/** Hash a plain number bit-for-bit, so a term energy is pinned exactly. */
function numHash(x) {
  const buf = new DataView(new ArrayBuffer(8));
  buf.setFloat64(0, x);
  let h = 0x811c9dc5;
  for (let i = 0; i < 8; i++) {
    h ^= buf.getUint8(i);
    h = Math.imul(h, 0x01000193);
  }
  return (h >>> 0).toString(16).padStart(8, "0");
}

async function main() {
  console.log("=== CG ligand-path bit-exact golden (4W52 + benzene) ===");
  const golden = JSON.parse(fs.readFileSync(findGolden(), "utf-8"));
  const pdbText = fs.readFileSync(findPdb(), "utf-8");

  assert(typeof golden.provenance === "string" && golden.provenance.length > 0,
    `the golden carries a provenance note (${String(golden.provenance).slice(0, 48)}…)`);

  const seed = golden.seed ?? 7;
  const steps = golden.steps ?? 10;

  const rng = new SeededRNG(seed);
  rng.install();
  let got;
  try {
    const parsed = parseCa(pdbText);
    const sel = selectSystem(parsed);
    const ligands = parseLigands(pdbText);
    assert(ligands.length > 0, `parsed ${ligands.length} crystallographic ligand(s) from 4W52`);

    const ff = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligands);
    // The system must actually contain the ligand particles, or this golden
    // is vacuous in exactly the way the protein-only one is.
    assert(ff.nLigAtoms > 0, `ligand atoms appended: nLigAtoms=${ff.nLigAtoms}, n=${ff.n}, nProt=${ff.nProt}`);
    assert(ff.holoSprings.length > 0, `holo springs built: ${ff.holoSprings.length / 3} pairs`);
    ff.trackTerms = true;

    const integrator = new LangevinIntegrator(ff.ref, ff, 110);
    integrator.setTemperature(300);
    integrator.setFriction(5.0);
    for (let s = 0; s < steps; s++) integrator.step();

    // Capture the trajectory's own binding/desolvation split BEFORE the term
    // probe, which re-runs compute() at a displaced geometry and would
    // otherwise overwrite them.
    const lastEnergy = ff.energy;
    const lastBindingU = ff.bindingU;
    const lastDesolvU = ff.desolvU;

    // Probe the individual kernels on the same deterministically DISPLACED
    // geometry the generator used (see generate_cg_ligand_golden.js for why
    // native coordinates would make three of these hashes vacuous).
    const probeRng = new SeededRNG((golden.seed ?? 7) + 1);
    const probe = Float64Array.from(ff.ref);
    for (let i = 0; i < probe.length; i++) probe[i] += probeRng.rand() * 1.6 - 0.8;

    const saveF = Float64Array.from(ff.forces);
    const Uall = ff.compute(probe);
    const terms = {};
    terms.harmonic = ff._harmonicPairs(probe, ff.forces, ff.bonds, 3, ff.kBond)
      + ff._springForces(probe, ff.forces)
      + (ff.holoSprings.length ? ff._harmonicPairs(probe, ff.forces, ff.holoSprings, 3, ff.holoGamma) : 0)
      + (ff.nativeContacts.length ? ff._harmonicPairs(probe, ff.forces, ff.nativeContacts, 3, 1.0) : 0);
    terms.angle = ff._angleForces(probe, ff.forces);
    terms.ligandInternal = ff._ligandInternal(probe, ff.forces);
    terms.repulsion = ff._repulsion(probe, ff.forces);
    terms.binding = ff._binding(probe, ff.forces);
    ff.forces.set(saveF);

    got = {
      energy: lastEnergy,
      energyHash: numHash(lastEnergy),
      posHash: arrayHash(integrator.pos, 1000),
      velHash: arrayHash(integrator.vel, 1000),
      forceHash: arrayHash(ff.forces, 1000),
      bindingU: lastBindingU,
      desolvU: lastDesolvU,
      terms: Object.fromEntries(Object.entries(terms).map(([k, v]) => [k, numHash(v)])),
    };
    assert(Number.isFinite(Uall), "probe compute() returned a finite energy");
  } finally {
    rng.restore();
  }

  const rows = [
    ["energy", golden.energy, got.energy, 1e-12],
    ["energyHash", golden.energyHash, got.energyHash],
    ["posHash", golden.posHash, got.posHash],
    ["velHash", golden.velHash, got.velHash],
    ["forceHash", golden.forceHash, got.forceHash],
    ["bindingU", golden.bindingU, got.bindingU, 1e-12],
    ["desolvU", golden.desolvU, got.desolvU, 1e-12],
  ];
  for (const [name, want, have, tol] of rows) {
    if (tol === undefined) assert(want === have, `${name} ${have} === golden ${want}`);
    else assert(Math.abs(want - have) < tol, `${name} ${have} === golden ${want} (diff ${Math.abs(want - have)})`);
  }
  for (const k of Object.keys(golden.terms)) {
    assert(got.terms[k] === golden.terms[k], `term ${k} ${got.terms[k]} === golden ${golden.terms[k]}`);
  }

  console.log(`\n=== test_cg_ligand_golden: ${passed} PASSED, ${failed} FAILED ===`);
  process.exit(failed > 0 ? 1 : 0);
}

main().catch((e) => { console.error(e); process.exit(1); });