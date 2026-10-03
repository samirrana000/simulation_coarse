/**
 * tests/test_heavy_golden.js — BIT-EXACT regression net for heavy mode.
 *
 * Run: node tests/test_heavy_golden.js        (also in tests/suites.js, FAST)
 *
 * WHY THIS EXISTS
 * ---------------
 * tests/test_golden.js pins the CG (Cα ENM) trajectory and imports nothing from
 * src/heavy.js, so the heavy engine had NO bit-level regression net — only
 * finite-difference tests, which assert a 1e-3 RELATIVE tolerance on 10 atoms
 * at one pose. A 1e-3 tolerance cannot see a summation reorder that changes the
 * last bit of a force sum, which is exactly the failure mode a "just move this
 * kernel" refactor can introduce. This test compares 4-lane FNV-1a hashes over
 * raw IEEE-754 bytes: equal hash == bit-identical, nothing softer.
 *
 * WHAT IS PINNED (see tests/generate_heavy_golden.js for the full rationale)
 *   [0] the fixture exists, is well-formed, and was generated from THIS system
 *       (n / topology counts / atom bookkeeping)
 *   [1] topology bit-identity: bonds, angles, propers, impropers, coord,
 *       ligandBonds, ref, masses
 *   [2] 6 constructor configurations x 12 poses: energy + force bits + 24
 *       per-term accumulators each. The configs switch on disjoint kernels
 *       (harmonicFlat / *PerK / _nonBondedGrid / _nonBondedGridOBC2 /
 *       _nonBondedGridNoGB / SasaModel / lcpoSasa / weak / membrane), so the
 *       set covers every branch of compute().
 *   [3] a 10-step Langevin trajectory, seed 42 — the heavy analogue of
 *       tests/test_golden.js
 *
 * IT COMPARES; IT NEVER REGENERATES. tests/generate_heavy_golden.js is not in
 * suites.js and is not a test. Editing the fixture by hand defeats the point;
 * regenerate only when a physics change is intended, and say so in the commit.
 *
 * KNOWN PRE-EXISTING DEFECT, PINNED ON PURPOSE
 *   The `full` configuration has par.gbModel === "obc2" AND trackTerms === true.
 *   In that combination _nonBondedGridNoGB reads `ligStart` before its `const`
 *   declaration (temporal dead zone) and throws ReferenceError; compute() catches
 *   it and falls back to _nonBondedGrid, but the OBC2 forces were already merged
 *   into `f` before the throw, so the reported forces are OBC2-GB + HCT-everything
 *   while the reported energy is pure HCT. That is a real physics bug, present
 *   before this module split and deliberately NOT fixed by it (a fix is a
 *   physics change). It is pinned here so that fixing it turns this test RED
 *   and the fix cannot be mistaken for a refactor.
 *
 * FAST tier: 6 constructors + 72 compute() calls + 10 integrator steps, no DOM.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { parseHeavy, selectHeavy, appendHeavyLigands, HeavyForceField } from "../src/heavy.js";
import { parseMol2 } from "../src/mol2.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { SeededRNG } from "../src/seeded-rng.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const FIXTURE = path.join(ROOT, "tests", "golden", "heavy_4w52_fp.json");

/** 4-lane FNV-1a over raw bytes -> 8 hex chars per lane. Bit-exact. */
function hashBytes(u8) {
  const L = [0x811c9dc5, 0x01000193 ^ 0x9e3779b9, 0x811c9dc5 ^ 0x85ebca6b, 0x27d4eb2f];
  for (let i = 0; i < u8.length; i++) {
    const b = u8[i];
    for (let k = 0; k < 4; k++) {
      L[k] ^= (b + k * 61) & 0xff;
      L[k] = Math.imul(L[k], 0x01000193) >>> 0;
    }
  }
  return L.map((x) => (x >>> 0).toString(16).padStart(8, "0")).join("");
}
const _f64 = new Float64Array(1);
const _u8 = new Uint8Array(_f64.buffer);
const hashNum = (x) => { _f64[0] = x; return hashBytes(_u8); };
const hashArr = (a) => hashBytes(new Uint8Array(a.buffer, a.byteOffset, a.byteLength));

const TERMS = [
  "energy", "bindingU", "desolvU", "repU", "bondU", "angleU", "improperU",
  "properU", "coordU", "elecU", "gbU", "sasaU", "hbondU", "membraneU",
  "springU", "weakU", "piU", "cpiU", "xbU", "coordAngleU",
  "bindLJU", "bindCoulU", "bindHBU", "_nanStrikes",
];

const CONFIGS = {
  legacy: { gamma: 2.0, temp: 300 },
  amber14: { gamma: 2.0, temp: 300, useAmber14: true },
  obc2: { gamma: 2.0, temp: 300, gbModel: "obc2" },
  lcpo: { gamma: 2.0, temp: 300, sasaModel: "lcpo" },
  full: {
    gamma: 2.0, temp: 300, useAmber14: true, gbModel: "obc2", sasaModel: "lcpo",
    weak: "on", metalAngles: true, membrane: { on: true, thickness: 15, width: 2 },
    trackTerms: true, binding: { on: true, charges: true, hbMode: "directional" },
  },
};

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

function buildSel() {
  const parsed = parseHeavy(fs.readFileSync(path.join(ROOT, "4w52.pdb"), "utf-8"));
  const mols = parseMol2(fs.readFileSync(path.join(ROOT, "benzene.mol2"), "utf-8"));
  let sel = selectHeavy(parsed, {
    heteroSelection: { "A|200|BNZ": false, "A|201|EPE": false },
    includePdbLigands: true,
    hasExternalLigand: true,
  });
  sel = appendHeavyLigands(sel, mols);
  return sel;
}

function poses(ref, count = 12) {
  const out = [new Float64Array(ref)];
  const rng = new SeededRNG(1234);
  for (let p = 1; p < count; p++) {
    const q = new Float64Array(ref);
    for (let i = 0; i < q.length; i++) q[i] += (rng.rand() * 2 - 1) * 0.05;
    out.push(q);
  }
  return out;
}

function main() {
  console.log("=== heavy-mode bit-exact golden (4W52 + benzene MOL2) ===");

  // ---- [0] fixture -------------------------------------------------------
  const hasFix = fs.existsSync(FIXTURE);
  assert(hasFix, `fixture tests/golden/heavy_4w52_fp.json exists`);
  if (!hasFix) { finish(); return; }
  const golden = JSON.parse(fs.readFileSync(FIXTURE, "utf-8"));
  assert(Object.keys(CONFIGS).every((k) => golden.configs && golden.configs[k]),
    `fixture carries all ${Object.keys(CONFIGS).length} configurations`);
  assert(golden.traj && golden.topo, "fixture carries the topology + trajectory rows");

  const sel = buildSel();

  // ---- [1] topology ------------------------------------------------------
  console.log("\n[1] topology bit-identity...");
  {
    const ff = new HeavyForceField({ atoms: sel.atoms }, CONFIGS.full, []);
    const t = golden.topo;
    assert(ff.bonds.length === t.counts.bonds && ff.angles.length === t.counts.angles &&
      ff.propers.length === t.counts.propers && ff.impropers.length === t.counts.impropers,
    `topology sizes ${ff.bonds.length / 3} bonds / ${ff.angles.length / 4} angles / ` +
    `${ff.propers.length / 5} propers / ${ff.impropers.length / 5} impropers`);
    for (const k of ["bonds", "angles", "propers", "impropers", "coord", "ligandBonds", "ref", "masses"]) {
      const h = hashArr(ff[k]);
      assert(h === t[k], `${k} bit-identical (${h}${h === t[k] ? "" : ` != ${t[k]}`})`);
    }
  }

  // ---- [2] configurations x poses ---------------------------------------
  for (const [name, par] of Object.entries(CONFIGS)) {
    console.log(`\n[2] config \`${name}\` (${golden.configs[name].poseHashes.length} poses)...`);
    const g = golden.configs[name];
    const ff = new HeavyForceField({ atoms: sel.atoms }, { ...par }, []);
    ff.trackTerms = par.trackTerms === true;
    assert(ff.n === g.n && ff.nProt === g.nProt && ff.nLigAtoms === g.nLigAtoms &&
      ff.ligandStart === g.ligandStart,
    `${name}: atom bookkeeping identical (n=${ff.n} nProt=${ff.nProt} nLig=${ff.nLigAtoms} ligandStart=${ff.ligandStart})`);
    assert(ff.physicsLevel === g.physicsLevel && ff.chargesOn === g.chargesOn && ff.hbMode === g.hbMode,
      `${name}: physics-level mirror identical (${ff.physicsLevel}/${ff.chargesOn}/${ff.hbMode})`);
    assert(ff._excluded.size === g.excluded && ff._scale14.size === g.scale14 && ff.coord.length === g.coordLen,
      `${name}: exclusion/1-4/coord sets identical (excl=${ff._excluded.size} s14=${ff._scale14.size} coord=${ff.coord.length})`);
    assert(ff._weakRings.length === g.weakRings && ff._weakCations.length === g.weakCations &&
      ff._weakHalogens.length === g.weakHalogens,
    `${name}: weak lists identical (rings=${ff._weakRings.length} cations=${ff._weakCations.length} halogens=${ff._weakHalogens.length})`);
    const mEn = ff._metalEnforce;
    assert((mEn ? mEn.metals.length : 0) === g.metalEnforce && (mEn ? mEn.hasGeometry.size : 0) === g.metalGeom,
      `${name}: metal coordination upgrade identical (metals=${mEn ? mEn.metals.length : 0} geom=${mEn ? mEn.hasGeometry.size : 0})`);

    const ps = poses(ff.ref, g.poseHashes.length);
    let uBad = 0, fBad = 0, tBad = [];
    for (let p = 0; p < ps.length; p++) {
      const U = ff.compute(ps[p]);
      if (hashNum(U) !== g.poseHashes[p].U) uBad++;
      if (hashArr(ff.forces) !== g.poseHashes[p].forces) fBad++;
      for (const t of TERMS) {
        if (hashNum(ff[t]) !== g.poseHashes[p].terms[t]) tBad.push(`pose${p}.${t}`);
      }
    }
    assert(uBad === 0, `${name}: energy bit-identical at all ${ps.length} poses (${uBad} mismatch)`);
    assert(fBad === 0, `${name}: forces bit-identical at all ${ps.length} poses (${fBad} mismatch)`);
    assert(tBad.length === 0,
      `${name}: all ${TERMS.length} accumulators bit-identical at all poses` +
      (tBad.length ? ` (${tBad.length} mismatch: ${tBad.slice(0, 6).join(", ")})` : ""));
  }

  // ---- [3] trajectory ----------------------------------------------------
  console.log("\n[3] 10-step Langevin trajectory (seed 42)...");
  {
    const rng = new SeededRNG(42);
    rng.install();
    let h;
    try {
      const ff = new HeavyForceField({ atoms: sel.atoms }, CONFIGS.full, []);
      ff.trackTerms = true;
      const itg = new LangevinIntegrator(ff.ref, ff, 110);
      itg.setTemperature(300);
      itg.setFriction(5.0);
      for (let s = 0; s < golden.traj.steps; s++) itg.step();
      h = { energy: hashNum(ff.energy), forces: hashArr(ff.forces), pos: hashArr(itg.pos), vel: hashArr(itg.vel), dt: hashNum(itg.dt) };
    } finally { rng.restore(); }
    for (const k of ["pos", "vel", "energy", "forces", "dt"]) {
      assert(h[k] === golden.traj[k],
        `trajectory ${k} bit-identical (${h[k]}${h[k] === golden.traj[k] ? "" : ` != ${golden.traj[k]}`})`);
    }
  }

  finish();
}

function finish() {
  console.log(`\n=== test_heavy_golden: ${passed} PASSED, ${failed} FAILED ===`);
  if (failed > 0) process.exit(1);
  console.log("PASS: heavy-mode trajectory + kernel outputs bit-identical to the golden fixture");
}

main();