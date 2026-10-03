/**
 * tests/generate_heavy_golden.js — write tests/golden/heavy_4w52_fp.json.
 *
 * Run: node tests/generate_heavy_golden.js
 *
 * WHY A SECOND GOLDEN
 * -------------------
 * tests/test_golden.js pins the Cα CG trajectory (parseCa + ForceField +
 * LangevinIntegrator). It imports NOTHING from src/heavy.js, so it is a
 * regression net for the CG engine only. The heavy engine has finite-difference
 * tests (test_forces_fd.js, test_gb_fd.js) — good for gradient CORRECTNESS,
 * but they assert a tolerance (1e-3 rel) on 10 atoms, not bit-identity, and a
 * 1e-3 tolerance cannot see a reordering that changes the last bits of a
 * force sum. This fixture pins heavy mode bit-exactly, so a refactor that
 * reorders a summation is caught at the bit.
 *
 * WHAT IS PINNED (all headless, no DOM)
 *   • topology: bonds/angles/propers/impropers/coord element-for-element
 *   • 6 constructor configurations, each of which switches on a DIFFERENT
 *     kernel pair, so together they cover every branch of compute():
 *       legacy     harmonicFlat / angleFlat / _nonBondedGrid (HCT) / SasaModel
 *       amber14    harmonicFlatPerK / angleFlatPerK
 *       obc2       _nonBondedGridOBC2 + _nonBondedGridNoGB
 *       lcpo       lcpoSasa
 *       full       weak (pi/cpi/xb) + metalAngles + membrane + per-term trackers
 *   • 12 deterministic poses (native + SeededRNG(1234) jitter), so the kernels
 *     are compared away from the native pose too — the FD tests only look at
 *     the native one
 *   • a 10-step Langevin trajectory with SeededRNG(42) — the heavy analogue of
 *     tests/test_golden.js
 *
 * HASHING: 4 independent 32-bit FNV-1a lanes over the RAW IEEE-754 bytes, so
 * equality here is bit equality, not tolerance equality.
 *
 * REGENERATED 2026-10-03, AFTER A DELIBERATE PHYSICS FIX (not a refactor)
 *   The `full` row used to pin a BUG on purpose: with par.gbModel === "obc2"
 *   AND ff.trackTerms === true, _nonBondedGridNoGB evaluated `ligStart` before
 *   its `const` declaration (TDZ) and threw; compute() caught it and fell back
 *   to _nonBondedGrid — but the OBC2 forces had already been merged into `f`
 *   before the throw, so forces carried OBC2-GB + HCT-everything while the
 *   energy was pure HCT. Fixed in src/heavy/nonbonded.js + src/heavy/energy.js
 *   (transactional OBC2 kernel, snapshot/restore around every substitutable
 *   stage, loud ff.physicsFallbacks). This fixture was then regenerated.
 *
 *   Only `obc2` forces, `full` (energy + forces + elecU + gbU) and the
 *   trajectory differ; `legacy`, `amber14`, `lcpo` and the whole `topo` block
 *   are byte-identical. No tolerance exists in this fixture to loosen — it is
 *   bit-exact by construction. The independent correctness proof for the new
 *   numbers (per-term tracking ON vs OFF agree bit-for-bit) is in
 *   tests/test_heavy_obc2_trackterms.js.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { parseHeavy, selectHeavy, appendHeavyLigands, HeavyForceField } from "../src/heavy.js";
import { parseMol2 } from "../src/mol2.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { SeededRNG } from "../src/seeded-rng.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");

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
function hashF64(arr) {
  _f64[0] = arr;
  return hashBytes(_u8);
}
const hashArr = (a) => hashBytes(new Uint8Array(a.buffer ?? a, (a.byteOffset ?? 0), a.byteLength ?? a.length * 8));
/** number -> its exact IEEE-754 bit pattern hash, so 0 vs -0 is distinguishable. */
function hashNum(x) { return hashF64(x); }

/** Accumulator fields pinned per configuration, in a FIXED order. */
const TERMS = [
  "energy", "bindingU", "desolvU", "repU", "bondU", "angleU", "improperU",
  "properU", "coordU", "elecU", "gbU", "sasaU", "hbondU", "membraneU",
  "springU", "weakU", "piU", "cpiU", "xbU", "coordAngleU",
  "bindLJU", "bindCoulU", "bindHBU", "_nanStrikes",
];

export const CONFIGS = {
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

/** 12 poses: native + 11 SeededRNG(1234) jitters of +/- 0.05 A. */
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

async function main() {
  const sel = buildSel();
  const out = { _generated: "tests/generate_heavy_golden.js", configs: {}, poses: {}, traj: {}, topo: {} };

  for (const [name, par] of Object.entries(CONFIGS)) {
    const ff = new HeavyForceField({ atoms: sel.atoms }, { ...par }, []);
    ff.trackTerms = par.trackTerms === true;
    const ps = poses(ff.ref);
    const rows = [];
    for (const pos of ps) {
      const U = ff.compute(pos);
      const terms = {};
      for (const t of TERMS) terms[t] = hashNum(ff[t]);
      rows.push({ U: hashNum(U), forces: hashArr(ff.forces), terms });
    }
    out.configs[name] = {
      n: ff.n, nProt: ff.nProt, nLigAtoms: ff.nLigAtoms, nHetero: ff.nHetero,
      ligandStart: ff.ligandStart,
      physicsLevel: ff.physicsLevel, chargesOn: ff.chargesOn, hbMode: ff.hbMode,
      weakRings: ff._weakRings.length, weakCations: ff._weakCations.length,
      weakHalogens: ff._weakHalogens.length,
      metalEnforce: ff._metalEnforce ? ff._metalEnforce.metals.length : 0,
      metalGeom: ff._metalEnforce ? ff._metalEnforce.hasGeometry.size : 0,
      excluded: ff._excluded.size, scale14: ff._scale14.size,
      coordLen: ff.coord.length,
      poseHashes: rows,
    };
  }

  // Topology is configuration-independent; pin it once, element for element.
  {
    const ff = new HeavyForceField({ atoms: sel.atoms }, CONFIGS.full, []);
    out.topo = {
      bonds: hashArr(ff.bonds), angles: hashArr(ff.angles),
      propers: hashArr(ff.propers), impropers: hashArr(ff.impropers),
      coord: hashArr(ff.coord), ligandBonds: hashArr(ff.ligandBonds),
      ref: hashArr(ff.ref), masses: hashArr(ff.masses),
      counts: { bonds: ff.bonds.length, angles: ff.angles.length, propers: ff.propers.length, impropers: ff.impropers.length },
    };
  }

  // 10-step heavy Langevin, seed 42 (the heavy analogue of test_golden.js).
  {
    const rng = new SeededRNG(42);
    rng.install();
    try {
      const ff = new HeavyForceField({ atoms: sel.atoms }, CONFIGS.full, []);
      ff.trackTerms = true;
      const itg = new LangevinIntegrator(ff.ref, ff, 110);
      itg.setTemperature(300);
      itg.setFriction(5.0);
      for (let s = 0; s < 10; s++) itg.step();
      out.traj = {
        seed: 42, steps: 10,
        energy: hashNum(ff.energy),
        forces: hashArr(ff.forces),
        pos: hashArr(itg.pos),
        vel: hashArr(itg.vel),
        dt: hashNum(itg.dt),
      };
    } finally { rng.restore(); }
  }

  const p = path.join(ROOT, "tests", "golden", "heavy_4w52_fp.json");
  fs.writeFileSync(p, JSON.stringify(out, null, 1) + "\n");
  console.log(`wrote ${p}`);
  for (const [name, c] of Object.entries(out.configs)) {
    console.log(`  ${name.padEnd(8)} n=${c.n} excluded=${c.excluded} scale14=${c.scale14} ` +
      `pose0.U=${c.poseHashes[0].U} pose0.F=${c.poseHashes[0].forces}`);
  }
  console.log(`  traj     pos=${out.traj.pos} vel=${out.traj.vel} E=${out.traj.energy}`);
}

main().catch((e) => { console.error(e); process.exit(1); });