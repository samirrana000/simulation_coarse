/**
 * validate_1crn_null.mjs — FP6 negative control: 1CRN apo-vs-apo null.
 *
 * Why 1CRN (documented choice): crambin (1CRN, 46 residues) is ligand-free
 * (0 HETATM ligands parsed), so the holo−apo thermo pipeline has NO binding
 * signal to find. Both legs are apo CG (binding off, no ligands, same record
 * protocol as scripts/test_thermo.mjs: T 300 K, zeta 8.0, mass 110 Da,
 * charges+directional-HB irrelevant with no ligand). Leg A plays the
 * holo role, leg B the apo role; binding-energy vectors are zeros by
 * construction, so ΔH ≈ 0 is the null expectation and the Schlitter
 * pocket difference measures sampling noise only (never a forced zero —
 * every assert is finiteness/boundedness, sign never asserted).
 *
 * Pocket rule: residues within 8 Å of the PROTEIN COM in ref (no ligand COM
 * exists). Live: 16 residues @8Å (THR2…ILE34 patch).
 *
 * Protocol (seeded, deterministic, fast <60s; CG-thermo class):
 *   STEPS 500 / STRIDE 2 → 250 frames/leg; SEED_HOLO 701 / SEED_APO 1701
 *   (`NULL_SEED_BASE` env override mirrors THERMO_SEED_BASE).
 *
 * Expectation (measured 2026-09-13, seeds 701/1701, twice bit-identical):
 *   ΔH 0.00 ± 0.00 (null-scale bound |ΔH| < 0.50 vs 4W52 signal −3.6/−6.8);
 *   −TΔS_pocket +5.67 (bounded < 20, same bar as test_thermo — sampling
 *   noise between two apo ensembles, honestly nonzero);
 *   f/DOF 5.2 UNDER-SAMPLED by design (pipeline null-smoke only).
 *
 * Speed class: CG-thermo (≈ 2 legs × 500 CG steps, <1 s). SLOW tier.
 * Run: node scripts/validate_1crn_null.mjs
 */

import { readFileSync } from "node:fs";
import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { computeThermodynamics } from "../src/analysis/thermodynamics.js";

/** Seeded legs: leg-A SEED_A, leg-B SEED_B (distinct from thermo 101/1101). */
const SEED_BASE = Number.parseInt(process.env.NULL_SEED_BASE || "0", 10) || 0;
const SEED_A = 701 + SEED_BASE;
const SEED_B = 1701 + SEED_BASE;
/** MD protocol: short CG relaxation (flexlig class). Units Å/ps/kcal/mol/Da. */
const STEPS = 500, STRIDE = 2, TEMP = 300, ZETA = 8.0, MASS = 110;
/** Pocket rule (Å) shared with the thermo path (here: of the protein COM). */
const POCKET_RCUT = 8.0;
/** Null-scale bounds (documented): |ΔH| vs 4W52 signal; |−TΔS| = test_thermo bar. */
const NULL_DH_BOUND = 0.5, NULL_TDS_BOUND = 20;

let passed = 0, failed = 0;
/**
 * Counted check (feeds the SLOW-tier grand total via the results line).
 * @param {boolean} c condition
 * @param {string} m message
 */
const assert = (c, m) => { if (c) { passed++; console.log(`  ✓ ${m}`); } else { failed++; console.error(`  ✗ FAIL: ${m}`); } };

/**
 * Run one seeded apo leg (binding off, no ligands — pure protein diffusion).
 * @param {object} ff live force field (binding off)
 * @param {number} seed RNG seed
 * @returns {{frames: Float32Array[], energies: number[][]}} frames + zero binding vectors
 */
function runLeg(ff, seed) {
  const integ = new LangevinIntegrator(ff.ref, ff, MASS, { seed });
  integ.setTemperature(TEMP);
  integ.setFriction(ZETA);
  const frames = [], energies = [];
  for (let s = 0; s < STEPS; s++) {
    integ.step();
    if (s % STRIDE === 0) {
      frames.push(Float32Array.from(integ.pos));
      energies.push([0, 0, 0, 0, 0, 0, 0]); // no ligand ⇒ no binding energy, by construction
    }
  }
  return { frames, energies };
}

/**
 * Residue indices within POCKET_RCUT of the protein COM in ref (no ligand COM).
 * @param {object} ff force field with .ref/.nProt
 * @param {object} sel selectSystem() output (for labels only)
 * @returns {{idx: number[], labels: string[]}}
 */
function pocket8A(ff, sel) {
  const nProt = ff.nProt;
  const com = [0, 0, 0];
  for (let i = 0; i < nProt; i++) {
    com[0] += ff.ref[3 * i] / nProt;
    com[1] += ff.ref[3 * i + 1] / nProt;
    com[2] += ff.ref[3 * i + 2] / nProt;
  }
  const idx = [];
  for (let i = 0; i < nProt; i++) {
    const d = Math.hypot(ff.ref[3 * i] - com[0], ff.ref[3 * i + 1] - com[1], ff.ref[3 * i + 2] - com[2]);
    if (d < POCKET_RCUT) idx.push(i);
  }
  return { idx, labels: idx.map((i) => `${sel.beads[i].resName}${sel.beads[i].resSeq}`) };
}

const pdbText = readFileSync(new URL("../1crn.pdb", import.meta.url), "utf8");
const parsed = parseCa(pdbText);
const sel = selectSystem(parsed);
const mols = parseLigands(pdbText);

// ---- null legs (apo-vs-apo, seeded) ----
console.log("=== 1CRN null (apo-vs-apo, seeded 701/1701, 500 steps, stride 2) ===");
console.log(`1CRN: ${sel.beads.length} Cα, ${mols.length} ligands (ligand-free control)`);
const ffA = new ForceField(sel, { rc: 10, gamma: 2.0, binding: { on: false } }, []);
const poc = pocket8A(ffA, sel);
console.log(`pocket (${poc.idx.length}): ${poc.labels.join(" ")}`);
const legA = runLeg(ffA, SEED_A);
const ffB = new ForceField(sel, { rc: 10, gamma: 2.0, binding: { on: false } }, []);
const legB = runLeg(ffB, SEED_B);
const res = computeThermodynamics({
  holoFrames: legA.frames, apoFrames: legB.frames, holoEnergies: legA.energies,
  pocketIdx: poc.idx, nProt: ffA.nProt, mass: MASS, T: TEMP,
});
const minusTds = -TEMP * res.dS.pocket;

// ---- report ----
const f2 = (v) => (Number.isFinite(v) ? v.toFixed(2) : "NaN");
console.log("\n── 1CRN null (apo-vs-apo, kcal/mol) ──");
console.log(`ΔH ${f2(res.dH.total)} ± ${res.dH_se.toFixed(2)} (no ligand ⇒ no binding energy, by construction)`);
console.log(`−TΔS_pocket ${f2(minusTds)} (sampling noise only; sign never asserted)`);
console.log(`S_pocket holo ${res.S_pocket.holo.toFixed(3)} vs apo ${res.S_pocket.apo.toFixed(3)}; f/DOF ${res.meta.framesPerDof.toFixed(1)}; frames ${res.meta.holoFrames}+${res.meta.apoFrames}`);

// ---- asserts (8; SLOW tier; boundedness only, never forced zero) ----
console.log("\n=== 1crn-null asserts ===");
assert(mols.length === 0, `1CRN ligand-free: ${mols.length} HETATM ligands (null control valid)`);
assert(sel.beads.length === 46, `1CRN 46 Cα (got ${sel.beads.length})`);
assert(poc.idx.length >= 5 && poc.idx.length <= 30, `pocket COM-8Å count sane: ${poc.idx.length} (16 live)`);
assert(Number.isFinite(res.dH.total), `null ΔH finite: ${f2(res.dH.total)}`);
assert(Math.abs(res.dH.total) < NULL_DH_BOUND, `null |ΔH| < ${NULL_DH_BOUND} (got ${f2(res.dH.total)}; 4W52 signal −3.6/−6.8)`);
assert(Number.isFinite(res.dH_se), `bootstrap SE finite: ±${res.dH_se.toFixed(2)}`);
assert(Number.isFinite(minusTds) && Math.abs(minusTds) < NULL_TDS_BOUND, `|−TΔS_pocket| bounded < ${NULL_TDS_BOUND}: ${f2(minusTds)} (test_thermo bar; noise only)`);
assert(Number.isFinite(res.S_pocket.holo) && Number.isFinite(res.S_pocket.apo) && res.meta.holoFrames === STEPS / STRIDE, `S_pocket finite (holo ${res.S_pocket.holo.toFixed(3)} / apo ${res.S_pocket.apo.toFixed(3)}) + ${STEPS / STRIDE} frames/leg`);

console.log(`\n=== validate_1crn_null: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed ? 1 : 0);
