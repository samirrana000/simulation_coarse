/**
 * validate_flexlig.mjs — Stage-5 second validation system: 4W52 EPE flexible ligand.
 *
 * System choice + why (documented):
 *   - Chosen: 4W52 EPE pocket (same T4L protein, different ligand — cheapest).
 *     EPE (HEPES buffer, 15 atoms / 15 bonds) carries 4 rotatable bonds
 *     ([[0,9],[3,6],[6,7],[9,10]]; ring + S(=O)3 fan excluded) vs BNZ rigid 0.
 *     Same protein isolates the ligand-flexibility variable; crystal pose needs
 *     no placement; pocket differs honestly (EPE surface site 7 residues vs BNZ
 *     cavity 14) so this is a genuine second system, not a re-run.
 *   - Rejected: MOL2 library ligand with ≥2 rotatable bonds placed headlessly.
 *     Measured on src/ligandLib.js: benzene 0, phenol 0, toluene 0,
 *     chlorobenzene 0, indole 0, imidazole 0, acetate 0, ethanolamine 1,
 *     DMSO 0, caffeine 0 — max 1 (methyl/methoxy rotors are excluded by the
 *     degree>1 rule), so NO library ligand reaches ≥2. Placement would add
 *     moving parts for zero gain.
 *
 * Protocol (seeded, deterministic, fast <60s; CG-thermo class):
 *   EPE-only and BNZ-only CG systems (charges + directional-HB, same as the
 *   record path), 500 steps / stride 2 → 250 frames/leg, T 300 K, zeta 8.0,
 *   seeds 501/1501. Ligand torsion ΔS via the Stage-4 autoTorsions fallback
 *   (p.ligand graph + nProt offset — EPE-only layout, so holo/apo offsets
 *   agree). ΔSASA via thermo_sasa.sasaBurial (stride 10, CG bead radius).
 *   Locked control: synthetic single-dihedral ensemble (S≈0) vs sampled EPE.
 *
 * What it proves: the rotbond path (unit-tested in test_rotbonds.js) is live
 * end-to-end — sampled EPE ΔS_lig is NONZERO while rigid BNZ is exactly 0.
 *
 * Speed class: CG-thermo (≈ 4 legs × 500 CG steps + SASA, seconds). SLOW tier.
 * Run: node scripts/validate_flexlig.mjs
 */

/** Seeded replicas: holo SEED_HOLO, apo SEED_APO (distinct from thermo 101/1101). */
import { readFileSync } from "node:fs";
import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { computeThermodynamics, torsionEntropy } from "../src/analysis/thermodynamics.js";
import { findRotatableBonds } from "../src/analysis/rotbonds.js";
import { sasaBurial, cgBeadExtendedRadius, selectedLigandElements } from "../src/analysis/thermo_sasa.js";

const SEED_HOLO = 501;
const SEED_APO = 1501;
/** MD protocol: short CG relaxation (task: 200-500 steps). Units Å/ps/kcal/mol/Da. */
const STEPS = 500, STRIDE = 2, TEMP = 300, ZETA = 8.0, MASS = 110;
/** Pocket rule (Å) shared with the thermo path. */
const POCKET_RCUT = 8.0;
/** Real-burial subsample stride (matches THERMO_SASA_STRIDE). */
const SASA_STRIDE = 10;

let passed = 0, failed = 0;
/**
 * Counted check (feeds the SLOW-tier grand total via the results line).
 * @param {boolean} c condition
 * @param {string} m message
 */
const assert = (c, m) => { if (c) { passed++; console.log(`  ✓ ${m}`); } else { failed++; console.error(`  ✗ FAIL: ${m}`); } };

/**
 * Run one seeded Langevin leg, collecting frames + per-frame 7-term binding vector.
 * @param {object} ff live force field (trackTerms must be on for holo legs)
 * @param {number} seed RNG seed
 * @returns {{frames: Float32Array[], energies: number[][]}}
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
      energies.push([ff.bindLJU, ff.bindCoulU, ff.bindHBU, ff.desolvU, 0, 0, 0]);
    }
  }
  return { frames, energies };
}

/**
 * Residue indices within POCKET_RCUT of the ligand COM in the reference structure.
 * @param {object} ff force field with .ref/.nProt/.nLigAtoms
 * @param {object} sel selectSystem() output (for labels only)
 * @returns {{idx: number[], labels: string[]}}
 */
function pocket8A(ff, sel) {
  const nProt = ff.nProt;
  const com = [0, 0, 0];
  for (let a = 0; a < ff.nLigAtoms; a++) {
    com[0] += ff.ref[3 * (nProt + a)] / ff.nLigAtoms;
    com[1] += ff.ref[3 * (nProt + a) + 1] / ff.nLigAtoms;
    com[2] += ff.ref[3 * (nProt + a) + 2] / ff.nLigAtoms;
  }
  const idx = [];
  for (let i = 0; i < nProt; i++) {
    const d = Math.hypot(ff.ref[3 * i] - com[0], ff.ref[3 * i + 1] - com[1], ff.ref[3 * i + 2] - com[2]);
    if (d < POCKET_RCUT) idx.push(i);
  }
  return { idx, labels: idx.map((i) => `${sel.beads[i].resName}${sel.beads[i].resSeq}`) };
}

/**
 * Synthetic 4-atom frame with torsion a-b-c-d = phi (deg); same rig as
 * test_thermo/test_rotbonds (locked vs uniform rotor).
 * @param {number} phiDeg dihedral angle in degrees
 * @returns {Float32Array} 4-atom frame
 */
function frameWithDihedral(phiDeg) {
  const rad = (phiDeg * Math.PI) / 180;
  const d = [1 + Math.cos(rad), 1, Math.sin(rad)];
  return Float32Array.from([0, 0, 0, 1, 0, 0, 1, 1, 0, d[0], d[1], d[2]]);
}

const pdbText = readFileSync(new URL("../4w52.pdb", import.meta.url), "utf8");
const parsed = parseCa(pdbText);
const sel = selectSystem(parsed);
const allMols = parseLigands(pdbText);
const epeMols = allMols.filter((m) => m.resName === "EPE");
const bnzMols = allMols.filter((m) => m.resName === "BNZ");

// ---- rotatable-bond inventory (EPE flexible vs BNZ rigid) ----
console.log("=== flexlig rotatable inventory (4W52 EPE vs BNZ) ===");
const rEpe0 = findRotatableBonds(epeMols[0].atoms, epeMols[0].bonds);
const rBnz0 = findRotatableBonds(bnzMols[0].atoms, bnzMols[0].bonds);
console.log(`EPE: ${epeMols[0].atoms.length} atoms / ${epeMols[0].bonds.length} bonds → ${rEpe0.count} rotatable ${JSON.stringify(rEpe0.rotBonds)}`);
console.log(`BNZ: ${bnzMols[0].atoms.length} atoms / ${bnzMols[0].bonds.length} bonds → ${rBnz0.count} rotatable`);

// ---- EPE legs (flexible system) ----
console.log("\n=== EPE-only legs (seeded 501/1501, 500 steps, stride 2) ===");
const ffE = new ForceField(sel, { rc: 10, gamma: 2.0, binding: { charges: true, hbMode: "directional" } }, epeMols);
ffE.trackTerms = true;
const pocE = pocket8A(ffE, sel);
console.log(`pocket (${pocE.idx.length}): ${pocE.labels.join(" ")}`);
const holoE = runLeg(ffE, SEED_HOLO);
const ffApoE = new ForceField(sel, { rc: 10, gamma: 2.0, binding: { on: false } }, epeMols);
const apoE = runLeg(ffApoE, SEED_APO);
const resE = computeThermodynamics({
  holoFrames: holoE.frames, apoFrames: apoE.frames, holoEnergies: holoE.energies,
  pocketIdx: pocE.idx, nProt: ffE.nProt, mass: MASS, T: TEMP,
  ligand: { atoms: epeMols[0].atoms, bonds: epeMols[0].bonds, offset: ffE.nProt },
});

// ---- BNZ control (rigid system, same seeds/steps) ----
console.log("\n=== BNZ-only control (rigid, same seeds/steps) ===");
const ffB = new ForceField(sel, { rc: 10, gamma: 2.0, binding: { charges: true, hbMode: "directional" } }, bnzMols);
ffB.trackTerms = true;
const pocB = pocket8A(ffB, sel);
console.log(`pocket (${pocB.idx.length}): ${pocB.labels.join(" ")}`);
const holoB = runLeg(ffB, SEED_HOLO);
const ffApoB = new ForceField(sel, { rc: 10, gamma: 2.0, binding: { on: false } }, bnzMols);
const apoB = runLeg(ffApoB, SEED_APO);
const resB = computeThermodynamics({
  holoFrames: holoB.frames, apoFrames: apoB.frames, holoEnergies: holoB.energies,
  pocketIdx: pocB.idx, nProt: ffB.nProt, mass: MASS, T: TEMP,
  ligand: { atoms: bnzMols[0].atoms, bonds: bnzMols[0].bonds, offset: ffB.nProt },
});

// ---- locked-rotor control (synthetic, no MD) ----
const locked = Array.from({ length: 200 }, () => frameWithDihedral(30));
const sLocked = torsionEntropy(locked, [[0, 1, 2, 3]]);

// ---- real SASA burial on the EPE legs (cheap: stride 10, sub-second) ----
const sasaE = sasaBurial(holoE.frames, apoE.frames, {
  nProt: ffE.nProt,
  selElements: selectedLigandElements(epeMols, [0]),
  protRadius: cgBeadExtendedRadius(ffE),
  stride: SASA_STRIDE,
});

// ---- report ----
const f2 = (v) => (Number.isFinite(v) ? v.toFixed(2) : "NaN");
console.log("\n── flexlig anchor (EPE flexible vs BNZ rigid, kcal/mol) ──");
console.log(`EPE ΔH ${f2(resE.dH.total)} ± ${resE.dH_se.toFixed(2)}  (LJ ${f2(resE.dH.lj)} / Coul ${f2(resE.dH.coul)} / HB ${f2(resE.dH.hb)} / desolv ${f2(resE.dH.desolv)})`);
console.log(`EPE ΔS_lig ${resE.dS.ligand.toFixed(5)}  (${resE.meta.ligNote} [${resE.meta.rotatableBonds} rotatable]) vs locked ${sLocked.toFixed(5)} (≈0)`);
console.log(`BNZ ΔH ${f2(resB.dH.total)} ± ${resB.dH_se.toFixed(2)}  (LJ ${f2(resB.dH.lj)} / Coul ${f2(resB.dH.coul)} / HB ${f2(resB.dH.hb)} / desolv ${f2(resB.dH.desolv)})`);
console.log(`BNZ ΔS_lig ${resB.dS.ligand.toFixed(5)}  (${resB.meta.ligNote} [${resB.meta.rotatableBonds} rotatable]) — rigid zero control`);
console.log(`EPE ΔSASA ${sasaE.dsasa.toFixed(1)} ± ${sasaE.se.toFixed(1)} Å² (dLig ${sasaE.dLig.toFixed(1)} ± ${sasaE.dLigSE.toFixed(1)} + dProt ${sasaE.dProt.toFixed(1)} ± ${sasaE.dProtSE.toFixed(1)}; free-ligand ⟨S⟩ ${sasaE.ligFreeMean.toFixed(1)} Å²; ${sasaE.nHoloEval}+${sasaE.nApoEval} frames)`);

// ---- asserts (10; SLOW tier) ----
console.log("\n=== flexlig asserts ===");
assert(rEpe0.count === 4, `EPE 4 rotatable (got ${rEpe0.count}: ${JSON.stringify(rEpe0.rotBonds)})`);
assert(rEpe0.torsions.length === 4, `EPE one torsion per rotatable bond (got ${rEpe0.torsions.length})`);
assert(rBnz0.count === 0, `BNZ rigid → 0 rotatable (got ${rBnz0.count})`);
assert(sLocked < 1e-9, `locked rotor S≈0 (${sLocked.toFixed(5)})`);
assert(resE.dS.ligand > 0.0005, `EPE sampled ΔS_lig NONZERO end-to-end (${resE.dS.ligand.toFixed(5)} > 0.0005)`);
assert(resB.dS.ligand === 0 && resB.meta.rotatableBonds === 0, `BNZ rigid control ΔS_lig exactly 0 (got ${resB.dS.ligand})`);
assert(Number.isFinite(resE.dH.total), `EPE ΔH finite: ${f2(resE.dH.total)} (LJ ${f2(resE.dH.lj)} / Coul ${f2(resE.dH.coul)} / HB ${f2(resE.dH.hb)} / desolv ${f2(resE.dH.desolv)})`);
assert(pocE.idx.length === 7, `EPE pocket 7 residues @8Å (${pocE.labels.join(" ")})`);
assert(sasaE.dsasa > 100 && sasaE.dsasa < 700, `EPE real ΔSASA EPE-scale 100–700 Å²: ${sasaE.dsasa.toFixed(1)} ± ${sasaE.se.toFixed(1)} (dLig ${sasaE.dLig.toFixed(1)} + dProt ${sasaE.dProt.toFixed(1)})`);
assert(Number.isFinite(sasaE.se) && sasaE.se > 0 && sasaE.se < 30, `EPE solvent SE finite & small: ±${sasaE.se.toFixed(2)} (stride ${sasaE.stride}, ${sasaE.nHoloEval}+${sasaE.nApoEval} frames)`);

console.log(`\n=== validate_flexlig: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed ? 1 : 0);
