/**
 * test_rev1_issue5_rmsd_split.js — Revolution 1 / Issue 5 regression.
 *
 * Scope: src/heavy.js HeavyForceField.rmsd ONLY (no HUD edit needed).
 *
 * Problem: HeavyForceField.rmsd averaged over ALL atoms (protein + hetero +
 * ligand) while ForceField.rmsd is protein-only, making HUD RMSD incomparable
 * across CG/heavy and letting ligand drift masquerade as fold instability.
 *
 * Fix verified here (minimal, backward compatible):
 *   - rmsd(pos) restricted to the protein-only slice i < nProt*3
 *     (mirrors ForceField.rmsd: sqrt(sum/nProt)).
 *   - rmsdLig(pos) covers the ligand block [ligandStart, n) (same slice the
 *     HUD ligRMSD uses); rmsdAll(pos) preserves the legacy all-atom value.
 *   - HUD call-site needs no edit: main.js already renders protein RMSD via
 *     ff.rmsd(pos) plus a separate ligRMSD block (ligStart..n / nLigAtoms).
 *
 * Run: node tests/test_rev1_issue5_rmsd_split.js (fast, <2s; NOT wired into
 * tests/test_all.js FAST so the 352 gate is untouched).
 */

import { HeavyForceField } from "../src/heavy.js";
import { ForceField } from "../src/forcefield.js";

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

function mkAtom(x, y, z, opts = {}) {
  return {
    element: "C",
    x, y, z,
    atomName: opts.atomName ?? "CA",
    resName: opts.resName ?? "ALA",
    chain: opts.chain ?? "A",
    resSeq: opts.resSeq ?? 1,
    serial: opts.serial ?? 0,
    isProtein: opts.isProtein ?? true,
    isMetal: opts.isMetal ?? false,
    isLigand: opts.isLigand ?? false,
    isHetero: opts.isHetero ?? false,
    heteroKey: opts.heteroKey ?? null,
  };
}

function main() {
  console.log("=== Rev1/Issue5: heavy RMSD protein-only split ===");

  // Static protein (4 atoms) + 1 hetero + 2 ligand atoms, ordered
  // protein -> hetero -> ligand so ligandStart = nProt + nHetero.
  const atoms = [
    mkAtom(0, 0, 0, { resSeq: 1 }),
    mkAtom(1.5, 0, 0, { resSeq: 2 }),
    mkAtom(0, 1.5, 0, { resSeq: 3 }),
    mkAtom(0, 0, 1.5, { resSeq: 4 }),
    mkAtom(10, 10, 10, { isProtein: false, isHetero: true, resName: "ZN", atomName: "ZN" }),
    mkAtom(20, 0, 0, { isProtein: false, isLigand: true, resName: "LIG", atomName: "C1" }),
    mkAtom(21.5, 0, 0, { isProtein: false, isLigand: true, resName: "LIG", atomName: "C2" }),
  ];
  const ff = new HeavyForceField({ atoms }, { gamma: 1.0 }, []);
  assert(ff.n === 7, `n = 7 (got ${ff.n})`);
  assert(ff.nProt === 4, `nProt = 4 (got ${ff.nProt})`);
  assert(ff.nLigAtoms === 2, `nLigAtoms = 2 (got ${ff.nLigAtoms})`);
  assert(ff.ligandStart === 5, `ligandStart = 5 (got ${ff.ligandStart})`);

  // Displace the ligand block by +5 A in x; protein + hetero stay put.
  const pos = new Float64Array(ff.ref);
  for (let a = ff.ligandStart; a < ff.n; a++) pos[3 * a] += 5.0;

  const prot = ff.rmsd(pos);
  const lig = ff.rmsdLig(pos);
  const all = ff.rmsdAll(pos);
  assert(prot === 0, `protein rmsd = 0 with static protein (got ${prot})`);
  assert(Math.abs(lig - 5.0) < 1e-12, `ligand rmsd = 5 A (got ${lig})`);
  assert(all > 0 && all < lig, `legacy all-atom rmsd diluted by protein (got ${all.toFixed(4)} A)`);
  // Legacy value must equal the old formula sqrt(sum_all / n).
  let sAll = 0;
  for (let i = 0; i < pos.length; i++) {
    const d = pos[i] - ff.ref[i];
    sAll += d * d;
  }
  assert(Math.abs(all - Math.sqrt(sAll / ff.n)) < 1e-12, `rmsdAll matches legacy formula`);
  // Hetero static => excluded from both splits.
  assert(prot === 0 && Math.abs(lig - 5.0) < 1e-12, `hetero excluded from protein + ligand splits`);

  // Protein motion still registers 1:1 in rmsd().
  const pos2 = new Float64Array(ff.ref);
  for (let a = 0; a < ff.nProt; a++) pos2[3 * a + 1] += 2.0;
  assert(Math.abs(ff.rmsd(pos2) - 2.0) < 1e-12, `protein 2 A shift => rmsd 2 A (got ${ff.rmsd(pos2)})`);
  assert(ff.rmsdLig(pos2) === 0, `ligand rmsd 0 when ligand static (got ${ff.rmsdLig(pos2)})`);

  // No-ligand system: rmsdLig returns 0 (no NaN), rmsd still works.
  const protOnly = [
    mkAtom(0, 0, 0, { resSeq: 1 }),
    mkAtom(1.5, 0, 0, { resSeq: 2 }),
  ];
  const ffP = new HeavyForceField({ atoms: protOnly }, { gamma: 1.0 }, []);
  assert(ffP.rmsdLig(ffP.ref) === 0, `no-ligand rmsdLig = 0`);
  assert(ffP.rmsd(ffP.ref) === 0, `no-ligand rmsd at ref = 0`);

  // Parity spot-check: ForceField.rmsd is protein-only by construction.
  const selCG = { beads: [{ x: 0, y: 0, z: 0 }, { x: 3, y: 0, z: 0 }], segments: [] };
  const ffCG = new ForceField(selCG, { rc: 10, gamma: 1.0 }, []);
  const cgPos = new Float64Array(ffCG.ref);
  assert(ffCG.rmsd(cgPos) === 0, `CG protein-only rmsd at ref = 0`);

  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log("=================================================");
  if (failed > 0) process.exit(1);
  console.log("PASS: test_rev1_issue5_rmsd_split.js — protein-only RMSD split validated");
}

main();
