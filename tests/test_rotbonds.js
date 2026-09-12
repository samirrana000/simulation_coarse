/**
 * tests/test_rotbonds.js — Stage-4 rotatable-bond detection + flexible-ligand ΔS.
 *
 * Covers src/analysis/rotbonds.js (findRotatableBonds / autoTorsions) and the
 * autoTorsions fallback wired into computeThermodynamics ligand-ΔS path:
 *  - rigid ligands (benzene MOL2, 4W52 BNZ: all-ring bonds) → 0 rotatable;
 *  - flexible ligand (4W52 HEPES/EPE, 15 bonds) → 4 rotatable
 *    (piperazine tails + ethanesulfonate chain; ring + S(=O)₃ head excluded);
 *  - rule unit checks: 4-atom chain (1), amide C–N excluded (0),
 *    tert-butyl-like fan excluded (0);
 *  - flexible-ligand physics: torsionEntropy > 0 on a rotating ensemble,
 *    0 when locked; computeThermodynamics auto fallback yields ΔS_lig > 0
 *    with rotatable count reported in meta + formatThermoTable ligand line;
 *    rigid/legacy paths stay exactly 0 (bit-identical).
 *
 * Runnable: node tests/test_rotbonds.js
 * Prints PASS/FAIL per assertion; exits 1 on any failure.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { findRotatableBonds, autoTorsions } from "../src/analysis/rotbonds.js";
import { computeThermodynamics, torsionEntropy, formatThermoTable } from "../src/analysis/thermodynamics.js";
import { parseLigands } from "../src/pdb.js";
import { parseMol2 } from "../src/mol2.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

let passed = 0, failed = 0;
function assert(c, m) {
  if (c) { passed++; console.log(`  ✓ ${m}`); }
  else { failed++; console.error(`  ✗ FAIL: ${m}`); }
}

function findFile(name) {
  const candidates = [
    path.resolve(process.cwd(), name),
    path.resolve(__dirname, "..", name),
    path.resolve(__dirname, name),
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`${name} not found (tried ${candidates.join(", ")})`);
}

/** 4-atom frame with torsion a-b-c-d = phi (deg), same rig as test_thermo. */
function frameWithDihedral(phi) {
  const rad = (phi * Math.PI) / 180;
  const d = [1 + Math.cos(rad), 1, Math.sin(rad)];
  return Float32Array.from([0, 0, 0, 1, 0, 0, 1, 1, 0, d[0], d[1], d[2]]);
}

function main() {
  console.log("=== test_rotbonds (Stage-4 rotatable bonds + flexible-ligand ΔS) ===");

  // ---- (1) rigid: benzene MOL2 → 0 ----
  const mol2Text = fs.readFileSync(findFile("benzene.mol2"), "utf-8");
  const bmol = parseMol2(mol2Text)[0];
  const rBenz = findRotatableBonds(bmol.atoms, bmol.bonds);
  assert(rBenz.count === 0, `benzene.mol2 rigid → 0 rotatable (got ${rBenz.count})`);

  // ---- (2) rigid: 4W52 BNZ (HETATM, all-ring) → 0 ----
  const pdbText = fs.readFileSync(findFile("4w52.pdb"), "utf-8");
  const ligs = parseLigands(pdbText);
  const bnz = ligs.find((m) => m.resName === "BNZ");
  const rBnz = findRotatableBonds(bnz.atoms, bnz.bonds);
  assert(rBnz.count === 0, `4W52 BNZ rigid → 0 rotatable (got ${rBnz.count})`);

  // ---- (3) flexible: 4W52 EPE (HEPES-like, 15 bonds) → 4 ----
  const epe = ligs.find((m) => m.resName === "EPE");
  const rEpe = findRotatableBonds(epe.atoms, epe.bonds);
  assert(rEpe.count === 4, `4W52 EPE flexible → 4 rotatable (got ${rEpe.count}: ${JSON.stringify(rEpe.rotBonds)})`);
  assert(rEpe.torsions.length === 4, `EPE emits one torsion per rotatable bond (got ${rEpe.torsions.length})`);

  // ---- (4) synthetic 4-atom chain → 1 rotatable, torsion [0,1,2,3] ----
  const chain = [
    { x: 0, y: 0, z: 0, element: "C" },
    { x: 1.5, y: 0, z: 0, element: "C" },
    { x: 3.0, y: 0, z: 0, element: "C" },
    { x: 4.5, y: 0, z: 0, element: "C" },
  ];
  const cbonds = [[0, 1], [1, 2], [2, 3]];
  const rChain = findRotatableBonds(chain, cbonds);
  assert(rChain.count === 1, `4-atom chain → 1 rotatable (got ${rChain.count})`);
  assert(JSON.stringify(rChain.torsions) === "[[0,1,2,3]]",
    `chain torsion quadruplet [0,1,2,3] (got ${JSON.stringify(rChain.torsions)})`);

  // ---- (5) amide C–N excluded (N-methylacetamide, C=O 1.23 Å) ----
  const amide = [
    { x: 0, y: 0, z: 0, element: "C" },
    { x: 1.5, y: 0, z: 0, element: "C" },
    { x: 2.0, y: 1.2, z: 0, element: "O" },
    { x: 2.2, y: -1.2, z: 0, element: "N" },
    { x: 3.7, y: -1.2, z: 0, element: "C" },
  ];
  const rAm = findRotatableBonds(amide, [[0, 1], [1, 2], [1, 3], [3, 4]]);
  assert(rAm.count === 0, `amide C–N excluded → 0 rotatable (got ${rAm.count})`);

  // ---- (6) tert-butyl-like fan excluded ----
  const tbut = [
    { x: 0, y: 0, z: 0, element: "C" },
    { x: 1.5, y: 0, z: 0, element: "C" },
    { x: -0.5, y: 1.4, z: 0, element: "C" },
    { x: -0.5, y: -1.4, z: 0, element: "C" },
    { x: 3.0, y: 0, z: 0, element: "C" },
    { x: 4.5, y: 0, z: 0, element: "C" },
  ];
  const rTb = findRotatableBonds(tbut, [[0, 1], [0, 2], [0, 3], [0, 4], [4, 5]]);
  assert(rTb.count === 0, `tert-butyl-like fan excluded → 0 rotatable (got ${rTb.count})`);

  // ---- (7) flexible-ligand physics: rotating S > 0, locked S = 0 ----
  const locked = Array.from({ length: 500 }, () => frameWithDihedral(30));
  const flex = [];
  for (let i = 0; i < 600; i++) flex.push(frameWithDihedral((i % 12) * 30 + 15));
  const sLocked = torsionEntropy(locked, [[0, 1, 2, 3]]);
  const sFlex = torsionEntropy(flex, [[0, 1, 2, 3]]);
  assert(sLocked < 1e-9, `locked rotor S≈0 (${sLocked.toFixed(5)})`);
  assert(sFlex > 0, `flexible rotor S>0 on sampled ensemble (${sFlex.toFixed(5)})`);

  // ---- (8) autoTorsions offset ----
  const off = autoTorsions(chain, cbonds, 10);
  assert(JSON.stringify(off) === "[[10,11,12,13]]",
    `autoTorsions offset 10 → [[10,11,12,13]] (got ${JSON.stringify(off)})`);

  // ---- (9) computeThermodynamics auto fallback: flexible → ΔS_lig > 0 ----
  const res = computeThermodynamics({
    holoFrames: flex, apoFrames: locked,
    ligand: { atoms: chain, bonds: cbonds },
    pocketIdx: [],
  });
  assert(res.dS.ligand > 0, `auto flexible-ligand ΔS_lig > 0 (${res.dS.ligand.toFixed(5)})`);
  assert(res.meta.rotatableBonds === 1, `meta.rotatableBonds = 1 (got ${res.meta.rotatableBonds})`);
  assert(/rotatable/.test(res.meta.ligNote), `ligNote reports rotatable count ("${res.meta.ligNote}")`);

  // ---- (10) rigid + legacy paths stay exactly 0 ----
  const resRigid = computeThermodynamics({
    holoFrames: flex, apoFrames: locked,
    ligand: { atoms: bmol.atoms, bonds: bmol.bonds },
    pocketIdx: [],
  });
  assert(resRigid.dS.ligand === 0 && resRigid.meta.rotatableBonds === 0,
    `rigid benzene auto → ΔS_lig 0, 0 rotatable (got ${resRigid.dS.ligand})`);
  const resLegacy = computeThermodynamics({ holoFrames: flex, apoFrames: locked, pocketIdx: [] });
  assert(resLegacy.dS.ligand === 0,
    `legacy path (no ligand graph) → ΔS_lig exactly 0 (got ${resLegacy.dS.ligand})`);

  // ---- (11) formatThermoTable ligand line carries the count ----
  const table = formatThermoTable(res);
  const ligLine = table.split("\n").find((l) => l.includes("ligand"));
  assert(/\[1 rotatable\]/.test(ligLine), `thermo table ligand line reports count ("${ligLine.trim()}")`);

  console.log(`\n=== test_rotbonds: ${passed} PASSED, ${failed} FAILED ===`);
  process.exit(failed ? 1 : 0);
}

main();
