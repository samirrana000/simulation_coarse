/**
 * tests/test_ala_noise_floor.js — Stage-4 alanine-scan honest display (Option B).
 *
 * CG Cα-ENM ΔΔG cancels in the holo−apo cycle (live 4W52 BNZ magnitudes all
 * |ΔΔG| < 0.02), so the scan is demoted to ranking-only: rows below
 * ALA_DDG_NOISE_FLOOR (0.05 kcal/mol) are flagged "~noise" by
 * formatMutationTable with a ranking-only disclaimer, and scanPocket tags
 * rows ({noise, noiseFloor}) deterministically. No force-field retuning —
 * display only. Option A (sidechain-count scaling) was measured and rejected:
 * same top-3/order, arbitrary inflation (TYR88 0.014→0.112) with no new
 * discriminating power — see docs/BINDING_LOOP2_DONE.md §18.
 *
 * Covers src/analysis/alanine_scanning.js (ALA_DDG_NOISE_FLOOR /
 * isNoiseDdG / annotateScanNoise / scanPocket tagging / formatMutationTable
 * flags) on synthetic rows + the live 4W52 BNZ cavity (rCut 8 Å, maxN 12,
 * relax 80 — calibration protocol): ordering bit-stable across two runs,
 * magnitudes honestly below floor, top-3 still contains MET102.
 *
 * Runnable: node tests/test_ala_noise_floor.js
 * Prints PASS/FAIL per assertion; exits 1 on any failure.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import {
  ALA_DDG_NOISE_FLOOR, isNoiseDdG, annotateScanNoise,
  pocketResidues, scanPocket, formatMutationTable,
} from "../src/analysis/alanine_scanning.js";

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

// ---- pure-unit: floor + flagging (no FF, deterministic) ----
assert(ALA_DDG_NOISE_FLOOR === 0.05, `noise floor is 0.05 kcal/mol (got ${ALA_DDG_NOISE_FLOOR})`);
assert(isNoiseDdG(0.014) === true, "isNoiseDdG(0.014) → noise (TYR88 live scale)");
assert(isNoiseDdG(0.06) === false, "isNoiseDdG(0.06) → signal (above floor)");
assert(isNoiseDdG(-0.049) === true, "isNoiseDdG(−0.049) → noise (just below floor)");
assert(isNoiseDdG(NaN) === true, "isNoiseDdG(NaN) → noise (never claim signal on NaN)");
assert(isNoiseDdG(0.04, 0.01) === false, "custom floor honored: isNoiseDdG(0.04, 0.01) → signal");

const synth = [{ label: "A", ddG: 0.01, note: "" }, { label: "B", ddG: 0.2, note: "" }];
annotateScanNoise(synth);
assert(synth[0].noise === true && synth[1].noise === false, "annotateScanNoise tags {noise} per row (true/false)");
const synthTable = formatMutationTable(synth);
assert(synthTable.includes("~noise") && synthTable.includes("ranking only"), "formatMutationTable flags ~noise + ranking-only disclaimer");

// ---- live 4W52 BNZ cavity (calibration protocol, deterministic: no RNG) ----
const pdbText = fs.readFileSync(findFile("4w52.pdb"), "utf-8");
const sel = selectSystem(parseCa(pdbText));
assert(sel.beads.length === 164, `4W52 Cα beads = 164 (got ${sel.beads.length})`);
const bnzMols = parseLigands(pdbText).filter((m) => m.resName === "BNZ");
const ff = new ForceField(sel, { rc: 10, gamma: 2.0, binding: { charges: true, hbMode: "directional" } }, bnzMols);
const sys = { mode: "cg", sel, ff, ligands: bnzMols };
const pocket = pocketResidues(sys, { rCut: 8.0, maxN: 12 });
assert(pocket.length === 12, `BNZ pocket n = 12 (got ${pocket.length})`);
const scan1 = scanPocket(sys, pocket.map((p) => p.resId), { relaxSteps: 80 });
const scan2 = scanPocket(sys, pocket.map((p) => p.resId), { relaxSteps: 80 });
assert(scan1.rows.every((r) => Number.isFinite(r.ddG)), `scan table finite ×${scan1.n}`);
const maxAbs = Math.max(...scan1.rows.map((r) => Math.abs(r.ddG)));
assert(maxAbs < ALA_DDG_NOISE_FLOOR, `honest magnitudes: max|ΔΔG| ${maxAbs.toFixed(4)} < floor 0.05 (all noise)`);
assert(scan1.rows.every((r) => r.noise === true), "floor honored: every live row tagged noise");
const sameOrder = scan1.rows.map((r) => r.label).join(",") === scan2.rows.map((r) => r.label).join(",")
  && scan1.rows.every((r, i) => r.ddG === scan2.rows[i].ddG);
assert(sameOrder, `ordering bit-stable across two runs (top ${scan1.rows[0].label} ${scan1.rows[0].ddG.toFixed(4)})`);
const top3 = scan1.rows.slice(0, 3).map((r) => r.label).join(", ");
assert(scan1.rows.slice(0, 3).some((r) => `${r.wtRes}${r.resSeq}` === "MET102"), `top-3 still contains cavity liner MET102 (${top3})`);
const liveTable = formatMutationTable(scan1.rows);
assert(liveTable.includes("~noise") && liveTable.includes("ranking only"), "live table flagged ~noise + disclaimer");

console.log(`\nTEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
process.exit(failed ? 1 : 0);
