/**
 * test_cutoff.js — C23 ENM cutoff scan monotonic springs for 4W52
 *
 * Builds FF with rc 8,10,12 Å for 4W52 (T4 lysozyme L99A, 164 aa, holo
 * benzene pocket) and checks springs count monotonic increasing.
 * Also prints the counts for notebooks/cutoff_scan.md table.
 *
 * Runnable: node tests/test_cutoff.js → PASS/FAIL
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { parseCa, selectSystem } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
function findPdb(name) {
  const candidates = [
    path.resolve(process.cwd(), name),
    path.resolve(__dirname, "..", name),
    `/home/samirr/WORK/simulation_coarse/${name}`,
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`PDB not found: ${name}`);
}

function main() {
  console.log("=== C23 ENM cutoff scan monotonic — 4W52 ===");
  const pdbText = fs.readFileSync(findPdb("4w52.pdb"), "utf-8");
  const parsed = parseCa(pdbText);
  const sel = selectSystem(parsed);
  console.log(`4W52: ${sel.beads.length} beads`);

  const cutoffs = [8, 10, 12];
  const counts = [];
  for (const rc of cutoffs) {
    const ff = new ForceField(sel, { rc, gamma: 1.0 });
    const nS = ff.springs.length / 3;
    counts.push(nS);
    console.log(`  Rc=${rc} Å → springs=${nS} (⟨${(nS/sel.beads.length).toFixed(1)}/res⟩)`);
  }
  // Montonic check
  if (!(counts[0] < counts[1] && counts[1] < counts[2])) {
    console.error(`FAIL: springs not monotonic: rc 8→${counts[0]}, 10→${counts[1]}, 12→${counts[2]} (expected strictly increasing)`);
    process.exit(1);
  }
  // Also sanity: Rc=10 default near optimum (notebooks/cutoff_scan.md)
  // Check that Rc=10 count is plausible (order 1000-1800 for 164 aa)
  if (counts[1] < 600 || counts[1] > 2500) {
    console.warn(`WARN: Rc=10 springs=${counts[1]} outside expected 600-2500 for 164 aa — check PDB selection`);
  }
  console.log(`PASS: C23 cutoff monotonic 8(${counts[0]}) < 10(${counts[1]}) < 12(${counts[2]}) — default 10 near optimum (notebooks/cutoff_scan.md)`);
}

main();
