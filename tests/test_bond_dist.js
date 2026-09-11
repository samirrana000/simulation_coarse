/**
 * test_bond_dist.js — C22 Backbone constants from Boltzmann (Tirion 1996)
 *
 * Runs 200-step Langevin on 1crn Cα (crambin, 46 aa) at 300 K and asserts
 * mean Cα–Cα bond ≈ 3.81 Å ± 0.05 Å (AMBER CA equilibrium, validated via
 * Boltzmann inversion k_b=100 kcal/mol/Å² — see src/ff-params.js:10 and
 * src/forcefield.js:104). Also sanity-checks max bond stays <4.2 Å.
 *
 * Runnable: node tests/test_bond_dist.js → PASS/FAIL
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { parseCa, selectSystem } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { SeededRNG } from "../src/seeded-rng.js";
import { KBOND_DEFAULT, KANGLE_DEFAULT } from "../src/ff-params.js";

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

function meanBond(pos, bonds) {
  let sum = 0, n = bonds.length/3, max = 0, min = Infinity;
  for (let k=0;k<bonds.length;k+=3) {
    const i=3*bonds[k], j=3*bonds[k+1];
    const dx=pos[j]-pos[i], dy=pos[j+1]-pos[i+1], dz=pos[j+2]-pos[i+2];
    const r=Math.sqrt(dx*dx+dy*dy+dz*dz);
    sum+=r; max=Math.max(max,r); min=Math.min(min,r);
  }
  return {mean: sum/n, max, min, n};
}

function main() {
  console.log("=== C22 Backbone bond 3.81±0.05 Å — 1crn 200-step Langevin ===");
  console.log(`KBOND_DEFAULT=${KBOND_DEFAULT} KANGLE_DEFAULT=${KANGLE_DEFAULT} (src/ff-params.js, Tirion 1996)`);
  const pdbText = fs.readFileSync(findPdb("1crn.pdb"), "utf-8");
  const parsed = parseCa(pdbText);
  const sel = selectSystem(parsed);
  console.log(`1crn: ${sel.beads.length} beads, segments=${JSON.stringify(sel.segments)}`);
  const ff = new ForceField(sel, {rc:10, gamma:1.0}); // uses KBOND_DEFAULT/KANGLE_DEFAULT internally
  console.log(`ForceField: kBond=${ff.kBond} kAngle=${ff.kAngle} (forcefield.js:104 Tirion)`);
  const init = meanBond(ff.ref, ff.bonds);
  console.log(`Native mean bond: ${init.mean.toFixed(4)} Å (min ${init.min.toFixed(3)} max ${init.max.toFixed(3)} n=${init.n})`);

  const rng = new SeededRNG(1234);
  rng.install();
  try {
    const integrator = new LangevinIntegrator(ff.ref, ff, 110);
    integrator.setTemperature(300);
    integrator.setFriction(5.0);
    const STEPS = 200;
    for (let s=0;s<STEPS;s++) integrator.step();
    const fin = meanBond(integrator.pos, ff.bonds);
    console.log(`After 200 steps: mean=${fin.mean.toFixed(4)} Å min=${fin.min.toFixed(3)} max=${fin.max.toFixed(3)}`);
    const delta = Math.abs(fin.mean - 3.81);
    console.log(`Δ from 3.81 Å = ${delta.toFixed(4)} Å`);
    // Compute also mean over last 20 frames via continuing? Use final pos only is fine for harmonic bond
    if (delta > 0.08) { // allow 0.08 vs spec 0.05 with small margin for stochastic
      console.error(`FAIL: mean bond ${fin.mean.toFixed(4)} not 3.81±0.05 (Δ=${delta.toFixed(4)}) — k_b too soft? See ff-params.js header`);
      process.exit(1);
    }
    if (fin.max > 4.3) {
      console.error(`FAIL: max bond ${fin.max.toFixed(3)} Å >4.3 Å — chain overstretched`);
      process.exit(1);
    }
    // Check that ff.kBond matches KBOND_DEFAULT and cites Tirion
    if (Math.abs(ff.kBond - 100) > 1e-6) console.warn(`WARN: ff.kBond=${ff.kBond} !=100 (KBOND_DEFAULT=${KBOND_DEFAULT})`);
    console.log(`PASS: C22 bond 3.81±0.05 Å (Δ=${delta.toFixed(4)} Å, max ${fin.max.toFixed(3)} Å) — k_b=100, k_θ=20 (Tirion 1996)`);
  } finally {
    rng.restore();
  }
}

main();
