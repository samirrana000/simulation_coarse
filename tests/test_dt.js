/**
 * tests/test_dt.js — G67 dt auto-tuning honest
 *
 * Verifies ligand-aware dt selection:
 *   CG alone 4fs (0.004 ps) vs CG+ligand 1.7fs (0.0017 ps) vs heavy 1fs (0.001 ps)
 *
 * Runnable: node tests/test_dt.js
 */

import { LangevinIntegrator } from "../src/integrator.js";

function mockFF({ n, nProt, heavy, nLigAtoms, ligandBondsLen }) {
  const masses = new Float64Array(n).fill(110);
  // Light ligand atoms: 12 Da for ligand region if present
  if (nLigAtoms > 0) {
    for (let i = nProt; i < n; i++) masses[i] = 12;
  }
  return {
    n,
    nProt,
    heavy: !!heavy,
    nLigAtoms,
    masses,
    covalentBonds: nLigAtoms ? new Float64Array([nProt, nProt+1, 1.4]) : new Float64Array([0,1,3.8]),
    ligandBonds: ligandBondsLen ? new Float64Array(new Array(ligandBondsLen).fill(0).map((_,i)=> i%3===2?1.4:i)) : new Float64Array(0),
    compute: () => 0,
    forces: new Float64Array(n*3),
    ref: new Float64Array(n*3),
  };
}

console.log("=== G67 dt auto-tuning: CG alone 4fs vs CG+ligand 1.7fs vs heavy 1fs ===");

// CG alone: n=164, heavy=false, no ligand
const ffCG = mockFF({ n:164, nProt:164, heavy:false, nLigAtoms:0, ligandBondsLen:0 });
const integCG = new LangevinIntegrator(new Float64Array(164*3), ffCG, 110);
console.log(`CG alone dt=${integCG.dt} ps (${(integCG.dt*1000).toFixed(1)}fs) — expected 0.004 ps / 4fs`);
if (Math.abs(integCG.dt - 0.004) > 1e-6) {
  console.error(`FAIL: CG alone dt ${integCG.dt} != 0.004`);
  process.exit(1);
}

// CG+ligand: n=170, ligand 6 atoms (benzene-like), heavy=false
const ffLig = mockFF({ n:170, nProt:164, heavy:false, nLigAtoms:6, ligandBondsLen: 18 });
const integLig = new LangevinIntegrator(new Float64Array(170*3), ffLig, 110);
console.log(`CG+ligand dt=${integLig.dt} ps (${(integLig.dt*1000).toFixed(2)}fs) — expected 0.0017 ps / 1.7fs`);
if (Math.abs(integLig.dt - 0.0017) > 1e-6) {
  console.error(`FAIL: CG+ligand dt ${integLig.dt} != 0.0017`);
  process.exit(1);
}

// Heavy: n=1308, heavy=true
const ffHeavy = mockFF({ n:1308, nProt:1308, heavy:true, nLigAtoms:0, ligandBondsLen:0 });
const integHeavy = new LangevinIntegrator(new Float64Array(1308*3), ffHeavy, 110);
console.log(`heavy dt=${integHeavy.dt} ps (${(integHeavy.dt*1000).toFixed(1)}fs) — expected 0.001 ps / 1fs`);
if (Math.abs(integHeavy.dt - 0.001) > 1e-6) {
  console.error(`FAIL: heavy dt ${integHeavy.dt} != 0.001`);
  process.exit(1);
}

console.log("PASS: dt values — CG alone 4fs vs CG+ligand 1.7fs vs heavy 1fs (G67)");
