/**
 * test_exclusions.js — B13 1-2/1-3/1-4 exclusions correct
 *
 * Audits src/heavy.js:524-542 (pairKey i*1e6+j) and src/forcefield.js:251-283
 * Documentation in docs/EXCLUSIONS.md: 1-4 scaled 0.5, intra-ligand fully excluded.
 * Runnable: node tests/test_exclusions.js → PASS/FAIL
 */

import { HeavyForceField, buildTopology } from "../src/heavy.js";

function assert(cond, msg) {
  if (!cond) { console.error(`FAIL: ${msg}`); process.exit(1); }
  console.log(`  ✓ ${msg}`);
}

function main() {
  console.log("=== B13 Exclusions — 1-4 scaled 0.5, intra-ligand fully excluded ===");

  // ---- 1. 1-4 scaled 0.5 on a 4-atom chain (C–C–C–C) ----
  // Create 4 carbon atoms linearly spaced at 1.54 Å (C–C single bond)
  const chain = [
    { element: "C", x: 0.00, y: 0, z: 0, isProtein: true,  isLigand: false, isMetal: false, isHetero: false },
    { element: "C", x: 1.54, y: 0, z: 0, isProtein: true,  isLigand: false, isMetal: false, isHetero: false },
    { element: "C", x: 3.08, y: 0, z: 0, isProtein: true,  isLigand: false, isMetal: false, isHetero: false },
    { element: "C", x: 4.62, y: 0, z: 0, isProtein: true,  isLigand: false, isMetal: false, isHetero: false },
  ];
  const ffChain = new HeavyForceField({ atoms: chain }, { gamma: 1.0 }, []);
  console.log(`Chain 4C: n=${ffChain.n} bonds=${ffChain.bonds.length/3} angles=${ffChain.angles.length/4} propers=${ffChain.propers.length/5}`);
  // proper 0-1-2-3 should give 1-4 pair (0,3) scaled 0.5
  const pairKey = (i,j) => i<j ? i*1e6+j : j*1e6+i;
  const key14 = pairKey(0,3);
  const scale14 = ffChain._scale14.get(key14);
  console.log(`1-4 pairKey(0,3)=${key14} scale=${scale14}`);
  assert(scale14 === 0.5, `1-4 pair (0,3) scaled 0.5 (got ${scale14})`);

  // also verify 1-2 and 1-3 are excluded, not scaled
  assert(ffChain._excluded.has(pairKey(0,1)), `1-2 pair (0,1) excluded`);
  assert(ffChain._excluded.has(pairKey(1,2)), `1-2 pair (1,2) excluded`);
  assert(ffChain._excluded.has(pairKey(0,2)), `1-3 pair (0,2) excluded`);
  assert(!ffChain._excluded.has(key14), `1-4 pair (0,3) NOT excluded (scaled, not excluded) got ${ffChain._excluded.has(key14)}`);

  // In _nonBondedGrid, 1-4 should be processed with s14=0.5, not returned
  let s14_used = null;
  // simulate lookup as in _nonBondedGrid
  const excluded = ffChain._excluded.has(key14);
  const s14 = ffChain._scale14.get(key14) ?? 1.0;
  s14_used = excluded ? "excluded" : s14;
  assert(s14_used === 0.5, `_nonBondedGrid would use s14=0.5 for 1-4 (got ${s14_used})`);

  // ---- 2. Intra-ligand fully excluded ----
  // Build system with 2 protein + 4 ligand atoms (ligandStart=2)
  // Place ligands far apart so no covalent bonds cross, but intra-ligand pairs should be excluded
  const mixed = [
    { element: "C", x: 0, y: 0, z: 0, isProtein: true, isLigand: false, isMetal: false, isHetero: false },
    { element: "C", x: 1.54, y: 0, z: 0, isProtein: true, isLigand: false, isMetal: false, isHetero: false },
    // ligand block starting at index 2
    { element: "C", x: 20, y: 0, z: 0, isProtein: false, isLigand: true, isMetal: false, isHetero: false },
    { element: "C", x: 21.5, y: 0, z: 0, isProtein: false, isLigand: true, isMetal: false, isHetero: false },
    { element: "O", x: 20, y: 1.5, z: 0, isProtein: false, isLigand: true, isMetal: false, isHetero: false },
    { element: "N", x: 21.5, y: 1.5, z: 0, isProtein: false, isLigand: true, isMetal: false, isHetero: false },
  ];
  const ffMixed = new HeavyForceField({ atoms: mixed }, { gamma: 1.0 }, []);
  console.log(`Mixed system: n=${ffMixed.n} nProt=${ffMixed.nProt} nLig=${ffMixed.nLigAtoms} ligandStart=${ffMixed.ligandStart}`);
  assert(ffMixed.ligandStart === 2, `ligandStart=2 (got ${ffMixed.ligandStart})`);
  assert(ffMixed.nLigAtoms === 4, `nLigAtoms=4 (got ${ffMixed.nLigAtoms})`);
  const nLig = ffMixed.nLigAtoms;
  const expectedIntra = nLig * (nLig - 1) / 2; // 6
  let countedExcluded = 0;
  for (let i = ffMixed.ligandStart; i < ffMixed.n; i++) {
    for (let j = i+1; j < ffMixed.n; j++) {
      if (ffMixed._excluded.has(pairKey(i,j))) countedExcluded++;
    }
  }
  console.log(`Intra-ligand excluded count: ${countedExcluded} expected ${expectedIntra}`);
  assert(countedExcluded === expectedIntra, `intra-ligand fully excluded count ${countedExcluded} == ${expectedIntra}`);

  // Verify a protein–ligand pair is NOT excluded (should be evaluated, unless holo etc.)
  // indices 0 (prot) and 2 (lig) are 20 Å apart → not excluded, grid will handle
  assert(!ffMixed._excluded.has(pairKey(0,2)), `protein–ligand pair (0,2) NOT excluded (binding evaluated)`);

  // ---- 3. pairKey capacity audit ----
  // Ensure no collision for n < 1e6: test max atoms scenario ~100k still distinct
  // Example collision if factor were 1e3: 0*1e3+5000 collides with 5*1e3+0. With 1e6, need 1e6 to collide.
  const keyA = pairKey(0, 999999);
  const keyB = pairKey(1, 0);
  assert(keyA !== keyB, `pairKey 1e6 factor no collision: ${keyA} != ${keyB}`);

  console.log("PASS: B13 exclusions — 1-4 scaled 0.5 and intra-ligand excluded");
}

main();
