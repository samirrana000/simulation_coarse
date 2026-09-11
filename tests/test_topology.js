/**
 * test_topology.js — B12 covalent radii / BOND_SLACK / cap 2.2 Å fidelity
 *
 * Loads 1crn.pdb (crambin, 46 aa, 3 SSBOND: CYS3-40 2.00, CYS4-32 2.04, CYS16-26 2.05)
 * via src/heavy.js parseHeavy + buildTopology. Asserts:
 *   - bonds >=3 disulfides (S–S) detected (CSD S radius 1.02, BOND_SLACK 1.15, cap 2.2 captures 2.04)
 *   - no spurious Ca–N 2.9 Å bond (synthetic Ca+N at 2.9 Å correctly rejected)
 * Runnable: node tests/test_topology.js → PASS/FAIL
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { parseHeavy, buildTopology } from "../src/heavy.js";
import { COVALENT_RADIUS, BOND_SLACK } from "../src/ff-params.js";

function findPdb(name) {
  const candidates = [
    name,
    path.join("..", name),
    path.join(path.dirname(fileURLToPath(import.meta.url)), "..", name),
    `/home/samirr/WORK/simulation_coarse/${name}`,
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`PDB not found: ${name} tried ${candidates.join(", ")}`);
}

function assert(cond, msg) {
  if (!cond) { console.error(`FAIL: ${msg}`); process.exit(1); }
  console.log(`  ✓ ${msg}`);
}

function main() {
  console.log("=== B12 Topology — CSD radii / BOND_SLACK 1.15 / cap 2.2 Å ===");
  console.log(`COVALENT_RADIUS S=${COVALENT_RADIUS.S} CA=${COVALENT_RADIUS.CA} N=${COVALENT_RADIUS.N} BOND_SLACK=${BOND_SLACK}`);

  // --- crambin disulfides ---
  const pdbPath = findPdb("1crn.pdb");
  const pdbText = fs.readFileSync(pdbPath, "utf-8");
  const parsed = parseHeavy(pdbText);
  const topo = buildTopology(parsed.atoms);
  console.log(`Parsed ${parsed.atoms.length} heavy atoms, bonds=${topo.bonds.length}`);

  // count S–S bonds (SG SG pairs)
  let ssCount = 0;
  let ssDists = [];
  for (const [i, j] of topo.bonds) {
    const A = parsed.atoms[i], B = parsed.atoms[j];
    if (A.element === "S" && B.element === "S") {
      const dx = A.x - B.x, dy = A.y - B.y, dz = A.z - B.z;
      const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
      ssCount++;
      ssDists.push(r.toFixed(3));
    }
  }
  console.log(`S–S disulfide bonds detected: ${ssCount} distances [${ssDists.join(", ")}]`);
  assert(ssCount >= 3, `disulfides >=3 (CSD S–S 2.04 Å captured via S 1.02*2=2.04 *1.15=2.35 cap 2.20) got ${ssCount}`);

  // also verify known 1crn SG pairs are within cap
  // Expect S–S 2.04 < 2.20 and < BOND_SLACK*sum validated
  const sum = COVALENT_RADIUS.S + COVALENT_RADIUS.S;
  assert(BOND_SLACK * sum > 2.04 && 2.2 > 2.04, `BOND_SLACK*sum=${(BOND_SLACK*sum).toFixed(2)} cap 2.2 correctly captures S–S 2.04 Å`);

  // --- spurious Ca–N 2.9 Å rejection ---
  // Synthetic two-atom system: CA at origin, N at 2.9 Å on x
  const synthetic = [
    { element: "CA", x: 0, y: 0, z: 0, isMetal: false },
    { element: "N",  x: 2.9, y: 0, z: 0, isMetal: false },
  ];
  const topoSyn = buildTopology(synthetic);
  const caN_bonded = topoSyn.bonds.some(([a,b]) => (a===0 && b===1));
  console.log(`Synthetic Ca–N 2.9 Å bond detected? ${caN_bonded} (should be false; Ca 1.76+N 0.75=2.51*1.15=2.89→2.9 fails & cap 2.2)`);
  assert(!caN_bonded, `no spurious Ca–N 2.9 Å bond (Ca 1.76+N 0.75=2.51*1.15=2.89, 2.9>2.89 && >2.2 cap)`);

  // also check that same Ca–N at 2.1 Å WOULD bond? Optional sanity
  const close = [
    { element: "CA", x: 0, y: 0, z: 0, isMetal: false },
    { element: "N",  x: 2.0, y: 0, z: 0, isMetal: false },
  ];
  const topoClose = buildTopology(close);
  const closeBonded = topoClose.bonds.length > 0;
  console.log(`Sanity Ca–N 2.0 Å bond? ${closeBonded} (should be true, inside cap and slack)`);

  // Also verify no Ca–N bonds in real 1crn topology with r ~2.9 falsely typed?
  // Scan 1crn for any Ca–N heavy pair <2.2 that would be bonded — expect 0 spurious long
  // We already synthetic-tested; real check just that no Ca–N bond >2.2 exists (by construction)
  let longCaN = 0;
  for (const [i,j] of topo.bonds) {
    const A = parsed.atoms[i], B = parsed.atoms[j];
    if ((A.element==="CA" && B.element==="N") || (A.element==="N" && B.element==="CA")) {
      const dx=A.x-B.x, dy=A.y-B.y, dz=A.z-B.z;
      const r=Math.sqrt(dx*dx+dy*dy+dz*dz);
      if (r>2.2) longCaN++;
    }
  }
  assert(longCaN===0, `no Ca–N bonds >2.2 Å in 1crn topology (cap enforced) got ${longCaN}`);

  console.log("PASS: B12 topology — disulfides + Ca–N rejection");
}

main();
