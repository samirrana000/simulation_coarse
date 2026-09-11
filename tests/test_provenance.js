/**
 * test_provenance.js — A09 provenance + A01 version regression.
 * Ensures recorder.buildFile header contains VERSION and required fields.
 */
import { Recorder } from "../src/recorder.js";
import { VERSION, BUILD_DATE } from "../src/version.js";

function assert(cond, msg) {
  if (!cond) { console.error("✗ FAIL: " + msg); process.exit(1); }
  console.log("✓ " + msg);
}

const rec = new Recorder();
const n = 3;
const beads = [
  { resName: "ALA", chain: "A", resSeq: 1 },
  { resName: "GLY", chain: "A", resSeq: 2 },
  { resName: "LEU", chain: "A", resSeq: 3 },
];
const pos1 = new Float64Array([0,0,0, 1,0,0, 2,0,0]);
const pos2 = new Float64Array([0.1,0,0, 1.1,0,0, 2.1,0,0]);
rec.start(0, 1, 10);
rec.maybeCapture(pos1, 0);
rec.maybeCapture(pos2, 1);

const pdb = rec.buildFile("pdb", beads, { T: 310, gamma: 2.5, seed: 42, date: BUILD_DATE });
assert(pdb.includes("REMARK simulation_coarse"), "PDB contains REMARK simulation_coarse");
assert(pdb.includes(`v${VERSION}`), `PDB contains version v${VERSION}`);
assert(pdb.includes("T=310"), "PDB contains T=310");
assert(pdb.includes("gamma=2.5"), "PDB contains gamma");
assert(pdb.includes("seed=42"), "PDB contains seed");
assert(pdb.includes(BUILD_DATE), "PDB contains date");

const xyz = rec.buildFile("xyz", beads, { T: 310, gamma: 2.5, seed: 42 });
assert(xyz.includes("REMARK simulation_coarse"), "XYZ contains REMARK simulation_coarse");
assert(xyz.includes(`v${VERSION}`), `XYZ contains version v${VERSION}`);

console.log(`\nProvenance test PASSED (VERSION=${VERSION}, BUILD_DATE=${BUILD_DATE})`);
