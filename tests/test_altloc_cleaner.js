/**
 * tests/test_altloc_cleaner.js — Stage-4 altloc / duplicate-atom cleaner.
 *
 * Covers the shared cleaner rule (src/pdb_altloc.js) wired into parseCa,
 * parseHeavy (src/heavy.js) and parseLigands (src/pdb.js):
 *   highest occupancy wins; 'A' on tie (order-independent); otherwise first;
 *   zero-occupancy copies dropped when a non-zero copy exists; exact
 *   duplicates (same altLoc/coords, repeated serial) collapse to one atom.
 *
 * Also guards bit-identity on clean files: 4W52 (164 Cα / 1308 heavy /
 * 2 hetero groups) and 1CRN (46 Cα / 327 heavy) must parse to the exact
 * pre-Stage-4 counts (the old first-wins code kept 'A' on every 0.50/0.50
 * tie, and the new rule keeps the same record).
 *
 * Runnable: node tests/test_altloc_cleaner.js
 * Prints PASS/FAIL per assertion; exits 1 on any failure.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { parseCa, parseLigands } from "../src/pdb.js";
import { parseHeavy } from "../src/heavy.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

let passed = 0, failed = 0;
function assert(c, m) {
  if (c) { passed++; console.log(`  ✓ ${m}`); }
  else { failed++; console.error(`  ✗ FAIL: ${m}`); }
}

function findPdb(name) {
  const candidates = [
    path.resolve(process.cwd(), name),
    path.resolve(__dirname, "..", name),
    path.resolve(__dirname, name),
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`${name} not found (tried ${candidates.join(", ")})`);
}

/**
 * Build one fixed-column PDB ATOM/HETATM line.
 * Columns (1-indexed): 1–6 record, 7–11 serial, 13–16 name, 17 altLoc,
 * 18–20 resName, 22 chain, 23–26 resSeq, 27 iCode, 31–38/39–46/47–54 xyz,
 * 55–60 occupancy, 61–66 B-factor, 77–78 element.
 */
function pdbLine({ rec = "ATOM  ", serial = 1, name = "CA", alt = " ", res = "ALA", chain = "A", seq = 1, icode = " ", x = 0, y = 0, z = 0, occ = 1.0, b = 20.0, elem = "C" }) {
  const f8 = (v) => v.toFixed(3).padStart(8, " ").slice(-8);
  const f6 = (v) => v.toFixed(2).padStart(6, " ").slice(-6);
  let line = rec.slice(0, 6).padEnd(6, " ");
  line += String(serial).padStart(5, " ").slice(-5);
  line += " ";
  line += name.padStart(4, " ").slice(-4);
  line += alt.slice(0, 1) || " ";
  line += res.slice(0, 3).padEnd(3, " ");
  line += " ";
  line += chain.slice(0, 1) || " ";
  line += String(seq).padStart(4, " ").slice(-4);
  line += icode.slice(0, 1) || " ";
  line += "   ";
  line += f8(x) + f8(y) + f8(z);
  line += f6(occ) + f6(b);
  line += "          ";
  line += elem.slice(0, 2).padStart(2, " ");
  return line;
}

function main() {
  console.log("=== test_altloc_cleaner (Stage-4 input dedup) ===");

  // ---- (1) parseCa: highest occupancy wins ----
  const caOcc = [
    pdbLine({ serial: 1, name: "CA", alt: "A", res: "ALA", seq: 1, x: 0, y: 0, z: 0, occ: 0.30 }),
    pdbLine({ serial: 2, name: "CA", alt: "B", res: "ALA", seq: 1, x: 1, y: 0, z: 0, occ: 0.70 }),
    pdbLine({ serial: 3, name: "CA", alt: " ", res: "GLY", seq: 2, x: 5, y: 0, z: 0, occ: 1.0 }),
  ].join("\n");
  const r1 = parseCa(caOcc);
  assert(r1.beads.length === 2, `parseCa occupancy winner: 2 beads kept (got ${r1.beads.length})`);
  assert(Math.abs(r1.beads[0].x - 1) < 1e-9, `parseCa keeps highest-occupancy copy (B @x=1, got x=${r1.beads[0].x})`);

  // ---- (2) parseCa: tie keeps first 'A' ----
  const caTie = [
    pdbLine({ serial: 1, name: "CA", alt: "A", res: "ARG", seq: 14, x: 10, y: 0, z: 0, occ: 0.50 }),
    pdbLine({ serial: 2, name: "CA", alt: "B", res: "ARG", seq: 14, x: 11, y: 0, z: 0, occ: 0.50 }),
    pdbLine({ serial: 3, name: "CA", alt: " ", res: "GLY", seq: 15, x: 15, y: 0, z: 0 }),
    pdbLine({ serial: 4, name: "CA", alt: " ", res: "SER", seq: 16, x: 20, y: 0, z: 0 }),
  ].join("\n");
  const r2 = parseCa(caTie);
  assert(r2.beads.length === 3, `parseCa tie: 3 beads (got ${r2.beads.length})`);
  assert(Math.abs(r2.beads[0].x - 10) < 1e-9, `parseCa tie keeps first 'A' (x=10, got ${r2.beads[0].x})`);

  // ---- (3) parseCa: B-then-A tie still prefers 'A' ----
  const caTieBA = [
    pdbLine({ serial: 1, name: "CA", alt: "B", res: "ARG", seq: 14, x: 11, y: 0, z: 0, occ: 0.50 }),
    pdbLine({ serial: 2, name: "CA", alt: "A", res: "ARG", seq: 14, x: 10, y: 0, z: 0, occ: 0.50 }),
    pdbLine({ serial: 3, name: "CA", alt: " ", res: "GLY", seq: 15, x: 15, y: 0, z: 0 }),
    pdbLine({ serial: 4, name: "CA", alt: " ", res: "SER", seq: 16, x: 20, y: 0, z: 0 }),
  ].join("\n");
  const r3 = parseCa(caTieBA);
  assert(Math.abs(r3.beads[0].x - 10) < 1e-9, `parseCa B-then-A tie prefers 'A' (x=10, got ${r3.beads[0].x})`);

  // ---- (4) parseCa: zero-occupancy duplicate dropped ----
  const caZero = [
    pdbLine({ serial: 1, name: "CA", alt: "A", res: "SER", seq: 3, x: 9, y: 9, z: 9, occ: 0.00 }),
    pdbLine({ serial: 2, name: "CA", alt: "B", res: "SER", seq: 3, x: 4, y: 0, z: 0, occ: 1.00 }),
    pdbLine({ serial: 3, name: "CA", alt: " ", res: "GLY", seq: 4, x: 8, y: 0, z: 0 }),
    pdbLine({ serial: 4, name: "CA", alt: " ", res: "ALA", seq: 5, x: 12, y: 0, z: 0 }),
  ].join("\n");
  const r4 = parseCa(caZero);
  assert(r4.beads.length === 3 && Math.abs(r4.beads[0].x - 4) < 1e-9,
    `parseCa drops zero-occupancy copy (kept x=4, got ${r4.beads[0]?.x})`);

  // ---- (5) parseHeavy: occupancy rule on a sidechain atom ----
  const hvSide = [
    pdbLine({ serial: 1, name: "N", res: "MET", seq: 1, x: 0, y: 0, z: 0, elem: "N" }),
    pdbLine({ serial: 2, name: "CA", res: "MET", seq: 1, x: 1, y: 0, z: 0, elem: "C" }),
    pdbLine({ serial: 3, name: "CG", alt: "A", res: "MET", seq: 1, x: 2, y: 0, z: 0, occ: 0.30, elem: "C" }),
    pdbLine({ serial: 4, name: "CG", alt: "B", res: "MET", seq: 1, x: 2.5, y: 0, z: 0, occ: 0.70, elem: "C" }),
    pdbLine({ serial: 5, name: "N", res: "GLY", seq: 2, x: 5, y: 0, z: 0, elem: "N" }),
  ].join("\n");
  const h1 = parseHeavy(hvSide);
  assert(h1.atoms.length === 4, `parseHeavy occupancy winner: 4 atoms (got ${h1.atoms.length})`);
  const cg = h1.atoms.find((a) => a.atomName === "CG");
  assert(cg && Math.abs(cg.x - 2.5) < 1e-9, `parseHeavy keeps highest-occupancy sidechain copy (x=2.5, got ${cg?.x})`);

  // ---- (6) parseHeavy: exact-duplicate HETATM collapses to one ----
  const hvDup = [
    pdbLine({ rec: "HETATM", serial: 100, name: "C1", res: "UNK", seq: 200, x: 0, y: 0, z: 0, elem: "C" }),
    pdbLine({ rec: "HETATM", serial: 101, name: "C1", res: "UNK", seq: 200, x: 0, y: 0, z: 0, elem: "C" }),
    pdbLine({ rec: "HETATM", serial: 102, name: "C2", res: "UNK", seq: 200, x: 1.5, y: 0, z: 0, elem: "C" }),
  ].join("\n");
  const h2 = parseHeavy(hvDup);
  assert(h2.atoms.length === 2, `parseHeavy exact-duplicate HETATM → single atom (got ${h2.atoms.length})`);

  // ---- (7) parseLigands: exact-duplicate HETATM → single atom, sane bonds ----
  const ligDup = [
    pdbLine({ rec: "HETATM", serial: 100, name: "C1", res: "LIG", seq: 500, x: 0, y: 0, z: 0, elem: "C" }),
    pdbLine({ rec: "HETATM", serial: 101, name: "C1", res: "LIG", seq: 500, x: 0, y: 0, z: 0, elem: "C" }),
    pdbLine({ rec: "HETATM", serial: 102, name: "C2", res: "LIG", seq: 500, x: 1.5, y: 0, z: 0, elem: "C" }),
  ].join("\n");
  const mols = parseLigands(ligDup);
  assert(mols.length === 1 && mols[0].atoms.length === 2,
    `parseLigands exact duplicate → 2 atoms (got ${mols.length ? mols[0].atoms.length : "no mol"})`);
  assert(mols[0].bonds.length === 1, `parseLigands deduped pair bonds once (got ${mols[0].bonds.length})`);

  // ---- (8) parseLigands: altloc duplicate keeps occupancy winner, no coincidence ----
  const ligAlt = [
    pdbLine({ rec: "HETATM", serial: 100, name: "C1", alt: "A", res: "LIG", seq: 501, x: 0, y: 0, z: 0, occ: 0.30, elem: "C" }),
    pdbLine({ rec: "HETATM", serial: 101, name: "C1", alt: "B", res: "LIG", seq: 501, x: 0.5, y: 0, z: 0, occ: 0.70, elem: "C" }),
    pdbLine({ rec: "HETATM", serial: 102, name: "C2", res: "LIG", seq: 501, x: 1.5, y: 0, z: 0, elem: "C" }),
  ].join("\n");
  const mols2 = parseLigands(ligAlt);
  assert(mols2[0].atoms.length === 2 && Math.abs(mols2[0].atoms[0].x - 0.5) < 1e-9,
    `parseLigands altloc winner kept (x=0.5, got ${mols2[0]?.atoms[0]?.x})`);
  const dLig = Math.hypot(
    mols2[0].atoms[0].x - mols2[0].atoms[1].x,
    mols2[0].atoms[0].y - mols2[0].atoms[1].y,
    mols2[0].atoms[0].z - mols2[0].atoms[1].z);
  assert(dLig > 1e-6, `parseLigands deduped ligand has no coincident pair (d=${dLig.toFixed(3)} Å)`);

  // ---- (9) cleaner warnings logged (count skipped) ----
  assert((r1.warnings?.length ?? 0) >= 1 && (h2.warnings?.length ?? 0) >= 1,
    `cleaner logs skip counts (parseCa ${r1.warnings.length}, parseHeavy ${h2.warnings.length} warnings)`);

  // ---- (10) bit-identity on clean files ----
  const t4 = fs.readFileSync(findPdb("4w52.pdb"), "utf-8");
  const p4 = parseCa(t4);
  assert(p4.beads.length === 164, `4W52 Cα unchanged = 164 (got ${p4.beads.length})`);
  const h4 = parseHeavy(t4);
  assert(h4.atoms.length === 1308, `4W52 heavy unchanged = 1308 (got ${h4.atoms.length})`);
  assert(h4.heteroGroups.length === 2, `4W52 hetero groups unchanged = 2 (got ${h4.heteroGroups.length})`);
  const t1 = fs.readFileSync(findPdb("1crn.pdb"), "utf-8");
  const p1 = parseCa(t1);
  assert(p1.beads.length === 46, `1CRN Cα unchanged = 46 (got ${p1.beads.length})`);
  const h1c = parseHeavy(t1);
  assert(h1c.atoms.length === 327, `1CRN heavy unchanged = 327 (got ${h1c.atoms.length})`);

  console.log(`\n=== test_altloc_cleaner: ${passed} PASSED, ${failed} FAILED ===`);
  process.exit(failed ? 1 : 0);
}

main();
