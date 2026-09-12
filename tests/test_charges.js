/**
 * test_charges.js — Loop-2 S1: CG formal bead charges (R2 term b, salt bridges)
 *
 * Validates the CG_FORMAL_CHARGES opt-in (par.binding.charges === true) that
 * revives the dead screened-Coulomb path in ff-binding.js:
 *
 *   E_coul = 332.0637 · q_i · q_a / (ε(r) · r) · sw(r),   ε(r) = 4 + 76·tanh(r/8)
 *
 * [1] Charge table: 4W52 Cα beads get Asp/Glu −1, Lys/Arg +1, all others
 *     (incl. His — neutral default, HIP hookup deferred) 0 when charges are
 *     ON; default ForceField keeps _protQ all-zero (pre-S1 behavior).
 * [2] Salt bridge: acetate (ligandLib) placed with O⁻ at 2.66 Å from the
 *     Lys35 Cα bead — the native HEPES O1S⁻···H₃N⁺–Lys35 distance (R2 §0) —
 *     produces a nonzero, attractive electrostatic contribution to bindingU.
 * [3] Analytic pair check: with only Lys35 charged, the measured ΔU equals
 *     the closed-form screened-Coulomb sum over ligand atoms to <1e-9.
 * [4] Screening profile: |ΔU_coul| decays with r (ε(r) grows, sw(r) → 0).
 * [5] Default-off bit-identity: a 10-step 4W52 + native-ligand Langevin run
 *     with default parameters reproduces the pre-S1 captured reference
 *     (energy/posHash/velHash) exactly — the flag flips nothing silently.
 *
 * Runnable: node tests/test_charges.js
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { parseCa, parseLigands, selectSystem } from "../src/pdb.js";
import { parseMol2 } from "../src/mol2.js";
import { LIGAND_LIBRARY } from "../src/ligandLib.js";
import { CG_FORMAL_CHARGES } from "../src/ff-params.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { SeededRNG } from "../src/seeded-rng.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

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

function findPdb() {
  const candidates = [
    path.resolve(process.cwd(), "4w52.pdb"),
    path.resolve(__dirname, "..", "4w52.pdb"),
    path.resolve(__dirname, "4w52.pdb"),
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`4w52.pdb not found (tried ${candidates.join(", ")})`);
}

/** FNV-1a hash over an array quantized to 1/scale (same as test_golden.js). */
function arrayHash(arr, scale = 1000) {
  let h = 0x811c9dc5;
  for (let i = 0; i < arr.length; i++) {
    const q = Math.round(arr[i] * scale);
    h ^= (q & 0xff);
    h = Math.imul(h, 0x01000193);
    h ^= ((q >> 8) & 0xff);
    h = Math.imul(h, 0x01000193);
    h ^= ((q >> 16) & 0xff);
    h = Math.imul(h, 0x01000193);
  }
  return (h >>> 0).toString(16).padStart(8, "0");
}

// Screened-Coulomb closed form, replicated from ff-binding.js PASS 1 exactly
// (ε(r) = 4 + 76·tanh(r/8); smoothstep sw = 1 for r ≤ 0.85·rc, 0 at rc = 9).
const ELC = 332.0637;
const epsr = (r) => 4 + 76 * Math.tanh(r / 8);
const swf = (r) => {
  const rc = 9.0, rsw = 0.85 * rc, inv = 1 / (rc - rsw);
  if (r <= rsw) return 1;
  const t = (r - rsw) * inv, t2 = t * t;
  return 1 - t2 * t * (10 - 15 * t + 6 * t2);
};

/**
 * Build a 4W52 CG ForceField with the acetate ligand placed so its O1
 * (carboxylate O, q = −0.5 via LIG_ELEMENT) sits exactly `dist` Å from the
 * Lys35 Cα bead along the outward COM direction. Molecule is mirrored so the
 * carbon skeleton points away from the bead (O1 is the closest atom).
 */
function buildAcetateSystem(sel, lys35Idx, lys35Pos, com, dist, charges) {
  const lib = LIGAND_LIBRARY.find((m) => m.id === "acetate");
  const mol = parseMol2(lib.mol2)[0];
  // mirror x → O1 becomes the leftmost atom, then translate to the target
  const dir = [lys35Pos[0] - com[0], lys35Pos[1] - com[1], lys35Pos[2] - com[2]];
  const dn = Math.hypot(...dir);
  const u = dir.map((v) => v / dn);
  const o1 = mol.atoms.find((a) => a.atomName === "O1");
  const target = lys35Pos.map((v, k) => v + u[k] * dist);
  const shift = [target[0] - o1.x, target[1] - o1.y, target[2] - o1.z];
  for (const a of mol.atoms) {
    const x = -a.x + shift[0], y = a.y + shift[1], z = a.z + shift[2];
    a.x = x; a.y = y; a.z = z;
  }
  const par = { rc: 10, gamma: 1.0, binding: { charges, holo: false } };
  return { ff: new ForceField(sel, par, [mol]), mol };
}

function main() {
  console.log("=== Loop-2 S1: CG formal charges (R2 term b — salt bridges) ===");
  const pdbText = fs.readFileSync(findPdb(), "utf-8");

  // ---------------------------------------------------------------------
  // [1] Charge table on 4W52 (164 Cα beads)
  // ---------------------------------------------------------------------
  console.log("\n[1] Charge table: parseCa → selectSystem → ForceField(4W52)...");
  const parsed = parseCa(pdbText);
  const sel = selectSystem(parsed);
  const ligands = parseLigands(pdbText);
  const ffOn = new ForceField(sel, { rc: 10, gamma: 1.0, binding: { charges: true } }, ligands);
  const ffDef = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligands);

  let nNeg = 0, nPos = 0, nZero = 0, nHis = 0, bad = 0;
  sel.beads.forEach((b, i) => {
    const q = ffOn._protQ[i];
    if (b.resName === "ASP" || b.resName === "GLU") {
      if (q === -1) nNeg++; else bad++;
    } else if (b.resName === "LYS" || b.resName === "ARG") {
      if (q === 1) nPos++; else bad++;
    } else {
      if (b.resName === "HIS" && q === 0) nHis++;
      if (q === 0) nZero++; else bad++;
    }
  });
  assert(bad === 0, `every bead charge matches residue identity (${nNeg} Asp/Glu −1, ${nPos} Lys/Arg +1, ${nZero} others 0, ${bad} mismatches)`);
  assert(nHis > 0, `His beads stay neutral by default (${nHis} His, q = 0; HIP +1 hookup deferred)`);
  assert(nNeg === 18 && nPos === 26, `4W52 charged census: Asp/Glu = 18, Lys/Arg = 26 (got ${nNeg}/${nPos})`);
  let allZero = true;
  for (let i = 0; i < ffDef._protQ.length; i++) if (ffDef._protQ[i] !== 0) allZero = false;
  assert(allZero && ffDef._protQ.length === 164, `default ForceField keeps _protQ all-zero (pre-S1 behavior, 164 beads)`);
  assert(ffOn.chargesOn === true && ffDef.chargesOn === false, `opt-in flag par.binding.charges (on=${ffOn.chargesOn}, default=${ffDef.chargesOn})`);
  assert(CG_FORMAL_CHARGES.ASP === -1 && CG_FORMAL_CHARGES.GLU === -1
    && CG_FORMAL_CHARGES.LYS === 1 && CG_FORMAL_CHARGES.ARG === 1
    && CG_FORMAL_CHARGES.HIS === undefined, `CG_FORMAL_CHARGES map exported (ASP/GLU −1, LYS/ARG +1, HIS absent ⇒ 0)`);

  // ---------------------------------------------------------------------
  // [2] Salt bridge: acetate O⁻ at 2.66 Å from Lys35 Cα
  // ---------------------------------------------------------------------
  console.log("\n[2] Salt-bridge ΔU: acetate (ligandLib) O⁻ near Lys35 (+1)...");
  const lysIdx = sel.beads.findIndex((b) => b.resName === "LYS" && b.resSeq === 35 && b.chain === "A");
  assert(lysIdx >= 0 && ffOn._protQ[lysIdx] === 1, `Lys35 A bead found at index ${lysIdx} with q = +1`);
  const lysPos = [sel.beads[lysIdx].x, sel.beads[lysIdx].y, sel.beads[lysIdx].z];
  const com = [0, 0, 0];
  for (const b of sel.beads) { com[0] += b.x; com[1] += b.y; com[2] += b.z; }
  com[0] /= sel.beads.length; com[1] /= sel.beads.length; com[2] /= sel.beads.length;

  const R_SB = 2.66; // native HEPES O1S⁻···H₃N⁺–Lys35 distance (R2 §0)
  const { ff: ffSbOn, mol } = buildAcetateSystem(sel, lysIdx, lysPos, com, R_SB, true);
  const { ff: ffSbOff } = buildAcetateSystem(sel, lysIdx, lysPos, com, R_SB, false);
  assert(mol.atoms.length === 4, `acetate parsed from ligandLib (4 heavy atoms)`);

  // identical geometry → ΔU is exactly the revived Coulomb contribution
  const posOn = ffSbOn.ref, posOff = ffSbOff.ref;
  let geoIdentical = posOn.length === posOff.length;
  for (let i = 0; i < posOn.length && geoIdentical; i++) if (posOn[i] !== posOff[i]) geoIdentical = false;
  assert(geoIdentical, `charges-on / charges-off systems share bit-identical coordinates`);

  ffSbOn.compute(posOn);
  ffSbOff.compute(posOff);
  const dUSb = ffSbOn.bindingU - ffSbOff.bindingU;
  const dUSbTot = ffSbOn.energy - ffSbOff.energy;
  console.log(`    bindingU(on) = ${ffSbOn.bindingU.toFixed(4)}, bindingU(off) = ${ffSbOff.bindingU.toFixed(4)}`);
  console.log(`    ΔU_coul(salt bridge @ ${R_SB} Å) = ${dUSb.toFixed(4)} kcal/mol (total-energy Δ = ${dUSbTot.toFixed(4)})`);
  assert(dUSb < -0.5, `charged ligand near Lys yields nonzero attractive electrostatic contribution (ΔU_coul = ${dUSb.toFixed(3)} kcal/mol)`);
  assert(Math.abs(dUSb - dUSbTot) < 1e-12, `ΔU lives entirely in the binding pass (Δ = ΔbindingU)`);

  // ---------------------------------------------------------------------
  // [3] Analytic single-pair check (only Lys35 charged)
  // ---------------------------------------------------------------------
  console.log("\n[3] Analytic screened-Coulomb check (isolated Lys35 pair)...");
  const { ff: ffIso } = buildAcetateSystem(sel, lysIdx, lysPos, com, R_SB, true);
  for (let i = 0; i < ffIso._protQ.length; i++) ffIso._protQ[i] = 0;
  ffIso._protQ[lysIdx] = 1; // only the salt-bridge partner carries charge
  ffIso.compute(ffIso.ref);
  const { ff: ffIsoOff } = buildAcetateSystem(sel, lysIdx, lysPos, com, R_SB, false);
  ffIsoOff.compute(ffIsoOff.ref);
  const dUIso = ffIso.bindingU - ffIsoOff.bindingU;

  // closed form: Σ_a 332.0637·q_i·q_a/(ε(r)·r)·sw(r) over acetand atoms in cutoff
  let expect = 0;
  for (const a of ffIso.ligandAtoms) {
    const qa = ffIso._ligQ[a.idx - ffIso.nProt];
    if (qa === 0) continue;
    const dx = ffIso.ref[3 * a.idx] - lysPos[0];
    const dy = ffIso.ref[3 * a.idx + 1] - lysPos[1];
    const dz = ffIso.ref[3 * a.idx + 2] - lysPos[2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
    if (r >= 9.0) continue;
    expect += ELC * 1 * qa / (epsr(r) * r) * swf(r);
  }
  console.log(`    measured ΔU_coul = ${dUIso.toFixed(6)}, analytic = ${expect.toFixed(6)} kcal/mol`);
  assert(Math.abs(dUIso - expect) < 1e-9, `isolated-pair ΔU_coul matches the closed-form screened Coulomb to <1e-9`);
  console.log(`    per-pair number: Lys(+1)···O⁻(−0.5) @ ${R_SB} Å → ${dUIso.toFixed(2)} kcal/mol (R2 §D predicted −2.20)`);

  // ---------------------------------------------------------------------
  // [4] Screening profile with distance
  // ---------------------------------------------------------------------
  console.log("\n[4] Screening profile ε(r)/sw(r): ΔU_coul decays with r...");
  const dUs = [];
  const dUsExpect = [];
  for (const r of [2.66, 4.0, 6.0, 8.0]) {
    const { ff: fOn } = buildAcetateSystem(sel, lysIdx, lysPos, com, r, true);
    const { ff: fOff } = buildAcetateSystem(sel, lysIdx, lysPos, com, r, false);
    for (let i = 0; i < fOn._protQ.length; i++) fOn._protQ[i] = 0;
    fOn._protQ[lysIdx] = 1;
    fOn.compute(fOn.ref);
    fOff.compute(fOff.ref);
    const d = fOn.bindingU - fOff.bindingU;
    dUs.push(d);
    // closed form for the same geometry (all in-cutoff acetate atoms)
    let ex = 0;
    for (const a of fOn.ligandAtoms) {
      const qa = fOn._ligQ[a.idx - fOn.nProt];
      if (qa === 0) continue;
      const dx = fOn.ref[3 * a.idx] - lysPos[0];
      const dy = fOn.ref[3 * a.idx + 1] - lysPos[1];
      const dz = fOn.ref[3 * a.idx + 2] - lysPos[2];
      const rr = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (rr >= 9.0) continue;
      ex += ELC * 1 * qa / (epsr(rr) * rr) * swf(rr);
    }
    dUsExpect.push(ex);
    console.log(`    r = ${r} Å: ΔU_coul = ${d.toFixed(4)} kcal/mol (analytic ${ex.toFixed(4)}, ε(r) = ${epsr(r).toFixed(1)})`);
  }
  let monotone = true;
  for (let i = 1; i < dUs.length; i++) if (Math.abs(dUs[i]) > Math.abs(dUs[i - 1]) + 1e-12) monotone = false;
  assert(monotone, `|ΔU_coul| monotonically screened with distance (ε 28→69, sw →0 at 9 Å)`);
  let maxDev = 0;
  for (let i = 0; i < dUs.length; i++) maxDev = Math.max(maxDev, Math.abs(dUs[i] - dUsExpect[i]));
  assert(maxDev < 1e-9, `every r point matches the closed-form screened Coulomb (max dev ${maxDev.toExponential(2)})`);
  assert(Math.abs(dUs[0]) > Math.abs(dUs[3]) * 10, `switch dominates by the cutoff edge: |ΔU(2.66 Å)| = ${Math.abs(dUs[0]).toFixed(2)} ≫ |ΔU(8.0 Å)| = ${Math.abs(dUs[3]).toFixed(3)}`);

  // ---------------------------------------------------------------------
  // [5] Default-off bit-identity vs pre-S1 reference (10 steps, seed 42)
  // ---------------------------------------------------------------------
  console.log("\n[5] Default (charges off) = pre-S1 behavior, bit-identical 10-step run...");
  // Captured BEFORE the S1 change on this machine (node /tmp/opencode/
  // ref_default_binding.mjs at HEAD~): same parseCa+selectSystem+parseLigands,
  // ForceField(sel, {rc:10, gamma:1.0}, ligands) — default binding params,
  // LangevinIntegrator T=300, zeta=5, SeededRNG(42), 10 BAOAB steps.
  const REF = { energy: 12.421069514062914, bindingU: -7.10639466730052, posHash: "999b22c9", velHash: "dbf98762" };
  const rng = new SeededRNG(42);
  rng.install();
  let energy, posHash, velHash, bindingU;
  try {
    const ff = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligands);
    const it = new LangevinIntegrator(ff.ref, ff, 110);
    it.setTemperature(300);
    it.setFriction(5.0);
    for (let s = 0; s < 10; s++) it.step();
    energy = ff.energy;
    bindingU = ff.bindingU;
    posHash = arrayHash(it.pos, 1000);
    velHash = arrayHash(it.vel, 1000);
  } finally {
    rng.restore();
  }
  console.log(`    current: energy=${energy} posHash=${posHash} velHash=${velHash}`);
  console.log(`    pre-S1 : energy=${REF.energy} posHash=${REF.posHash} velHash=${REF.velHash}`);
  assert(Math.abs(energy - REF.energy) < 1e-9, `10-step energy matches pre-S1 reference exactly (Δ = ${Math.abs(energy - REF.energy)})`);
  assert(Math.abs(bindingU - REF.bindingU) < 1e-9, `10-step bindingU matches pre-S1 reference (${bindingU.toFixed(6)})`);
  assert(posHash === REF.posHash, `posHash bit-identical to pre-S1 (${posHash} vs ${REF.posHash})`);
  assert(velHash === REF.velHash, `velHash bit-identical to pre-S1 (${velHash} vs ${REF.velHash})`);

  // explicit charges:false must equal the default (same opt-in plumbing)
  const ffExplicit = new ForceField(sel, { rc: 10, gamma: 1.0, binding: { charges: false } }, ligands);
  let qSame = true;
  for (let i = 0; i < ffExplicit._protQ.length; i++) if (ffExplicit._protQ[i] !== ffDef._protQ[i]) qSame = false;
  assert(qSame, `binding.charges:false is bit-identical to the default (off)`);

  // ---------------------------------------------------------------------
  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log("=================================================");
  if (failed > 0) process.exit(1);
  console.log("PASS: test_charges.js — Loop-2 S1 salt-bridge charges validated");
}

main();
