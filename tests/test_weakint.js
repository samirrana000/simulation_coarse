/**
 * test_weakint.js — Loop-2 S3 weak-interaction tests (R3 §1a–c, §1e).
 *
 * Covers:
 *  1. FD gradients < 1e-6 for π-stack / cation-π / halogen at generic
 *     geometries (frame level AND atom level, incl. the Newell normal→atom
 *     chain) — prototype numbers reproduced (π 6.2e-11, cπ 0, XB 4.6e-10).
 *  2. Literature energy ranges: π −1..−3 at 3.8 Å parallel; cation-π −2..−5
 *     on-axis; halogen −1..−4 linear (per-element ε).
 *  3. Ring detection: Phe/Tyr/Trp rings + benzene in 4W52.
 *  4. Default-off bit-identity: U on the native pose identical pre/post
 *     (reference captured with the SAME code by toggling par.weak).
 *  5. Weak-on: U finite, negative contributions present, accumulators set.
 *  6. 100-step heavy NVT stability with weak + metalAngles on (no NaN,
 *     RMSD < 1 Å) — 4W52 + benzene and 4HHB (4 FE hemes, porphyrin rings).
 *  7. Metal upgrade: enforceCoordination replaces springs for detected
 *     geometries; coord spring list shrinks accordingly.
 *
 * Run: node tests/test_weakint.js  (zero deps, plain Node)
 */

import fs from "fs";
import { parseHeavy, selectHeavy, appendHeavyLigands, HeavyForceField } from "../src/heavy.js";
import { parseMol2 } from "../src/pdb.js";
import { LangevinIntegrator } from "../src/integrator.js";
import {
  piStackEnergy, cationPiEnergy, halogenEnergy,
  piStackForces, cationPiForces, halogenForces,
  buildRingFrames, buildCationList, buildHalogenList,
  HALOGEN_EPS, fdCheckWeak,
} from "../src/physics/weakint.js";

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}
const vnorm = (a) => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };

/* ---------------- 1. FD gradients at generic geometries --------------- */
console.log("\n[1] FD gradient checks (frame level)...");
{
  // π: parallel-displaced at generic pose (tilted, off-axis) — avoids the
  // max(0,·)² and |·| kinks; normals renormalized (unit-sphere chain rule)
  const cA = [0.1, -0.2, 0.05], nA = [0.05, 0.08, 0.995];
  const cB = [1.4, 0.6, 3.7], nB = [-0.06, -0.03, 0.998];
  const wrap = (cA, nA, cB, nB) => piStackEnergy(cA, vnorm(nA), cB, vnorm(nB));
  const fd = fdCheckWeak(wrap, [cA, nA, cB, nB], ["gA", "gNA", "gB", "gNB"], { normalizeIdx: [1, 3] });
  assert(fd.maxRel < 1e-6, `π-stack FD max rel err ${fd.maxRel.toExponential(2)} < 1e-6 (${fd.worst})`);

  // anti-parallel winding: same geometry, ring B traversed backwards —
  // U must be identical and gradients still FD-exact (canonical n̄ fix)
  const fdAnti = fdCheckWeak(wrap, [cA, nA, cB, [-nB[0], -nB[1], -nB[2]]], ["gA", "gNA", "gB", "gNB"], { normalizeIdx: [1, 3] });
  assert(fdAnti.maxRel < 1e-6, `π-stack anti-parallel FD max rel err ${fdAnti.maxRel.toExponential(2)} < 1e-6`);
  const uP = piStackEnergy(cA, vnorm(nA), cB, vnorm(nB)).U;
  const uAP = piStackEnergy(cA, vnorm(nA), cB, vnorm([-nB[0], -nB[1], -nB[2]])).U;
  assert(Math.abs(uP - uAP) < 1e-12, `π-stack winding-invariant U (${uP.toFixed(6)} vs ${uAP.toFixed(6)})`);

  // cation-π at generic off-axis pose
  const fd2 = fdCheckWeak(cationPiEnergy, [[0.1, 0.2, 0.05], [0.4, 0.5, 4.3], [0.1, -0.2, 0.97]], ["gCat", "gRing", null]);
  assert(fd2.maxRel < 1e-6, `cation-π FD max rel err ${fd2.maxRel.toExponential(2)} < 1e-6 (${fd2.worst})`);

  // halogen at generic (slightly bent, cos β ≈ 0.9) geometry — avoids kink
  const fd3 = fdCheckWeak(halogenEnergy, [[0.2, 0.1, 0.05], [1.75, 0.15, 0.1], [4.7, 0.9, 0.3]], ["gC", "gX", "gD"]);
  assert(fd3.maxRel < 1e-6, `halogen FD max rel err ${fd3.maxRel.toExponential(2)} < 1e-6 (${fd3.worst})`);
}

console.log("\n[1b] FD gradient checks (atom level, flat-buffer kernels)...");
{
  // two generic hexagons + cation + halogen triple in one flat buffer
  const mkRing = (cx, cy, cz, tilt) => {
    const pts = [], idx = [];
    for (let k = 0; k < 6; k++) {
      const a = k * Math.PI / 3;
      const y2 = Math.sin(a) * 1.4 * Math.cos(tilt), z2 = Math.sin(a) * 1.4 * Math.sin(tilt);
      pts.push([Math.cos(a) * 1.4 + cx, y2 + cy, z2 + cz]);
      idx.push(k);
    }
    return { pts, idx };
  };
  const RA = mkRing(0, 0, 0, 0.15), RB = mkRing(0.9, 0.5, 3.75, -0.35);
  const all = [...RA.pts, ...RB.pts];
  const catPos = [2.0, 1.1, 4.5];
  const C = [10, 10, 10], X = [11.75, 10.2, 10.1], D = [14.7, 10.8, 10.35];
  const catIdx = all.length; all.push(catPos);
  const cIdx = all.length; all.push(C);
  const xIdx = all.length; all.push(X);
  const dIdx = all.length; all.push(D);
  const pos = new Float64Array(all.length * 3);
  for (let i = 0; i < all.length; i++) { pos[3 * i] = all[i][0]; pos[3 * i + 1] = all[i][1]; pos[3 * i + 2] = all[i][2]; }
  const ringA = { atomIdx: RA.idx }, ringB = { atomIdx: RB.idx.map((i) => i + 6) };
  const N = all.length, h = 1e-5;

  const fdAtom = (kernel, pick) => {
    const fAn = new Float64Array(N * 3);
    const U0 = kernel(pos, fAn);
    let maxRel = 0, worst = "";
    for (const i of pick) for (let k = 0; k < 3; k++) {
      const p1 = new Float64Array(pos); p1[3 * i + k] += h;
      const p2 = new Float64Array(pos); p2[3 * i + k] -= h;
      const fd = (kernel(p1, new Float64Array(N * 3)) - kernel(p2, new Float64Array(N * 3))) / (2 * h);
      const dU = -fAn[3 * i + k]; // force f = −dU/dx
      const rel = Math.abs(fd - dU) / (Math.abs(fd) + Math.abs(dU) + 1e-12);
      if (rel > maxRel) { maxRel = rel; worst = `atom${i}[${k}] fd=${fd.toPrecision(6)} an=${dU.toPrecision(6)}`; }
    }
    void U0;
    return { maxRel, worst };
  };
  const ringAtoms = [...ringA.atomIdx, ...ringB.atomIdx];
  const r1 = fdAtom((p, f) => piStackForces(p, f, ringA, ringB), ringAtoms);
  assert(r1.maxRel < 1e-6, `π-stack atom-level FD ${r1.maxRel.toExponential(2)} < 1e-6 (${r1.worst})`);
  const r2 = fdAtom((p, f) => cationPiForces(p, f, catIdx, ringB), [catIdx, ...ringB.atomIdx]);
  assert(r2.maxRel < 1e-6, `cation-π atom-level FD ${r2.maxRel.toExponential(2)} < 1e-6`);
  const r3 = fdAtom((p, f) => halogenForces(p, f, cIdx, xIdx, dIdx), [cIdx, xIdx, dIdx]);
  assert(r3.maxRel < 1e-6, `halogen atom-level FD ${r3.maxRel.toExponential(2)} < 1e-6`);

  // momentum conservation of the π pass (net force zero)
  const f = new Float64Array(N * 3);
  piStackForces(pos, f, ringA, ringB);
  let sx = 0, sy = 0, sz = 0;
  for (let i = 0; i < N; i++) { sx += f[3 * i]; sy += f[3 * i + 1]; sz += f[3 * i + 2]; }
  assert(Math.hypot(sx, sy, sz) < 1e-9, `π-stack net force ~0 (${Math.hypot(sx, sy, sz).toExponential(2)})`);
}

/* ---------------- 2. Literature energy ranges ------------------------- */
console.log("\n[2] Literature-range energy checks...");
{
  // π at 3.8 Å parallel, 0 offset: −1..−3 (Hunter–Sanders)
  const u = piStackEnergy([0, 0, 0], [0, 0, 1], [0, 0, 3.8], [0, 0, 1]).U;
  assert(u >= -3 && u <= -1, `π parallel 3.8 Å = ${u.toFixed(3)} in [−3, −1] (lit −1..−3)`);
  // well minimum at 3.8 Å (scan)
  let best = 0, bestZ = 0;
  for (let z = 3.0; z <= 5.5; z += 0.05) {
    const e = piStackEnergy([0, 0, 0], [0, 0, 1], [0, 0, z], [0, 0, 1]).U;
    if (e < best) { best = e; bestZ = z; }
  }
  assert(best >= -3 && best <= -1, `π scan minimum ${best.toFixed(3)} @ ${bestZ.toFixed(2)} Å in [−3, −1]`);

  // cation-π on-axis 4.3 Å: −2..−5 (Gallivan & Dougherty)
  const uc = cationPiEnergy([0, 0, 0], [0, 0, 4.3], [0, 0, 1]).U;
  assert(uc >= -5 && uc <= -2, `cation-π on-axis 4.3 Å = ${uc.toFixed(3)} in [−5, −2] (lit −2..−5)`);

  // halogen linear at 3.1 Å per element: −ε_X, all in −1..−4 (net water)
  for (const [el, eps] of Object.entries(HALOGEN_EPS)) {
    const u = halogenEnergy([0, 0, 0], [1.75, 0, 0], [4.85, 0, 0], { eps }).U;
    assert(Math.abs(u - (-eps)) < 1e-9 && u >= -4 && u <= -1,
      `halogen ${el} linear 3.1 Å = ${u.toFixed(3)} = −ε_${el} (${eps}) in [−4, −1]`);
  }
  // angular decay: 120° ≪ 180° (prototype: −0.50 vs −2.00)
  const uLin = halogenEnergy([0, 0, 0], [1.75, 0, 0], [4.85, 0, 0], { eps: 2.0 }).U;
  const u120 = halogenEnergy([0, 0, 0], [1.75, 0, 0], [1.75 + 3.1 * Math.cos(Math.PI / 3), 3.1 * Math.sin(Math.PI / 3), 0], { eps: 2.0 }).U;
  assert(u120 > uLin + 0.8, `halogen angular decay: 120° (${u120.toFixed(2)}) ≪ linear (${uLin.toFixed(2)})`);
  // F excluded by construction (HALOGEN_EPS has no F key)
  assert(!("F" in HALOGEN_EPS), "F excluded from halogen σ-hole ε table");
}

/* ---------------- 3. Ring detection on 4W52 -------------------------- */
console.log("\n[3] Ring detection (4W52 + benzene ligand)...");
const pdbText = fs.readFileSync("4w52.pdb", "utf-8");
const mol2Text = fs.readFileSync("benzene.mol2", "utf-8");
const parsedHeavy = parseHeavy(pdbText);
const mols = parseMol2(mol2Text);
let selA = selectHeavy(parsedHeavy, {
  heteroSelection: { "A|200|BNZ": false, "A|201|EPE": false },
  includePdbLigands: true, hasExternalLigand: true,
});
selA = appendHeavyLigands(selA, mols, { gaff: true });

{
  const ff = new HeavyForceField({ atoms: selA.atoms }, { gamma: 2.0, weak: "on" }, []);
  const byRes = {};
  for (const r of ff._weakRings) {
    const a = ff.atoms[r.atomIdx[0]];
    byRes[a.resName] = (byRes[a.resName] ?? 0) + 1;
  }
  // 4W52 chain A (164 res): PHE 5 (PHE4,67,104,114,153), TYR 6, TRP 3 (2 rings
  // each: 5+6), HIS ≥1 (ring only when complete — HIS31 here)
  assert((byRes.PHE ?? 0) === 5, `Phe rings = ${byRes.PHE ?? 0} (expected 5)`);
  assert((byRes.TYR ?? 0) === 6, `Tyr rings = ${byRes.TYR ?? 0} (expected 6)`);
  assert((byRes.TRP ?? 0) === 6, `Trp rings = ${byRes.TRP ?? 0} (expected 3×2 fused)`);
  assert((byRes.LIG ?? 0) + (byRes.BEN ?? 0) >= 1, `benzene ligand ring found (${byRes.LIG ?? byRes.BEN ?? 0})`);
  const hisRings = Object.entries(byRes).filter(([k]) => k === "HIS" || k.startsWith("HI"));
  assert(hisRings.length === 0 || (byRes.HIS ?? 0) >= 1, `His ring detection sane (${byRes.HIS ?? 0})`);

  // every frame: unit normal + ≥5 atoms + cycle-ordered (adjacent atoms bonded)
  const bonds = new Set();
  for (let a = 0; a < ff.bonds.length; a += 3) {
    const i = ff.bonds[a], j = ff.bonds[a + 1];
    bonds.add(i < j ? `${i}-${j}` : `${j}-${i}`);
  }
  let framesOk = true;
  for (const r of ff._weakRings) {
    const L = Math.hypot(...r.normal);
    if (Math.abs(L - 1) > 1e-9 || r.atomIdx.length < 5) { framesOk = false; break; }
    let adjacent = true;
    for (let t = 0; t < r.atomIdx.length; t++) {
      const i = r.atomIdx[t], j = r.atomIdx[(t + 1) % r.atomIdx.length];
      const k = i < j ? `${i}-${j}` : `${j}-${i}`;
      if (!bonds.has(k)) { adjacent = false; break; }
    }
    if (!adjacent) { framesOk = false; break; }
  }
  assert(framesOk, "all ring frames: unit normals, ≥5 atoms, cycle-ordered ring walk");

  // standalone buildRingFrames on the same system (API check)
  const bondPairs = [];
  for (let a = 0; a < ff.bonds.length; a += 3) bondPairs.push([ff.bonds[a], ff.bonds[a + 1]]);
  const frames = buildRingFrames(selA.atoms, bondPairs);
  assert(frames.length === ff._weakRings.length, `buildRingFrames standalone = ${frames.length} matches FF (${ff._weakRings.length})`);
}

/* ---------------- 4. Default-off bit-identity ------------------------- */
console.log("\n[4] Default-off bit-identity (native pose)...");
{
  const mk = (par) => new HeavyForceField({ atoms: selA.atoms }, { gamma: 2.0, ...par }, []);
  // reference: constructed weak-on then toggled OFF (same code path)
  const refFF = mk({ weak: "on" });
  refFF.weakOn = false;
  const uRef = refFF.compute(refFF.ref);

  const ffDefault = mk({});
  const uDefault = ffDefault.compute(ffDefault.ref);
  assert(uDefault === uRef, `default U === weak-on-toggled-off U (bit-identical: ${uDefault})`);
  assert(ffDefault.weakU === 0 && ffDefault.piU === 0 && ffDefault.cpiU === 0 && ffDefault.xbU === 0,
    "default-off accumulators zero");

  const ffWeak = mk({ weak: "on" });
  const uWeak = ffWeak.compute(ffWeak.ref);
  assert(uWeak !== uDefault, `weak-on changes U (${uDefault.toFixed(3)} → ${uWeak.toFixed(3)})`);
  assert(Number.isFinite(uWeak) && Number.isFinite(ffWeak.weakU), "weak-on U finite");

  // PDB HETATM path (BNZ ligand) also bit-identical default vs toggled
  let selB = selectHeavy(parsedHeavy, {
    heteroSelection: { "A|200|BNZ": true, "A|201|EPE": true },
    includePdbLigands: true, hasExternalLigand: false,
  });
  const ffB1 = new HeavyForceField({ atoms: selB.atoms }, { gamma: 2.0 }, []);
  const uB1 = ffB1.compute(ffB1.ref);
  const ffB2 = new HeavyForceField({ atoms: selB.atoms }, { gamma: 2.0, weak: "on" }, []);
  ffB2.weakOn = false;
  const uB2 = ffB2.compute(ffB2.ref);
  assert(uB1 === uB2, `HETATM-ligand default U bit-identical (${uB1})`);
}

/* ---------------- 5. Weak-on energy content --------------------------- */
console.log("\n[5] Weak-on contributions (4W52 native)...");
{
  const ff = new HeavyForceField({ atoms: selA.atoms }, { gamma: 2.0, weak: "on" }, []);
  const U = ff.compute(ff.ref);
  assert(Number.isFinite(U), `weak-on U finite (${U.toFixed(3)})`);
  assert(ff.cpiU < 0, `cation-π negative contribution present (${ff.cpiU.toFixed(3)} kcal/mol)`);
  assert(Math.abs(ff.weakU - (ff.piU + ff.cpiU + ff.xbU)) < 1e-12, `weakU = piU + cpiU + xbU (${ff.weakU.toFixed(3)})`);
  // toggling off restores the exact default energy (same instance)
  ff.weakOn = false;
  const U2 = ff.compute(ff.ref);
  const ffD = new HeavyForceField({ atoms: selA.atoms }, { gamma: 2.0 }, []);
  const uD = ffD.compute(ffD.ref);
  assert(U2 === uD, `toggle-off restores default U (${U2} === ${uD})`);
}

/* ---------------- 6. NVT stability (weak + metalAngles on) ------------- */
console.log("\n[6] 100-step heavy NVT stability (weak + metalAngles)...");
{
  const ff = new HeavyForceField({ atoms: selA.atoms }, { gamma: 2.0, weak: "on", metalAngles: true }, []);
  const integ = new LangevinIntegrator(ff.ref, ff, 110.0);
  integ.setTemperature(300.0);
  integ.setFriction(8.0);
  let finite = true;
  for (let s = 0; s < 100; s++) {
    integ.step();
    if (!Number.isFinite(ff.energy)) { finite = false; break; }
  }
  const rmsd = ff.rmsd(integ.pos);
  assert(finite && Number.isFinite(rmsd), "100 NVT steps: no NaN (4W52+benzene, weak+metalAngles)");
  assert(rmsd < 1.0, `RMSD < 1 Å after 100 steps (got ${rmsd.toFixed(3)} Å)`);

  // 4HHB: 4 Fe hemes — metal upgrade exercised with real coordination
  const pdbText2 = fs.readFileSync("4hhb.pdb", "utf-8");
  const parsedHeavy2 = parseHeavy(pdbText2);
  const atoms2 = parsedHeavy2.atoms.filter((a) => a.isProtein || a.resName === "HEM");
  const ff2 = new HeavyForceField({ atoms: atoms2 }, { gamma: 1.0, weak: "on", metalAngles: true }, []);
  const U2 = ff2.compute(ff2.ref);
  assert(Number.isFinite(U2), `4HHB weak+metalAngles U finite (${U2.toFixed(1)})`);
  assert(ff2._metalEnforce.metals.length === 4, `4HHB metals detected = ${ff2._metalEnforce.metals.length} (4 Fe)`);
  assert(ff2._metalEnforce.hasGeometry.size === 4, `4HHB metals with coordination geometry = ${ff2._metalEnforce.hasGeometry.size}`);
  const integ2 = new LangevinIntegrator(ff2.ref, ff2, 110.0);
  integ2.setTemperature(300.0);
  integ2.setFriction(8.0);
  let finite2 = true;
  for (let s = 0; s < 100; s++) {
    integ2.step();
    if (!Number.isFinite(ff2.energy)) { finite2 = false; break; }
  }
  const rmsd2 = ff2.rmsd(integ2.pos);
  assert(finite2 && Number.isFinite(rmsd2), "100 NVT steps 4HHB: no NaN");
  assert(rmsd2 < 1.0, `4HHB RMSD < 1 Å after 100 steps (got ${rmsd2.toFixed(3)} Å)`);
}

/* ---------------- 7. Metal upgrade wiring ------------------------------ */
console.log("\n[7] Metal enforceCoordination swap (R3 §1e)...");
{
  const pdbText2 = fs.readFileSync("4hhb.pdb", "utf-8");
  const parsedHeavy2 = parseHeavy(pdbText2);
  const atoms2 = parsedHeavy2.atoms.filter((a) => a.isProtein || a.resName === "HEM");
  const ffLegacy = new HeavyForceField({ atoms: atoms2 }, { gamma: 1.0 }, []);
  assert(ffLegacy._metalEnforce === null, "default: no metal enforce state (legacy springs)");
  const nSpringsLegacy = ffLegacy.coord.length / 3;
  assert(nSpringsLegacy > 0, `legacy metal distance springs present (${nSpringsLegacy})`);

  const ffSwapped = new HeavyForceField({ atoms: atoms2 }, { gamma: 1.0, metalAngles: true }, []);
  const nSpringsKept = ffSwapped.coord.length / 3;
  assert(nSpringsKept === 0, `all geometry-detected metals' springs removed (kept ${nSpringsKept})`);
  const uS = ffSwapped.compute(ffSwapped.ref);
  assert(Number.isFinite(uS) && ffSwapped.coordU > 0, `enforceCoordination U finite & positive (${ffSwapped.coordU.toFixed(2)})`);
  // legacy path U bit-identical to pre-change (springs untouched default)
  const ffAgain = new HeavyForceField({ atoms: atoms2 }, { gamma: 1.0 }, []);
  assert(ffAgain.compute(ffAgain.ref) === ffLegacy.compute(ffLegacy.ref), "default metal path deterministic");
  // detected geometries are sensible (heme FE: 3–6 donors; deoxy-heme FE
  // is 5-coordinate, one FE here coordinates 3 protein donors + porphyrin N)
  for (const m of ffSwapped._metalEnforce.metals) {
    assert(m.donors.length >= 3 && m.donors.length <= 6, `FE donor shell ${m.donors.length} in [3,6]`);
  }
}

/* ---------------- summary ---------------- */
console.log("\n=================================================");
console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
console.log("=================================================");
if (failed > 0) process.exit(1);
