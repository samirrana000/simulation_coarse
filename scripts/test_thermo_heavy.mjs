/**
 * test_thermo_heavy.mjs — Stage-1: heavy all-atom pocket-entropy sign test (4W52/BNZ).
 *
 * Loop-2 closed with an honest finding (docs/BINDING_PHYSICS_R4.md §5): at Cα-only
 * resolution pocket-ΔS is model-dependent/systematically negative (−TΔS < 0) because
 * the ligand bath rattles pocket Cα while sidechain rotors are invisible. S7
 * hypothesized heavy all-atom mode restores the physical sign (holo < apo entropy,
 * i.e. positive −TΔS_pocket). This script tests that hypothesis headlessly.
 *
 * Design choices (documented per task):
 *  - System built the same way src/main.js buildSystem() builds heavy systems:
 *    parseHeavy + selectHeavy + new HeavyForceField({atoms}, {gamma, weak:"on"}, []).
 *  - Ligand = BNZ only (heteroSelection BNZ:true, EPE:false — same selection as
 *    tests/test_negative.js). EPE is a crystallization buffer; keeping it would
 *    corrupt the ligand-COM pocket definition with a second diffusing body.
 *  - Holo = heavy FF with binding on (charges inherent + weak on), trackTerms on.
 *    Apo = ligand-free clone (protein atoms only, same indices 0..nProt-1), so
 *    bindingU ≡ 0 by construction — the heavy analogue of CG binding:{on:false}
 *    (HeavyForceField has no binding toggle; cf. apoSystemOf() heavy branch in
 *    src/analysis/alanine_scanning.js).
 *  - Protocol per replica: 300-step equilibration (discarded) + 1500 production
 *    steps, stride 2 → 750 frames/leg. LangevinIntegrator reads ff.masses
 *    per-particle (heavy-compatible; dt auto-tunes to 1 fs).
 *  - Pocket = heavy protein atoms within 8 Å of BNZ COM in ref (atom indices,
 *    protein-only so they are valid in both holo and apo frames). Split reported:
 *    backbone (N/CA/C/O) vs sidechain (rest) — the split carries the verdict's
 *    mechanistic content.
 *  - Masses: per-atom ff.masses via computeThermodynamics({masses}) — uniform-110
 *    would mis-weight C/N/O/S in Schlitter's mass-weighted covariance.
 *  - ΔH from the 4 exact protein↔ligand cross terms (lj/coul/hb/desolv); the S3
 *    weak totals (pi/cpi/xb) are whole-system and reported separately for context,
 *    NOT folded into ΔH (they contain protein–protein stacking present in apo too).
 *
 * Verdict rule: the script's own 3-replica mean at 750 frames / 318 DOF is
 * UNDER-SAMPLED by construction, so the script reports its numbers + a Cα control
 * and records the Stage-1 verdict from the deep-sampling pilot (NOT-RESTORED —
 * see the VERDICT block and docs/BINDING_PHYSICS_R4.md §5). The sign itself is
 * never asserted (either outcome exits 0 when the pipeline asserts pass).
 *
 * Run: node scripts/test_thermo_heavy.mjs   (~3-4 min: 6 legs × 1800 × ~20 ms)
 *   SMOKE (Stage-6 fast CI): node scripts/test_thermo_heavy.mjs --smoke
 *   or HEAVY_SMOKE=1 node scripts/test_thermo_heavy.mjs — 1 replica,
 *   50 equil + 300 production steps stride 2 → 150 frames/leg
 *   (0.30 ps production + 0.05 ps equil at heavy dt 0.001 ps = 1 fs),
 *   ~15-25 s. Same asserts (finiteness/boundedness only — sign never
 *   asserted). FULL (default) stays 3 replicas × 300+1500 stride 2.
 */

import { readFileSync } from "node:fs";
import { computeThermodynamics, formatThermoTable } from "../src/analysis/thermodynamics.js";
import { parseHeavy, selectHeavy, HeavyForceField } from "../src/heavy.js";
import { LangevinIntegrator } from "../src/integrator.js";

let passed = 0, failed = 0;
const assert = (c, m) => { if (c) { passed++; console.log(`  ✓ ${m}`); } else { failed++; console.error(`  ✗ FAIL: ${m}`); } };

const SMOKE = process.argv.includes("--smoke") || process.env.HEAVY_SMOKE === "1";
const EQUIL = SMOKE ? 50 : 300, STEPS = SMOKE ? 300 : 1500, STRIDE = 2, REPS = SMOKE ? 1 : 3, T = 300;
const DT_PS = 0.001; // measured heavy dt on this system (integrator._pickDt); 1 fs
const prodPs = STEPS * DT_PS, equilPs = EQUIL * DT_PS;
// Stage-2 determinism: fixed seeds per replica (holo SEEDS[rep], apo SEEDS[rep]+1000
// to decorrelate legs). Seeded via LangevinIntegrator opts.seed (mulberry32);
// replica-mean reporting unchanged (sign never asserted — see VERDICT block).
// Optional override: THERMO_HEAVY_SEED_BASE env shifts all seeds (default 0).
const SEED_BASE = Number.parseInt(process.env.THERMO_HEAVY_SEED_BASE || "0", 10) || 0;
const SEEDS = [1001, 2002, 3003].map((s) => s + SEED_BASE);

// ---- build heavy 4W52/BNZ system (main.js buildSystem heavy path) ----
console.log(`=== Stage-1 setup: 4W52 heavy (BNZ-only ligand) [${SMOKE ? "SMOKE" : "FULL"}] ===`);
console.log(`  protocol: ${REPS} replica(s) × (holo ${EQUIL}+${STEPS} + apo ${EQUIL}+${STEPS} steps, stride ${STRIDE})`);
console.log(`  time: ${prodPs.toFixed(2)} ps production + ${equilPs.toFixed(2)} ps equil per leg (dt ${DT_PS * 1000} fs)`);
const pdbText = readFileSync(new URL("../4w52.pdb", import.meta.url), "utf8");
const parsedHeavy = parseHeavy(pdbText);
const selH = selectHeavy(parsedHeavy, {
  heteroSelection: { "A|200|BNZ": true, "A|201|EPE": false },
  includePdbLigands: true,
});
const nProt = selH.atoms.filter((a) => a.isProtein).length;
const nLig = selH.atoms.filter((a) => a.isLigand).length;
console.log(`  atoms ${selH.atoms.length} = ${nProt} protein + ${nLig} ligand (BNZ)`);
assert(nLig === 6, `BNZ-only ligand selection (${nLig} atoms)`);

function makeHoloFF() {
  const ff = new HeavyForceField({ atoms: selH.atoms }, { gamma: 2.0, weak: "on" }, []);
  ff.trackTerms = true;
  return ff;
}
// Apo: ligand-free clone — protein atom order/indices preserved (ligands are last).
const apoAtoms = selH.atoms.filter((a) => !a.isLigand);
function makeApoFF() {
  return new HeavyForceField({ atoms: apoAtoms }, { gamma: 2.0, weak: "on" }, []);
}
assert(apoAtoms.length === nProt, `apo keeps ${apoAtoms.length} protein atoms`);

// ---- pocket: heavy protein atoms within 8 Å of BNZ COM (ref frame) ----
const probeFF = makeHoloFF();
const ligStart = probeFF.ligandStart;
let lcom = [0, 0, 0];
for (let a = 0; a < probeFF.nLigAtoms; a++) {
  lcom[0] += probeFF.ref[3 * (ligStart + a)] / probeFF.nLigAtoms;
  lcom[1] += probeFF.ref[3 * (ligStart + a) + 1] / probeFF.nLigAtoms;
  lcom[2] += probeFF.ref[3 * (ligStart + a) + 2] / probeFF.nLigAtoms;
}
const BB = new Set(["N", "CA", "C", "O"]);
const pocketIdx = [], bbIdx = [], scIdx = [], caIdx = [];
for (let i = 0; i < probeFF.nProt; i++) {
  const dx = probeFF.ref[3 * i] - lcom[0], dy = probeFF.ref[3 * i + 1] - lcom[1], dz = probeFF.ref[3 * i + 2] - lcom[2];
  if (Math.hypot(dx, dy, dz) < 8.0) {
    pocketIdx.push(i);
    (BB.has(selH.atoms[i].atomName) ? bbIdx : scIdx).push(i);
    if (selH.atoms[i].atomName === "CA") caIdx.push(i);
  }
}
console.log(`  pocket: ${pocketIdx.length} heavy atoms (${bbIdx.length} backbone N/CA/C/O + ${scIdx.length} sidechain, of which ${caIdx.length} Cα), DOF ${pocketIdx.length * 3}`);
assert(pocketIdx.length > 0 && bbIdx.length > 0 && scIdx.length > 0, "pocket + backbone/sidechain subsets nonempty");
assert(caIdx.length >= 5, `Cα subset sizable (${caIdx.length} atoms → ${caIdx.length * 3} DOF, frames/DOF ≈ ${((STEPS / STRIDE) / (caIdx.length * 3)).toFixed(1)})`);
// Apo sanity: bindingU ≡ 0 with no ligand present.
probeFF.compute(probeFF.ref);
const probeApo = makeApoFF();
probeApo.compute(probeApo.ref);
assert(probeApo.bindingU === 0, `apo bindingU ≡ 0 (got ${probeApo.bindingU})`);
console.log(`  holo native U ${probeFF.energy.toFixed(1)} (bindU ${probeFF.bindingU.toFixed(2)}) · apo native U ${probeApo.energy.toFixed(1)}`);

// ---- run one leg: equilibrate (discard) then collect frames + exact cross-term energies ----
/**
 * Run one holo/apo leg deterministically.
 * @param {object} ff force field instance
 * @param {boolean} collectE whether to collect exact cross-term energies + weak means
 * @param {number|null} [seed=null] optional mulberry32 seed (null → Math.random path)
 */
function runLeg(ff, collectE, seed = null) {
  const integ = seed === null || seed === undefined
    ? new LangevinIntegrator(ff.ref, ff, 110.0)
    : new LangevinIntegrator(ff.ref, ff, 110.0, { seed });
  integ.setTemperature(T); integ.setFriction(8.0);
  for (let s = 0; s < EQUIL; s++) integ.step();
  const frames = [], energies = [];
  let wPi = 0, wCpi = 0, wXb = 0, nW = 0;
  for (let s = 0; s < STEPS; s++) {
    integ.step();
    if (s % STRIDE === 0) {
      frames.push(Float32Array.from(integ.pos));
      if (collectE) {
        energies.push([ff.bindLJU, ff.bindCoulU, ff.bindHBU, ff.desolvU, 0, 0, 0]);
        wPi += ff.piU; wCpi += ff.cpiU; wXb += ff.xbU; nW++;
      }
    }
  }
  return { frames, energies, weakMean: nW ? [wPi / nW, wCpi / nW, wXb / nW] : [0, 0, 0] };
}

// ---- replicas (FULL 3, SMOKE 1) ----
console.log(`=== ${REPS} replica(s) × (holo ${EQUIL}+${STEPS} + apo ${EQUIL}+${STEPS} steps, stride ${STRIDE}) ===`);
const t0 = Date.now();
const reps = [];
for (let rep = 0; rep < REPS; rep++) {
  const holo = runLeg(makeHoloFF(), true, SEEDS[rep]);
  const apo = runLeg(makeApoFF(), false, SEEDS[rep] + 1000);
  const masses = makeHoloFF().masses; // per-atom Da (identical protein block in apo)
  const res = computeThermodynamics({
    holoFrames: holo.frames, apoFrames: apo.frames, holoEnergies: holo.energies,
    pocketIdx, nProt, T, masses,
  });
  // Backbone / sidechain / Cα splits (same frames, subset DOFs; no energies → ΔS only).
  // The Cα split is the well-conditioned control (frames/DOF ≥ 10): it tests whether
  // the CG negative sign persists for backbone atoms under heavy dynamics.
  const resBB = computeThermodynamics({ holoFrames: holo.frames, apoFrames: apo.frames, pocketIdx: bbIdx, nProt, T, masses });
  const resSC = computeThermodynamics({ holoFrames: holo.frames, apoFrames: apo.frames, pocketIdx: scIdx, nProt, T, masses });
  const resCA = computeThermodynamics({ holoFrames: holo.frames, apoFrames: apo.frames, pocketIdx: caIdx, nProt, T, masses });
  reps.push({ res, resBB, resSC, resCA, weakMean: holo.weakMean });
  console.log(`  rep${rep}: −TΔS full ${(-T * res.dS.pocket).toFixed(2)} | bb ${(-T * resBB.dS.pocket).toFixed(2)} | sc ${(-T * resSC.dS.pocket).toFixed(2)} | Cα ${(-T * resCA.dS.pocket).toFixed(2)} | ΔH ${res.dH.total.toFixed(2)} kcal/mol`);
  console.log(`         S_holo ${res.S_pocket.holo.toFixed(4)} vs S_apo ${res.S_pocket.apo.toFixed(4)} kcal/mol/K`);
}
console.log(`  wall ${(Date.now() - t0) / 1000 | 0}s`);

// ---- report rep 0 table + means ----
console.log("=== rep0 full table ===");
console.log(formatThermoTable(reps[0].res));
console.log(`  mass model: ${reps[0].res.meta.massModel} · frames/DOF ${reps[0].res.meta.framesPerDof.toFixed(1)}` +
  (reps[0].res.meta.warning ? ` · ⚠ ${reps[0].res.meta.warning}` : ""));
const mean = (f) => reps.reduce((s, r) => s + f(r), 0) / reps.length;
const mTdS = mean((r) => -T * r.res.dS.pocket);
const mTdSbb = mean((r) => -T * r.resBB.dS.pocket);
const mTdSsc = mean((r) => -T * r.resSC.dS.pocket);
const mTdSca = mean((r) => -T * r.resCA.dS.pocket);
const mDH = mean((r) => r.res.dH.total);
const perRun = reps.map((r) => (-T * r.res.dS.pocket).toFixed(2)).join(", ");
console.log(`=== means over ${REPS} replicas ===`);
console.log(`  −TΔS_pocket full [${perRun}] → mean ${mTdS.toFixed(2)} kcal/mol`);
console.log(`  −TΔS backbone mean ${mTdSbb.toFixed(2)} · sidechain mean ${mTdSsc.toFixed(2)} · Cα-control mean ${mTdSca.toFixed(2)} kcal/mol`);
console.log(`  ΔH mean ${mDH.toFixed(2)} (LJ ${mean((r) => r.res.dH.lj).toFixed(2)} / Coul ${mean((r) => r.res.dH.coul).toFixed(2)} / HB ${mean((r) => r.res.dH.hb).toFixed(2)} / desolv ${mean((r) => r.res.dH.desolv).toFixed(2)})`);
console.log(`  weak totals (whole-system, context only) pi ${mean((r) => r.weakMean[0]).toFixed(2)} / cpi ${mean((r) => r.weakMean[1]).toFixed(2)} / xb ${mean((r) => r.weakMean[2]).toFixed(2)}`);

// ---- non-flaky asserts (finiteness + bounds; sign is reported, NOT asserted) ----
console.log("=== asserts ===");
assert(reps.every((r) => Number.isFinite(r.res.dH.total)), "ΔH finite all replicas");
assert(reps.every((r) => Number.isFinite(-T * r.res.dS.pocket)), "−TΔS_pocket finite all replicas");
assert(reps.every((r) => Number.isFinite(-T * r.resBB.dS.pocket) && Number.isFinite(-T * r.resSC.dS.pocket)), "backbone/sidechain −TΔS finite");
assert(reps.every((r) => Number.isFinite(-T * r.resCA.dS.pocket)), "Cα-control −TΔS finite");
assert(Math.abs(mTdS) < 100, `|mean −TΔS| sane < 100 (${mTdS.toFixed(2)})`);
assert(reps.every((r) => r.res.dH_se >= 0 && Number.isFinite(r.res.dH_se)), "bootstrap SE finite");
assert(reps[0].res.meta.massModel === "per-atom", "per-atom mass model active");
assert(reps.every((r) => r.res.meta.holoFrames === STEPS / STRIDE && r.res.meta.apoFrames === STEPS / STRIDE), `frame counts ${STEPS / STRIDE}/leg`);

// ---- verdict (published either way; does not affect exit code) ----
// The 750-frame / 318-DOF protocol is UNDER-SAMPLED by construction (task spec);
// its 3-rep mean is variance-dominated and NOT the verdict by itself. The recorded
// Stage-1 verdict (deep-sampling pilot 2026-09-12, single pair, 500 equil + 6000
// production steps stride 1 → 6000 frames/leg, f/DOF 18.9: full −4.78, bb −5.82,
// sc +1.04, Cα −1.26; block curve full −36.20/−0.95/+9.64/−4.78 at
// 750/1500/3000/6000 frames) is NOT-RESTORED — see docs/BINDING_PHYSICS_R4.md §5.
console.log("=== STAGE-1 VERDICT ===");
console.log(`  this run: full-pocket mean −TΔS = ${mTdS.toFixed(2)} [${perRun}] (UNDER-SAMPLED f/DOF ${reps[0].res.meta.framesPerDof.toFixed(1)} — sign undetermined from this protocol alone)`);
console.log(`  this run: Cα control (well-conditioned) mean −TΔS = ${mTdSca.toFixed(2)} kcal/mol`);
console.log("  recorded Stage-1 verdict: NOT-RESTORED — deep-sampling pilot stays negative (full −4.78); ligand-bath effect persists, sidechain restriction (+1.04) too small to flip the net sign.");

console.log(`\n=== test_thermo_heavy (${SMOKE ? "SMOKE" : "FULL"}): ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed ? 1 : 0);
