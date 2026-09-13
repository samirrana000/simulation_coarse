/**
 * test_pocket_entropy.mjs — Stage-1 follow-up: deep-sampling pocket-entropy protocol (4W52/BNZ heavy).
 *
 * Prior routine closed (commit b3cf17d): pocket-ΔS sign NOT-restored — 750f/leg
 * variance-dominated (span −13.7…+59.9); at 6000f depth: full −4.78, backbone
 * −5.82, sidechain +1.04, Cα −1.26. Prescription: (a) 10+ ps × ≥3 replicas,
 * (b) explicit pocket-χ torsion-Shannon term (~100 frames/rotor),
 * (c) sidechain-only pocket entropy.
 *
 * Protocol (documented per task):
 *  - System built exactly like scripts/test_thermo_heavy.mjs: parseHeavy +
 *    selectHeavy + HeavyForceField({gamma:2.0, weak:"on"}), BNZ-only ligand
 *    (EPE excluded — same heteroSelection as tests/test_negative.js).
 *  - Holo = binding on + trackTerms; apo = ligand-free clone (protein indices
 *    0..nProt-1 preserved), bindingU ≡ 0 by construction.
 *  - Integrator dt: heavy auto-tunes to dt = 0.001 ps = 1 fs (see
 *    src/integrator.js:_pickDt — heavy maxDt 0.001; measured 0.001 on this
 *    system). Hence 10 ps production = 10000 steps. Equilibration is extra
 *    (discarded) and documented in ps below.
 *  - PILOT (default, fast): EQUIL 50 + STEPS 600, stride 2 → 300 frames/leg
 *    (0.60 ps production + 0.05 ps equil), 1 replica. Exercises every code
 *    path (full/bb/sc/Cα + χ + blocks) in ~35 s.
 *  - LONG (--long / --slow, opt-in, never default): EQUIL 500 + STEPS 10000,
 *    stride 2 → 5000 frames/leg (10.0 ps production + 0.5 ps equil),
 *    3 seeded replicas (~27 min at ~26 ms/step + ~2-4 min analysis).
 *  - Seeds fixed per replica (holo SEEDS[rep], apo SEEDS[rep]+1000;
 *    THERMO_CHI_SEED_BASE env override) — two pilot runs are bit-identical.
 *  - Pocket = heavy protein atoms within 8 Å of BNZ COM (same rule);
 *    backbone (N/CA/C/O) vs sidechain (CB outward) via
 *    splitPocketByBackbone(); Cα control alongside.
 *  - χ torsions: pocketChiTorsions() — residues with ≥1 sidechain atom in the
 *    pocket; ALA (CB but no rotor) skipped, GLY none; all quadruplets resolved
 *    against the full atom table (protein-absolute ⇒ valid holo + apo).
 *    On 4W52: 25 pocket residues → 16 sidechain residues → 14 χ residues
 *    (2 ALA skipped, 0 unresolved) → 26 χ torsions (14 χ1 + 11 χ2 + 1 χ3);
 *    LONG frames/rotor = 5000/26 ≈ 192 (≥100 target).
 *  - Block-convergence curve slices the FIRST N production frames
 *    (contiguous, deterministic): LONG [750,1500,3000,5000], PILOT [150,300].
 *  - Sign is reported, NEVER asserted (either outcome exits 0 when the
 *    pipeline asserts pass) — no forced sign.
 *
 * Run: node scripts/test_pocket_entropy.mjs            (pilot, ~35 s)
 *      node scripts/test_pocket_entropy.mjs --long     (full, ~30 min)
 */

import { readFileSync } from "node:fs";
import {
  computeThermodynamics, formatThermoTable, pocketChiTorsions,
  splitPocketByBackbone, chiEntropyDelta,
} from "../src/analysis/thermodynamics.js";
import { parseHeavy, selectHeavy, HeavyForceField } from "../src/heavy.js";
import { LangevinIntegrator } from "../src/integrator.js";

let passed = 0, failed = 0;
const assert = (c, m) => { if (c) { passed++; console.log(`  ✓ ${m}`); } else { failed++; console.error(`  ✗ FAIL: ${m}`); } };

const LONG = process.argv.includes("--long") || process.argv.includes("--slow");
const EQUIL = LONG ? 500 : 50, STEPS = LONG ? 10000 : 600, STRIDE = 2, REPS = LONG ? 3 : 1, T = 300;
const BLOCKS = LONG ? [750, 1500, 3000, 5000] : [150, 300];
const DT_PS = 0.001; // measured heavy dt on this system (integrator._pickDt); 1 fs
const prodPs = STEPS * DT_PS, equilPs = EQUIL * DT_PS;
// Stage-2 determinism pattern (same seeds as test_thermo_heavy for comparability).
const SEED_BASE = Number.parseInt(process.env.THERMO_CHI_SEED_BASE || process.env.THERMO_HEAVY_SEED_BASE || "0", 10) || 0;
const SEEDS = [1001, 2002, 3003].map((s) => s + SEED_BASE);

// ---- build heavy 4W52/BNZ system (same as test_thermo_heavy) ----
console.log(`=== Stage-1 pocket-entropy ${LONG ? "LONG" : "PILOT"}: 4W52 heavy (BNZ-only) ===`);
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
const apoAtoms = selH.atoms.filter((a) => !a.isLigand);
function makeApoFF() {
  return new HeavyForceField({ atoms: apoAtoms }, { gamma: 2.0, weak: "on" }, []);
}
assert(apoAtoms.length === nProt, `apo keeps ${apoAtoms.length} protein atoms`);

// ---- pocket + splits + χ ----
const probeFF = makeHoloFF();
const ligStart = probeFF.ligandStart;
let lcom = [0, 0, 0];
for (let a = 0; a < probeFF.nLigAtoms; a++) {
  lcom[0] += probeFF.ref[3 * (ligStart + a)] / probeFF.nLigAtoms;
  lcom[1] += probeFF.ref[3 * (ligStart + a) + 1] / probeFF.nLigAtoms;
  lcom[2] += probeFF.ref[3 * (ligStart + a) + 2] / probeFF.nLigAtoms;
}
const pocketIdx = [];
for (let i = 0; i < probeFF.nProt; i++) {
  const dx = probeFF.ref[3 * i] - lcom[0], dy = probeFF.ref[3 * i + 1] - lcom[1], dz = probeFF.ref[3 * i + 2] - lcom[2];
  if (Math.hypot(dx, dy, dz) < 8.0) pocketIdx.push(i);
}
const { bb: bbIdx, sc: scIdx } = splitPocketByBackbone(selH.atoms, pocketIdx);
const caIdx = pocketIdx.filter((i) => selH.atoms[i].atomName === "CA");
console.log(`  pocket: ${pocketIdx.length} heavy atoms (${bbIdx.length} backbone + ${scIdx.length} sidechain CB-outward, of which ${caIdx.length} Cα), DOF ${pocketIdx.length * 3}`);
assert(pocketIdx.length > 0 && bbIdx.length > 0 && scIdx.length > 0, "pocket + backbone/sidechain subsets nonempty");
assert(caIdx.length >= 5, `Cα subset sizable (${caIdx.length} atoms)`);
const chiBuilt = pocketChiTorsions(selH.atoms, pocketIdx);
const chiTorsions = chiBuilt.torsions;
console.log(`  χ: ${chiTorsions.length} torsions over ${chiBuilt.coverage.nChiResidues} residues (${chiBuilt.coverage.nSidechainResidues} sidechain residues in ${chiBuilt.coverage.nPocketResidues} pocket residues; ALA skipped ${chiBuilt.coverage.alaSkipped}, unresolved ${chiBuilt.coverage.unresolved})`);
console.log(`     ${chiBuilt.details.map((d) => `${d.resKey.split("|")[1]}${d.resKey.split("|")[0].split("|")[1]}:${d.label}`).join(" ")}`);
assert(chiTorsions.length >= 20, `χ coverage sizable (${chiTorsions.length} torsions, ALA-free CB+ residues)`);
assert(chiBuilt.coverage.unresolved === 0, "all χ quadruplets resolved in full atom table");
probeFF.compute(probeFF.ref);
const probeApo = makeApoFF();
probeApo.compute(probeApo.ref);
assert(probeApo.bindingU === 0, `apo bindingU ≡ 0 (got ${probeApo.bindingU})`);

/**
 * Run one holo/apo leg deterministically.
 * @param {object} ff force field instance
 * @param {boolean} collectE whether to collect exact cross-term energies
 * @param {number|null} [seed=null] mulberry32 seed (null → Math.random path)
 */
function runLeg(ff, collectE, seed = null) {
  const integ = seed === null || seed === undefined
    ? new LangevinIntegrator(ff.ref, ff, 110.0)
    : new LangevinIntegrator(ff.ref, ff, 110.0, { seed });
  integ.setTemperature(T); integ.setFriction(8.0);
  for (let s = 0; s < EQUIL; s++) integ.step();
  const frames = [], energies = [];
  for (let s = 0; s < STEPS; s++) {
    integ.step();
    if (s % STRIDE === 0) {
      frames.push(Float32Array.from(integ.pos));
      if (collectE) energies.push([ff.bindLJU, ff.bindCoulU, ff.bindHBU, ff.desolvU, 0, 0, 0]);
    }
  }
  return { frames, energies };
}

/** ΔS-only split result on a frame slice (first N frames, contiguous). */
function splitRes(holoFrames, apoFrames, idx, masses, N, withChi) {
  return computeThermodynamics({
    holoFrames: holoFrames.slice(0, N), apoFrames: apoFrames.slice(0, N),
    pocketIdx: idx, nProt, T, masses,
    ...(withChi ? { chiTorsions } : {}),
  });
}

// ---- replicas ----
console.log(`=== ${REPS} replica(s) ===`);
const t0 = Date.now();
const reps = [];
const FRAMES_PER_LEG = STEPS / STRIDE;
for (let rep = 0; rep < REPS; rep++) {
  const holo = runLeg(makeHoloFF(), true, SEEDS[rep]);
  const apo = runLeg(makeApoFF(), false, SEEDS[rep] + 1000);
  const masses = makeHoloFF().masses;
  const res = computeThermodynamics({
    holoFrames: holo.frames, apoFrames: apo.frames, holoEnergies: holo.energies,
    pocketIdx, nProt, T, masses, chiTorsions,
  });
  const resBB = computeThermodynamics({ holoFrames: holo.frames, apoFrames: apo.frames, pocketIdx: bbIdx, nProt, T, masses });
  const resSC = computeThermodynamics({ holoFrames: holo.frames, apoFrames: apo.frames, pocketIdx: scIdx, nProt, T, masses });
  const resCA = computeThermodynamics({ holoFrames: holo.frames, apoFrames: apo.frames, pocketIdx: caIdx, nProt, T, masses });
  reps.push({ res, resBB, resSC, resCA, holoFrames: holo.frames, apoFrames: apo.frames, masses });
  const chi = res.chi;
  console.log(`  rep${rep}: −TΔS full ${(-T * res.dS.pocket).toFixed(2)} | bb ${(-T * resBB.dS.pocket).toFixed(2)} | sc ${(-T * resSC.dS.pocket).toFixed(2)} | Cα ${(-T * resCA.dS.pocket).toFixed(2)} | χ −TΔS ${(-T * chi.dS).toFixed(2)} ± ${(T * chi.se).toFixed(2)} | ΔH ${res.dH.total.toFixed(2)}`);
}
console.log(`  wall ${(Date.now() - t0) / 1000 | 0}s`);

// ---- per-replica tables ----
console.log("=== rep0 full table (with χ line) ===");
console.log(formatThermoTable(reps[0].res));
console.log(`  mass model: ${reps[0].res.meta.massModel} · frames/DOF ${reps[0].res.meta.framesPerDof.toFixed(1)} · χ ${reps[0].res.chi.nChi} rotors × ${reps[0].res.chi.framesPerRotor.toFixed(0)} frames/rotor` +
  (reps[0].res.meta.warning ? ` · ⚠ ${reps[0].res.meta.warning}` : ""));
const mean = (f) => reps.reduce((s, r) => s + f(r), 0) / reps.length;
const sd = (f) => {
  if (reps.length < 2) return 0;
  const m = mean(f);
  return Math.sqrt(reps.reduce((s, r) => s + (f(r) - m) ** 2, 0) / (reps.length - 1));
};
const mTdS = mean((r) => -T * r.res.dS.pocket);
const mTdSbb = mean((r) => -T * r.resBB.dS.pocket);
const mTdSsc = mean((r) => -T * r.resSC.dS.pocket);
const mTdSca = mean((r) => -T * r.resCA.dS.pocket);
const mTdSchi = mean((r) => -T * r.res.chi.dS);
const mDH = mean((r) => r.res.dH.total);
console.log(`=== means over ${REPS} replica(s) ===`);
console.log(`  −TΔS full ${mTdS.toFixed(2)} (SD ${sd((r) => -T * r.res.dS.pocket).toFixed(2)}) · bb ${mTdSbb.toFixed(2)} · sc ${mTdSsc.toFixed(2)} · Cα ${mTdSca.toFixed(2)} · χ ${mTdSchi.toFixed(2)} kcal/mol`);
console.log(`  ΔH mean ${mDH.toFixed(2)} (LJ ${mean((r) => r.res.dH.lj).toFixed(2)} / Coul ${mean((r) => r.res.dH.coul).toFixed(2)} / HB ${mean((r) => r.res.dH.hb).toFixed(2)} / desolv ${mean((r) => r.res.dH.desolv).toFixed(2)})`);

// ---- block-convergence curve (replica-mean −TΔS at first-N frames) ----
console.log(`=== block convergence (first-N frames; f/DOF full @N/${pocketIdx.length * 3} DOF, χ fr/rot @N/${chiTorsions.length}) ===`);
const blockRows = [];
for (const N of BLOCKS) {
  if (N > FRAMES_PER_LEG) continue;
  const row = { N };
  for (const [key, idx, withChi] of [["full", pocketIdx, true], ["bb", bbIdx, false], ["sc", scIdx, false], ["ca", caIdx, false]]) {
    const vals = reps.map((r) => -T * splitRes(r.holoFrames, r.apoFrames, idx, r.masses, N, withChi).dS.pocket);
    row[key] = vals.reduce((s, v) => s + v, 0) / vals.length;
  }
  const chiVals = reps.map((r) => -T * chiEntropyDelta(r.holoFrames.slice(0, N), r.apoFrames.slice(0, N), chiTorsions).dS);
  row.chi = chiVals.reduce((s, v) => s + v, 0) / chiVals.length;
  row.fDof = (N / (pocketIdx.length * 3)).toFixed(1);
  row.frRot = (N / chiTorsions.length).toFixed(0);
  blockRows.push(row);
  console.log(`  N=${N} (f/DOF ${row.fDof}, fr/rot ${row.frRot}): full ${row.full.toFixed(2)} | bb ${row.bb.toFixed(2)} | sc ${row.sc.toFixed(2)} | ca ${row.ca.toFixed(2)} | χ ${row.chi.toFixed(2)}`);
}

// ---- non-flaky asserts (finiteness + coverage; sign reported, NOT asserted) ----
console.log("=== asserts ===");
assert(reps.every((r) => r.holoFrames.length === FRAMES_PER_LEG && r.apoFrames.length === FRAMES_PER_LEG), `frame counts ${FRAMES_PER_LEG}/leg`);
assert(reps.every((r) => Number.isFinite(r.res.dH.total)), "ΔH finite all replicas");
assert(reps.every((r) => Number.isFinite(-T * r.res.dS.pocket)), "−TΔS full finite all replicas");
assert(reps.every((r) => Number.isFinite(-T * r.resBB.dS.pocket) && Number.isFinite(-T * r.resSC.dS.pocket)), "backbone/sidechain-only −TΔS finite");
assert(reps.every((r) => Number.isFinite(-T * r.resCA.dS.pocket)), "Cα −TΔS finite");
assert(reps.every((r) => Number.isFinite(r.res.chi.dS) && Number.isFinite(r.res.chi.se) && r.res.chi.framesPerRotor > 0), "χ ΔS + jackknife SE finite, frames/rotor > 0");
assert(reps[0].res.meta.massModel === "per-atom", "per-atom mass model active");
assert(reps.every((r) => r.res.dH_se >= 0 && Number.isFinite(r.res.dH_se)), "bootstrap SE finite");
assert(blockRows.length === BLOCKS.filter((N) => N <= FRAMES_PER_LEG).length && blockRows.every((row) => ["full", "bb", "sc", "ca", "chi"].every((k) => Number.isFinite(row[k]))), "block curve finite all N × splits");

// ---- verdict (published either way; does not affect exit code) ----
console.log("=== STAGE-1 VERDICT ===");
{
  const repSDsc = sd((r) => -T * r.resSC.dS.pocket);
  const lastTwo = blockRows.slice(-2);
  const scStablePos = lastTwo.length === 2 && lastTwo.every((row) => row.sc > 0);
  const chiStablePos = lastTwo.length === 2 && lastTwo.every((row) => row.chi > 0);
  const chiSE = mean((r) => T * r.res.chi.se);
  const samplingOK = reps[0].res.meta.framesPerDof >= 10 && reps[0].res.chi.framesPerRotor >= 100;
  console.log(`  sampling: frames/DOF ${reps[0].res.meta.framesPerDof.toFixed(1)} (need ≥10) · χ frames/rotor ${reps[0].res.chi.framesPerRotor.toFixed(0)} (target ~100) → ${samplingOK ? "ADEQUATE" : "UNDER-SAMPLED"}`);
  console.log(`  full-depth means: sc ${mTdSsc.toFixed(2)} (repSD ${repSDsc.toFixed(2)}) · χ ${mTdSchi.toFixed(2)} (jackSE ~${chiSE.toFixed(2)}) · bb ${mTdSbb.toFixed(2)} · full ${mTdS.toFixed(2)}`);
  console.log(`  block stability (last two N): sidechain ${scStablePos ? "POSITIVE×2" : "not-stable-positive"} · χ ${chiStablePos ? "POSITIVE×2" : "not-stable-positive"}`);
  const separable = mTdSsc > 0 && mTdSchi > 0 && (mTdSsc - mTdSbb) > 0 && scStablePos && chiStablePos;
  if (samplingOK && separable) {
    console.log(`  verdict: YES — converged positive-restriction signal (sidechain/χ) separable from backbone rattle: sc ${mTdSsc.toFixed(2)} vs bb ${mTdSbb.toFixed(2)} (Δ ${(mTdSsc - mTdSbb).toFixed(2)}), χ ${mTdSchi.toFixed(2)} ± ${chiSE.toFixed(2)}.`);
  } else if (!samplingOK) {
    console.log("  verdict: NO — sampling floor: Schlitter and/or χ below their conditioning targets; deepen sampling before judging the sign.");
  } else {
    console.log(`  verdict: NO — no converged separable restriction: net full ${mTdS.toFixed(2)} with backbone ${mTdSbb.toFixed(2)} vs sidechain ${mTdSsc.toFixed(2)}, χ ${mTdSchi.toFixed(2)}; ` +
      (Math.abs(mTdSsc) < repSDsc ? "sidechain signal within replica noise (sampling floor dominates). " : "sidechain resolved but overwhelmed by backbone bath (force-field ceiling at 8 Å/ps scale). ") +
      "See §15.");
  }
}

console.log(`\n=== test_pocket_entropy (${LONG ? "LONG" : "PILOT"}): ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed ? 1 : 0);
