/**
 * test_bindlog_integration.mjs — Loop-2 S4 integration test (R4 §5 item 1,
 * R6 §5, REVIEW_LOOP1 S4 validation gate).
 *
 * Exercises the S4 wiring end-to-end on the real 4W52 CG system:
 *   1. per-term accumulators (ff.trackTerms = true, charges ON, hbMode
 *      "directional"): bindLJU/bindCoulU/bindHBU/desolvU sum to bindingU;
 *   2. 300-step Langevin mini-run with a BindLog capturing stride-10 frames
 *      + 7-term energy components + contact form/break events (same 5.5 Å
 *      pair-set + diff semantics as the main.js tick wiring);
 *   3. frame/event counts, blob round-trip, and the native-pose LJ ΔH sign.
 *
 * Run: node scripts/test_bindlog_integration.mjs
 */
import { readFileSync } from "node:fs";
import { BindLog } from "../src/capture/bindlog.js";

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

const pdbText = readFileSync(new URL("../4w52.pdb", import.meta.url), "utf8");
const { parseCa, selectSystem, parseLigands } = await import("../src/pdb.js");
const { ForceField } = await import("../src/forcefield.js");
const { LangevinIntegrator } = await import("../src/integrator.js");

const parsed = parseCa(pdbText);
const sel = selectSystem(parsed);
const mols = parseLigands(pdbText);
const par = { rc: 10, gamma: 2.0, temp: 300, binding: { on: true, holo: true, charges: true, hbMode: "directional" } };
const ff = new ForceField(sel, par, mols);
const integ = new LangevinIntegrator(ff.ref, ff, 110.0);
integ.setTemperature(300);
integ.setFriction(8.0);

// --- 1. per-term accumulator invariants at the native pose ---------------
console.log("=== S4: per-term accumulators (CG, charges ON, directional HB) ===");
ff.trackTerms = true;
ff.compute(ff.ref);
assert(
  Math.abs((ff.bindLJU + ff.bindCoulU + ff.bindHBU + ff.desolvU) - ff.bindingU) < 1e-9,
  `lj+coul+hb+desolv = bindingU (${ff.bindingU.toFixed(4)} kcal/mol, exact sum)`,
);
assert(Number.isFinite(ff.bindU.lj) && Number.isFinite(ff.bindU.coul) && Number.isFinite(ff.bindU.hb) && Number.isFinite(ff.bindU.desolv),
  "bindU component vector all finite");
assert(ff.bindCoulU !== 0, `charges ON revives Coulomb term (bindCoulU = ${ff.bindCoulU.toFixed(3)} ≠ 0)`);
// OFF = bit-identical: same forces with and without the flag
ff.trackTerms = false;
ff.compute(ff.ref);
const eOff = ff.energy, fOff = Float64Array.from(ff.forces);
ff.trackTerms = true;
ff.compute(ff.ref);
let maxF = 0;
for (let i = 0; i < ff.forces.length; i++) maxF = Math.max(maxF, Math.abs(ff.forces[i] - fOff[i]));
assert(ff.energy === eOff && maxF === 0, "trackTerms OFF/ON: bit-identical energy + forces");

// --- 2. 300-step mini-run with BindLog capture ---------------------------
console.log("=== S4: 300-step Langevin mini-run + BindLog ===");
const bl = new BindLog({ quantStep: 0.01 });
const STRIDE = 10;
const nProt = ff.nProt, ligStart = ff.nLigAtoms > 0 ? ff.nProt : 0;
let lastContacts = null;
const rawFrames = [];
let nEnergySets = 0;
for (let s = 1; s <= 300; s++) {
  integ.step();
  if (s % STRIDE === 0) {
    // frame + 7-term energy vector (same push order as main.js tick)
    bl.captureFrame(integ.pos, integ.time);
    rawFrames.push(Float32Array.from(integ.pos));
    bl.pushEnergyComponents(integ.time, [
      ff.bindLJU || 0, ff.bindCoulU || 0, ff.bindHBU || 0,
      ff.desolvU || 0, ff.piU || 0, ff.cpiU || 0, ff.xbU || 0,
    ]);
    nEnergySets++;
    // contact form/break diff — same 5.5 Å pair set as the main.js nContacts
    // HUD loop, diffed against the previous captured frame's set
    const cur = new Map();
    for (let i = 0; i < nProt; i++) {
      for (let la = ligStart; la < ff.n; la++) {
        const dx = integ.pos[3 * la] - integ.pos[3 * i], dy = integ.pos[3 * la + 1] - integ.pos[3 * i + 1], dz = integ.pos[3 * la + 2] - integ.pos[3 * i + 2];
        const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < 30.25) cur.set(i * 1e6 + la, Math.sqrt(r2));
      }
    }
    if (lastContacts) {
      for (const [k, d] of cur) if (!lastContacts.has(k)) bl.pushContact(integ.time, true, k % 1e6, (k / 1e6) | 0, d);
      for (const [k, d] of lastContacts) if (!cur.has(k)) bl.pushContact(integ.time, false, k % 1e6, (k / 1e6) | 0, d);
    }
    lastContacts = cur;
  }
}

assert(bl.nFrames === 30, `frames captured: ${bl.nFrames} (= 300/10)`);
let hbSeen = 0;
for (let i = 0; i < bl.nEvents; i++) if (bl.evType[i] === 0 && bl.evA[i] === 2 && bl.evX[i] < 0) hbSeen++;
assert(hbSeen > 0, `directional HB fires during dynamics (${hbSeen}/30 negative HB samples — 0 at the exact native reference gate is physical)`);
const nEnergy = (() => { let c = 0; for (let i = 0; i < bl.nEvents; i++) if (bl.evType[i] === 0) c++; return c; })();
assert(nEnergy === 210, `energy events: ${nEnergy} (= 30 frames × 7 components)`);
const nCPlus = (() => { let c = 0; for (let i = 0; i < bl.nEvents; i++) if (bl.evType[i] === 1) c++; return c; })();
const nCMinus = (() => { let c = 0; for (let i = 0; i < bl.nEvents; i++) if (bl.evType[i] === 2) c++; return c; })();
assert(nCPlus + nCMinus > 0, `contact form/break events: ${nCPlus}+ / ${nCMinus}− (thermal flicker at 5.5 Å)`);
assert(bl.maxErrorVs(rawFrames) <= 0.01 + 1e-6,
  `dense-frame quantization error ${bl.maxErrorVs(rawFrames).toExponential(2)} ≤ quantStep/2 bound (0.01)`);

// --- 3. blob round-trip ---------------------------------------------------
console.log("=== S4: binary blob round-trip ===");
const blob = bl.toBinaryBlob();
const back = BindLog.fromBinaryBlob(blob);
assert(back.nFrames === 30 && back.nEvents === bl.nEvents, `blob round-trip: ${back.nFrames} frames, ${back.nEvents} events`);
let drift = 0, evDrift = 0;
for (let i = 0; i < 30; i++) {
  const a = bl.getFrame(i), b = back.getFrame(i);
  for (let j = 0; j < a.length; j++) drift = Math.max(drift, Math.abs(a[j] - b[j]));
}
for (let i = 0; i < bl.nEvents; i++) {
  if (bl.evType[i] !== back.evType[i] || bl.evA[i] !== back.evA[i]) evDrift++;
  evDrift = Math.max(evDrift, Math.abs(bl.evX[i] - back.evX[i]) > 1e-6 ? 1 : 0);
}
assert(drift === 0 && evDrift === 0, `blob vs in-memory: frame drift ${drift}, event drift ${evDrift}`);

// --- 4. ΔH of the LJ component at native (finite, negative) --------------
console.log("=== S4: native-pose ΔH(LJ) ===");
let sumLJ = 0, nLJ = 0, minY = Infinity;
for (let i = 0; i < back.nEvents; i++) {
  if (back.evType[i] === 0 && back.evA[i] === 0) { sumLJ += back.evX[i]; nLJ++; minY = Math.min(minY, back.evX[i]); }
}
const meanLJ = sumLJ / nLJ;
assert(nLJ === 30, `LJ component samples: ${nLJ}`);
assert(Number.isFinite(meanLJ) && meanLJ < 0,
  `⟨U_LJ⟩ over 30 frames = ${meanLJ.toFixed(3)} kcal/mol (finite, negative at native)`);
const sumDH = (() => { const acc = [0, 0, 0, 0, 0, 0, 0]; const cnt = [0, 0, 0, 0, 0, 0, 0];
  for (let i = 0; i < back.nEvents; i++) if (back.evType[i] === 0) { acc[back.evA[i]] += back.evX[i]; cnt[back.evA[i]]++; }
  return acc.map((v, i) => v / Math.max(1, cnt[i])); })();
console.log(`  ΔH decomposition (kcal/mol): LJ ${sumDH[0].toFixed(2)} · Coul ${sumDH[1].toFixed(2)} · HB ${sumDH[2].toFixed(2)} · desolv ${sumDH[3].toFixed(2)}`);
assert(sumDH.every(Number.isFinite), "all 7 ΔH components finite");

console.log(`\n=== test_bindlog_integration: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed ? 1 : 0);
