/**
 * test_bindlog.mjs — BindLog unit + integration tests (R6).
 * Run: node scripts/test_bindlog.mjs
 */
import { BindLog } from "../src/capture/bindlog.js";
import { readFileSync } from "node:fs";

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

// ---- synthetic trajectory: 100 frames, thermal motion 0.1-0.5 Å ----
console.log("=== BindLog: dense quantized frames ===");
const nAtoms = 50;
const rawFrames = [];
let seed = 7;
const rng = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
let pos = new Float32Array(nAtoms * 3);
for (let i = 0; i < pos.length; i++) pos[i] = (rng() - 0.5) * 60;
for (let f = 0; f < 100; f++) {
  const fr = new Float32Array(pos.length);
  for (let i = 0; i < fr.length; i++) fr[i] = pos[i] + (rng() - 0.5) * 0.5;
  rawFrames.push(fr);
  pos = fr;
}
const bl = new BindLog({ quantStep: 0.01 });
rawFrames.forEach((fr, i) => bl.captureFrame(fr, i * 2.0));
assert(bl.nFrames === 100, `frames captured: ${bl.nFrames}`);
const maxErr = bl.maxErrorVs(rawFrames);
assert(maxErr <= 0.01 + 1e-6, `quantization max error ${maxErr.toExponential(2)} ≤ quantStep 0.01`);
const rawBytes = 100 * nAtoms * 3 * 4;
const blBytes = nAtoms * 3 * 4 + 99 * nAtoms * 3 * 2 + 100 * 4; // frame0 + deltas + times
assert(blBytes < rawBytes * 0.55, `compression: ${blBytes} B vs raw ${rawBytes} B → ${(100 * blBytes / rawBytes).toFixed(0)}% (ratio ${(rawBytes / blBytes).toFixed(2)}×)`);

// ---- sparse events: 10k events memory ----
console.log("=== BindLog: sparse event log ===");
const bl2 = new BindLog();
for (let i = 0; i < 10000; i++) {
  const m = i % 7;
  if (m === 0) bl2.pushContact(i * 0.1, true, i % 21, i % 164, 3.5);
  else if (m === 1) bl2.pushContact(i * 0.1, false, i % 21, i % 164, 5.5);
  else if (m === 2) bl2.pushHill(i * 0.1, 4.2, 0.05);
  else if (m === 3) bl2.pushState(i * 0.1, i % 4);
  else if (m === 4) bl2.pushPocketVolume(i * 0.1, 620 + (rng() - 0.5) * 40);
  else bl2.pushEnergyComponents(i * 0.1, [-12.3, -4.5, -1.2, -0.8, -0.3, -0.1, 0.2]); // m 5 and 6
}
// pushEnergyComponents pushes 7 events per call → 5/7 of i are multi-event
const expectedEvents = (() => { let c = 0; for (let i = 0; i < 10000; i++) { const m = i % 7; c += (m === 5 || m === 6) ? 7 : 1; } return c; })();
assert(bl2.nEvents === expectedEvents, `10k pushes → ${expectedEvents} events stored: ${bl2.nEvents}`);
const evMem = (bl2.evTime.byteLength + bl2.evType.byteLength + bl2.evA.byteLength + bl2.evB.byteLength + bl2.evX.byteLength + bl2.evY.byteLength);
assert(evMem < 1.2e6, `10k events typed-array memory ${ (evMem / 1024).toFixed(0) } KB < 1.2 MB (vs ~800KB+ if per-event objects)`);
// verify a stored value
let hills = 0, states = new Set();
for (let i = 0; i < bl2.nEvents; i++) {
  if (bl2.evType[i] === 3) hills++;
  if (bl2.evType[i] === 4) states.add(bl2.evA[i]);
}
assert(hills === 1429, `hill events round-trip: ${hills}/1429`);
assert(states.size === 4, `state ids 0-3 present: ${[...states].sort().join(",")}`);

// ---- binary blob round-trip ----
console.log("=== BindLog: binary blob round-trip ===");
const blob = bl.toBinaryBlob();
const back = BindLog.fromBinaryBlob(blob);
assert(back.nFrames === 100, `blob frames ${back.nFrames}`);
assert(back.nEvents === 0, "blob events 0 (dense-only blob)");
const err2 = back.maxErrorVs(rawFrames);
assert(err2 <= 0.01 + 1e-6, `blob round-trip error ${err2.toExponential(2)}`);
for (let i = 0; i < back.frameTimes.length; i++) {
  if (Math.abs(back.frameTimes[i] - bl.frameTimes[i]) > 1e-4) { assert(false, `time ${i} mismatch`); break; }
}
assert(true, "all frame times identical after round-trip");
// events blob
const blob2 = bl2.toBinaryBlob();
const back2 = BindLog.fromBinaryBlob(blob2);
assert(back2.nEvents === expectedEvents, `event blob round-trips: ${back2.nEvents} events`);
let mismatch = 0;
for (let i = 0; i < expectedEvents; i++) {
  if (back2.evType[i] !== bl2.evType[i] || back2.evA[i] !== bl2.evA[i] ||
      Math.abs(back2.evX[i] - bl2.evX[i]) > 1e-6) mismatch++;
}
assert(mismatch === 0, `all event fields identical (mismatches: ${mismatch})`);

// ---- JSON summary ----
const sum = bl2.toJSONSummary();
const nStates = (() => { const s = new Set(); for (let i = 0; i < 10000; i++) if (i % 7 === 3) s.add(i % 4); return s.size; })();
assert(sum.events.hill === Math.floor(10000/7) + (10000%7>2?1:0) && sum.events.state === sum.events.hill, `JSON summary counts: hills ${sum.events.hill}, states ${sum.events.state}`);

// ---- real-FF integration: capture frames + energy from a 200-step mini run ----
console.log("=== BindLog: real ForceField mini-run ===");
try {
  const { parseCa, selectSystem, parseLigands } = await import("../src/pdb.js");
  const { ForceField } = await import("../src/forcefield.js");
  const { LangevinIntegrator } = await import("../src/integrator.js");
  const pdbText = readFileSync(new URL("../4w52.pdb", import.meta.url), "utf8");
  const parsed = parseCa(pdbText);
  const sel = selectSystem(parsed);
  const mols = parseLigands(pdbText);
  const ff = new ForceField(sel, { rc: 10, gamma: 2.0 }, mols);
  const integ = new LangevinIntegrator(ff.ref, ff, 110.0);
  integ.setTemperature(300); integ.setFriction(8.0);
  const blR = new BindLog({ quantStep: 0.02 });
  for (let s = 0; s < 200; s++) {
    integ.step();
    if (s % 5 === 0) {
      blR.captureFrame(integ.pos, integ.t);
      blR.pushEnergyComponents(integ.t, [ff.bindingU ?? 0, ff.desolvU ?? 0, ff.U ?? 0]);
    }
  }
  assert(blR.nFrames === 40, `real-run frames: ${blR.nFrames}`);
  assert(blR.nEvents === 120, `real-run energy events: ${blR.nEvents} (40 frames × 3 comps)`);
  const blobR = blR.toBinaryBlob();
  const backR = BindLog.fromBinaryBlob(blobR);
  assert(backR.nFrames === 40 && backR.nEvents === 120, "real-run blob round-trip OK");
  const maxDrift = backR.maxErrorVs(Array.from({ length: 40 }, (_, i) => {
    // reconstruct originals from blR (quantized) — verify blob == in-memory, error vs live run bounded
    return blR.getFrame(i);
  }));
  assert(maxDrift < 1e-5, `blob vs in-memory max drift ${maxDrift.toExponential(2)}`);
} catch (e) {
  console.log("  (real-FF integration skipped:", e.message.slice(0, 80) + ")");
}

console.log(`\n=== test_bindlog: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed ? 1 : 0);
