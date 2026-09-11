/**
 * tests/test_recorder_cap.js — G69 memory leak guard maxFrames cap
 *
 * Verifies recorder caps at maxFrames and auto-stops.
 *
 * Runnable: node tests/test_recorder_cap.js
 */

import { Recorder } from "../src/recorder.js";

console.log("=== G69 Recorder maxFrames cap ===");

const rec = new Recorder();
const maxFrames = 5;
rec.start(0, 1.0, maxFrames); // stride 1ps, cap 5
console.log(`Started recorder maxFrames=${maxFrames} (G69 cap)`);

const pos = new Float64Array([0,0,0, 1,1,1]);
for (let i = 0; i < 10; i++) {
  // advance time by 1ps each frame
  rec.maybeCapture(pos, i * 1.0);
  // mutate pos slightly to ensure copy
  pos[0] += 0.01;
}

console.log(`Recorder count after 10 capture attempts: ${rec.count} (maxFrames=${maxFrames})`);
if (rec.count !== maxFrames) {
  console.error(`FAIL: recorder caps at maxFrames — expected ${maxFrames} got ${rec.count}`);
  process.exit(1);
}

if (rec.recording) {
  console.error(`FAIL: recorder should auto-stop after reaching maxFrames`);
  process.exit(1);
}

// Ensure no leak on extra captures beyond cap
rec.maybeCapture(pos, 20);
if (rec.count !== maxFrames) {
  console.error(`FAIL: extra capture beyond cap changed count`);
  process.exit(1);
}

console.log(`PASS: recorder caps at maxFrames=${maxFrames} — G69 memory leak guard verified`);
