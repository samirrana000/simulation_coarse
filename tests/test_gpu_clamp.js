/**
 * tests/test_gpu_clamp.js — G64 GPU LJ clamp sr<5
 *
 * Checks that sigma/r is clamped to 5: sr=200 → 5
 * Uses Math.min(sigma/r,5) clamp as in src/gpu.js
 *
 * Runnable: node tests/test_gpu_clamp.js
 */

import { gpuClampSr } from "../src/gpu.js";

function clampSr(sigma, r) {
  return Math.min(sigma / r, 5); // clamp sigma/r to 5 — G64
}

console.log("=== G64 GPU clamp sr<5 — sr=200 → 5 ===");

// Direct Math.min test
const srRaw = 200;
const srClamped = Math.min(srRaw, 5);
console.log(`clamp: Math.min(${srRaw},5) = ${srClamped}`);
if (srClamped !== 5) {
  console.error(`FAIL: clamp ${srRaw} -> ${srClamped} expected 5`);
  process.exit(1);
}

// Via gpuClampSr helper: sigma=20, r=0.1 => sigma/r=200 -> 5
const sigma = 20;
const r = 0.1;
const sr = gpuClampSr(sigma, r);
console.log(`gpuClampSr(sigma=${sigma}, r=${r}) = ${sr} (sigma/r=${sigma/r})`);
if (sr !== 5) {
  console.error(`FAIL: gpuClampSr ${sigma}/${r}=200 should clamp to 5, got ${sr}`);
  process.exit(1);
}

// Also test non-clamped case: sigma=3.4, r=4 => 0.85 -> 0.85
const sr2 = gpuClampSr(3.4, 4);
console.log(`gpuClampSr(3.4,4) = ${sr2} (expected 0.85)`);
if (Math.abs(sr2 - 0.85) > 1e-9) {
  console.error(`FAIL: expected 0.85 got ${sr2}`);
  process.exit(1);
}

console.log("PASS: G64 clamp sr=200→5 verified (Math.min(sigma/r,5) + clamp)");
