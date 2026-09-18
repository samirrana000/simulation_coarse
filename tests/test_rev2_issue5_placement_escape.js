/**
 * test_rev2_issue5_placement_escape.js — Revolution 2 / Issue 5 regression.
 *
 * Scope: src/placement.js coincident-escape branch ONLY (clashGrad r2 < 1e-6,
 * src/placement.js:89). No energy-kernel or caller change.
 *
 * Problem: the coincident degenerate case pushed a fixed +x+y+z diagonal
 * (Fx += 1, Fy += 1, Fz += 1, minRatio = 0) — brittle when +x+y+z is walled
 * in a dense hetero shell (net push walks into the wall, step halves to
 * <1e-10, pose stalls), and minRatio = 0 reported residualClash = Infinity
 * (non-finite degenerate pose) even though the clash energy itself is finite.
 *
 * Fix verified here (minimal, backward compatible):
 *   - Deterministic index-hashed escape direction per (ligand a, collider i)
 *     plus the matching torque arm (rotation participates in the escape).
 *     Fixed given indices, so the full placement stays deterministic given
 *     seed (seed still owns the initial SO(3) rotation).
 *   - minRatio floor (1e-3/rE) instead of 0, so the residual (1/minRatio)
 *     stays finite even for a still-coincident pose.
 *
 * Coverage (headless, no fixtures):
 *   [0] single-atom coincident escape leaves along a non-diagonal direction
 *       (cos with +x+y+z < 0.9; the old code gives exactly 1.0), converges,
 *       finite residual, bit-identical twice.
 *   [1] still-coincident pose (maxIters:0) reports a FINITE residual
 *       (≈1000·rE, old code: Infinity) and honestly converged:false.
 *   [2] dense hetero cage (octahedral ±2.5 Å shell + coincident center):
 *       finite residual, escape attempted (iters > 0), bit-identical twice.
 *   [3] placeLigand seed determinism in the dense cage (same seed twice
 *       bit-identical, finite residual).
 *
 * Run: node tests/test_rev2_issue5_placement_escape.js (fast, <2s; NOT wired
 * into tests/test_all.js FAST so the 352 gate is untouched).
 */

import { placeLigand, relaxClash } from "../src/placement.js";

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

function maxAbsDiff(a, b) {
  let m = 0;
  const n = Math.min(a.length, b.length);
  for (let i = 0; i < n; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
  return m;
}

function main() {
  console.log("=== Rev2/Issue5: deterministic coincident escape (no fixed diagonal) ===");

  // -----------------------------------------------------------------
  // [0] Single-atom coincident escape is NOT the fixed +x+y+z diagonal
  // -----------------------------------------------------------------
  console.log("\n[0] single-atom coincident escape direction...");
  const mol1 = { atoms: [{ element: "C", x: 0, y: 0, z: 0 }] };
  const single = (extra = {}) =>
    relaxClash(new Float64Array([0, 0, 0]), mol1,
      { pos: new Float64Array([0, 0, 0]), sigma: new Float64Array([3.4]) }, extra);
  const s1 = single();
  const disp = [s1.pos[0], s1.pos[1], s1.pos[2]];
  const len = Math.hypot(...disp);
  const cosDiag = (disp[0] + disp[1] + disp[2]) / (len * Math.sqrt(3));
  assert(s1.converged === true, `single coincident escapes to converged (iters=${s1.iterations})`);
  assert(Number.isFinite(s1.residualClash), `single escape residual finite (${s1.residualClash.toFixed(3)})`);
  assert(len > 1.0, `single escape moves clear of coincidence (|d|=${len.toFixed(3)} Å)`);
  assert(cosDiag < 0.9, `escape NOT along fixed +x+y+z diagonal (cos=${cosDiag.toFixed(3)} < 0.9; old code = 1.0)`);
  const s2 = single();
  assert(maxAbsDiff(s1.pos, s2.pos) === 0 && s1.residualClash === s2.residualClash,
    `single escape bit-identical twice (deterministic given geometry)`);

  // -----------------------------------------------------------------
  // [1] Still-coincident pose reports a FINITE residual (no Infinity)
  // -----------------------------------------------------------------
  console.log("\n[1] still-coincident residual finite...");
  const mol3 = {
    atoms: [
      { element: "C", x: 0, y: 0, z: 0 },
      { element: "C", x: 1.4, y: 0, z: 0 },
      { element: "O", x: -1.2, y: 0.5, z: 0 },
    ],
  };
  const coincidentPos = () => new Float64Array([0, 0, 0, 1.4, 0, 0, -1.2, 0.5, 0]);
  const stuck = relaxClash(coincidentPos(), mol3,
    { pos: new Float64Array([0, 0, 0]), sigma: new Float64Array([3.4]) }, { maxIters: 0 });
  assert(stuck.converged === false, `unrelaxed coincident honestly unconverged`);
  assert(Number.isFinite(stuck.residualClash),
    `unrelaxed coincident residual finite (${stuck.residualClash.toFixed(1)}; old code: Infinity)`);
  assert(stuck.residualClash > 100,
    `unrelaxed coincident residual still signals clash (${stuck.residualClash.toFixed(1)} > 100)`);

  // -----------------------------------------------------------------
  // [2] Dense hetero cage + coincident center: finite + deterministic
  // -----------------------------------------------------------------
  console.log("\n[2] dense hetero cage with coincident center...");
  const cageXYZ = [2.5, 0, 0, -2.5, 0, 0, 0, 2.5, 0, 0, -2.5, 0, 0, 0, 2.5, 0, 0, -2.5, 0, 0, 0];
  const cagePos = new Float64Array(cageXYZ);
  const cageSig = new Float64Array(cageXYZ.length / 3).fill(3.4);
  const c1 = relaxClash(coincidentPos(), mol3, { pos: cagePos, sigma: cageSig });
  assert(Number.isFinite(c1.residualClash), `cage+coincident residual finite (${c1.residualClash.toFixed(3)})`);
  assert(c1.iterations > 0, `cage+coincident escape attempted (${c1.iterations} iters)`);
  assert(c1.converged === true, `cage+coincident converges off the degenerate point`);
  const c2 = relaxClash(coincidentPos(), mol3, { pos: cagePos, sigma: cageSig });
  assert(maxAbsDiff(c1.pos, c2.pos) === 0 && c1.residualClash === c2.residualClash,
    `cage+coincident bit-identical twice (deterministic)`);

  // -----------------------------------------------------------------
  // [3] placeLigand seed determinism in the dense cage
  // -----------------------------------------------------------------
  console.log("\n[3] placeLigand seed determinism (dense cage)...");
  const target = [0, 0, 0];
  const p1 = placeLigand(mol3, target, { protein: { pos: cagePos, sigma: cageSig }, seed: 7 });
  const p2 = placeLigand(mol3, target, { protein: { pos: cagePos, sigma: cageSig }, seed: 7 });
  assert(Number.isFinite(p1.residualClash), `dense-cage placement residual finite (${p1.residualClash.toFixed(3)})`);
  assert(maxAbsDiff(p1.pos, p2.pos) === 0 && p1.residualClash === p2.residualClash,
    `same seed twice bit-identical (deterministic given seed)`);

  // -----------------------------------------------------------------
  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log("=================================================");
  if (failed > 0) process.exit(1);
  console.log("PASS: test_rev2_issue5_placement_escape.js — Rev2/Issue5 escape validated");
}

main();
