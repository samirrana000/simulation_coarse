/**
 * test_hbond.js — B17 H-bond directional
 *
 * In src/physics/hbond.js:71 evaluateDirectional(i,j,dx,dy,dz,r, donorAngle)
 * returns 0 at 90° vs -2.5 at 180° (linear). Keep existing evaluatePair as wrapper.
 * Runnable: node tests/test_hbond.js → PASS/FAIL
 */

import { DirectionalHBond, HBOND_EQ_DIST } from "../src/physics/hbond.js";

function assert(cond, msg) {
  if (!cond) { console.error(`FAIL: ${msg}`); process.exit(1); }
  console.log(`  ✓ ${msg}`);
}

function main() {
  console.log("=== B17 H-bond directional — 180°=-2.5, 90°~0 ===");
  const hb = new DirectionalHBond({ epsHB: 2.5 });
  const r = HBOND_EQ_DIST; // 2.9 Å equilibrium
  const dx = r, dy = 0, dz = 0;

  // Linear 180° should give -2.5 (full strength)
  const e180 = hb.evaluateDirectional(0, 1, dx, dy, dz, r, 180);
  console.log(`180° r=2.9 energy=${e180.energy.toFixed(6)} expected -2.5`);
  assert(Math.abs(e180.energy + 2.5) < 1e-6, `linear 180° energy -2.5 (got ${e180.energy})`);

  // Also test radians: π rad should also give -2.5
  const ePi = hb.evaluateDirectional(0, 1, dx, dy, dz, r, Math.PI);
  console.log(`π rad energy=${ePi.energy.toFixed(6)} expected -2.5`);
  assert(Math.abs(ePi.energy + 2.5) < 1e-6, `π rad energy -2.5 (got ${ePi.energy})`);

  // 90° should be ~0 (cos² gating)
  const e90 = hb.evaluateDirectional(0, 1, dx, dy, dz, r, 90);
  console.log(`90° energy=${e90.energy.toExponential(6)} expected ~0`);
  assert(Math.abs(e90.energy) < 1e-6, `90° energy ~0 (got ${e90.energy})`);

  const e90rad = hb.evaluateDirectional(0, 1, dx, dy, dz, r, Math.PI/2);
  console.log(`π/2 rad energy=${e90rad.energy.toExponential(6)} expected ~0`);
  assert(Math.abs(e90rad.energy) < 1e-6, `π/2 rad energy ~0 (got ${e90rad.energy})`);

  // evaluatePair wrapper should still be radial and equal to 180° (ideal)
  const ePair = hb.evaluatePair(0, 1, dx, dy, dz, r);
  console.log(`evaluatePair r=2.9 energy=${ePair.energy.toFixed(6)} (wrapper → 180°)`);
  assert(Math.abs(ePair.energy + 2.5) < 1e-6, `evaluatePair wrapper 180° == -2.5 (got ${ePair.energy})`);

  // Off-equilibrium r still scales but angular factor preserved
  const r2 = 3.5;
  const g = Math.exp(-((r2 - HBOND_EQ_DIST)**2)/(2*0.5*0.5));
  const expected180_r2 = -2.5 * g;
  const e180_r2 = hb.evaluateDirectional(0,1, r2,0,0, r2, 180);
  console.log(`180° r=3.5 energy=${e180_r2.energy.toFixed(6)} expected ${expected180_r2.toFixed(6)} g=${g.toFixed(4)}`);
  assert(Math.abs(e180_r2.energy - expected180_r2) < 1e-6, `180° r=3.5 scaled by g (got ${e180_r2.energy})`);
  const e90_r2 = hb.evaluateDirectional(0,1, r2,0,0, r2, 90);
  assert(Math.abs(e90_r2.energy) < 1e-6, `90° r=3.5 still ~0 (got ${e90_r2.energy})`);

  // 0° and 180° both max due to cos² symmetry (optional check)
  const e0 = hb.evaluateDirectional(0,1, dx,dy,dz,r, 0);
  assert(Math.abs(e0.energy + 2.5) < 1e-6, `0° energy -2.5 (cos² symmetric, got ${e0.energy})`);

  console.log("PASS: B17 hbond directional");
}

main();
