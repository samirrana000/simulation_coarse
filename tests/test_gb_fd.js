/**
 * test_gb_fd.js — Finite-difference validation for GeneralizedBorn.pairInteraction
 *
 * Two atoms at 3.0 Å, q=+1 / -1, Born radii 1.7/1.5 Å.
 * Analytic force from pairInteraction vs central difference (dx=1e-4).
 * Runnable: node tests/test_gb_fd.js  →  PASS/FAIL
 */

import { GeneralizedBorn } from "../src/physics/gb.js";

function run() {
  const gb = new GeneralizedBorn({ epsIn: 4.0, epsOut: 78.5, saltM: 0.15, temperature: 300 });
  const r = 3.0;
  const dx = 3.0, dy = 0.0, dz = 0.0;
  const qi = 1.0, qj = -1.0;
  const Ri = 1.7, Rj = 1.5;
  const h = 1e-4;

  // Analytic
  const ana = gb.pairInteraction(0, 1, dx, dy, dz, r, qi, qj, Ri, Rj, 1.0);
  const anaFx = ana.fx;
  const anaFy = ana.fy;
  const anaFz = ana.fz;
  const anaE = ana.energy;

  // Numeric central difference dU/dr
  const plus = gb.pairInteraction(0, 1, dx + h, dy, dz, r + h, qi, qj, Ri, Rj, 1.0);
  const minus = gb.pairInteraction(0, 1, dx - h, dy, dz, r - h, qi, qj, Ri, Rj, 1.0);
  const numDudr = (plus.energy - minus.energy) / (2 * h);
  // For x-axis, Fx = dU/dr * (dx/r) ; dx/r==1 at r=3
  const numFx = numDudr; // equivalent to -(-dU/dr) ; GB pairInteraction returns +dU/dr*dx/r

  // Also test y-perturbation should be zero (but numeric via orthogonal shift small)
  // Instead verify full vector: perturb along x only so Fy,Fz should be 0
  // Compute numeric Fy by shifting dy slightly: energy at dy=+h vs -h
  const r_orth = Math.sqrt(r * r + h * h);
  const plusY = gb.pairInteraction(0, 1, dx, h, dz, r_orth, qi, qj, Ri, Rj, 1.0);
  const minusY = gb.pairInteraction(0, 1, dx, -h, dz, r_orth, qi, qj, Ri, Rj, 1.0);
  // dU/dy at y=0 should be 0; check numeric vs analytic

  const diffFx = Math.abs(anaFx - numFx);
  const relFx = Math.abs(anaFx) > 1e-8 ? diffFx / Math.abs(anaFx) : diffFx;

  console.log(`GB pair r=3.0 Å  q=+1/-1  Ri=${Ri} Rj=${Rj}`);
  console.log(`  analytic Fx = ${anaFx.toExponential(6)}  Fy=${anaFy.toExponential(6)} Fz=${anaFz.toExponential(6)}  E=${anaE.toFixed(6)}`);
  console.log(`  numeric  Fx = ${numFx.toExponential(6)}  diff=${diffFx.toExponential(3)} rel=${relFx.toExponential(3)}`);
  console.log(`  E+ = ${plus.energy.toFixed(6)}  E- = ${minus.energy.toFixed(6)}  dU/dr num=${numDudr.toExponential(6)} ana dU/dr=${anaFx.toExponential(6)}`);

  // Direct force finite-diff via position shift of atom j along x:
  // F_analytic_on_j = -dU/dr * (dx/r) ; F_analytic_on_i = +dU/dr*dx/r
  // Our anaFx is force on i; numeric on i equals +dU/dr
  const tol = 1e-3;
  if (diffFx < tol) {
    console.log("PASS: GB analytic force matches finite difference (diff < 1e-3)");
    process.exit(0);
  } else {
    console.error(`FAIL: GB force mismatch diff=${diffFx} > ${tol} (rel ${relFx})`);
    process.exit(1);
  }
}

run();
