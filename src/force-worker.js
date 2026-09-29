/**
 * force-worker.js — Web Worker script for parallel force evaluation.
 *
 * Runs non-bonded and bonded force evaluations on a separate CPU core/thread.
 *
 * The Coulomb constant is imported from the single src/units.js contract.
 * It previously carried a bare `332.0` inline in the Coulomb branch, which is
 * 1.9e-4 below the CODATA value every other non-bonded kernel uses.
 *
 * DELIBERATE PER-MODULE EXCEPTION — `R_CUT = 8.5` below is NOT a duplicate of
 * src/heavy.js's `R_CUT`, and must not be consolidated with it:
 *   • this kernel's electrostatics is exponential-screened, exp(-r/8)/r^2,
 *     whereas heavy.js uses a smooth C2 switch over [6.5, 8.5] plus GB. They
 *     are different potentials that happen to share a truncation radius, so
 *     binding them would assert a physical identity that does not hold;
 *   • the value is already parameterised upstream for the GPU path
 *     (compute/webgpu_backend.js RCUT_DEFAULT, its own default), i.e. it is a
 *     per-backend tuning constant, not a unit;
 *   • importing heavy.js to reach one number would add a 26-module / 377 KB
 *     static-import closure to EVERY worker thread (measured), for a constant
 *     that is not shared. units.js, by contrast, is a zero-import leaf.
 * See tests/test_constant_ledger.js ALLOW_MULTI_SITE for the guard entry.
 */

/* global self */

// Loaded as `new Worker(url, { type: "module" })` (see worker-pool.js), so
// static ES imports are available. units.js has zero imports (it is the leaf
// of the dependency graph), so this edge cannot create a cycle.
import { COULOMB_CONST } from "./units.js";

let system = null;
let elem = null;
let excluded = null;
let scale14 = null;
let bornRadii = null;

self.onmessage = function (e) {
  const msg = e.data;
  if (!msg) return;

  if (msg.type === "init") {
    system = msg.system;
    elem = msg.elem;
    excluded = new Set(msg.excluded);
    scale14 = new Map(msg.scale14);
    bornRadii = msg.bornRadii;
    self.postMessage({ type: "init_ok" });
  } else if (msg.type === "compute_batch") {
    const pos = new Float64Array(msg.posBuffer);
    const atomStart = msg.atomStart;
    const atomEnd = msg.atomEnd;
    const n = msg.n;
    const forces = new Float64Array(n * 3);

    let ljTot = 0, elecTot = 0;
    // Per-backend tuning constant — see the header note on why this is not
    // src/heavy.js's R_CUT despite the identical value.
    const R_CUT = 8.5;
    const cut2 = R_CUT * R_CUT;

    for (let i = atomStart; i < atomEnd; i++) {
      const xi = 3 * i;
      const ei = elem[i];
      const qi = ei.q;

      for (let j = i + 1; j < n; j++) {
        const k = i < j ? i * 1e6 + j : j * 1e6 + i;
        if (excluded.has(k)) continue;
        const s14 = scale14.get(k) ?? 1.0;

        const xj = 3 * j;
        const dx = pos[xj] - pos[xi];
        const dy = pos[xj + 1] - pos[xi + 1];
        const dz = pos[xj + 2] - pos[xi + 2];
        const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 >= cut2 || r2 <= 1e-4) continue;

        const r = Math.sqrt(r2);
        const ej = elem[j];
        const qj = ej.q;

        // LJ
        const s = 0.5 * (ei.sigma + ej.sigma);
        const eps = Math.sqrt(ei.eps * ej.eps);
        const sr = s / r, sr6 = sr * sr * sr * sr * sr * sr;
        const ljE = 4 * eps * (sr6 * sr6 - sr6);
        const ljF = 4 * eps * (12 * sr6 * sr6 - 6 * sr6) / r;

        // Coulomb
        let ee = 0, ef = 0;
        if (qi !== 0 && qj !== 0) {
          const qq = qi * qj;
          const fac = Math.exp(-r / 8.0) / (r * r);
          ee = COULOMB_CONST * qq * fac;
          ef = ee * (-1 / 8.0 - 2 / r);
        }

        const totE = s14 * (ljE + ee);
        const totF = s14 * (-ljF + ef);

        ljTot += s14 * ljE;
        elecTot += s14 * ee;

        const fx = totF * dx / r;
        const fy = totF * dy / r;
        const fz = totF * dz / r;

        forces[xi] += fx;
        forces[xi + 1] += fy;
        forces[xi + 2] += fz;
        forces[xj] -= fx;
        forces[xj + 1] -= fy;
        forces[xj + 2] -= fz;
      }
    }

    self.postMessage(
      {
        type: "batch_result",
        batchId: msg.batchId,
        forcesBuffer: forces.buffer,
        lj: ljTot,
        elec: elecTot,
      },
      [forces.buffer]
    );
  }
};
