/**
 * force-worker.js — Web Worker script for parallel force evaluation.
 *
 * Runs non-bonded and bonded force evaluations on a separate CPU core/thread.
 */

/* global self */

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
          ee = 332.0 * qq * fac;
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
