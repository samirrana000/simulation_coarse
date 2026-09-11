/**
 * worker-pool.js — Worker Pool manager for CPU multi-processing.
 *
 * Spawns and coordinates Web Workers across available CPU cores.
 */

/**
 * Guard: ensure an ArrayBuffer-backed typed array was not detached
 * by an accidental transfer of its underlying buffer. Transferring
 * `pos.buffer` directly would detach the original `pos`; we clone
 * via `pos.slice()` and transfer the clone so the source stays live.
 * @param {Float64Array} pos
 * @param {number} beforeLen - pos.byteLength before transfer
 */
function assertNotDetached(pos, beforeLen) {
  if (pos.byteLength !== beforeLen) {
    throw new Error(`WorkerPool detach guard: pos.byteLength changed ${beforeLen} -> ${pos.byteLength} (buffer was detached)`);
  }
  if (pos.buffer.byteLength !== beforeLen) {
    throw new Error(`WorkerPool detach guard: underlying buffer detached`);
  }
}

export class WorkerPool {
  /**
   * @param {number} [numWorkers] Number of worker threads (defaults to navigator.hardwareConcurrency - 1 or 2)
   */
  constructor(numWorkers) {
    const defaultWorkers = typeof navigator !== "undefined" && navigator.hardwareConcurrency
      ? Math.max(1, Math.min(8, navigator.hardwareConcurrency - 1))
      : 2;
    this.numWorkers = numWorkers ?? defaultWorkers;
    this.workers = [];
    this.ready = false;
    this._busy = false; // guard against concurrent computeParallel
    this.initWorkers();
  }

  initWorkers() {
    if (typeof Worker === "undefined") {
      this.ready = false;
      return;
    }

    try {
      this.terminate();
      this.workers = [];
      for (let i = 0; i < this.numWorkers; i++) {
        const w = new Worker(new URL("./force-worker.js", import.meta.url), { type: "module" });
        this.workers.push(w);
      }
      this.ready = true;
    } catch (e) {
      console.warn("WorkerPool initialization fallback to single-thread:", e);
      this.ready = false;
    }
  }

  setNumWorkers(n) {
    const clamped = Math.max(1, Math.min(16, n));
    if (clamped !== this.numWorkers) {
      this.numWorkers = clamped;
      this.initWorkers();
    }
  }

  initSystem(ff) {
    if (!this.ready || this.workers.length === 0) return;

    const excludedList = Array.from(ff._excluded || []);
    const scale14List = Array.from((ff._scale14 || new Map()).entries());

    for (const w of this.workers) {
      w.postMessage({
        type: "init",
        elem: ff._elem,
        excluded: excludedList,
        scale14: scale14List,
        bornRadii: ff._bornRadii,
      });
    }
  }

  /**
   * Dispatches parallel force evaluation batches across all workers.
   *
   * @param {Float64Array} pos
   * @param {number} n
   * @returns {Promise<{forces: Float64Array, lj: number, elec: number}>}
   */
  async computeParallel(pos, n) {
    if (this._busy) throw new Error("WorkerPool: computeParallel already in progress (concurrent calls not allowed)");
    if (!this.ready || this.workers.length <= 1) {
      return null;
    }
    this._busy = true;
    try {
      const nWorkers = this.workers.length;
      const chunkSize = Math.ceil(n / nWorkers);
      const promises = [];

      for (let w = 0; w < nWorkers; w++) {
        const atomStart = w * chunkSize;
        const atomEnd = Math.min(n, (w + 1) * chunkSize);
        if (atomStart >= n) break;

        const worker = this.workers[w];
        // Guard: clone pos via pos.slice() so the transfer does not detach the original pos buffer.
        // We transfer posCopy.buffer, not pos.buffer, and assert pos.byteLength unchanged before/after.
        const beforeLen = pos.byteLength;
        const posCopy = pos.slice(); // Float64Array clone — owns its own ArrayBuffer
        // Verify original pos not detached by the slice itself
        assertNotDetached(pos, beforeLen);

        const p = new Promise((resolve, reject) => {
          const handler = (e) => {
            if (e.data.type === "batch_result") {
              cleanup();
              resolve(e.data);
            }
          };
          const errHandler = (err) => {
            cleanup();
            reject(err);
          };
          const cleanup = () => {
            worker.removeEventListener("message", handler);
            worker.removeEventListener("error", errHandler);
          };
          worker.addEventListener("message", handler);
          worker.addEventListener("error", errHandler);
          worker.postMessage(
            {
              type: "compute_batch",
              batchId: w,
              atomStart,
              atomEnd,
              n,
              posBuffer: posCopy.buffer,
            },
            [posCopy.buffer]
          );
          // After postMessage with transfer of posCopy.buffer, the original pos must remain attached
          assertNotDetached(pos, beforeLen);
        });

        promises.push(p);
      }

      const results = await Promise.all(promises);
      const totalForces = new Float64Array(n * 3);
      let totalLJ = 0;
      let totalElec = 0;

      for (const res of results) {
        const f = new Float64Array(res.forcesBuffer);
        for (let i = 0; i < totalForces.length; i++) totalForces[i] += f[i];
        totalLJ += res.lj;
        totalElec += res.elec;
      }

      return { forces: totalForces, lj: totalLJ, elec: totalElec };
    } finally {
      this._busy = false;
    }
  }

  terminate() {
    for (const w of this.workers) {
      try { w.terminate(); } catch (_) {}
    }
    this.workers = [];
  }
}
