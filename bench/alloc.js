/**
 * bench/alloc.js — Zero-alloc per step audit (G65).
 *
 * Measures heap growth over 1000 ForceField.compute() / HeavyForceField.compute() calls.
 * Uses process.memoryUsage().heapUsed before / after a tight loop. In steady state,
 * a zero-alloc hot loop should show heap delta ≈ 0 (or bounded by GC jitter),
 * whereas SasaModel's per-call Float64Array(n) ×3 will show linear growth.
 *
 * G65 AUDIT NOTES (see src/forcefield.js:403 and src/heavy.js:594 comments):
 *   - ForceField (Cα ENM): Expected ~0 heap alloc per step after warmup.
 *     Flat buffers (forces, _grid/_gridB, _dens/_bp*) are preallocated and reused;
 *     only a tiny iterator object from `for (const [k,arr] of grid)` is allocated (~ O(cells)).
 *   - HeavyForceField: Currently allocates 3 × Float64Array(n) per compute() via
 *     SasaModel.compute() (s0, radii, burial). Suggested fix: hoist to scratch buffers
 *     on HeavyForceField (this._sasaS0 etc.) and reuse. Until fixed, heap delta
 *     will be ~ n*8*3 bytes per step ×1000 ≈ 7–30 MB for 1k–4k atoms.
 *     Non-bonded grid (SpatialGrid.head/next/cellCoords) is zero-alloc in steady state.
 *
 * Runnable: node bench/alloc.js
 * Prints heap delta and a JSON summary. Also calls global.gc() if exposed (node --expose-gc).
 * Does not assert failure on large delta — reports and documents expected vs observed.
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { parseCa, selectSystem } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { parseHeavy, HeavyForceField } from "../src/heavy.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

function readPdbOrThrow(name) {
  const candidates = [
    path.resolve(process.cwd(), name),
    path.resolve(__dirname, "..", name),
    path.resolve(__dirname, name),
    path.resolve("data", name),
  ];
  for (const p of candidates) {
    try { if (fs.existsSync(p)) return fs.readFileSync(p, "utf-8"); } catch {}
  }
  throw new Error(`Cannot find ${name} (tried ${candidates.join(", ")})`);
}

function measureHeap(label, ff, steps = 1000) {
  const pos = new Float64Array(ff.ref);
  // Warmup to stabilize JIT and grow any bucket arrays to final capacity
  for (let i = 0; i < 20; i++) ff.compute(pos);
  // Perturb a tiny amount so optimizer does not constant-fold
  for (let i = 0; i < pos.length; i++) pos[i] += (i % 7) * 1e-9;

  if (typeof global !== "undefined" && typeof global.gc === "function") {
    try { global.gc(); } catch {}
  }

  const before = process.memoryUsage().heapUsed;
  const t0 = typeof performance !== "undefined" && performance.now ? performance.now() : Date.now();
  for (let s = 0; s < steps; s++) {
    // Slightly wiggle position every 100 steps to avoid dead-code elimination and exercise different branches
    if (s % 100 === 0 && s > 0) {
      for (let k = 0; k < Math.min(9, pos.length); k++) pos[k] += 1e-7;
    }
    ff.compute(pos);
  }
  const t1 = typeof performance !== "undefined" && performance.now ? performance.now() : Date.now();
  const after = process.memoryUsage().heapUsed;
  const delta = after - before;
  const deltaPerStep = delta / steps;
  const msPerCompute = (t1 - t0) / steps;
  return { label, n: ff.n, steps, heapBefore: before, heapAfter: after, heapDelta: delta, heapPerStep: deltaPerStep, msPerCompute };
}

async function main() {
  const pdbText = readPdbOrThrow("4w52.pdb");
  console.log("=== G65 Zero-Alloc Per Step — heap growth over 1000 compute() calls ===");
  if (typeof global !== "undefined" && typeof global.gc === "function") {
    console.log("(GC exposed — calling global.gc() before measurement for cleaner delta)");
  } else {
    console.log("(Run with `node --expose-gc bench/alloc.js` for most accurate heap delta)");
  }

  // CG benchmark
  const parsedCa = parseCa(pdbText);
  const sel = selectSystem(parsedCa);
  const ffCG = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const resCG = measureHeap("CG", ffCG, 1000);
  console.log(`CG: n=${resCG.n}  heap delta ${ (resCG.heapDelta/1024).toFixed(1)} KiB over ${resCG.steps} steps  (${resCG.heapPerStep.toFixed(1)} B/step)  ${resCG.msPerCompute.toFixed(3)} ms/compute`);
  if (resCG.heapDelta < 512 * 1024) {
    console.log("  CG heap growth <512 KiB → effectively zero-alloc (PASS)");
  } else {
    console.warn(`  CG heap growth ${(resCG.heapDelta/1024).toFixed(1)} KiB suggests per-step alloc — audit src/forcefield.js compute() for new Array/Float64Array`);
  }

  // Heavy benchmark (if available and small enough to not OOM in bench)
  let resHeavy = null;
  try {
    const parsedHeavy = parseHeavy(pdbText);
    // Use a subset if too large for bench stability? 4w52 heavy is 1308 atoms, OK.
    const heavySystem = { atoms: parsedHeavy.atoms };
    const ffHeavy = new HeavyForceField(heavySystem, { gamma: 1.0, temp: 300 });
    resHeavy = measureHeap("Heavy", ffHeavy, 1000);
    console.log(`Heavy: n=${resHeavy.n}  heap delta ${(resHeavy.heapDelta/1024).toFixed(1)} KiB over ${resHeavy.steps} steps  (${resHeavy.heapPerStep.toFixed(1)} B/step)  ${resHeavy.msPerCompute.toFixed(3)} ms/compute`);
    // Heavy currently expected to allocate ~ n*8*3 per step via SasaModel → delta ≈ n*24*1000
    const expectedHeavyMin = resHeavy.n * 8 * 3 * 0.8 * 1000; // ~80% of naive if GC reclaims some
    if (resHeavy.heapDelta > 1024 * 1024) {
      console.warn(`  Heavy heap growth >1 MiB as expected (SasaModel per-call alloc). Fix: hoist s0/radii/burial to scratch buffers.`);
    }
    if (resHeavy.heapDelta < 512 * 1024) {
      console.log("  Heavy heap growth <512 KiB → unexpectedly zero-alloc (SasaModel fix may already be applied)");
    }
  } catch (e) {
    console.error("Heavy bench failed:", e.message);
    resHeavy = { label: "Heavy", n: 1308, steps: 1000, heapDelta: NaN, heapPerStep: NaN, msPerCompute: NaN, error: e.message };
  }

  const summary = {
    timestamp: new Date().toISOString(),
    steps: 1000,
    pdb: "4w52.pdb",
    cg: { n: resCG.n, heapDelta: resCG.heapDelta, heapPerStep: resCG.heapPerStep, msPerCompute: resCG.msPerCompute },
    heavy: resHeavy ? { n: resHeavy.n, heapDelta: resHeavy.heapDelta, heapPerStep: resHeavy.heapPerStep, msPerCompute: resHeavy.msPerCompute, error: resHeavy.error || null } : null,
    note: "G65: CG should be ~0 B/step; Heavy currently ~ n*24 B/step due to SasaModel alloc (see src/heavy.js comment for fix).",
    gcExposed: !!(typeof global !== "undefined" && typeof global.gc === "function"),
  };
  console.log(JSON.stringify(summary, null, 2));
}

main().catch(e => { console.error(e); process.exit(1); });
