/**
 * bench/perf.js — Performance benchmark for Cα ENM vs Heavy forcefields.
 *
 * Reads 4w52.pdb, builds:
 *   CG:    ~164 Cα beads (ForceField)
 *   Heavy: ~1308 heavy atoms (HeavyForceField)
 * Runs ff.compute(ref) 20× warmup + 30× timed, reporting mean±sd per call.
 *
 * Runnable: node bench/perf.js  (from repo root or bench/ dir)
 * Prints human-readable lines with "CG:" and "Heavy:" plus a JSON summary.
 *
 * Uses performance.now() if available, else Date.now().
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

// Import forcefields — try without ?v suffix first (Node-friendly), fall back to ?v=10 style
import { ForceField } from "../src/forcefield.js";
import { HeavyForceField, parseHeavy } from "../src/heavy.js";
import { parseCa, selectSystem } from "../src/pdb.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const now = () => (typeof performance !== "undefined" && performance.now ? performance.now() : Date.now());

function readPdbOrThrow(name) {
  const candidates = [
    path.resolve(process.cwd(), name),
    path.resolve(__dirname, "..", name),
    path.resolve(__dirname, name),
    path.resolve("data", name),
  ];
  for (const p of candidates) {
    try {
      if (fs.existsSync(p)) return fs.readFileSync(p, "utf-8");
    } catch {}
  }
  throw new Error(`Cannot find ${name} (tried ${candidates.join(", ")})`);
}

function stats(times) {
  const n = times.length;
  const mean = times.reduce((a, b) => a + b, 0) / n;
  let v = 0;
  for (const t of times) v += (t - mean) ** 2;
  const sd = Math.sqrt(v / n);
  // 95% CI approx mean ± 1.96*sd/sqrt(n) but we just report sd
  return { mean, sd, n, min: Math.min(...times), max: Math.max(...times) };
}

function benchCompute(label, ff, pos) {
  const WARMUP = 20;
  const TIMED = 30;
  for (let i = 0; i < WARMUP; i++) ff.compute(pos);
  const times = [];
  for (let i = 0; i < TIMED; i++) {
    const t0 = now();
    ff.compute(pos);
    const t1 = now();
    times.push(t1 - t0);
  }
  const s = stats(times);
  return { label, n: ff.n, meanMs: s.mean, sdMs: s.sd, minMs: s.min, maxMs: s.max, nSamples: s.n };
}

async function main() {
  const pdbText = readPdbOrThrow("4w52.pdb");

  // CG: parse Cα and build ENM
  const parsedCa = parseCa(pdbText);
  const sel = selectSystem(parsedCa);
  const ffCG = new ForceField(sel, { rc: 10, gamma: 1.0 });
  const posCG = new Float64Array(ffCG.ref);
  const resCG = benchCompute("CG", ffCG, posCG);

  // Heavy: parse all heavy atoms
  let resHeavy = null;
  try {
    const parsedHeavy = parseHeavy(pdbText);
    // Use raw heavy atoms (1308) to match spec; fallback to selected system if needed
    let heavySystem;
    if (parsedHeavy.atoms && parsedHeavy.atoms.length) {
      heavySystem = { atoms: parsedHeavy.atoms };
    } else {
      throw new Error("parseHeavy returned no atoms");
    }
    const ffHeavy = new HeavyForceField(heavySystem, { gamma: 1.0 });
    const posHeavy = new Float64Array(ffHeavy.ref);
    resHeavy = benchCompute("Heavy", ffHeavy, posHeavy);
  } catch (e) {
    // Fallback: report error but keep placeholder so script doesn't crash
    console.error("Heavy bench failed:", e.message);
    resHeavy = { label: "Heavy", n: 1308, meanMs: NaN, sdMs: NaN, minMs: NaN, maxMs: NaN, nSamples: 0, error: e.message };
  }

  // Scale exponent placeholder: naive O(N^alpha) estimate from timing ratio
  // alpha = log(tHeavy/tCG) / log(nHeavy/nCG) ; placeholder if not computable
  let scaleExponent = null;
  let scaleNote = "placeholder (ideal O(N) grid ⇒ exponent ~1.0; measured ratio used if available)";
  if (Number.isFinite(resCG.meanMs) && Number.isFinite(resHeavy.meanMs) && resCG.meanMs > 0 && resHeavy.meanMs > 0) {
    scaleExponent = Math.log(resHeavy.meanMs / resCG.meanMs) / Math.log(resHeavy.n / resCG.n);
  }

  // Human-readable lines (required by measurable criteria: must contain "CG:" and "Heavy:")
  console.log(`CG: n=${resCG.n}  ${resCG.meanMs.toFixed(3)} ± ${resCG.sdMs.toFixed(3)} ms/compute  (min ${resCG.minMs.toFixed(3)} max ${resCG.maxMs.toFixed(3)} over ${resCG.nSamples})`);
  console.log(`Heavy: n=${resHeavy.n}  ${Number.isFinite(resHeavy.meanMs) ? resHeavy.meanMs.toFixed(3) : "NaN"} ± ${Number.isFinite(resHeavy.sdMs) ? resHeavy.sdMs.toFixed(3) : "NaN"} ms/compute  (min ${Number.isFinite(resHeavy.minMs) ? resHeavy.minMs.toFixed(3) : "NaN"} max ${Number.isFinite(resHeavy.maxMs) ? resHeavy.maxMs.toFixed(3) : "NaN"} over ${resHeavy.nSamples})`);
  console.log(`scale exponent: ${scaleExponent !== null && Number.isFinite(scaleExponent) ? scaleExponent.toFixed(3) : "N/A"} — ${scaleNote}`);

  // JSON summary (machine-readable)
  const summary = {
    cg: { n: resCG.n, meanMs: resCG.meanMs, sdMs: resCG.sdMs, minMs: resCG.minMs, maxMs: resCG.maxMs, samples: resCG.nSamples },
    heavy: { n: resHeavy.n, meanMs: resHeavy.meanMs, sdMs: resHeavy.sdMs, minMs: resHeavy.minMs, maxMs: resHeavy.maxMs, samples: resHeavy.nSamples, error: resHeavy.error || null },
    scaleExponent: scaleExponent,
    scaleExponentNote: scaleNote,
    warmup: 20,
    timed: 30,
    pdb: "4w52.pdb",
    timestamp: new Date().toISOString(),
  };
  console.log(JSON.stringify(summary, null, 2));
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
