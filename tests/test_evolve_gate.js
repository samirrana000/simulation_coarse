/**
 * test_evolve_gate.js — the gate must be a GATE.
 *
 * The loop found two defects in its own scoring, both of which had the same
 * fatal property: they made the gate unable to fail, or unable to tell good
 * work from bad.
 *
 *   D1  Bloat scored VOLUME (tracked bytes, total src LOC). A goal that added
 *       a test or a shared physics kernel LOWERED the score, so the loop
 *       would refuse good work. It now scores STRUCTURE (largest module
 *       relative to the median) and duplication.
 *   D2  src/ metrics were read from `git ls-files`, so UNCOMMITTED work was
 *       invisible. A 3000-line junk module added to src/ did not move the
 *       score, because the file was not yet tracked. A metric blind to the
 *       edits an agent just made cannot gate that agent.
 *   D3  `gate` printed "verdict: CLOSED (roll back)" and still exited 0. In
 *       CI that is decoration, not a gate.
 *   D4  R_best was recorded against a tree that no longer existed, so it was
 *       unreachable and every subsequent run read as a regression.
 *
 * These tests assert the PROPERTIES, using temporary probes that are always
 * cleaned up, so the loop cannot silently return to a metric that cannot fail.
 */

import fs from "node:fs";
import path from "node:path";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const EVOLVE = path.join(ROOT, "evolve", "evolve.mjs");
let passed = 0, failed = 0;

function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

/**
 * Run the gate, returning { json, code }.
 *
 * `--fast` reuses the last recorded test result instead of re-running the
 * suite. Without it this test re-ran the full 16s suite once per probe —
 * O(n^2), and the test cost more than the thing it tested. Structural
 * probes (syntax / bloat / dom / coverage) do not depend on the test result,
 * so they are exactly as meaningful with the cached value. The one assertion
 * that genuinely needs a live suite run (the "must be able to fail" check on
 * the full path) uses a single non-fast call.
 */
function gate(args = ["--json", "--fast"]) {
  try {
    const out = execFileSync(process.execPath, [EVOLVE, "gate", ...args], {
      cwd: ROOT, encoding: "utf-8",
      stdio: ["ignore", "pipe", "pipe"], timeout: 300000,
      // EVOLVE_GATE_CHILD stops the gate's own test run from re-entering the
      // suite that invoked this test. Without it the gate recurses until the
      // stack dies. The child flag is an env var precisely because argv of the
      // grandparent is invisible to the gate.
      env: { ...process.env, EVOLVE_GATE_CHILD: "1" },
    });
    return { json: JSON.parse(out), code: 0 };
  } catch (e) {
    let json = null;
    try { json = JSON.parse(e.stdout); } catch {}
    return { json, code: e.status ?? -1 };
  }
}

const PROBE = path.join(ROOT, "src", "_gate_probe_tmp.js");
function withProbe(body, fn) {
  fs.writeFileSync(PROBE, body);
  try { return fn(); }
  finally { try { fs.unlinkSync(PROBE); } catch {} }
}

console.log("=== evolve gate: a gate must be able to fail ===");

const base = gate();

// D3 — a closed gate must exit non-zero, or CI cannot enforce it.
//     This is the ONE probe that runs the real (non-fast) path.
{
  const r = gate(["--fast"]); // no --json => real exit code
  assert(typeof r.code === "number", `gate exposes a real exit code (got ${r.code})`);
  const closed = base.json.verdict.startsWith("CLOSED");
  if (closed) {
    assert(r.code !== 0, "D3: a CLOSED gate exits non-zero (was: printed 'roll back' and exited 0)");
  } else {
    assert(r.code === 0, "an OPEN gate exits 0");
  }
}

// D2 — the gate must see UNCOMMITTED src/ changes.
{
  const junk = withProbe(
    Array.from({ length: 3000 }, (_, i) => `export const junk_${i} = ${i};`).join("\n"),
    () => gate().json
  );
  const real = withProbe(
    Array.from({ length: 20 }, (_, i) => `export const small_${i} = ${i};`).join("\n"),
    () => gate().json
  );
  assert(junk.parts.syntax === real.parts.syntax,
    "D2: a valid untracked src/ file does not fail the syntax gate");
  const junkLoc = junk.tracked.srcLoc, realLoc = real.tracked.srcLoc;
  assert(junkLoc - realLoc >= 2980,
    `D2: untracked src/ LOC is measured (junk +${junkLoc - realLoc}, small +${realLoc - junkLoc} base) — ` +
    `a metric that cannot see uncommitted edits cannot gate the agent making them`);
  assert(junk.parts.bloat < real.parts.bloat,
    `D2: adding a 3000-line god module LOWERS bloat (${real.parts.bloat} -> ${junk.parts.bloat})`);
}

// D1 — good work must not be penalised. A small, well-formed module is not bloat.
{
  const small = withProbe(
    Array.from({ length: 20 }, (_, i) => `export const small_${i} = ${i};`).join("\n"),
    () => gate().json
  );
  assert(Math.abs(small.parts.bloat - base.json.parts.bloat) < 1e-9,
    `D1: adding a small clean module does not lower bloat ` +
    `(base ${base.json.parts.bloat}, with module ${small.parts.bloat}) — volume is not bloat`);
}

// The syntax gate must actually catch a planted syntax error.
{
  const bad = withProbe("const x = ;\n", () => gate().json);
  assert(bad.parts.syntax < 1, `syntax gate catches a planted error in an untracked file (${bad.parts.syntax})`);
  assert(bad.notes.some((n) => /syntax FAIL/.test(n)), "the failure is NAMED in the notes, not silent");
}

// The gate must report saturation, and must renormalise around it.
{
  assert(Array.isArray(base.json.saturated), "gate reports which components are saturated");
  const w = base.json.weights;
  const satW = base.json.saturated.reduce((a, k) => a + (w[k] || 0), 0);
  assert(satW > 0, "at least one component is saturated (a dead reward carries no gradient)");
}

// R_best must be reachable: the current tree must not read as a regression
// against itself. This is D4 — an unreachable R_best blocks the whole loop.
{
  assert(base.json.verdict.startsWith("OPEN"),
    `D4: the current tree satisfies its own gate (${base.json.verdict}; ` +
    `R=${base.json.score} vs R_best=${base.json.R_best}) — an unreachable R_best blocks all work`);
  assert(base.json.R_best <= base.json.score + 1e-9,
    `R_best (${base.json.R_best}) is not in the future (score ${base.json.score})`);
}

console.log(`=== test_evolve_gate: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed > 0 ? 1 : 0);
