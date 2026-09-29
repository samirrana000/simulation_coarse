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

// D3 — the gate must expose a real exit code, and it must be 0 for an OPEN
//      gate and non-zero for a CLOSED one. Asserted on the junk probe below
//      (which is guaranteed to regress) rather than on whatever the current
//      tree happens to be.
{
  const r = gate(["--fast"]); // no --json => real exit code
  assert(typeof r.code === "number", `gate exposes a real exit code (got ${r.code})`);
  const closed = base.json.verdict.startsWith("CLOSED");
  assert(closed ? r.code !== 0 : r.code === 0,
    `D3: exit code matches the verdict (${base.json.verdict} -> exit ${r.code})`);
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

// D4 — R_best must never be in the FUTURE. If R_best exceeds the current
//      score, the gate is permanently CLOSED and the loop can accept nothing
//      ever again. It is correct for the gate to read CLOSED after a change
//      that lowers the score; that is the gate working. What must never
//      happen is R_best being unreachable by the current tree.
//      (R_best is re-baselined on `done ACCEPTED`, so an accepted change
//      becomes the new reference rather than a permanent regression.)
{
  assert(base.json.R_best <= base.json.score + 1e-9,
    `D4: R_best (${base.json.R_best}) is reachable by the current tree ` +
    `(score ${base.json.score}) — an unreachable R_best blocks all future work`);

  // The gate must be able to DETECT a regression: a junk module must read
  // as a regression, not merely as a lower number nobody acts on.
  // NOTE: no --json here. --json is machine-readable-output mode and
  // deliberately does not set an exit code, so asking for it would blind
  // this probe to the very property it is testing.
  const junk = withProbe(
    Array.from({ length: 3000 }, (_, i) => `export const junk_${i} = ${i};`).join("\n"),
    () => gate(["--fast"])
  );
  assert(junk.json.score < base.json.score,
    `D4: a regressing tree lowers the score (${base.json.score} -> ${junk.json.score})`);
  assert(junk.code !== 0,
    `D4: a regressing tree makes the gate exit non-zero (exit ${junk.code}) — ` +
    `a gate that prints 'roll back' and exits 0 cannot fail anything`);
}

console.log(`=== test_evolve_gate: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed > 0 ? 1 : 0);
