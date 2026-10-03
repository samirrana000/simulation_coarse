/**
 * tests/test_main_module_size.js — main.js must stay a composition root.
 *
 * WHY THIS EXISTS
 * ---------------
 * src/main.js was 1686 LOC: 30 top-level functions and 23 imports covering
 * structure loading, system construction, parameter hot reload, the live-term
 * mirror, dock sparklines, BindViz, session save/restore, the recording
 * panel, hotkeys and the RAF loop. It is the root of the ES-module graph that
 * index.html loads, so one syntax error anywhere in its transitive closure
 * killed the whole app with ZERO terminal output — the server kept serving,
 * the canvas stayed blank, and the only evidence was a SyntaxError in a
 * browser console. That is exactly how this project shipped a blank-canvas
 * failure once (AUDIT_REPORT.md §1.1 / §2.1).
 *
 * The split (2026-10) moved the implementation into src/controllers/, one
 * module per responsibility, leaving main.js as a table of contents: import
 * the modules, call their init functions in the original startup order, start
 * the loop. 1686 -> 159 LOC.
 *
 * A refactor without a guard regrows. Nothing in the toolchain, the physics or
 * the UI pushes code back up; people do. So this test is the whole point.
 *
 * WHAT IT ASSERTS
 *   1. SIZE: main.js <= MAX_MAIN_LOC lines. The threshold is 400 — roughly
 *      2.5x the current file — because the failure mode is "someone pasted a
 *      whole subsystem back in", not "someone added one line". A file at 400
 *      lines is still unmistakably a table of contents.
 *   2. FAN-IN: main.js references at most MAX_MAIN_IMPORTS module specifiers.
 *      This is the metric that grows when a controller is inlined back into
 *      the root instead of being split out. Pre-split value: 23.
 *   3. SHAPE: main.js holds no physics and no canvas painting — it must not
 *      import forcefield.js / heavy.js / integrator.js / funnel.js /
 *      viewer.js / capture/*, and it must contain no `new Funnel(`,
 *      `new ForceField(`, `advance(` or `getContext(`. Each controller
 *      declared in the header comment must exist and be importable, so
 *      "main.js is small" cannot be achieved by deleting the app.
 *
 * FAILING DELIBERATELY (the proof the gate bites)
 *   printf '\n// bloat probe\n'.repeat(120) >> src/main.js  -> FAILS
 *   git checkout src/main.js (or remove the block)         -> passes
 *   tests/test_main_module_size.js prints the measured LOC,
 *   the import count and the FIRST offending line, so a failure
 *   points at the block that was pasted back.
 *
 * FAST tier: reads 13 files, no DOM, no physics, ~30 ms.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const MAIN = path.join(ROOT, "src", "main.js");

/** Hard ceiling on src/main.js. Pre-split: 1686. Post-split: 159. */
const MAX_MAIN_LOC = 400;
/** Hard ceiling on module specifiers referenced by main.js. Pre-split: 23. */
const MAX_MAIN_IMPORTS = 20;

/** Modules the composition root must never import (physics / rendering). */
const FORBIDDEN_IN_MAIN = [
  "forcefield.js", "heavy.js", "integrator.js", "funnel.js", "viewer.js",
  "recorder.js", "capture/bindlog.js", "heavy_progress.js",
];
/** Syntax that only belongs inside a controller, never in the root. */
const FORBIDDEN_TOKENS = [
  "new ForceField(", "new HeavyForceField(", "new Funnel(", "new LangevinIntegrator(",
  ".advance(", "getContext(", "addEventListener(",
];
/** Every controller the header table promises must exist and parse. */
const REQUIRED_CONTROLLERS = [
  "physics-tier.js", "param-binding.js", "live-terms.js", "dock.js",
  "binding-insights.js", "accelerate.js", "structure-input.js",
  "system-build.js", "transport.js", "recording.js", "guide.js", "tick.js",
];

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

/** Every module specifier main.js names, via `import` or `export ... from`. */
const specs0 = (text) => [...text.matchAll(/(?:^|\n)\s*(?:import|export)\s+(?:[^;"']*?\s+from\s+)?["']([^"']+)["']/g)].map((m) => m[1]);

const src = fs.readFileSync(MAIN, "utf-8");
const lines = src.split("\n");
// evolve/evolve.mjs counts `split("\n").length`, which is one more than
// `wc -l` on a newline-terminated file. Both are printed so a number quoted in
// a doc can be matched to the counter that produced it.
const loc = lines.length;
const wcLoc = loc - (lines[loc - 1] === "" ? 1 : 0);

console.log("=== main.js composition-root guard ===");
console.log(`      measured: ${loc} LOC by the gate's counter (${wcLoc} by wc -l), limit ${MAX_MAIN_LOC}` +
  `; ${specs0(src).length} module specifiers (limit ${MAX_MAIN_IMPORTS})`);

// ---- 1. SIZE ------------------------------------------------------------
assert(loc <= MAX_MAIN_LOC,
  `main.js is ${loc} LOC (<= ${MAX_MAIN_LOC}; was 1686 before the controllers/ split)`);

// ---- 2. FAN-IN ----------------------------------------------------------
const specs = specs0(src);
const uniqueSpecs = [...new Set(specs)];
assert(specs.length <= MAX_MAIN_IMPORTS,
  `main.js references ${specs.length} module specifiers (<= ${MAX_MAIN_IMPORTS}; was 23)`);

// ---- 3. SHAPE -----------------------------------------------------------
for (const bad of FORBIDDEN_IN_MAIN) {
  assert(!specs.includes(`./${bad}`) && !specs.includes(`./${bad}`),
    `main.js does not import ${bad} (the root must not own physics or rendering)`);
}
for (const tok of FORBIDDEN_TOKENS) {
  // Strip comments so prose about `advance(` in the header cannot trip it.
  const code = src.replace(/\/\*[\s\S]*?\*\//g, "").replace(/^\s*\/\/.*$/gm, "");
  const line = code.split("\n").findIndex((l) => l.includes(tok)) + 1;
  assert(line === 0, `main.js contains no \`${tok}\`` + (line ? ` (line ${line})` : ""));
}
for (const c of REQUIRED_CONTROLLERS) {
  const p = path.join(ROOT, "src", "controllers", c);
  const exists = fs.existsSync(p) && fs.readFileSync(p, "utf-8").trim().length > 0;
  assert(exists, `src/controllers/${c} exists and is non-empty`);
  if (exists) {
    const rel = `./controllers/${c}`;
    assert(uniqueSpecs.includes(rel),
      `main.js imports ./controllers/${c} (the header table is not a lie)`);
  }
}

// ---- 4. the composition-root contract is textual and load-bearing --------
assert(/startTick\(\)/.test(src), "main.js starts the animation loop (startTick)");
const inits = [...src.matchAll(/^init[A-Z]\w*\(/gm)].map((m) => m[0]);
assert(inits.length >= 8,
  `main.js calls ${inits.length} controller init functions (${inits.join(" ")})`);
assert(/import\s+"\.\/analysis-panel\.js"/.test(src),
  "main.js still imports analysis-panel.js for its side effects");
assert(!/^\s*(function|const)\s+(tick|buildSystem|onParamChange)\b/m.test(src),
  "main.js defines none of tick/buildSystem/onParamChange itself");

console.log(`\n${passed} PASSED, ${failed} FAILED`);
process.exit(failed > 0 ? 1 : 0);