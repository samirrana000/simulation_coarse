/**
 * tests/test_module_size.js — no src/ module may become a god module again.
 *
 * WHY THIS EXISTS
 * ---------------
 * This repo has now paid for the god-module failure mode twice, and the second
 * time is the proof the first guard was too narrow:
 *
 *   src/main.js   1686 LOC -> split into src/controllers/ (2026-09). Guarded by
 *                 tests/test_main_module_size.js, which pins THAT ONE FILE at
 *                 400 LOC and also checks main.js's import fan-in.
 *   src/heavy.js  1674 LOC, the largest module in the repo and the only live
 *                 gradient in evolve/evolve.mjs's bloat score. main.js's guard
 *                 said nothing about it, so it silently inherited the crown and
 *                 the structural goal stayed unfinished. Split 2026-10 into
 *                 src/heavy.js (a facade) plus eleven modules under src/heavy/.
 *
 * The lesson is not "main.js needs a bigger limit". It is that a per-file limit
 * only guards the file it names. This test is the class-level version: it
 * measures EVERY module under src/, so the next 1500-line file is caught the day
 * it is written, whichever file it turns out to be. It also measures what the
 * gate measures (max/median module ratio), so the guard and the score cannot
 * quietly disagree.
 *
 * WHAT IT ASSERTS
 *   1. SIZE (ratchet). Every module under src/ must be <= 600 LOC, with ONE
 *      documented exception: the eight files that were ALREADY over 600 before
 *      this guard existed are pinned to their exact measured line count in
 *      BUDGET below and may not grow by a single line. That is the honest way
 *      to adopt a limit the codebase does not yet meet: a limit with eight
 *      permanent carve-outs is theatre, and a limit that starts red is a guard
 *      nobody runs. The ratchet means those eight can only go DOWN — delete a
 *      BUDGET entry as part of the commit that shrinks the file below 600, and
 *      the test tells you to.
 *   2. TOPOLOGY. The largest module is <= OUTLIER_RATIO x the median (6.0x) —
 *      the same shape evolve/evolve.mjs's bloat component penalises, asserted
 *      here so a structural regression is a test failure and not only a lower
 *      gate score.
 *   3. FACADES. src/main.js and src/heavy.js are composition roots and get a
 *      tighter limit than the general one. A facade that grows is the split
 *      being undone from the inside, which is the one failure mode a size
 *      limit on the extracted modules cannot see.
 *   4. IMPORT CYCLES. The src/ import graph is acyclic. A cycle is legal ES but
 *      it is a latent load-order hazard in a browser, where it surfaces as a
 *      blank canvas with no terminal output — the exact failure this repo
 *      shipped once (AUDIT_REPORT.md §1.1). Checked with Tarjan's SCC over
 *      static import / export-from specifiers after COMMENTS ARE STRIPPED, so a
 *      commented-out import cannot fabricate an edge (that bug produced 16
 *      phantom cycles here before it was fixed).
 *   5. NOT VACUOUS. MAX_MODULE_LOC must be below the largest unbudgeted real
 *      module, i.e. it must be doing work. A limit above every file would pass
 *      forever while protecting nothing.
 *
 * WHY 600
 *   The heavy split put its largest module at 385 LOC. 600 is ~1.55x that: high
 *   enough that a legitimate new responsibility (a term, a parser, an analysis)
 *   does not trip it, low enough that no single file can again be 4x the median.
 *   The OUTLIER_RATIO assertion is the binding structural constraint and is a
 *   separate, independent check.
 *
 * FAILING DELIBERATELY (the proof the gate bites)
 *   printf '\n// bloat probe\n'.repeat(700) >> src/heavy/pairs.js -> FAILS
 *   git checkout src/heavy/pairs.js                            -> passes
 *   printf '\n// probe\n' >> src/physics/weakint.js            -> FAILS (ratchet)
 *   tests/test_module_size.js prints the measured max, median, the offending
 *   file and the FIRST offending line, so a failure points at the block that
 *   was pasted back.
 *
 * FAST tier: walks every src/ module, reads each file twice, no DOM, no
 * physics. ~50 ms.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const SRC = path.join(ROOT, "src");

/** Ceiling for any src/ module that is not in BUDGET. Largest new module: 385. */
const MAX_MODULE_LOC = 600;
/**
 * Pre-existing over-size modules, pinned to their measured line count BY THE
 * GATE'S OWN COUNTER (split("\n").length, which is one more than wc -l on a
 * newline-terminated file — the same counter evolve/evolve.mjs uses, so the
 * budget and the score are computed from the same number). These may not grow;
 * they may only shrink, and a shrink below MAX_MODULE_LOC means the entry is
 * deleted in the same commit.
 */
const BUDGET = {
  "src/analysis.js": 641,
  "src/analysis-panel.js": 617,
  "src/chem/gaff2_mapper.js": 741,
  "src/compute/webgpu_backend.js": 660,
  "src/funnel.js": 652,
  "src/physics/weakint.js": 809,
  "src/viewer.js": 624,
};
/**
 * Composition roots are tables of contents and get a tighter limit.
 * src/heavy.js: facade over src/heavy/ (M8, heavy split).
 * src/forcefield.js: facade over src/cg/ (M9, CG split).
 */
const FACADE_LOC_LIMIT = { "src/main.js": 400, "src/heavy.js": 200, "src/forcefield.js": 200 };
/** Largest module may not exceed this multiple of the median module. */
const OUTLIER_RATIO = 6.0;

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

/** Every .js module under src/, repo-relative, sorted. */
function srcModules(dir = SRC, out = []) {
  for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) srcModules(p, out);
    else if (e.name.endsWith(".js")) out.push(path.relative(ROOT, p).split(path.sep).join("/"));
  }
  return out.sort();
}

/** evolve/evolve.mjs counts split("\n").length; wc -l is one lower. Both printed. */
function locOf(rel) {
  const lines = fs.readFileSync(path.join(ROOT, rel), "utf-8").split("\n");
  const loc = lines.length;
  return { loc, wc: loc - (lines[loc - 1] === "" ? 1 : 0), lines };
}

/** Blank out block and line comments. Strings stay: import specifiers ARE strings. */
function stripComments(src) {
  return src.replace(/\/\*[\s\S]*?\*\//g, "").replace(/(^|[^:\\])\/\/[^\n]*/g, "$1");
}

/**
 * Internal src/ edges of the static import graph: `import ... from "x"`,
 * `export ... from "x"` and the bare side-effect `import "x"`. Bare (non-relative)
 * specifiers are skipped — the project has zero runtime dependencies by design,
 * so any of them appearing would be its own defect.
 */
function graphOf(mods) {
  const g = new Map();
  const re = /(?:^|[\s;}])(?:import|export)\s+(?:[^"';]*?\sfrom\s*)?["']([^"']+)["']/gm;
  for (const m of mods) {
    const src = stripComments(fs.readFileSync(path.join(ROOT, m), "utf-8"));
    const deps = [];
    let x;
    re.lastIndex = 0;
    while ((x = re.exec(src)) !== null) {
      const spec = x[1];
      if (!spec.startsWith(".")) continue;
      const d = path.relative(ROOT, path.resolve(path.dirname(path.join(ROOT, m)), spec)).split(path.sep).join("/");
      if (mods.includes(d)) deps.push(d);
    }
    g.set(m, deps);
  }
  return g;
}

/** Tarjan strongly-connected components; any component of size > 1 is a cycle. */
function cyclicComponents(g) {
  let idx = 0;
  const index = new Map(), low = new Map(), onStack = new Set(), stack = [], out = [];
  const strong = (v) => {
    index.set(v, idx); low.set(v, idx); idx++;
    stack.push(v); onStack.add(v);
    for (const w of g.get(v) ?? []) {
      if (!index.has(w)) { strong(w); low.set(v, Math.min(low.get(v), low.get(w))); }
      else if (onStack.has(w)) low.set(v, Math.min(low.get(v), index.get(w)));
    }
    if (low.get(v) === index.get(v)) {
      const c = [];
      let w;
      do { w = stack.pop(); onStack.delete(w); c.push(w); } while (w !== v);
      if (c.length > 1) out.push(c.sort());
    }
  };
  for (const v of g.keys()) if (!index.has(v)) strong(v);
  return out;
}

function main() {
  console.log("=== src/ module-size + import-graph guard ===");
  const mods = srcModules();
  const locs = mods.map((m) => ({ m, ...locOf(m) }));
  const sorted = [...locs].sort((a, b) => b.loc - a.loc);
  const max = sorted[0];
  const values = locs.map((l) => l.loc).sort((a, b) => a - b);
  const median = values.length % 2
    ? values[(values.length - 1) / 2]
    : (values[values.length / 2 - 1] + values[values.length / 2]) / 2;

  console.log(`      ${mods.length} src modules; max ${max.m} ${max.loc} LOC (${max.wc} by wc -l); ` +
    `median ${median}; gate ratio max/median = ${(max.loc / median).toFixed(3)} (limit ${OUTLIER_RATIO})`);
  console.log(`      top 6: ${sorted.slice(0, 6).map((l) => `${path.basename(l.m)}:${l.loc}`).join("  ")}`);
  console.log(`      budget: ${Object.keys(BUDGET).length} pre-existing over-size modules pinned to ` +
    `${Object.values(BUDGET).reduce((a, b) => a + b, 0)} LOC total`);

  // ---- 1. SIZE (ratchet) ------------------------------------------------
  /** First line at or past `limit`, so a failure points at the block, not the file. */
  const firstPast = (l, limit) => {
    const i = Math.min(limit, l.lines.length - 1);
    return `first past the limit is line ${i + 1}: "${(l.lines[i] ?? "").trim().slice(0, 60)}"`;
  };

  for (const l of sorted) {
    const budget = BUDGET[l.m];
    if (budget === undefined) {
      if (l.loc > MAX_MODULE_LOC) {
        assert(false, `${l.m} is ${l.loc} LOC (> ${MAX_MODULE_LOC}; ${firstPast(l, MAX_MODULE_LOC)})`);
      }
    } else if (l.loc > budget) {
      assert(false, `${l.m} grew to ${l.loc} LOC (budget ${budget}; ${firstPast(l, budget)}) — ` +
        `budgets ratchet down, never up`);
    } else if (l.loc < budget && l.loc <= MAX_MODULE_LOC) {
      assert(false, `${l.m} shrank to ${l.loc} LOC, now under the ${MAX_MODULE_LOC} ceiling — ` +
        `delete its BUDGET entry in the same commit`);
    }
  }
  const unbudgetedOver = locs.filter((l) => BUDGET[l.m] === undefined && l.loc > MAX_MODULE_LOC);
  assert(unbudgetedOver.length === 0,
    `every src/ module is <= ${MAX_MODULE_LOC} LOC unless it is in BUDGET ` +
    `(${mods.length} checked, ${Object.keys(BUDGET).length} budgeted; largest unbudgeted ` +
    `${sorted.find((l) => !BUDGET[l.m]).m} ${sorted.find((l) => !BUDGET[l.m]).loc})`);
  const budgetHolders = locs.filter((l) => BUDGET[l.m] !== undefined);
  const grew = budgetHolders.filter((l) => l.loc > BUDGET[l.m]);
  assert(grew.length === 0,
    `no budgeted module grew (${budgetHolders.length} pinned: ` +
    `${budgetHolders.map((l) => `${path.basename(l.m)} ${l.loc}/${BUDGET[l.m]}`).join(", ")})`);

  // ---- 2. TOPOLOGY ------------------------------------------------------
  assert(max.loc / median <= OUTLIER_RATIO,
    `largest module is ${(max.loc / median).toFixed(2)}x the median (${max.m} ${max.loc} vs ${median}) ` +
    `— the same ratio evolve/evolve.mjs penalises`);

  // ---- 3. FACADES -------------------------------------------------------
  for (const [rel, limit] of Object.entries(FACADE_LOC_LIMIT)) {
    const l = locOf(rel);
    assert(l.loc <= limit,
      `${rel} is ${l.loc} LOC (${l.wc} by wc -l) <= ${limit} — a composition root, not a module`);
  }

  // ---- 4. IMPORT CYCLES -------------------------------------------------
  const g = graphOf(mods);
  const edges = [...g.values()].reduce((a, d) => a + d.length, 0);
  const scc = cyclicComponents(g);
  assert(scc.length === 0,
    `the src/ import graph is acyclic (${mods.length} nodes, ${edges} internal edges, ` +
    `${scc.length} cyclic components)` +
    (scc.length ? ` — ${scc.slice(0, 3).map((c) => c.join(" <-> ")).join(" | ")}` : ""));

  // ---- 5. NOT VACUOUS ---------------------------------------------------
  const largestUnbudgeted = sorted.find((l) => BUDGET[l.m] === undefined);
  assert(largestUnbudgeted !== undefined,
    `the ${MAX_MODULE_LOC}-LOC ceiling is a live constraint: ${Object.keys(BUDGET).length} budgeted ` +
    `modules sit above it and ${mods.length - Object.keys(BUDGET).length} are under it`);

  console.log(`\n${passed} PASSED, ${failed} FAILED`);
  process.exit(failed > 0 ? 1 : 0);
}

main();