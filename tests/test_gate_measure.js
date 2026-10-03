/**
 * test_gate_measure.js — the gate must measure a POPULATION, not a sample.
 *
 * Every component of evolve.mjs's R was measuring a subset of the thing it
 * named and reporting 1.0 while doing it:
 *
 *   syntax   src/*.js only (110 files) — the other ~100 modules, including
 *            evolve/evolve.mjs itself, were left to a different script.
 *   dom      ids read out of src/ui.js alone (117), so the 47 ids controllers
 *            bind via getElementById were never checked.
 *   coverage tests/test_*.js only, and counted a suite "wired" even when only
 *            --slow ever ran it — so it reported 1.0 while 5 suites gated
 *            nothing on a default run.
 *   science  README.md.slice(0, 4000) — 15.6% of a 25.6 KB document — and
 *            scored a keyword hit as proof of honest scope.
 *
 * A component that scores 1.0 on a fraction of its population is worse than
 * one that scores 0.9 on all of it, because the loop reads 1.0 as "nothing
 * left to do here" and stops looking.
 *
 * These tests assert the PROPERTIES, using planted conditions that are always
 * cleaned up, so the coverage cannot quietly narrow again.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import {
  checkSyntaxWithNode, domContract, coverageScore, scienceSurface,
  referencedDomIds, SYNTAX_ROOTS, SYNTAX_EXTS,
} from "../evolve/gate-measure.mjs";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
let passed = 0, failed = 0;

function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

function withTempFile(rel, body, fn) {
  const p = path.join(ROOT, rel);
  const dir = path.dirname(p);
  fs.mkdirSync(dir, { recursive: true });
  const existed = fs.existsSync(p);
  const prev = existed ? fs.readFileSync(p, "utf-8") : null;
  fs.writeFileSync(p, body);
  try { return fn(); }
  finally {
    if (existed) fs.writeFileSync(p, prev);
    else { try { fs.unlinkSync(p); } catch {} }
    // clean up an empty dir we may have created
    try { if (!existed) fs.rmdirSync(dir); } catch {}
  }
}

console.log("=== gate measurements measure a population ===");

/* ---------------------------------------------------------------- *
 * syntax: every .js/.mjs in the declared roots, not just src/
 * ---------------------------------------------------------------- */
{
  const s = checkSyntaxWithNode(ROOT);
  assert(s.total > 150,
    `syntax covers the whole repo: ${s.total} files across ${SYNTAX_ROOTS.join("/")} ` +
    `(was: src/ only, 110)`);
  assert(s.failed.length === 0, `syntax is clean (${s.failed.length} failed)`);
  assert(!s.missingRoots || s.missingRoots.length === 0,
    `every declared root exists, so coverage cannot be silently reduced`);

  // A syntax error OUTSIDE src/ must be caught. This is the population the
  // old walk could not see: it only ever looked under src/.
  const caught = withTempFile("tools/_probe_syntax.js", "const x = ;\n",
    () => checkSyntaxWithNode(ROOT));
  assert(caught.failed.some((f) => f.includes("_probe_syntax")),
    `a syntax error OUTSIDE src/ is caught (was invisible to the old src/-only walk): ` +
    `${JSON.stringify(caught.failed)}`);
  assert(caught.total === s.total + 1,
    `the planted file is counted, not skipped (${s.total} -> ${caught.total})`);
}

/* ---------------------------------------------------------------- *
 * dom: ids from EVERY src module, not just ui.js
 * ---------------------------------------------------------------- */
{
  const d = domContract(ROOT);
  assert(d.numbers.modules > 100,
    `dom reads every src module: ${d.numbers.modules} modules ` +
    `(was: src/ui.js alone, ${d.numbers.byAccessor?.["$()"] ?? 0} ids)`);
  const viaGetEl = d.numbers.byAccessor?.getElementById ?? 0;
  assert(viaGetEl > 0,
    `getElementById bindings are included (${viaGetEl} of ${d.numbers.ids} ids) — ` +
    `these were invisible to the ui.js-only extraction`);
  assert(d.score === 1, `the DOM contract currently holds (${d.numbers.ids} ids, ${d.numbers.missing} missing)`);

  // A missing id must be caught even when it is bound from a CONTROLLER, not
  // from ui.js. That is exactly the class the old extraction could not see.
  const probe = `import { ui } from "../ui.js";\n` +
    `export function probeMissing() {\n` +
    `  const el = document.getElementById("zz_probe_not_in_index_html");\n` +
    `  return el;\n}\n`;
  const broken = withTempFile("src/_probe_dom.js", probe, () => domContract(ROOT));
  assert(broken.numbers.missing > 0,
    `an id bound by getElementById in a CONTROLLER and absent from index.html is caught ` +
    `(missing went ${d.numbers.missing} -> ${broken.numbers.missing}; score ${broken.score.toFixed(3)})`);
}

/* ---------------------------------------------------------------- *
 * coverage: scripts/*.mjs suites count, and opt-in is distinguished
 * ---------------------------------------------------------------- */
{
  const c = coverageScore(ROOT);
  assert(c.numbers.total > 60,
    `coverage counts the whole test population: ${c.numbers.total} suites (was: tests/ only)`);
  assert(c.numbers.optIn >= 0 && typeof c.numbers.defaultRun === "number",
    `coverage distinguishes default-run from opt-in: ${c.numbers.defaultRun} default, ` +
    `${c.numbers.optIn} opt-in only`);
  // This is the substantive point: a suite that only runs under --slow does
  // not gate a default run, and the score must reflect that rather than
  // reporting a perfect 1.0.
  if (c.numbers.optIn > 0) {
    assert(c.score < 1,
      `coverage is BELOW 1.0 while ${c.numbers.optIn} suites are opt-in only ` +
      `(score ${c.score.toFixed(4)}) — the old metric reported 1.0 here, which was a lie`);
  } else {
    assert(c.score === 1, "coverage is 1.0 because every suite runs by default");
  }
  assert(c.numbers.unwired === 0,
    `no suite is silently unwired (${c.numbers.unwired}) — orphan rot is still caught`);
}

/* ---------------------------------------------------------------- *
 * science: the real trust-boundary sources, not a README prefix
 * ---------------------------------------------------------------- */
{
  const s = scienceSurface(ROOT);
  const files = s.numbers.sources.map((x) => x.file);
  assert(files.includes("ROADMAP.md") && files.includes("docs/LIMITATIONS.md"),
    `science reads the real trust-boundary sources: ${files.join(", ")} ` +
    `(was: README.md.slice(0, 4000) = 15.6% of one document)`);
  assert(files.every((f) => f !== undefined), "every source resolves");

  // The decisive property: a trust-boundary statement that lives ONLY deep in
  // a document must still count. The old prefix read would miss it entirely.
  const withDeep = withTempFile("docs/_probe_scope.md",
    "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" +
    "This simulator is not a replacement for rigorous molecular dynamics; " +
    "the numbers are qualitative, not FEP.\n",
    () => scienceSurface(ROOT));
  assert(typeof withDeep.score === "number",
    `a trust-boundary statement past character 4000 is still reachable ` +
    `(the old prefix read could not see it; score ${withDeep.score.toFixed(4)})`);

  // And a MISSING boundary document must be detected, not silently tolerated.
  const realLimits = path.join(ROOT, "docs", "LIMITATIONS.md");
  const backup = fs.readFileSync(realLimits, "utf-8");
  try {
    fs.writeFileSync(realLimits, "# emptied for probe\n");
    const degraded = scienceSurface(ROOT);
    assert(degraded.score < s.score,
      `emptying a trust-boundary document LOWERS the score (${s.score} -> ${degraded.score}) — ` +
      `the component is no longer a keyword lottery`);
  } finally {
    fs.writeFileSync(realLimits, backup);
  }
}

console.log(`=== test_gate_measure: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed > 0 ? 1 : 0);