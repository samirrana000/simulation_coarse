#!/usr/bin/env node
/**
 * tests/test_suite_registry.js — anti-rot self-check for the test registry.
 *
 * THE DEFECT THIS GUARDS AGAINST
 * -------------------------------
 * A test file that exists but never runs is a false claim of coverage: it
 * looks rigorous and gates nothing. Before this file existed, tests/ held
 * ~30 test_*.js scripts that no tier, CI step, or gate ever executed, so
 * they could rot silently and forever (two of them had already been deleted
 * out from under the repo while still "present" in the working tree).
 *
 * THE CONTRACT
 * ------------
 * Every tests/test_*.js file must be EITHER
 *   (a) registered in tests/suites.js under some tier, OR
 *   (b) present under tests/manual/ (with a README explaining how to run it),
 *   (c) matched by an EXEMPT entry that carries a written reason.
 * The mirror rule holds for scripts/test_*.mjs.
 *
 * The registry is also checked for its own integrity: no `expect:` assertion
 * counts (those are derived at run time now), no duplicate file entries, no
 * entries pointing at files that do not exist, and only known tier names.
 *
 * Runnable standalone: node tests/test_suite_registry.js
 * Runs as part of the FAST tier: it is itself registered in tests/suites.js.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { SUITES, EXEMPT, TIERS, FAST_BUDGET_S, suitesIn } from "./suites.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..");

let passed = 0;
let failed = 0;

function assert(condition, message) {
  if (condition) {
    passed++;
    console.log(`  ✓ ${message}`);
  } else {
    failed++;
    console.error(`  ✗ FAIL: ${message}`);
  }
}

/** Repo-relative, "/" separators. */
function rel(p) {
  return path.relative(ROOT, p).split(path.sep).join("/");
}

function listDir(dir, re) {
  if (!fs.existsSync(dir)) return [];
  return fs.readdirSync(dir, { withFileTypes: true })
    .filter((d) => d.isFile() && re.test(d.name))
    .map((d) => rel(path.join(dir, d.name)))
    .sort();
}

/**
 * Recursive relative paths of TEST files under tests/manual/. The README is
 * infrastructure, not a test, so it is excluded here (it is checked for
 * separately in section 2).
 */
function manualFiles(dir = path.join(HERE, "manual")) {
  if (!fs.existsSync(dir)) return [];
  const out = [];
  for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) out.push(...manualFiles(p));
    else if (/\.m?js$/.test(e.name) && e.name !== "README.md") out.push(rel(p));
  }
  return out.sort();
}

function exemptFor(relPath) {
  return EXEMPT.find((x) => x.match.test(relPath));
}

console.log("=== Suite registry anti-rot self-check ===");

// ---------------------------------------------------------------------
// 1. Every test file on disk is accounted for.
// ---------------------------------------------------------------------
const registered = new Set(SUITES.map((s) => s.file));
const manual = new Set(manualFiles());
const manualBasenames = new Set([...manual].map((f) => path.basename(f)));

const subjectDirs = [
  { dir: path.join(HERE), re: /^test_.*\.js$/, label: "tests/test_*.js" },
  { dir: path.join(ROOT, "scripts"), re: /^test_.*\.mjs$/, label: "scripts/test_*.mjs" },
];

for (const { dir, re, label } of subjectDirs) {
  const files = listDir(dir, re).filter((f) => f !== rel(path.join(HERE, "test_all.js")));
  const orphans = files.filter((f) => !registered.has(f) && !exemptFor(f) && !manualBasenames.has(path.basename(f)));
  assert(
    orphans.length === 0,
    `${label}: no unwired files (${files.length} found, ${SUITES.length} registered, ${manualBasenames.size} in tests/manual/)` +
    (orphans.length ? ` — ORPHANS: ${orphans.join(", ")}` : "")
  );
}

// ---------------------------------------------------------------------
// 2. tests/manual/ has a README and a non-empty justification per file.
// ---------------------------------------------------------------------
const manualDir = path.join(HERE, "manual");
if (fs.existsSync(manualDir) && manualFiles().length) {
  const readme = path.join(manualDir, "README.md");
  assert(fs.existsSync(readme), "tests/manual/ has a README.md explaining what lives there and how to run it");
  const body = fs.existsSync(readme) ? fs.readFileSync(readme, "utf-8") : "";
  for (const f of manualFiles()) {
    const base = path.basename(f);
    assert(body.includes(base), `tests/manual/README.md documents ${base}`);
  }
} else {
  passed++;
  console.log("  ✓ tests/manual/ is empty or absent (nothing to justify)");
}

// ---------------------------------------------------------------------
// 3. Registry integrity: no stale assertion counts, no dupes, no ghosts.
// ---------------------------------------------------------------------
const countFields = SUITES.filter((s) => "expect" in s || "asserts" in s);
assert(
  countFields.length === 0,
  `registry declares no hardcoded assertion counts (found ${countFields.length}; counts are derived from child output)` +
  (countFields.length ? ` — ${countFields.map((s) => s.file).join(", ")}` : "")
);

const seen = new Map();
const dupes = [];
for (const s of SUITES) seen.set(s.file, (seen.get(s.file) || 0) + 1);
for (const [f, n] of seen) if (n > 1) dupes.push(`${f} x${n}`);
assert(dupes.length === 0, `no duplicate registry entries (${seen.size} unique of ${SUITES.length})` + (dupes.length ? ` — ${dupes.join(", ")}` : ""));

const ghosts = SUITES.filter((s) => !fs.existsSync(path.join(ROOT, s.file)));
assert(ghosts.length === 0, `every registered suite exists on disk (${SUITES.length})` + (ghosts.length ? ` — MISSING: ${ghosts.map((s) => s.file).join(", ")}` : ""));

const badTier = SUITES.filter((s) => !TIERS.includes(s.tier));
assert(badTier.length === 0, `every suite declares a known tier (${TIERS.join("/")})` + (badTier.length ? ` — ${badTier.map((s) => `${s.file}:${s.tier}`).join(", ")}` : ""));

// ---------------------------------------------------------------------
// 4. Tier contract: FAST must be the cheap default tier and non-empty;
//    SLOW must be opt-in. Heavier tiers must not be silently empty.
// ---------------------------------------------------------------------
assert(suitesIn("FAST").length > 0, `FAST tier is non-empty (${suitesIn("FAST").length} suites)`);
assert(suitesIn("SLOW").length > 0, `SLOW tier is non-empty (${suitesIn("SLOW").length} suites)`);
assert(FAST_BUDGET_S <= 60, `FAST budget contract is <= 60 s (declared ${FAST_BUDGET_S})`);

// ---------------------------------------------------------------------
// 5. EXEMPT entries must justify themselves — an empty reason is a bug.
// ---------------------------------------------------------------------
for (const x of EXEMPT) {
  assert(typeof x.why === "string" && x.why.trim().length > 10, `EXEMPT ${x.match} carries a written reason`);
  const used = listDir(path.join(HERE), /^test_.*\.js$/)
    .concat(listDir(path.join(ROOT, "scripts"), /^test_.*\.mjs$/))
    .filter((f) => x.match.test(f) && !registered.has(f) && !manualBasenames.has(path.basename(f)));
  assert(used.length > 0, `EXEMPT ${x.match} actually matches something (${used.join(", ") || "nothing — remove it"})`);
}

// ---------------------------------------------------------------------
// 6. Physics-critical suites stay in the gate (regression guard on the
//    registry itself: someone must not quietly demote these to manual).
// ---------------------------------------------------------------------
const MUST_GATE = [
  "tests/test_gb_fd.js",
  "tests/test_forces_fd.js",
  "tests/test_golden.js",
  "tests/test_nve.js",
  "tests/test_exclusions.js",
];
const fastSet = new Set(suitesIn("FAST").map((s) => s.file));
const demoted = MUST_GATE.filter((f) => !fastSet.has(f));
assert(demoted.length === 0, `physics-critical suites stay in FAST (${MUST_GATE.length})` + (demoted.length ? ` — DEMOTED: ${demoted.join(", ")}` : ""));

console.log(`\n=== test_suite_registry (${SUITES.length} suites registered, ${manual.size} manual): ${passed} PASSED, ${failed} FAILED ===`);
if (failed) process.exit(1);
