/**
 * bench/budget_check.js — the ONE place the performance contract is evaluated.
 *
 * WHY THIS FILE IS ITSELF PART OF THE FIX
 * ---------------------------------------
 * Before this existed, the same ~15 lines of python were copy-pasted into two
 * places: `.github/workflows/check.yml` and `manuscript/reproduce.sh`. They had
 * already drifted — reproduce.sh had a double-typo `if hv and hv and`, and BOTH
 * printed `fps budget {n}` for a quantity neither of them measured. One
 * implementation, called by both, is the only way a budget stays a contract
 * rather than a slogan.
 *
 * WHAT IT DOES
 * ------------
 *   1. runs bench/perf.js and parses its JSON summary,
 *   2. for every key in bench/budget.json with `_meta[key].measured === true`,
 *      compares the produced metric against the budget and reports OK / OVER,
 *   3. for every key with `measured === false`, prints NOT MEASURED plus the
 *      stated reason and compares NOTHING. An unmeasured key can never be
 *      silently reported as passing, and can never be silently dropped.
 *
 * EXIT CODE
 * ---------
 *   0  every measured key is within budget
 *   1  at least one measured key is over budget
 *
 * The caller decides whether exit 1 is a warning or an error; CI chooses warn
 * (this gate has never been a blocker and turning it into one is a policy
 * change, not a bug fix), reproduce.sh propagates it under `set -e`.
 *
 * Usage:  node bench/budget_check.js [--perf-log FILE]
 *           --perf-log reuses an existing bench/perf.js output instead of
 *           re-running the benchmark (used by reproduce.sh, which already ran it).
 */

import { execFileSync } from "child_process";
import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const budget = JSON.parse(fs.readFileSync(path.join(ROOT, "bench", "budget.json"), "utf-8"));

const argv = process.argv.slice(2);
const logIdx = argv.indexOf("--perf-log");
let perfText;
if (logIdx >= 0) {
  perfText = fs.readFileSync(argv[logIdx + 1], "utf-8");
} else {
  perfText = execFileSync(process.execPath, [path.join(ROOT, "bench", "perf.js")], {
    cwd: ROOT,
    encoding: "utf-8",
    stdio: ["ignore", "pipe", "inherit"],
  });
}

// The JSON summary is the last top-level {...} block perf.js prints. Take the
// first line that parses to an object carrying both metrics — matching on shape
// rather than on "the last brace" survives reformatting of perf.js.
let metrics = null;
for (const chunk of perfText.split("\n")) {
  const t = chunk.trim();
  if (!t.startsWith("{")) continue;
  try {
    const o = JSON.parse(t);
    if (o && typeof o.heavy?.meanMs === "number" && typeof o.cg?.meanMs === "number") {
      metrics = o;
      break;
    }
  } catch {}
}
if (!metrics) {
  // perf.js pretty-prints multi-line JSON. Re-scan with a brace-depth walk.
  const start = perfText.indexOf('{\n  "cg"');
  if (start >= 0) {
    let depth = 0;
    for (let i = start; i < perfText.length; i++) {
      if (perfText[i] === "{") depth++;
      else if (perfText[i] === "}") {
        depth--;
        if (depth === 0) {
          try {
            metrics = JSON.parse(perfText.slice(start, i + 1));
          } catch {}
          break;
        }
      }
    }
  }
}
if (!metrics) {
  console.error("budget_check: could not parse the bench/perf.js JSON summary. Refusing to report a budget result from a benchmark that did not run — an absent measurement is not a passing measurement.");
  process.exit(1);
}

/** Read `a.b.c` out of the perf summary. */
function dig(obj, dotted) {
  return dotted.split(".").reduce((o, k) => (o == null ? undefined : o[k]), obj);
}

const keys = Object.keys(budget).filter((k) => !k.startsWith("_"));
const meta = budget._meta || {};
let over = 0;
let unmeasured = 0;

console.log("performance contract — bench/budget.json");
for (const key of keys) {
  const m = meta[key];
  if (!m) {
    console.error(`  ERROR    ${key}: declared with no _meta entry. A budget key must say whether it is measured and by what.`);
    over++;
    continue;
  }
  if (m.measured !== true) {
    unmeasured++;
    const firstLine = String(m.reason || "(NO REASON GIVEN)").split("\n")[0];
    console.log(`  UNMEASURED ${key} — no number is claimed. ${firstLine}`);
    continue;
  }
  const got = dig(metrics, m.metric);
  if (typeof got !== "number" || !Number.isFinite(got)) {
    console.error(`  ERROR    ${key}: declared measured by ${m.producer} as ${m.metric}, but perf.js did not produce a finite number there. The producer and the budget disagree.`);
    over++;
    continue;
  }
  const limit = budget[key];
  if (typeof limit !== "number") {
    console.error(`  ERROR    ${key}: declared measured but its budget value is not a number (${JSON.stringify(limit)}).`);
    over++;
    continue;
  }
  const ratio = got / limit;
  if (ratio > 1) {
    over++;
    console.log(`  OVER      ${key}: ${got.toFixed(3)} ms > ${limit} ms budget (${ratio.toFixed(2)}x) — ::warning :: ${key} ${got.toFixed(3)} exceeds budget ${limit}`);
  } else {
    console.log(`  OK        ${key}: ${got.toFixed(3)} ms <= ${limit} ms budget (${ratio.toFixed(2)}x of budget; measured baseline ${m.measured_value})`);
  }
}

console.log(
  over === 0
    ? `budget: all ${keys.length - unmeasured} measured key(s) within budget; ${unmeasured} key(s) declared unmeasured and skipped.`
    : `budget: ${over} key(s) over budget.`
);

process.exit(over > 0 ? 1 : 0);
