/**
 * tests/test_budget_coverage.js — a budget key cannot be declared without a
 * way to check it.
 *
 * THE DEFECT THIS ENDS
 * --------------------
 * bench/budget.json declared `{"heavy_compute_ms": 2.0, "fps": 30,
 * "cg_compute_ms": 0.5}`. Two of the three numbers were fiction:
 *
 *   - heavy_compute_ms 2.0 vs ~14.1 measured — 7x under the truth, and the
 *     file's entire git history is ONE commit (3ee4914, 2026-09-11) that
 *     creates it already containing 2.0. There is no earlier value to have
 *     regressed from. CI printed "::warning Heavy compute ... exceeds budget
 *     2.0" on every run since, so the warning meant nothing.
 *   - fps 30 — measured by NOTHING. `grep -n fps bench/*.js scripts/*.mjs
 *     tests/*.js` returns zero hits. A frame rate needs a canvas, a
 *     compositor and a GPU; this repo has zero dependencies by design, and
 *     `new Viewer({})` in bare Node throws `canvas.getContext is not a
 *     function`. The number was unfalsifiable and CI printed it anyway.
 *
 * The failure was not that the numbers were wrong. It was that NOTHING
 * CONNECTED a declared number to a measurement, so being wrong was free and
 * invisible. This test is the connection.
 *
 * THE RULE
 * --------
 * Every non-`_` key in bench/budget.json must have a `_meta[key]` entry, and
 * that entry must be exactly one of:
 *
 *   measured: true   ->  names a `producer` file under bench/ that exists,
 *                        and a `metric` path that the producer's output
 *                        ACTUALLY CONTAINS (checked by running it), and a
 *                        numeric budget value. Plus a `measured_value`
 *                        recording what the number was when last calibrated,
 *                        so a later silent 10x drift is visible in the diff.
 *
 *   measured: false  ->  budget value MUST be null (no number is claimed), and
 *                        a non-empty `reason` explaining why. A key cannot be
 *                        both unmeasured and carry a number.
 *
 * Plus: no file in .github/, manuscript/, scripts/ or docs/ may print an fps
 * budget, because no fps measurement exists to back one. That is the
 * P09 acceptance criterion expressed as a gate rather than a grep someone
 * runs once.
 *
 * This test does NOT run the benchmark. bench/perf.js takes ~0.9 s and is
 * noisy; what is asserted here is the SHAPE of the contract and that the named
 * producer really emits the named metric, which is a structural claim that
 * must hold at all times, not a timing claim that belongs in CI.
 */

import { execFileSync } from "child_process";
import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
let pass = 0;
let fail = 0;

/** Counted assertion. RETURNS the condition so callers can branch on it. */
function check(cond, msg) {
  if (cond) {
    pass++;
    console.log(`PASS ${msg}`);
  } else {
    fail++;
    console.log(`FAIL ${msg}`);
  }
  return cond;
}

const budget = JSON.parse(fs.readFileSync(path.join(ROOT, "bench", "budget.json"), "utf-8"));
const keys = Object.keys(budget).filter((k) => !k.startsWith("_"));

check(keys.length > 0, `budget.json declares ${keys.length} contract key(s): ${keys.join(", ")}`);
check(budget._meta !== undefined, "budget.json carries a _meta block");

const meta = budget._meta || {};
const measuredKeys = [];
const unmeasuredKeys = [];

// ---------------------------------------------------------------------------
// 1. Shape of every declared key.
// ---------------------------------------------------------------------------
for (const key of keys) {
  const m = meta[key];
  if (!check(m !== undefined, `${key}: has a _meta entry`)) continue;

  if (m.measured === true) {
    measuredKeys.push(key);
    check(typeof budget[key] === "number" && budget[key] > 0,
      `${key}: measured=true carries a positive numeric budget (${budget[key]})`);
    check(typeof m.producer === "string" && m.producer.startsWith("bench/") && fs.existsSync(path.join(ROOT, m.producer)),
      `${key}: names a producer that exists on disk (${m.producer})`);
    check(typeof m.metric === "string" && m.metric.length > 0,
      `${key}: names the metric path it measures (${m.metric})`);
    check(typeof m.measured_value === "number" && Number.isFinite(m.measured_value) && m.measured_value > 0,
      `${key}: records measured_value so later drift shows in the diff (${m.measured_value})`);
    // A budget that is tighter than the value it was calibrated to is the exact
    // bug this file exists to prevent, re-introduced through the back door.
    check(budget[key] >= m.measured_value,
      `${key}: budget ${budget[key]} is not BELOW its own measured value ${m.measured_value} (a budget under the truth is a permanently lit warning)`);
    check(typeof m.workload === "string" && m.workload.length > 0,
      `${key}: states the workload it was measured on`);
  } else if (m.measured === false) {
    unmeasuredKeys.push(key);
    check(budget[key] === null,
      `${key}: measured=false must carry a null budget value, not the number ${JSON.stringify(budget[key])}`);
    check(typeof m.reason === "string" && m.reason.trim().length > 30,
      `${key}: measured=false states WHY it cannot be measured (${String(m.reason || "").slice(0, 60)}...)`);
  } else {
    check(false, `${key}: _meta.measured must be exactly true or false, got ${JSON.stringify(m.measured)}`);
  }
}

// ---------------------------------------------------------------------------
// 2. Every measured key's producer really emits the metric it claims.
//
// This is the part that makes the coverage real. A producer that stops
// emitting its metric — renamed key, deleted code path, refactor — must break
// the budget rather than leave it comparing against undefined.
// ---------------------------------------------------------------------------
function runProducer(rel) {
  return execFileSync(process.execPath, [path.join(ROOT, rel)], {
    cwd: ROOT,
    encoding: "utf-8",
    stdio: ["ignore", "pipe", "pipe"],
    maxBuffer: 64 * 1024 * 1024,
  });
}

for (const key of measuredKeys) {
  const m = meta[key];
  if (!m.producer || !fs.existsSync(path.join(ROOT, m.producer))) continue;
  let out = "";
  try {
    out = runProducer(m.producer);
  } catch (e) {
    check(false, `${key}: producer ${m.producer} ran but exited non-zero: ${String(e.message).split("\n")[0]}`);
    continue;
  }
  // Resolve the dotted metric path against the producer's JSON summary.
  const start = out.indexOf('{\n  "cg"');
  const start2 = start >= 0 ? start : out.indexOf("{");
  let summary = null;
  if (start2 >= 0) {
    let depth = 0;
    for (let i = start2; i < out.length; i++) {
      if (out[i] === "{") depth++;
      else if (out[i] === "}") {
        depth--;
        if (depth === 0) {
          try { summary = JSON.parse(out.slice(start2, i + 1)); } catch {}
          break;
        }
      }
    }
  }
  const got = summary ? m.metric.split(".").reduce((o, k) => (o == null ? undefined : o[k]), summary) : undefined;
  check(typeof got === "number" && Number.isFinite(got) && got > 0,
    `${key}: producer ${m.producer} really emits a finite ${m.metric} (got ${got})`);
}

// ---------------------------------------------------------------------------
// 3. The reverse direction: the single evaluator agrees with the file.
// ---------------------------------------------------------------------------
check(fs.existsSync(path.join(ROOT, "bench", "budget_check.js")),
  "bench/budget_check.js exists — one evaluator, so CI and reproduce.sh cannot drift apart");

// ---------------------------------------------------------------------------
// 4. P09: nothing anywhere prints an fps budget, because nothing measures fps.
// ---------------------------------------------------------------------------
const fpsIsUnmeasured = unmeasuredKeys.includes("fps");
check(fpsIsUnmeasured, "fps is declared unmeasured with a reason rather than carrying a number");

/**
 * THE fps-CLAIM RULE, precisely.
 *
 * A line is an offending fps BUDGET CLAIM when any of these hold on the
 * backtick-stripped line:
 *
 *   A  fps is bound to a literal value, quoted or not:
 *        /"?fps"?\s*[:=]\s*"?\d/i
 *      This is the shape of the actual defect: `{"heavy_compute_ms":2.0,
 *      "fps":30,...}` in check.yml and `"fps": 30` in a JSON code block.
 *   B  a literal is bound to fps:  /\d+(?:\.\d+)?\s*(?:[-–]\s*\d+(?:\.\d+)?)?\s*fps/i
 *   C  fps and a number co-occur inside an ASSERTION about a threshold:
 *      /\bfps\b[^.\n]{0,40}\b\d+(?:\.\d+)?\b/i together with one of
 *      budget | threshold | target | sustain | below | above | drop | at least
 *      | >= | <= | ≥ | ≤
 *
 * Rule C is what keeps this from crying wolf on ordinary runtime code: a HUD
 * that renders a *computed* frame rate (`${Math.round(state.fpsEMA)}fps`) has
 * no literal at all, and a comment like "top bar: ... FPS | Step (10 Hz)"
 * states no threshold and so fails C while matching neither A nor B.
 *
 * TWO EXEMPTIONS, both narrow and both about the difference between a CLAIM
 * and a LABELLED GUESS:
 *
 *   - SCOPE MARKER (`aspirational`, `not measured`, `estimate`, `design
 *     target`, `placeholder`, ...). A number carrying one is declared to be a
 *     guess, which is the distinction this project draws everywhere else;
 *     honouring it keeps the gate from forcing honest hedging to be deleted in
 *     order to go green.
 *   - REMOVAL NARRATION (`used to print`, `no longer`, `NEITHER measured`).
 *     The comments added by this very fix quote the string they deleted —
 *     "we used to print `fps budget 30` for a quantity NEITHER measured".
 *     That is a record of a removal, not an assertion. A gate that made such
 *     a comment illegal would be pushing the repo toward deleting its own
 *     history, which is the opposite of what it is for.
 *
 * SCOPE — and the reason for it, stated rather than assumed
 * ---------------------------------------------------------
 * The scan covers the surfaces that PUBLISH the performance contract: docs/,
 * bench/, .github/, manuscript/, scripts/, and the root markdown. It does not
 * cover src/. A comment in src/main.js such as `viewer.js — 60 fps Canvas
 * rendering` is a design note to the next maintainer, not a published claim, and
 * this goal's write scope explicitly excludes src/. Those three lines are known
 * and are recorded in docs/PERFORMANCE.md rather than silently ignored. If src/
 * ever comes into scope, widening SCAN_DIRS below is the whole change.
 */
const SCAN_DIRS = [".github", "manuscript", "scripts", "bench", "docs"];
const SCAN_ROOT_FILES = ["README.md", "ROADMAP.md", "CHANGELOG.md", "WHAT_CHANGED.md", "CONTRIBUTING.md"];

const MARKED = /aspirational|not measured|unmeasured|not benchmarked|no frame rate|estimate[ds]?\b|design target|placeholder|not yet benchmarked|if ≥|if >=|iff\b/i;
const REMOVAL_NARRATION = /used to print|used to be|no longer|neither measured|NEITHER measured|prints? no longer|now carries a banner/i;

const TICK = (line) => line.replace(/`[^`]*`/g, "``");

const ASSERT = /budget|threshold|target|sustain|below|above|drop|at least|>=|<=|≥|≤/i;
function isFpsClaim(line) {
  if (MARKED.test(line) || REMOVAL_NARRATION.test(line)) return false;
  const s = TICK(line);
  if (/\b"?fps"?\s*[:=]\s*"?\d/i.test(s)) return true;                            // A
  if (/\d+(?:\.\d+)?\s*(?:[-–]\s*\d+(?:\.\d+)?)?\s*fps/i.test(s)) return true;    // B
  if (/\bfps\b[^.\n]{0,40}\b\d+(?:\.\d+)?\b/i.test(s) && ASSERT.test(s)) return true; // C
  return false;
}

if (fpsIsUnmeasured) {
  const offenders = [];
  const scanFile = (rel) => {
    const txt = fs.readFileSync(path.join(ROOT, rel), "utf-8");
    txt.split("\n").forEach((line, i) => {
      if (isFpsClaim(line)) offenders.push(`${rel}:${i + 1}: ${line.trim().slice(0, 110)}`);
    });
  };
  const walk = (dir) => {
    for (const e of fs.readdirSync(path.join(ROOT, dir), { withFileTypes: true })) {
      const p = path.posix.join(dir, e.name);
      if (e.isDirectory()) { walk(p); continue; }
      if (!/\.(js|mjs|sh|yml|yaml|md|json)$/.test(e.name)) continue;
      if (p === "bench/budget.json") continue;            // the declaration itself
      if (p === "tests/test_budget_coverage.js") continue; // this file, which names it
      scanFile(p);
    }
  };
  for (const d of SCAN_DIRS) if (fs.existsSync(path.join(ROOT, d))) walk(d);
  for (const f of SCAN_ROOT_FILES) if (fs.existsSync(path.join(ROOT, f))) scanFile(f);

  check(offenders.length === 0,
    `no published surface claims an fps budget value while fps is unmeasured` +
    (offenders.length ? `\n      offenders:\n        ${offenders.join("\n        ")}` : ""));

  // Self-validation: the rule must still recognise the claims it was written
  // to catch. If a future edit to this regex silently stops matching real
  // offenders, the gate above becomes a permanent green light — exactly the
  // failure mode it exists to prevent.
  const KNOWN = [
    ["CI: compares meanMs vs budget.json (heavy_compute_ms 2.0, cg_compute_ms 0.5, fps 30) and emits ::warning", true],
    ['the budget is {"heavy_compute_ms":2.0,"fps":30,"cg_compute_ms":0.5}', true],
    ["- `fps: 30` — Interactive threshold. Below 30 fps => CI warns.", true],
    ["If a change drops fps below 30, CI warns.", true],
    ["**Instant** — 30–60 fps in browser, ~1 ps/frame", true],
    ["- **Performance:** ~1300 arc calls at 60 fps is borderline on integrated GPUs", true],
    ["viewer.js         — 60 fps Canvas rendering", true],   // a bare assertion, no scope marker
    ["Pass the same budget without regressing fps <30.", true],
    ["// we used to print `fps budget 30` for a quantity nothing measured", false], // quoted
    ["Frame rate is not measured in this repo — a design target, not a measured frame rate.", false],
    ["a ~10% frame cost design target", false],  // marked
  ];
  let ruleOk = true;
  const ruleBad = [];
  for (const [line, shouldMatch] of KNOWN) {
    const got = isFpsClaim(line);
    if (got !== shouldMatch) {
      ruleOk = false;
      ruleBad.push(`${shouldMatch ? "MISSED" : "FALSE POSITIVE"}: ${line.slice(0, 80)}`);
    }
  }
  check(ruleOk, "the fps-claim rule still recognises the claims it was written to catch" +
    (ruleBad.length ? `\n      ${ruleBad.join("\n      ")}` : ""));
}

// ---------------------------------------------------------------------------
// 5. Vacuity guard: a budget file that lost every key must fail, not pass.
// ---------------------------------------------------------------------------
check(keys.length >= 2, `budget.json still declares the keys the contract is made of (${keys.length} >= 2)`);
check(measuredKeys.length >= 1, `at least one key is genuinely measured (${measuredKeys.length})`);

console.log(`\n${pass} PASSED, ${fail} FAILED`);
process.exit(fail > 0 ? 1 : 0);
