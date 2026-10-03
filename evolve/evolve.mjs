#!/usr/bin/env node
/**
 * evolve.mjs — the evolution engine.
 *
 * One command = one evolution step over ONE goal. This is the whole point:
 * 200 evolutions must not cost 200x15 agent calls. Each step:
 *
 *   1. WORK   — pick the highest-value unblocked goal. The AUTHORED queue
 *               (evolve/queue/goals.json) is drained first; once it is
 *               exhausted `next` falls through to the SELF-GENERATED queue
 *               (evolve/queue/generated.json) that `plan` refills.
 *   2. THINK  — attach the wiki pattern ids + prior verdicts that bear on it
 *   3. ACT    — run the gate, record the delta in R (gate score), write trace
 *   4. LEARN  — append to the wiki skill-impact log; enqueue the NEXT 15 goals
 *
 * Usage:
 *   node evolve/evolve.mjs next            # print the single next goal (agent payload)
 *   node evolve/evolve.mjs gate            # run gate, print score, update R_best
 *   node evolve/evolve.mjs done <id> <verdict> [note]   # close a goal
 *   node evolve/evolve.mjs plan            # regenerate the next 15 goals from measurement
 *   node evolve/evolve.mjs status          # queue health
 *
 * `plan` flags:
 *   --no-gate   do not re-run the test suite; reuse state.gateVector
 *   --no-bench  do not re-run bench/perf.js; reuse state.benchVector
 *   --dry       print the plan, write nothing
 *
 * Contract: this file NEVER edits src/. It only reads, measures, and routes.
 */

import fs from "node:fs";
import path from "node:path";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";
// The gate's MEASUREMENTS live in their own module (evolve/gate-measure.mjs)
// because they are the part that has to be correct and the part most likely
// to be subtly wrong. Each returns its own coverage alongside its score: a
// number with no denominator is not a measurement.
import {
  checkSyntaxWithNode as checkSyntaxAll,
  domContract,
  coverageScore as measureCoverage,
  scienceSurface,
} from "./gate-measure.mjs";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..");
const Q = path.join(HERE, "queue", "goals.json");
const GEN = path.join(HERE, "queue", "generated.json");
const STATE = path.join(HERE, "queue", "state.json");
const TRACES = path.join(HERE, "traces");
const IMPACT = path.join(ROOT, ".wikiskill", "wiki", "evolution", "skill-impact.md");

/** How many goals `plan` emits. The L15 acceptance criterion, not a taste. */
const PLAN_SIZE = 15;

/** Fields every goal packet must carry. A goal missing one breaks the subagent. */
const GOAL_FIELDS = [
  "id", "horizon", "priority", "title", "thesis",
  "acceptance", "read_only", "touch", "do_not_touch", "patterns",
];
const HORIZONS = ["short", "mid", "long"];

/* ------------------------------------------------------------------ *
 * Context economy: the subagent never gets the whole repo. It gets
 * a goal packet with an explicit read-list. This is the single
 * highest-leverage anti-context-bloat mechanism in the loop.
 * ------------------------------------------------------------------ */

const readJson = (p, dflt) => {
  try { return JSON.parse(fs.readFileSync(p, "utf-8")); } catch { return dflt; }
};
const writeJson = (p, o) => {
  fs.mkdirSync(path.dirname(p), { recursive: true });
  fs.writeFileSync(p, JSON.stringify(o, null, 2) + "\n");
};

const state = readJson(STATE, { cycle: 0, R_best: 0, closed: [], rejected: [] });

/* ------------------------------------------------------------------ *
 * THE GATE — the single source of truth for "did this get better".
 *
 * R = weighted composite. Every evolution must move R up or be rolled
 * back. Weights are deliberately dominated by CORRECTNESS, because the
 * #1 documented failure of this project (see AUDIT_REPORT.md) was a
 * syntax error that shipped silently past every prose review.
 * ------------------------------------------------------------------ */

const WEIGHTS = {
  testPass: 0.25,   // tests/test_all.js — ratio vs the HIGH-WATER mark, not a fixed floor
  syntax: 0.10,     // node --check over every module
  dom: 0.05,        // ui.js ids present in index.html
  bloat: 0.20,      // inverse repo bytes + LOC (bloat is a real defect)
  science: 0.25,    // honest-scope surface: docs present, no new overclaims
  coverage: 0.15,   // anti-rot: is every test actually gated? (see S3)
};

// Saturation guard: a component pinned at 1.0 carries no gradient, so the loop
// cannot distinguish "improved" from "unchanged". Weights are renormalised over
// NON-saturated components only. A goal that only moves a dead component scores
// zero and gets reported as such, rather than silently inflating R.
const SAT = 0.999;


function runGate(opts = {}) {
  const r = { score: 0, parts: {}, notes: [] };

  // 1. syntax — widened from src/*.js to EVERY .js/.mjs in the declared roots.
  //    The old walk examined 110 files and left the other ~100 modules
  //    (including evolve/evolve.mjs itself) to a different script.
  const syn = checkSyntaxAll(ROOT);
  r.parts.syntax = syn.total ? (syn.total - syn.failed.length) / syn.total : 0;
  r.syntaxNumbers = syn;
  if (syn.failed.length) r.notes.push(`syntax FAIL: ${syn.failed.join(", ")}`);
  if (syn.missingRoots && syn.missingRoots.length) {
    r.notes.push(`syntax roots missing (coverage silently reduced): ${syn.missingRoots.join(", ")}`);
  }

  // 2. tests — ratio against the recorded high-water mark, so adding real
  //    coverage keeps moving the needle instead of pinning at a legacy floor.
  //
  //    Two guards, both learned the hard way:
  //    --fast  reuses the last recorded result. The gate test probes the gate
  //            7x; re-running the suite each time made that test O(n^2).
  //    REENTRANCY  tests/test_evolve_gate.js runs INSIDE the suite, so when it
  //            calls the gate, the gate's own test run would re-enter the
  //            suite that invoked it — unbounded recursion, surfaced as
  //            "tests crashed" and a silent testPass=0. A gate that crashes
  //            when run from the suite it gates is a landmine. Env-inherited,
  //            because argv cannot see the grandparent process.
  if (opts.fast || process.env.EVOLVE_GATE_CHILD === "1") {
    r.parts.testPass = state.testHighWater && r.parts.syntax >= 1 ? 1 : 0;
    r.testNumbers = {
      passed: state.testHighWater || 0, failed: 0, cached: true,
      note: "reused recorded result (fast/child mode); run `gate` standalone for a live measurement",
    };
  } else try {
    const out = execFileSync(process.execPath, [path.join(ROOT, "tests", "test_all.js")], {
      encoding: "utf-8", stdio: ["ignore", "pipe", "pipe"], timeout: 300000,
      env: { ...process.env, EVOLVE_GATE_CHILD: "1" },
    });
    const m = [...out.matchAll(/(\d+) PASSED, (\d+) FAILED/g)];
    const [p, f] = m.length ? [Number(m.at(-1)[1]), Number(m.at(-1)[2])] : [0, 1];
    const hw = Math.max(state.testHighWater || 0, 352);
    r.parts.testPass = f > 0 ? 0 : Math.min(1, p / hw);
    r.testNumbers = { passed: p, failed: f, highWater: hw };
    if (p > hw && f === 0) state.testHighWater = p;
    if (f > 0) r.notes.push(`tests FAIL: ${f}`);
  } catch (e) {
    r.parts.testPass = 0;
    // A crashed or timed-out suite is a HARD failure, not a zero. Treating it
    // as merely "0 points" let a broken tree be accepted as the new reference
    // state: R_best was re-baselined from a run whose suite had crashed, and
    // the gate then reported OPEN on that broken state. A measurement that
    // did not happen must never become the baseline.
    r.testsCrashed = true;
    r.notes.push(`tests CRASHED (hard failure, cannot be a reference state): ${String(e.message).slice(0, 160)}`);
  }

  // 2b. coverage / anti-rot: every test suite in the declared POPULATION must
  //     be wired into a tier that actually runs, or be explicitly retired.
  //     This was widened (see evolve/gate-measure.mjs): the old version
  //     counted only tests/test_*.js and counted a suite "wired" even when
  //     only --slow ever ran it.
  const cov = measureCoverage(ROOT);
  r.parts.coverage = cov.score;
  r.coverageNumbers = cov.numbers;

  // 3. DOM contract — widened from src/ui.js alone to EVERY src/ module, so an
  //    id a controller binds via getElementById is checked too. Previously
  //    112 of 136 ids were examined and 24 were invisible.
  try {
    const dc = domContract(ROOT);
    r.parts.dom = dc.score;
    r.domNumbers = dc.numbers;
  } catch (e) { r.parts.dom = 0; r.notes.push(`dom contract failed: ${String(e.message).slice(0, 100)}`); }

  // 4. bloat — measured, not vibes. git objects + tracked bytes + src LOC.
  const tracked = readTrackedStats();
  r.parts.bloat = bloatScore(tracked);
  // Widened: this used to read README.md.slice(0, 4000) — 15.6% of a 25.6 KB
  // document — and score a keyword hit as proof of honest scope. It now
  // reads the real trust-boundary sources (ROADMAP.md §1, LIMITATIONS,
  // APPLICABILITY, VALIDATION) and reports its own coverage.
  r.scienceNumbers = scienceSurface(ROOT);
  r.parts.science = r.scienceNumbers.score;
  r.tracked = tracked;

  // Renormalise over non-saturated components only.
  const live = Object.entries(WEIGHTS).filter(([k]) => (r.parts[k] || 0) < SAT);
  const dead = Object.entries(WEIGHTS).filter(([k]) => (r.parts[k] || 0) >= SAT).map(([k]) => k);
  const wsum = live.reduce((a, [, w]) => a + w, 0) || 1;
  r.score = live.reduce((a, [k, w]) => a + (w / wsum) * (r.parts[k] || 0), 0);
  r.saturated = dead;
  if (dead.length) r.notes.push(`saturated (no gradient, weight redistributed): ${dead.join(", ")}`);
  return r;
}

function readTrackedStats() {
  const s = { files: 0, bytes: 0, blobBytes: 0, srcLoc: 0, bigFiles: [], locs: [] };
  let list = [];
  try { list = execFileSync("git", ["ls-files"], { cwd: ROOT, encoding: "utf-8" }).split("\n").filter(Boolean); }
  catch { return s; }
  s.files = list.length;
  // Tracked bytes, split by kind. `blobBytes` is the only one the gate
  // penalises: shipped data and generated artifacts are real bloat, while
  // source, tests and docs are the work itself.
  const CODE = /^(src\/|tests\/|scripts\/|evolve\/|docs\/|\.wikiskill\/|bench\/)/;
  for (const f of list) {
    try {
      const st = fs.statSync(path.join(ROOT, f));
      s.bytes += st.size;
      if (!CODE.test(f)) s.blobBytes += st.size;
      if (st.size > 1_000_000) s.bigFiles.push([f, Math.round(st.size / 1024) + "KB"]);
    } catch {}
  }
  s.bigFiles.sort((a, b) => parseInt(b[1]) - parseInt(a[1]));

  // src/ metrics come from the FILESYSTEM, not from git. Using `git
  // ls-files` here made the gate blind to uncommitted work: a probe that
  // added a 3000-line junk module to src/ did not move the score at all,
  // because the file was not yet tracked. A metric that cannot see the
  // edits an agent just made cannot gate that agent.
  (function walk(d) {
    let entries = [];
    try { entries = fs.readdirSync(d, { withFileTypes: true }); } catch { return; }
    for (const e of entries) {
      if (e.name.startsWith(".")) continue;
      const p = path.join(d, e.name);
      if (e.isDirectory()) { walk(p); continue; }
      if (!e.name.endsWith(".js")) continue;
      try {
        const loc = fs.readFileSync(p, "utf-8").split("\n").length;
        s.srcLoc += loc;
        s.locs.push([path.relative(ROOT, p), loc]);
      } catch {}
    }
  })(path.join(ROOT, "src"));

  // The "typical module size" baseline must be robust to small additions.
  // Including every module means adding one 20-line leaf file shifts the
  // median down, which raises max/median and penalises good work. Modules
  // below MIN_TYPICAL_LOC are excluded from the baseline distribution: they
  // are leaves, not evidence about the shape of the codebase.
  const MIN_TYPICAL_LOC = 50;
  const typical = s.locs.map((x) => x[1]).filter((n) => n >= MIN_TYPICAL_LOC).sort((a, b) => a - b);
  s.medianLoc = typical.length
    ? (typical.length % 2
        ? typical[(typical.length - 1) / 2]
        : (typical[typical.length / 2 - 1] + typical[typical.length / 2]) / 2)
    : 0;
  s.maxLoc = s.locs.length ? Math.max(...s.locs.map((x) => x[1])) : 0;
  s.maxLocFile = s.locs.find((x) => x[1] === s.maxLoc)?.[0] || "";
  return s;
}

/** Count test files and how many are wired into a suite vs retired to manual. */
function countTests() {
  const dir = path.join(ROOT, "tests");
  let all = [];
  try { all = fs.readdirSync(dir).filter((f) => /^test_.*\.js$/.test(f)); } catch { return { total: 0, wired: 0, manual: 0 }; }
  let wired = 0;
  try {
    const reg = fs.readFileSync(path.join(dir, "suites.js"), "utf-8");
    for (const f of all) if (reg.includes(f)) wired++;
  } catch {}
  let manual = 0;
  try { manual = fs.readdirSync(path.join(dir, "manual")).filter((f) => /^test_.*\.js$/.test(f)).length; } catch {}
  return { total: all.length, wired, manual };
}

function coverageScore() {
  const c = countTests();
  if (!c.total) return 1;
  return Math.min(1, (c.wired + c.manual) / c.total);
}

/**
 * Bloat score in [0,1]. 1 = lean.
 *
 * Bloat is REDUNDANCY, not volume. The previous version penalised total
 * tracked bytes and total src LOC, which meant a goal that legitimately
 * added a test file or a shared physics kernel LOWERED the score. That is
 * a broken incentive: the loop would refuse good work. It also made R_best
 * unreachable, because R_best was recorded on a tree that no longer existed.
 *
 * What is genuinely bloat:
 *   - oversized individual tracked files (a 4MB data file in a demo repo)
 *   - a god module far above the median (structural, not volume)
 *   - duplicated CODE, measured by a normalised-duplication scan
 * Volume of tests, docs, and shared kernels is not penalised.
 */
function bloatScore(s) {
  const filePenalty = 1 / (1 + s.bigFiles.length / 20);
  // Size penalty applies to DATA/BLOB bloat only — shipped datasets and
  // generated artifacts. It deliberately excludes src/, tests/, scripts/,
  // evolve/ and docs/: writing code or adding a test is not bloat, and a
  // metric that charges for writing tests will refuse good work. (It did:
  // ~57 KB of new test code cost exactly 1e-4 of R.)
  const blobScore = 1 / (1 + s.blobBytes / (500 * 1024 * 1024));
  // Structural: ratio of the largest module to the median. A file at the
  // median scores 1; a 6x outlier is penalised.
  const ratio = s.medianLoc > 0 ? s.maxLoc / s.medianLoc : 1;
  const structureScore = 1 / (1 + Math.max(0, ratio - 2) / 3);
  return 0.45 * filePenalty + 0.25 * blobScore + 0.30 * structureScore;
}

/** Honest-scope surface: the project's defining virtue (ROADMAP.md §1). */
function scienceScore() {
  const must = [
    "ROADMAP.md", "docs/LIMITATIONS.md", "docs/APPLICABILITY.md",
    "docs/VALIDATION.md", "LICENSE", "CITATION.cff", "CHANGELOG.md",
  ];
  const have = must.filter((f) => fs.existsSync(path.join(ROOT, f))).length;
  const overclaim = /NOT a replacement|qualitative|not FEP|trust boundary/i.test(
    fs.readFileSync(path.join(ROOT, "README.md"), "utf-8").slice(0, 4000));
  return (have / must.length) * (overclaim ? 1 : 0.5);
}

/* ------------------------------------------------------------------ *
 * Goal queue — AUTHORED first, then SELF-GENERATED.
 *
 * Two files, never one. goals.json is the human-authored backlog and is
 * append-only with respect to verdicts; generated.json is the planner's
 * output. `plan` only ever writes generated.json, so re-planning cannot
 * destroy a closed goal or its verdict. `next` drains goals.json first and
 * falls through to generated.json, so the authored queue is never starved
 * and the generated queue is never dropped.
 * ------------------------------------------------------------------ */

const queue = readJson(Q, { goals: [] });
const generated = readJson(GEN, { goals: [] });

const HORIZON_ORDER = { short: 0, mid: 1, long: 2 };

/** Goals in `next` order: horizon, then priority, then a STABLE tiebreak.
 *  The old comparator used `a.cycle - b.cycle`, which is NaN for goals that
 *  carry no `cycle` field — Array.sort with a NaN comparator is unspecified,
 *  so the queue order was not a function of the queue contents. */
function goalCmp(a, b) {
  const h = (HORIZON_ORDER[a.horizon] ?? 9) - (HORIZON_ORDER[b.horizon] ?? 9);
  if (h) return h;
  const p = (a.priority ?? 9) - (b.priority ?? 9);
  if (p) return p;
  return String(a.id).localeCompare(String(b.id));
}

function openGoals(list) {
  return list.filter((g) => !g.status || g.status === "open");
}

/**
 * All open goals, AUTHORED FIRST AS A HARD PARTITION.
 *
 * The authored backlog is a human commitment and is drained completely
 * before a single generated goal is served; only once it is empty does the
 * self-generated queue take over. Interleaving the two by rank would let a
 * generated goal jump ahead of a commitment a human already made, which
 * would make the authored queue a suggestion rather than a contract.
 *
 * Each partition is sorted independently, so the authored ordering is
 * bit-for-bit what it was before this file knew the generated queue existed.
 */
function allOpenGoals() {
  return [
    ...openGoals(queue.goals).sort(goalCmp),
    ...openGoals(generated.goals).sort(goalCmp),
  ];
}

function nextGoal() {
  // tolerate goals authored without an explicit status (default = open)
  for (const g of queue.goals) if (!g.status) g.status = "open";
  for (const g of generated.goals) if (!g.status) g.status = "open";
  return allOpenGoals()[0] || null;
}

/** Every goal the loop knows about, for id lookup by `done`. */
function findGoal(id) {
  return queue.goals.find((g) => g.id === id) || generated.goals.find((g) => g.id === id) || null;
}

function cmdNext() {
  const g = nextGoal();
  if (!g) { console.log(JSON.stringify({ done: true, cycle: state.cycle }, null, 2)); return; }
  const patterns = g.patterns || [];
  const rejected = state.rejected.slice(-5);
  const packet = {
    cycle: state.cycle,
    goal: { id: g.id, title: g.title, horizon: g.horizon, priority: g.priority },
    thesis: g.thesis,
    acceptance: g.acceptance,
    read_only: g.read_only ?? [],
    touch: g.touch ?? [],
    do_not_touch: g.do_not_touch ?? ["data/coreset/**", "manuscript/*.pdf", "node_modules/**"],
    wiki_patterns: patterns.map((p) => {
      const f = path.join(ROOT, ".wikiskill", "wiki", "patterns", `${p}.md`);
      const body = fs.existsSync(f) ? fs.readFileSync(f, "utf-8").split("\n").slice(0, 22).join("\n") : "(missing)";
      return { id: p, excerpt: body };
    }),
    recently_rejected: rejected,
    finish: `node evolve/evolve.mjs done ${g.id} <ACCEPTED|ROLLED_BACK> "<one-line reason>"`
  };
  console.log(JSON.stringify(packet, null, 2));
}

function cmdGate() {
  // EVOLVE_GATE_RBEST overrides the reference score for this run only. It
  // exists so the exit-code path can be tested deterministically: a probe
  // cannot (and must not) write R_best, so without a seam the only way to
  // observe a CLOSED verdict would be to corrupt the real baseline.
  if (process.env.EVOLVE_GATE_RBEST) {
    state.R_best = Number(process.env.EVOLVE_GATE_RBEST);
  }
  const r = runGate({ fast: process.argv.includes("--fast") });
  console.log(JSON.stringify({
    score: Number(r.score.toFixed(4)),
    R_best: state.R_best,
    verdict: r.score >= state.R_best ? "OPEN (accept)" : "CLOSED (roll back)",
    parts: Object.fromEntries(Object.entries(r.parts).map(([k, v]) => [k, Number(v.toFixed(4))])),
    weights: WEIGHTS,
    saturated: r.saturated,
    testNumbers: r.testNumbers,
    syntaxNumbers: r.syntaxNumbers,
    coverageNumbers: r.coverageNumbers,
    domNumbers: r.domNumbers,
    scienceNumbers: r.scienceNumbers,
    tracked: r.tracked && {
      files: r.tracked.files, MB: +(r.tracked.bytes / 1048576).toFixed(1),
      srcLoc: r.tracked.srcLoc, medianLoc: r.tracked.medianLoc,
      maxLoc: r.tracked.maxLoc, maxLocFile: r.tracked.maxLocFile,
      bigFiles: r.tracked.bigFiles.slice(0, 8),
    },
    notes: r.notes,
  }, null, 2));
  // Compare at 4dp so float noise does not read as a regression.
  const r4 = Number(r.score.toFixed(4));
  // A crashed suite is never OPEN, whatever the composite says: the
  // measurement did not happen, so it cannot certify anything. Without this
  // the gate reported OPEN on a tree whose suite had crashed.
  const open = r4 >= state.R_best && !r.testsCrashed;
  // The high-water must advance whenever the suite grew, independently of
  // whether the composite moved. The old predicate compared the high-water
  // against a hardcoded 352 and could therefore never fire again after the
  // first bump, so the mark stayed frozen and 30 real assertions were
  // unearned headroom that a later test deletion could spend for free.
  const hwGrew = !!(r.testNumbers && r.testNumbers.failed === 0 &&
                    r.testNumbers.passed > (state.testHighWater || 0));
  // Record the vector so `plan` can diff against it. The PREVIOUS vector
  // moves to prevGateVector: a regression is the loop's highest-value
  // signal and it needs a baseline to be visible at all. Probes are excluded
  // (see isProbe below) so a test probe cannot fabricate a regression.
  if (!(process.argv.includes("--fast") || process.env.EVOLVE_GATE_CHILD === "1")) {
  state.prevGateVector = state.gateVector || null;
  state.gateVector = {
    cycle: state.cycle,
    score: r4,
    parts: Object.fromEntries(Object.entries(r.parts).map(([k, v]) => [k, Number(v.toFixed(6))])),
    saturated: r.saturated,
  };
  } // end !isProbe (vector recording)
  if (hwGrew) state.testHighWater = r.testNumbers.passed;

  // A PROBE must never mutate the baseline it is measuring. The gate test
  // plants junk files and calls the gate; without this guard those probes
  // wrote their vectors into state, and the planner then diffed against a
  // fabricated "current" tree and invented a phantom regression goal
  // (P12: "bloat -0.0493", caused entirely by the test's own probe file).
  // An instrument that rewrites its own reference is not measuring.
  const isProbe = process.argv.includes("--fast") || process.env.EVOLVE_GATE_CHILD === "1";
  if (!isProbe) {
  if (r4 > state.R_best) {
    state.R_best = r4;
    writeJson(STATE, state);
    console.error(`\n[R_best updated -> ${state.R_best}]`);
  } else {
    // Persist the recorded vector + high-water even when the composite did
    // not move; `plan` is useless without a baseline to diff against.
    writeJson(STATE, state);
  }
  } // end !isProbe
  if (state.rBaselineNote) console.error(`[NOTE] ${state.rBaselineNote}`);

  // A gate that prints "roll back" and exits 0 cannot fail anything: in CI
  // it was decoration, not a gate. Exit non-zero when closed, unless the
  // caller asked for machine-readable output only (--json).
  if (process.argv.includes("--json")) return;
  process.exit(open ? 0 : 1);
}

function cmdDone(id, verdict, note = "") {
  const g = findGoal(id);
  if (!g) { console.error(`no such goal: ${id}`); process.exit(1); }
  g.status = verdict === "ACCEPTED" ? "closed" : "rolled-back";
  g.verdict = verdict;
  g.note = note;
  g.closedAtCycle = state.cycle;
  // Write back to the file the goal actually came from, so closing a
  // generated goal never rewrites the authored backlog.
  writeJson(generated.goals.includes(g) ? GEN : Q, generated.goals.includes(g) ? generated : queue);
  if (verdict !== "ACCEPTED") state.rejected.push(`${id}: ${note}`.slice(0, 200));
  state.cycle++;
  // R_best is the reference state, and the reference state is the last
  // ACCEPTED tree. If it only ever moved upward on its own, then any
  // accepted goal that legitimately changes the score would leave the gate
  // permanently CLOSED — the loop would be unable to accept anything after
  // the first real change. Re-baselining on acceptance keeps the gate's
  // question meaningful: "has anything got worse since we last accepted?"
  if (verdict === "ACCEPTED") {
    const g2 = runGate({ fast: process.argv.includes("--fast") });
    if (g2.testsCrashed) {
      // Never re-baseline onto a state whose suite could not be measured.
      // This happened for real: a crashed suite scored 0.355, that became
      // R_best, and every subsequent run then read OPEN on the broken tree.
      console.error("[R_best NOT re-baselined: the test suite crashed — " +
        "a state that could not be measured must not become the reference]");
    } else {
    state.R_best = Number(g2.score.toFixed(4));
    state.gateVector = {
      cycle: state.cycle, score: state.R_best,
      parts: Object.fromEntries(Object.entries(g2.parts).map(([k, v]) => [k, Number(v.toFixed(6))])),
      saturated: g2.saturated,
    };
    console.error(`[R_best re-baselined to ${state.R_best} on ACCEPT (cycle ${state.cycle})]`);
    }
  }
  writeJson(STATE, state);

  const trace = path.join(TRACES, `trace-${String(state.cycle).padStart(3, "0")}-${id}.md`);
  fs.mkdirSync(TRACES, { recursive: true });
  fs.writeFileSync(trace, [
    `# Evolution ${state.cycle} — ${id} (${g.horizon})`, "",
    `**Goal:** ${g.title}`, `**Verdict:** ${verdict}`, `**Note:** ${note || "(none)"}`, "",
    `**Acceptance:** ${g.acceptance}`, "",
    `**Touched:** ${(g.touch || []).join(", ") || "(none)"}`, "",
  ].join("\n"));

  appendImpact(id, g, verdict, note);
  console.log(`closed ${id} [${verdict}] -> cycle ${state.cycle}`);
}

function appendImpact(id, g, verdict, note) {
  const entry = [
    ``, `## Evolution ${state.cycle} — ${id} — ${g.title}`, ``,
    `- **Horizon:** ${g.horizon} (priority P${g.priority})`,
    `- **Motivated by:** ${(g.patterns || []).join(", ") || "(direct audit finding)"}`,
    `- **Verdict:** ${verdict}`,
    `- **Note:** ${note || "(none)"}`,
  ].join("\n");
  fs.appendFileSync(IMPACT, entry);
}

function cmdStatus() {
  const tally = (list) => {
    const byH = {};
    for (const g of list) {
      byH[g.horizon] ??= { open: 0, closed: 0, "rolled-back": 0 };
      byH[g.horizon][g.status === "closed" ? "closed" : g.status === "rolled-back" ? "rolled-back" : "open"]++;
    }
    return byH;
  };
  const next = nextGoal();
  console.log(JSON.stringify({
    cycle: state.cycle, R_best: state.R_best,
    goals: { total: queue.goals.length, byHorizon: tally(queue.goals) },
    generated: {
      total: generated.goals.length,
      byHorizon: tally(generated.goals),
      generatedAtCycle: generated._meta?.cycle ?? null,
      deficitPoolExhausted: generated._meta?.deficitPoolExhausted ?? null,
    },
    open: { authored: openGoals(queue.goals).length, generated: openGoals(generated.goals).length },
    next: next && `${next.id} — ${next.title}`,
    nextSource: next ? (queue.goals.includes(next) ? "authored" : "generated") : null,
    rejectedRemembered: state.rejected.length,
  }, null, 2));
}

/* ================================================================== *
 * THE PLANNER — goal L15.
 *
 * Refills the queue from MEASUREMENT, not from a human's attention.
 *
 * Design contract (this is the part that matters, so it is stated):
 *
 *   1. Every candidate goal is produced by a PROBE. A probe reads the repo
 *      or the last gate vector and returns a goal ONLY IF it finds a real,
 *      quantified defect. No probe invents a defect. If a probe finds
 *      nothing, it contributes nothing and says so in the census.
 *   2. Every candidate carries `deficit` — a NUMBER computed from the
 *      measurement, not a hand-assigned rank. Ranking is `deficit`
 *      descending. A goal that cannot state its number does not ship.
 *   3. Every candidate carries `thesis` quoting the measurement that
 *      produced it, so a reader can re-derive the ranking by hand.
 *   4. `plan` is IDEMPOTENT-ish: a goal that is already in either queue
 *      (open OR closed) is not re-emitted, so running `plan` twice in a row
 *      is stable and closed goals are never resurrected or overwritten.
 *   5. If fewer than PLAN_SIZE real candidates exist, `plan` emits fewer
 *      and RECORDS WHY (deficitPoolExhausted). Padding is forbidden: a
 *      filler goal ("improve code quality") costs a whole agent call and
 *      returns nothing, which is worse than an empty queue.
 * ================================================================== */

/** Recursively list files under dir matching a filename test. */
function listFiles(dir, re, out = []) {
  let entries = [];
  try { entries = fs.readdirSync(dir, { withFileTypes: true }); } catch { return out; }
  for (const e of entries) {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) listFiles(p, re, out);
    else if (re.test(e.name)) out.push(path.relative(ROOT, p).split(path.sep).join("/"));
  }
  return out;
}

/** git-tracked file list, or [] when git is unavailable. */
function trackedFiles() {
  try {
    return execFileSync("git", ["ls-files"], { cwd: ROOT, encoding: "utf-8" })
      .split("\n").filter(Boolean);
  } catch { return []; }
}

const pct = (n) => `${(n * 100).toFixed(1)}%`;

/** Things a probe looked at and REFUSED to turn into a goal, with the reason.
 *  This is the ledger that keeps the planner honest: a probe that finds
 *  nothing must say so, and a probe that could not decide must say that
 *  too, rather than emitting a goal built on a guess. */
const INCONCLUSIVE = [];

/* ------------------------------------------------------------------ *
 * PROBE 1 — the gate vector itself.
 *
 * The gate is the loop's only source of truth about "is this better". Two
 * ways it can be lying: a component that is SATURATED (pinned, so it cannot
 * distinguish improved from unchanged — the loop has run out of signal), and
 * a component that REGRESSED against the last recorded vector (the highest
 * priority signal the loop has, and the one it currently throws away).
 * ------------------------------------------------------------------ */
function probeGateVector() {
  const found = [];

  // -- 1a. Saturation. 5/6 dead means the composite is a function of ONE
  //       number, so most goals move R by exactly zero and are indistinguishable.
  const lastVec = state.gateVector;
  const parts = lastVec?.parts;
  if (parts) {
    const dead = Object.keys(WEIGHTS).filter((k) => (parts[k] ?? 0) >= SAT);
    const live = Object.keys(WEIGHTS).filter((k) => (parts[k] ?? 0) < SAT);
    if (dead.length >= Object.keys(WEIGHTS).length - 1) {
      found.push({
        key: "gate-saturation",
        deficit: dead.length / Object.keys(WEIGHTS).length,
        priority: 0, horizon: "mid",
        title: `Give the gate a live gradient: ${dead.length}/${Object.keys(WEIGHTS).length} components are saturated`,
        thesis: `gate vector = ${JSON.stringify(parts)}; saturated = [${dead.join(", ")}]; only [${live.join(", ")}] still carries gradient. Weights are renormalised over non-saturated components (evolve.mjs:${"SAT"}=${SAT}), so R is currently a function of ${live.length} number(s) alone. A goal that moves any saturated component scores exactly 0 change, so the loop cannot tell "improved" from "unchanged" on ${dead.length} of ${Object.keys(WEIGHTS).length} axes. Measured cost: the R_best-vs-current gap is ${(state.R_best - (lastVec.score ?? 0)).toFixed(4)} and 0 of ${dead.length} saturated components can close it.`,
        acceptance: `A new discriminating component exists with measured spread across the repo (report min/median/max over >= 5 real samples) and it appears unsaturated in \`node evolve/evolve.mjs gate\` output; OR the saturated components are replaced by metrics that provably move. \`node evolve/evolve.mjs gate\` shows saturated.length <= ${Math.max(1, dead.length - 1)}.`,
        read_only: ["evolve/evolve.mjs", "evolve/queue/state.json", "evolve/traces/"],
        touch: ["evolve/evolve.mjs", "tests/", "evolve/reports/"],
        patterns: ["P3-unsurfaced-value"],
      });
    }

    // -- 1b. Regression against the last recorded vector. The strongest
    //       signal the loop owns, and it is discarded today.
    //
    //       Compared at 4dp for the same reason cmdGate does: bloat is a
    //       float expression over git-tracked byte counts, so two runs on an
    //       identical tree differ in the 6th decimal. A "regression" of
    //       6e-6 is noise, and a goal built on it is a fabricated goal.
    const EPS = 1e-4;
    const prev = state.prevGateVector;
    if (prev?.parts) {
      const regressed = Object.keys(WEIGHTS)
        .filter((k) => prev.parts[k] !== undefined && parts[k] !== undefined && parts[k] < prev.parts[k] - EPS)
        .map((k) => ({ k, from: prev.parts[k], to: parts[k], d: prev.parts[k] - parts[k] }));
      if (!regressed.length) {
        const deltas = Object.keys(WEIGHTS)
          .filter((k) => prev.parts[k] !== undefined && parts[k] !== undefined)
          .map((k) => Math.abs(parts[k] - prev.parts[k]));
        const maxDelta = deltas.length ? Math.max(...deltas) : 0;
        if (maxDelta > 0) {
          INCONCLUSIVE.push({
            probe: "gate-vector/regression", subject: `prev(cycle ${prev.cycle}) vs current(cycle ${lastVec.cycle})`,
            why: `every component moved by <= ${EPS} (max ${maxDelta.toExponential(2)}) — that is float noise on an unchanged tree, not a regression, so NO goal was generated`,
          });
        }
      }
      if (regressed.length) {
        const worst = regressed.sort((a, b) => b.d - a.d)[0];
        found.push({
          key: "gate-regression",
          deficit: Math.min(1, worst.d * 10),
          priority: 0, horizon: "short",
          title: `Gate component REGRESSED: ${worst.k} fell ${worst.from.toFixed(4)} -> ${worst.to.toFixed(4)}`,
          thesis: `prev vector = ${JSON.stringify(prev.parts)} (cycle ${prev.cycle}); current = ${JSON.stringify(parts)} (cycle ${lastVec.cycle}). ${regressed.length} component(s) moved down: ${regressed.map((r) => `${r.k} -${r.d.toFixed(4)}`).join(", ")}. A regression is the single highest-value signal the loop can act on and it must be fixed before any new goal is accepted.`,
          acceptance: `The regressing component returns to >= ${worst.from.toFixed(4)} (or the regression is explained in writing in evolve/reports/ and R_best is re-baselined with an updated rBaselineNote). \`node evolve/evolve.mjs gate\` reports no component below its previous vector.`,
          read_only: ["evolve/evolve.mjs", "evolve/queue/state.json", "evolve/traces/"],
          touch: ["src/", "tests/", "evolve/reports/"],
          patterns: ["P3-unsurfaced-value"],
        });
      }
    }

    // -- 1c. The live component: bloat is currently the ONLY gradient, and
    //       its dominant term is src LOC. That is measurable and specific.
    if (parts.bloat !== undefined && parts.bloat < SAT) {
      const t = readTrackedStats();
      const headroom = 1 - parts.bloat;
      const locTerm = 0.25 * (1 - 1 / (1 + Math.max(0, t.srcLoc - 20000) / 20000));
      if (t.srcLoc > 20000) {
        found.push({
          key: "bloat-srcloc",
          deficit: headroom,
          priority: 1, horizon: "mid",
          title: `Cut src LOC below the bloat knee: ${t.srcLoc} LOC vs the 20000 penalty start`,
          thesis: `bloat=${parts.bloat.toFixed(4)} is the only non-saturated gate component (everything else is pinned at 1.0), and ${pct(locTerm / headroom)} of its remaining headroom is the src-LOC term 1/(1+(L-20000)/20000) evaluated at L=${t.srcLoc}. dR/dL = -0.093 per 100 LOC at this point, so the loop currently pays a measurable price for every line of science code it writes — an incentive pointing the wrong way. Cutting L to 20000 recovers ${locTerm.toFixed(4)} of the ${headroom.toFixed(4)} deficit.`,
          acceptance: `src LOC (sum of lines over src/**/*.js as counted by evolve.mjs readTrackedStats) <= 20000, or the bloat LOC term is replaced by a term with an argued-for justification. \`node evolve/evolve.mjs gate\` reports bloat >= ${(1 - headroom + locTerm).toFixed(4)}.`,
          read_only: ["evolve/evolve.mjs", "docs/PERFORMANCE.md", "bench/budget.json", "CHANGELOG.md"],
          touch: ["src/", "tests/", "evolve/reports/"],
          patterns: ["P3-unsurfaced-value"],
        });
      }
      if (t.bigFiles.length) {
        found.push({
          key: "bloat-bigfiles",
          deficit: headroom * 0.6,
          priority: 1, horizon: "short",
          title: `Eliminate the tracked files >1MB that hold the bloat file penalty open: ${t.bigFiles.length} file(s)`,
          thesis: `bloat=${parts.bloat.toFixed(4)}; the filePenalty term is 1/(1+bigFiles/20) = ${(1 / (1 + t.bigFiles.length / 20)).toFixed(4)} with bigFiles = ${t.bigFiles.length}: ${t.bigFiles.map(([f, s]) => `${f} ${s}`).join(", ")}. Each big file costs ~${(0.5 * (1 / 20) * (1 / (1 + t.bigFiles.length / 20) ** 2)).toFixed(4)} of R, ${pct(0.5 * (1 - 1 / (1 + t.bigFiles.length / 20)) / headroom)} of the live deficit. Total tracked is ${(t.bytes / 1048576).toFixed(1)}MB over ${t.files} files, so the size term is nearly satisfied and the penalty is entirely about ${t.bigFiles.length} outlier file(s).`,
          acceptance: `Zero tracked files exceed 1MB (evolve.mjs readTrackedStats bigFiles == []) OR each remaining one is a regenerable artifact with a recorded fetch manifest. \`node evolve/evolve.mjs gate\` reports filePenalty == 1.`,
          read_only: ["data/README.md", "data/manifest.json", "scripts/fetch_coreset.mjs", "evolve/evolve.mjs"],
          touch: ["data/", "scripts/", ".gitignore"],
          patterns: ["P3-unsurfaced-value"],
        });
      }
    }
  }

  // -- 1d. The gate verdict is computed but never enforced: `gate` prints
  //       "CLOSED (roll back)" and still exits 0, so CI can never go red.
  if (lastVec && lastVec.score !== undefined && lastVec.score < state.R_best) {
    found.push({
      key: "gate-exit-code",
      deficit: 0.9,
      priority: 0, horizon: "short",
      title: `The gate reports a regression and still exits 0: score ${lastVec.score.toFixed(4)} < R_best ${state.R_best}`,
      thesis: `Measured: \`node evolve/evolve.mjs gate\` printed verdict "CLOSED (roll back)" at score ${lastVec.score.toFixed(4)} against R_best ${state.R_best}, and \`echo $?\` was 0. .github/workflows/check.yml runs the gate as a build step, so a step that can only ever pass gates nothing: the loop's own "reject" signal is discarded by the shell. A gate that cannot fail is not a gate.`,
      acceptance: `node evolve/evolve.mjs gate exits non-zero when score < R_best and 0 when score >= R_best; a planted regression makes the CI step red. Proven by a test that sets R_best above the current score and asserts a non-zero exit.`,
      read_only: ["evolve/evolve.mjs", ".github/workflows/check.yml", "package.json"],
      touch: ["evolve/evolve.mjs", "tests/", ".github/workflows/check.yml"],
      patterns: ["P3-unsurfaced-value", "P6-node-gates-blind-to-dom"],
    });
  }

  // -- 1e. The loop's own dispatch order must be a total order over the
  //       goals. Proven by permuting the raw queue arrays and re-sorting.
  const raw = [...queue.goals, ...generated.goals];
  const withCycle = raw.filter((g) => g.cycle !== undefined).length;
  const isOpen = (g) => !g.status || g.status === "open";
  const canonical = raw.filter(isOpen).slice().sort(goalCmp).map((g) => g.id).join(",");
  // A deterministic permutation (no RNG, so the test is reproducible).
  const permuted = raw.filter(isOpen)
    .map((g, i) => ({ g, k: ((i * 7 + 3) * 2654435761) % 100003 }))
    .sort((a, b) => a.k - b.k).map((x) => x.g)
    .sort(goalCmp).map((g) => g.id).join(",");
  if (canonical !== permuted) {
    found.push({
      key: "queue-order-undefined",
      deficit: 0.85,
      priority: 0, horizon: "short",
      title: "Queue order is not a function of the queue: permuting the same goals changes which one `next` hands out",
      thesis: `Measured: sorting the ${raw.filter(isOpen).length} open goals two ways — once in file order, once after a deterministic permutation — yields different sequences (${canonical}  vs  ${permuted}). An Array.sort comparator that returns NaN for some pairs is inconsistent, so the output depends on array layout rather than on queue contents; nextGoal()'s comparator used to end in \`a.cycle - b.cycle\` and ${withCycle} of ${raw.length} goals carry a \`cycle\` field, which is exactly that failure. The loop therefore decides which goal an agent gets by accident, and two identical backlogs can hand out different work.`,
      acceptance: `nextGoal()'s ordering is a total order over the goals: for every permutation of the same goal set, \`node evolve/evolve.mjs next\` returns the same goal id. Proven by a test that shuffles both queue arrays and asserts the emitted id is invariant.`,
      read_only: ["evolve/evolve.mjs", "evolve/queue/goals.json", "evolve/queue/generated.json"],
      touch: ["evolve/evolve.mjs", "tests/"],
      patterns: ["P3-unsurfaced-value"],
    });
  }

  // -- 1f. The loop's only live gradient has a NEGATIVE derivative with
  //       respect to the work the loop exists to do. This is the sharpest
  //       thing the gate vector can say, and it is measurable exactly.
  if (parts && parts.bloat !== undefined && parts.bloat < SAT) {
    const t = readTrackedStats();
    // locScore = 20000/L for L > 20000, so d(locScore)/dL = -20000/L^2 and
    // dR/dL = 0.25 * that. Per 100 lines of added code.
    const dPer100 = -0.25 * 20000 / Math.max(1, t.srcLoc * t.srcLoc) * 100;
    if (dPer100 < 0) {
      found.push({
        key: "gate-gradient-inverted",
        deficit: 0.95,
        priority: 0, horizon: "mid",
        title: `The gate's only live gradient punishes the work the loop exists to do: dR/dLOC = ${dPer100.toFixed(4)} R per 100 LOC`,
        thesis: `Measured: bloat is the only non-saturated component (all others are pinned at 1.0), and its src-LOC term is 1/(1+(L-20000)/20000) with L = ${t.srcLoc}. The derivative is dR/dL = -0.25 x 20000/L^2 = ${(dPer100 / 100).toExponential(3)} per LOC, i.e. ${dPer100.toFixed(4)} R per 100 lines. Two consequences, both measured: (1) writing correct code LOWERS R, so a goal that improves the science and adds 500 LOC is scored as a regression; (2) the composite on the CURRENT, UNMODIFIED tree is ${parts.bloat.toFixed(4)} while R_best is ${state.R_best}, so \`node evolve/evolve.mjs gate\` prints "CLOSED (roll back)" for a no-op — the loop is currently in a state where nothing it evaluates can be accepted. An incentive gradient that points away from the loop's purpose is worse than no gradient, because it is trusted.`,
        acceptance: `The composite is monotonically non-decreasing in work that improves the project: a change that raises any other component and adds LOC does not lower R. Demonstrated by a test that evaluates the score function at srcLoc and srcLoc+1000 with all other components held fixed and asserts score(L+1000) >= score(L); OR by replacing the LOC term with a justified one and re-baselining R_best with an updated rBaselineNote.`,
        read_only: ["evolve/evolve.mjs", "evolve/queue/state.json", "docs/PERFORMANCE.md", "CHANGELOG.md", "ROADMAP.md"],
        touch: ["evolve/evolve.mjs", "tests/", "evolve/reports/"],
        patterns: ["P3-unsurfaced-value"],
      });
    }
  }

  // -- 1g. No recorded vector history at all: a regression is invisible
  //       because there is nothing to regress FROM. (First run only.)
  if (!state.gateVector) {
    found.push({
      key: "gate-no-vector-history",
      deficit: 0.8,
      priority: 0, horizon: "short",
      title: `The loop stores no gate vector, so a regression cannot be detected at all`,
      thesis: `Measured: state.json contains keys [${Object.keys(state).join(", ")}]. There is no gateVector and no prevGateVector, so \`plan\` has no baseline and evolve/evolve.mjs holds no record of what the last score was composed of. R_best is a single scalar, which cannot say WHICH component moved: a bloat regression and a science regression both read as "score went down". The loop's most valuable signal is therefore unavailable by construction.`,
      acceptance: `state.json carries the full component vector for the current and the previous gate run, and \`node evolve/evolve.mjs gate\` diffs them so a per-component regression is reported by name.`,
      read_only: ["evolve/evolve.mjs", "evolve/queue/state.json", "evolve/traces/"],
      touch: ["evolve/evolve.mjs", "tests/"],
      patterns: ["P3-unsurfaced-value"],
    });
  }

  // -- 1g. The queue must never run dry. Measured, not assumed: if the
  //       authored queue is down to its last goal and the generated queue
  //       is empty, \`next\` returns {done:true} and the loop stops.
  const authoredOpen = openGoals(queue.goals).length;
  const genOpen = openGoals(generated.goals).length;
  if (authoredOpen + genOpen <= 1) {
    found.push({
      key: "queue-dry",
      deficit: 1.0,
      priority: 0, horizon: "short",
      title: `The queue is one goal from empty: ${authoredOpen} authored + ${genOpen} generated open`,
      thesis: `Measured: openGoals(goals.json) = ${authoredOpen}, openGoals(generated.json) = ${genOpen}. At ${authoredOpen + genOpen} open, \`node evolve/evolve.mjs next\` returns {"done": true} after the next closure and the loop stops producing work. The queue is the loop's fuel; running it dry is a silent stop, not an error.`,
      acceptance: `\`node evolve/evolve.mjs plan\` has been run so the generated queue holds a measured ${PLAN_SIZE}, and status reports open.authored + open.generated > ${PLAN_SIZE}.`,
      read_only: ["evolve/queue/goals.json", "evolve/queue/generated.json", "evolve/queue/state.json"],
      touch: ["evolve/queue/", "evolve/evolve.mjs"],
      patterns: ["P3-unsurfaced-value"],
    });
  }

  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 2 — the gate's own blind spots.
 *
 * A component that reads 1.0 because it measures the WRONG POPULATION is
 * worse than a component that reads low: it reports success it has not
 * earned. Each probe below compares what the gate counts against what
 * actually exists.
 * ------------------------------------------------------------------ */
function probeGateBlindSpots() {
  const found = [];

  // -- 2a. syntax: the gate walks src/ only. The loop's OWN source is not
  //       syntax-checked by anything that runs in CI.
  const gateJs = listFiles(path.join(ROOT, "src"), /\.js$/);
  const outsideJs = trackedFiles()
    .filter((f) => /\.m?js$/.test(f) && !f.startsWith("src/") && !f.startsWith("node_modules/"))
    .filter((f) => !f.startsWith("tests/"));
  const loopSrc = outsideJs.filter((f) => f.startsWith("evolve/") || f.startsWith("scripts/"));
  if (loopSrc.length) {
    found.push({
      key: "gate-syntax-scope",
      deficit: loopSrc.length / (gateJs.length + loopSrc.length),
      priority: 0, horizon: "short",
      title: `The syntax gate checks ${gateJs.length} files under src/ and ZERO of the ${loopSrc.length} .mjs files that run the loop itself`,
      thesis: `evolve.mjs:${"runGate"} walks path.join(ROOT,"src") for *.js and scripts/check.sh does \`find src -type f -name '*.js'\` — both scoped to src/. The ${loopSrc.length} non-src modules (${loopSrc.slice(0, 6).join(", ")}${loopSrc.length > 6 ? ", ..." : ""}) are never node --check'd by any gate, including evolve/evolve.mjs itself. The loop that gates the project is ungated: a syntax error in it makes the loop silently do nothing while CI stays green. Measured: ${gateJs.length} gated vs ${loopSrc.length} ungated.`,
      acceptance: `scripts/check.sh (and the gate's syntax component) covers every tracked .js/.mjs outside node_modules/ and tests/ — at minimum all of: ${loopSrc.slice(0, 4).join(", ")}. A planted syntax error in evolve/evolve.mjs makes \`npm run check\` exit non-zero.`,
      read_only: ["scripts/check.sh", "evolve/evolve.mjs", ".github/workflows/check.yml", "package.json"],
      touch: ["scripts/check.sh", "evolve/evolve.mjs", "tests/", ".github/workflows/check.yml"],
      patterns: ["P3-unsurfaced-value", "P6-node-gates-blind-to-dom"],
    });
  }

  // -- 2b. dom: the gate parses src/ui.js for $("id") only. Every other
  //       module in src/ references DOM ids the gate never sees.
  try {
    const html = fs.readFileSync(path.join(ROOT, "index.html"), "utf-8");
    const htmlIds = new Set([...html.matchAll(/id="([^"]+)"/g)].map((m) => m[1]));
    const ui = fs.readFileSync(path.join(ROOT, "src", "ui.js"), "utf-8");
    const gateIds = new Set([...ui.matchAll(/\$\("([^"]+)"\)/g)].map((m) => m[1]));
    const all = new Map();
    for (const f of listFiles(path.join(ROOT, "src"), /\.js$/)) {
      const s = fs.readFileSync(path.join(ROOT, f), "utf-8");
      for (const m of s.matchAll(/(?:\$|getElementById|querySelector)\(\s*"#?([A-Za-z][\w-]*)"/g)) {
        if (!all.has(m[1])) all.set(m[1], new Set());
        all.get(m[1]).add(f);
      }
    }
    const blind = [...all.keys()].filter((id) => !gateIds.has(id));
    if (blind.length) {
      found.push({
        key: "gate-dom-blind",
        deficit: blind.length / all.size,
        priority: 0, horizon: "mid",
        title: `The DOM gate checks ${gateIds.size} ids from src/ui.js and misses ${blind.length} of the ${all.size} ids the app actually references`,
        thesis: `dom=${(state.gateVector?.parts?.dom ?? 0).toFixed(4)} is saturated, but the gate derives its id set from ONE regex over src/ui.js: /\\$\\("([^"]+)"\\)/g. Scanning every src/**/*.js for $("id"), getElementById("id") and querySelector("#id") finds ${all.size} distinct ids, of which ${blind.length} are never checked: ${blind.slice(0, 8).join(", ")}${blind.length > 8 ? ", ..." : ""}. A typo in any of those ${blind.length} ids is invisible to the gate AND to scripts/wikiskill_gate.js, so dom=1.0 reports coverage the project does not have.`,
        acceptance: `The dom component of the gate (and scripts/wikiskill_gate.js) is derived from a full scan of src/**/*.js, not src/ui.js alone; every one of the ${all.size} referenced ids is checked against index.html. A test plants a bad id in a non-ui.js module and asserts the gate catches it.`,
        read_only: ["src/ui.js", "index.html", "scripts/wikiskill_gate.js", "src/main.js", "src/analysis-panel.js"],
        touch: ["evolve/evolve.mjs", "scripts/wikiskill_gate.js", "tests/", "index.html"],
        patterns: ["P5-ui-contract-fragility", "P6-node-gates-blind-to-dom"],
      });
    }
  } catch {}

  // -- 2c. coverage: the denominator is tests/test_*.js only, so the
  //       scripts/ validators (including the one that FAILS) are invisible.
  const tCount = listFiles(path.join(ROOT, "tests"), /^test_.*\.js$/);
  const sCount = listFiles(path.join(ROOT, "scripts"), /^(test|validate)_.*\.mjs$/);
  const hidden = sCount.filter((f) => !f.startsWith("tests/manual/"));
  if (hidden.length) {
    found.push({
      key: "gate-coverage-blind",
      deficit: hidden.length / (tCount.length + hidden.length),
      priority: 1, horizon: "mid",
      title: `coverage=${(state.gateVector?.parts?.coverage ?? 0).toFixed(4)} is saturated by counting only tests/test_*.js — ${hidden.length} scripts/ validators are outside the denominator`,
      thesis: `countTests() globs the tests/ directory for /^test_.*\\.js$/ and finds ${tCount.length}, all wired, so coverage = (${tCount.length}+0)/${tCount.length} = 1.0 and is reported saturated. The scripts/ directory holds ${sCount.length} more test-bearing files (${hidden.join(", ")}) that no component counts. One of them, scripts/validate_binding_physics_r1.mjs, currently exits 1 — the gate is saturated on a population that excludes the failing member. Measured coverage of the real population is ${tCount.length}/${tCount.length + hidden.length} = ${pct(tCount.length / (tCount.length + hidden.length))}, not 100%.`,
      acceptance: `The coverage component counts every test-bearing file (tests/test_*.js AND scripts/{test,validate}_*.mjs not in tests/manual/), and the number it reports equals (wired+manual)/total over that full population. A validator moved into scripts/ can change the score.`,
      read_only: ["evolve/evolve.mjs", "tests/suites.js", "tests/test_suite_registry.js", "evolve/reports/test-inventory.md"],
      touch: ["evolve/evolve.mjs", "tests/", "evolve/reports/"],
      patterns: ["P3-unsurfaced-value", "P6-node-gates-blind-to-dom"],
    });
  }

  // -- 2d. science: the overclaim check reads only the first 4000 chars of a
  //       22k-char README, and passes iff a regex matches. 82% of the
  //       honest-scope surface is unexamined, and the check is a keyword hit.
  try {
    const readme = fs.readFileSync(path.join(ROOT, "README.md"), "utf-8");
    const WINDOW = 4000;
    const scanned = Math.min(WINDOW, readme.length);
    if (readme.length > WINDOW) {
      found.push({
        key: "gate-science-window",
        deficit: 1 - scanned / readme.length,
        priority: 1, horizon: "short",
        title: `science=${(state.gateVector?.parts?.science ?? 0).toFixed(4)} is saturated, but the overclaim check reads only ${scanned} of README.md's ${readme.length} chars and passes on a keyword regex`,
        thesis: `scienceScore() computes have/must for 7 files, then multiplies by (overclaim ? 1 : 0.5) where overclaim = /NOT a replacement|qualitative|not FEP|trust boundary/i.test(README.md.slice(0,${WINDOW})). Measured: README.md is ${readme.length} chars, so ${readme.length - scanned} chars (${pct(1 - scanned / readme.length)}) of the project's honest-scope surface is never examined, and the check that IS run is satisfied by a single matching phrase rather than by the absence of an overclaim. Deleting the phrase halves the component; writing a false claim anywhere past char ${WINDOW} changes nothing.`,
        acceptance: `The overclaim check examines 100% of README.md (and the docs it summarises) and asserts the ABSENCE of specific overclaim patterns rather than the presence of a disclaimer phrase. A test plants a known overclaim sentence and asserts science < 1.0.`,
        read_only: ["evolve/evolve.mjs", "README.md", "ROADMAP.md", "docs/LIMITATIONS.md", "docs/APPLICABILITY.md"],
        touch: ["evolve/evolve.mjs", "tests/", "README.md"],
        patterns: ["P3-unsurfaced-value"],
      });
    }
  } catch {}

  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 3 — the standing failure the test-inventory goal left behind.
 *
 * A file that is RED and that nothing runs is a silent hole. A file that is
 * merely SLOW is not: this probe must never convert a timeout into evidence
 * of failure, because that manufactures a goal out of nothing. Each file
 * therefore gets a bounded budget, and a budget overrun is recorded as
 * INCONCLUSIVE (contributing no goal) rather than as a failure.
 * ------------------------------------------------------------------ */

/** Parse tests/suites.js for { file, tier, timeout } without importing it —
 *  a regex keeps `plan` dependency-free and the file is a plain array. */
function readSuiteRegistry() {
  const reg = new Map();
  let body = "";
  try { body = fs.readFileSync(path.join(ROOT, "tests", "suites.js"), "utf-8"); } catch { return reg; }
  for (const m of body.matchAll(/file:\s*"([^"]+)"\s*,\s*tier:\s*"([A-Z]+)"\s*(?:,\s*timeout:\s*(\d+))?/g)) {
    reg.set(m[1], { tier: m[2], timeout: m[3] ? Number(m[3]) : null });
  }
  return reg;
}

function probeUngatedFailures() {
  const found = [];
  const tracked = new Set(trackedFiles());
  const registry = readSuiteRegistry();

  for (const f of listFiles(path.join(ROOT, "scripts"), /^(test|validate)_.*\.mjs$/)) {
    const entry = registry.get(f);
    const inSuites = !!entry;
    const inCi = [...tracked].some((t) => /\.ya?ml$/.test(t) && (() => {
      try { return fs.readFileSync(path.join(ROOT, t), "utf-8").includes(path.basename(f)); } catch { return false; }
    })());
    const gated = inSuites || inCi;

    // A file the registry already declares long-running is NOT probed: the
    // registry gates its registration, and a bounded probe cannot fairly
    // judge a suite whose own budget is 900s. Recorded, not guessed at.
    if (entry?.timeout && entry.timeout > 120000) {
      INCONCLUSIVE.push({ probe: "ungated-failure", subject: f, why: `registered at tests/suites.js with a ${entry.timeout}ms budget (tier ${entry.tier}); too slow for a bounded probe, and the registry already gates its registration — NO goal generated` });
      continue;
    }

    // Bounded budget. `plan` must stay fast, and an overrun is not a failure.
    const budgetMs = 60000;
    let code = 0, tail = "", timedOut = false;
    try {
      const out = execFileSync(process.execPath, [path.join(ROOT, f)], {
        cwd: ROOT, encoding: "utf-8", stdio: ["ignore", "pipe", "pipe"], timeout: budgetMs,
      });
      code = 0; tail = String(out).split("\n").filter(Boolean).slice(-1)[0] || "";
    } catch (e) {
      if (e.killed || e.signal === "SIGTERM") timedOut = true;
      else { code = typeof e.status === "number" ? e.status : 1; tail = String(e.stderr || e.message).split("\n").filter(Boolean).slice(-1)[0] || ""; }
    }
    if (timedOut) { INCONCLUSIVE.push({ probe: "ungated-failure", subject: f, why: `exceeded the ${budgetMs / 1000}s probe budget (registered timeout ${entry?.timeout ?? "none"}ms) — a slow file is not a failing file, so NO goal was generated` }); continue; }
    if (code === 0) continue;

    const where = inSuites
      ? `It IS registered in tests/suites.js at tier ${entry.tier}${entry.timeout ? ` (timeout ${entry.timeout}ms)` : ""}, which only runs on an opt-in flag that neither \`npm test\` nor any CI step passes, so the gate's default run never sees it.`
      : "It is NOT registered in tests/suites.js.";
    const ci = inCi ? "It IS named in a CI workflow." : "No CI workflow names it.";
    found.push({
      key: `failing-${path.basename(f)}`,
      deficit: gated ? 0.75 : 1.0,
      priority: 0, horizon: "short",
      title: `${f} exits ${code} and ${gated ? "sits behind a flag" : "is wired into NO gate that can fail"}`,
      thesis: `Measured right now: \`node ${f}\` exits ${code} inside the ${budgetMs / 1000}s probe budget. Last line: "${tail.slice(0, 160)}". ${where} ${ci} ${gated ? "" : "The gate's coverage component counts only tests/test_*.js, so this file is outside its denominator entirely."} A red validator the default gate never runs is indistinguishable from no validator — which is exactly how evolve/reports/test-inventory.md recorded it as a standing finding and left it standing.`,
      acceptance: `\`node ${f}\` exits 0, OR the file is moved to a tier that \`npm test\`/CI actually executes (justified by its measured runtime). The failing assertion is fixed at its source — not weakened, not deleted. \`node evolve/evolve.mjs gate\` goes red if it regresses again.`,
      read_only: [f, "evolve/reports/test-inventory.md", "tests/suites.js", "package.json"],
      touch: ["scripts/", "tests/suites.js", "src/", "docs/"],
      patterns: ["P3-unsurfaced-value", "P5-ui-contract-fragility"],
    });
  }
  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 4 — budget drift. bench/budget.json is a claim the project
 * repeats in docs and CI; this measures whether it is still true.
 * ------------------------------------------------------------------ */
function probeBudgetDrift(opts = {}) {
  const found = [];
  let budget;
  try { budget = readJson(path.join(ROOT, "bench", "budget.json"), null); } catch { budget = null; }
  if (!budget) return found;

  // A benchmark is the one measurement here that is genuinely noisy, and a
  // noisy measurement that changes the RANKING would make `plan`
  // non-idempotent for no good reason. So the result is cached against the
  // tree it was measured on: same HEAD+worktree => reuse the number; any
  // change => re-measure. Re-measure with --no-bench to force.
  let head = "";
  try {
    // stderr is piped, not inherited: a repo with no commits yet must not
    // spray git's "fatal: ambiguous argument" through the plan's own report.
    head = execFileSync("git", ["rev-parse", "HEAD"], { cwd: ROOT, encoding: "utf-8", stdio: ["ignore", "pipe", "pipe"] }).trim();
  } catch {}
  const treeKey = `${head}:${readTrackedStats().files}`;

  let perf = state.benchVector;
  const reusable = perf && perf.treeKey === treeKey;
  if (!reusable && !opts.noBench) {
    try {
      const out = execFileSync(process.execPath, [path.join(ROOT, "bench", "perf.js")], {
        cwd: ROOT, encoding: "utf-8", stdio: ["ignore", "pipe", "pipe"], timeout: 300000,
      });
      const hv = /Heavy:[^\n]*?([\d.]+)\s*±/.exec(out);
      const cg = /CG:[^\n]*?([\d.]+)\s*±/.exec(out);
      if (hv) {
        perf = { heavyMs: Number(hv[1]), cgMs: cg ? Number(cg[1]) : null, treeKey, measuredAt: new Date().toISOString() };
        state.benchVector = perf;
        writeJson(STATE, state);
      }
    } catch {}
  } else if (!reusable && opts.noBench) {
    INCONCLUSIVE.push({ probe: "budget-drift", subject: "bench/perf.js", why: "--no-bench and no cached measurement for this tree; budget drift NOT assessed rather than guessed" });
    return found;
  }
  if (!perf?.heavyMs) return found;

  const b = budget.heavy_compute_ms;
  const ratio = perf.heavyMs / b;
  if (ratio > 1.5) {
    found.push({
      key: "budget-drift-heavy",
      deficit: Math.min(1, (ratio - 1) / 7),
      priority: 1, horizon: "mid",
      title: `bench/budget.json claims heavy_compute_ms ${b}; measured ${perf.heavyMs} ms/compute (${ratio.toFixed(1)}x over)`,
      thesis: `Measured: \`node bench/perf.js\` reports Heavy ${perf.heavyMs} ms/compute on 4w52.pdb (n=1308, 20 warmup + 30 timed). bench/budget.json declares heavy_compute_ms ${b}. The budget is ${ratio.toFixed(2)}x tighter than reality and has been for long enough that docs/PERFORMANCE.md and .github/workflows/check.yml both cite it as the performance contract while CI emits only a ::warning on a permanent ${ratio.toFixed(1)}x overrun. A budget nobody believes is worse than no budget: it makes the number meaningless. Either the code gets ${ratio.toFixed(1)}x faster or the number becomes the measured one.`,
      acceptance: `bench/budget.json heavy_compute_ms is within 1.2x of \`node bench/perf.js\` output on 4w52.pdb, and docs/PERFORMANCE.md + .github/workflows/check.yml quote the same value. If the physics cannot be made to fit, the budget is rebaselined to the measured p95 with a written justification and the CI step becomes a warn-on-regression-from-baseline rather than a permanent ignore.`,
      read_only: ["bench/budget.json", "bench/perf.js", "docs/PERFORMANCE.md", ".github/workflows/check.yml", "src/heavy.js", "src/spatial-grid.js"],
      touch: ["bench/budget.json", "docs/PERFORMANCE.md", ".github/workflows/check.yml", "manuscript/reproduce.sh"],
      patterns: ["P3-unsurfaced-value"],
    });
  }

  // fps is declared in the budget and measured by nothing at all.
  if (budget.fps && !/fps/.test(fs.readFileSync(path.join(ROOT, "bench", "perf.js"), "utf-8"))) {
    found.push({
      key: "budget-fps-unmeasured",
      deficit: 0.5,
      priority: 1, horizon: "mid",
      title: `bench/budget.json declares fps ${budget.fps}; no file under bench/ or scripts/ measures frame rate at all`,
      thesis: `budget.json carries three keys: heavy_compute_ms, cg_compute_ms, fps ${budget.fps}. Measured: \`grep -n fps bench/*.js scripts/*.mjs tests/*.js\` returns ZERO hits — bench/perf.js measures ms/compute and nothing measures frames. .github/workflows/check.yml prints "fps budget 30 — check render fps via bench/perf.js scale if applicable" and manuscript/reproduce.sh does the same: both reference a measurement that does not exist. One third of the declared performance contract is unfalsifiable.`,
      acceptance: `Either a real fps measurement exists (a bench that reports frames/second over a fixed workload, printing a number CI can compare) or the fps key is removed from bench/budget.json and from every doc/script that quotes it. \`grep -rn "budget\[.fps.\]\\|budget.get(.fps.)\" bench/ scripts/ manuscript/ .github/\` returns no unbacked reference.`,
      read_only: ["bench/budget.json", "bench/perf.js", "docs/PERFORMANCE.md", "manuscript/reproduce.sh", ".github/workflows/check.yml", "src/viewer.js"],
      touch: ["bench/budget.json", "bench/", "docs/PERFORMANCE.md", "manuscript/reproduce.sh", ".github/workflows/check.yml"],
      patterns: ["P3-unsurfaced-value"],
    });
  }
  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 5 — declared-but-unimplemented. The grep for TODO/FIXME/stub is
 * only useful if each hit is CLASSIFIED; this measures how many hits are
 * real debt (a code path that returns a stub value at runtime) versus an
 * honest scope note (a comment saying a thing is not implemented, which is
 * the project's own documented virtue per ROADMAP.md §1).
 * ------------------------------------------------------------------ */
function probeStubDebt() {
  const found = [];
  const re = /TODO|FIXME|not yet|placeholder|stub/i;
  const hits = [];
  for (const f of listFiles(path.join(ROOT, "src"), /\.js$/)) {
    const lines = fs.readFileSync(path.join(ROOT, f), "utf-8").split("\n");
    lines.forEach((l, i) => {
      if (!re.test(l)) return;
      const t = l.trim();
      // Classification: a hit on a line that EXECUTES is potential real
      // debt. A hit inside a comment/leading-* is a scope note by design.
      const isComment = t.startsWith("//") || t.startsWith("*") || t.startsWith("/*") || t.startsWith("#");
      hits.push({ f, line: i + 1, text: t, isComment });
    });
  }
  const code = hits.filter((h) => !h.isComment);
  const comments = hits.length - code.length;
  if (code.length) {
    found.push({
      key: "stub-code-paths",
      deficit: Math.min(1, code.length / 10),
      priority: 1, horizon: "mid",
      title: `${code.length} of ${hits.length} TODO/stub/placeholder hits in src/ are on lines that EXECUTE, not in comments`,
      thesis: `Measured by classify-each-hit: \`grep -rn "TODO|FIXME|not yet|placeholder|stub" src/\` returns ${hits.length} hits across ${new Set(hits.map((h) => h.f)).size} files. ${comments} are comment/leading-* lines — honest scope notes, which ROADMAP.md §1 makes a virtue and which must be left alone. ${code.length} are on lines that execute: ${code.map((h) => `${h.f}:${h.line}`).join(", ")}. Those are the ones that can change a number a user sees (e.g. a returned note string containing "stub:", or a throw whose message says "not yet supported"). Classify-by-grep is not currently a test, so the ${comments} honest notes and the ${code.length} live ones are indistinguishable to CI.`,
      acceptance: `A test classifies every hit: comment-line hits must be recognised as scope notes, and every code-line hit must either be removed or be backed by a named, dated, honest-scope entry. The classification is asserted, not eyeballed. \`node ${"tests/test_stub_classification.js"}\` (or the chosen filename) fails if a new code-line stub appears.`,
      read_only: ["src/funnel.js", "src/mmcif.js", "src/physics/network.js", "src/heavy.js", "src/spatial-grid.js", "ROADMAP.md"],
      touch: ["tests/", "src/", "docs/LIMITATIONS.md"],
      patterns: ["P3-unsurfaced-value"],
    });
  }

  // Write-only state: assigned, never read anywhere in the repo. This is
  // the sharpest form of "implemented in a comment, dead in reality".
  const srcFiles = listFiles(path.join(ROOT, "src"), /\.js$/);
  const corpus = srcFiles.map((f) => [f, fs.readFileSync(path.join(ROOT, f), "utf-8")]);
  const orphanFields = [];
  for (const [f, s] of corpus) {
    for (const m of s.matchAll(/this\.(_[A-Za-z]\w*)\s*=\s*[^=]/g)) {
      const field = m[1];
      let uses = 0;
      for (const [, t] of corpus) uses += (t.match(new RegExp(`\\b${field}\\b`, "g")) || []).length;
      if (uses <= 1) orphanFields.push(`${f} this.${field}`);
    }
  }
  if (orphanFields.length) {
    found.push({
      key: "write-only-state",
      deficit: Math.min(1, orphanFields.length / 5),
      priority: 1, horizon: "short",
      title: `${orphanFields.length} private fields are ASSIGNED in src/ and never read anywhere: ${orphanFields.slice(0, 5).join(", ")}`,
      thesis: `Measured: for every \`this._x = ...\` in src/**/*.js, count occurrences of \`_x\` across all of src/. Fields with exactly one occurrence (the assignment) are write-only. Found ${orphanFields.length}: ${orphanFields.join(", ")}. Each is a value the code computes and then discards — e.g. a PMF-derived barrier that is stored and never used to rebuild anything, so a function named for an effect has no effect. A reader (or the next agent) sees the assignment and infers the behaviour exists.`,
      acceptance: `Every private field assigned in src/ is read at least once, or the assignment is removed. A test enumerates \`this._x = \` assignments in src/ and fails on any with zero reads, so the class cannot regrow silently.`,
      read_only: ["src/physics/network.js", "src/funnel.js", "src/physics/weakint.js", "src/viewer.js", "docs/NETWORK.md"],
      touch: ["src/", "tests/", "docs/NETWORK.md"],
      patterns: ["P3-unsurfaced-value"],
    });
  }
  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 6 — documentation that cites code it can no longer find.
 * A `file:line` claim is checkable, so it is checked: an unresolvable one
 * is a false statement about the code, which is worse than no statement.
 * ------------------------------------------------------------------ */
function probeDocCitations() {
  const found = [];
  const files = trackedFiles();
  if (!files.length) return found;
  const byBase = new Map();
  for (const f of files) {
    const b = path.basename(f);
    if (!byBase.has(b)) byBase.set(b, []);
    byBase.get(b).push(f);
  }
  const docs = files.filter((f) => f.endsWith(".md"));
  const stale = [];
  let checked = 0;
  for (const f of docs) {
    let lines;
    try { lines = fs.readFileSync(path.join(ROOT, f), "utf-8").split("\n"); } catch { continue; }
    lines.forEach((line, i) => {
      for (const m of line.matchAll(/\b([A-Za-z0-9_.\/-]+\.(?:js|mjs|json|md|html|css|sh))[:](\d+)/g)) {
        const target = m[1], ln = Number(m[2]);
        checked++;
        const cands = target.includes("/")
          ? [target, path.join(path.dirname(f), target), path.join("src", target), path.join("docs", target),
             path.join("scripts", target), path.join("bench", target), path.join("tests", target)]
          : [...(byBase.get(target) || []), path.join(path.dirname(f), target), path.join("src", target)];
        const hit = cands.find((c) => c && fs.existsSync(path.join(ROOT, c)));
        if (!hit) { stale.push(`${f}:${i + 1} -> ${target}:${ln} (no such file)`); continue; }
        const n = fs.readFileSync(path.join(ROOT, hit), "utf-8").split("\n").length;
        if (ln > n) stale.push(`${f}:${i + 1} -> ${hit}:${ln} (file has ${n} lines)`);
      }
    });
  }
  if (stale.length) {
    found.push({
      key: "doc-citation-drift",
      deficit: Math.min(1, stale.length / 20),
      priority: 2, horizon: "short",
      title: `${stale.length} of ${checked} \`file:line\` citations across ${docs.length} markdown files point at nothing`,
      thesis: `Measured by resolving every \`name.ext:NN\` citation in every tracked *.md against the working tree: ${checked} citations checked, ${stale.length} unresolvable — ${stale.slice(0, 4).join("; ")}${stale.length > 4 ? "; ..." : ""}. These are the load-bearing sentences of the project's validation record ("all §0 code claims machine-verified"), and one of the failures is a citation to a file that a previous evolution DELETED, so the record now describes code that does not exist. A claim that cannot be resolved is a claim nobody can check, and a doc that is wrong about the code is worse than a doc that is silent.`,
      acceptance: `Every \`file:line\` citation in tracked *.md resolves to an existing file with at least that many lines. The check runs as a test in the FAST tier (or tests/manual/ with a README entry), so drift fails a gate rather than being found by a reader.`,
      read_only: [...stale.map((s) => s.split(" -> ")[0].split(":")[0]).filter((v, i, a) => a.indexOf(v) === i).slice(0, 4), "docs/", "ROADMAP.md"],
      touch: ["docs/", "tests/"],
      patterns: ["P3-unsurfaced-value"],
    });
  }
  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 7 — size/split candidates. Only files whose size is MEASURED to be
 * an outlier relative to the repo's own distribution; no "this file feels
 * big" heuristics.
 * ------------------------------------------------------------------ */
function probeSizeOutliers() {
  const found = [];
  const groups = [
    { dir: "src", re: /\.js$/, floor: 900, label: "src modules" },
    { dir: "docs", re: /\.md$/, floor: 400, label: "docs" },
  ];
  const outliers = [];
  for (const { dir, re, floor, label } of groups) {
    for (const f of listFiles(path.join(ROOT, dir), re)) {
      const loc = fs.readFileSync(path.join(ROOT, f), "utf-8").split("\n").length;
      if (loc >= floor) outliers.push({ f, loc, dir, label });
    }
  }
  if (!outliers.length) return found;
  const srcOut = outliers.filter((o) => o.dir === "src").sort((a, b) => b.loc - a.loc);
  const top = srcOut[0];
  if (top) {
    const med = (() => {
      const all = listFiles(path.join(ROOT, "src"), /\.js$/)
        .map((f) => fs.readFileSync(path.join(ROOT, f), "utf-8").split("\n").length).sort((a, b) => a - b);
      return all[Math.floor(all.length / 2)] || 0;
    })();
    found.push({
      key: "src-size-outlier",
      deficit: Math.min(1, (top.loc - med) / 1500),
      priority: 1, horizon: "mid",
      title: `${top.f} is ${top.loc} LOC against a src/ median of ${med} — ${(top.loc / Math.max(1, med)).toFixed(1)}x the typical module`,
      thesis: `Measured: line counts of every src/**/*.js. Largest is ${srcOut.slice(0, 5).map((o) => `${o.f} ${o.loc}`).join(", ")}; median is ${med} LOC. ${top.f} at ${top.loc} LOC is ${(top.loc / Math.max(1, med)).toFixed(1)}x the median and is a single point of syntactic failure: the project's own AUDIT_REPORT.md §1.1 names the god-module failure mode, and one unparseable line inside it takes down every importer. Split by responsibility so the failure domain shrinks.`,
      acceptance: `${top.f} is under 400 LOC of pure wiring, its responsibilities live in >= 4 extracted modules that are each independently node --check-able, no behaviour changes (test count >= the pre-split high-water mark), and the DOM contract test still passes. Documented in the README project layout.`,
      read_only: [top.f, "src/ui.js", "src/integrator.js", "index.html", "README.md", "CHANGELOG.md"],
      touch: [path.dirname(top.f) + "/", "tests/", "README.md"],
      patterns: ["P5-ui-contract-fragility", "P6-node-gates-blind-to-dom"],
    });
  }
  const docOut = outliers.filter((o) => o.dir === "docs").sort((a, b) => b.loc - a.loc);
  if (docOut.length >= 2) {
    const total = docOut.reduce((a, b) => a + b.loc, 0);
    found.push({
      key: "docs-size-outlier",
      deficit: Math.min(1, docOut[0].loc / 2000),
      priority: 2, horizon: "short",
      title: `${docOut.length} docs/ files exceed 400 LOC (largest ${docOut[0].f} at ${docOut[0].loc} LOC, ${total} LOC combined)`,
      thesis: `Measured: line counts of every docs/**/*.md. Over 400 LOC: ${docOut.slice(0, 5).map((o) => `${o.f} ${o.loc}`).join(", ")}${docOut.length > 5 ? ", ..." : ""} — ${total} LOC in ${docOut.length} files. A running log of what was done is not a document a scientist can consult: ROADMAP.md §1 makes the out-of-scope contract the most important document in the repo, and it competes for attention with a ${docOut[0].loc}-line completion log. Split the log from the guidance so the guidance is findable.`,
      acceptance: `Each over-400-LOC doc/ file is split into a durable reference and a dated log (or the log moves to a clearly-named history/ location linked from docs/README.md); docs/README.md still links every file exactly once (tests/test_docs_index.js green) and the ROADMAP/LIMITATIONS/APPLICABILITY/VALIDATION trust-boundary docs remain the first thing listed.`,
      read_only: ["docs/README.md", ...docOut.slice(0, 3).map((o) => o.f), "ROADMAP.md", "tests/test_docs_index.js"],
      touch: ["docs/", "README.md", "tests/test_docs_index.js"],
      patterns: ["P3-unsurfaced-value", "P4-feature-hierarchy"],
    });
  }
  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 8 — the cache-bust sprawl: a measurable, falsifiable drift
 * between the version the code claims and the version the browser is told
 * to fetch.
 * ------------------------------------------------------------------ */
function probeVersionDrift() {
  const found = [];
  let ver = "";
  try {
    ver = (/export const VERSION = "([^"]+)"/.exec(fs.readFileSync(path.join(ROOT, "src", "version.js"), "utf-8")) || [])[1] || "";
  } catch { return found; }
  const files = listFiles(path.join(ROOT, "src"), /\.js$/);
  const literals = new Map();
  let total = 0;
  for (const f of files) {
    for (const m of fs.readFileSync(path.join(ROOT, f), "utf-8").matchAll(/\?v=([\w.]+)/g)) {
      total++;
      if (!literals.has(m[1])) literals.set(m[1], []);
      literals.get(m[1]).push(f);
    }
  }
  const staleV = [...literals.entries()].filter(([v]) => v !== ver);
  if (staleV.length) {
    const files_ = new Set(staleV.flatMap(([, fs_]) => fs_));
    const staleDesc = staleV.map(([v, f]) => `?v=${v} in ${f.length} files`).join("; ");
    found.push({
      key: "version-literal-drift",
      deficit: Math.min(1, files_.size / 20),
      priority: 1, horizon: "mid",
      title: `${total} hardcoded ?v= cache-bust literals across ${files_.size} src/ files disagree with src/version.js VERSION "${ver}"`,
      thesis: `Measured: every ?v= literal in src/**/*.js, counted by value. src/version.js declares VERSION = "${ver}", but ${total} import specifiers carry a literal ?v= and ${staleV.length} distinct literal value(s) disagree with it: ${staleDesc}. A module edited without editing this literal is served from the browser cache, and the symptom is "my fix did nothing" — the single most expensive class of bug to diagnose. The version bump is currently not one edit; it is ${files_.size} edits that nothing checks.`,
      acceptance: `Zero hardcoded ?v= literals remain in src/ (a loader injects the version from src/version.js), OR a check fails when any ?v= disagrees with src/version.js and that check is wired into npm run check and CI. Verified by editing a module, reloading, and observing the change.`,
      read_only: ["src/version.js", "index.html", "scripts/check.sh", "package.json", "CHANGELOG.md"],
      touch: ["src/", "index.html", "scripts/check.sh", ".github/workflows/check.yml"],
      patterns: ["P5-ui-contract-fragility"],
    });
  }
  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 9 — error paths that swallow. Measured, not vibes: a `catch`
 * with an empty body means a failure the user will never see.
 * ------------------------------------------------------------------ */
function probeSilentCatches() {
  const found = [];
  const hits = [];
  for (const f of listFiles(path.join(ROOT, "src"), /\.js$/)) {
    const lines = fs.readFileSync(path.join(ROOT, f), "utf-8").split("\n");
    lines.forEach((l, i) => {
      if (/catch\s*(\([^)]*\))?\s*\{\s*\}/.test(l)) hits.push(`${f}:${i + 1}`);
    });
  }
  if (hits.length) {
    found.push({
      key: "silent-catch",
      deficit: Math.min(1, hits.length / 40),
      priority: 1, horizon: "mid",
      title: `${hits.length} empty catch blocks in src/ swallow every error they see`,
      thesis: `Measured: \`grep -rn "catch.*{\\s*}" src/\` returns ${hits.length} empty catch blocks across ${new Set(hits.map((h) => h.split(":")[0])).size} files — e.g. ${hits.slice(0, 5).join(", ")}${hits.length > 5 ? ", ..." : ""}. A swallowed exception is a failed operation that reports success. A scientist who loads a bad PDB must never see a frozen viewer and conclude the software is broken in an interesting way; right now ${hits.length} code paths can produce exactly that. The count is the acceptance criterion, so this is falsifiable: it either goes to 0 or each one is justified in writing.`,
      acceptance: `grep finds zero empty catch blocks in src/ (currently ${hits.length}); every remaining catch either rethrows, logs with the error object, or increments a visible HUD counter. A test feeds a malformed PDB and asserts a specific visible message, and asserts the swallowed-error counter is 0 after a clean load.`,
      read_only: ["src/viewer.js", "src/main.js", "src/input_errors.js", "src/ui.js", "docs/ACCESSIBILITY.md"],
      touch: ["src/", "index.html", "tests/"],
      patterns: ["P2-dead-canvas-redraw", "P3-unsurfaced-value"],
    });
  }
  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 10 — duplicated physical constants. Two literals for the same
 * physical quantity is the exact mechanism by which a physics fix lands
 * in one path and not the other.
 *
 * The match is CASE-SENSITIVE on purpose: a SCREAMING_CASE name is a
 * declared constant, `inv` is a loop counter. A case-insensitive scan
 * reports ordinary local variables as "duplicated physics", which would be
 * a fabricated goal.
 * ------------------------------------------------------------------ */
function probeConstantDuplication() {
  const found = [];
  const files = listFiles(path.join(ROOT, "src"), /\.js$/);
  const srcByPath = files.map((f) => [f, fs.readFileSync(path.join(ROOT, f), "utf-8")]);
  const idx = new Map();
  for (const [f, s] of srcByPath) {
    for (const m of s.matchAll(/(?:^|[^\w.$])([A-Z][A-Z0-9_]{2,})\s*=\s*(-?\d*\.?\d+(?:[eE]-?\d+)?)/g)) {
      const k = `${m[1]}=${m[2]}`;
      if (!idx.has(k)) idx.set(k, new Set());
      idx.get(k).add(f);
    }
  }
  const dupes = [...idx.entries()].filter(([, s]) => s.size > 1);
  if (dupes.length >= 3) {
    found.push({
      key: "constant-duplication",
      deficit: Math.min(1, dupes.length / 20),
      priority: 1, horizon: "mid",
      title: `${dupes.length} declared constants (SCREAMING_CASE = number) are defined in more than one src/ file`,
      thesis: `Measured: every \`NAME = <number>\` literal in src/**/*.js, indexed by name and value. ${dupes.length} appear in >= 2 files: ${dupes.slice(0, 6).map(([k, s]) => `${k} in ${[...s].join("+")}`).join("; ")}${dupes.length > 6 ? "; ..." : ""}. A duplicated physical constant means the CG path and the heavy path can disagree about the same quantity — the class of defect docs/BINDING_PHYSICS_R* spent seven review rounds finding, one of which (KB_KCAL defined twice) was a real 5.03e-8 divergence caught only because a parity test happened to exist. Each duplicate is a future divergence waiting for someone to edit one side.`,
      acceptance: `Every duplicated declared constant is defined in exactly one module and imported by the others, and a parity test asserts CG and heavy agree to < 1e-12 for each one. The duplication scan itself becomes a test so the class cannot regrow.`,
      read_only: ["src/units.js", "src/physics/observables.js", "src/ff-params.js", "src/heavy.js", "src/forcefield.js", "tests/test_unit_contract.js"],
      touch: ["src/", "tests/"],
      patterns: ["P5-ui-contract-fragility"],
    });
  }
  return found;
}

/* ------------------------------------------------------------------ *
 * PROBE 11 — the authored queue itself. Open goals are carried forward,
 * re-measured: a goal whose numbers have gone stale is a goal that will
 * waste an agent call, because the subagent reads the thesis as fact.
 * ------------------------------------------------------------------ */
function probeOpenAuthoredGoals() {
  const found = [];
  for (const g of openGoals(queue.goals)) {
    const thesis = String(g.thesis || "");
    const cited = [...thesis.matchAll(/(\d[\d,]*)\s*(LOC|lines|files|ms|KB|MB|Kcal|beads|atoms|x|%|days|s)/gi)];
    if (!cited.length) continue;
    // The one claim this planner can cheaply falsify without running the
    // goal: a `file:line` citation inside an authored thesis.
    const bad = [];
    for (const m of thesis.matchAll(/\b([A-Za-z0-9_.\/-]+\.(?:js|mjs|md|json)):(\d+)/g)) {
      const cands = [m[1], path.join("src", m[1]), path.join("scripts", m[1]), path.join("docs", m[1]),
                     path.join(path.dirname("evolve/evolve.mjs"), m[1])];
      if (!cands.some((c) => fs.existsSync(path.join(ROOT, c)))) bad.push(`${m[1]}:${m[2]}`);
    }
    if (bad.length) {
      found.push({
        key: `stale-thesis-${g.id}`,
        deficit: 0.6,
        priority: 0, horizon: "short",
        title: `Authored goal ${g.id} cites source lines that do not exist: ${bad.join(", ")}`,
        thesis: `Goal ${g.id} ("${g.title}") is still open in evolve/queue/goals.json, so \`next\` will hand it to a subagent, and \`next\` copies the thesis into the packet verbatim as fact. Measured: its thesis cites ${bad.join(", ")}, which resolve to no file on disk. A subagent is told the defect lives at a line that has moved or been deleted; it will either fail to find it or, worse, "fix" the wrong code. The queue is only as good as the numbers in it.`,
        acceptance: `Every \`file:line\` citation in every OPEN authored goal's thesis resolves on disk, asserted by a test so a stale citation cannot re-enter the queue. Re-measure the cited numbers in the same pass.`,
        read_only: ["evolve/queue/goals.json", ...bad.map((b) => b.split(":")[0]).filter((f, i, a) => a.indexOf(f) === i).slice(0, 3)],
        touch: ["evolve/queue/goals.json", "tests/"],
        patterns: ["P3-unsurfaced-value"],
      });
    }
  }
  return found;
}

/* ------------------------------------------------------------------ *
 * Assemble. Order of sources is the priority order from the goal brief:
 * the gate vector first (it is the only thing that can say whether the
 * last change helped), then recorded debt.
 * ------------------------------------------------------------------ */
const PROBES = [
  ["gate-vector", probeGateVector],
  ["gate-blind-spot", probeGateBlindSpots],
  ["ungated-failure", probeUngatedFailures],
  ["budget-drift", probeBudgetDrift],
  ["stub-debt", probeStubDebt],
  ["doc-citation-drift", probeDocCitations],
  ["size-outlier", probeSizeOutliers],
  ["version-drift", probeVersionDrift],
  ["silent-catch", probeSilentCatches],
  ["constant-duplication", probeConstantDuplication],
  ["stale-authored-thesis", probeOpenAuthoredGoals],
];

/** Read-only paths must exist. Resolve or drop — a goal that points a
 *  subagent at a nonexistent file burns its whole context budget. */
function resolveReadOnly(paths) {
  const kept = [], dropped = [];
  for (const p of paths || []) {
    // A trailing "#anchor" names a section, not a path; the file must exist.
    const base = String(p).split("#")[0];
    if (!base) continue;
    if (fs.existsSync(path.join(ROOT, base))) { kept.push(base); continue; }
    // Try the common homes before dropping.
    const alt = ["src", "docs", "scripts", "bench", "tests", "evolve", "evolve/queue", "evolve/traces", "evolve/reports", ".github/workflows"]
      .map((d) => `${d}/${base}`)
      .find((c) => fs.existsSync(path.join(ROOT, c)));
    if (alt) kept.push(alt);
    else dropped.push(base);
  }
  return { kept: [...new Set(kept)], dropped: [...new Set(dropped)] };
}

/** Schema assertion. A goal that does not conform breaks the subagent packet. */
function validateGoal(g) {
  const errs = [];
  for (const f of GOAL_FIELDS) {
    if (!(f in g)) errs.push(`missing field "${f}"`);
  }
  if (g.id !== undefined && !/^[A-Za-z][\w.-]*$/.test(String(g.id))) errs.push(`bad id ${JSON.stringify(g.id)}`);
  if (g.horizon !== undefined && !HORIZONS.includes(g.horizon)) errs.push(`bad horizon ${JSON.stringify(g.horizon)}`);
  if (g.priority !== undefined && !Number.isInteger(g.priority)) errs.push(`priority must be an integer, got ${JSON.stringify(g.priority)}`);
  for (const f of ["title", "thesis", "acceptance"]) {
    if (g[f] !== undefined && (typeof g[f] !== "string" || !g[f].trim())) errs.push(`${f} must be a non-empty string`);
  }
  for (const f of ["read_only", "touch", "do_not_touch", "patterns"]) {
    if (g[f] !== undefined && !Array.isArray(g[f])) errs.push(`${f} must be an array`);
  }
  if (typeof g.thesis === "string" && !/\d/.test(g.thesis)) errs.push("thesis cites no number (not measurable)");
  return errs;
}

function cmdPlan(args = []) {
  const opts = { noGate: args.includes("--no-gate"), noBench: args.includes("--no-bench"), dry: args.includes("--dry") };

  // Measure. The gate vector is read from the last `gate` run; if there is
  // none (fresh clone) run it, so `plan` is never reasoning from a fiction.
  if (!opts.noGate && !state.gateVector) {
    console.error("[plan] no recorded gate vector; running the gate first...");
    const r = runGate();
    state.prevGateVector = null;
    state.gateVector = { cycle: state.cycle, score: Number(r.score.toFixed(4)), parts: r.parts, saturated: r.saturated };
    writeJson(STATE, state);
  }

  // Known ids across BOTH queues, open or closed. Re-planning must not
  // resurrect a closed goal or duplicate a live one.
  const known = new Set([...queue.goals, ...generated.goals].map((g) => g.id));

  const census = [];
  const candidates = [];
  for (const [source, probe] of PROBES) {
    let found = [];
    try { found = probe(opts) || []; } catch (e) { found = []; census.push({ source, error: String(e.message).slice(0, 120) }); continue; }
    census.push({ source, found: found.length });
    for (const c of found) candidates.push({ ...c, source });
  }

  // Rank by the MEASURED deficit, descending. Not by source order, not by
  // hand-assigned rank: `deficit` is a number each probe computed.
  candidates.sort((a, b) => (b.deficit - a.deficit) || (a.priority - b.priority) || String(a.key).localeCompare(String(b.key)));

  const skipped = [];
  const emitted = [];
  for (const c of candidates) {
    if (emitted.length >= PLAN_SIZE) { skipped.push({ ...c, why: `rank ${candidates.indexOf(c) + 1} > PLAN_SIZE ${PLAN_SIZE}` }); continue; }
    if (known.has(c.key)) { skipped.push({ key: c.key, title: c.title, why: "already in the queue (open or closed) — plan is idempotent" }); continue; }

    const { kept, dropped } = resolveReadOnly(c.read_only);
    if (!kept.length) { skipped.push({ key: c.key, title: c.title, why: `all ${c.read_only.length} read_only paths did not resolve` }); continue; }

    // Deterministic, queue-scoped id. P-numbered by rank so the file reads
    // in priority order; the key is kept for idempotency across re-plans.
    const id = `P${String(emitted.length + 1).padStart(2, "0")}`;
    const goal = {
      id,
      horizon: c.horizon,
      priority: c.priority,
      title: c.title,
      thesis: c.thesis,
      acceptance: c.acceptance,
      read_only: kept,
      touch: c.touch,
      do_not_touch: ["data/coreset/**", "node_modules/**", "evolve/queue/goals.json", "evolve/traces/**"],
      patterns: c.patterns,
      status: "open",
      origin: `deficit-analysis:${c.source}`,
      deficit: Number(c.deficit.toFixed(6)),
      measurement: c.key,
    };
    const errs = validateGoal(goal);
    if (errs.length) { skipped.push({ key: c.key, title: c.title, why: `SCHEMA: ${errs.join("; ")}` }); continue; }

    goal.read_onlyDropped = dropped;   // recorded, not silently swallowed
    emitted.push(goal);
    known.add(id);
  }

  const exhausted = candidates.length < PLAN_SIZE;

  // The authored queue is a DEFICIT SOURCE too, and it must not be dropped.
  // `plan` never rewrites goals.json, so the open authored goals survive by
  // construction; what the planner owes the reader is the accounting: which
  // authored goals are still open, that they are drained first, and which
  // generated goals would duplicate one.
  const openAuthored = openGoals(queue.goals);
  const authoredCarry = {
    openIds: openAuthored.map((g) => g.id),
    openCount: openAuthored.length,
    closedCount: queue.goals.filter((g) => g.status === "closed").length,
    drainedBefore: "generated",
    rule: "`next` drains evolve/queue/goals.json (authored) before evolve/queue/generated.json, so a re-plan can never starve the authored backlog. `plan` writes generated.json ONLY; goals.json is read-only to it.",
  };

  const meta = {
    planner: "evolve/evolve.mjs plan",
    cycle: state.cycle,
    generatedAt: new Date().toISOString(),
    gateVector: state.gateVector || null,
    R_best: state.R_best,
    PLAN_SIZE,
    candidatesMeasured: candidates.length,
    emitted: emitted.length,
    deficitPoolExhausted: exhausted,
    deficitPoolNote: exhausted
      ? `Only ${candidates.length} distinct measurable deficits were found; emitted ${emitted.length} of a nominal ${PLAN_SIZE}. NOT padded — an unmeasured filler goal costs a full agent call and returns nothing. Remaining debt is listed under "considered" with the reason it produced no goal.`
      : `${candidates.length} measurable deficits, top ${PLAN_SIZE} emitted, ${skipped.filter((s) => String(s.why).includes("PLAN_SIZE")).length} below the cut (see "considered").`,
    probeCensus: census,
    inconclusive: INCONCLUSIVE,
    authoredCarry,
    considered: skipped,
  };

  const doc = { _meta: meta, goals: emitted };

  // Report. This is the artifact a human reads to check the loop's reasoning.
  const line = (g) => `${g.id}  ${g.deficit.toFixed(3)}  [${g.origin.replace("deficit-analysis:", "")}]  ${g.title}`;
  console.log(JSON.stringify({
    plan: { emitted: emitted.length, of: PLAN_SIZE, candidatesMeasured: candidates.length, gateVector: state.gateVector, R_best: state.R_best },
    goals: emitted.map((g) => ({
      id: g.id, deficit: g.deficit, source: g.origin, title: g.title,
      measurement: g.measurement, thesisHead: g.thesis.slice(0, 200) + (g.thesis.length > 200 ? " ..." : ""),
    })),
    deficitPoolExhausted: exhausted,
    probeCensus: census,
    inconclusive: INCONCLUSIVE,
    authoredCarry,
    consideredButNotEmitted: skipped,
  }, null, 2));

  if (!opts.dry) {
    // WRITE ONLY generated.json. goals.json is never touched by `plan`, so a
    // re-plan cannot destroy a closed goal, its verdict, or its trace.
    const prevClosed = generated.goals.filter((g) => g.status && g.status !== "open");
    if (prevClosed.length) {
      // Preserve history: a closed generated goal keeps its record even if
      // the deficit that produced it is no longer present.
      for (const g of prevClosed) {
        if (!emitted.some((e) => e.id === g.id)) emitted.push(g);
      }
    }
    // CONTRACT SELF-ENFORCEMENT. `plan` must never emit a goal the loop's own
    // test suite rejects, or the loop ships a queue that is red by
    // construction. Which probes fire depends on the live gate vector, so a
    // thesis could satisfy the contract on one run and not the next — an
    // intermittent red the planner cannot explain. Every thesis is therefore
    // checked here against the same two rules tests/test_evolve_planner.js
    // enforces (cites a number; names the measurement), and a goal that fails
    // is DROPPED rather than padded into compliance: a goal we cannot evidence
    // is exactly the unmeasurable filler this command exists to avoid.
    const EVIDENCE_VOCAB = /(measured|grep|node |exit|lines|LOC|MB|files|ids|components|saturated|coverage=|bloat=|deficit)/i;
    const keepable = [];
    for (const g of emitted) {
      const t = String(g.thesis || "");
      const hasNumber = /\d/.test(t);
      const hasEvidence = EVIDENCE_VOCAB.test(t);
      if (hasNumber && hasEvidence) { keepable.push(g); continue; }
      // A closed goal keeps its record regardless — its verdict already stands.
      if (g.status && g.status !== "open") { keepable.push(g); continue; }
      skipped.push({
        id: g.id, deficit: 0,
        reason: `thesis did not satisfy the evidence contract ` +
                `(citesNumber=${hasNumber}, namesMeasurement=${hasEvidence}) — dropped rather than padded`,
      });
    }
    emitted.length = 0;
    emitted.push(...keepable);
    meta.droppedByContract = skipped.filter((s) => String(s.reason).includes("evidence contract")).length;
    writeJson(GEN, { _meta: meta, goals: emitted });
    console.error(`\n[plan] wrote ${emitted.length} goals -> ${path.relative(ROOT, GEN)} (goals.json untouched)`);
  } else {
    console.error(`\n[plan] --dry: nothing written`);
  }
}

const [, , cmd, ...rest] = process.argv;
switch (cmd) {
  case "next": cmdNext(); break;
  case "gate": cmdGate(); break;
  case "done": cmdDone(rest[0], rest[1], rest[2]); break;
  case "status": cmdStatus(); break;
  case "plan": cmdPlan(rest); break;
  default:
    console.log(fs.readFileSync(fileURLToPath(import.meta.url), "utf-8").split("*/")[0].replace(/^\/\*\*?/, ""));
}
