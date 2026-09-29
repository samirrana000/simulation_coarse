#!/usr/bin/env node
/**
 * evolve.mjs — the evolution engine.
 *
 * One command = one evolution step over ONE goal. This is the whole point:
 * 200 evolutions must not cost 200x15 agent calls. Each step:
 *
 *   1. WORK   — pick the highest-value unblocked goal from evolve/queue/goals.json
 *   2. THINK  — attach the wiki pattern ids + prior verdicts that bear on it
 *   3. ACT    — run the gate, record the delta in R (gate score), write trace
 *   4. LEARN  — append to the wiki skill-impact log; enqueue the NEXT 15 goals
 *
 * Usage:
 *   node evolve/evolve.mjs next            # print the single next goal (agent payload)
 *   node evolve/evolve.mjs gate            # run gate, print score, update R_best
 *   node evolve/evolve.mjs done <id> <verdict> [note]   # close a goal
 *   node evolve/evolve.mjs plan            # regenerate the next 15 goals
 *   node evolve/evolve.mjs status          # queue health
 *
 * Contract: this file NEVER edits src/. It only reads, measures, and routes.
 */

import fs from "node:fs";
import path from "node:path";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..");
const Q = path.join(HERE, "queue", "goals.json");
const STATE = path.join(HERE, "queue", "state.json");
const TRACES = path.join(HERE, "traces");
const IMPACT = path.join(ROOT, ".wikiskill", "wiki", "evolution", "skill-impact.md");

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


function runGate() {
  const r = { score: 0, parts: {}, notes: [] };

  // 1. syntax
  const srcFiles = [];
  (function walk(d) {
    for (const e of fs.readdirSync(d, { withFileTypes: true })) {
      const p = path.join(d, e.name);
      if (e.isDirectory()) walk(p);
      else if (e.name.endsWith(".js")) srcFiles.push(p);
    }
  })(path.join(ROOT, "src"));
  let syntaxOk = 0;
  const syntaxBad = [];
  for (const f of srcFiles) {
    try { execFileSync(process.execPath, ["--check", f], { stdio: "pipe" }); syntaxOk++; }
    catch (e) { syntaxBad.push(path.relative(ROOT, f)); }
  }
  r.parts.syntax = srcFiles.length ? syntaxOk / srcFiles.length : 0;
  r.syntaxNumbers = { ok: syntaxOk, total: srcFiles.length };
  if (syntaxBad.length) r.notes.push(`syntax FAIL: ${syntaxBad.join(", ")}`);

  // 2. tests — ratio against the recorded high-water mark, so adding real
  //    coverage keeps moving the needle instead of pinning at a legacy floor.
  try {
    const out = execFileSync(process.execPath, [path.join(ROOT, "tests", "test_all.js")], {
      encoding: "utf-8", stdio: ["ignore", "pipe", "pipe"], timeout: 300000,
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
    r.notes.push(`tests crashed: ${String(e.message).slice(0, 120)}`);
  }

  // 2b. coverage / anti-rot: every tests/test_*.js must be wired or in tests/manual/.
  r.parts.coverage = coverageScore();
  r.coverageNumbers = countTests();


  // 3. DOM contract
  try {
    const uiSrc = fs.readFileSync(path.join(ROOT, "src", "ui.js"), "utf-8");
    const html = fs.readFileSync(path.join(ROOT, "index.html"), "utf-8");
    const ids = [...uiSrc.matchAll(/\$\("([^"]+)"\)/g)].map((m) => m[1]);
    const missing = ids.filter((id) => !new RegExp(`id="${id}"`).test(html));
    r.parts.dom = ids.length ? (ids.length - missing.length) / ids.length : 1;
    r.domNumbers = { ids: ids.length, missing: missing.length };
  } catch { r.parts.dom = 0; }

  // 4. bloat — measured, not vibes. git objects + tracked bytes + src LOC.
  const tracked = readTrackedStats();
  r.parts.bloat = bloatScore(tracked);
  r.parts.science = scienceScore();
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
  const s = { files: 0, bytes: 0, srcLoc: 0, bigFiles: [] };
  let list = [];
  try { list = execFileSync("git", ["ls-files"], { cwd: ROOT, encoding: "utf-8" }).split("\n").filter(Boolean); }
  catch { return s; }
  s.files = list.length;
  for (const f of list) {
    try {
      const st = fs.statSync(path.join(ROOT, f));
      s.bytes += st.size;
      if (st.size > 1_000_000) s.bigFiles.push([f, Math.round(st.size / 1024) + "KB"]);
      if (f.startsWith("src/") && f.endsWith(".js")) {
        s.srcLoc += fs.readFileSync(path.join(ROOT, f), "utf-8").split("\n").length;
      }
    } catch {}
  }
  s.bigFiles.sort((a, b) => parseInt(b[1]) - parseInt(a[1]));
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

/** Bloat score in [0,1]. 1 = lean. Penalises >1MB tracked files hard. */function bloatScore(s) {
  const filePenalty = 1 / (1 + s.bigFiles.length / 20);
  const sizeScore = 1 / (1 + s.bytes / (500 * 1024 * 1024)); // 500MB -> 0.5
  const locScore = 1 / (1 + Math.max(0, s.srcLoc - 20000) / 20000);
  return 0.5 * filePenalty + 0.25 * sizeScore + 0.25 * locScore;
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
 * Goal queue
 * ------------------------------------------------------------------ */

const queue = readJson(Q, { goals: [] });

const HORIZON_ORDER = { short: 0, mid: 1, long: 2 };

function nextGoal() {
  // tolerate goals authored without an explicit status (default = open)
  for (const g of queue.goals) if (!g.status) g.status = "open";
  const open = queue.goals.filter((g) => g.status === "open");
  open.sort((a, b) => {
    if (a.horizon !== b.horizon) return HORIZON_ORDER[a.horizon] - HORIZON_ORDER[b.horizon];
    if (a.priority !== b.priority) return a.priority - b.priority;
    return a.cycle - b.cycle;
  });
  return open[0] || null;
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
  const r = runGate();
  console.log(JSON.stringify({
    score: Number(r.score.toFixed(4)),
    R_best: state.R_best,
    verdict: r.score >= state.R_best ? "OPEN (accept)" : "CLOSED (roll back)",
    parts: Object.fromEntries(Object.entries(r.parts).map(([k, v]) => [k, Number(v.toFixed(4))])),
    weights: WEIGHTS,
    saturated: r.saturated,
    testNumbers: r.testNumbers,
    coverageNumbers: r.coverageNumbers,
    domNumbers: r.domNumbers,
    tracked: r.tracked && { files: r.tracked.files, MB: +(r.tracked.bytes / 1048576).toFixed(1), srcLoc: r.tracked.srcLoc, bigFiles: r.tracked.bigFiles.slice(0, 8) },
    notes: r.notes,
  }, null, 2));
  // Compare at 4dp so float noise does not read as a regression.
  const r4 = Number(r.score.toFixed(4));
  const hwMoved = (state.testHighWater || 0) < 352 + 1 && r.testNumbers && r.testNumbers.passed > 352;
  if (r4 > state.R_best) {
    state.R_best = r4;
    writeJson(STATE, state);
    console.error(`\n[R_best updated -> ${state.R_best}]`);
  } else if (hwMoved) {
    // Persist the coverage high-water even when the composite did not move.
    writeJson(STATE, state);
  }
  if (state.rBaselineNote) console.error(`[NOTE] ${state.rBaselineNote}`);
}

function cmdDone(id, verdict, note = "") {
  const g = queue.goals.find((x) => x.id === id);
  if (!g) { console.error(`no such goal: ${id}`); process.exit(1); }
  g.status = verdict === "ACCEPTED" ? "closed" : "rolled-back";
  g.verdict = verdict;
  g.note = note;
  g.closedAtCycle = state.cycle;
  writeJson(Q, queue);
  if (verdict !== "ACCEPTED") state.rejected.push(`${id}: ${note}`.slice(0, 200));
  state.cycle++;
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
  const byH = {};
  for (const g of queue.goals) {
    byH[g.horizon] ??= { open: 0, closed: 0, "rolled-back": 0 };
    byH[g.horizon][g.status === "closed" ? "closed" : g.status === "rolled-back" ? "rolled-back" : "open"]++;
  }
  const next = nextGoal();
  console.log(JSON.stringify({
    cycle: state.cycle, R_best: state.R_best,
    goals: { total: queue.goals.length, byHorizon: byH },
    next: next && `${next.id} — ${next.title}`,
    rejectedRemembered: state.rejected.length,
  }, null, 2));
}

const [, , cmd, ...rest] = process.argv;
switch (cmd) {
  case "next": cmdNext(); break;
  case "gate": cmdGate(); break;
  case "done": cmdDone(rest[0], rest[1], rest[2]); break;
  case "status": cmdStatus(); break;
  default:
    console.log(fs.readFileSync(fileURLToPath(import.meta.url), "utf-8").split("*/")[0].replace(/^\/\*\*?/, ""));
}
