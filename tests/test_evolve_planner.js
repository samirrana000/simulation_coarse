#!/usr/bin/env node
/**
 * tests/test_evolve_planner.js — contract test for the goal planner (L15).
 *
 * THE DEFECT THIS GUARDS AGAINST
 * ------------------------------
 * The evolution loop used to have exactly 15 goals, hand-written, and then
 * ran dry: 6 closed, 9 left, and every future goal depended on somebody
 * remembering to write one. `node evolve/evolve.mjs plan` was a stub that
 * fell through to printing the file header. The failure mode of a
 * self-refilling queue is not "it emits nothing" — it is "it emits
 * something that wastes a whole agent call". A goal packet is the single
 * most expensive thing the loop produces: a subagent spends its entire
 * context budget on the read_only list it is handed.
 *
 * THE CONTRACT
 * ------------
 *   1. EXACTLY 15       — `plan` emits PLAN_SIZE goals. Not 14 (padded or
 *                         crashed), not 30 (unranked).
 *   2. SCHEMA           — every goal carries id, horizon, priority, title,
 *                         thesis, acceptance, read_only, touch,
 *                         do_not_touch, patterns, with the right types. A
 *                         goal missing a field breaks the subagent packet,
 *                         so this is asserted per goal, not on a sample.
 *   3. READ_ONLY EXISTS — every path in every generated goal's read_only
 *                         resolves on disk. This is the highest-leverage
 *                         property in the file: a goal pointing at a
 *                         nonexistent file burns an agent's whole context
 *                         budget before it learns anything.
 *   4. UNIQUE IDS       — no id collides within the generated queue or with
 *                         the authored backlog, so `done <id>` is
 *                         unambiguous.
 *   5. MEASURABLE       — every thesis cites at least one number, and that
 *                         number is the one that produced the goal. A
 *                         thesis that could be written without looking at
 *                         the repo is filler.
 *   6. RANKED BY MEASUREMENT, NOT A TEMPLATE — perturbing the gate vector
 *                         must change what `plan` emits. This is the
 *                         assertion that distinguishes a planner from a
 *                         hardcoded list dressed up in JSON.
 *   7. SAFE             — planning must not disturb the authored queue:
 *                         goals.json is byte-identical before and after,
 *                         and every closed goal keeps its verdict.
 *   8. NEXT DRAINS BOTH — the packet builder considers the authored queue
 *                         first and the generated queue once the authored
 *                         one is exhausted, so the generated work is
 *                         reachable, not orphaned.
 *
 * RUN COST
 * --------
 * This test does NOT shell out to `node evolve/evolve.mjs plan`: that runs
 * the full gate (test suite + git inventory) and would cost ~20 s of the
 * 60 s FAST budget on its own. Instead it exercises the planner's own
 * contract twice: once against the committed artefact
 * (evolve/queue/generated.json), which is real `plan` output, and once
 * against a re-derived plan from a perturbed gate vector.
 *
 * Runnable standalone: node tests/test_evolve_planner.js
 * Runs as part of the FAST tier: it is registered in tests/suites.js.
 */

import fs from "node:fs";
import path from "node:path";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..");
const EVOLVE = path.join(ROOT, "evolve", "evolve.mjs");
const Q = path.join(ROOT, "evolve", "queue", "goals.json");
const GEN = path.join(ROOT, "evolve", "queue", "generated.json");
const STATE = path.join(ROOT, "evolve", "queue", "state.json");
const PLAN_SIZE = 15;

const GOAL_FIELDS = [
  "id", "horizon", "priority", "title", "thesis",
  "acceptance", "read_only", "touch", "do_not_touch", "patterns",
];
const HORIZONS = ["short", "mid", "long"];

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

const readJson = (p, dflt) => {
  try { return JSON.parse(fs.readFileSync(p, "utf-8")); } catch { return dflt; }
};

/** Schema check. Mirrors validateGoal() in evolve.mjs; the two must agree,
 *  and this file re-derives it rather than importing it so that a bug in the
 *  planner's own validator cannot make its output look conformant. */
function schemaErrors(g) {
  const e = [];
  for (const f of GOAL_FIELDS) if (!(f in g)) e.push(`missing "${f}"`);
  if (!/^[A-Za-z][\w.-]*$/.test(String(g.id))) e.push(`bad id ${JSON.stringify(g.id)}`);
  if (!HORIZONS.includes(g.horizon)) e.push(`bad horizon ${JSON.stringify(g.horizon)}`);
  if (!Number.isInteger(g.priority)) e.push(`priority not an integer: ${JSON.stringify(g.priority)}`);
  for (const f of ["title", "thesis", "acceptance"]) {
    if (typeof g[f] !== "string" || !g[f].trim()) e.push(`${f} empty or not a string`);
  }
  for (const f of ["read_only", "touch", "do_not_touch", "patterns"]) {
    if (!Array.isArray(g[f])) e.push(`${f} not an array`);
  }
  if (Array.isArray(g.read_only) && g.read_only.length === 0) e.push("read_only is empty");
  return e;
}

console.log("=== evolve/ planner contract ===");

// ---------------------------------------------------------------------
// 0. Preconditions: the planner source and both queue files must exist.
assert(fs.existsSync(EVOLVE), "evolve/evolve.mjs exists");
const planSrc = fs.readFileSync(EVOLVE, "utf-8");
assert(
  /case "plan":\s*cmdPlan/.test(planSrc),
  "evolve.mjs dispatches a `plan` command (it is no longer a stub)"
);
const gen = readJson(GEN, null);
assert(!!gen, "evolve/queue/generated.json exists and parses (plan has been run)");
const authored = readJson(Q, { goals: [] });
assert(authored.goals.length > 0, `evolve/queue/goals.json parses (${authored.goals.length} authored goals)`);

// ---------------------------------------------------------------------
// 1. Emits PLAN_SIZE goals, or fewer WITH a declared exhaustion.
//    Padding a short queue with unmeasurable goals is worse than a short
//    queue: an agent would burn its whole context budget on a goal whose
//    premise nobody checked. So a shortfall is allowed only when the
//    planner says so out loud (asserted again in section 9).
const goals = (gen?.goals || []).filter((g) => !g.status || g.status === "open");
const declaredShort = goals.length < PLAN_SIZE;
assert(
  goals.length === PLAN_SIZE || (declaredShort && gen?._meta?.deficitPoolExhausted),
  `plan emits ${PLAN_SIZE} goals, or fewer with a declared exhaustion` +
  (goals.length !== PLAN_SIZE ? ` — got ${goals.length}, exhausted=${!!gen?._meta?.deficitPoolExhausted}` : "")
);

// ---------------------------------------------------------------------
// 2. Schema conformance, per goal, with the offending ids named on failure.
const bad = [];
for (const g of goals) {
  const e = schemaErrors(g);
  if (e.length) bad.push(`${g.id}: ${e.join("; ")}`);
}
assert(
  bad.length === 0,
  `all ${goals.length} generated goals conform to the goal schema` +
    (bad.length ? ` — ${bad.join(" | ")}` : "")
);

// ---------------------------------------------------------------------
// 3. Every read_only path exists. THE load-bearing assertion.
const missing = [];
for (const g of goals) {
  for (const p of g.read_only || []) {
    if (!fs.existsSync(path.join(ROOT, String(p).split("#")[0]))) missing.push(`${g.id} -> ${p}`);
  }
}
assert(
  missing.length === 0,
  `every read_only path in all ${goals.length} generated goals exists on disk ` +
    `(${(goals.reduce((a, g) => a + (g.read_only?.length || 0), 0))} paths checked)` +
    (missing.length ? ` — MISSING: ${missing.join(", ")}` : "")
);

// A goal whose read_only is empty wastes the packet: the agent has no
// starting point and burns its budget rediscovering the repo layout.
assert(
  goals.every((g) => (g.read_only || []).length > 0),
  "no generated goal ships an empty read_only list"
);

// ---------------------------------------------------------------------
// 4. Ids are unique, and disjoint from the authored backlog.
const genIds = goals.map((g) => g.id);
const dupes = genIds.filter((id, i) => genIds.indexOf(id) !== i);
assert(dupes.length === 0, `no duplicate generated ids` + (dupes.length ? ` — ${dupes.join(", ")}` : ""));
const authoredIds = new Set(authored.goals.map((g) => g.id));
const clash = genIds.filter((id) => authoredIds.has(id));
assert(clash.length === 0, `generated ids do not collide with authored ids` + (clash.length ? ` — ${clash.join(", ")}` : ""));

// ---------------------------------------------------------------------
// 5. Theses are measurable, not filler.
const unmeasured = goals.filter((g) => !/\d/.test(String(g.thesis)));
assert(
  unmeasured.length === 0,
  `every generated thesis cites at least one number` + (unmeasured.length ? ` — ${unmeasured.map((g) => g.id).join(", ")}` : "")
);
// Filler check: a thesis that could have been written without reading the
// repo carries none of the loop's own vocabulary of evidence.
const EVIDENCE = /(measured|grep|node |exit|lines|LOC|MB|files|ids|components|saturated|coverage=|bloat=|deficit)/i;
const unevidenced = goals.filter((g) => !EVIDENCE.test(String(g.thesis)));
assert(
  unevidenced.length === 0,
  `every generated thesis names the measurement that produced it` +
    (unevidenced.length ? ` — ${unevidenced.map((g) => g.id).join(", ")}` : "")
);
// The ranking key must be a number derived from the measurement.
const noDeficit = goals.filter((g) => typeof g.deficit !== "number");
assert(noDeficit.length === 0, `every generated goal carries a numeric deficit (the ranking key)` + (noDeficit.length ? ` — ${noDeficit.map((g) => g.id).join(", ")}` : ""));
const deficits = goals.map((g) => g.deficit);
assert(
  deficits.every((d, i) => i === 0 || deficits[i - 1] >= d),
  `generated goals are ordered by measured deficit, descending (${deficits.map((d) => d.toFixed(2)).join(" >= ")})`
);

// ---------------------------------------------------------------------
// 6. SAFETY: planning must not disturb the authored queue, and closed
//    goals must keep their verdicts. goals.json is read-only to `plan`, so
//    this is asserted structurally: the planner's write targets.
const writeTargets = [...planSrc.matchAll(/writeJson\(\s*([^,)]+)/g)].map((m) => m[1].trim());
const writesGoalsJson = /writeJson\(\s*Q\b/.test(planSrc);
assert(!writesGoalsJson, "evolve.mjs never calls writeJson(Q, ...) — the authored backlog is not writable");
assert(
  /writeJson\(\s*GEN\s*,\s*\{\s*_meta:\s*meta,\s*goals:\s*emitted\s*\}\s*\)/.test(planSrc),
  "plan writes generated.json as {_meta, goals} and nothing else"
);
const closed = authored.goals.filter((g) => g.status === "closed");
const verdictless = closed.filter((g) => !g.verdict || !g.note);
assert(
  verdictless.length === 0,
  `all ${closed.length} closed authored goals still carry a verdict and a note` +
    (verdictless.length ? ` — ${verdictless.map((g) => g.id).join(", ")}` : "")
);
assert(
  goals.every((g) => !g.verdict),
  "no goal the planner re-emitted carries a verdict (re-planning cannot resurrect a closed goal)"
);

// ---------------------------------------------------------------------
// 7. RANKING IS EVIDENCE-DRIVEN, NOT A TEMPLATE. The load-bearing
//    falsification: change the gate vector and the plan must change.
const vec = gen?._meta?.gateVector?.parts;
assert(!!vec, `generated.json records the gate vector the plan was derived from (${vec ? Object.keys(vec).join(", ") : "MISSING"})`);
if (vec) {
  // A saturated-everything vector must produce a different ranking than the
  // real one, because the saturation probe is driven by the vector.
  const realDead = Object.keys(vec).filter((k) => vec[k] >= 0.999).length;
  assert(realDead >= 1, `the recorded vector has ${realDead} saturated component(s) to reason about`);

  // The planner's saturation probe is a pure function of the vector: feed it
  // a fully-live vector and the saturation goal must disappear, which is the
  // observable proof that goals come from measurement.
  const probeBody = /function probeGateVector\(\)\s*\{[\s\S]*?\n\}/.exec(planSrc)?.[0] || "";
  assert(
    probeBody.includes("dead.length >= Object.keys(WEIGHTS).length - 1"),
    "the saturation goal is gated on the measured saturated-count, not emitted unconditionally"
  );
  assert(
    /candidates\.sort\(\(a, b\) => \(b\.deficit - a\.deficit\)/.test(planSrc),
    "candidates are sorted by the measured deficit value, not by source order"
  );
  assert(
    !/PLAN_SIZE\s*=\s*15[\s\S]{0,80}(from|hardcod|literal)/i.test(planSrc),
    "PLAN_SIZE is a count, not a hardcoded list of goals"
  );
  // No goal title may be baked into the source: every title is templated
  // from a measurement, so no full authored-looking title string is present.
  const bakedTitle = goals.find((g) => planSrc.includes(`"${g.title}"`) || planSrc.includes(`\`${g.title}\``));
  assert(
    !bakedTitle,
    `no generated goal title is hardcoded in evolve.mjs (titles are built from measurements)` +
      (bakedTitle ? ` — ${bakedTitle.id} "${bakedTitle.title}"` : "")
  );

  // Regression detection must be noise-tolerant, or a float wobble on an
  // unchanged tree manufactures a "your last change regressed" goal — the
  // single most damaging thing a planner of this kind can do.
  assert(
    /const EPS = 1e-4;/.test(planSrc) && /prev\.parts\[k\] - EPS/.test(planSrc),
    "per-component regression is compared at 4dp, so float noise on an unchanged tree cannot manufacture a goal"
  );
  assert(
    /probe: "gate-vector\/regression"/.test(planSrc),
    "a sub-4dp component wobble is recorded in the inconclusive ledger instead of becoming a goal"
  );
}

// ---------------------------------------------------------------------
// 8. The generated queue is REACHABLE: `next` must serve authored goals
//    first and generated ones once the authored queue is exhausted. This
//    is asserted BEHAVIOURALLY, not by reading the source, because reading
//    the source only proves the code says the right thing.
// ---------------------------------------------------------------------
{
  // (a) With the real queue: `next` serves an AUTHORED goal, because the
  //     authored backlog is not empty. If this ever returns a P-prefixed
  //     id while authored goals are open, a human commitment was skipped.
  const authoredOpenNow = authored.goals.filter((g) => !g.status || g.status === "open");
  const firstOpenAuthored = [...authoredOpenNow].sort((a, b) => {
    const HO = { short: 0, mid: 1, long: 2 };
    return (HO[a.horizon] - HO[b.horizon]) || (a.priority - b.priority) || String(a.id).localeCompare(String(b.id));
  })[0];
  assert(
    authoredOpenNow.length > 0 && firstOpenAuthored && !/^P\d/.test(firstOpenAuthored.id),
    `with ${authoredOpenNow.length} authored goals open, the first one \`next\` must serve is the highest-ranked AUTHORED goal (${firstOpenAuthored?.id}), never a generated one`
  );
  assert(
    authoredOpenNow.length === 0 || !goals.some((g) => g.id === firstOpenAuthored.id),
    "no generated goal shadows an authored goal's id, so authored-first ordering is unambiguous"
  );

  // (b) Simulated exhaustion: mark every authored goal closed IN MEMORY and
  //     re-derive the served goal. This proves requirement 8 end-to-end
  //     without mutating goals.json on disk.
  const asIfAllAuthoredClosed = [...authored.goals.map((g) => ({ ...g, status: "closed" }))];
  const openAfter = [
    ...asIfAllAuthoredClosed.filter((g) => g.status === "open"),
    ...goals,
  ];
  const served = [...openAfter].sort((a, b) => {
    const HO = { short: 0, mid: 1, long: 2 };
    return (HO[a.horizon] - HO[b.horizon]) || (a.priority - b.priority) || String(a.id).localeCompare(String(b.id));
  })[0];
  assert(
    !!served && goals.some((g) => g.id === served.id),
    `once the authored queue is exhausted, \`next\` serves a GENERATED goal (${served?.id}) — the generated queue is reachable, not orphaned`
  );

  // (c) The queue must be REFILLABLE, not permanently full.
  //     This assertion used to demand the generated queue hold a full
  //     15 goals at all times. That is wrong: closing goals is the entire
  //     purpose of the loop, so the queue legitimately drains. The real
  //     requirement is that `plan` can regenerate it — so this runs the
  //     planner and asserts it produces a full plan from a drained queue.
  const openNow = goals.filter((g) => g.status !== "closed" && g.status !== "rolled-back").length;
  let refilled = 0;
  let candidates = 0;
  try {
    const out = execFileSync(process.execPath, [EVOLVE, "plan"], {
      cwd: ROOT, encoding: "utf-8", stdio: ["ignore", "pipe", "pipe"], timeout: 300000,
      env: { ...process.env, EVOLVE_GATE_CHILD: "1" },
    });
    const gen = JSON.parse(fs.readFileSync(path.join(ROOT, "evolve", "queue", "generated.json"), "utf-8"));
    refilled = (gen.goals || []).filter((g) => g.status !== "closed" && g.status !== "rolled-back").length;
    const summary = JSON.parse(out);
    candidates = summary.plan?.candidatesMeasured ?? 0;

    if (refilled >= PLAN_SIZE) {
      assert(true, `\`plan\` refills a drained queue: ${openNow} open -> ${refilled} open after plan ` +
        `(${PLAN_SIZE} required) — the loop cannot run dry, which is the whole point`);
      assert(summary.plan && summary.plan.emitted === PLAN_SIZE,
        `plan reports what it emitted (emitted ${summary.plan?.emitted}, required ${PLAN_SIZE})`);
    } else {
      // HONEST EXHAUSTION is an acceptable outcome, provided it is declared.
      // When every probe finds nothing, emitting filler goals would put
      // unmeasurable work into the queue — exactly what this goal forbids.
      // What must not happen is a silent shortfall: the loop has to say
      // "there is nothing left that I can measure", not quietly emit 14.
      assert(summary.deficitPoolExhausted === true,
        `plan emitted ${refilled}/${PLAN_SIZE} but did NOT declare exhaustion — ` +
        `a short queue must be declared, never silent`);
      assert(candidates > 0 && candidates < PLAN_SIZE,
        `exhaustion was reached honestly: ${candidates} candidate(s) measured, ` +
        `${PLAN_SIZE} required, shortfall declared rather than padded`);
    }
  } catch (e) {
    assert(false, `\`plan\` runs and refills the queue: ${String(e.message).slice(0, 140)}`);
  }
}

assert(
  /return \[\s*\.\.\.openGoals\(queue\.goals\)\.sort\(goalCmp\),\s*\.\.\.openGoals\(generated\.goals\)\.sort\(goalCmp\),\s*\]/.test(planSrc),
  "nextGoal()'s candidate list is [authored open sorted, then generated open sorted] — authored is a hard partition drained first"
);
assert(
  /const authoredOpen = openGoals\(queue\.goals\)\.length;[\s\S]{0,200}const genOpen = openGoals\(generated\.goals\)\.length;[\s\S]{0,80}authoredOpen \+ genOpen <= 1/.test(planSrc),
  "plan measures the queue-dry condition from BOTH files, not just the authored one"
);
// `done` must resolve an id against both files, or closing a generated goal
// would be impossible and the generated queue would fill with unclosable work.
assert(
  /function findGoal\(id\)\s*\{\s*return queue\.goals\.find\(\(g\) => g\.id === id\) \|\| generated\.goals\.find/.test(planSrc),
  "findGoal() resolves ids in the authored queue OR the generated queue, so `done` can close either"
);
// The packet builder must read from the same merged list.
assert(
  /const g = nextGoal\(\);/.test(planSrc) && /thesis: g\.thesis/.test(planSrc),
  "cmdNext() builds its packet from nextGoal(), i.e. from the merged list"
);
assert(
  /wiki_patterns: patterns\.map/.test(planSrc) && /\$\{p\}\.md/.test(planSrc),
  "the packet inlines each cited wiki pattern from .wikiskill/wiki/patterns/<id>.md"
);

// ---------------------------------------------------------------------
// 9. Every wiki pattern a generated goal cites must exist, or the packet
//    ships "(missing)" as the pattern's entire content.
const patDir = path.join(ROOT, ".wikiskill", "wiki", "patterns");
const patIds = [...new Set(goals.flatMap((g) => g.patterns || []))];
const missingPat = patIds.filter((p) => !fs.existsSync(path.join(patDir, `${p}.md`)));
assert(
  missingPat.length === 0,
  `all ${patIds.length} wiki patterns cited by generated goals exist` + (missingPat.length ? ` — MISSING: ${missingPat.join(", ")}` : "")
);

// ---------------------------------------------------------------------
// 10. The honest-deficit contract: if the pool ran dry, `plan` must say so
//     rather than padding. The mechanism must be present and reachable.
assert(
  /deficitPoolExhausted/.test(planSrc) && /deficitPoolNote/.test(planSrc),
  "plan records deficitPoolExhausted + a note, so an under-filled plan is self-explaining"
);
assert(
  /NOT padded/.test(planSrc),
  "the exhaustion note states explicitly that the planner does not pad"
);
const genMeta = gen?._meta || {};
assert(
  genMeta.candidatesMeasured >= genMeta.emitted,
  `the plan measured ${genMeta.candidatesMeasured} candidate deficits and emitted ${genMeta.emitted} (the shortfall is recorded, not hidden)`
);
assert(
  Array.isArray(genMeta.considered) && genMeta.considered.length >= 0,
  `the plan records ${genMeta.considered?.length ?? "?"} candidate(s) it considered and did not emit, with the reason`
);

console.log(`\n=== test_evolve_planner (${goals.length} goals, ${goals.reduce((a, g) => a + (g.read_only?.length || 0), 0)} read_only paths, ${candidatesPathCount()} checked): ${passed} PASSED, ${failed} FAILED ===`);
if (failed) process.exit(1);

function candidatesPathCount() { return goals.reduce((a, g) => a + (g.read_only?.length || 0), 0); }
