/**
 * tests/test_doc_citations.js — every `path:line` claim in every tracked
 * markdown file must resolve to a real, non-blank line.
 *
 * WHY THIS EXISTS
 * ---------------
 * A `src/heavy.js:415` in a doc is a falsifiable claim: the reader is told
 * "the thing I just described is at that exact line". It is the cheapest and
 * most numerous kind of claim this repo makes, and until now nothing checked
 * it. The evolution loop's planner measured "3 of 612 broken"; measurement
 * with the rule below finds 37 (1 file deleted, 36 pointers into blank
 * lines). Either way, a silently-wrong rate on a category the project treats
 * as evidence is the defect class ROADMAP section 1 exists to prevent.
 *
 * THE RULE (decide once, document it, then hold to it)
 * ---------------------------------------------------
 * 1. SCOPE: every `*.md` file tracked by git (`git ls-files '*.md'`). Not the
 *    working tree — an untracked scratch file must not be able to fail CI.
 * 2. PATTERN: a repo-relative path with a source extension, followed by `:`
 *    and a line number or an inclusive `a-b` range:
 *
 *        src/heavy.js:415
 *        src/placement.js:89-104
 *        `src/viewer.js:110`          (backticks are not part of the citation)
 *
 * 3. RESOLUTION — a four-step deterministic ladder, first hit wins:
 *      (a) the path as written, relative to the REPO ROOT;
 *      (b) the path as written, relative to the DOC'S OWN directory;
 *      (c) `src/` + the path as written (docs/ writes `physics/hbond.js`
 *          meaning `src/physics/hbond.js`);
 *      (d) a git-tracked file whose BASENAME equals the path's basename.
 *    Step (d) exists because the survey docs (BINDING_PHYSICS_R*) use
 *    `heavy.js:383-404` as shorthand throughout; it is allowed ONLY when
 *    exactly one tracked file has that basename, so an ambiguous pointer
 *    (four files are named README.md) FAILS as ambiguous rather than silently
 *    picking one. Ladder order is fixed so the result cannot depend on
 *    filesystem ordering.
 *    The resolved file must be git-tracked.
 * 4. THE LINE MUST EXIST and MUST NOT BE BLANK. A pointer past EOF is
 *    obviously broken. A pointer into a blank line is equally broken and has
 *    zero false positives: no answer to "where is this?" is ever whitespace.
 *    A range `a-b` requires a <= b and both endpoints to satisfy 3 and 4.
 * 5. `:??N` / `:?` is a deliberately-unknown line, not a literal citation, and
 *    is skipped by rule. It is rare enough that no doc currently uses it.
 * 6. Anything genuinely exempt goes in ALLOWLIST below, keyed on the exact
 *    parsed citation in a named doc, with a reason. An ALLOWLIST entry whose
 *    text no longer appears is reported as STALE and must be deleted, so the
 *    escape hatch cannot outlive its justification.
 *
 * WHY THE ALLOWLIST IS AN EXPLICIT LIST AND NOT A REGEX
 * ---------------------------------------------------
 * The tempting version is `if (/:(\d+)\?*[-–]?(\d+)?/` and hope the shape is
 * obvious. That is how a checker rots: a real broken citation whose line
 * happens to be followed by a question mark stops being checked and nobody
 * notices. An explicit `doc + exact-citation -> reason` list is auditable.
 *
 * WHAT THIS TEST DELIBERATELY DOES NOT DO
 * --------------------------------------
 * It does not verify that line 415 CONTAINS what the sentence claims. That is
 * semantic drift, not mechanical drift, and pretending otherwise would be the
 * exact over-claiming this test exists to remove. It proves the pointer
 * resolves to a real line; a human still owns the claim.
 */

import { execFileSync } from "child_process";
import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(__dirname, "..");

/**
 * Citations intentionally NOT checked, each with its reason.
 *
 * Every entry is keyed on the EXACT parsed citation (`path:spec`) in a named
 * document, so the exemption is countable and cannot silently widen: a new
 * broken citation in the same document is still checked, and an entry whose
 * text no longer appears is reported STALE and must be deleted.
 *
 * All of these live in DATED ARTIFACTS whose `file:line` values are a snapshot
 * of the tree on the day the artifact was written. Repointing them at today's
 * code would be the worse lie: the surrounding snippets are labelled
 * "verbatim" and the entire value of a dated audit/log is that it records what
 * was true then. Both files now carry a banner saying exactly that.
 */
export const ALLOWLIST = [
  // --- AUDIT_REPORT.md, dated 2026-08-29, header carries a STATUS banner ---
  { doc: "AUDIT_REPORT.md", cite: "src/settings-panel.js:73-79", why: "2026-08-29 line numbers of the buggy snapshot (section heading naming the defect site)" },
  { doc: "AUDIT_REPORT.md", cite: "src/settings-panel.js:73-85", why: "2026-08-29 line numbers inside a block labelled 'verbatim'" },
  { doc: "AUDIT_REPORT.md", cite: "src/settings-panel.js:78", why: "2026-08-29 line number in the P0 fix checklist" },
  { doc: "AUDIT_REPORT.md", cite: "viewer.js:226-234", why: "2026-08-29 line numbers of the dead-code _resize() guard" },
  { doc: "AUDIT_REPORT.md", cite: "src/viewer.js:226-234", why: "2026-08-29 line numbers, same snippet, path-prefixed form" },
  { doc: "AUDIT_REPORT.md", cite: "src/viewer.js:226", why: "2026-08-29 line number in the P1 fix checklist" },
  { doc: "AUDIT_REPORT.md", cite: "viewer.js:150-167", why: "2026-08-29 line numbers of the centroid/radius sizing bug" },
  // The file's own STATUS banner already declares ALL of its `file:line`
  // citations to be the 2026-08-29 snapshot and "deliberately NOT repointed at
  // today's code — repointing them would falsify the record". When src/main.js
  // was split into src/controllers/ (2026-10) its main.js pointers went past
  // EOF; the banner's stated policy is allowlisting, so these follow it.
  { doc: "AUDIT_REPORT.md", cite: "main.js:31", why: "2026-08-29 line number of the settings-panel import in the graph diagram" },
  { doc: "AUDIT_REPORT.md", cite: "main.js:620", why: "2026-08-29 line number of requestAnimationFrame(tick) in the browser-consequence list" },
  { doc: "AUDIT_REPORT.md", cite: "main.js:538-616", why: "2026-08-29 line range of the viewer-guarded tick body" },
  { doc: "AUDIT_REPORT.md", cite: "main.js:206-282", why: "2026-08-29 line range of the early-return analysis in §3.2" },
  { doc: "AUDIT_REPORT.md", cite: "main.js:206-234", why: "2026-08-29 line range of buildSystem in §3.2" },
  { doc: "AUDIT_REPORT.md", cite: "main.js:230-234", why: "2026-08-29 line range quoted in §3.2" },
  { doc: "AUDIT_REPORT.md", cite: "main.js:251", why: "2026-08-29 line number of the workerPool.initSystem call in §3.4" },
  { doc: "AUDIT_REPORT.md", cite: "src/main.js:230", why: "2026-08-29 line number in the §5 fix checklist", },

  // --- TRANSFORMATION_PLAN_100.md: a dated 2026-08-29 delivery log ----------
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "settings-panel.js:73", why: "2026-08-29 rationale pointer in the A02 row of a dated delivery log" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "viewer.js:226", why: "2026-08-29 rationale pointer in the same A02 row" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "pdb.js:29", why: "2026-08-29 pointer to the pre-REST-fetch import block" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "src/heavy.js:87", why: "2026-08-29 pointer to the warnings/countWarnings helper, before the header block grew" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "src/recorder.js:56", why: "2026-08-29 pointer to the getFrame() H78 placeholder" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "README.md:29", why: "2026-08-29 pointer to the README no-PME bullet, before the doc-index section landed" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "main.js:591", why: "2026-08-29 line number of the live-tracking classifyPose call in the E49 row" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "main.js:525", why: "2026-08-29 line number of integ.advance(steps,14) in the G68 row" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "main.js:602", why: "2026-08-29 line number of the HUD debounce guard in the H76 row" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "src/main.js:551", why: "2026-08-29 line number of lastHudUpdate in the H76 verification table" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "src/main.js:595", why: "2026-08-29 line number of the funnel nHills HUD row in the D33 verification table" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "src/main.js:665", why: "2026-08-29 line number of the HUD version prefix in the A01 verification table" },
  { doc: "TRANSFORMATION_PLAN_100.md", cite: "src/main.js:232", why: "2026-08-29 line number of console.error in buildSystem in the A03 verification table" },
];

let pass = 0;
let fail = 0;

function ok(name, detail) {
  pass++;
  console.log(`PASS ${name}${detail ? ` — ${detail}` : ""}`);
}

function bad(name, detail) {
  fail++;
  console.log(`FAIL ${name} — ${detail}`);
}

function gitLines(args) {
  return execFileSync("git", args, { cwd: ROOT, encoding: "utf-8" })
    .split("\n")
    .filter((s) => s.length > 0);
}

/** Split a file into lines; a trailing newline terminates the last line. */
function linesOf(abs) {
  const txt = fs.readFileSync(abs, "utf-8");
  if (txt.length === 0) return [];
  const l = txt.split("\n");
  if (l[l.length - 1] === "") l.pop();
  return l;
}

const trackedMd = gitLines(["ls-files", "*.md"]);
if (trackedMd.length === 0) {
  bad("citation scan is not vacuous", "git ls-files '*.md' returned nothing");
} else {
  ok("citation scan is not vacuous", `${trackedMd.length} tracked markdown files`);

  const allTracked = gitLines(["ls-files"]);
  const tracked = new Set(allTracked);
  const byBasename = new Map();
  for (const p of allTracked) {
    const b = path.posix.basename(p);
    if (!byBasename.has(b)) byBasename.set(b, []);
    byBasename.get(b).push(p);
  }

  // Extensions a citation may point at. Anything else (.csv, .pdf, .pdb, .ipynb)
  // is a data file or a build artifact; a line number into one is noise.
  const CITE =
    /`?((?:\.{0,2}\/)?[A-Za-z0-9_.\-\/]+\.(?:js|mjs|sh|json|html|css|py|yml|yaml|md|ts|glsl|wgsl)):(\?\d*|\d+(?:\s*[-–]\s*\d+)?)`?/g;

  /** The four-step ladder. Returns {rel, how} or {why}. */
  function resolve(target, doc) {
    const docDir = path.posix.dirname(doc);
    const steps = [
      path.posix.normalize(target),
      path.posix.normalize(path.posix.join(docDir, target)),
      path.posix.normalize(path.posix.join("src", target)),
    ];
    for (const rel of steps) {
      if (rel.startsWith("..")) continue;
      if (tracked.has(rel) && fs.existsSync(path.join(ROOT, rel))) return { rel, how: "as-written" };
    }
    const b = path.posix.basename(target);
    const hits = byBasename.get(b) || [];
    if (hits.length === 1) return { rel: hits[0], how: "unique-basename" };
    if (hits.length > 1) return { why: `ambiguous basename ${b} (${hits.length} candidates: ${hits.join(", ")})` };
    return { why: `no such file: ${target}` };
  }

  let total = 0;
  const brokenLine = [];
  const brokenTarget = [];
  const usedAllow = new Set();
  const ladderUse = new Map();

  for (const doc of trackedMd) {
    const abs = path.join(ROOT, doc);
    if (!fs.existsSync(abs)) {
      bad("doc exists", `${doc} is tracked but absent from the working tree`);
      continue;
    }
    const docLines = fs.readFileSync(abs, "utf-8").split("\n");
    for (let i = 0; i < docLines.length; i++) {
      CITE.lastIndex = 0;
      let m;
      while ((m = CITE.exec(docLines[i])) !== null) {
        const [, target, spec] = m;
        const key = `${target}:${spec}`;
        total++;
        const where = `${doc}:${i + 1}`;

        const exempt = ALLOWLIST.find((a) => a.doc === doc && a.cite === key);
        if (exempt) {
          usedAllow.add(exempt);
          continue;
        }
        if (spec.startsWith("?")) continue; // rule 5: deliberately unknown

        const r = resolve(target, doc);
        if (r.why) {
          brokenTarget.push({ where, key, why: r.why });
          continue;
        }
        ladderUse.set(r.how, (ladderUse.get(r.how) || 0) + 1);

        const L = linesOf(path.join(ROOT, r.rel));
        const nums = spec.split(/[-–]/).map((s) => parseInt(s.trim(), 10));
        let why = null;
        for (const k of nums) {
          if (!Number.isFinite(k) || k < 1) { why = `line ${k} is not >= 1`; break; }
          if (k > L.length) { why = `line ${k} is past EOF (${r.rel} has ${L.length} lines)`; break; }
          if (L[k - 1].trim() === "") { why = `line ${k} is BLANK in ${r.rel} (a pointer into whitespace points at nothing)`; break; }
        }
        if (!why && nums.length === 2 && nums[0] > nums[1]) why = `range ${spec} is inverted`;
        if (why) brokenLine.push({ where, key, why });
      }
    }
  }

  ok("citation census", `${total} \`path:line\` citations parsed from ${trackedMd.length} markdown files`);
  console.log(`      resolution ladder: ${[...ladderUse].map(([k, v]) => `${k}=${v}`).join(", ")}`);
  console.log(`      allowlist: ${ALLOWLIST.length} entries, all in dated artifacts`);

  for (const b of brokenTarget) bad("citation target resolves", `${b.where}: \`${b.key}\` — ${b.why}`);
  for (const b of brokenLine) bad("citation line is real and non-blank", `${b.where}: \`${b.key}\` — ${b.why}`);

  if (!brokenTarget.length && !brokenLine.length) {
    ok("every file:line citation resolves", `${total - usedAllow.size} checked (${total} total, ${usedAllow.size} allowlisted as dated-snapshot)`);
  }

  for (const a of ALLOWLIST) {
    if (!usedAllow.has(a)) bad("allowlist has no stale entry", `${a.doc}: \`${a.cite}\` is allowlisted but no longer appears in that doc`);
  }
  ok("allowlist is exercised", `${usedAllow.size}/${ALLOWLIST.length} exemptions matched`);
}

console.log(`\n${pass} PASSED, ${fail} FAILED`);
process.exit(fail > 0 ? 1 : 0);
