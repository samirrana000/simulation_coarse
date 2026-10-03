#!/usr/bin/env node
/**
 * evolve/gate-measure.mjs — the gate's MEASUREMENTS, extracted.
 *
 * WHY THIS FILE EXISTS
 * --------------------
 * Every component of evolve.mjs's R was measuring a SUBSET of the thing it
 * named, and every one of them reported 1.0 while doing it. A component that
 * scores 1.0 on 18% of its population is worse than one that scores 0.9 on all
 * of it, because the loop reads 1.0 as "nothing left to do here" and stops
 * looking. These four measurements were each fixed; the fix is the same in
 * every case — measure the POPULATION, not a convenient sample of it:
 *
 *   syntax   walked src/ for *.js (110 files) and left the other 100 modules
 *            — including evolve/evolve.mjs itself — to a different script.
 *   dom      read ids out of src/ui.js alone (117), so the 24 ids bound by
 *            controllers through getElementById/querySelector were unchecked.
 *   coverage counted tests/test_*.js (67) and excluded the 9 scripts/*.mjs
 *            suites, and counted a suite "wired" even when only --slow runs it.
 *   science  read README.md.slice(0, 4000) — 15.6% of a 25.6 KB document — and
 *            scored a keyword hit as proof of honest scope.
 *
 * Each function below reports, alongside its score, WHAT IT LOOKED AT. A
 * number with no denominator is not a measurement.
 *
 * DEPENDENCIES: none. No build step. `node >= 20`.
 *
 * CLI (used by the gate's own fast syntax path):
 *   node --experimental-vm-modules evolve/gate-measure.mjs --syntax-scope <root>
 */

import fs from "node:fs";
import path from "node:path";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import { createRequire } from "node:module";

const SELF = fileURLToPath(import.meta.url);
/** vm.SourceTextModule only exists under --experimental-vm-modules; in the
 *  gate's own process it is absent, which is precisely what makes the helper
 *  CLI below the fast path rather than this file. */
const require = createRequire(import.meta.url);

/* ================================================================== *
 * SYNTAX — every .js/.mjs in the repo, not just src/.
 * ================================================================== */

/**
 * Roots and extensions are DECLARED, not globbed implicitly, and both are
 * reported in the gate output. scripts/check.sh has carried that discipline
 * (plus a per-extension vacuity guard) since it was written; the gate adopts
 * it so the two agree about what "every module" means.
 */
export const SYNTAX_ROOTS = ["src", "scripts", "tests", "bench", "evolve", "tools"];
export const SYNTAX_EXTS = [".js", ".mjs"];
export const SYNTAX_ROOT_FILES = true;

/**
 * Parse every declared file in THIS process and report which do not parse.
 * Returns null when the process was not started with --experimental-vm-modules,
 * which is the signal for checkSyntax() to fall back rather than to pass.
 */
export function checkSyntaxFast(root) {
  let vm = null;
  try { vm = require("node:vm"); } catch { return null; }
  if (!vm || typeof vm.SourceTextModule !== "function") return null;
  const { files } = syntaxScope(root);
  const failed = [];
  for (const f of files) {
    const src = readText(path.join(root, f));
    if (src === null) { failed.push(`${f}: unreadable`); continue; }
    try { new vm.SourceTextModule(src, { identifier: f }); }
    catch (e) { failed.push(`${f}: ${String(e.message).split("\n")[0].slice(0, 90)}`); }
  }
  return { failed, total: files.length };
}

/**
 * The set of files the syntax component is responsible for: every declared
 * extension under every declared root, plus the repo root's own entry points.
 * Returned repo-relative and sorted, so the number is reproducible.
 */
export function syntaxScope(root) {
  const out = [];
  const missingRoots = [];
  for (const r of SYNTAX_ROOTS) {
    const abs = path.join(root, r);
    if (!fs.existsSync(abs) || !fs.statSync(abs).isDirectory()) { missingRoots.push(r); continue; }
    walk(abs, root, out);
  }
  if (SYNTAX_ROOT_FILES) {
    let entries = [];
    try { entries = fs.readdirSync(root, { withFileTypes: true }); } catch {}
    for (const e of entries) {
      if (e.isFile() && SYNTAX_EXTS.some((x) => e.name.endsWith(x))) out.push(e.name);
    }
  }
  out.sort();
  return { files: out, missingRoots };
}

/** Authoritative fallback: one `node --check` per file. Slow, but the same
 *  judgement scripts/check.sh makes, used when the in-process parser is
 *  unavailable so the component can never silently pass. */
export function checkSyntaxWithNode(root) {
  const { files, missingRoots } = syntaxScope(root);
  const failed = [];
  for (const f of files) {
    try { execFileSync(process.execPath, ["--check", path.join(root, f)], { stdio: "pipe" }); }
    catch (e) { failed.push(`${f}: ${String(e.stderr || e.message).split("\n")[0].slice(0, 90)}`); }
  }
  return { failed, total: files.length, missingRoots };
}

/**
 * checkSyntax(root) — the gate's syntax component.
 *
 * WHY A HELPER PROCESS for the fast path: the previous implementation shelled
 * out to `node --check <file>` once per file — 110 spawns for src/ alone,
 * 4.2 s for the whole repo, and this function is called ~14 times inside
 * tests/test_evolve_gate.js. Widening the population without changing the
 * cost would have taken that suite from 18 s to 34 s and pushed the FAST tier
 * past its 60 s budget — which is how the CI gate began failing
 * intermittently. `vm.SourceTextModule` parses (and does NOT execute) a module
 * in ~65 ms for all 210 files.
 *
 * EQUIVALENCE IS NOT ASSUMED: tests/test_gate_measurement.js asserts that this
 * parser and `node --check` agree over a sample of the real corpus, and that
 * both reject the same hand-planted malformed sources. If the fast path ever
 * becomes unavailable, checkSyntax() falls back to spawn-per-file `node
 * --check` rather than silently reporting 1.0 — a gate that cannot run must
 * say so, not pass.
 */
export function checkSyntax(root, { fast = true } = {}) {
  const { files, missingRoots } = syntaxScope(root);
  const byExt = Object.fromEntries(SYNTAX_EXTS.map((e) => [e, 0]));
  for (const f of files) {
    const e = SYNTAX_EXTS.find((x) => f.endsWith(x));
    if (e) byExt[e]++;
  }
  let failed = null, parser = "node --check (spawn per file)";
  if (fast && files.length) {
    try {
      const out = execFileSync(
        process.execPath,
        ["--experimental-vm-modules", "--no-warnings", SELF, "--syntax-scope", root],
        { encoding: "utf-8", stdio: ["ignore", "pipe", "pipe"], timeout: 120000 },
      );
      const res = JSON.parse(out);
      failed = res.failed;
      parser = "vm.SourceTextModule (parse-only, no execution)";
    } catch { failed = null; }
  }
  if (failed === null) {
    const slow = checkSyntaxWithNode(root);
    failed = slow.failed;
    parser = "node --check (fallback: fast parser unavailable)";
  }
  const ok = files.length - failed.length;
  return {
    score: files.length ? ok / files.length : 0,
    numbers: {
      ok, total: files.length, failed, byExt,
      roots: SYNTAX_ROOTS, rootFiles: SYNTAX_ROOT_FILES,
      missingRoots, parser,
    },
  };
}

/* ================================================================== *
 * DOM — every id any src/ module binds, not just src/ui.js's $("...").
 * ================================================================== */

/**
 * Strip comments before extracting ids. Measured necessity, not tidiness:
 * src/controllers/dg-chain-ui.js:21 contains the literal text
 * `$("id")` inside a block comment describing this very gate, and a naive
 * scan reads it as a reference to an element called `id` — one guaranteed
 * false "missing id", i.e. a score depressed for a reason that does not exist.
 */
export function stripJsComments(src) {
  return src
    .replace(/\/\*[\s\S]*?\*\//g, " ")
    .replace(/(^|[^:])\/\/.*$/gm, "$1 ");
}

/** The three ways this codebase reaches an element by id. A bare-id argument
 *  only: `"#controls > .panel"` and `"[data-ex]"` are CSS selectors, not ids,
 *  and treating them as ids would manufacture missing-element reports. */
const ID_ACCESSORS = [
  { name: "$()", re: /\$\(\s*["'`]#?([A-Za-z][\w-]*)["'`]/g },
  { name: "getElementById", re: /getElementById\(\s*["'`]#?([A-Za-z][\w-]*)["'`]/g },
  { name: "querySelector", re: /querySelector(?:All)?\(\s*["'`]#([A-Za-z][\w-]*)["'`]/g },
];

/** id -> Set of modules that bind it. */
export function referencedDomIds(root) {
  const refs = new Map();
  const files = listFiles(path.join(root, "src"), /\.js$/, root);
  const byAccessor = Object.fromEntries(ID_ACCESSORS.map((a) => [a.name, 0]));
  for (const f of files) {
    const src = stripJsComments(readText(path.join(root, f)) || "");
    for (const acc of ID_ACCESSORS) {
      acc.re.lastIndex = 0;
      for (const m of src.matchAll(acc.re)) {
        byAccessor[acc.name]++;
        if (!refs.has(m[1])) refs.set(m[1], new Set());
        refs.get(m[1]).add(f);
      }
    }
  }
  return { refs, files, byAccessor };
}

export function domContract(root) {
  let html = "";
  try { html = fs.readFileSync(path.join(root, "index.html"), "utf-8"); }
  catch { return { score: 0, numbers: { ids: 0, missing: 0, error: "index.html unreadable" } }; }
  const htmlIds = new Set([...html.matchAll(/\bid="([^"]+)"/g)].map((m) => m[1]));
  const { refs, files, byAccessor } = referencedDomIds(root);
  const ids = [...refs.keys()];
  const missing = ids.filter((id) => !htmlIds.has(id));
  const htmlOnly = [...htmlIds].filter((id) => !refs.has(id));
  return {
    score: ids.length ? (ids.length - missing.length) / ids.length : 1,
    numbers: {
      ids: ids.length,
      missing: missing.length,
      missingIds: missing.slice(0, 12),
      modules: files.length,
      byAccessor,
      htmlIds: htmlIds.size,
      htmlOnlyCount: htmlOnly.length,
    },
  };
}

/* ================================================================== *
 * COVERAGE — every test-bearing file, and only the ones a gate runs.
 * ================================================================== */

/**
 * The population is DECLARED per directory, and each directory's pattern is
 * reported, so adding a test to a directory that is not listed here is a
 * visible omission rather than a silent one.
 *
 *   tests/test_*.js        — the main suites, registered in tests/suites.js
 *   scripts/{test,validate}_*.mjs  — nine real suites that the old coverage
 *                            metric did not count AT ALL
 *   bench/{test,validate}_*.js     — same convention, currently empty; listed
 *                            so that adding one is covered by construction
 */
export const TEST_POPULATION = [
  { dir: "tests", re: /^test_.*\.js$/, label: "tests/test_*.js" },
  { dir: "scripts", re: /^(test|validate)_.*\.mjs$/, label: "scripts/{test,validate}_*.mjs" },
  { dir: "bench", re: /^(test|validate)_.*\.js$/, label: "bench/{test,validate}_*.js" },
];
/** The tier `npm test` and .github/workflows/check.yml actually execute. */
export const DEFAULT_TIER = "FAST";

/**
 * The runner is not a suite. tests/test_all.js matches tests/test_*.js, and
 * the old countTests() counted it — as WIRED, because tests/suites.js mentions
 * "test_all.js" in a comment on line 2 and the wiring test was a substring
 * search. So one of the 67 files that made coverage read 67/67 = 1.0 was the
 * harness that RUNS the other 66, "covered" by a sentence of prose. A
 * registry ENTRY is now required, and the harness is excluded by name.
 */
export const NOT_SUITES = ["tests/test_all.js"];

export function testPopulation(root) {
  const suites = [];
  const bySource = {};
  for (const { dir, re, label } of TEST_POPULATION) {
    const files = listFiles(path.join(root, dir), re, root).filter((f) => !NOT_SUITES.includes(f));
    bySource[label] = files.length;
    for (const f of files) suites.push({ file: f, label });
  }
  const manual = listFiles(path.join(root, "tests", "manual"), /^test_.*\.js$/, root);

  let reg = "";
  try { reg = fs.readFileSync(path.join(root, "tests", "suites.js"), "utf-8"); } catch {}
  const tiers = new Map();
  for (const m of reg.matchAll(/file:\s*"([^"]+)"\s*,\s*tier:\s*"([A-Z]+)"\s*(?:,\s*timeout:\s*(\d+))?/g)) {
    tiers.set(m[1], { tier: m[2], timeout: m[3] ? Number(m[3]) : null });
  }

  let wired = 0, defaultRun = 0, optIn = 0, unwired = 0;
  const unwiredFiles = [], optInFiles = [];
  for (const s of suites) {
    const t = tiers.get(s.file);
    if (!t) { unwired++; unwiredFiles.push(s.file); continue; }
    wired++;
    if (t.tier === DEFAULT_TIER) defaultRun++;
    else { optIn++; optInFiles.push(`${s.file} [${t.tier}]`); }
  }
  const total = suites.length + manual.length;
  return {
    total, wired, defaultRun, optIn, unwired,
    manual: manual.length,
    unwiredFiles: unwiredFiles.slice(0, 12),
    optInFiles: optInFiles.slice(0, 12),
    bySource,
    excludedHarness: NOT_SUITES,
    definition: `total = every file matching ${TEST_POPULATION.map((p) => p.label).join(" / ")} plus tests/manual/, minus the runner (${NOT_SUITES.join(", ")}); a suite counts as gated only when tests/suites.js registers it at tier ${DEFAULT_TIER}, the tier \`npm test\` and CI execute`,
  };
}

export function coverageScore(root) {
  const c = testPopulation(root);
  if (!c.total) return { score: 1, numbers: c };
  return { score: Math.min(1, (c.defaultRun + c.manual) / c.total), numbers: c };
}

/* ================================================================== *
 * SCIENCE — the whole honest-scope surface, and the absence of an overclaim.
 * ================================================================== */

export const MUST_FILES = [
  "ROADMAP.md", "docs/LIMITATIONS.md", "docs/APPLICABILITY.md",
  "docs/VALIDATION.md", "LICENSE", "CITATION.cff", "CHANGELOG.md",
];

/** The documents that actually carry the trust boundary. ROADMAP.md §1 is the
 *  out-of-scope contract; the other three are what a reader is sent to. */
export const TRUST_SOURCES = [
  "README.md", "ROADMAP.md", "docs/LIMITATIONS.md",
  "docs/APPLICABILITY.md", "docs/VALIDATION.md",
];

/** Language that states a boundary. Counted per document, not per file: a
 *  single phrase anywhere in the README used to be enough to score the whole
 *  honest-scope surface 1.0. */
export const TRUST_PHRASE =
  /trust boundary|not a replacement|qualitative|not FEP|out[- ]of[- ]scope|out of scope|will NOT do/gi;

/**
 * Claims that are false for THIS model whatever the surrounding prose says.
 * Each is a sentence that cannot appear in an honest document: ROADMAP.md §1
 * rules out FEP/TI/MBAR, QM/MM, membranes and PME, and forbids calling the
 * toy 4-state kinetics publishable. The check is for the ABSENCE of these —
 * the old component could only ever be made to pass by adding a phrase, never
 * to fail by writing a wrong claim, which is the wrong direction for a gate
 * guarding honest scope.
 */
export const OVERCLAIM_PATTERNS = [
  { name: "production-ready", re: /\bproduction[- ]ready\b/gi },
  { name: "AMBER-quality", re: /\bAMBER[- ]quality\b/gi },
  { name: "GROMACS-equivalent", re: /\bGROMACS[- ]equivalent\b/gi },
  { name: "quantitative-free-energy", re: /quantitative(?:ly)? (?:binding )?free energ/gi },
  { name: "converged-affinity", re: /converged (?:binding )?affinit/gi },
  { name: "trust-boundary-met", re: /trust boundary (?:is )?(?:satisfied|met|closed)/gi },
  { name: "experimentally-validated", re: /\bexperimentally validated\b/gi },
  { name: "exact-binding-affinity", re: /\bexact\b[^.\n]{0,24}binding affinity/gi },
];

export function scienceSurface(root) {
  const have = MUST_FILES.filter((f) => fs.existsSync(path.join(root, f)));
  const missingMust = MUST_FILES.filter((f) => !fs.existsSync(path.join(root, f)));

  const sources = [];
  let charsScanned = 0;
  for (const f of TRUST_SOURCES) {
    const text = readText(path.join(root, f));
    if (text === null) { sources.push({ file: f, present: false, chars: 0, phrases: 0 }); continue; }
    charsScanned += text.length;
    const phrases = [...text.matchAll(TRUST_PHRASE)].length;
    sources.push({ file: f, present: true, chars: text.length, phrases });
  }
  const withBoundary = sources.filter((s) => s.present && s.phrases > 0).length;

  // Reachability: a trust-boundary document nobody links is not part of the
  // honest-scope surface, it is a file. Measured on README.md, whose first
  // link to docs/LIMITATIONS.md sits at char 5038 — past the 4000-char window
  // the old component read, i.e. entirely invisible to it.
  const readme = readText(path.join(root, "README.md")) || "";
  const linkable = TRUST_SOURCES.filter((f) => f !== "README.md");
  const linked = linkable.filter((f) => readme.includes(f));

  const overclaims = [];
  for (const file of TRUST_SOURCES) {
    const text = readText(path.join(root, file));
    if (text === null) continue;
    for (const p of OVERCLAIM_PATTERNS) {
      p.re.lastIndex = 0;
      const hits = [...text.matchAll(p.re)];
      if (hits.length) {
        overclaims.push({ file, pattern: p.name, count: hits.length, sample: String(hits[0][0]).slice(0, 60) });
      }
    }
  }

  let score = (have.length / MUST_FILES.length) * (withBoundary / TRUST_SOURCES.length);
  score *= linked.length / linkable.length;            // the boundary must be reachable
  for (let i = 0; i < overclaims.length; i++) score *= 0.5;   // each overclaim halves it

  return {
    score: Math.max(0, Math.min(1, score)),
    numbers: {
      must: MUST_FILES.length, have: have.length, missingMust,
      sources, charsScanned,
      charsTotal: sources.reduce((a, s) => a + s.chars, 0),
      withBoundary, trustSources: TRUST_SOURCES.length,
      linked: linked.length, linkable: linkable.length,
      linkedDocs: linkable,
      unlinkedDocs: linkable.filter((f) => !readme.includes(f)),
      overclaims,
      definition: "score = (must files present) x (trust-boundary sources carrying boundary language) x (those docs linked from README.md) x 0.5^(overclaim hits); the whole of every source is scanned",
    },
  };
}

/* ================================================================== *
 * shared helpers
 * ================================================================== */

function readText(p) { try { return fs.readFileSync(p, "utf-8"); } catch { return null; } }

function walk(dir, root, out) {
  let entries = [];
  try { entries = fs.readdirSync(dir, { withFileTypes: true }); } catch { return; }
  for (const e of entries) {
    if (e.name === "node_modules" || e.name.startsWith(".")) continue;
    const p = path.join(dir, e.name);
    if (e.isDirectory()) walk(p, root, out);
    else if (SYNTAX_EXTS.some((x) => e.name.endsWith(x))) {
      out.push(path.relative(root, p).split(path.sep).join("/"));
    }
  }
}

export function listFiles(dir, re, root = process.cwd(), out = []) {
  let entries = [];
  try { entries = fs.readdirSync(dir, { withFileTypes: true }); } catch { return out; }
  for (const e of entries) {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) listFiles(p, re, root, out);
    else if (re.test(e.name)) out.push(path.relative(root, p).split(path.sep).join("/"));
  }
  return out;
}

/* CLI: the fast syntax child. Prints one JSON object and nothing else. */
if (process.argv[2] === "--syntax-scope") {
  const root = path.resolve(process.argv[3] || process.cwd());
  const res = checkSyntaxFast(root);
  process.stdout.write(JSON.stringify(res || { failed: ["__no_fast_parser__"], total: 0 }));
}
