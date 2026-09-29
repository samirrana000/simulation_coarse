#!/usr/bin/env node
/**
 * tests/test_cache_contract.js — the app must never serve a stale module.
 *
 * THE DEFECT THIS GUARDS AGAINST
 * -------------------------------
 * Until 2026-09-29 every module edge carried a hand-typed cache-bust
 * suffix: `import { ForceField } from "./forcefield.js?v=10"`. Measured at
 * the time: 121 literals across 42 files under src/ plus one on index.html's
 * module tag, all frozen at the value `10` while src/version.js said
 * `1.1.0-fp7`. A module edited without also editing its neighbours' literals
 * kept the same url, the browser served the cached bytes, and the symptom was
 * "my fix did nothing" — indistinguishable from a physics bug, on a tool
 * whose whole premise is interactive exploration.
 *
 * WHAT IS GUARDED, AND HOW
 * ------------------------
 * 1. No query string on any module specifier, anywhere in shipped code
 *    (src/**, index.html, sw.js) or in node-side code that imports the app
 *    (tests/**, bench/**). Comments are stripped with a real tokenizer
 *    first, so prose about `?v=` cannot satisfy or trip the rule.
 * 2. The mechanism that replaced it exists and is wired up: sw.js, whose
 *    only response path is `fetch(request, { cache: "no-store" })` and which
 *    holds no cache of its own (no `caches.open`, no `.match(`), registered
 *    from index.html with a query-free url, `updateViaCache: "none"`, an
 *    availability guard and a bounded one-shot reload.
 * 3. That mechanism actually does what it claims — sw.js is executed in a
 *    node:vm sandbox with a recording `fetch`/`caches`, and the recorded
 *    calls are asserted. This is a mechanism test, not a grep: a sw.js that
 *    stopped passing `no-store` goes red here even though every literal is
 *    gone and the file still looks plausible.
 * 4. The single-instance invariant the old literals quietly broke. Two
 *    test files had to import `../src/ui.js?v=10` to reach the same `state`
 *    object the panels read, because `ui.js` and `ui.js?v=10` are two
 *    module instances. With the suffixes gone there is exactly one, and
 *    that is asserted by importing both ends for real.
 *
 * Vacuity guards: the file walk must find at least MIN_SHIPPED_MODULES
 * shipped modules and MIN_APP_MODULES import edges, so a broken glob or a
 * renamed directory reads as FAIL, never as "nothing to complain about".
 *
 * Deliberately NOT scanned: docs/, *.md, AUDIT_REPORT.md and
 * TRANSFORMATION_PLAN_100.md. Those are prose and history — an audit report
 * that quotes the `?v=10` it diagnosed is a record, not a defect — and
 * evolve/evolve.mjs, whose probe 8 exists precisely to re-detect literals
 * in src/ if they ever come back.
 *
 * Run: node tests/test_cache_contract.js
 * Tier: FAST (registered in tests/suites.js)
 */

import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..");

/** src/ must hold at least this many modules for the walk to be believed. */
const MIN_SHIPPED_MODULES = 60;
/** The app must still have at least this many module edges to be scanned. */
const MIN_APP_MODULES = 100;

let passed = 0;
let failed = 0;
function assert(condition, message, detail) {
  if (condition) {
    passed++;
    console.log(`  ✓ ${message}`);
  } else {
    failed++;
    console.error(`  ✗ FAIL: ${message}${detail ? `\n      ${detail}` : ""}`);
  }
}

function rel(p) {
  return path.relative(ROOT, p).split(path.sep).join("/");
}

function walk(dir, re, out = []) {
  if (!fs.existsSync(dir)) return out;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) walk(full, re, out);
    else if (re.test(entry.name)) out.push(full);
  }
  return out;
}

/**
 * Strip comments while collecting the source's string literals.
 *
 * A regex that deletes `//.*$` would mangle `"https://esm.sh"` and, worse,
 * could truncate a line that holds a real specifier after a URL in a
 * comment — a guard that can be blinded is not a guard. So this is a small
 * state machine over the three quote kinds, \ escapes, both comment forms,
 * and REGEX LITERALS.
 *
 * Regex literals are the part that is easy to get wrong, and the first
 * version of this file got it wrong: `/!\(\s*["']serviceWorker["']\s+in\s+navigator\s*\)/`
 * contains a quote character, so a scanner that only knows about strings
 * opens a string at the first `"`, runs on to the end of the file, and then
 * reports a clean scan of code it never really read (that is how
 * `"../src/ui.js?v=10"` inside a *comment* of this very file was reported as
 * a live specifier). So: a `/` opens a regex literal when the previous
 * significant character is one of `( , = : [ ! & | ? { } ;` (or the previous
 * token is `return`, `typeof`, `case`, `in`, `of`, `new`, …) — the standard
 * disambiguation — and the regex is consumed to its closing unescaped `/`,
 * honouring `[...]` character classes.
 *
 * Second safety net: an unterminated string or comment at EOF throws. A
 * scanner that could not read a file must fail the test loudly; it must
 * never report that file as clean.
 *
 * @returns {{code: string, strings: {value: string, line: number}[]}}
 */
const REGEX_PRECEDER_CHARS = new Set([..."(,=:[!&|?{};+-*%~^<>"]);
const REGEX_PRECEDER_WORDS = new Set([
  "return", "typeof", "instanceof", "in", "of", "new", "delete", "void",
  "case", "do", "else", "yield", "await", "throw",
]);
const WORD_CHAR = /[A-Za-z0-9_$]/;

export function stripComments(src) {
  let code = "";
  const strings = [];
  let i = 0;
  let line = 1;
  const n = src.length;
  let lastSig = "";        // last non-whitespace character emitted
  let lastWord = "";       // last identifier token emitted
  const emit = (text) => {
    code += text;
    const trimmed = text.replace(/\s+$/, "");
    if (!trimmed) return;
    lastSig = trimmed[trimmed.length - 1];
    if (WORD_CHAR.test(lastSig)) {
      const m = /[A-Za-z0-9_$]+$/.exec(trimmed);
      lastWord = m ? m[0] : "";
    } else {
      lastWord = "";
    }
  };
  const regexAllowed = () =>
    lastSig === "" || REGEX_PRECEDER_CHARS.has(lastSig) || REGEX_PRECEDER_WORDS.has(lastWord);

  while (i < n) {
    const c = src[i];
    const c2 = src[i + 1];
    if (c === "\n") line++;
    if (c === "/" && c2 === "/") {
      while (i < n && src[i] !== "\n") i++;
      continue;
    }
    if (c === "/" && c2 === "*") {
      const startLine = line;
      i += 2;
      while (i < n && !(src[i] === "*" && src[i + 1] === "/")) {
        if (src[i] === "\n") line++;
        i++;
      }
      if (i >= n) throw new Error(`unterminated block comment opened before line ${startLine}`);
      i += 2;
      emit(" ");
      continue;
    }
    if (c === "/" && regexAllowed()) {
      const start = i;
      i++;
      let inClass = false;
      let closed = false;
      while (i < n) {
        if (src[i] === "\\") { i += 2; continue; }
        if (src[i] === "\n") break;              // not a regex after all
        if (src[i] === "[") inClass = true;
        else if (src[i] === "]") inClass = false;
        else if (src[i] === "/" && !inClass) { i++; closed = true; break; }
        i++;
      }
      if (closed) {
        while (i < n && WORD_CHAR.test(src[i])) i++;   // flags
        emit(src.slice(start, i));
        continue;
      }
      i = start;                                // fall through as a plain '/'
    }
    if (c === '"' || c === "'" || c === "`") {
      const quote = c;
      const startLine = line;
      const start = i;
      let value = "";
      i++;
      while (i < n && src[i] !== quote) {
        if (src[i] === "\\") { value += src[i + 1] ?? ""; i += 2; continue; }
        if (src[i] === "\n") line++;
        value += src[i];
        i++;
      }
      if (i >= n) throw new Error(`unterminated ${quote} string opened before line ${startLine}`);
      i++;
      strings.push({ value, line: startLine });
      emit(src.slice(start, i));   // keep the literal itself: only comments are blanked
      continue;
    }
    emit(c);
    i++;
  }
  return { code, strings };
}

/** A string literal that looks like a module url carrying a query string. */
function isQueriedModuleSpec(value) {
  return value.includes("?") && /\.js(\?|$)/.test(value);
}

const shipped = [
  ...walk(path.join(ROOT, "src"), /\.js$/),
  path.join(ROOT, "sw.js"),
];
const nodeSide = [
  ...walk(path.join(ROOT, "tests"), /\.js$/),
  ...walk(path.join(ROOT, "bench"), /\.js$/),
];
const appModules = walk(path.join(ROOT, "src"), /\.js$/);

/* ---------------------------------------------------------------- 0. */
/* The scanner's own self-test, then the vacuity guards.                 */
/* ---------------------------------------------------------------- */
{
  // Assembled, not written out, so this file does not trip its own scan.
  const Q = "?v" + "=10";
  const probe = [
    'const re = /["\']serviceWorker["\']/;',
    '// a comment mentioning "../src/ui.js' + Q + '" must be ignored',
    'import { X } from "./x.js' + Q + '";',
  ].join("\n");
  let got;
  try {
    got = stripComments(probe).strings.filter((s) => isQueriedModuleSpec(s.value)).map((s) => s.value);
  } catch (e) {
    got = `threw: ${e.message}`;
  }
  assert(
    Array.isArray(got) && got.length === 1 && got[0] === `./x.js${Q}`,
    "the scanner is not blind: a quote inside a regex literal does not desynchronise it, and prose in a comment is ignored",
    `scanner returned ${JSON.stringify(got)} — if this fails, every other scan in this file is untrustworthy`
  );
}

/* Vacuity: a walk that finds nothing must not read as "all clear".     */
const srcModules = appModules.filter((f) => f !== path.join(ROOT, "sw.js"));
assert(
  srcModules.length >= MIN_SHIPPED_MODULES,
  `walk found ${srcModules.length} shipped modules (>= ${MIN_SHIPPED_MODULES} required, so this scan cannot be vacuous)`,
  `expected >= ${MIN_SHIPPED_MODULES} files under src/; a wrong ROOT or a renamed directory reads as PASS otherwise`
);

let importEdges = 0;
let scanned = 0;
const unreadable = [];
for (const f of appModules) {
  try {
    const { code } = stripComments(fs.readFileSync(f, "utf-8"));
    scanned++;
    importEdges += (code.match(/\bfrom\s*["']/g) || []).length +
      (code.match(/\bimport\s*\(\s*["']/g) || []).length;
  } catch (e) {
    unreadable.push(`${rel(f)}: ${e.message}`);
  }
}
assert(
  unreadable.length === 0,
  `every src/ module parsed cleanly (${scanned}/${appModules.length}) — an unparseable file is a blind spot, not a pass`,
  unreadable.slice(0, 5).join("; ")
);
assert(
  importEdges >= MIN_APP_MODULES,
  `walk resolved ${importEdges} module edges (>= ${MIN_APP_MODULES} required)`,
  "if the import-edge scan finds nothing, assertion 1 below is vacuous"
);

/* ---------------------------------------------------------------- 1. */
/* No query string on any module specifier, in shipped or node-side code.*/
/* ---------------------------------------------------------------- */
{
  // This scanner has to name the token it forbids, so it is the one file
  // exempt from the raw-token rule below — a guard cannot ban the string it
  // must search for. The specifier rule still applies to it, and the exemption
  // is asserted to be exactly this one file rather than a growing list.
  const SELF = fileURLToPath(import.meta.url);
  const files = [...shipped, ...nodeSide];
  const offenders = [];
  let tokenScanned = 0;
  let parsed = 0;
  for (const f of files) {
    let code, strings;
    try {
      ({ code, strings } = stripComments(fs.readFileSync(f, "utf-8")));
      parsed++;
    } catch (e) {
      offenders.push(`${rel(f)}  ${e.code === "ENOENT" ? "MISSING" : "UNPARSEABLE"} (${e.message}) — a file the rules cannot read is reported, never passed`);
      continue;
    }
    for (const s of strings) {
      if (isQueriedModuleSpec(s.value)) {
        offenders.push(`${rel(f)}:${s.line}  specifier ${JSON.stringify(s.value)}`);
      }
    }
    if (f === SELF) continue;
    tokenScanned++;
    for (const m of code.matchAll(/\?v=/g)) {
      const ln = code.slice(0, m.index).split("\n").length;
      offenders.push(`${rel(f)}:${ln}  literal ?v= in code`);
    }
  }
  assert(
    parsed === files.length,
    `all ${files.length} shipped+node files parsed cleanly by the scanner`,
    "an unparseable file is a blind spot in every rule below"
  );
  assert(
    tokenScanned === files.length - 1,
    `the version-token rule is exempt for exactly one file (this scanner); ${tokenScanned} of ${files.length} files were token-scanned, all ${files.length} specifier-scanned`,
    "a second exemption means the guard is being widened instead of the defect being fixed"
  );
  assert(
    offenders.length === 0,
    `zero version query literals in ${files.length} shipped+node files (${importEdges} edges scanned)`,
    offenders.slice(0, 10).join("\n      ") + (offenders.length > 10 ? `\n      … and ${offenders.length - 10} more` : "")
  );
}

{
  const rawHtml = fs.readFileSync(path.join(ROOT, "index.html"), "utf-8");
  // Comments are prose, not code: strip them before the literal scan so a
  // sentence that *describes* the rule cannot satisfy or trip it.
  const html = rawHtml.replace(/<!--[\s\S]*?-->/g, " ");
  const offenders = [];
  for (const m of html.matchAll(/(?:src|href)\s*=\s*"([^"]*)"/g)) {
    if (m[1].includes("?")) offenders.push(`${m[1]}  (attribute src/href)`);
  }
  for (const m of html.matchAll(/\?v=/g)) {
    offenders.push(`line ${html.slice(0, m.index).split("\n").length}: literal ?v=`);
  }
  assert(
    offenders.length === 0,
    "index.html carries no query string on any src/href and no ?v= literal",
    offenders.join("; ")
  );
  assert(
    /<script\s+type="module"\s+src="src\/main\.js"><\/script>/.test(html),
    "index.html loads src/main.js as a query-free module script",
    "expected exactly <script type=\"module\" src=\"src/main.js\"></script>"
  );
}

/* ---------------------------------------------------------------- 2. */
/* The mechanism exists, and it cannot serve anything it did not fetch.   */
/* ---------------------------------------------------------------- */
const swPath = path.join(ROOT, "sw.js");
const swSrc = fs.existsSync(swPath) ? fs.readFileSync(swPath, "utf-8") : "";
const swCode = stripComments(swSrc).code;
{
  assert(swSrc.length > 0, "sw.js exists and is non-empty", "sw.js is the freshness mechanism; an empty file leaves every reload at the server's mercy");
  assert(
    /cache:\s*["']no-store["']/.test(swCode),
    "sw.js re-fetches with cache:\"no-store\" (bypasses the HTTP cache entirely, so heuristic freshness cannot apply)",
    "no `cache: \"no-store\"` in sw.js code — a bare fetch() would let the browser reuse a cached copy of an edited module"
  );
  assert(
    !/caches\.open\s*\(/.test(swCode) && !/\.match\s*\(\s*request|\.match\s*\(\s*req\b|cache\.match\s*\(/.test(swCode),
    "sw.js has no cache-read path (no caches.open, no cache.match) — its only response is a fresh network response",
    "a cache read here would be a stored copy of the app, i.e. the original bug back"
  );
  assert(
    !/["'`][^"'`]*\?[^"'`]*\.js/.test(swCode),
    "sw.js adds no query string to any url it handles"
  );
  assert(
    /self\.addEventListener\(\s*["']install["']/.test(swCode) && /skipWaiting/.test(swCode),
    "sw.js installs without waiting (covers the first visit too)"
  );
  assert(
    /self\.addEventListener\(\s*["']activate["']/.test(swCode) && /clients\.claim/.test(swCode),
    "sw.js claims open clients on activate"
  );
}

{
  const html = fs.readFileSync(path.join(ROOT, "index.html"), "utf-8");
  const reg = /serviceWorker\.register\(\s*["']([^"']+)["']/.exec(html);
  assert(!!reg, "index.html registers the worker");
  assert(
    reg && !reg[1].includes("?"),
    `the registered worker url is query-free (${reg ? reg[1] : "n/a"}) — the browser's own byte comparison keeps the mechanism itself fresh`,
    "a query on the worker url would need a version literal, i.e. the defect, one level up"
  );
  assert(
    /updateViaCache:\s*["']none["']/.test(html),
    "index.html registers with updateViaCache:\"none\" (worker script always revalidated)"
  );
  assert(
    /serviceWorker\.register[\s\S]{0,200}?\.catch\s*\(/.test(html),
    "registration failure is caught (insecure origin / file:// must not break the app)"
  );
  assert(
    /!\(\s*["']serviceWorker["']\s+in\s+navigator\s*\)/.test(html),
    "availability is guarded (`serviceWorker` in navigator) before use"
  );
  assert(
    /sessionStorage/.test(html) && /location\.reload\(\)/.test(html),
    "the one-shot post-activation reload is bounded by sessionStorage (no reload loop)"
  );
}

/* ---------------------------------------------------------------- 3. */
/* Execute sw.js and assert what it actually does.                      */
/* ---------------------------------------------------------------- */
if (swSrc.trim().length === 0) {
  // Reported, not crashed on: a missing worker must produce a counted
  // failure and a summary line, because a bare stack trace here would be
  // indistinguishable from a broken harness.
  assert(false, "sw.js could not be executed — the freshness mechanism is absent (reported again above, and by scripts/check.sh)");
} else {
  const listeners = new Map();
  const calls = { fetch: [], cache: [] };
  const ORIGIN = "http://127.0.0.1:8471";
  const NET = new Response("fresh-bytes", { status: 200 });
  const sandbox = {
    console,
    URL,
    Response,
    fetch: (request, init) => {
      calls.fetch.push({ request, init });
      return Promise.resolve(NET);
    },
    caches: {
      keys: () => Promise.resolve([...calls.cache]),
      delete: (k) => { calls.cache.delete(k); return Promise.resolve(true); },
    },
    self: {
      location: { origin: ORIGIN },
      addEventListener: (type, fn) => listeners.set(type, fn),
      skipWaiting: () => { calls.skipped = true; return Promise.resolve(); },
      clients: { claim: () => { calls.claimed = true; return Promise.resolve(); } },
    },
  };
  vm.createContext(sandbox);
  vm.runInContext(swSrc, sandbox, { filename: "sw.js" });

  const on = (type) => {
    const fn = listeners.get(type);
    assert(typeof fn === "function", `sw.js registers a "${type}" handler`);
    return fn;
  };
  on("install"); on("activate"); on("fetch");

  // install / activate
  sandbox.self.skipWaiting = () => { calls.skipped = true; return Promise.resolve(); };
  calls.cache = new Set(["simcoarse-1.0.0-transform"]);
  const actEvent = { waitUntil: (p) => { calls.waitUntil = p; } };
  listeners.get("install")({ waitUntil: (p) => { calls.waitUntil = p; } });
  listeners.get("activate")(actEvent);
  assert(calls.skipped === true, "install handler calls skipWaiting()");
  if (calls.waitUntil) await calls.waitUntil;
  assert(calls.claimed === true, "activate handler claims clients");
  assert(calls.cache.size === 0, "activate handler deletes leftover Cache Storage from an older build", `left: ${[...calls.cache].join(", ")}`);

  // same-origin GET -> intercepted, re-fetched with no-store
  const get = (url) => ({ method: "GET", url });
  const ev1 = { request: get(`${ORIGIN}/src/main.js`), respondWith: (p) => { ev1.p = p; }, waitUntil() {} };
  listeners.get("fetch")(ev1);
  assert(ev1.p !== undefined, "same-origin GET is intercepted (respondWith called)");
  assert(calls.fetch.length === 1 && calls.fetch[0].init && calls.fetch[0].init.cache === "no-store",
    `the intercepted module fetch is issued with cache:"no-store" (saw ${JSON.stringify(calls.fetch[0] && calls.fetch[0].init)})`,
    "this is the whole mechanism: without no-store the browser may answer from its HTTP cache");
  assert(calls.fetch[0] && calls.fetch[0].request.url === `${ORIGIN}/src/main.js`,
    "the original request (url included) is what gets re-fetched, not a rewritten one");
  if (ev1.p) assert((await ev1.p) === NET, "the network response is returned verbatim (no stored copy can shadow it)");

  // cross-origin + non-GET -> untouched
  const ev2 = { request: get("https://esm.sh/x.js"), respondWith: () => { ev2.hit = true; }, waitUntil() {} };
  const ev3 = { request: { method: "POST", url: `${ORIGIN}/x` }, respondWith: () => { ev3.hit = true; }, waitUntil() {} };
  listeners.get("fetch")(ev2);
  listeners.get("fetch")(ev3);
  assert(ev2.p === undefined && !ev2.hit, "cross-origin GET is left to the default network path");
  assert(ev3.p === undefined && !ev3.hit, "non-GET request is left to the default network path");
  assert(calls.fetch.length === 1, "only the same-origin GET reached the network (2 further events ignored)", `fetch was called ${calls.fetch.length} times`);
}

/* ---------------------------------------------------------------- 4. */
/* One ui.js instance, not two.                                          */
/* ---------------------------------------------------------------- */
{
  // Before the suffixes were removed, tests/test_placement_hetero.js and
  // tests/test_rev1_issue4_live_terms.js had to import "../src/ui.js?v=10"
  // to reach the same `state` object the panels read: ui.js and ui.js?v=10
  // are two module instances, and a panel writing one while the app read the
  // other is a silent, physics-shaped bug. The suffix is gone; these
  // assertions exist so a second instance cannot come back unnoticed.
  const a = await import("../src/ui.js");
  const b = await import("../src/ui.js");
  assert(a === b, "two imports of src/ui.js resolve to the identical module namespace (one instance, one `state`)");
  assert(typeof a.state === "object" && a.state !== null, "src/ui.js exports a live `state` object under Node");

  const panel = stripComments(fs.readFileSync(path.join(ROOT, "src", "ligand-panel.js"), "utf-8"));
  const uiEdges = panel.strings.filter((s) => /(^|\/)ui\.js$/.test(s.value));
  assert(
    uiEdges.length > 0 && uiEdges.every((s) => s.value === "./ui.js"),
    `ligand-panel.js reaches ui.js through the plain specifier only (${uiEdges.map((s) => s.value).join(", ") || "none"})`,
    "a queried twin specifier here is a second module instance with its own state"
  );
}

/* ---------------------------------------------------------------- */
console.log(`\n${passed} passed, ${failed} failed`);
if (failed > 0) {
  console.error("FAILED");
  process.exitCode = 1;
} else {
  console.log("PASS");
}
