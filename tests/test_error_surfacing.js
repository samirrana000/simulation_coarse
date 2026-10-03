#!/usr/bin/env node
/**
 * tests/test_error_surfacing.js — no error may be silently swallowed.
 *
 * THE DEFECT THIS GUARDS AGAINST
 * ------------------------------
 * Measured before this suite existed: 181 `catch` blocks under src/, of which
 * 92 had an empty body and 56 of those had no annotation either — a bare
 * `catch (_) {}`. On a simulator whose entire premise is that its numbers mean
 * something, that is not a style problem. A user drops a structure, a step of
 * the pipeline throws, the code swallows it and keeps drawing, and the user is
 * watching a plausible animation that is not the simulation they asked for —
 * indistinguishable from a physics bug, and corrosive to trust in every number
 * the tool reports.
 *
 * THE RULE (src/errors.js header states it; this suite enforces it)
 * -----------------------------------------------------------------
 * Every `catch` must be handled in exactly one of four ways:
 *   (a) rethrown
 *   (b) a specific, user-visible recovery
 *   (c) recorded via recordError() — counted, and surfaced on the top bar
 *   (d) a documented no-op via ignore() — counted as "ignored", never silent
 * A body that is empty (comment-only or not) is (d) minus the documentation and
 * the count, which is the thing this suite exists to delete permanently.
 *
 * WHAT IS ASSERTED
 * ----------------
 * 1. GUARD (the rot-stopper): zero bare catches in src/. A catch counts as
 *    handled when its body is non-empty after comment stripping, OR it carries
 *    an explicit `throw`. Comments and doc comments are stripped with a real
 *    scanner first, so prose about `catch {}` cannot satisfy or trip the rule.
 * 2. THE RECORDER ACTUALLY RECORDS AND COUNTS: same error twice → one slot with
 *    count 2 (a per-frame failure must not become 60 slots/second); different
 *    errors → separate slots; the ring is bounded; the counters never decrement.
 * 3. A DELIBERATELY FAILING OPERATION INCREMENTS THE COUNTER AND SURFACES
 *    THROUGH THE STATUS PATH: the recorder is driven from a fake DOM and from a
 *    real global-error handler, and the resulting text is asserted against what
 *    index.html + css/style.css + src/controllers/tick.js will actually render.
 * 4. THE SURFACE EXISTS IN THE MARKUP: #topErrors is in index.html, is wired
 *    into src/ui.js, is inside the EXISTING #topStatus bar (no new panel, so
 *    the Digit1-7 hotkey contract is untouched), and is hidden by CSS when empty.
 * 5. NO EXECUTING STUB: a TODO/FIXME/not yet/placeholder/stub keyword surviving
 *    a comment strip on a line that executes is a value a user can see. Zero
 *    allowed. Comment-line hits are the project's honest scope notes
 *    (ROADMAP.md §1) and are explicitly NOT this suite's business.
 *
 * 6. THE RECORDER CANNOT BECOME A FAILURE MODE: a hostile value, a null
 *    context, a 5000-char message and a DOM whose setter throws must all be
 *    survivable, because this code runs ON the failure path. src/errors.js's own
 *    terminal guards are the one legitimate place a catch body is empty — they
 *    are annotated as such and asserted separately below, because a recorder
 *    that recursed into itself would be worse than the silence it replaces.
 *
 * VACUITY GUARDS
 * --------------
 * MIN_CATCHES / MIN_FILES: a walk that silently finds nothing would report a
 * clean bill of health for a src/ tree that does not exist. A broken glob or a
 * renamed directory must read as FAIL.
 *
 * Deliberately NOT scanned: docs/, *.md, AUDIT_REPORT.md and
 * TRANSFORMATION_PLAN_100.md — those are prose and dated history.
 *
 * Run: node tests/test_error_surfacing.js
 * Tier: FAST (registered in tests/suites.js)
 */

import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";
import { fileURLToPath } from "node:url";

import {
  MAX_ERROR_SLOTS,
  ERROR_PAINT_MIN_MS,
  recordError,
  ignore,
  errorSummary,
  errorCounterText,
  errorHudPrefix,
  errorReport,
  resetErrors,
  installGlobalErrorHandlers,
  drainBootstrapQueue,
  setErrorSink,
} from "../src/errors.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..");

/** src/ must hold at least this many modules for the walk to be believed. */
const MIN_FILES = 60;
/** …and at least this many catch blocks, or the guard below is vacuous. */
const MIN_CATCHES = 150;

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

/** Every .js under src/, recursively, sorted for stable output. */
function srcModules() {
  const out = [];
  (function walk(dir) {
    for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
      const p = path.join(dir, e.name);
      if (e.isDirectory()) walk(p);
      else if (e.name.endsWith(".js")) out.push(p);
    }
  })(path.join(ROOT, "src"));
  return out.sort();
}

/**
 * Blank out comments, preserving offsets and line structure, so a comment can
 * neither satisfy nor trip a body-content rule. Handles line comments to EOL
 * and block comments (the doc-comment form). String literals are NOT stripped:
 * a string that says "stub:" is a value a user can see, and that is the whole
 * point of check 5.
 */
function stripComments(src) {
  let out = "";
  let i = 0;
  const n = src.length;
  // A hand-rolled scanner, because the obvious regexes all lie here. Three
  // things bit this suite on its first run and each is recorded, because a
  // comment stripper that is wrong is WORSE than none: it decides which source
  // text the guard below actually sees.
  //   1. An apostrophe inside a comment ("the browser's cache") opened a fake
  //      single-quoted string that never closed, so from that line to EOF every
  //      real comment survived — two doc comments in src/pdb.js and
  //      src/forcefield.js then read as "executing stub" hits. The quote char
  //      is now remembered per string so it can only be closed by itself.
  //   2. `${` in a template put the scanner into code mode for the rest of the
  //      FILE, with no way back. A stack of brace depths finds the closing `}`.
  //   3. A `/` in code is either a regex opener or a DIVISION operator, and a
  //      regex literal may contain quote characters that are not string
  //      delimiters. The standard "what precedes the slash" heuristic is used;
  //      it has to recognise an identifier or number as an OPERAND, not just
  //      `)]}`. Getting that wrong in the permissive direction turned
  //      `integ.time / 1000` inside a template in src/controllers/tick.js into
  //      a regex that ran to end-of-line, so the following `catch` was read as
  //      pattern text and one catch looked bare when it was not.
  let state = "code"; // code | line | block | str | tmpl | regex
  let quote = ""; // the delimiter that opened `str`
  const tmplStack = [];
  let depth = 0;
  let inClass = false; // inside a regex character class [...] — `/` does not end it
  /** True when a `/` in code position starts a regex rather than dividing. */
  const regexAllowed = () => {
    for (let k = out.length - 1; k >= 0; k--) {
      const ch = out[k];
      if (ch === " " || ch === "\t" || ch === "\n") continue;
      // After an operand (identifier char, digit, `)`, `]`, `}`) a `/` divides.
      // After anything else — an operator, a keyword, `(`, `,`, `=` — it opens
      // a regex literal.
      if (/[A-Za-z0-9_$)\]}]/.test(ch)) return false;
      return true;
    }
    return true; // start of file
  };
  while (i < n) {
    const c = src[i], d = src[i + 1];
    if (state === "code") {
      if (c === "/" && d === "/") { state = "line"; out += "  "; i += 2; continue; }
      if (c === "/" && d === "*") { state = "block"; out += "  "; i += 2; continue; }
      if (c === '"' || c === "'") { state = "str"; quote = c; out += c; i++; continue; }
      if (c === "`") { state = "tmpl"; tmplStack.push(-1); out += c; i++; continue; }
      if (c === "/" && regexAllowed()) { state = "regex"; inClass = false; out += c; i++; continue; }
      if (c === "{") depth++;
      if (c === "}") {
        depth--;
        // Return to the enclosing template WITHOUT popping the stack: the stack
        // is popped by the closing backtick, not by an interpolation. Popping
        // here made the SECOND `${...}` in a template like
        // `parseCa: … resSeq=${resSeq} x=${x} y=${y} …` find an empty stack,
        // stay in code mode, and then treat every comment from there to EOF as
        // source — which is how two doc comments in src/pdb.js and
        // src/forcefield.js first read as "executing stub" hits.
        if (tmplStack.length && tmplStack[tmplStack.length - 1] === depth) {
          state = "tmpl";
          out += c;
          i++;
          continue;
        }
      }
      out += c; i++; continue;
    }
    if (state === "line") {
      if (c === "\n") { state = "code"; out += "\n"; } else out += " ";
      i++; continue;
    }
    if (state === "block") {
      if (c === "*" && d === "/") { state = "code"; out += "  "; i += 2; continue; }
      out += c === "\n" ? "\n" : " ";
      i++; continue;
    }
    if (state === "str") {
      out += c;
      if (c === "\\") { out += d ?? ""; i += 2; continue; }
      if (c === quote) state = "code";
      i++; continue;
    }
    if (state === "regex") {
      out += c;
      if (c === "\\") { out += d ?? ""; i += 2; continue; }
      if (c === "[") inClass = true;
      else if (c === "]") inClass = false;
      else if (c === "/" && !inClass) state = "code";
      else if (c === "\n") state = "code"; // unterminated: do not run away
      i++; continue;
    }
    // tmpl
    out += c;
    if (c === "\\") { out += d ?? ""; i += 2; continue; }
    if (c === "`") { state = "code"; tmplStack.pop(); i++; continue; }
    if (c === "$" && d === "{") {
      state = "code";
      // Store the depth the closing `}` will restore (the `}` branch below
      // decrements BEFORE comparing), then count the `{` like any other brace so
      // nested braces inside the interpolation balance out.
      tmplStack[tmplStack.length - 1] = depth;
      depth++;
      out += "${";
      i += 2;
      continue;
    }
    i++; continue;
  }
  return out;
}

/** Index of the matching close for the bracket at `open`, or -1. */
function matchBracket(s, open, o, c) {
  let depth = 0;
  for (let i = open; i < s.length; i++) {
    if (s[i] === o) depth++;
    else if (s[i] === c) { depth--; if (depth === 0) return i; }
  }
  return -1;
}

/** Every block-form catch: {line, body, hasThrow} over comment-stripped source. */
function findCatches(stripped) {
  const out = [];
  const re = /\bcatch\b/g;
  let m;
  while ((m = re.exec(stripped)) !== null) {
    let i = m.index + 5;
    while (stripped[i] === " ") i++;
    if (stripped[i] === "(") {
      const e = matchBracket(stripped, i, "(", ")");
      if (e < 0) continue;
      i = e + 1;
    }
    while (stripped[i] === " ") i++;
    if (stripped[i] !== "{") continue; // `catch` used as an identifier
    const end = matchBracket(stripped, i, "{", "}");
    if (end < 0) continue;
    out.push({
      line: stripped.slice(0, m.index).split("\n").length,
      body: stripped.slice(i + 1, end).trim(),
      hasThrow: /\bthrow\b/.test(stripped.slice(i + 1, end)),
      index: m.index,
    });
    re.lastIndex = m.index + 5;
  }
  return out;
}

const files = srcModules();
const corpus = files.map((f) => ({ f, rel: path.relative(ROOT, f), src: fs.readFileSync(f, "utf8") }));

// ---- 0. VACUITY: the walk must actually see the tree (2) ----
assert(files.length >= MIN_FILES,
  `vacuity: walked ${files.length} src modules (need >= ${MIN_FILES})`);
let totalCatches = 0;
const bare = [];
for (const { rel, src } of corpus) {
  const stripped = stripComments(src);
  for (const c of findCatches(stripped)) {
    totalCatches++;
    const handled = c.body.length > 0 || c.hasThrow;
    if (!handled) bare.push(`${rel}:${c.line}`);
  }
}
assert(totalCatches >= MIN_CATCHES,
  `vacuity: found ${totalCatches} catch blocks in src/ (need >= ${MIN_CATCHES})`);

// ---- 1. GUARD: no bare catch anywhere in src/ (1 + 1 per offender) ----
// src/errors.js is exempted BY PATH and only for its own terminal guards: a
// recorder that called recordError on its own failure would recurse. The
// exemption is narrow (path, not pattern) and asserted non-empty below, so a
// future bare catch in errors.js still turns this red — the exemption only ever
// *narrows*, it can never hide a new offender.
const RECORDER = "src/errors.js";
const bareOutside = bare.filter((b) => !b.startsWith(`${RECORDER}:`));
const bareInRecorder = bare.filter((b) => b.startsWith(`${RECORDER}:`));
assert(bareOutside.length === 0,
  bareOutside.length === 0
    ? `no bare catch outside ${RECORDER} — all ${totalCatches - bareInRecorder.length} catches there rethrow, recover, record, or document`
    : `bare catch (empty body, no throw) at ${bareOutside.length} site(s): ${bareOutside.join(", ")}`);
// The bound was `<= 6` with no measurement behind it, so the exemption had two
// slots of silent headroom: two more bare catches could be planted in
// src/errors.js and this would still pass. Measured 2026-10-03: there are 4, and
// all 4 are annotated (asserted immediately below). Pinning the ceiling to the
// measured count closes those slots — the exemption can now only ever shrink,
// which is the property the comment above already claimed it had.
assert(bareInRecorder.length > 0 && bareInRecorder.length <= 4,
  `${RECORDER}'s terminal guards are the only empty bodies, and there are no more of them ` +
  `than the 4 this ceiling was pinned to (${bareInRecorder.length}: ${bareInRecorder.join(", ")})`);
// …and each one must SAY why a terminal guard is the right answer there, so the
// exemption is a documented decision rather than a quiet hole. Counted on the
// raw source, since the comment is exactly what the scanner blanks out.
const unannotated = bareInRecorder.filter((b) => {
  const ln = Number(b.split(":")[1]);
  const line = fs.readFileSync(path.join(ROOT, RECORDER), "utf-8").split("\n")[ln - 1] ?? "";
  return !/TERMINAL GUARD/.test(line);
});
assert(unannotated.length === 0,
  `every terminal guard in ${RECORDER} carries a TERMINAL GUARD annotation (${bareInRecorder.length - unannotated.length}/${bareInRecorder.length} annotated)`);

// ---- 2. THE RECORDER RECORDS AND COUNTS (9) ----
resetErrors();
assert(errorSummary().total === 0, "resetErrors(): counters start at zero");
assert(errorCounterText() === "", "clean recorder renders an empty status string");
assert(errorHudPrefix() === "", "clean recorder contributes nothing to #hud");

recordError(new Error("force eval exploded"), "ff.compute@tick");
recordError(new Error("force eval exploded"), "ff.compute@tick");
let s = errorSummary();
assert(s.total === 2 && s.errors === 2, `repeat of the same error counts twice (total ${s.total})`);
assert(s.slots === 1, `repeat of the same error occupies ONE slot (slots ${s.slots})`);
assert(s.last.count === 2, `the slot carries the occurrence count (${s.last?.count})`);
assert(s.last.message === "force eval exploded" && s.last.context === "ff.compute@tick",
  "the slot keeps message + context (a number alone cannot be acted on)");
assert(Number.isFinite(s.last.timestamp), `the slot is timestamped (${s.last?.timestamp})`);

ignore(new TypeError("x is not a function"), "headless", "ui ref null under Node");
s = errorSummary();
assert(s.errors === 2 && s.ignored === 1 && s.total === 3,
  `ignore() is a separate kind: errors ${s.errors} / ignored ${s.ignored} / total ${s.total}`);
assert(s.last.why === "ui ref null under Node", "a documented no-op keeps its justification on the slot");

// Ring bound: MAX_ERROR_SLOTS distinct keys must not grow the map without bound.
for (let i = 0; i < MAX_ERROR_SLOTS + 20; i++) recordError(`synthetic ${i}`, `flood@${i}`);
s = errorSummary();
assert(s.slots <= MAX_ERROR_SLOTS,
  `ring is bounded: ${MAX_ERROR_SLOTS + 20} distinct errors occupy ${s.slots} slots (cap ${MAX_ERROR_SLOTS})`);
assert(s.total === 3 + MAX_ERROR_SLOTS + 20,
  `eviction loses detail, never the total (${s.total} occurrences still counted)`);

// ---- 3. A FAILING OPERATION SURFACES THROUGH THE STATUS PATH (7) ----
// Driven through a real installGlobalErrorHandlers on a fake window, so the
// assertion covers the wiring, not just the counter arithmetic.
resetErrors();
const painted = [];
setErrorSink((t) => painted.push(t));
const fakeWin = {
  listeners: new Map(),
  addEventListener(type, fn) { (this.listeners.get(type) ?? this.listeners.set(type, []).get(type)).push(fn); },
  fire(type, ev) { for (const fn of this.listeners.get(type) ?? []) fn(ev); },
};
const installed = installGlobalErrorHandlers(fakeWin);
assert(installed === true, "installGlobalErrorHandlers() installs on a fresh target");
assert(installGlobalErrorHandlers(fakeWin) === false, "install is idempotent (no double reporting)");

// The deliberately failing operation: a force evaluation that throws.
try { throw new Error("force eval exploded: bad sigma"); } catch (e) { recordError(e, "ff.compute@step"); }
s = errorSummary();
assert(s.errors === 1, `a failing operation increments the counter (errors ${s.errors})`);
assert(painted.length === 1 && /Errors: 1/.test(painted[painted.length - 1]),
  `the top-bar surface is painted with the failure: "${painted[painted.length - 1]}"`);

// Uncaught error + unhandled rejection must land on the same counter.
fakeWin.fire("error", { error: new Error("uncaught boom"), message: "uncaught boom" });
fakeWin.fire("unhandledrejection", { reason: new Error("rejected") });
s = errorSummary();
assert(s.errors === 3 && s.total === 3,
  `global handlers feed the SAME counter (${s.errors} errors: forced + uncaught + unhandled)`);

// The pre-module net: an error before src/errors.js loaded (AUDIT_REPORT §2.1).
const boot = { __appErrorQueue: [{ message: "SyntaxError: bad import", context: "pre-module:src/main.js:5" }] };
const drained = drainBootstrapQueue(boot);
assert(drained === 1 && errorSummary().errors === 4,
  `a pre-module failure is drained, not lost (${drained} drained, ${errorSummary().errors} total)`);

// The rendered strings must match what the UI will actually show.
const bar = errorCounterText();
assert(/^Errors: \d+ ⚠/.test(bar), `top-bar text has the house shape: "${bar}"`);
const prefix = errorHudPrefix();
assert(prefix.startsWith("⚠ ") && prefix.endsWith("  ·  "),
  `#hud prefix matches the NON-FINITE ENERGY warning's shape (${JSON.stringify(prefix.slice(0, 56))}…)`);
assert(/ERRORS RECORDED/.test(prefix) || /ERROR RECORDED/.test(prefix),
  `#hud prefix names the failure ("${prefix.slice(0, 40)}…")`);
assert(errorReport().split("\n").length === errorSummary().slots,
  "errorReport() emits one line per slot, newest first");

// A suppressed-only state must not claim an error.
resetErrors();
ignore(new Error("cosmetic"), "paint@caption", "caption only");
assert(!/ERRORS RECORDED/.test(errorHudPrefix()) && /DOCUMENTED NO-OP/.test(errorHudPrefix()),
  "a documented no-op alone is reported as such, never as an error");
assert(/\(\+1 n\/o\)/.test(errorCounterText()),
  `the top bar shows the no-op count as secondary: "${errorCounterText()}"`);
resetErrors();
setErrorSink(null);

// The throttle must not make the bar LIE. A bare `if (too soon) return` drops
// the update, so if the last error of a session lands inside the window and
// nothing follows, the count on screen is permanently wrong. Caught by running
// the real recorder against a stub DOM in one burst.
const burstEls = new Map([["topErrors", { id: "topErrors", textContent: "", dataset: {}, title: "" }]]);
const prevDoc = globalThis.document;
globalThis.document = { getElementById: (id) => burstEls.get(id) ?? null };
resetErrors();
for (let i = 0; i < 5; i++) recordError(new Error(`burst ${i}`), `burst@${i}`);
const barEl = burstEls.get("topErrors");
assert(barEl.textContent === errorCounterText(),
  `the top bar shows the true count after a burst inside the throttle window (bar "${barEl.textContent}" vs summary "${errorCounterText()}")`);
assert(Number(barEl.dataset.errors) === errorSummary().errors,
  `the bar's data-errors counter agrees with the summary (${barEl.dataset.errors})`);
globalThis.document = prevDoc;
resetErrors();

// ---- 4. THE SURFACE EXISTS IN THE MARKUP (6) ----
const html = fs.readFileSync(path.join(ROOT, "index.html"), "utf-8");
const css = fs.readFileSync(path.join(ROOT, "css", "style.css"), "utf-8");
const uiSrc = fs.readFileSync(path.join(ROOT, "src", "ui.js"), "utf-8");
const tickSrc = fs.readFileSync(path.join(ROOT, "src", "controllers", "tick.js"), "utf-8");

assert(/id="topErrors"/.test(html), "index.html carries the #topErrors span");
// Inside the EXISTING top status bar: a new <details class="panel"> would move
// the Digit1-7 hotkey targets, which scripts/wikiskill_gate.js checks.
const inTopBar = /id="topStatus"[\s\S]*?id="topErrors"[\s\S]*?<\/div>/.test(html);
assert(inTopBar, "#topErrors lives INSIDE the existing #topStatus bar (no new panel → hotkey contract intact)");
assert(/\$?\("topErrors"\)|topErrors: \$\("topErrors"\)/.test(uiSrc),
  "src/ui.js wires topErrors into the ui map (so the DOM-id contract covers it)");
assert(/#topErrors:empty\s*\{\s*display:\s*none/.test(css),
  "css/style.css hides #topErrors when empty (zero visual change on a clean run)");
assert(/errorHudPrefix\(\)/.test(tickSrc),
  "src/controllers/tick.js puts the recorder's prefix on #hud, beside the NON-FINITE warning");
assert(/NON-FINITE ENERGY/.test(tickSrc) && /import \{[^}]*errorHudPrefix[^}]*\} from "\.\.\/errors\.js"/.test(tickSrc),
  "the #hud prefix sits on the SAME line as the NON-FINITE ENERGY precedent (one house style, not two)");
assert(!/window\.onerror\s*=.*getElementById\("hud"\)/.test(html),
  "index.html no longer clobbers #hud wholesale (the tick loop rewrites it 100 ms later)");

// ---- 4b. THE SCANNER IS PROVEN ON ITS OWN TRICKY INPUT (4) ----
// A guard that cannot see is not a guard. Each case below is a real shape that
// appears in src/ and that broke an earlier version of this scanner, which then
// reported a plausible-looking line while actually seeing nothing.
const probe = (code) => {
  const s = stripComments(code);
  return findCatches(s).filter((c) => !(c.body.length > 0 || c.hasThrow)).length;
};
assert(probe('try { a(); } catch (_) { /* x */ }') === 1,
  "scanner: a comment-only body is BARE (prose does not satisfy the rule)");
assert(probe('try { a(); } catch (e) { ignore(e, "c"); }') === 0,
  "scanner: a recorded body is handled");
assert(probe('try { a(); } catch (e) { throw e; }') === 0, "scanner: a rethrow is handled");
assert(probe('// try { a(); } catch (_) {}\ntry { b(); } catch (_) {}') === 1,
  "scanner: a catch inside a comment is not counted (only the real one is)");
assert(probe('const s = "the browser\'s cache"; // note\ntry { a(); } catch (_) {}') === 1,
  "scanner: an apostrophe in a COMMENT cannot open a string and swallow the rest of the file");
assert(probe('const r = /ab\'c/;\ntry { a(); } catch (_) {}') === 1,
  "scanner: a quote inside a REGEX literal is not a string delimiter");
assert(probe('const q = `a ${x / 2} b`;\ntry { a(); } catch (_) {}') === 1,
  "scanner: division inside a template interpolation is not a regex opener");
assert(probe('try { a(); } catch (_) {}') === 1, "scanner: the plain case is caught");

// ---- 5. NO EXECUTING STUB (1 + 1 per offender) ----
// evolve/evolve.mjs probe 5 counts a TODO/FIXME/not yet/placeholder/stub
// keyword on a line that EXECUTES as live debt: it is a value a user can see.
// Comment-line hits are honest scope notes (ROADMAP.md §1) and are left alone.
const STUB_RE = /\b(TODO|FIXME|not yet|placeholder|stub)\b/i;
const execStubs = [];
for (const { rel, src } of corpus) {
  const stripped = stripComments(src).split("\n");
  stripped.forEach((line, i) => {
    if (STUB_RE.test(line)) execStubs.push(`${rel}:${i + 1}`);
  });
}
assert(execStubs.length === 0,
  execStubs.length === 0
    ? "no executing TODO/stub/placeholder in src/ (comment-line scope notes are out of scope, per ROADMAP.md §1)"
    : `executing stub marker at ${execStubs.length} site(s): ${execStubs.join(", ")}`);

// ---- 6. THE RECORDER ITSELF CANNOT BECOME A FAILURE MODE (3) ----
assert(ERROR_PAINT_MIN_MS > 0, `the DOM repaint is throttled (${ERROR_PAINT_MIN_MS} ms), matching the 10 Hz #hud debounce`);
// A pathological context/message must not throw or leak unbounded text.
let ok = true;
try {
  recordError({ toString() { throw new Error("hostile toString"); } }, "hostile@toString");
  recordError(undefined, "undefined@err");
  recordError(null, null);
  recordError("x".repeat(5000), "huge@message");
} catch (_) { ok = false; }
assert(ok, "recordError survives a hostile value, a null context and a 5000-char message");
// A cross-realm Error (Web Worker, iframe) is NOT `instanceof Error` here. Using
// instanceof rendered "[object Object]" — the one message that tells the user
// nothing — and this project posts from Workers (src/worker-pool.js).
resetErrors();
const foreign = vm.runInNewContext("new Error('worker kernel failed')", {});
recordError(foreign, "worker.onerror");
assert(errorSummary().last.message === "worker kernel failed",
  `a cross-realm Error still yields its message (got "${errorSummary().last.message}")`);
recordError({ type: "unhandledrejection" }, "bareEvent");
assert(/unhandledrejection event/.test(errorSummary().last.message),
  `an ErrorEvent with no reason still renders something (got "${errorSummary().last.message}")`);
resetErrors();
// A DOM stub that throws on every property must not propagate either.
try {
  globalThis.document = { getElementById() { return { set textContent(_) { throw new Error("nope"); } }; } };
  recordError(new Error("dom hostile"), "hostile@dom");
  delete globalThis.document;
  ok = true;
} catch (_) { ok = false; }
assert(ok, "a throwing DOM does not propagate out of the recorder (it runs ON the failure path)");

console.log(`\n${passed} PASSED, ${failed} FAILED`);
process.exit(failed === 0 ? 0 : 1);