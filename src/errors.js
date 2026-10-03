/**
 * errors.js — the ONE place a swallowed failure becomes countable and visible.
 *
 * THE DEFECT THIS EXISTS TO KILL
 * -----------------------------
 * Measured before this module: 181 `catch` blocks under src/, of which 92 had
 * an empty body and 56 of those had no annotation either — a bare `catch (_) {}`
 * that swallows everything and continues. On a scientific simulator that is not
 * a style problem, it is a correctness problem: a user loads a structure, a step
 * of the pipeline fails, the code swallows it and keeps drawing, and the user is
 * left watching a plausible animation that is not the simulation they think they
 * are running. That is indistinguishable from a physics bug and it destroys trust
 * in every number the tool reports.
 *
 * THE RULE THIS ENFORCES
 * ----------------------
 * Every `catch` in src/ must now be handled in exactly one of four ways:
 *   (a) rethrown                                     — nothing here does this today
 *   (b) a specific, user-visible recovery            — the call sites that already
 *                                                       had one (formatInputError
 *                                                       into #hud, setThermoCaption,
 *                                                       the RESPA BAOAB fallback)
 *   (c) RECORDED here, counted, and surfaced          — `recordError()`
 *   (d) an annotated, documented no-op that is STILL  — `ignore()`
 *       counted, so "it never fires" is a measurement
 * `catch (_) {}` with no body is banned outright; tests/test_error_surfacing.js
 * plants one and reads red.
 *
 * WHY (d) IS COUNTED AND NOT JUST ANNOTATED
 * ------------------------------------------
 * Most of the empty catches are legitimately defensive (a browser API absent in
 * Node, a DOM ref captured as null by a headless import, a GPU buffer that is
 * already gone). Rethrowing those would abort module initialisation or the
 * requestAnimationFrame loop, so they keep the catch. But "defensive" is a claim
 * about what happens in the field, and until now nothing measured it. Routing
 * them through `ignore()` keeps the no-op AND makes the first time one fires a
 * number on the top status bar instead of a silence.
 *
 * TWO KINDS, AND WHY THEY ARE SEPARATE
 * ------------------------------------
 *   recordError()  kind "error"     — the user must know. Shown on the top bar
 *                                    and prefixed onto #hud next to the
 *                                    NON-FINITE ENERGY warning (same house style:
 *                                    src/controllers/tick.js:246).
 *   ignore()       kind "ignored"   — a documented no-op. Counted and shown as a
 *                                    secondary figure so it can never be mistaken
 *                                    for "nothing happened", but it does not
 *                                    shout over a real failure.
 *
 * WHY NO IMPORTS
 * --------------
 * Zero deps, and specifically no import of ./ui.js: ui.js is imported by nearly
 * every controller, so importing ui back from here would close a cycle (the
 * import graph is asserted acyclic by tests/test_module_size.js). The one DOM
 * node this module paints is looked up lazily by id, which also keeps the module
 * importable in Node with no document at all.
 *
 * HEADLESS-SAFE AND THROW-PROOF BY CONSTRUCTION
 * ---------------------------------------------
 * `recordError` runs *on* the failure path. If it threw, it would replace the
 * error the caller was already handling with a second, less useful one. Every
 * entry point here is wrapped so the recorder cannot throw, and the surface
 * lookup can only ever miss.
 *
 * See docs/ERRORS.md for the classification table of all 181 catches.
 */

/**
 * Distinct `context|message` slots retained in memory. Beyond this the oldest
 * slot is evicted, so a run that fails on every frame in a dozen places cannot
 * grow this map without bound. Occurrence counts inside a slot are NOT capped,
 * so eviction loses the detail, never the total.
 */
export const MAX_ERROR_SLOTS = 64;

/**
 * Minimum milliseconds between DOM repaints of the error surface. The recorder
 * counts on every call but repaints the top bar at most this often, for the same
 * reason src/controllers/tick.js:227 debounces #hud to 10 Hz: a per-frame caller
 * (three of them in the idle branch of tick) would otherwise cause a DOM write
 * per frame. The FIRST error always paints immediately (the clock starts at
 * -Infinity), so the user is never left waiting on a throttle.
 */
export const ERROR_PAINT_MIN_MS = 250;

/** Message truncation for the summary strings (the full text stays in slots). */
export const ERROR_MESSAGE_MAX = 200;

/** key -> { message, context, timestamp, count, kind, seq } in insertion order. */
const slots = new Map();
/** Total occurrences of every kind (never decremented, never evicted). */
const counts = { error: 0, ignored: 0 };
/** Monotonic sequence: two records in the same millisecond still order stably. */
let seq = 0;
let lastPaint = { text: "", at: -Infinity };
/** Pending trailing repaint (see paint()); null when none is scheduled. */
let pending = null;
let sink = null;

/** Test seam: a function called with the same text the DOM surface receives. */
function setErrorSink(fn) {
  sink = typeof fn === "function" ? fn : null;
}

function clip(s) {
  const t = String(s == null ? "" : s).replace(/\s+/g, " ").trim();
  return t.length > ERROR_MESSAGE_MAX ? `${t.slice(0, ERROR_MESSAGE_MAX - 1)}…` : t;
}

function now() {
  return typeof Date !== "undefined" && typeof Date.now === "function" ? Date.now() : 0;
}

/**
 * Extract a displayable message from a caught value.
 *
 * NOT `err instanceof Error`: an Error thrown inside a Web Worker or an iframe
 * belongs to that realm, so `instanceof Error` is FALSE for it in this one — and
 * the recorder rendered `[object Object]`, which is the one message that tells
 * the user nothing at all. This project posts from Workers (src/worker-pool.js)
 * and can host an iframe, so duck-typing is the correct check, not a nicety.
 * @param {unknown} err the caught value
 * @returns {string} a message, never empty
 */
function messageOf(err) {
  if (err && typeof err === "object" && typeof err.message === "string" && err.message) {
    return err.message;
  }
  if (err && typeof err === "object" && typeof err.type === "string") {
    return `(${err.type} event)`; // ErrorEvent / unhandledrejection with no reason
  }
  const s = clip(err);
  return s || String(err);
}

/**
 * The one write path for every recorded failure. Never throws.
 * @param {unknown} err the caught value
 * @param {string} context where it happened ("parseLigands@buildSystem")
 * @param {"error"|"ignored"} kind "error" = user must know; "ignored" = documented no-op
 * @param {string} [why] justification text for the "ignored" kind
 * @returns {{message:string, context:string, timestamp:number, count:number, kind:string, seq:number}|null}
 */
function record(err, context, kind, why = "") {
  try {
    const message = messageOf(err);
    const ctx = clip(context) || "unknown";
    const key = `${ctx} :: ${message}`;
    const existing = slots.get(key);
    if (existing) {
      existing.count += 1;
      existing.timestamp = now();
    } else {
      if (slots.size >= MAX_ERROR_SLOTS) {
        const oldest = slots.keys().next();
        if (!oldest.done) slots.delete(oldest.value);
      }
      slots.set(key, { message, context: ctx, timestamp: now(), count: 1, kind, seq: ++seq, ...(why ? { why: clip(why) } : {}) });
    }
    counts[kind] += 1;
    paint();
    return slots.get(key);
  } catch (_) {
    // A recorder that throws while handling an error would mask it. Give up
    // quietly: the counters are best-effort telemetry, never a control path.
    return null;
  }
}

/**
 * Record a failure the user must know about. Surfaces on the top status bar and
 * (via errorHudPrefix) on #hud. Never throws.
 * @param {unknown} err the caught value
 * @param {string} [context] where it happened
 * @returns {object|null} the slot, or null if the recorder gave up
 */
export function recordError(err, context = "unknown") {
  return record(err, context, "error");
}

/**
 * Record a DOCUMENTED no-op: the catch stays (rethrow would break startup or
 * the render loop) and the justification is kept on the slot, but the occurrence
 * is counted and shown, so "this never fires" stops being an assumption.
 * Never throws.
 * @param {unknown} err the caught value
 * @param {string} [context] where it happened
 * @param {string} [why] why swallowing is safe here
 * @returns {object|null} the slot, or null if the recorder gave up
 */
export function ignore(err, context = "unknown", why = "") {
  return record(err, context, "ignored", why);
}

/**
 * Aggregate view for the UI and for tests. A plain snapshot: mutating it cannot
 * corrupt the recorder.
 * @returns {{total:number, errors:number, ignored:number, slots:number, last:object|null}}
 */
export function errorSummary() {
  let last = null;
  for (const s of slots.values()) if (!last || s.seq > last.seq) last = s;
  return {
    total: counts.error + counts.ignored,
    errors: counts.error,
    ignored: counts.ignored,
    slots: slots.size,
    last: last ? { ...last } : null,
  };
}

/**
 * Counter text for the top status bar, or "" when nothing has ever failed (so
 * the span stays empty and CSS hides it — no visual change on a clean run).
 * @returns {string} e.g. "Errors: 3 ⚠ (+12 n/o)" or ""
 */
export function errorCounterText() {
  const s = errorSummary();
  if (s.total === 0) return "";
  const noop = s.ignored > 0 ? ` (+${s.ignored} n/o)` : "";
  return `Errors: ${s.errors} ⚠${noop}`;
}

/**
 * #hud prefix, matching the existing NON-FINITE ENERGY warning verbatim in
 * shape (see src/controllers/tick.js:246): "⚠ <what> — <what to do>.  ·  ".
 * Returns "" when clean, so the caller can concatenate it unconditionally.
 * @returns {string} prefix ending in "  ·  " or ""
 */
export function errorHudPrefix() {
  const s = errorSummary();
  if (s.total === 0) return "";
  const one = s.total === 1 ? "" : "S";
  const what = s.errors > 0
    ? `${s.errors} ERROR${one} RECORDED`
    : `${s.ignored} DOCUMENTED NO-OP${one}`;
  const last = s.last ? `last: ${s.last.context}: ${clip(s.last.message).slice(0, 80)}` : "";
  return `⚠ ${what}${last ? ` — ${last}` : ""}  ·  `;
}

/**
 * Multi-line detail for the console (and for a copy/paste bug report).
 * @returns {string} one line per slot, newest first
 */
export function errorReport() {
  const rows = [...slots.values()].sort((a, b) => b.seq - a.seq);
  if (rows.length === 0) return "no recorded errors";
  return rows
    .map((r) => `[${r.kind}] ${r.context}: ${r.message} (x${r.count}${r.why ? `; why swallowed: ${r.why}` : ""})`)
    .join("\n");
}

/**
 * Repaint the error surface. The DOM node is resolved lazily by id so this
 * module stays importable under Node. Throttled by ERROR_PAINT_MIN_MS and a
 * no-op when the text has not changed.
 * @returns {void}
 */
function paint() {
  try {
    const text = errorCounterText();
    const t = now();
    if (text === lastPaint.text && slots.size > 0) return;
    const due = t - lastPaint.at >= ERROR_PAINT_MIN_MS;
    if (!due && lastPaint.text !== "") {
      // THROTTLED — but a throttle that just drops the update is how the bar
      // ends up permanently stale: if the LAST error of a session lands inside
      // the window and nothing further arrives, the count the user is looking
      // at is simply wrong forever. So a trailing repaint is scheduled instead,
      // and the DOM text is written NOW with the current (possibly not-yet-
      // due) value rather than left at the old one. Proof of the old bug:
      // #topErrors read "Errors: 1 ⚠" while errorSummary() said 2.
      writeSurface(text);
      if (pending === null && typeof setTimeout === "function") {
        pending = setTimeout(() => { pending = null; paint(); }, ERROR_PAINT_MIN_MS);
      }
      return;
    }
    if (pending !== null && typeof clearTimeout === "function") { clearTimeout(pending); pending = null; }
    lastPaint = { text, at: t };
    writeSurface(text);
  } catch (_) { /* TERMINAL GUARD: the surface must never become a failure mode; recording here would recurse */ }
}

/**
 * The actual DOM write, split out so a throttled paint can still refresh the
 * visible text without disturbing the throttle bookkeeping. Never throws.
 * @param {string} text the bar text to render ("" clears it)
 * @returns {void}
 */
function writeSurface(text) {
  try {
    if (sink) sink(text);
    if (typeof document === "undefined" || !document.getElementById) return;
    const el = document.getElementById("topErrors");
    if (!el) return;
    el.textContent = text;
    const s = errorSummary();
    el.dataset.errors = String(s.errors);
    el.dataset.ignored = String(s.ignored);
    // The tooltip carries the detail the 11px bar has no room for.
    el.title = text ? errorReport() : "";
  } catch (_) { /* TERMINAL GUARD: the surface must never become a failure mode; recording here would recurse */ }
}

/**
 * Reset every counter and slot, and clear the surface. Test-only seam — nothing
 * in src/ calls this.
 * @returns {void}
 */
export function resetErrors() {
  slots.clear();
  counts.error = 0;
  counts.ignored = 0;
  seq = 0;
  lastPaint = { text: "", at: -Infinity };
  // A trailing repaint left over from before the reset would repaint a stale
  // count onto the bar, so cancel it.
  if (pending !== null && typeof clearTimeout === "function") { clearTimeout(pending); pending = null; }
  try {
    if (typeof document !== "undefined" && document.getElementById) {
      const el = document.getElementById("topErrors");
      if (el) { el.textContent = ""; el.title = ""; delete el.dataset.errors; delete el.dataset.ignored; }
    }
  } catch (_) { /* TERMINAL GUARD: test-only reset; a DOM that throws must not fail the suite */ }
}

/**
 * Drain what the pre-module bootstrap in index.html queued before this module
 * loaded. A classic script runs before the deferred module graph, so it is the
 * only thing that can catch a failure to PARSE src/main.js — the exact defect in
 * AUDIT_REPORT.md §2.1, where a SyntaxError showed the user a half-rendered page
 * and no message.
 * @param {object} win the window that carries __appErrorQueue
 * @returns {number} how many queued entries were recorded
 */
export function drainBootstrapQueue(win) {
  let n = 0;
  try {
    const q = win && win.__appErrorQueue;
    if (!Array.isArray(q)) return 0;
    while (q.length) {
      const e = q.shift();
      recordError(e, e?.context || "pre-module");
      n++;
    }
  } catch (_) { /* TERMINAL GUARD: pre-module net, best effort by definition; recursing would re-enter the recorder */ }
  return n;
}

/**
 * Route uncaught errors, unhandled rejections and ErrorEvents into the recorder
 * so they land on the same counter as everything else. Idempotent.
 *
 * This replaces three inline handlers in index.html that (i) reported every
 * window error TWICE (window.onerror and an "error" listener both fire for the
 * same event), (ii) overwrote #hud wholesale — the line the tick loop rewrites
 * 100 ms later, so the message lasted a tenth of a second — and (iii) counted
 * nothing, so "how many times did this fail" had no answer.
 *
 * @param {object} [win] defaults to globalThis
 * @returns {boolean} true when this call installed the handlers
 */
export function installGlobalErrorHandlers(win = (typeof globalThis !== "undefined" ? globalThis : undefined)) {
  if (!win || win.__appErrorHandlersInstalled) return false;
  try {
    win.__appErrorDrain = () => { drainBootstrapQueue(win); };
    win.onerror = (message, source, line, col, err) => {
      recordError(err || message, `window.onerror@${String(source || "inline").split("/").pop()}:${line || 0}:${col || 0}`);
    };
    if (typeof win.addEventListener === "function") {
      win.addEventListener("error", (e) => {
        recordError(e?.error || e?.message || "ErrorEvent with no error", "window.error");
      });
      win.addEventListener("unhandledrejection", (e) => {
        recordError(e?.reason ?? "unhandled rejection (no reason)", "window.unhandledrejection");
      });
    }
    win.__appErrorHandlersInstalled = true;
    drainBootstrapQueue(win);
    return true;
  } catch (_) {
    return false;
  }
}

export { setErrorSink };