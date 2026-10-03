/**
 * physics/terms/registry.js — THE pluggable physics-term contract.
 *
 * WHAT THIS IS
 * ------------
 * Before this module, adding a force-field term meant editing an inlined
 * `U += …` summation sequence inside a class method in src/cg/compute.js or
 * src/heavy/energy.js. That is a single point of syntactic failure for every
 * consumer of the engine, and it puts the numeric summation order — the one
 * thing the bit-exact goldens exist to protect — one careless edit away from
 * moving. This module is the seam: a TERM is a declarative descriptor, and a
 * compute() assembly walks an ordered PLAN of descriptors instead of an
 * inlined statement sequence.
 *
 * THE INTERFACE (concrete; the whole thing)
 * ----------------------------------------
 *   id       string   "engine.name", unique, stable forever. Also the key the
 *                      UI, the docs and the analysis modules address a term by.
 *   engine   string   "cg" | "heavy" — which compute() assembly owns it.
 *   order    integer  position in the running-U accumulation, ascending. The
 *                      ORDER IS THE PHYSICS: U is a float64 sum, so reordering
 *                      two terms changes the last bits. See SUMMATION ORDER below.
 *   label    string   one line of human text, for the HUD and for docs.
 *   enabled  fn|null  `function () { return …this… }`. null ⇒ always on. A
 *                      term whose predicate is false is SKIPPED ENTIRELY: its
 *                      energy function is never called and U is not touched.
 *   energy   fn       `function (pos, f) { … ; return dU; }` with `this` bound
 *                      to the force field. MUST accumulate −∇U into `f` (the
 *                      caller's shared force buffer) and return its own ΔU as
 *                      ONE float64. Per-atom / per-term detail goes on the
 *                      force field (ff.repU etc.), not in the return value.
 *   reports  string[] OPTIONAL, documentation only: the names of the ff
 *                      accumulators this term publishes. Never read at runtime.
 *   report   fn       OPTIONAL, documentation only: `function () { return
 *                      this.repU; }` — how describeTerms() reads that term's
 *                      energy back for the UI. Never called in the hot path.
 *
 * That is five required fields and two optional ones. It is deliberately small:
 * the interface exists so a term can be IMPLEMENTED, not so every conceivable
 * term can be modelled.
 *
 * WHY `this` AND NOT AN EXPLICIT `ff` ARGUMENT
 * -------------------------------------------
 * Every force kernel in this repo already takes the field as `this`
 * (`_forceSnapshot()`, `_weakInteractions()`, `_nonBondedGrid()`), so binding it
 * keeps a third-party term indistinguishable from the built-in ones. The
 * dispatcher therefore calls `d.energy.call(ff, pos, f)`, and registerTerm()
 * REJECTS an arrow function (`energy.prototype === undefined`) rather than
 * letting it silently receive `this === undefined`. That check is the whole
 * reason the mistake cannot happen quietly.
 *
 * WHY THIS DOES NOT DEOPT THE HOT PATH (the hard requirement)
 * -----------------------------------------------------------
 * A naive registry — an array of closures, a megamorphic dispatch, a result
 * object per term per step — costs 2-5x on a 13 ms heavy compute. What makes
 * this one free is that it holds no per-step work at all:
 *
 *   • The PLAN is a plain frozen array, built ONCE per (engine, generation).
 *     The per-call cost is a single module-variable read plus an integer
 *     compare, done by the caller (src/cg/terms.js / src/heavy/terms.js), not
 *     by a Map lookup here.
 *   • `enabled === null` is the marker for "always on", so the common case is a
 *     null compare rather than a call.
 *   • Nothing is allocated per call: a term returns a number, never a record.
 *   • The only new cost is the dispatch itself — one `.call` per term, plus the
 *     loop that replaced the straight-line `U += …` sequence.
 *
 * MEASURED, not asserted (interleaved A/B, `bench/perf.js` and a block-timed
 * 9 x 4000-call CG benchmark; A/B blocks alternated so machine drift hits both
 * arms equally):
 *   • Heavy (1308 atoms, ~13.9 ms/compute): 8 pairs, 13.66-14.16 ms with the
 *     registry vs 13.76-14.16 ms without. Mean difference -0.3 %, inside the
 *     run-to-run spread. NO CHANGE.
 *   • CG (164 beads, ~0.0907 ms/compute block-timed): 15 pairs, mean 0.09118 ms
 *     with the registry vs 0.09065 ms without: +0.53 us = +0.58 %, positive in
 *     12 of 15 pairs. Reproducible, and NOT dispatch: three alternatives were
 *     implemented and measured, and all three landed on the same +0.55 us —
 *     (a) a `switch` on an integer `kind` with nine DIRECT monomorphic call
 *     sites, (b) an explicit `energy(ff, pos, f)` signature instead of `.call`,
 *     (c) unfrozen descriptors and plan array. So it is the cost of the loop
 *     itself, not of the indirection, and the simplest design is the one
 *     shipped. Against a benchmark whose own repeated-run spread is +-1 %, that
 *     is at the resolution limit of the harness.
 *
 * SUMMATION ORDER IS PART OF THE CONTRACT
 * ---------------------------------------
 * `U += d.energy(…)` happens once per term, in ascending `order`, at exactly
 * the position the term occupied in the inlined sequence. Two rules follow, and
 * both are load-bearing for the bit-exact goldens:
 *
 *   1. N terms the original code summed in ONE `+=` (`U += a + b`,
 *      `U += nb.lj + nb.elec + nb.gb + nb.hbond`) stay inside ONE descriptor,
 *      which returns that same expression. Splitting them into N descriptors
 *      would insert N−1 extra roundings — a physics change.
 *   2. The CG assembly has one hard-coded term (the inlined holo compression
 *      floor) that CANNOT be a descriptor: it accumulates U into the running
 *      total pair-by-pair, so giving it its own accumulator reassociates the
 *      sum. src/cg/compute.js keeps it inline and splits the plan around it.
 *      See HOLO_FLOOR_BARRIER_ORDER and docs/PHYSICS_TERMS.md.
 *
 * WHAT THIS IS NOT
 * ----------------
 * Not a general plugin system. There is no sandbox, no dependency injection, no
 * per-field registry (terms are process-global), no versioning of a term's
 * physics, and no way to register a term whose energy is not a single float64.
 * It is a fixed SHAPE — "an enabled predicate and a (pos, f) → ΔU function
 * that writes forces" — that happens to be data rather than an inlined
 * statement sequence, plus a validated registry and an introspection surface.
 * docs/PHYSICS_TERMS.md §"What this does not support" states the same list for
 * readers.
 *
 * Zero static imports, so it sits at the leaf of the dependency graph exactly
 * like src/physics/params.js and src/units.js. Node-importable, no DOM globals.
 */

/** Engines that own a compute() assembly and therefore a term plan. */
export const TERM_ENGINES = ["cg", "heavy"];

/**
 * Monotone counter bumped by every registry mutation. Exported as a LIVE
 * BINDING (`let`, not a function) so a plan cache in another module reads the
 * current value with a plain variable load and no call. This is the whole
 * invalidation mechanism — there is no other.
 */
let _generation = 0;
export let TERM_GENERATION = 0;

/** engine -> Map(id -> term). The live registry. */
const LIVE = new Map();
/** engine -> Map(id -> term). The immutable seed installBuiltins() installed. */
const BUILTIN = new Map();

function bump() {
  _generation += 1;
  TERM_GENERATION = _generation;
}

/** engine -> its Map, created on demand. Throws for an unknown engine. */
function slot(map, engine, what) {
  if (engine !== "cg" && engine !== "heavy") {
    throw new RangeError(
      `${what}: unknown engine ${JSON.stringify(engine)} — expected one of ${TERM_ENGINES.join(", ")}`,
    );
  }
  let m = map.get(engine);
  if (m === undefined) { m = new Map(); map.set(engine, m); }
  return m;
}

/**
 * `engine.name` — a lowercase engine segment, a dot, then a camelCase name.
 * Deliberately narrow: an id is a stable public name that the HUD, the docs
 * and the analysis modules all address, so it must not be a free-form string.
 */
const ID_RE = /^[a-z][a-z0-9_]*\.[A-Za-z][A-Za-z0-9_]*$/;

/**
 * Validate a descriptor against the interface documented above.
 *
 * The BUILT-IN terms go through this too (installBuiltins calls it), so the
 * registry is not a contract that only third parties are held to — a malformed
 * built-in descriptor is a load-time throw, not a silently mis-ordered term.
 *
 * @param {object} t
 * @returns {object} the same object, normalized (`enabled` → null, unit set)
 * @throws {TypeError|RangeError} with a message naming the offending field
 */
export function validateTerm(t) {
  if (t === null || typeof t !== "object" || Array.isArray(t)) {
    throw new TypeError("registerTerm: a term must be a plain object");
  }
  const { id, engine, order, label, energy, enabled } = t;
  if (typeof id !== "string" || !ID_RE.test(id)) {
    throw new TypeError(
      `registerTerm: id must match ${ID_RE} (got ${JSON.stringify(id)}) — ` +
      `use "engine.name", e.g. "cg.myWell" or "heavy.myWell"`,
    );
  }
  if (typeof engine !== "string") {
    throw new TypeError(`registerTerm(${id}): engine must be a string`);
  }
  if (typeof order !== "number" || !Number.isInteger(order) || !Number.isFinite(order)) {
    throw new TypeError(
      `registerTerm(${id}): order must be an integer — it is the position in the ` +
      `running-U accumulation, and a fractional position has no meaning`,
    );
  }
  if (typeof label !== "string" || label.length === 0) {
    throw new TypeError(`registerTerm(${id}): label must be a non-empty string (UI/docs text)`);
  }
  if (typeof energy !== "function") {
    throw new TypeError(`registerTerm(${id}): energy must be a function (pos, f) -> dU`);
  }
  // The `this === force field` contract. An arrow function has no `prototype`,
  // so this rejects the exact mistake that would otherwise fail silently with
  // `this === undefined` at step 1 of a 100-step trajectory.
  if (energy.prototype === undefined) {
    throw new TypeError(
      `registerTerm(${id}): energy must be a \`function (pos, f) {}\` — not an arrow ` +
      `function. The force field arrives as \`this\`; an arrow would silently ` +
      `receive undefined.`,
    );
  }
  if (enabled !== null && enabled !== undefined && typeof enabled !== "function") {
    throw new TypeError(
      `registerTerm(${id}): enabled must be null (always on) or a function () -> boolean`,
    );
  }
  if (t.reports !== undefined && !Array.isArray(t.reports)) {
    throw new TypeError(`registerTerm(${id}): reports, when present, must be an array of names`);
  }
  if (t.report !== undefined && typeof t.report !== "function") {
    throw new TypeError(`registerTerm(${id}): report, when present, must be a function () -> number`);
  }
  return {
    id, engine, order, label, energy,
    enabled: enabled ?? null,
    reports: t.reports ?? null,
    report: t.report ?? null,
    unit: typeof t.unit === "string" ? t.unit : "kcal/mol",
  };
}

/**
 * Seed an engine's built-in terms. Idempotent per engine: a second call
 * REPLACES the seed (and drops any runtime overrides), because the built-in
 * list is source code and source code wins. Called at module load by
 * src/cg/terms.js and src/heavy/terms.js — never by a caller.
 * @param {string} engine
 * @param {object[]} terms
 * @returns {object[]} the validated, order-sorted built-ins
 */
export function installBuiltins(engine, terms) {
  if (!Array.isArray(terms)) throw new TypeError("installBuiltins: terms must be an array");
  const seed = new Map();
  for (const t of terms) {
    const v = validateTerm(t);
    if (v.engine !== engine) {
      throw new TypeError(
        `installBuiltins(${engine}): term ${v.id} declares engine ${JSON.stringify(v.engine)}`,
      );
    }
    if (seed.has(v.id)) throw new TypeError(`installBuiltins(${engine}): duplicate term id ${v.id}`);
    seed.set(v.id, v);
  }
  BUILTIN.set(engine, seed);
  slot(LIVE, engine, "installBuiltins");
  LIVE.set(engine, new Map(seed));
  bump();
  return sorted(seed);
}

/**
 * Register (or replace) one term. Re-registering an existing id is the
 * supported way to override a built-in — e.g. an analysis module that wants a
 * term's energy zeroed, or a site-specific parameterisation that reuses a
 * built-in's slot with different constants.
 *
 * The registry is PROCESS-GLOBAL: there is no per-force-field registry. A term
 * registered here is live for every ForceField / HeavyForceField in the
 * process from the next compute() on. Registering after a compute() is fine
 * (the plan cache invalidates on TERM_GENERATION); unregistering in a `finally`
 * is what a test should do.
 *
 * @param {object} term  a descriptor, validated by validateTerm()
 * @returns {object} the normalized term actually stored
 */
export function registerTerm(term) {
  const v = validateTerm(term);
  const m = slot(LIVE, v.engine, "registerTerm");
  m.set(v.id, v);
  bump();
  return v;
}

/**
 * Remove a term by id — including a built-in, which is the supported way to
 * turn a term off for a whole process (an analysis module perturbing the model
 * does not have to reimplement the term to leave it out).
 * @param {string} id
 * @returns {boolean} true if a term was removed
 */
export function unregisterTerm(id) {
  for (const engine of TERM_ENGINES) {
    const m = LIVE.get(engine);
    if (m !== undefined && m.delete(id)) { bump(); return true; }
  }
  return false;
}

/**
 * Drop every runtime override for an engine and restore its built-in set.
 * @param {string} engine
 */
export function resetTerms(engine) {
  const seed = BUILTIN.get(engine);
  if (seed === undefined) { slot(LIVE, engine, "resetTerms"); bump(); return; }
  LIVE.set(engine, new Map(seed));
  bump();
}

/** One descriptor by id, or undefined. Searches both engines. */
export function getTerm(id) {
  for (const engine of TERM_ENGINES) {
    const m = LIVE.get(engine);
    const t = m === undefined ? undefined : m.get(id);
    if (t !== undefined) return t;
  }
  return undefined;
}

/**
 * The engine's terms, ordered by `order` ascending. THIS is the plan the
 * assemblies walk. Rebuilt per call — so callers must memoise it against
 * TERM_GENERATION rather than calling this inside compute(). The memo lives in
 * src/cg/terms.js and src/heavy/terms.js for exactly that reason.
 *
 * Ordering is total and deterministic: `order` ascending, then `id` ascending
 * as a tiebreak, so two terms that claim the same slot get a defined order
 * instead of depending on registration sequence.
 *
 * @param {string} engine
 * @returns {ReadonlyArray<object>} a frozen array of frozen descriptors
 */
export function listTerms(engine) {
  return sorted(slot(LIVE, engine, "listTerms"));
}

function sorted(m) {
  const out = [...m.values()];
  out.sort((a, b) => (a.order - b.order) || (a.id < b.id ? -1 : a.id > b.id ? 1 : 0));
  for (const d of out) Object.freeze(d);
  return Object.freeze(out);
}

/**
 * Introspection for the HUD, the docs and the analysis modules. NOT a hot-path
 * function: it allocates one record per term and calls each term's `report`,
 * so it belongs in a UI refresh or a test, never in compute().
 *
 * @param {string} engine
 * @param {object} ff  a force field of that engine, for `enabled` / `report`
 * @returns {Array<{id:string,label:string,order:number,unit:string,
 *   active:boolean,energy:(number|null),reports:(string[]|null)}>}
 */
export function describeTerms(engine, ff) {
  const out = [];
  for (const d of listTerms(engine)) {
    const active = d.enabled === null ? true : d.enabled.call(ff) === true;
    let energy = null;
    if (active && d.report !== null) {
      const v = d.report.call(ff);
      energy = typeof v === "number" && Number.isFinite(v) ? v : null;
    }
    out.push({
      id: d.id, label: d.label, order: d.order, unit: d.unit,
      active, energy, reports: d.reports,
    });
  }
  return out;
}
