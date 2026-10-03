/**
 * results-record.js — ONE small JSON file that carries a run's provenance, its
 * input identity, its computed values WITH their real uncertainties, its
 * run health, and the machine-readable statement of what it does not mean.
 *
 * THE USER PROBLEM
 * ----------------
 * Before this, an answer from this tool lived in three artefacts: a trajectory
 * file whose header carried a REMARK line, a PMF CSV, and the operator's own
 * notes. The parameters, the seed, the code version, which bytes were loaded,
 * the values, and the uncertainty were never in one place — so reproducing or
 * citing a number was guesswork. That is a defect for a project whose entire
 * thesis is honest, reproducible, interactive biophysics.
 *
 * WHAT IS DELIBERATELY NOT HERE
 * -----------------------------
 * No frames, no coordinates, no PDB text. This is a summary with provenance,
 * not a trajectory dump — the same line src/session.js draws at
 * `SESSION_MAX_BYTES`. Frames go to the recorder export, the numbers go here.
 * Results are capped at RESULTS_MAX_VALUES entries, one row each.
 *
 * THE THREE RULES THAT MAKE IT WORTH CITING
 * -----------------------------------------
 * 1. AN UNCERTAINTY IS NEVER ZERO. `uncertainty` is either an object with a
 *    finite POSITIVE value, or null with a machine-readable
 *    `uncertaintyReason`. This matters because the analysis modules DO return
 *    0 for "I could not compute this" — `computeThermodynamics` leaves `dH_se`
 *    at 0 when the block bootstrap has fewer than 2 blocks, and
 *    `p.sasa.se ?? 0` does the same. Writing that 0 into a record as an error
 *    bar would claim a precision nobody measured. `normalizeUncertainty()`
 *    converts it to null + reason instead.
 * 2. NO VALUE EXPORTS WITHOUT A TRUST BOUNDARY. Every result id must have a
 *    `VALUE_TRUST_BOUNDARY` entry in src/scope.js, or `resultsExportGate()`
 *    fails closed. "We forgot to say what ΔG means" can never ship as a
 *    clean-looking JSON file.
 * 3. A RUN WITH RECORDED ERRORS IS NOT EXPORTABLE AS IF IT WERE CLEAN.
 *    `src/errors.js` is the project's single failure recorder; the record
 *    carries its counts, and the gate refuses the export while
 *    `runHealth.errors > 0`. Documented no-ops (`ignore()`) are counted and
 *    carried too, but they do not block: they are the defensive catches, and
 *    the record says how many there were rather than pretending there were
 *    none.
 * 4. A ROW WITH NO VALUE IS ABSENT, NOT ZERO. `resultRow()` returns null when
 *    the value is not finite and the row is dropped, so a quantity the
 *    analysis could not compute (RMSIP with no elastic network, ΔG with no
 *    funnel) never appears as `0`. This is rule 1 applied to the value instead
 *    of the error bar, and it is why `coverage.included` is worth reading: it
 *    lists what the record actually holds.
 *
 * DOM-FREE AND NODE-IMPORTABLE
 * ----------------------------
 * No `document`, no `window`, no `fetch`, no `node:` imports — the only browser
 * API used is `TextEncoder` (in src/sha256.js), which is also a Node global, so
 * plain `node tests/test_results_record.js` exercises the same code the browser
 * runs. Downloads go through the existing `downloadText` pattern in the
 * controller, not from here.
 *
 * The honesty block is BUILT from src/scope.js at call time, and scope.js is
 * held to ROADMAP.md §1 / docs/LIMITATIONS.md by tests/test_results_record.js.
 * See src/scope.js's header for what is and is not derived at runtime.
 */

import { VERSION, BUILD_DATE } from "./version.js";
import { sha256Text } from "./sha256.js";
import { scopeBlock, VALUE_TRUST_BOUNDARY } from "./scope.js";
import { createInputError } from "./input_errors.js";

/** Results-record schema version (frozen contract — tests assert === 1). */
export const RESULTS_RECORD_VERSION = 1;

/** Refuse to serialize a record larger than this (bytes). A record is ~2-6 KiB. */
export const RESULTS_MAX_BYTES = 64 * 1024;

/** Refuse to load a record file larger than this (bytes). */
export const RESULTS_FILE_MAX_BYTES = 1024 * 1024;

/** Hard cap on rows in `results` — keeps the file small and the shape stable. */
export const RESULTS_MAX_VALUES = 24;

/**
 * The result ids the shipped UI wiring produces. Declared here (not in the
 * controller) so a test can assert every one of them has a trust-boundary
 * entry — the gate would refuse them at runtime, which is the right behaviour,
 * but a test failure points at the cause instead of the symptom.
 */
export const RESULTS_VALUE_IDS = [
  "b_factor_pearson_r",
  "b_factor_mean_sim_ang2",
  "rmsf_mean_ang",
  "var_top_k",
  "rmsip_pca_enm",
  "occupancy_bound_fraction",
  "contact_lifetime_mean_ps",
  "pocket_volume_mean_ang3",
  "pocket_volume_sd_ang3",
  "dg_bind",
];

/** Model modes accepted by the record (mirrors src/session.js). */
export const RESULTS_MODEL_MODES = ["cg", "heavy"];

/** Physics levels accepted by the record (mirrors src/session.js). */
export const RESULTS_PHYSICS_LEVELS = ["L0", "L1", "L2"];

/**
 * Strict numeric coercion: null/undefined/"" are NOT zero.
 *
 * `Number(null) === 0` and `Number("") === 0`, so the obvious
 * `Number.isFinite(Number(v))` test silently turns a missing seed into a
 * reproducible-looking `seed: 0`, which is the single most damaging thing this
 * module could do to a provenance field: the record would claim a seed for a run
 * that had none. Rejecting nullish and empty-string here is what keeps
 * "unseeded" distinguishable from "seed zero" all the way to the JSON.
 */
const finite = (v, dflt = null) => (
  v === null || v === undefined || v === "" || !Number.isFinite(Number(v)) ? dflt : Number(v)
);

/**
 * Normalize one uncertainty into the record's shape.
 *
 * Returns `{uncertainty, uncertaintyAvailable, uncertaintyReason}`. The ONE
 * rule: a zero, a NaN, an Infinity or an absent error bar becomes
 * `uncertainty: null` plus a reason that says WHY it is unknown — never a 0,
 * never a silent omission. `kind`/`method` name the estimator so a reader knows
 * whether it is a bootstrap, a jackknife or a hill-count convergence bound.
 *
 * @param {object|null|undefined} u `{value, kind, method, n}` or null
 * @param {string} fallbackReason sentence used when u carries no value
 * @param {number} [n] sample count when known
 */
function normalizeUncertainty(u, fallbackReason, n = null) {
  const v = u && typeof u === "object" ? finite(u.value, null) : null;
  if (v === null || v <= 0) {
    return {
      uncertainty: null,
      uncertaintyAvailable: false,
      uncertaintyReason: v === 0
        ? `the analysis code returned 0 for this error bar, which means "not computed", not "no error": ${fallbackReason}`
        : fallbackReason,
    };
  }
  return {
    uncertainty: {
      value: v,
      kind: typeof u.kind === "string" ? u.kind : "unspecified",
      method: typeof u.method === "string" ? u.method : "",
      ...(Number.isFinite(Number(n)) ? { n: Number(n) } : (Number.isFinite(Number(u.n)) ? { n: Number(u.n) } : {})),
    },
    uncertaintyAvailable: true,
    uncertaintyReason: null,
  };
}

/**
 * Build one `results[]` row.
 * @param {object} spec `{id, label, value, unit, uncertainty, noUncertaintyReason, n, dispersion}`
 * @returns {object|null} the row, or null when there is no value to report
 */
function resultRow(spec) {
  const v = finite(spec?.value, null);
  if (v === null) return null; // nothing computed — represented by absence, not by 0
  const id = String(spec.id ?? "");
  const boundary = VALUE_TRUST_BOUNDARY[id] ?? null;
  const u = normalizeUncertainty(spec.uncertainty, spec.noUncertaintyReason ?? "this quantity is a single deterministic statistic of the recorded frames; the analysis code computes no error bar for it", spec.n);
  const row = {
    id,
    label: String(spec.label ?? id),
    value: v,
    unit: String(spec.unit ?? ""),
    ...u,
    ...(Number.isFinite(Number(spec.n)) ? { n: Number(spec.n) } : {}),
  };
  // A dispersion (e.g. the SD of a per-frame track) is NOT an uncertainty on the
  // mean. Keeping it in its own field is what stops a track SD being read as an
  // error bar by whatever consumes the CSV.
  if (spec.dispersion && finite(spec.dispersion.value, null) !== null) {
    row.dispersion = {
      value: Number(spec.dispersion.value),
      kind: String(spec.dispersion.kind ?? "unspecified"),
      note: "a dispersion of the underlying sample, NOT an uncertainty on the mean above",
    };
  }
  row.meaning = boundary ? boundary.meaning : null;
  row.notA = boundary ? boundary.notA : null;
  row.notComputed = boundary ? boundary.notComputed.slice() : null;
  row.substitute = boundary ? boundary.substitute : null;
  row.documentedIn = boundary ? "src/scope.js VALUE_TRUST_BOUNDARY" : null;
  row.honestyMissing = !boundary;
  return row;
}

/**
 * Build a results record from a plain snapshot (pure; never throws on missing
 * fields, never touches the DOM, never stores frames or PDB text).
 *
 * @param {object} [snap]
 * @param {object} [snap.input] `{pdbId, fileName, origin, pdbText}`
 * @param {object} [snap.selection] `{chains, resFrom, resTo, includeLigand}`
 * @param {object} [snap.physics] the live force field / integrator values (see header schema)
 * @param {object} [snap.run] `{seed, seedSource, simulatedPs, steps, dtPs, startedAt}`
 * @param {object} [snap.trajectory] `{nFrames, spanPs, stridePs}`
 * @param {object[]} [snap.results] specs for `resultRow`
 * @param {object} [snap.runHealth] `{errors, ignored, lastError}` from src/errors.js
 * @param {object[]} [snap.coverageNotIncluded] families this record omits, with reasons
 * @param {string} [snap.generatedAt] ISO timestamp (defaults to now; tests inject)
 * @returns {object} the record (version === RESULTS_RECORD_VERSION)
 */
export function buildResultsRecord(snap = {}) {
  const s = snap && typeof snap === "object" ? snap : {};
  const input = s.input && typeof s.input === "object" ? s.input : {};
  const sel = s.selection && typeof s.selection === "object" ? s.selection : {};
  const ph = s.physics && typeof s.physics === "object" ? s.physics : {};
  const run = s.run && typeof s.run === "object" ? s.run : {};
  const traj = s.trajectory && typeof s.trajectory === "object" ? s.trajectory : {};
  const rh = s.runHealth && typeof s.runHealth === "object" ? s.runHealth : {};
  const tier = ph.tierFlags && typeof ph.tierFlags === "object" ? ph.tierFlags : {};

  const rows = [];
  for (const spec of Array.isArray(s.results) ? s.results.slice(0, RESULTS_MAX_VALUES) : []) {
    const row = resultRow(spec);
    if (row) rows.push(row);
  }

  const errors = Math.max(0, Math.round(finite(rh.errors, 0)));
  const ignored = Math.max(0, Math.round(finite(rh.ignored, 0)));
  const seed = finite(run.seed, null);

  return {
    version: RESULTS_RECORD_VERSION,
    format: "results-record-v1",
    app: "simulation_coarse",
    generatedAt: typeof s.generatedAt === "string" ? s.generatedAt : new Date().toISOString(),
    provenance: {
      code: { version: VERSION, buildDate: BUILD_DATE },
      run: {
        seed,
        // null seed is not a missing value: it is the honest statement that the
        // run used the unseeded Math.random path (src/integrator.js getSeed()
        // returns null there), which is NOT reproducible.
        seedSource: seed === null
          ? "unseeded-Math.random (src/integrator.js getSeed() === null) — this run is NOT reproducible from the record"
          : `SeededRNG (mulberry32, src/seeded-rng.js) seed ${seed}`,
        startedAt: typeof run.startedAt === "string" ? run.startedAt : null,
        simulatedPs: finite(run.simulatedPs, 0),
        steps: Math.max(0, Math.round(finite(run.steps, 0))),
        dtPs: finite(run.dtPs, null),
      },
      physics: {
        modelMode: RESULTS_MODEL_MODES.includes(ph.modelMode) ? ph.modelMode : "cg",
        physicsLevel: RESULTS_PHYSICS_LEVELS.includes(ph.physicsLevel) ? ph.physicsLevel : "L0",
        tierFlags: { charges: tier.charges === true, hbMode: String(tier.hbMode ?? "off"), weak: String(tier.weak ?? "off"), bindLog: tier.bindLog === true },
        temperatureK: finite(ph.temperatureK, null),
        frictionPerPs: finite(ph.frictionPerPs, null),
        beadMassDa: finite(ph.beadMassDa, null),
        dtPs: finite(ph.dtPs ?? run.dtPs, null),
        cutoffAng: finite(ph.cutoffAng, null),
        springGamma: finite(ph.springGamma, null),
        solventModel: typeof ph.solventModel === "string" ? ph.solventModel : null,
        saltM: finite(ph.saltM, null),
        epsIn: finite(ph.epsIn, null),
        epsOut: finite(ph.epsOut, null),
        sasaGamma: finite(ph.sasaGamma, null),
        electrostatics: "screened GB/Coulomb, switching cutoff 6.5→8.5 Å, NO PME (ROADMAP.md §1)",
        respa: { on: ph.respaOn === true, innerFs: finite(ph.respaInnerFs, null), outerFs: finite(ph.respaOuterFs, null) },
        binding: { on: ph.bindingOn !== false, holoSprings: ph.holoSprings !== false },
        system: {
          nParticles: Math.max(0, Math.round(finite(ph.nParticles, 0))),
          nProtein: Math.max(0, Math.round(finite(ph.nProtein, 0))),
          nLigandAtoms: Math.max(0, Math.round(finite(ph.nLigandAtoms, 0))),
          nLigandMolecules: Math.max(0, Math.round(finite(ph.nLigandMolecules, 0))),
        },
      },
      notAPhysicsParameter: {
        motionGain: finite(ph.motionGain, null),
        why: "the #motionGain slider is a Canvas2D display amplification (src/controllers/param-binding.js viewer.setMotionGain) — it scales pixels, not forces",
      },
    },
    input: {
      pdbId: typeof input.pdbId === "string" ? input.pdbId : "",
      fileName: typeof input.fileName === "string" ? input.fileName : "",
      origin: typeof input.origin === "string" ? input.origin : "unknown",
      contentHash: typeof input.pdbText === "string" && input.pdbText.length
        ? sha256Text(input.pdbText)
        : null,
      hashAlgorithm: "sha256",
      hashOf: "the raw PDB text the run parsed, UTF-8 bytes — the data/manifest.json convention, so the two are comparable with ===",
      selection: {
        chains: typeof sel.chains === "string" ? sel.chains : "",
        resFrom: finite(sel.resFrom, null),
        resTo: finite(sel.resTo, null),
        includeLigand: sel.includeLigand !== false,
      },
    },
    trajectory: {
      nFrames: Math.max(0, Math.round(finite(traj.nFrames, 0))),
      spanPs: finite(traj.spanPs, 0),
      stridePs: finite(traj.stridePs, null),
      framesOmitted: "all — this is a summary; the trajectory is a separate export (docs/EXPORT.md)",
    },
    results: rows,
    coverage: {
      included: rows.map((r) => r.id),
      notIncluded: (Array.isArray(s.coverageNotIncluded) ? s.coverageNotIncluded : []).map((c) => ({
        family: String(c?.family ?? "unknown"),
        why: String(c?.why ?? ""),
        howToGet: String(c?.howToGet ?? ""),
      })),
    },
    runHealth: {
      errors,
      documentedNoOps: ignored,
      // `clean` is derived, never supplied: a caller cannot mark a run with
      // recorded errors as clean by passing clean: true.
      clean: errors === 0,
      lastError: rh.lastError && typeof rh.lastError === "object"
        ? { context: String(rh.lastError.context ?? "unknown"), message: String(rh.lastError.message ?? ""), count: Math.max(1, Math.round(finite(rh.lastError.count, 1))) }
        : null,
      recordedBy: "src/errors.js recordError()/ignore() counters (see the top status bar and #topErrors)",
    },
    scope: scopeBlock(),
  };
}

/**
 * Serialize a record to pretty JSON (size-guarded, newline-terminated).
 * @param {object} rec
 * @returns {string}
 * @throws {Error} mapped SYSTEM_TOO_LARGE past RESULTS_MAX_BYTES
 */
export function serializeResultsRecord(rec) {
  const text = JSON.stringify(rec, null, 2) + "\n";
  if (text.length > RESULTS_MAX_BYTES) {
    throw createInputError(
      "SYSTEM_TOO_LARGE",
      `results record ${text.length} B over the ${RESULTS_MAX_BYTES} B limit (results must stay a summary — never embed frames)`
    );
  }
  return text;
}

/**
 * Validate a parsed record (pure, never throws).
 * @param {*} obj
 * @returns {{ok:boolean, error:Error|null}}
 */
export function validateResultsRecord(obj) {
  try {
    if (!obj || typeof obj !== "object" || Array.isArray(obj)) {
      return { ok: false, error: createInputError("GENERIC", "results record is not a JSON object") };
    }
    if (obj.version !== RESULTS_RECORD_VERSION) {
      return { ok: false, error: createInputError("GENERIC", `unsupported results-record version ${String(obj.version)} (want ${RESULTS_RECORD_VERSION})`) };
    }
    if (!Array.isArray(obj.results)) {
      return { ok: false, error: createInputError("GENERIC", "results record has no results array") };
    }
    return { ok: true, error: null };
  } catch (_) {
    return { ok: false, error: createInputError("GENERIC", null) };
  }
}

/**
 * Parse a results-record file back (pure, never throws uncaught).
 * @param {string} text
 * @returns {{ok:boolean, data:object|null, error:Error|null}}
 */
export function parseResultsRecord(text) {
  try {
    if (typeof text !== "string" || !text.trim()) {
      return { ok: false, data: null, error: createInputError("GENERIC", "empty results record") };
    }
    if (text.length > RESULTS_FILE_MAX_BYTES) {
      return { ok: false, data: null, error: createInputError("SYSTEM_TOO_LARGE", `results record file ${text.length} B over the ${RESULTS_FILE_MAX_BYTES} B load limit`) };
    }
    let obj;
    try { obj = JSON.parse(text); } catch (e) {
      return { ok: false, data: null, error: createInputError("GENERIC", e instanceof Error ? e.message : String(e)) };
    }
    const v = validateResultsRecord(obj);
    if (!v.ok) return { ok: false, data: null, error: v.error };
    return { ok: true, data: obj, error: null };
  } catch (_) {
    return { ok: false, data: null, error: createInputError("GENERIC", null) };
  }
}

/**
 * Characters that force a CSV cell to be quoted. A named constant rather than
 * an inline literal because tests/test_cache_contract.js's tokenizer tracks only
 * the LAST emitted CHARACTER (it calls `emit(c)` per character, so its
 * `REGEX_PRECEDER_WORDS` allowance never fires), which means an inline
 * `return /[",\n]/` reads as a `/` followed by an unterminated string to that
 * scanner. After `=` the scanner does recognise a regex, so this parses — and a
 * hoisted constant is the better shape for a table of CSV metacharacters anyway.
 */
const CSV_NEEDS_QUOTE = /[",\n]/;

/**
 * The tabular half of the record as CSV — the part a downstream script joins on.
 * `uncertainty` is the EMPTY STRING when there is no error bar (never `0`), and
 * `uncertainty_reason` carries the sentence explaining why.
 * @param {object} rec
 * @returns {string} CSV text
 */
export function resultsCsv(rec) {
  const rows = (rec && Array.isArray(rec.results) ? rec.results : []).map((r) => r || {});
  const q = (v) => {
    const s = String(v == null ? "" : v);
    return CSV_NEEDS_QUOTE.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  const out = [
    ["id", "label", "value", "unit", "uncertainty", "uncertainty_kind", "n", "uncertainty_reason", "meaning", "not_a"].join(","),
  ];
  for (const r of rows) {
    out.push([
      r.id, r.label, r.value, r.unit,
      r.uncertainty ? r.uncertainty.value : "",
      r.uncertainty ? r.uncertainty.kind : "",
      Number.isFinite(Number(r.n)) ? r.n : "",
      r.uncertaintyReason == null ? "" : r.uncertaintyReason,
      r.meaning, r.notA,
    ].map(q).join(","));
  }
  return out.join("\n") + "\n";
}

/** Split one CSV line, honouring the doubled-quote escape that resultsCsv writes. */
function splitCsvLine(line) {
  const cells = [];
  let cur = "", inQ = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (inQ) {
      if (ch === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else inQ = false; }
      else cur += ch;
    } else if (ch === '"') inQ = true;
    else if (ch === ",") { cells.push(cur); cur = ""; }
    else cur += ch;
  }
  cells.push(cur);
  return cells;
}

/**
 * Parse the CSV half back into rows (pure, never throws uncaught).
 *
 * `uncertainty` comes back as `null` for an empty cell — which is the point:
 * "no error bar" survives the round trip as "no error bar", and a consumer that
 * did `Number(cells[uncertainty]) || 0` gets 0 only by its own choice.
 * @param {string} text
 * @returns {{ok:boolean, rows:object[], error:Error|null}}
 */
export function parseResultsCsv(text) {
  try {
    if (typeof text !== "string" || !text.trim()) {
      return { ok: false, rows: [], error: createInputError("GENERIC", "empty results CSV") };
    }
    const lines = text.trim().split(/\r?\n/);
    const head = splitCsvLine(lines[0]);
    if (head[0] !== "id" || head[2] !== "value" || head[4] !== "uncertainty") {
      return { ok: false, rows: [], error: createInputError("GENERIC", `bad results CSV header ${JSON.stringify(lines[0].slice(0, 60))}`) };
    }
    const rows = [];
    for (let k = 1; k < lines.length; k++) {
      if (!lines[k].trim()) continue;
      const c = splitCsvLine(lines[k]);
      const value = Number(c[2]);
      if (!c[0] || !Number.isFinite(value)) {
        return { ok: false, rows: [], error: createInputError("GENERIC", `bad results CSV row ${k}: ${JSON.stringify(lines[k].slice(0, 60))}`) };
      }
      const u = c[4] === "" ? null : Number(c[4]);
      rows.push({
        id: c[0],
        value,
        unit: c[3],
        uncertainty: Number.isFinite(u) ? { value: u, kind: c[5] } : null,
        n: c[6] === "" ? null : Number(c[6]),
        uncertaintyReason: c[7] === "" ? null : c[7],
        meaning: c[8],
        notA: c[9],
      });
    }
    return { ok: true, rows, error: null };
  } catch (_) {
    return { ok: false, rows: [], error: createInputError("GENERIC", null) };
  }
}

/**
 * Should this record be allowed out of the browser?
 *
 * Fails CLOSED for the three things that would make the file a lie, and returns
 * WARNINGS for the things that make it merely limited:
 *
 *   BLOCKS — any recorded error (src/errors.js), because a run with recorded
 *     errors must not be exported as if it were clean. The counter is cumulative
 *     from page load (src/errors.js is not per-run), so the remedy says reload
 *     rather than pretending a fresh run would clear it.
 *     Any result with no `VALUE_TRUST_BOUNDARY` entry, because a value with no
 *     stated meaning must not ship as a clean-looking JSON.
 *     No content hash, because an un-identified input cannot be reproduced.
 *
 *   WARNS — an unseeded run (the integrator's default Math.random path, which is
 *     what most browser runs actually are: `getSeed()` returns null there). This
 *     is deliberately NOT a block: blocking it would make the export unusable for
 *     every default session while adding no honesty, because the record already
 *     says in `provenance.run.seedSource` that the run is not reproducible.
 *
 * @param {object} rec
 * @returns {{ok:boolean, reason:string, remedy:string, warnings:string[]}}
 */
export function resultsExportGate(rec) {
  const deny = (reason, remedy) => ({ ok: false, reason, remedy, warnings: [] });
  if (!rec || typeof rec !== "object") return deny("no record built", "build the system and run a trajectory first");
  if (rec.version !== RESULTS_RECORD_VERSION) return deny(`record version ${String(rec.version)} is not ${RESULTS_RECORD_VERSION}`, "re-export with the current build");
  const rh = rec.runHealth || {};
  if (Number(rh.errors) > 0) {
    const last = rh.lastError ? ` (last: ${rh.lastError.context}: ${rh.lastError.message})` : "";
    return deny(
      `${rh.errors} error(s) recorded since page load${last}`,
      "a run with recorded errors is not exportable as if it were clean — reload the page to clear the counter, re-run, then export"
    );
  }
  const missing = (Array.isArray(rec.results) ? rec.results : []).filter((r) => r && r.honestyMissing).map((r) => r.id);
  if (missing.length) {
    return deny(
      `no trust-boundary statement for ${missing.join(", ")} (src/scope.js VALUE_TRUST_BOUNDARY)`,
      "add the entry, then re-export — a value with no stated meaning must not ship as a clean JSON"
    );
  }
  if (!rec.input || !rec.input.contentHash) {
    return deny("no content hash for the input structure", "load a structure (PDB id, preset or file) so its bytes can be hashed");
  }
  const warnings = [];
  if (rec.provenance?.run?.seed == null) {
    warnings.push("run is UNSEEDED (integrator on the Math.random path) — the record names this, but the run cannot be reproduced from it");
  }
  if (Number(rh.documentedNoOps) > 0) {
    warnings.push(`${rh.documentedNoOps} documented no-op catch(es) fired (recorded, not errors) — counted in runHealth.documentedNoOps`);
  }
  const nF = Number(rec.trajectory?.nFrames ?? 0);
  if (nF < 2) warnings.push(`only ${nF} recorded frame(s) — the analysis values below are from a single frame`);
  const noBar = (Array.isArray(rec.results) ? rec.results : []).filter((r) => r && r.uncertaintyAvailable === false).map((r) => r.id);
  if (noBar.length) warnings.push(`${noBar.length} of ${rec.results.length} value(s) have NO error bar: ${noBar.join(", ")}`);
  return { ok: true, reason: "", remedy: "", warnings };
}

/**
 * Conventional download filename: `results_<pdbId|input>_v<VERSION>.json`.
 * @param {object} rec
 * @returns {string}
 */
export function resultsFilename(rec) {
  const id = String(rec?.input?.pdbId || rec?.input?.fileName || "run").replace(/\W+/g, "_") || "run";
  return `results_${id}_v${rec?.provenance?.code?.version ?? VERSION}.json`;
}
