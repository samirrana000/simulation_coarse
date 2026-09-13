/**
 * session.js — FP4 export/persistence: session save/load + export-matrix builders.
 *
 * Gap closed: BindLog BLG1 blob existed (`src/capture/bindlog.js:166`
 * `toBinaryBlob`), PMF CSV + thermo table existed ad hoc
 * (`src/analysis.js:619` `pmfCsv`, `src/analysis/thermodynamics.js:548`
 * `formatThermoTable`), but no session save/load roundtrip and no documented
 * export matrix. This module is the headless-testable core (zero deps beyond
 * the existing `src/input_errors.js` validators; no DOM access so plain Node
 * can import it):
 *
 *   Session file (JSON, `{version: 1, ...}`):
 *     {version, app, savedAt, pdbId, modelMode, chains, resFrom, resTo,
 *      includeLig, physicsLevel, ligand:{selected}, thermoLig,
 *      settings:{...sim.* keys}, dynamics:{...}, recording:{...},
 *      recorderMeta:{count, spanPs, includeFrames:false, framesOmitted, note}}
 *   Limits (documented, enforced): full PDB text and full trajectory frames
 *   are NEVER persisted — `recorderMeta` carries counts only
 *   (`includeFrames` is always false on save; the `frames` key is never
 *   written). Serialized sessions over `SESSION_MAX_BYTES` (256 KiB) refuse
 *   to serialize; session files over `SESSION_FILE_MAX_BYTES` (1 MiB) refuse
 *   to load. No multi-MB localStorage: callers must not stash sessions in
 *   localStorage (download/upload files only).
 *
 *   Export builders (pure, all downloads go through the existing
 *   `downloadText`/blob pattern in `src/recorder.js:168`):
 *     dccmCsv(matrix, n) → "i,j,C" CSV for the DCCM subpanel
 *     parseDccmCsv(text) → {n, matrix} (roundtrip parse-back)
 *     trajectoryJson(frames, times, meta) → JSON frames export
 *     settingsJson(st) → sim.* settings JSON export
 *     downloadBlob(buffer, filename, mime) → binary download (BLG1)
 *
 * Load validation reuses `src/input_errors.js` failure classes: malformed
 * JSON / bad version / bad enum → GENERIC; oversized file → SYSTEM_TOO_LARGE;
 * empty input → EMPTY_PDB-family via `validatePdbText`-style check mapped to
 * GENERIC with session technical detail (sessions are not PDB text, so the
 * code is GENERIC, never a PDB class).
 */

/* eslint-disable no-unused-vars */
import { createInputError } from "./input_errors.js?v=10";

/** Session schema version (frozen contract — tests assert === 1). */
export const SESSION_VERSION = 1;

/** Refuse to serialize sessions larger than this (bytes, UTF-8 approx). */
export const SESSION_MAX_BYTES = 256 * 1024; // 256 KiB — sessions are ~1 KiB (counts only, never frames)

/** Refuse to load session files larger than this (bytes). */
export const SESSION_FILE_MAX_BYTES = 1024 * 1024; // 1 MiB

/** Physics levels accepted by the session schema (mirrors settings-panel). */
export const SESSION_PHYSICS_LEVELS = ["L0", "L1", "L2"];

/** Model modes accepted by the session schema. */
export const SESSION_MODEL_MODES = ["cg", "heavy"];

/**
 * Estimate raw trajectory bytes for a frame count (Float32 acquisition).
 * @param {number} nFrames frame count
 * @param {number} nAtoms atom/bead count
 * @returns {number} estimated bytes (nFrames × nAtoms × 3 × 4)
 */
export function estimateFramesBytes(nFrames, nAtoms) {
  const f = Number(nFrames) || 0, a = Number(nAtoms) || 0;
  return Math.max(0, Math.round(f)) * Math.max(0, Math.round(a)) * 12;
}

/**
 * Build a session object from a UI snapshot (pure, never throws on missing
 * fields — absent keys fall back to app defaults, which are unchanged).
 * NEVER stores full frames or PDB text (counts/meta only).
 * @param {object} [snap] flat snapshot (see header schema)
 * @returns {object} session object (version === SESSION_VERSION)
 */
export function buildSession(snap = {}) {
  const s = snap && typeof snap === "object" ? snap : {};
  const settings = s.settings && typeof s.settings === "object" ? s.settings : {};
  const dynamics = s.dynamics && typeof s.dynamics === "object" ? s.dynamics : {};
  const recording = s.recording && typeof s.recording === "object" ? s.recording : {};
  const recMeta = s.recorderMeta && typeof s.recorderMeta === "object" ? s.recorderMeta : {};
  const ligand = s.ligand && typeof s.ligand === "object" ? s.ligand : {};
  const count = Number(recMeta.count) || 0;
  return {
    version: SESSION_VERSION,
    app: "simulation_coarse",
    savedAt: new Date().toISOString(),
    pdbId: typeof s.pdbId === "string" ? s.pdbId : "",
    modelMode: SESSION_MODEL_MODES.includes(s.modelMode) ? s.modelMode : "cg",
    chains: typeof s.chains === "string" ? s.chains : "",
    resFrom: s.resFrom ?? null,
    resTo: s.resTo ?? null,
    includeLig: s.includeLig !== false,
    physicsLevel: SESSION_PHYSICS_LEVELS.includes(s.physicsLevel) ? s.physicsLevel : "L0",
    ligand: { selected: typeof ligand.selected === "string" ? ligand.selected : "" },
    thermoLig: typeof s.thermoLig === "string" ? s.thermoLig : "auto",
    settings: {
      backend: typeof settings.backend === "string" ? settings.backend : "auto",
      numThreads: Number.isFinite(Number(settings.numThreads)) ? Number(settings.numThreads) : 2,
      solventModel: typeof settings.solventModel === "string" ? settings.solventModel : "gb",
      saltM: Number.isFinite(Number(settings.saltM)) ? Number(settings.saltM) : 0.15,
      epsIn: Number.isFinite(Number(settings.epsIn)) ? Number(settings.epsIn) : 4.0,
      epsOut: Number.isFinite(Number(settings.epsOut)) ? Number(settings.epsOut) : 78.5,
      sasaGamma: Number.isFinite(Number(settings.sasaGamma)) ? Number(settings.sasaGamma) : 0.0072,
      respaOn: settings.respaOn === true,
      respaOuterFs: Number.isFinite(Number(settings.respaOuterFs)) ? Number(settings.respaOuterFs) : 4,
      chemicalNetworkOn: settings.chemicalNetworkOn !== false,
    },
    dynamics: {
      rc: Number.isFinite(Number(dynamics.rc)) ? Number(dynamics.rc) : 10,
      gamma: Number.isFinite(Number(dynamics.gamma)) ? Number(dynamics.gamma) : 2.0,
      temp: Number.isFinite(Number(dynamics.temp)) ? Number(dynamics.temp) : 300,
      fric: Number.isFinite(Number(dynamics.fric)) ? Number(dynamics.fric) : 8.0,
      mass: Number.isFinite(Number(dynamics.mass)) ? Number(dynamics.mass) : 110,
      motionGain: Number.isFinite(Number(dynamics.motionGain)) ? Number(dynamics.motionGain) : 1.3,
      bindPot: dynamics.bindPot !== false,
      holoSprings: dynamics.holoSprings !== false,
    },
    recording: {
      stridePs: Number.isFinite(Number(recording.stridePs)) ? Number(recording.stridePs) : 2.0,
      maxFrames: Number.isFinite(Number(recording.maxFrames)) ? Number(recording.maxFrames) : 500,
      exportFmt: typeof recording.exportFmt === "string" ? recording.exportFmt : "xyz",
    },
    recorderMeta: {
      count,
      spanPs: Number.isFinite(Number(recMeta.spanPs)) ? Number(recMeta.spanPs) : 0,
      includeFrames: false, // cutoff: full frames are never persisted (see header)
      framesOmitted: count > 0,
      note: count > 0
        ? `counts only — ${count} frame(s) stay in memory; re-record after load`
        : "no frames recorded at save time",
    },
  };
}

/**
 * Serialize a session to pretty JSON (size-guarded).
 * @param {object} sess session object
 * @returns {string} JSON text
 * @throws {Error} mapped SYSTEM_TOO_LARGE when over SESSION_MAX_BYTES
 */
export function serializeSession(sess) {
  const text = JSON.stringify(sess, null, 2) + "\n";
  if (text.length > SESSION_MAX_BYTES) {
    throw createInputError(
      "SYSTEM_TOO_LARGE",
      `session ${text.length} B over the ${SESSION_MAX_BYTES} B persist limit (frames/PDB text must not be embedded)`
    );
  }
  return text;
}

/**
 * Validate a parsed session object (pure, never throws).
 * @param {*} obj parsed JSON value
 * @returns {{ok:boolean, error:Error|null}}
 */
export function validateSession(obj) {
  try {
    if (!obj || typeof obj !== "object" || Array.isArray(obj)) {
      return { ok: false, error: createInputError("GENERIC", "session is not a JSON object") };
    }
    if (obj.version !== SESSION_VERSION) {
      return {
        ok: false,
        error: createInputError("GENERIC", `unsupported session version ${String(obj.version)} (want ${SESSION_VERSION})`),
      };
    }
    if (obj.physicsLevel !== undefined && !SESSION_PHYSICS_LEVELS.includes(obj.physicsLevel)) {
      return { ok: false, error: createInputError("GENERIC", `bad physicsLevel ${String(obj.physicsLevel)}`) };
    }
    if (obj.modelMode !== undefined && !SESSION_MODEL_MODES.includes(obj.modelMode)) {
      return { ok: false, error: createInputError("GENERIC", `bad modelMode ${String(obj.modelMode)}`) };
    }
    return { ok: true, error: null };
  } catch (_) {
    return { ok: false, error: createInputError("GENERIC", null) };
  }
}

/**
 * Parse session file text (pure, never throws uncaught).
 * @param {*} text raw file text
 * @returns {{ok:boolean, data:object|null, error:Error|null}}
 */
export function parseSession(text) {
  try {
    if (typeof text !== "string" || !text.trim()) {
      return { ok: false, data: null, error: createInputError("GENERIC", "empty session text") };
    }
    if (text.length > SESSION_FILE_MAX_BYTES) {
      return {
        ok: false,
        data: null,
        error: createInputError(
          "SYSTEM_TOO_LARGE",
          `session file ${text.length} B over the ${SESSION_FILE_MAX_BYTES} B load limit`
        ),
      };
    }
    let obj;
    try {
      obj = JSON.parse(text);
    } catch (e) {
      return { ok: false, data: null, error: createInputError("GENERIC", e instanceof Error ? e.message : String(e)) };
    }
    const v = validateSession(obj);
    if (!v.ok) return { ok: false, data: null, error: v.error };
    return { ok: true, data: obj, error: null };
  } catch (_) {
    return { ok: false, data: null, error: createInputError("GENERIC", null) };
  }
}

/**
 * Serialize a DCCM matrix as CSV (pure).
 * @param {Float64Array|Array} matrix row-major n×n correlations
 * @param {number} n residue count
 * @returns {string} CSV text ("i,j,C" header + n² rows)
 */
export function dccmCsv(matrix, n) {
  if (!Number.isInteger(n) || n < 2) throw new Error("dccmCsv: n (≥ 2) is required.");
  if (!matrix || matrix.length < n * n) throw new Error(`dccmCsv: matrix too short (${matrix ? matrix.length : 0} < ${n * n}).`);
  const out = ["i,j,C"];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const c = Number(matrix[i * n + j]);
      if (!Number.isFinite(c)) throw new Error(`dccmCsv: non-finite C at (${i},${j}).`);
      out.push(`${i},${j},${c.toFixed(6)}`);
    }
  }
  return out.join("\n") + "\n";
}

/**
 * Parse a DCCM CSV back (pure, never throws uncaught).
 * @param {*} text CSV text
 * @returns {{ok:boolean, n:number, matrix:Float64Array|null, error:Error|null}}
 */
export function parseDccmCsv(text) {
  try {
    if (typeof text !== "string" || !text.trim()) {
      return { ok: false, n: 0, matrix: null, error: createInputError("GENERIC", "empty DCCM CSV") };
    }
    const lines = text.trim().split(/\r?\n/);
    if (lines[0].trim() !== "i,j,C") {
      return { ok: false, n: 0, matrix: null, error: createInputError("GENERIC", `bad DCCM header ${JSON.stringify(lines[0])}`) };
    }
    const rows = lines.slice(1).filter((l) => l.trim() !== "");
    const n = Math.round(Math.sqrt(rows.length));
    if (n < 2 || n * n !== rows.length) {
      return { ok: false, n: 0, matrix: null, error: createInputError("GENERIC", `DCCM row count ${rows.length} is not a square ≥ 2×2`) };
    }
    const m = new Float64Array(n * n);
    for (let k = 0; k < rows.length; k++) {
      const parts = rows[k].split(",");
      if (parts.length !== 3) {
        return { ok: false, n: 0, matrix: null, error: createInputError("GENERIC", `bad DCCM row ${k}`) };
      }
      const i = Number(parts[0]), j = Number(parts[1]), c = Number(parts[2]);
      if (!Number.isInteger(i) || !Number.isInteger(j) || !Number.isFinite(c)) {
        return { ok: false, n: 0, matrix: null, error: createInputError("GENERIC", `bad DCCM row ${k}`) };
      }
      m[i * n + j] = c;
    }
    return { ok: true, n, matrix: m, error: null };
  } catch (_) {
    return { ok: false, n: 0, matrix: null, error: createInputError("GENERIC", null) };
  }
}

/**
 * Serialize recorded frames as JSON (pure; export download only — never used
 * for session persist, so no size cutoff here beyond available memory).
 * @param {Array<Float32Array|Float64Array>} frames snapshots
 * @param {Array<number>} times ps per frame
 * @param {object} [meta] provenance ({T, gamma, seed, date, nProt})
 * @returns {string} JSON text
 */
export function trajectoryJson(frames, times, meta = {}) {
  if (!frames || frames.length === 0) throw new Error("No frames recorded yet.");
  const fr = frames.map((f) => Array.from(f));
  return JSON.stringify({
    version: 1,
    format: "trajectory-json-v1",
    count: frames.length,
    nAtoms: frames[0].length / 3,
    times: Array.from(times ?? []),
    meta: { T: meta.T ?? 300, gamma: meta.gamma ?? 2, seed: meta.seed ?? 0, date: meta.date ?? "" },
    frames: fr,
  }) + "\n";
}

/**
 * Serialize the sim.* settings subset as JSON (pure).
 * @param {object} st settingsState-like object
 * @returns {string} JSON text
 */
export function settingsJson(st = {}) {
  const s = st && typeof st === "object" ? st : {};
  return JSON.stringify({
    version: 1,
    format: "settings-json-v1",
    "sim.physicsLevel": SESSION_PHYSICS_LEVELS.includes(s.physicsLevel) ? s.physicsLevel : "L0",
    backend: s.backend ?? "auto",
    numThreads: s.numThreads ?? 2,
    solventModel: s.solventModel ?? "gb",
    saltM: s.saltM ?? 0.15,
    epsIn: s.epsIn ?? 4.0,
    epsOut: s.epsOut ?? 78.5,
    sasaGamma: s.sasaGamma ?? 0.0072,
    respaOn: s.respaOn === true,
    respaOuterFs: s.respaOuterFs ?? 4,
    chemicalNetworkOn: s.chemicalNetworkOn !== false,
  }, null, 2) + "\n";
}

/**
 * Trigger a browser download of binary `buffer` (BLG1 pattern: mirrors
 * `downloadText` in `src/recorder.js:168`). Headless-safe (returns false).
 * @param {ArrayBuffer|Uint8Array} buffer binary payload
 * @param {string} filename download name
 * @param {string} [mime] MIME type
 * @returns {boolean} true when a download was triggered
 */
export function downloadBlob(buffer, filename, mime = "application/octet-stream") {
  try {
    if (typeof document === "undefined" || typeof URL === "undefined") return false;
    const u8 = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
    const blob = new Blob([u8.slice().buffer], { type: mime });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 5000);
    return true;
  } catch (_) {
    return false;
  }
}
