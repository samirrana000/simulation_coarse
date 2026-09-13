/**
 * test_session_roundtrip.js — FP4 session save/load + export parse-back.
 *
 * Fast (<10s), deterministic (no Math.random, no wall-clock asserts).
 * Run: node tests/test_session_roundtrip.js
 */

import {
  SESSION_VERSION, SESSION_MAX_BYTES, SESSION_FILE_MAX_BYTES,
  buildSession, serializeSession, validateSession, parseSession,
  dccmCsv, parseDccmCsv, trajectoryJson, settingsJson,
  estimateFramesBytes,
} from "../src/session.js";
import { BindLog } from "../src/capture/bindlog.js";
import { Funnel } from "../src/funnel.js";
import { pmfCsv } from "../src/analysis.js";

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

// ---------------------------------------------------------------- session --
console.log("=== session save→load roundtrip ===");
const snap = {
  pdbId: "4W52", modelMode: "heavy", chains: "A", resFrom: 1, resTo: 100,
  includeLig: true, physicsLevel: "L2",
  ligand: { selected: "BNZ" }, thermoLig: "all",
  settings: {
    backend: "cpu", numThreads: 4, solventModel: "gb", saltM: 0.15,
    epsIn: 4.0, epsOut: 78.5, sasaGamma: 0.0072,
    respaOn: true, respaOuterFs: 4, chemicalNetworkOn: true,
  },
  dynamics: {
    rc: 10, gamma: 2.0, temp: 300, fric: 8.0, mass: 110,
    motionGain: 1.3, bindPot: true, holoSprings: true,
  },
  recording: { stridePs: 2.0, maxFrames: 500, exportFmt: "json" },
  recorderMeta: { count: 42, spanPs: 84.0 },
};
const sess = buildSession(snap);
assert(sess.version === SESSION_VERSION && SESSION_VERSION === 1, `session version === 1 (got ${sess.version})`);
assert(sess.physicsLevel === "L2", `physicsLevel round-trips L2 (got ${sess.physicsLevel})`);
assert(sess.thermoLig === "all", `picker thermoLig round-trips "all" (got ${sess.thermoLig})`);
assert(sess.ligand.selected === "BNZ", `ligand selection round-trips BNZ (got ${sess.ligand.selected})`);
assert(sess.settings.backend === "cpu" && sess.settings.solventModel === "gb", "settings backend/solvent round-trip");
assert(sess.dynamics.temp === 300 && sess.dynamics.gamma === 2.0, "dynamics temp/gamma round-trip");
assert(sess.recording.exportFmt === "json", `recording exportFmt round-trips json (got ${sess.recording.exportFmt})`);
assert(sess.recorderMeta.count === 42, `recorderMeta count round-trips 42 (got ${sess.recorderMeta.count})`);
assert(sess.recorderMeta.includeFrames === false, "recorderMeta includeFrames is false (counts only)");
assert(!("frames" in sess), "no full-frame persist (no `frames` key on the session)");

const text = serializeSession(sess);
assert(text.length < SESSION_MAX_BYTES, `session ${text.length} B < ${SESSION_MAX_BYTES} B persist limit`);
assert(text.includes('"version": 1'), "serialized session carries version 1");

const parsed = parseSession(text);
assert(parsed.ok === true && parsed.data !== null, "parseSession(save) ok");
assert(parsed.data.physicsLevel === "L2", "load restores physicsLevel L2");
assert(parsed.data.thermoLig === "all", "load restores picker thermoLig all");
assert(parsed.data.settings.backend === "cpu", "load restores settings.backend cpu");
assert(parsed.data.dynamics.temp === 300, "load restores dynamics.temp 300");
assert(parsed.data.recorderMeta.count === 42, "load restores recorderMeta.count 42");

// ---------------------------------------------------------------- guards ---
console.log("=== session guards (input_errors codes) ===");
const empty = parseSession("");
assert(empty.ok === false && empty.error && empty.error.code === "GENERIC", `empty session → GENERIC (got ${empty.error && empty.error.code})`);
const badJson = parseSession("{not json");
assert(badJson.ok === false && badJson.error.code === "GENERIC", "malformed JSON → GENERIC");
const badVer = validateSession({ version: 99 });
assert(badVer.ok === false && badVer.error.code === "GENERIC", "bad version → GENERIC");
const badLvl = validateSession({ version: 1, physicsLevel: "L9" });
assert(badLvl.ok === false, "bad physicsLevel fails validation");
const fallback = buildSession({ physicsLevel: "L9" });
assert(fallback.physicsLevel === "L0", `bad physicsLevel falls back to L0 default (got ${fallback.physicsLevel})`);
const huge = parseSession("x".repeat(SESSION_FILE_MAX_BYTES + 1));
assert(huge.ok === false && huge.error.code === "SYSTEM_TOO_LARGE", "oversized session file → SYSTEM_TOO_LARGE");
assert(estimateFramesBytes(500, 164) === 500 * 164 * 12, `estimateFramesBytes(500,164) = ${estimateFramesBytes(500, 164)}`);

// ------------------------------------------------------------------ BLG1 ---
console.log("=== BLG1 blob parse-back ===");
const bl = new BindLog({ quantStep: 0.01 });
const nAtoms = 4;
for (let f = 0; f < 3; f++) {
  const pos = new Float32Array(nAtoms * 3);
  for (let i = 0; i < pos.length; i++) pos[i] = i * 0.5 + f * 0.05;
  bl.captureFrame(pos, f * 2.0);
}
bl.pushEnergyComponents(0, [-1.9, 0.5, -0.1, -5.7, 0, 0, 0]);
bl.pushHill(2.0, 4.2, 0.05);
bl.pushContact(2.0, true, 1, 2, 3.5);
bl.pushState(4.0, 3);
bl.pushPocketVolume(4.0, 620);
const buf = bl.toBinaryBlob();
const dv = new DataView(buf);
assert(dv.getUint32(0, true) === 0x424c4731, "BLG1 magic ok (\"BLG1\")");
assert(dv.getUint32(4, true) === 1, "BLG1 version === 1");
const hLen = dv.getUint32(8, true);
assert(hLen > 0, `BLG1 header length ${hLen} > 0`);
let headerStr = "";
const u8 = new Uint8Array(buf);
for (let i = 0; i < hLen; i++) headerStr += String.fromCharCode(u8[12 + i]);
const header = JSON.parse(headerStr);
assert(header.nF === 3 && header.n === nAtoms, `BLG1 header nF=3 n=4 (got nF=${header.nF} n=${header.n})`);
assert(header.nE === bl.nEvents, `BLG1 header nE matches (${header.nE} events)`);
const back = BindLog.fromBinaryBlob(buf);
assert(back.nFrames === 3, `BLG1 parse-back frames ${back.nFrames}`);
assert(back.nEvents === bl.nEvents, `BLG1 parse-back events ${back.nEvents}`);
assert(back.frameTimes.length === 3, "BLG1 parse-back frame times ×3");

// --------------------------------------------------------------- PMF CSV ---
console.log("=== PMF CSV parse-back ===");
const nProt = 6, nLig = 2;
const ref = new Float64Array(3 * (nProt + nLig));
for (let i = 0; i < nProt; i++) { ref[3 * i] = i * 3.0; ref[3 * i + 1] = 0; ref[3 * i + 2] = 0; }
for (let a = 0; a < nLig; a++) { const c = 3 * (nProt + a); ref[c] = 2.0; ref[c + 1] = 0.2 * a; ref[c + 2] = 0; }
const funnel = new Funnel({ nProt, n: nProt + nLig, ref, biasFactor: 6, rMax: 24, bins: 96 });
funnel.T = 300;
for (let i = 0; i < 5; i++) funnel.deposit(3 + i * 0.4);
const csv = pmfCsv(funnel);
const lines = csv.trim().split("\n");
assert(/# T=300, gamma=6, hills=5, V0=1660\.54/.test(lines[0]), `PMF header provenance (got ${lines[0].slice(0, 44)}…)`);
assert(lines[2] === "r_Ang,pMF_kcal_per_mol", "PMF column header row ok");
const { r } = funnel.getPMF();
assert(lines.slice(3, 3 + r.length).length === r.length, `PMF CSV parse-back row count ${r.length}`);
assert(csv.includes("# hills,5"), "PMF footer hills=5");

// ---------------------------------------------------------------- DCCM ----
console.log("=== DCCM CSV parse-back ===");
const dn = 4;
const dm = new Float64Array(dn * dn);
for (let i = 0; i < dn; i++) for (let j = 0; j < dn; j++) dm[i * dn + j] = i === j ? 1 : (i - j) * 0.1;
const dcsv = dccmCsv(dm, dn);
const dback = parseDccmCsv(dcsv);
assert(dback.ok === true && dback.n === dn, `DCCM CSV parse-back n=${dback.n}`);
assert(dback.matrix.length === dn * dn, `DCCM CSV parse-back rows ${dback.matrix.length / dn}×${dn}`);
assert(dback.matrix[0] === 1 && dback.matrix[5] === 1, "DCCM diagonal parses back to 1");
assert(Math.abs(dback.matrix[1] - dm[1]) < 1e-6, "DCCM off-diagonal parses back");

// ------------------------------------------------------- trajectory JSON ---
console.log("=== trajectory JSON + settings JSON ===");
const fr = [new Float32Array([1, 2, 3, 4, 5, 6]), new Float32Array([1.1, 2.1, 3.1, 4.1, 5.1, 6.1])];
const tj = trajectoryJson(fr, [0, 2.0], { T: 300, gamma: 2, seed: 7, date: "2026-09-13" });
const tjBack = JSON.parse(tj);
assert(tjBack.count === 2 && tjBack.nAtoms === 2, `trajectory JSON count=2 nAtoms=2 (got ${tjBack.count}/${tjBack.nAtoms})`);
assert(tjBack.frames.length === 2 && tjBack.times.length === 2, "trajectory JSON frames+times parse back");
let threw = false;
try { trajectoryJson([], []); } catch (_) { threw = true; }
assert(threw === true, "trajectoryJson([]) throws (no silent empty export)");
const sjBack = JSON.parse(settingsJson({ physicsLevel: "L2", backend: "workers", numThreads: 8 }));
assert(sjBack["sim.physicsLevel"] === "L2", "settings JSON carries sim.physicsLevel L2");
assert(sjBack.backend === "workers" && sjBack.numThreads === 8, "settings JSON backend/threads parse back");

console.log(`\n${passed} PASSED, ${failed} FAILED`);
process.exit(failed ? 1 : 0);
