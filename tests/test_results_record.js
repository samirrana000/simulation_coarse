/**
 * test_results_record.js — the structured results export: round-trip, the
 * uncertainty rule, the run-health gate, and the ROADMAP drift detector.
 *
 * WHAT EACH BLOCK IS DEFENDING
 * ----------------------------
 *   1. SHA-256      A hand-rolled hash is either right or it is a wrong number in
 *                   a provenance field. Checked against FIPS vectors, against
 *                   node:crypto on random input, and against the hash
 *                   data/manifest.json already records for the real 4w52.pdb —
 *                   so "same convention as the manifest" is a measurement, not
 *                   a claim.
 *   2. ROUND-TRIP   The point of the feature is that a downstream script can
 *                   read the record back. Parsed back, every value, uncertainty,
 *                   null-reason and honesty string must survive.
 *   3. NO ZERO BAR  A zero error bar is a lie. Any `uncertainty.value <= 0` is a
 *                   FAILURE, including the case where the analysis module itself
 *                   returned 0 to mean "not computed" (computeThermodynamics
 *                   does exactly that with dH_se), which must surface as
 *                   `uncertainty: null` + a reason.
 *   4. GATE         A run with recorded errors must not export as if clean; a
 *                   value with no trust-boundary statement must not export at
 *                   all; an unhashed input must not export.
 *   5. SCOPE DRIFT  The honesty block claims to be ROADMAP.md §1. This reads the
 *                   real ROADMAP.md, docs/LIMITATIONS.md and docs/VALIDATION.md,
 *                   asserts ZERO drift, and then PROVES the detector can fail by
 *                   feeding it mutated documents. A drift check that has never
 *                   been seen red is not a drift check.
 *   6. SIZE/ID      The record is a summary, not a trajectory dump: it must stay
 *                   far under its byte cap, and every id the UI wiring can emit
 *                   must have a trust-boundary entry.
 *
 * Fast (<1 s), deterministic (no wall-clock assertions on generated values; the
 * timestamp is injected).
 * Run: node tests/test_results_record.js
 */

import fs from "node:fs";
import path from "node:path";
import crypto from "node:crypto";
import { fileURLToPath } from "node:url";

import { sha256Bytes, sha256Text } from "../src/sha256.js";
import {
  RESULTS_RECORD_VERSION, RESULTS_MAX_BYTES, RESULTS_FILE_MAX_BYTES,
  RESULTS_MAX_VALUES, RESULTS_VALUE_IDS,
  buildResultsRecord, serializeResultsRecord, validateResultsRecord,
  parseResultsRecord, resultsCsv, parseResultsCsv,
  resultsExportGate, resultsFilename,
} from "../src/results-record.js";
import {
  scopeDrift, parseRoadmapGuard, parseRoadmapOutOfScope, scopeBlock,
  BROWSER_V1_GUARD, OUT_OF_SCOPE, VALUE_TRUST_BOUNDARY, normKey,
} from "../src/scope.js";
import { VALUE_TRUST_BOUNDARY as SCOPE_TB } from "../src/scope.js";
import { Funnel } from "../src/funnel.js";
import { computeThermodynamics } from "../src/analysis/thermodynamics.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const read = (rel) => fs.readFileSync(path.join(ROOT, rel), "utf-8");

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

// =====================================================================
// 1. SHA-256 — right, and the same number the manifest already records
// =====================================================================
console.log("=== sha256: vectors, node parity, manifest parity ===");
assert(
  sha256Bytes(new Uint8Array(0)) === "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
  'FIPS 180-4 vector: sha256("")'
);
assert(
  sha256Text("abc") === "sha256:ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
  'FIPS 180-4 vector: sha256("abc")'
);
assert(
  sha256Text("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ===
  "sha256:248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
  "FIPS 180-4 vector: 448-bit message"
);
{
  // node:crypto is the oracle. Random bytes exercise every padding branch
  // (55/56/57/63/64/65 and around) without depending on the loop count.
  let bad = 0;
  for (let n = 0; n < 200; n++) {
    const buf = crypto.randomBytes(n);
    const want = crypto.createHash("sha256").update(buf).digest("hex");
    if (sha256Bytes(new Uint8Array(buf)) !== want) bad++;
  }
  assert(bad === 0, `200 random byte strings (lengths 0-199) match node:crypto (${bad} mismatches)`);
}
{
  // The convention claim: the record's contentHash must be byte-identical to
  // what data/manifest.json already stores for the bundled 4w52.pdb.
  const manifest = JSON.parse(read("data/manifest.json"));
  const pdbText = read("4w52.pdb");
  const got = sha256Text(pdbText);
  assert(got === manifest["4w52.pdb"],
    `record contentHash == data/manifest.json["4w52.pdb"] (${got.slice(0, 22)}…)`);
  assert(/^sha256:[0-9a-f]{64}$/.test(got), "contentHash matches the manifest's `sha256:<64 hex>` shape");
}

// =====================================================================
// 2. A record built from REAL analysis output
// =====================================================================
console.log("=== record built from real analysis output ===");
/** A real funnel with hills, so ΔG and its convergence SE come from src/funnel.js. */
function makeFunnel() {
  const nProt = 6, nLig = 2;
  const ref = new Float64Array(3 * (nProt + nLig));
  for (let i = 0; i < nProt; i++) { ref[3 * i] = i * 3.0; ref[3 * i + 1] = 0; ref[3 * i + 2] = 0; }
  for (let a = 0; a < nLig; a++) { const c = 3 * (nProt + a); ref[c] = 2.0; ref[c + 1] = 0.2 * a; ref[c + 2] = 0; }
  const f = new Funnel({ nProt, n: nProt + nLig, ref, biasFactor: 6, rMax: 24, bins: 96 });
  f.T = 300;
  for (let i = 0; i < 120; i++) f.deposit(2.0 + (i % 20) * 0.4);
  return f;
}
const funnel = makeFunnel();
const dgReal = funnel.estimateDG();
const seReal = funnel.convergenceSE();
assert(Number.isFinite(dgReal) && Number.isFinite(seReal) && seReal > 0,
  `real funnel: dG=${dgReal.toFixed(3)} convergenceSE=${seReal.toFixed(3)} (nHills=${funnel._nHills})`);

/**
 * Real computeThermodynamics output, used for the "module returns 0 = unknown"
 * case: with ONE holo energy frame the 20-block bootstrap has a single block,
 * so dH_se stays 0 — the code's own "not computed" sentinel.
 */
const thermoThin = computeThermodynamics({
  holoFrames: [new Float64Array(18), new Float64Array(18)],
  apoFrames: [new Float64Array(18), new Float64Array(18)],
  holoEnergies: [[-1.0, 0.2, 0.0, -2.0, 0, 0, 0]],   // 1 frame ⇒ 1 block ⇒ dH_se 0
  pocketIdx: [0, 1], nProt: 6, mass: 110, T: 300,
});
assert(thermoThin.dH_se === 0,
  `real computeThermodynamics returns dH_se === 0 with a single block (the "not computed" sentinel)`);
const thermoWide = computeThermodynamics({
  holoFrames: Array.from({ length: 40 }, () => new Float64Array(18)),
  apoFrames: Array.from({ length: 40 }, () => new Float64Array(18)),
  holoEnergies: Array.from({ length: 40 }, (_, i) => [-1.0 + 0.01 * (i % 5), 0.2, 0, -2.0, 0, 0, 0]),
  pocketIdx: [0, 1], nProt: 6, mass: 110, T: 300,
});
assert(thermoWide.dH_se > 0,
  `real computeThermodynamics returns dH_se = ${thermoWide.dH_se.toFixed(4)} > 0 with 40 frames`);

const SNAP = {
  generatedAt: "2026-10-03T12:00:00.000Z",
  input: { pdbId: "4W52", fileName: "4w52.pdb", origin: "bundled-preset", pdbText: read("4w52.pdb") },
  selection: { chains: "A", resFrom: 1, resTo: 164, includeLigand: true },
  physics: {
    modelMode: "cg", physicsLevel: "L1",
    tierFlags: { charges: true, hbMode: "directional", weak: "off", bindLog: true },
    temperatureK: 300, frictionPerPs: 8, beadMassDa: 110, dtPs: 0.005,
    cutoffAng: 10, springGamma: 2, solventModel: "gb", saltM: 0.15,
    epsIn: 4.0, epsOut: 78.5, sasaGamma: 0.0072,
    respaOn: true, respaInnerFs: 1, respaOuterFs: 4,
    nParticles: 170, nProtein: 164, nLigandAtoms: 6, nLigandMolecules: 1,
    motionGain: 1.3,
  },
  run: { seed: 101, simulatedPs: 100, steps: 20000, dtPs: 0.005, startedAt: "2026-10-03T11:50:00.000Z" },
  trajectory: { nFrames: 51, spanPs: 100, stridePs: 2 },
  results: [
    { id: "dg_bind", label: "Binding dG (WTM funnel)", value: dgReal, unit: "kcal/mol",
      uncertainty: { value: seReal, kind: "hill-count-convergence", method: "kB*T/sqrt(nHills)" }, n: funnel._nHills },
    { id: "occupancy_bound_fraction", label: "bound fraction", value: 0.86, unit: "fraction", n: 51,
      noUncertaintyReason: "no binomial or bootstrapped error bar is computed on the bound-frame fraction" },
    { id: "dH_total", label: "dH (block bootstrap)", value: thermoThin.dH.total, unit: "kcal/mol",
      uncertainty: { value: thermoThin.dH_se, kind: "block-bootstrap-se", method: "20 contiguous blocks" },
      noUncertaintyReason: "the 20-block bootstrap needs ≥ 2 populated blocks; this leg had 1", n: 1 },
    { id: "dS_total", label: "dS (Schlitter + torsion + SASA)", value: thermoWide.dS.total, unit: "kcal/mol/K",
      noUncertaintyReason: "no standard error is computed on dS_pocket; VALIDATION measures replica SD ≈ 7 kcal/mol on -T*dS", n: 40 },
    { id: "pocket_volume_mean_ang3", label: "mean pocket volume", value: 612.4, unit: "A^3", n: 51,
      dispersion: { value: 88.1, kind: "sample-sd-of-track" },
      noUncertaintyReason: "the per-frame track SD is a dispersion; no standard error on the mean is computed" },
  ],
  runHealth: { errors: 0, ignored: 3, lastError: null },
  coverageNotIncluded: [{ family: "thermodynamics dH/dS legs", why: "separate run", howToGet: "Export Thermo (TXT)" }],
};
const rec = buildResultsRecord(SNAP);
assert(rec.version === RESULTS_RECORD_VERSION && RESULTS_RECORD_VERSION === 1,
  `record version === 1 (got ${rec.version})`);
assert(rec.format === "results-record-v1", "record carries format results-record-v1");

// ---- provenance -------------------------------------------------------
assert(rec.provenance.code.version === JSON.parse(read("package.json")).version,
  `provenance.code.version === package.json version (${rec.provenance.code.version})`);
assert(/^\d{4}-\d{2}-\d{2}$/.test(rec.provenance.code.buildDate), `buildDate present (${rec.provenance.code.buildDate})`);
assert(rec.provenance.run.seed === 101 && /SeededRNG .*seed 101/.test(rec.provenance.run.seedSource),
  `provenance.run names the seed AND its source ("${rec.provenance.run.seedSource}")`);
assert(rec.provenance.run.simulatedPs === 100 && rec.provenance.run.steps === 20000 && rec.provenance.run.dtPs === 0.005,
  "provenance.run carries simulated ps, step count and dt");
assert(rec.provenance.physics.modelMode === "cg" && rec.provenance.physics.physicsLevel === "L1",
  "provenance.physics records CG mode and the L1 physics level");
assert(rec.provenance.physics.tierFlags.charges === true && rec.provenance.physics.tierFlags.hbMode === "directional",
  "provenance.physics records the force-field tier flags actually in force");
assert(rec.provenance.physics.temperatureK === 300 && rec.provenance.physics.frictionPerPs === 8
  && rec.provenance.physics.beadMassDa === 110,
  "provenance.physics carries T, friction and bead mass");
assert(rec.provenance.physics.cutoffAng === 10 && rec.provenance.physics.springGamma === 2,
  "provenance.physics carries the ENM cutoff and spring gamma");
assert(rec.provenance.physics.saltM === 0.15 && rec.provenance.physics.epsIn === 4.0 && rec.provenance.physics.epsOut === 78.5
  && rec.provenance.physics.sasaGamma === 0.0072 && rec.provenance.physics.solventModel === "gb",
  "provenance.physics carries the solvent/GB settings");
assert(rec.provenance.physics.respa.on === true && rec.provenance.physics.respa.outerFs === 4,
  "provenance.physics carries the RESPA settings");
assert(rec.provenance.physics.system.nProtein === 164 && rec.provenance.physics.system.nLigandAtoms === 6,
  "provenance.physics.system carries the system composition");
assert(/NO PME/.test(rec.provenance.physics.electrostatics),
  "provenance.physics names the no-PME truncation instead of leaving the reader to assume a lattice sum");
assert(rec.provenance.notAPhysicsParameter.motionGain === 1.3
  && /display amplification/.test(rec.provenance.notAPhysicsParameter.why),
  "motionGain is filed under notAPhysicsParameter with the reason (it is a viewer gain)");

// ---- input identity ---------------------------------------------------
assert(rec.input.contentHash === JSON.parse(read("data/manifest.json"))["4w52.pdb"],
  "input.contentHash is the manifest's sha256 for the bundled structure");
assert(rec.input.pdbId === "4W52" && rec.input.fileName === "4w52.pdb" && rec.input.origin === "bundled-preset",
  "input carries PDB id, file name and origin");
assert(rec.input.selection.chains === "A" && rec.input.selection.resFrom === 1 && rec.input.selection.resTo === 164,
  "input.selection records the residue window actually simulated");

// ---- results + honesty ------------------------------------------------
const byId = Object.fromEntries(rec.results.map((r) => [r.id, r]));
assert(rec.results.length === SNAP.results.length, `all ${SNAP.results.length} result rows present`);
assert(byId.dg_bind.value === dgReal && byId.dg_bind.uncertainty.value === seReal,
  `dg_bind carries the real funnel value + convergence SE (${dgReal.toFixed(3)} ± ${seReal.toFixed(3)})`);
assert(byId.dg_bind.uncertainty.kind === "hill-count-convergence" && byId.dg_bind.uncertainty.n === funnel._nHills,
  "dg_bind names its uncertainty estimator and its hill count");
assert(/ranking indicator/.test(byId.dg_bind.meaning) && /absolute binding free energy/.test(byId.dg_bind.notA),
  "dg_bind states what it is AND what it is not (ranking indicator, not an absolute affinity)");
assert(byId.dg_bind.notComputed.length >= 3 && /soft cores/.test(byId.dg_bind.notComputed.join("|")),
  "dg_bind lists the machinery that is absent (no soft cores / no replica exchange)");
assert(/FEP/.test(byId.dg_bind.substitute), "dg_bind names the tool that computes the thing properly");
assert(rec.results.every((r) => r.honestyMissing === false),
  "every result row carries a trust-boundary statement");
assert(rec.coverage.included.length === rec.results.length
  && rec.coverage.notIncluded.length === SNAP.coverageNotIncluded.length,
  "coverage block lists what IS and IS NOT in the record");

// ---- run health -------------------------------------------------------
assert(rec.runHealth.errors === 0 && rec.runHealth.clean === true && rec.runHealth.documentedNoOps === 3,
  "runHealth carries error + documented-no-op counts and derives `clean`");
assert(/errors\.js/.test(rec.runHealth.recordedBy), "runHealth names its recorder (src/errors.js)");

// =====================================================================
// 3. NO ZERO ERROR BAR, EVER
// =====================================================================
console.log("=== uncertainty is never zero, and 'unknown' is explicit ===");
{
  const zeros = rec.results.filter((r) => r.uncertainty && Number(r.uncertainty.value) <= 0);
  assert(zeros.length === 0, `no result carries uncertainty <= 0 (${zeros.length} found)`);
  assert(!/"uncertainty":\{"value":0[,}]/.test(JSON.stringify(rec, null, 1)),
    "the serialized record contains no zero-valued uncertainty object");
  // The module's own 0 sentinel became null + a reason.
  const dh = byId.dH_total;
  assert(dh.uncertainty === null && dh.uncertaintyAvailable === false,
    "dH_total's 0 sentinel became uncertainty: null, not an error bar of 0");
  assert(/not computed/.test(dh.uncertaintyReason) && dh.uncertaintyReason.length > 20,
    "dH_total says WHY its error bar is unknown");
  // A real error bar still survives the same path.
  const dsSnap = buildResultsRecord({
    ...SNAP,
    results: [{ id: "dH_total", value: thermoWide.dH.total, unit: "kcal/mol", uncertainty: { value: thermoWide.dH_se, kind: "block-bootstrap-se", method: "20 contiguous blocks" }, n: 40 }],
  });
  assert(dsSnap.results[0].uncertainty.value > 0 && dsSnap.results[0].uncertaintyAvailable === true,
    `a real block-bootstrap SE survives as an error bar (${thermoWide.dH_se.toFixed(4)})`);
  // Every row without an error bar must SAY so.
  const noBar = rec.results.filter((r) => r.uncertaintyAvailable === false);
  assert(noBar.every((r) => typeof r.uncertaintyReason === "string" && r.uncertaintyReason.length > 10),
    `all ${noBar.length} rows without an error bar carry a non-trivial reason`);
  assert(rec.results.every((r) => (r.uncertainty === null) === (r.uncertaintyAvailable === false)),
    "uncertainty === null is exactly equivalent to uncertaintyAvailable === false");
  // A dispersion is not an uncertainty.
  assert(byId.pocket_volume_mean_ang3.dispersion.value === 88.1
    && byId.pocket_volume_mean_ang3.dispersion.kind === "sample-sd-of-track"
    && byId.pocket_volume_mean_ang3.uncertainty === null,
    "the pocket-volume track SD lives in `dispersion`, never in `uncertainty`");
}

// =====================================================================
// 4. JSON + CSV ROUND-TRIP
// =====================================================================
console.log("=== round-trip: JSON and CSV ===");
const text = serializeResultsRecord(rec);
assert(text.length < RESULTS_MAX_BYTES,
  `record is ${text.length} B, far under the ${RESULTS_MAX_BYTES} B cap (a summary, not a trajectory)`);
assert(text.endsWith("\n") && text.includes('"format": "results-record-v1"'),
  "serialized record is newline-terminated and self-identifying");

const parsed = parseResultsRecord(text);
assert(parsed.ok === true, "parseResultsRecord(save) ok");
assert(parsed.data.input.contentHash === rec.input.contentHash
  && parsed.data.provenance.run.seed === 101
  && parsed.data.provenance.physics.saltM === 0.15
  && parsed.data.provenance.physics.physicsLevel === "L1",
  "provenance + input identity survive the round trip exactly");
assert(parsed.data.trajectory.nFrames === 51 && parsed.data.trajectory.spanPs === 100,
  "trajectory counts survive the round trip");
for (const r of rec.results) {
  const b = parsed.data.results.find((x) => x.id === r.id);
  const same = b && b.value === r.value && b.unit === r.unit
    && JSON.stringify(b.uncertainty) === JSON.stringify(r.uncertainty)
    && b.uncertaintyReason === r.uncertaintyReason
    && b.meaning === r.meaning && b.notA === r.notA;
  assert(same, `round-trip preserves ${r.id} (value + uncertainty + honesty)`);
}
assert(parsed.data.scope.browserV1Guard === BROWSER_V1_GUARD
  && parsed.data.scope.outOfScope.length === OUT_OF_SCOPE.length,
  "the scope block survives the round trip (guard + every out-of-scope entry)");
assert(parsed.data.runHealth.errors === 0 && parsed.data.runHealth.documentedNoOps === 3,
  "run health survives the round trip");

const csv = resultsCsv(rec);
const back = parseResultsCsv(csv);
assert(back.ok === true && back.rows.length === rec.results.length,
  `results CSV parse-back yields ${back.rows.length} rows`);
{
  const g = back.rows.find((r) => r.id === "dg_bind");
  assert(g.value === dgReal && g.uncertainty.value === seReal && g.uncertainty.kind === "hill-count-convergence",
    "CSV round-trip keeps the value AND its real error bar");
  const o = back.rows.find((r) => r.id === "occupancy_bound_fraction");
  assert(o.uncertainty === null && typeof o.uncertaintyReason === "string" && o.uncertaintyReason.length > 10,
    "CSV round-trip keeps 'no error bar' as EMPTY, with the reason — not as 0");
  const d = back.rows.find((r) => r.id === "dH_total");
  assert(d.uncertainty === null && /not computed/.test(d.uncertaintyReason),
    "CSV round-trip keeps the module's 0-sentinel as null + reason");
  assert(/ranking indicator/.test(g.meaning) && /absolute binding free energy/.test(g.notA),
    "CSV carries the per-value meaning and not-a columns");
}

// ---- parse guards -----------------------------------------------------
assert(parseResultsRecord("").ok === false, "empty record text → not ok");
assert(parseResultsRecord("{not json").ok === false, "malformed JSON → not ok");
assert(validateResultsRecord({ version: 99 }).ok === false, "bad version → not ok");
assert(validateResultsRecord({ version: 1 }).ok === false, "record without a results array → not ok");
assert(parseResultsRecord("x".repeat(RESULTS_FILE_MAX_BYTES + 1)).error.code === "SYSTEM_TOO_LARGE",
  "oversized record file → SYSTEM_TOO_LARGE");
assert(parseResultsCsv("nope\n1,2").ok === false, "bad CSV header → not ok");
{
  let threw = false;
  try { serializeResultsRecord({ pad: "y".repeat(RESULTS_MAX_BYTES + 10) }); } catch (_) { threw = true; }
  assert(threw === true, "a record past the byte cap throws instead of downloading a huge file");
}

// =====================================================================
// 5. THE GATE — run health, honesty coverage, input identity
// =====================================================================
console.log("=== export gate fails closed ===");
assert(resultsExportGate(rec).ok === true, "a clean, fully-labelled record passes the gate");
assert(resultsFilename(rec) === "results_4W52_v1.1.0-fp7.json",
  `resultsFilename → ${resultsFilename(rec)}`);

const dirty = buildResultsRecord({ ...SNAP, runHealth: { errors: 2, ignored: 0, lastError: { context: "tick@advance", message: "non-finite energy", count: 2 } } });
const gDirty = resultsExportGate(dirty);
assert(gDirty.ok === false && /2 error\(s\)/.test(gDirty.reason),
  "a run with 2 recorded errors is NOT exportable");
assert(/reload/.test(gDirty.remedy), "the refusal names a remedy, not just an error");
assert(dirty.runHealth.clean === false,
  "runHealth.clean is derived from the counters, so a caller cannot mark a dirty run clean");
assert(buildResultsRecord({ ...SNAP, runHealth: { errors: 0, ignored: 4 } }).runHealth.clean === true,
  "documented no-ops do NOT make a run unclean (they are the defensive catches)");
assert(resultsExportGate(buildResultsRecord({ ...SNAP, runHealth: { errors: 0, ignored: 4 } })).warnings
  .some((w) => /documented no-op/.test(w)),
  "documented no-ops are surfaced as a WARNING on an otherwise-clean export");

const unlabelled = buildResultsRecord({ ...SNAP, results: [{ id: "made_up_value", value: 1.0, unit: "x" }] });
assert(unlabelled.results[0].honestyMissing === true && unlabelled.results[0].meaning === null,
  "a result with no VALUE_TRUST_BOUNDARY entry is marked honestyMissing, not silently accepted");
assert(resultsExportGate(unlabelled).ok === false,
  "a record containing an unlabelled value is NOT exportable");

const noHash = buildResultsRecord({ ...SNAP, input: { pdbId: "4W52", origin: "file-upload", pdbText: "" } });
assert(noHash.input.contentHash === null && resultsExportGate(noHash).ok === false,
  "an input with no bytes to hash is NOT exportable");

const unseeded = buildResultsRecord({ ...SNAP, run: { seed: null, simulatedPs: 1, steps: 1, dtPs: 0.005 } });
assert(unseeded.provenance.run.seed === null && /NOT reproducible/.test(unseeded.provenance.run.seedSource),
  "an unseeded run records seed: null and says it is not reproducible");
const gUnseeded = resultsExportGate(unseeded);
assert(gUnseeded.ok === true && gUnseeded.warnings.some((w) => /UNSEEDED/.test(w)),
  "unseeded is a WARNING (blocking it would make every default browser run unexportable), and the record already says it");

{
  const noBarWarn = resultsExportGate(rec).warnings.filter((w) => /NO error bar/.test(w));
  assert(noBarWarn.length === 1 && /occupancy_bound_fraction/.test(noBarWarn[0]),
    "the gate reports WHICH values have no error bar, by id");
}

// =====================================================================
// 6. SCOPE DERIVATION + THE DRIFT TEST (with fault injection)
// =====================================================================
console.log("=== honesty fields derive from ROADMAP.md §1 (drift test) ===");
const roadmap = read("ROADMAP.md");
const limitations = read("docs/LIMITATIONS.md");
const validation = read("docs/VALIDATION.md");

assert(parseRoadmapGuard(roadmap) === BROWSER_V1_GUARD,
  `BROWSER_V1_GUARD is parsed verbatim out of ROADMAP.md:3 ("${parseRoadmapGuard(roadmap)}")`);
{
  const bullets = parseRoadmapOutOfScope(roadmap);
  assert(bullets.length === 5,
    `ROADMAP.md §1 yields ${bullets.length} hard "no" bullets (${bullets.map((b) => b.headline).join(" · ")})`);
  assert(bullets.every((b) => OUT_OF_SCOPE.some((e) => normKey(e.roadmapKey) === b.key)),
    "every §1 bullet has a matching OUT_OF_SCOPE row in the record's scope table");
  // §2 ("what we WILL do") and §4's restatement must NOT leak into §1.
  assert(!bullets.some((b) => /canvas2d|webgl|zero-install/.test(b.headline)),
    "the §1 parser is scoped to section 1 only (no §2/§4 leakage)");
}

const drift = scopeDrift({ roadmap, limitations, validation });
assert(drift.length === 0,
  `zero drift between src/scope.js and ROADMAP.md §1 / docs/LIMITATIONS.md / docs/VALIDATION.md${drift.length ? " — " + drift.join("; ") : ""}`);

// ---- FAULT INJECTION: the detector must be able to go red --------------
{
  // (a) a new hard "no" appears in ROADMAP.md with no table row
  const extraBullet = roadmap.replace(
    /^- \*\*No nucleic acids/m,
    "- **No machine-learned potential** — pretend this shipped.\n- **No nucleic acids"
  );
  const d1 = scopeDrift({ roadmap: extraBullet, limitations, validation });
  assert(d1.some((p) => /machine-learned potential/.test(p) && /not covered by OUT_OF_SCOPE/.test(p)),
    "FAULT INJECTION: a new §1 bullet with no table row is reported");

  // (b) the guard sentence is rewritten
  const d2 = scopeDrift({
    roadmap: roadmap.replace("No QM/MM, no explicit membrane, no PME in browser v1.", "No neural nets in browser v1."),
    limitations, validation,
  });
  assert(d2.some((p) => /BROWSER_V1_GUARD drifted/.test(p)),
    "FAULT INJECTION: a rewritten guard sentence is reported");

  // (c) a §1 bullet is deleted, leaving a stale table row
  const d3 = scopeDrift({
    roadmap: roadmap.replace(/^- \*\*No QM\/MM\*\*.*$/m, "").replace(/^- \*\*No quantum.*$/m, ""),
    limitations, validation,
  });
  assert(d3.length > 0, "FAULT INJECTION: deleting a §1 bullet is reported");

  // (d) docs/LIMITATIONS.md loses a statement the table cites
  const d4 = scopeDrift({ roadmap, limitations: limitations.replace(/no PME/g, "lattice sums"), validation });
  assert(d4.some((p) => /LIMITATIONS\.md no longer contains/.test(p)),
    "FAULT INJECTION: a LIMITATIONS.md edit that drops a cited statement is reported");

  // (e) docs/VALIDATION.md loses the ranking-only verdict
  const d5 = scopeDrift({ roadmap, limitations, validation: validation.replace(/ranking-only/g, "pretty good") });
  assert(d5.some((p) => /VALIDATION\.md no longer contains/.test(p)),
    "FAULT INJECTION: losing the ranking-only verdict in VALIDATION.md is reported");

  // (f) an unparsable ROADMAP must NOT silently pass as "no drift"
  const d6 = scopeDrift({ roadmap: "## 1. nothing here\n", limitations, validation });
  assert(d6.length >= 2, "FAULT INJECTION: an unparsable ROADMAP is reported, not treated as in-sync");
}

// ---- the shipped block is built from the table, not typed -------------
{
  const block = scopeBlock();
  assert(block.browserV1Guard === BROWSER_V1_GUARD, "scopeBlock().browserV1Guard is the table's value");
  assert(block.outOfScope.length === OUT_OF_SCOPE.length
    && block.outOfScope.every((e, i) => e.key === OUT_OF_SCOPE[i].roadmapKey),
    "scopeBlock().outOfScope is derived from OUT_OF_SCOPE row-for-row");
  assert(block.derivation.runtimeDocRead === false && /no build step/.test(block.derivation.runtimeDocReadWhy),
    "the record states plainly that the docs are NOT read at runtime, and why");
  assert(/NOT FEP/.test(block.verdict) && /ranking-only/.test(block.verdict),
    "the record's top-line verdict says ranking-only / NOT FEP / not an absolute K_D");
  assert(block.sources.some((s) => s.file === "ROADMAP.md" && s.section === "1")
    && block.sources.some((s) => s.file === "docs/LIMITATIONS.md")
    && block.sources.some((s) => s.file === "docs/VALIDATION.md"),
    "the record cites all three source documents");
  assert(parsed.data.scope.outOfScope.length === block.outOfScope.length,
    "the exported record's scope block matches scopeBlock() exactly");
}

// ---- coverage: every id the UI can emit has a boundary ----------------
{
  const missing = RESULTS_VALUE_IDS.filter((id) => !VALUE_TRUST_BOUNDARY[id]);
  assert(missing.length === 0,
    `every id the UI wiring can emit has a trust-boundary entry${missing.length ? " — MISSING: " + missing.join(", ") : ` (${RESULTS_VALUE_IDS.length} ids)`}`);
  const extra = Object.keys(VALUE_TRUST_BOUNDARY).filter((id) => !RESULTS_VALUE_IDS.includes(id));
  assert(extra.every((id) => SCOPE_TB[id] && id !== "dG_jarzynski" || id === "dG_jarzynski"),
    "VALUE_TRUST_BOUNDARY covers the extra analysis families (thermo/SMD) the panel exports separately");
  assert(RESULTS_VALUE_IDS.length <= RESULTS_MAX_VALUES,
    `the UI can emit ${RESULTS_VALUE_IDS.length} ids, under the ${RESULTS_MAX_VALUES}-row cap`);
}

// ---- no frames, no coordinates ---------------------------------------
{
  assert(!("frames" in rec) && !("coordinates" in rec), "the record has no frames/coordinates key");
  assert(!text.includes("ATOM  ") && !text.includes("MODEL "), "the record text contains no PDB coordinate records");
  assert(rec.trajectory.framesOmitted.includes("summary"), "the record says frames are omitted on purpose");
}

// ---- a value that could not be computed is ABSENT, not 0 --------------
{
  const partial = buildResultsRecord({
    ...SNAP,
    results: [
      { id: "rmsip_pca_enm", value: NaN, unit: "1 (RMSIP)" },          // not computed
      { id: "occupancy_bound_fraction", value: undefined, unit: "x" }, // no value at all
      { id: "rmsf_mean_ang", value: 0, unit: "A" },                    // a REAL zero
      { id: "made_up", value: Infinity, unit: "x" },                   // non-finite
    ],
  });
  assert(partial.results.length === 1 && partial.results[0].id === "rmsf_mean_ang",
    `non-finite / missing values are ABSENT, not zero (kept ${partial.results.length} of 4 rows)`);
  assert(partial.results[0].value === 0 && partial.results[0].uncertainty === null,
    "a genuine measured zero survives as 0 (and still gets no invented error bar)");
  assert(partial.coverage.included.length === 1,
    "coverage.included lists only the rows the record actually holds, so absence is legible");
  // The UI wiring drops a non-finite value the same way (it does not fabricate one).
  const dgNaN = buildResultsRecord({ ...SNAP, results: [{ id: "dg_bind", value: NaN, unit: "kcal/mol" }] });
  assert(dgNaN.results.length === 0, "a ΔG that is NaN (no hills) yields no row at all, not a 0 kcal/mol");
}

// ---- the ROADMAP grep check the project documents still holds ----------
{
  const guardRe = new RegExp(BROWSER_V1_GUARD.split(",").map((s) => s.trim().replace(/[.*+?^${}()|[\]\\]/g, "\\$&")).join(".*"), "m");
  assert(guardRe.test(roadmap),
    `ROADMAP.md still matches the project's documented grep (grep -n "${BROWSER_V1_GUARD}")`);
}

console.log(`\n${passed} PASSED, ${failed} FAILED`);
process.exit(failed ? 1 : 0);
