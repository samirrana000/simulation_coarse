/**
 * tests/test_heavy_obc2_trackterms.js — OBC2 x per-term-tracking regression.
 *
 * Run: node tests/test_heavy_obc2_trackterms.js   (also in tests/suites.js, FAST)
 *
 * THE BUG THIS PINS
 * -----------------
 * src/heavy/nonbonded.js `_nonBondedGridNoGB` read `ligStart` in its `trk`
 * initializer five lines BEFORE the `const nProt = ..., ligStart = ...`
 * declaration that introduced it — a temporal-dead-zone reference. Because
 *   const trk = this.trackTerms === true && ligStart > 0
 * short-circuits when trackTerms is false, the ReferenceError was latent: the
 * default configuration never noticed. Turn on `gbModel: "obc2"` AND
 * `trackTerms: true` and it threw, compute() caught it, and the catch ran the
 * HCT kernel on TOP of forces the OBC2 kernel had already merged. The field
 * then reported pure-HCT ENERGIES with OBC2-GB + HCT-everything FORCES: a
 * force field that was neither model, with energy and gradient disagreeing.
 *
 * WHY THE HEADLINE IS A PROPERTY, NOT A GOLDEN BLOB
 * -------------------------------------------------
 * The reproduction configuration cannot be compared against a golden, because
 * the golden was captured from the buggy tree. So the real evidence is
 * structural: per-term tracking is an accounting switch, not a physics switch.
 * With tracking ON and OFF, the SAME system at the SAME poses must produce
 * BIT-IDENTICAL energy and BIT-IDENTICAL forces — tracking only decides whether
 * bindLJU/bindCoulU/bindHBU/desolvU get filled in, and it must not perturb a
 * single bit of anything else. Before the fix they DISAGREED (27.566792 vs
 * 27.530422 max|F|; one of the two was an OBC2+HCT hybrid). After the fix they
 * agree exactly. That is a property of the engine, checkable by anyone, and it
 * cannot be satisfied by a golden blob.
 *
 * COVERAGE (headless, 4W52 + benzene, no DOM)
 *   [0] the reproducing configuration completes with NO substitution:
 *       ff.physicsFallbacks === 0, ff.lastPhysicsFallback === null
 *   [1] HEADLINE: trackTerms ON vs OFF agree bit-for-bit (energy, forces and
 *       every non-tracking accumulator) across 12 deterministic poses
 *   [2] the answer is really OBC2, not HCT-in-disguise: U and forces both
 *       differ from a pure-HCT field of the same configuration
 *   [3] the substitute path is transactional and LOUD, proven by fault
 *       injection: a kernel that mutates `f` and then throws must leave forces
 *       bit-identical to the clean substitute, and must announce itself
 *   [4] every substitutable catch in src/heavy/ is bracketed by
 *       _forceSnapshot()/_forceRestore() and announced via
 *       _notePhysicsFallback() — a structural pin, so the defect class cannot
 *       be reintroduced in a new stage without turning this red
 *
 * FAST tier: 3 fields x 12 poses x 2 configurations + 2 fault-injected computes.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { parseHeavy, selectHeavy, appendHeavyLigands, HeavyForceField } from "../src/heavy.js";
import { parseMol2 } from "../src/mol2.js";
import { SeededRNG } from "../src/seeded-rng.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");

/** The exact configuration that reproduced the bug. */
const OBC2_PAR = { gamma: 2.0, temp: 300, gbModel: "obc2" };
const HCT_PAR = { gamma: 2.0, temp: 300 };
/** 12 poses: native + 11 SeededRNG(1234) jitters of +/- 0.05 A (same recipe as the golden). */
const N_POSES = 12;

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

function buildSel() {
  const parsed = parseHeavy(fs.readFileSync(path.join(ROOT, "4w52.pdb"), "utf-8"));
  const mols = parseMol2(fs.readFileSync(path.join(ROOT, "benzene.mol2"), "utf-8"));
  let sel = selectHeavy(parsed, {
    heteroSelection: { "A|200|BNZ": false, "A|201|EPE": false },
    includePdbLigands: true,
    hasExternalLigand: true,
  });
  sel = appendHeavyLigands(sel, mols);
  return sel;
}

function poses(ref, count = N_POSES) {
  const out = [new Float64Array(ref)];
  const rng = new SeededRNG(1234);
  for (let p = 1; p < count; p++) {
    const q = new Float64Array(ref);
    for (let i = 0; i < q.length; i++) q[i] += (rng.rand() * 2 - 1) * 0.05;
    out.push(q);
  }
  return out;
}

/** Accumulators that must be identical whether or not tracking is on. */
const PHYSICS_TERMS = [
  "repU", "bondU", "angleU", "improperU", "properU", "coordU", "elecU", "gbU",
  "sasaU", "hbondU", "membraneU", "springU", "weakU", "bindingU",
];

/** Captures console.warn so a fallback can be asserted rather than merely printed. */
function withWarnCapture(fn) {
  const seen = [];
  const real = console.warn;
  console.warn = (...a) => { seen.push(a.join(" ")); };
  try { return { value: fn(), warns: seen }; }
  finally { console.warn = real; }
}

const maxAbs = (a) => { let m = 0; for (const x of a) m = Math.max(m, Math.abs(x)); return m; };
const sameArr = (a, b) => a.length === b.length && a.every((x, i) => Object.is(x, b[i]));

function main() {
  console.log("=== OBC2 x per-term tracking: no silent HCT substitution (4W52 heavy) ===");
  const sel = buildSel();
  const mk = (par, track) => {
    const ff = new HeavyForceField({ atoms: sel.atoms }, { ...par }, []);
    ff.trackTerms = track === true;
    return ff;
  };

  // -----------------------------------------------------------------
  // [0] the reproducing configuration no longer substitutes a force field
  // -----------------------------------------------------------------
  console.log("\n[0] the reproducing configuration completes with NO substitution...");
  const ffOn = mk(OBC2_PAR, true);
  const ps = poses(ffOn.ref);
  const cap = withWarnCapture(() => {
    for (const p of ps) ffOn.compute(p);
    return null;
  });
  assert(cap.warns.length === 0,
    `obc2 + trackTerms:true emits no fallback warning over ${ps.length} poses (got ${cap.warns.length})`);
  assert((ffOn.physicsFallbacks ?? 0) === 0,
    `ff.physicsFallbacks === 0 (got ${ffOn.physicsFallbacks ?? 0})`);
  assert(ffOn.lastPhysicsFallback === null || ffOn.lastPhysicsFallback === undefined,
    `ff.lastPhysicsFallback unset (got ${JSON.stringify(ffOn.lastPhysicsFallback ?? null)})`);
  assert(ffOn._obc2Radii !== null && ffOn._obc2Radii.length === ffOn.n,
    `OBC-II Born radii published (${ffOn._obc2Radii ? ffOn._obc2Radii.length : "null"})`);
  const sumBind = ffOn.bindLJU + ffOn.bindCoulU + ffOn.bindHBU + ffOn.desolvU;
  assert(Number.isFinite(sumBind),
    `per-term accumulators populated while tracking ON (lj=${ffOn.bindLJU} coul=${ffOn.bindCoulU} hb=${ffOn.bindHBU})`);

  // -----------------------------------------------------------------
  // [1] HEADLINE: tracking ON vs OFF agree bit-for-bit
  // -----------------------------------------------------------------
  console.log("\n[1] HEADLINE: per-term tracking ON vs OFF must not change the physics...");
  const ffOff = mk(OBC2_PAR, false);
  let uDiff = 0, fDiff = 0, tDiff = [];
  for (let p = 0; p < ps.length; p++) {
    const Uon = ffOn.compute(ps[p]);
    const fOn = Float64Array.from(ffOn.forces);
    const termsOn = Object.fromEntries(PHYSICS_TERMS.map((t) => [t, ffOn[t]]));
    const Uoff = ffOff.compute(ps[p]);
    if (!Object.is(Uon, Uoff)) uDiff++;
    if (!sameArr(fOn, ffOff.forces)) fDiff++;
    for (const t of PHYSICS_TERMS) if (!Object.is(termsOn[t], ffOff[t])) tDiff.push(`p${p}.${t}`);
  }
  assert(uDiff === 0,
    `energy bit-identical ON vs OFF at all ${ps.length} poses (${uDiff} differ)`);
  assert(fDiff === 0,
    `FORCES bit-identical ON vs OFF at all ${ps.length} poses (${fDiff} differ)`);
  assert(tDiff.length === 0,
    `all ${PHYSICS_TERMS.length} physics accumulators bit-identical ON vs OFF` +
    (tDiff.length ? ` (${tDiff.length} differ: ${tDiff.slice(0, 5).join(", ")})` : ""));
  console.log(`    max|F| tracking ON  = ${maxAbs(ffOn.forces).toFixed(6)}`);
  console.log(`    max|F| tracking OFF = ${maxAbs(ffOff.forces).toFixed(6)}`);

  // -----------------------------------------------------------------
  // [2] ...and the shared answer really is OBC2, not HCT in disguise
  // -----------------------------------------------------------------
  console.log("\n[2] the shared answer is OBC2, not the HCT fallback...");
  const ffHct = mk(HCT_PAR, false);
  const Uhct = ffHct.compute(ffHct.ref);
  const Uobc2 = ffOff.compute(ffOff.ref);
  const hctF = Float64Array.from(ffHct.forces);
  const obcF = Float64Array.from(ffOff.forces);
  assert(Object.is(Uobc2, Uhct) === false,
    `OBC2 energy differs from the HCT energy (E_obc2=${Uobc2.toFixed(6)} E_hct=${Uhct.toFixed(6)})`);
  assert(sameArr(obcF, hctF) === false,
    `OBC2 forces differ from the HCT forces (maxF_obc2=${maxAbs(obcF).toFixed(6)} maxF_hct=${maxAbs(hctF).toFixed(6)})`);
  assert(ffOff.gbU !== ffHct.gbU,
    `GB accumulators differ (gbU_obc2=${ffOff.gbU.toFixed(6)} gbU_hct=${ffHct.gbU.toFixed(6)})`);

  // -----------------------------------------------------------------
  // [3] the substitute path is transactional AND loud (fault injection)
  // -----------------------------------------------------------------
  console.log("\n[3] fault injection: a throwing kernel must not leave its forces behind...");
  // Reference: what a clean pure-HCT field produces on the same system/pose.
  const refHct = mk(HCT_PAR, false);
  const Uref = refHct.compute(refHct.ref);
  const fRef = Float64Array.from(refHct.forces);

  // Injected fault: _nonBondedGridNoGB writes into the shared force buffer and
  // THEN throws — precisely the shape of the original TDZ failure, which is why
  // the original fallback produced OBC2-GB + HCT-everything forces.
  const ffBad = mk(OBC2_PAR, true);
  ffBad._nonBondedGridNoGB = function injectPartialThenThrow(pos, f) {
    for (let i = 0; i < f.length; i++) f[i] += 1e3; // partial, garbage contribution
    throw new Error("injected: _nonBondedGridNoGB died mid-pass");
  };
  const capBad = withWarnCapture(() => ffBad.compute(ffBad.ref));
  assert(Number.isFinite(capBad.value), `compute() still returns a finite energy after the fault (U=${capBad.value})`);
  assert(capBad.warns.length === 1, `exactly one warning announces the substitution (got ${capBad.warns.length})`);
  assert(/substituted/i.test(capBad.warns[0] ?? ""),
    `warning names the substituted model ("${(capBad.warns[0] ?? "").match(/substituted \w+/)?.[0] ?? "?"}")`);
  assert(/physicsFallbacks=1/.test(capBad.warns[0] ?? ""),
    `warning reports the counter (physicsFallbacks=1)`);
  assert(ffBad.physicsFallbacks === 1, `ff.physicsFallbacks === 1 (got ${ffBad.physicsFallbacks})`);
  assert(ffBad.lastPhysicsFallback?.stage === "OBC2 non-bonded"
    && ffBad.lastPhysicsFallback?.substitutedTo === "HCT non-bonded"
    && /injected/.test(ffBad.lastPhysicsFallback?.error ?? ""),
    `ff.lastPhysicsFallback records stage/model/error (${JSON.stringify(ffBad.lastPhysicsFallback)})`);
  assert(Object.is(capBad.value, Uref),
    `substituted ENERGY is exactly the clean HCT energy (Δ=${Math.abs(capBad.value - Uref).toExponential(2)})`);
  assert(sameArr(ffBad.forces, fRef),
    `substituted FORCES are bit-identical to clean HCT — partial 1e3 contribution rolled back ` +
    `(maxΔ=${maxAbs(ffBad.forces.map((x, i) => x - fRef[i])).toExponential(2)})`);
  assert(ffBad._obc2Radii === null,
    `a failed pass does not publish Born radii (got ${ffBad._obc2Radii === null ? "null" : "set"})`);
  // A second fault must count, not overwrite.
  withWarnCapture(() => ffBad.compute(ffBad.ref));
  assert(ffBad.physicsFallbacks === 2 && ffBad.lastPhysicsFallback.count === 2,
    `repeated faults accumulate (physicsFallbacks=${ffBad.physicsFallbacks}, count=${ffBad.lastPhysicsFallback?.count})`);

  // -----------------------------------------------------------------
  // [4] structural pin: every substitutable catch in src/heavy/ is guarded
  // -----------------------------------------------------------------
  console.log("\n[4] structural pin over src/heavy/ catch blocks...");
  const energySrc = fs.readFileSync(path.resolve(ROOT, "src", "heavy", "energy.js"), "utf-8");
  const catchIdx = [];
  for (let i = 0; i < energySrc.length; i++) if (energySrc.startsWith("catch (", i)) catchIdx.push(i);
  assert(catchIdx.length === 3, `src/heavy/energy.js has 3 substitutable catches (found ${catchIdx.length})`);
  const nonBondedSrc = energySrc;
  const STAGE_MARKERS = [/this\.gbModel === "obc2"/, /this\.sasaModel === "lcpo"/, /this\.membraneOpts\?\.on/];
  for (const i of catchIdx) {
    const before = nonBondedSrc.slice(Math.max(0, i - 1600), i);
    const after = nonBondedSrc.slice(i, i + 700);
    const stage = STAGE_MARKERS.map((re, k) => (re.test(before) ? k : -1)).filter((k) => k >= 0).pop();
    const label = stage === undefined ? `catch@${i}` : ["OBC2 non-bonded", "LCPO SASA", "membrane slab"][stage];
    assert(/_forceSnapshot\(\)/.test(before),
      `catch (${label}) is preceded by _forceSnapshot()`);
    assert(/_forceRestore\(/.test(after) && /_notePhysicsFallback\(/.test(after),
      `catch (${label}) restores forces AND calls _notePhysicsFallback`);
  }
  const nbSrc = fs.readFileSync(path.resolve(ROOT, "src", "heavy", "nonbonded.js"), "utf-8");
  assert(/const base = this\._nonBondedGridNoGB\(pos, f\);\s*\n\s*for \(let i = 0; i < f\.length; i\+\+\) f\[i\] \+= gbF\[i\];/.test(nbSrc),
    "_nonBondedGridOBC2 runs the short-range pass BEFORE merging gbF (transactional order)");
  const trkIdx = nbSrc.indexOf("const trk = this.trackTerms === true && ligStart > 0;", nbSrc.indexOf("_nonBondedGridNoGB"));
  const declIdx = nbSrc.indexOf("const nProt = this.nProt, ligStart = this.ligandStart;", nbSrc.indexOf("_nonBondedGridNoGB"));
  assert(declIdx > 0 && trkIdx > declIdx,
    "in _nonBondedGridNoGB, ligStart is DECLARED BEFORE the trk initializer reads it (TDZ regression)");

  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log("=================================================");
  if (failed > 0) process.exit(1);
  console.log("PASS: OBC2 x trackTerms no longer silently substitutes the HCT force field");
}

main();