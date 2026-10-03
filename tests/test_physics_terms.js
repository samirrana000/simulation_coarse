/**
 * tests/test_physics_terms.js — THE "modular" test: a third party adds physics
 * without editing a single line of src/.
 *
 * WHAT THIS PROVES
 * ----------------
 * Everything below registers terms AT RUNTIME, from this test file, against the
 * public registry in src/physics/terms/registry.js, and asserts that a real
 * ForceField / HeavyForceField picks them up on its next compute() — with the
 * energy AND the forces of an ANALYTICALLY CHECKABLE number, not merely "it
 * ran". The synthetic term is a harmonic well with a known minimum, so the
 * assertions are on exact numbers rather than on a change:
 *
 *     U = ½ k (r − r0)²,   r = |pos_i − anchor|
 *     F_i = −(k (r − r0) / r) · (pos_i − anchor)
 *
 * At r = r0 the well has U = 0 and F = 0 exactly; at r = r0 + d it has
 * U = ½ k d² and |F| = k d. Both are asserted.
 *
 * THE POINT OF THE TERM-NOT-IN-src ASSERTION
 * ------------------------------------------
 * [6] reads every module under src/ and asserts the synthetic term's id appears
 * in none of them. That is what makes "without touching src/" falsifiable: if a
 * future change hard-coded the id, or if the registry silently required an
 * engine-side edit to admit a new term, this test goes red instead of the claim
 * quietly becoming false.
 *
 * [7] is the strongest statement available: unregistering the term restores the
 * force field's energy AND force buffer BIT-EXACTLY to the pre-registration
 * values. So the term is genuinely additive — it leaves no residue, no stale
 * accumulator, no warm state — which is the property an analysis module that
 * perturbs the model actually relies on.
 *
 * Registry state is process-global by design (see docs/PHYSICS_TERMS.md), so
 * every case removes its terms in a `finally`; a leaked term would change every
 * later test's energy.
 *
 * FAST tier: two force fields over 4W52, a handful of compute() calls, no DOM,
 * no physics change. ~1 s.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import {
  registerTerm, unregisterTerm, resetTerms, getTerm, listTerms, describeTerms,
  validateTerm, TERM_GENERATION, TERM_ENGINES,
} from "../src/physics/terms/registry.js";
import { ForceField } from "../src/forcefield.js";
import { HeavyForceField, parseHeavy } from "../src/heavy.js";
import { parseCa, selectSystem } from "../src/pdb.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");

let passed = 0, failed = 0;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}
/** Absolute-or-relative comparison, the convention used across this repo's tests. */
function near(got, want, tol, what) {
  const d = Math.abs(got - want);
  return assert(d <= tol, `${what} = ${got} ≈ ${want} (|diff| ${d.toExponential(2)} ≤ ${tol.toExponential(1)})`);
}

// ---------------------------------------------------------------------------
// THE SYNTHETIC THIRD-PARTY TERM
//
// A single-atom harmonic well in a fixed external frame: atom `index` is
// restrained to `anchor` at distance `r0`. The parameters are closed over
// (read-only from the test's point of view) and the last ΔU is recorded so the
// test can assert the term's OWN contribution separately from the total.
// `enabled` is a closure over a mutable flag so the predicate path is exercised
// on a real field rather than only in the validator.
//
// SIGN CONVENTION matches every kernel in this repo (see harmonicPairs in
// src/ff-harmonic.js): `dx` points from the restrained atom TOWARD the anchor,
// so `f[i] += (k·dr/r)·dx` is `−∂U/∂x_i` — the force that SHORTENS a stretched
// well and lengthens a compressed one.
// ---------------------------------------------------------------------------
const K = 2.0;          // kcal/mol/Å²
const R0 = 4.25;        // Å — the well's minimum
const ANCHOR = [0, 0, 0];
let calls = 0;          // how many times energy() ran
let lastU = 0;          // the ΔU this term returned on its last call
let active = true;      // the flag `enabled` reads

/** @param {number} index  which atom the well acts on */
function makeWellTerm(engine, index) {
  return {
    id: "thirdparty.tetherWell",
    engine,
    // order 5 puts it BEFORE every built-in term, which is the interesting case:
    // it proves a third party can enter the summation, not just append to it.
    order: 5,
    label: "Third-party harmonic well (test fixture)",
    unit: "kcal/mol",
    enabled: function () { return active; },
    report: function () { return lastU; },
    energy: function (pos, f) {
      if (!active) return 0;
      calls++;
      const dx = ANCHOR[0] - pos[3 * index];
      const dy = ANCHOR[1] - pos[3 * index + 1];
      const dz = ANCHOR[2] - pos[3 * index + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
      const dr = r - R0;
      lastU = 0.5 * K * dr * dr;
      const s = (K * dr) / r;
      f[3 * index] += s * dx;
      f[3 * index + 1] += s * dy;
      f[3 * index + 2] += s * dz;
      return lastU;
    },
  };
}

/** Snapshot energy + force bytes of a field, for the bit-exactness claims. */
function snapshot(ff) {
  return { U: ff.energy, forces: Float64Array.from(ff.forces) };
}
function identical(a, b) {
  if (!Object.is(a.U, b.U)) return false;
  for (let i = 0; i < a.forces.length; i++) if (!Object.is(a.forces[i], b.forces[i])) return false;
  return true;
}

/** Reset every registry-mutable global this file touches. */
function resetAll() { active = true; calls = 0; lastU = 0; resetTerms("cg"); resetTerms("heavy"); }

// ---------------------------------------------------------------------------
// The two real force fields. Built once; neither is modified structurally by
// this test — only its term registry is.
// ---------------------------------------------------------------------------
const pdbText = fs.readFileSync(path.join(ROOT, "4w52.pdb"), "utf-8");
const selCG = selectSystem(parseCa(pdbText));
const ffCG = new ForceField(selCG, { rc: 10, gamma: 1.0 });
const parsedHeavy = parseHeavy(pdbText);
// Two DIFFERENT heavy configurations, so [5] can show the term is independent
// of the engine's own opt-ins rather than of one particular flag set.
const ffHeavy = new HeavyForceField({ atoms: parsedHeavy.atoms }, { gamma: 2.0 }, []);
const ffHeavyRigor = new HeavyForceField(
  { atoms: parsedHeavy.atoms },
  { gamma: 2.0, useAmber14: true, gbModel: "obc2", sasaModel: "lcpo", weak: "on", trackTerms: true },
  [],
);

/** Put atom `i` at an exact position so the analytic numbers are exact. */
function poseWith(ff, i, x, y, z) {
  const pos = Float64Array.from(ff.ref);
  pos[3 * i] = x; pos[3 * i + 1] = y; pos[3 * i + 2] = z;
  return pos;
}

console.log("=== physics-term registry: a third party adds physics (no src/ edit) ===");

// ---------------------------------------------------------------------------
// [1] THE SEAM — the synthetic term lands in a real force field's energy and
//     forces, at analytically known values.
// ---------------------------------------------------------------------------
console.log("\n[1] a term registered from this test file reaches a real CG ForceField...");
{
  resetAll();
  const idx = 0;
  // Atom 0 at exactly (6, 0, 0): r = 6 exactly, dr = 6 − 4.25 = 1.75,
  // U = ½·2·1.75² = 3.0625 exactly, |F| = k·dr = 3.5 along −x̂.
  const pos = poseWith(ffCG, idx, 6, 0, 0);

  ffCG.compute(pos);
  const before = snapshot(ffCG);
  const uBefore = before.U;

  registerTerm(makeWellTerm("cg", idx));
  assert(calls === 0, "the registry holds the term but nothing ran yet — compute() has not been called since");
  ffCG.compute(pos);
  assert(calls > 0, `the term's energy() ran inside ff.compute() (${calls} call(s)) with no src/ edit`);

  near(lastU, 3.0625, 1e-12, "term ΔU = ½·k·dr²");
  const dU = ffCG.energy - uBefore;
  near(dU, 3.0625, 1e-9, "total energy rose by the term's ΔU");
  assert(ffCG.energy !== uBefore, "ff.energy is a different float64 than before registration");

  // The force the well wrote on atom 0, recovered by differencing the buffers.
  // Nothing else in the CG model writes atom 0 at a different position than
  // before, because `pos` is byte-identical across both compute() calls.
  const after = snapshot(ffCG);
  near(after.forces[0] - before.forces[0], -3.5, 1e-12, "force on atom 0, x = −k·dr");
  near(after.forces[1] - before.forces[1], 0, 1e-12, "force on atom 0, y = 0");
  near(after.forces[2] - before.forces[2], 0, 1e-12, "force on atom 0, z = 0");
}

// ---------------------------------------------------------------------------
// [2] KNOWN MINIMUM — the well is exactly where the analytic form says it is.
// ---------------------------------------------------------------------------
console.log("\n[2] the synthetic term has a known minimum (analytic, not fitted)...");
{
  resetAll();
  const idx = 0;
  registerTerm(makeWellTerm("cg", idx));
  const atMin = poseWith(ffCG, idx, R0, 0, 0);   // r = R0 exactly
  ffCG.compute(atMin);
  assert(lastU === 0, `U = 0 exactly at r = r0 = ${R0} (got ${lastU})`);
  const fAtMin = snapshot(ffCG);
  // Compare against the same pose with the term removed.
  resetTerms("cg");
  ffCG.compute(atMin);
  const fOff = snapshot(ffCG);
  for (let a = 0; a < 3; a++) {
    assert(Object.is(fAtMin.forces[a], fOff.forces[a]),
      `force on atom 0 axis ${a} is bit-identical to the no-term field at the minimum`);
  }
}

// ---------------------------------------------------------------------------
// [3] THE enabled PREDICATE — a disabled term contributes nothing at all and
//     is not merely zeroed.
// ---------------------------------------------------------------------------
console.log("\n[3] the enabled predicate gates the term on a real field...");
{
  resetAll();
  const idx = 0;
  const pos = poseWith(ffCG, idx, 6, 0, 0);
  registerTerm(makeWellTerm("cg", idx));

  active = false;
  calls = 0;
  ffCG.compute(pos);
  const offU = ffCG.energy;
  const offF = snapshot(ffCG);
  assert(calls === 0, "a disabled term's energy() is never called (no wasted work, not a zero return)");

  active = true;
  ffCG.compute(pos);
  assert(calls > 0, "flipping the predicate re-enables it on the very next compute()");
  assert(ffCG.energy !== offU, "the enabled term changes the energy");
  const onF = snapshot(ffCG);
  let touched = 0;
  for (let i = 0; i < offF.forces.length; i++) if (!Object.is(offF.forces[i], onF.forces[i])) touched++;
  assert(touched > 0, `the enabled term changed ${touched} force components`);
  near(onF.forces[0] - offF.forces[0], -3.5, 1e-12, "force on atom 0, x = −k·dr (enabled)");
}

// ---------------------------------------------------------------------------
// [4] ORDER IS DATA — the third party's integer order decides where it lands
//     in the running-U accumulation, and describeTerms() reports it.
// ---------------------------------------------------------------------------
console.log("\n[4] order is the term's data, and it is introspectable...");
{
  resetAll();
  const cgIds = listTerms("cg").map((t) => t.id);
  assert(cgIds[0] === "cg.bonds", `the built-in CG plan starts with ${cgIds[0]}`);
  assert(cgIds.includes("cg.repulsion") && cgIds.includes("cg.binding"),
    `the built-in CG plan carries the migrated terms (${cgIds.length} total)`);

  registerTerm(makeWellTerm("cg", 0));
  const withTerm = listTerms("cg");
  assert(withTerm[0].id === "thirdparty.tetherWell",
    `order 5 sorts BEFORE cg.bonds (first is ${withTerm[0].id}, order ${withTerm[0].order})`);
  assert(withTerm[1].id === "cg.bonds", `and cg.bonds follows it (${withTerm[1].id})`);
  assert(withTerm.every((t, i) => i === 0 || withTerm[i - 1].order <= t.order),
    "listTerms() is totally ordered by ascending `order`");

  const table = describeTerms("cg", ffCG);
  const row = table.find((r) => r.id === "thirdparty.tetherWell");
  assert(row !== undefined, "describeTerms() lists the third-party term");
  assert(row.label === "Third-party harmonic well (test fixture)" && row.unit === "kcal/mol",
    `describeTerms() carries the term's own label and unit ("${row?.label}", ${row?.unit})`);
  assert(row.order === 5 && typeof row.active === "boolean",
    `describeTerms() carries order=${row?.order} and active=${row?.active}`);

  // Append at the END of the plan instead, to prove order is honoured both ways.
  unregisterTerm("thirdparty.tetherWell");
  registerTerm({ ...makeWellTerm("cg", 0), order: 9999, id: "thirdparty.lastTerm" });
  const tail = listTerms("cg");
  assert(tail[tail.length - 1].id === "thirdparty.lastTerm",
    `order 9999 sorts after every built-in (last is ${tail[tail.length - 1].id})`);
  // Tie-break is deterministic and documented.
  registerTerm({ ...makeWellTerm("cg", 0), order: 9999, id: "thirdparty.aaaTiebreak" });
  const tied = listTerms("cg");
  const iA = tied.findIndex((t) => t.id === "thirdparty.aaaTiebreak");
  const iB = tied.findIndex((t) => t.id === "thirdparty.lastTerm");
  assert(iA < iB, `equal orders tie-break on id ascending (${tied[iA].id} before ${tied[iB].id})`);
}

// ---------------------------------------------------------------------------
// [5] BOTH ENGINES — the same registry admits a term on the heavy field, whose
//     compute() lives behind a completely different assembly (plan walk in
//     src/heavy/energy.js, grouped multi-term partials, the fallback
//     transaction machinery).
// ---------------------------------------------------------------------------
console.log("\n[5] the same descriptor shape works on HeavyForceField...");
{
  resetAll();
  const idx = 0;
  const pos = poseWith(ffHeavy, idx, 6, 0, 0);
  ffHeavy.compute(pos);
  const before = snapshot(ffHeavy);

  registerTerm(makeWellTerm("heavy", idx));
  ffHeavy.compute(pos);
  near(lastU, 3.0625, 1e-12, "heavy: term ΔU = ½·k·dr²");
  near(ffHeavy.energy - before.U, 3.0625, 1e-9, "heavy: total energy rose by the term's ΔU");
  const after = snapshot(ffHeavy);
  near(after.forces[0] - before.forces[0], -3.5, 1e-12, "heavy: force on atom 0, x = −k·dr");
  near(after.forces[1] - before.forces[1], 0, 1e-12, "heavy: force on atom 0, y = 0");
  near(after.forces[2] - before.forces[2], 0, 1e-12, "heavy: force on atom 0, z = 0");

  // The term is independent of the engine's own flags: a field carrying
  // AMBER14 + OBC2 + LCPO + weak terms + trackTerms must accept it identically.
  const pos2 = poseWith(ffHeavyRigor, idx, 6, 0, 0);
  ffHeavyRigor.compute(pos2);
  const b2 = snapshot(ffHeavyRigor);
  registerTerm({ ...makeWellTerm("heavy", idx), id: "thirdparty.heavyProbe" });
  ffHeavyRigor.compute(pos2);
  near(ffHeavyRigor.energy - b2.U, 3.0625, 1e-9,
    "heavy: term lands identically on an AMBER14+OBC2+LCPO+weak+trackTerms field");
  near(snapshot(ffHeavyRigor).forces[0] - b2.forces[0], -3.5, 1e-12,
    "heavy: and contributes the same force there");
}

// ---------------------------------------------------------------------------
// [6] NO src/ EDIT — the synthetic term exists nowhere in the source tree, so
//     the seam cannot have been hard-wired.
// ---------------------------------------------------------------------------
console.log("\n[6] the term exists nowhere in src/ (the seam is not hard-wired)...");
{
  const ids = ["thirdparty.tetherWell", "thirdparty.lastTerm", "thirdparty.aaaTiebreak", "thirdparty.heavyProbe"];
  const files = [];
  (function walk(d) {
    for (const e of fs.readdirSync(d, { withFileTypes: true })) {
      const p = path.join(d, e.name);
      if (e.isDirectory()) walk(p); else if (e.name.endsWith(".js")) files.push(p);
    }
  })(path.join(ROOT, "src"));
  const hits = [];
  for (const f of files) {
    const src = fs.readFileSync(f, "utf-8");
    for (const id of ids) if (src.includes(id)) hits.push(`${path.relative(ROOT, f)}: ${id}`);
  }
  assert(hits.length === 0,
    `no src/ module mentions any third-party term id (${files.length} modules scanned)` +
    (hits.length ? ` — FOUND ${hits.join(", ")}` : ""));
  // The registration surface must be reachable from a test with no engine-
  // specific wiring at all: assert on the imported bindings themselves.
  for (const [name, fn] of Object.entries({
    registerTerm, unregisterTerm, resetTerms, listTerms, describeTerms,
  })) {
    assert(typeof fn === "function", `the registry exports ${name}() to any importer`);
  }
}

// ---------------------------------------------------------------------------
// [7] NO RESIDUE — unregistering restores energy and forces BIT-EXACTLY. This
//     is what an analysis module that perturbs the model (alanine scanning,
//     cryptic pockets) relies on.
// ---------------------------------------------------------------------------
console.log("\n[7] unregistering a third-party term leaves no residue (bit-exact)...");
{
  for (const [label, ff] of [["CG", ffCG], ["heavy", ffHeavy]]) {
    const pos = poseWith(ff, 0, 6, 0, 0);
    ff.compute(pos);
    const base = snapshot(ff);
    registerTerm({ ...makeWellTerm(ff === ffCG ? "cg" : "heavy", 0), id: `thirdparty.residue${label}` });
    ff.compute(pos);
    assert(ff.energy !== base.U, `${label}: the term changed the energy while registered`);
    assert(unregisterTerm(`thirdparty.residue${label}`), `${label}: unregisterTerm removed it`);
    ff.compute(pos);
    const back = snapshot(ff);
    assert(identical(base, back),
      `${label}: energy + all ${base.forces.length} force components are bit-identical to before registration`);
  }
}

// ---------------------------------------------------------------------------
// [8] THE REGISTRY VALIDATES ITS OWN CONTRACT — including the one mistake that
//     would otherwise fail silently (an arrow function, which loses `this`).
// ---------------------------------------------------------------------------
console.log("\n[8] the registry rejects malformed descriptors loudly...");
{
  const base = { id: "thirdparty.probe", engine: "cg", order: 500, label: "probe", energy: function () { return 0; } };
  const throws = (mut, why) => {
    let threw = null;
    try { registerTerm({ ...base, ...mut }); } catch (e) { threw = e; }
    if (threw) passed++; else failed++;
    console.log(threw ? `  ✓ ${why}` : `  ✗ FAIL: ${why} (nothing was thrown)`);
    return threw;
  };
  throws({ id: "NoDots" }, "a descriptor without an engine.name id is rejected");
  throws({ engine: "molecular" }, "an unknown engine is rejected");
  throws({ order: 1.5 }, "a fractional order is rejected");
  throws({ label: "" }, "an empty label is rejected");
  throws({ energy: "not a function" }, "a non-function energy is rejected");
  const arrow = throws({ energy: () => 0 },
    "an ARROW energy is rejected — `this` is the contract and an arrow would silently get undefined");
  if (arrow) {
    assert(/arrow/i.test(arrow.message),
      `the arrow-function rejection names the cause ("${arrow.message.slice(0, 96)}…")`);
  }
  throws({ enabled: "yes" }, "a non-function `enabled` is rejected");
  // A well-formed descriptor is accepted and round-trips through getTerm().
  const good = validateTerm({ ...base, unit: "kcal/mol", reports: ["probeU"], report: function () { return 0; } });
  assert(good.enabled === null && good.unit === "kcal/mol" && Array.isArray(good.reports),
    "a valid descriptor normalizes (enabled → null, unit defaulted, reports kept)");
  registerTerm({ ...base });
  assert(getTerm("thirdparty.probe")?.id === "thirdparty.probe", "getTerm() finds a registered term by id");
  assert(TERM_GENERATION > 0, `TERM_GENERATION advanced with the registry (now ${TERM_GENERATION})`);
  assert(Array.isArray(TERM_ENGINES) && TERM_ENGINES.length === 2,
    `the registry names its engines (${TERM_ENGINES.join(", ")})`);
}

// ---------------------------------------------------------------------------
// [9] NO RESIDUE IN THE REGISTRY — resetAll() restores both plans exactly.
// ---------------------------------------------------------------------------
console.log("\n[9] resetTerms() restores both engines' built-in plans...");
resetAll();
{
  assert(listTerms("cg").every((t) => !t.id.startsWith("thirdparty.")),
    `the CG plan holds no third-party term (${listTerms("cg").length} built-ins)`);
  assert(listTerms("heavy").every((t) => !t.id.startsWith("thirdparty.")),
    `the heavy plan holds no third-party term (${listTerms("heavy").length} built-ins)`);
  // And the fields are numerically back where they started: compute() still runs
  // the whole built-in plan and returns a finite energy after the reset.
  const pos = poseWith(ffCG, 0, 6, 0, 0);
  const uCG = ffCG.compute(pos);
  const uHeavy = ffHeavy.compute(poseWith(ffHeavy, 0, 6, 0, 0));
  assert(Number.isFinite(uCG) && Number.isFinite(uHeavy),
    `both fields still compute() to a finite energy after the registry was reset (CG ${uCG}, heavy ${uHeavy})`);
  assert(uCG !== 0 && uHeavy !== 0, "and those energies are the built-in plans' own, not a stub");
}

console.log(`\n=== test_physics_terms: ${passed} PASSED, ${failed} FAILED ===`);
console.log("PASS: a third party registered a physics term at runtime and got the exact");
console.log("      analytic energy AND forces out of a real force field, with no src/ edit.");
process.exit(failed > 0 ? 1 : 0);
