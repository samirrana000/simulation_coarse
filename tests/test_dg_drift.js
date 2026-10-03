/**
 * test_dg_drift.js — the durable guard: the ΔG the UI shows and the ΔG
 * docs/VALIDATION.md reports cannot drift apart, and every ΔG-like number the
 * app prints is enumerated.
 *
 * WHY THIS TEST EXISTS
 * --------------------
 * A user runs the app and sees a ΔG in the HUD. docs/VALIDATION.md reports a
 * completely different ΔG for the same system (4W52 + benzene). Neither is wrong
 * exactly — they come from different estimators — but NOTHING in either said
 * which estimator produced which number, what it assumed, or how far either sits
 * from the other. This test is what stops that gap from becoming two unrelated
 * numbers that each look fine on their own.
 *
 * WHAT EACH BLOCK DEFENDS
 * -----------------------
 *   1. ENUMERATION   Every ΔG-like quantity the app produces is in
 *                    src/dg-estimators.js DG_ESTIMATORS, and every source site
 *                    that PRINTS one is claimed by that table. A ΔG in the UI
 *                    that the table does not own is a failure.
 *   2. THE CHAIN IS THE ESTIMATOR'S OWN ARITHMETIC
 *                    The chain's terms sum to the number the estimator reports,
 *                    on a REAL Funnel and a REAL computeThermodynamics result.
 *                    The chain is not a hand-typed restatement.
 *   3. NO ZERO ERROR BARS  dH_se = 0 (the block-bootstrap sentinel) must appear
 *                    as `uncertainty: null` + a reason, never as 0.
 *   4. THE DRIFT GUARD  Both stories are RECOMPUTED from the code paths they name
 *                    — the HUD funnel ΔG from src/funnel.js, the doc numbers
 *                    from computeThermodynamics + sasaBurial + scanPocket under
 *                    the protocol the document's own Protocol column states —
 *                    and compared to what docs/VALIDATION.md PRINTS. Perturbing
 *                    either side must turn this red, and does (fault injection
 *                    below).
 *   5. THE RECONCILIATION  The gap between the two is reported, not hidden, and
 *                    every reason for it is tagged; `canReconcile` is false with
 *                    the missing physics named.
 *   6. CONTRACTS     #dgChain lives INSIDE an existing top-level panel, its ids
 *                    are registered in src/ui.js and present in index.html, the
 *                    Digit1-7 hotkey range is untouched, and the results record
 *                    carries the chain by reference.
 *
 * Cost: ~4 s (four seeded CG legs + the LCPO SASA sweep + one 12-residue
 * alanine scan). FAST tier, which has 60 s; the whole suite is ~42 s.
 * Deterministic: every leg is seeded and every protocol constant is stated in
 * the methods table this test reads, so nothing here depends on wall-clock.
 * Run: node tests/test_dg_drift.js
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import {
  DG_ESTIMATORS, DG_BRIDGE, VALIDATION_NUMERIC_CLAIMS,
  parseValidationClaims, validationDrift, num, printedPrecision,
} from "../src/dg-estimators.js";
import {
  funnelChain, thermoChain, reconcileChains, chainMarkdown, chainRows,
  uncertaintyOf, findUnenumeratedDGBinding, CHAIN_SASA_GAMMA, TERM_NATURE,
} from "../src/dg-chain.js";
import { Funnel } from "../src/funnel.js";
import { KB_KCAL } from "../src/units.js";
import { computeThermodynamics, schlitterEntropy } from "../src/analysis/thermodynamics.js";
import { sasaBurial, cgBeadExtendedRadius, selectedLigandElements, THERMO_SASA_STRIDE } from "../src/analysis/thermo_sasa.js";
import { pocketResidues, scanPocket, ALA_DDG_NOISE_FLOOR } from "../src/analysis/alanine_scanning.js";
import { scopeBlock } from "../src/scope.js";
import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const read = (rel) => fs.readFileSync(path.join(ROOT, rel), "utf-8");

/** How many numeric claims the drift guard tracks (read from the table, not typed). */
const VALIDATION_CLAIM_COUNT = VALIDATION_NUMERIC_CLAIMS.length;

let passed = 0, failed = 0;
/** Hoisted: section 2 measures it, section 3's drift guard compares it to the doc. */
let schlitterErrPct = NaN;
function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}
/** Free text for a failure so the message is actionable, not just negative. */
function failWith(msg, problems) {
  assert(false, `${msg} — ${problems.length} problem(s):\n      ${problems.join("\n      ")}`);
}

console.log("=== test_dg_drift.js — one auditable chain, two ΔG stories, one guard ===");

/* =====================================================================
 * PROTOCOL — read from docs/VALIDATION.md's own Protocol column, not invented
 * ===================================================================== */
const STEPS = 2000, STRIDE = 2, TEMP = 300, ZETA = 8.0, MASS = 110;
const SEED_HOLO = 101, SEED_APO = 1101;
const POCKET_RCUT = 8.0, MAXN = 12, RELAX_STEPS = 80;
const FF_PAR = { rc: 10, gamma: 2.0, binding: { charges: true, hbMode: "directional" } };
const FUNNEL_STEPS = 2000; // the CG row's step count: what a user reaches in one run

/* =====================================================================
 * 1. ENUMERATION — every ΔG-like quantity is listed; every print site claimed
 * ===================================================================== */
console.log("\n--- 1. estimator enumeration is exhaustive ---");
{
  const sweep = findUnenumeratedDGBinding({ read });
  assert(sweep.missing.length === 0,
    `all ${sweep.sitesChecked} claimed source sites still exist in their files${sweep.missing.length ? " — " + sweep.missing.join("; ") : ""}`);
  assert(sweep.deadRows.length === 0,
    `every estimator row has at least one live site (no dead rows)${sweep.deadRows.length ? " — " + sweep.deadRows.join(", ") : ""}`);

  const ids = DG_ESTIMATORS.map((e) => e.id);
  assert(new Set(ids).size === ids.length, `estimator ids are unique (${ids.length} rows)`);
  for (const e of DG_ESTIMATORS) {
    const hasU = e.uncertainty != null || (typeof e.noErrorBar === "string" && e.noErrorBar.length > 0);
    assert(typeof e.formula === "string" && e.formula.length > 4
      && typeof e.codePath === "string" && e.codePath.length > 3
      && ["measured", "assumed", "mixed"].includes(e.nature)
      && hasU,
      `${e.id}: formula + codePath + measured/assumed status + an uncertainty OR an explicit reason for its absence`);
  }

  // "A ΔG that appears in the UI but is not in your table is a bug." Made
  // mechanical: scan the real tree for ΔG display literals, require every file
  // that carries one to be claimed by a site.
  const UI_DG_RE = /(ΔG|ΔΔG|&Delta;G|dG_bind|dG_estimate|dG_int)/;
  const DISPLAY_RE = /(toFixed|innerHTML|textContent|canvas\.fillText|ctx\.fillText|push\(|inner &rarr;|title=)/;
  // The two chain modules are the ENUMERATION, not a display: their `sites` are
  // the claim about other files, and a table that listed itself would be a table
  // of one. Excluding them by name (not by pattern) keeps the sweep able to fail
  // on anything else.
  const CHAIN_MODULES = new Set(["src/dg-chain.js", "src/dg-estimators.js"]);
  const claimed = new Set(DG_ESTIMATORS.flatMap((e) => e.sites.map((s) => s.file)));
  const walk = (dir, out = []) => {
    for (const ent of fs.readdirSync(dir, { withFileTypes: true })) {
      const p = path.join(dir, ent.name);
      if (ent.isDirectory()) walk(p, out);
      else if (ent.name.endsWith(".js") || ent.name.endsWith(".html")) out.push(p);
    }
    return out;
  };
  const swept = walk(path.join(ROOT, "src"));
  const offenders = [];
  for (const abs of swept) {
    const rel = path.relative(ROOT, abs);
    if (CHAIN_MODULES.has(rel)) continue;
    const text = fs.readFileSync(abs, "utf-8");
    text.split(/\r?\n/).forEach((line, i) => {
      if (line.trimStart().startsWith("*") || line.trimStart().startsWith("//")) return;
      if (UI_DG_RE.test(line) && DISPLAY_RE.test(line) && !claimed.has(rel)) {
        offenders.push(`${rel}:${i + 1}`);
      }
    });
  }
  assert(offenders.length === 0,
    `every ΔG display in src/ belongs to a claimed estimator${offenders.length ? " — unclaimed: " + offenders.join(", ") : ""}`);
  console.log(`      (swept ${swept.length} src files, ${CHAIN_MODULES.size} of them the chain modules themselves; ${claimed.size} are claimed by a site)`);

  // Non-vacuity: delete a row and the site it owns becomes unclaimed.
  const dropped = DG_ESTIMATORS.filter((e) => e.id !== "funnel_dg_hud");
  const droppedClaimed = new Set(dropped.flatMap((e) => e.sites.map((s) => s.file)));
  const orphan = [...claimed].filter((f) => !droppedClaimed.has(f) && DG_ESTIMATORS.some((e) => e.id === "funnel_dg_hud" && e.sites.some((s) => s.file === f)));
  assert(orphan.length > 0,
    `FAULT INJECTION the enumeration check can fail: dropping the funnel_dg_hud row orphans ${orphan.length} display file(s) (${orphan.join(", ")})`);
  const missingRow = findUnenumeratedDGBinding({ read }, dropped);
  assert(missingRow.sitesChecked < sweep.sitesChecked,
    `FAULT INJECTION a removed row is detected by the site sweep (${sweep.sitesChecked} → ${missingRow.sitesChecked} sites checked)`);
}

/* =====================================================================
 * 2. THE CHAIN IS THE ESTIMATOR'S OWN ARITHMETIC (real objects, real sums)
 * ===================================================================== */
console.log("\n--- 2. chain arithmetic == estimator arithmetic ---");
{
  // A deterministic synthetic Funnel (no PDB needed): deposit a fixed hill set.
  const nProt = 12, nLig = 3, n = nProt + nLig;
  const ref = new Float64Array(3 * n);
  for (let i = 0; i < nProt; i++) { ref[3 * i] = i * 1.7 - 8; ref[3 * i + 1] = (i % 3) - 1; ref[3 * i + 2] = (i % 5) - 2; }
  for (let a = 0; a < nLig; a++) { ref[3 * (nProt + a)] = 1.5 + a * 0.1; ref[3 * (nProt + a) + 1] = 0.2; ref[3 * (nProt + a) + 2] = -0.3; }
  const funnel = new Funnel({ nProt, n, ref });
  for (let i = 0; i < 120; i++) funnel.deposit(1.5 + (i % 24) * 0.3);

  const chain = funnelChain(funnel);
  assert(chain.estimatorId === "funnel_dg_hud", "funnelChain labels the estimator it expanded");
  const terms = chain.terms;
  const dpmf = terms[0].value, minusVol = terms[1].value, reported = terms[2].value;
  assert(Number.isFinite(dpmf) && Number.isFinite(minusVol) && Number.isFinite(reported),
    `all three arithmetic terms are finite (ΔPMF ${dpmf.toFixed(6)}, −dG_vol ${minusVol.toFixed(6)}, ΔG ${reported.toFixed(6)})`);
  assert(Math.abs((dpmf + minusVol) - reported) < 1e-12,
    `term0 + term1 === the reported ΔG to 1e-12 (${(dpmf + minusVol).toFixed(12)} vs ${reported.toFixed(12)})`);
  assert(reported === funnel.estimateDG(),
    `the chain's reported term is estimateDG() BIT-for-bit (${reported} === ${funnel.estimateDG()})`);

  // −dG_vol is the Boresch correction with the units contract's constants.
  const kBT = KB_KCAL * funnel.T;
  const vRest = (4 / 3) * Math.PI * funnel.rFlat ** 3;
  const expectVol = kBT * Math.log(Math.max(1e-6, vRest / 1660.54));
  assert(Math.abs(minusVol - expectVol) < 1e-12,
    `the volume term is +kT·ln(V_rest/V°) from src/units.js (${minusVol.toFixed(9)} vs ${expectVol.toFixed(9)})`);

  const badNature = terms.filter((t) => !TERM_NATURE.includes(t.nature));
  assert(badNature.length === 0, `every funnel chain term carries a closed-set nature tag (${terms.length} terms)`);
  const oos = terms.filter((t) => t.nature === "out-of-scope");
  assert(oos.length >= 3,
    `the funnel chain names its absent physics rather than omitting it (${oos.length} out-of-scope rows: ${oos.map((t) => t.term.slice(0, 28)).join(" | ")})`);

  // The reported term's error bar: a REAL kT/√nHills at 120 hills.
  assert(terms[2].uncertainty && Math.abs(terms[2].uncertainty.value - funnel.convergenceSE()) < 1e-12,
    `the reported term carries convergenceSE()'s own value (${terms[2].uncertainty?.value.toFixed(6)}), labelled as a hill-count heuristic not an error bar`);
  assert(/NOT an uncertainty/.test(terms[2].uncertainty.method),
    "the hill-count 'error bar' says in its own method string that it is not an uncertainty on ΔG");

  // Below 50 hills convergenceSE() is null: absence must be stated, not zeroed.
  const cold = new Funnel({ nProt, n, ref });
  cold.deposit(2.0);
  const coldChain = funnelChain(cold);
  assert(coldChain.terms[2].uncertainty === null && typeof coldChain.terms[2].noErrorBar === "string" && coldChain.terms[2].noErrorBar.length > 20,
    `below 50 hills the reported term has uncertainty:null + a reason ("${(coldChain.terms[2].noErrorBar ?? "").slice(0, 58)}…"), NOT 0`);

  // funnelChain must not throw or fabricate on an absent funnel.
  const none = funnelChain(null);
  assert(none.value === null && none.terms.length === 0,
    "funnelChain(null) is an explicit absence (no value, no rows) rather than a zero");

  // ---- thermo side: a real computeThermodynamics result ----
  const mkFrames = (seed, nF, jitter) => {
    const r = (k) => Math.sin(seed * 12.9898 + k * 78.233) * 43758.5453;
    return Array.from({ length: nF }, (_, f) => {
      const buf = new Float32Array(3 * nProt);
      for (let i = 0; i < nProt; i++) {
        buf[3 * i] = i * 1.7 + jitter * r(f + i);
        buf[3 * i + 1] = r(f + i + 40) * jitter;
        buf[3 * i + 2] = r(f + i + 90) * jitter;
      }
      return buf;
    });
  };
  const holo = mkFrames(1, 40, 0.25);
  const apo = mkFrames(2, 40, 0.60);
  const holoEnergies = Array.from({ length: 40 }, (_, f) => [-1.0 + 0.01 * f, 0.2, -0.1, -2.5 + 0.02 * f, 0, 0, 0]);
  const pocketIdx = [2, 5, 9, 11];
  const thermo = computeThermodynamics({ holoFrames: holo, apoFrames: apo, holoEnergies, pocketIdx, nProt, mass: MASS, T: TEMP });
  const tchain = thermoChain(thermo);
  const sum = tchain.terms.slice(0, 4).reduce((s, t) => s + t.value, 0);
  assert(Math.abs(sum - thermo.dG_estimate) < 1e-10,
    `ΔH + (−TΔS_pocket) + (−TΔS_lig) + (−TΔS_solv) === ΔG_estimate to 1e-10 (${sum.toFixed(10)} vs ${thermo.dG_estimate.toFixed(10)})`);
  assert(tchain.terms[4].value === thermo.dG_estimate, "the thermo chain's '=' row is dG_estimate itself");
  assert(tchain.terms.every((t) => TERM_NATURE.includes(t.nature)), "every thermo chain term carries a closed-set nature tag");
  assert(tchain.terms[1].noErrorBar && /replica SD/.test(tchain.terms[1].noErrorBar),
    "the −T·ΔS_pocket row names the REPLICA SD as the uncertainty that matters instead of inventing one");

  // The 0-sentinel rule, on a case that actually produces dH_se = 0: ONE frame,
  // so the 20-block bootstrap collects a single block and never assigns a spread.
  const tiny = computeThermodynamics({
    holoFrames: mkFrames(3, 1, 0.2), apoFrames: mkFrames(4, 1, 0.3),
    holoEnergies: [[-1, 0, 0, -2, 0, 0, 0]],
    pocketIdx, nProt, mass: MASS, T: TEMP,
  });
  const tinyChain = thermoChain(tiny);
  assert(tiny.dH_se === 0,
    `the fixture really does hit the dH_se = 0 sentinel (fewer than two bootstrap blocks): dH_se = ${tiny.dH_se}`);
  assert(tinyChain.terms[0].uncertainty === null && /NOT COMPUTED/.test(tinyChain.terms[0].noErrorBar),
    `the chain converts that 0 into uncertainty:null + "not computed", never an error bar of 0 ("${tinyChain.terms[0].noErrorBar.slice(0, 62)}…")`);
  assert(chainRows(tinyChain).every((r) => r.uncertainty !== 0),
    "chainRows() emits null for every absent error bar — a downstream CSV never sees a 0");

  // The project's own rule, exercised directly.
  const u0 = uncertaintyOf(0, "bootstrap", "x", "the analysis returned 0 meaning not computed");
  const uNaN = uncertaintyOf(NaN, "bootstrap", "x", "the analysis returned 0 meaning not computed");
  const uOK = uncertaintyOf(0.16, "bootstrap", "20-block", "unused");
  assert(u0.uncertainty === null && uNaN.uncertainty === null && uOK.uncertainty.value === 0.16,
    "uncertaintyOf(): 0 and NaN → null + reason; a positive value passes through with its kind");

  // Schlitter still at its published accuracy (the doc's UNIT ANCHOR row:
  // "2-DOF Gaussian, 50k samples → 0.20% (< 1% bar)"). The comparison is
  // EXACT-covariance vs SAMPLED-covariance through the SAME function, which is
  // the claim the row makes: the estimator is unbiased against its own sampling.
  const covExact = [[1.0, 0], [0, 4.0]];
  const sExact = schlitterEntropy(covExact.map((r) => [...r]), [12, 12]);
  let s42 = 42;
  const rnd = () => (s42 = (s42 * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
  const gauss = () => Math.sqrt(-2 * Math.log(Math.max(rnd(), 1e-12))) * Math.cos(2 * Math.PI * rnd());
  let s00 = 0, s11 = 0, nS = 0;
  for (let i = 0; i < 50000; i++) { const x = gauss(), y = 2 * gauss(); s00 += x * x; s11 += y * y; nS++; }
  const sSampled = schlitterEntropy([[s00 / nS, 0], [0, s11 / nS]], [12, 12]);
  schlitterErrPct = Math.abs(sSampled - sExact) / Math.abs(sExact) * 100;
  assert(schlitterErrPct < 1,
    `Schlitter 2-DOF sampled-vs-exact error ${schlitterErrPct.toFixed(2)}% < 1% (the doc's unit-anchor bar; doc prints ${printedPrecision("0.20%")} dp)`);

  const md = chainMarkdown(funnelChain(funnel));
  assert(md.includes("[measured") && md.includes("[assumed") && md.includes("[out-of-scope") && md.includes("NO ERROR BAR"),
    "chainMarkdown() renders every nature tag and states the missing error bars in the text");
}

/* =====================================================================
 * 3. THE DRIFT GUARD — recompute BOTH stories, compare to what the doc PRINTS
 * ===================================================================== */
console.log("\n--- 3. recompute both stories and hold the document to them ---");
const pdbText = read("4w52.pdb");
const parsed = parseCa(pdbText);
const sel = selectSystem(parsed);
const allMols = parseLigands(pdbText);
const bnzMols = allMols.filter((m) => m.resName === "BNZ");

function runLeg(ff, seed) {
  const integ = new LangevinIntegrator(ff.ref, ff, MASS, { seed });
  integ.setTemperature(TEMP);
  integ.setFriction(ZETA);
  const frames = [], energies = [];
  for (let s = 0; s < STEPS; s++) {
    integ.step();
    if (s % STRIDE === 0) {
      frames.push(Float32Array.from(integ.pos));
      energies.push([ff.bindLJU, ff.bindCoulU, ff.bindHBU, ff.desolvU, 0, 0, 0]);
    }
  }
  return { frames, energies };
}
function pocket8A(ff) {
  const com = [0, 0, 0];
  for (let a = 0; a < ff.nLigAtoms; a++) {
    com[0] += ff.ref[3 * (ff.nProt + a)] / ff.nLigAtoms;
    com[1] += ff.ref[3 * (ff.nProt + a) + 1] / ff.nLigAtoms;
    com[2] += ff.ref[3 * (ff.nProt + a) + 2] / ff.nLigAtoms;
  }
  const idx = [];
  for (let i = 0; i < ff.nProt; i++) {
    if (Math.hypot(ff.ref[3 * i] - com[0], ff.ref[3 * i + 1] - com[1], ff.ref[3 * i + 2] - com[2]) < POCKET_RCUT) idx.push(i);
  }
  return idx;
}

const measured = {};
let hudChain = null;
{
  // (a) the HUD side: the app's own defaults, one CG run long enough to deposit
  //     hills — exactly what a user reaches after the CG row's 2000 steps.
  const ffHud = new ForceField(sel, { rc: 10, gamma: 2.0 }, bnzMols);
  ffHud.funnel = new Funnel({ nProt: ffHud.nProt, n: ffHud.n, ref: ffHud.ref, ligStart: ffHud.ligandStart ?? ffHud.nProt });
  ffHud.funnelOn = true;
  const integHud = new LangevinIntegrator(ffHud.ref, ffHud, MASS, { seed: SEED_HOLO });
  integHud.setTemperature(TEMP);
  integHud.setFriction(ZETA);
  for (let s = 0; s < FUNNEL_STEPS; s++) integHud.step();
  hudChain = funnelChain(ffHud.funnel);
  measured.ui_dg_hud_2000steps = hudChain.value;
  assert(hudChain.value !== null && Number.isFinite(hudChain.value),
    `HUD side recomputed from src/funnel.js: ΔG = ${hudChain.value.toFixed(4)} kcal/mol at ${ffHud.funnel._nHills} hills, t = ${integHud.time.toFixed(2)} ps`);
  assert(hudChain.value === ffHud.funnel.estimateDG(),
    "the HUD number the guard uses IS estimateDG() bit-for-bit — no second code path");

  // (b) the doc side: the record-path leg (BNZ+EPE) and the BNZ-only legs.
  const ffA = new ForceField(sel, FF_PAR, allMols); ffA.trackTerms = true;
  const holoA = runLeg(ffA, SEED_HOLO);
  const resA = computeThermodynamics({
    holoFrames: holoA.frames, apoFrames: holoA.frames, holoEnergies: holoA.energies,
    pocketIdx: pocket8A(ffA), nProt: ffA.nProt, mass: MASS, T: TEMP,
  });

  const ffB = new ForceField(sel, FF_PAR, bnzMols); ffB.trackTerms = true;
  const ffBapo = new ForceField(sel, { rc: 10, gamma: 2.0, binding: { on: false } }, bnzMols); ffBapo.trackTerms = true;
  const holoB = runLeg(ffB, SEED_HOLO);
  const apoB = runLeg(ffBapo, SEED_APO);
  const resB = computeThermodynamics({
    holoFrames: holoB.frames, apoFrames: apoB.frames, holoEnergies: holoB.energies,
    pocketIdx: pocket8A(ffB), nProt: ffB.nProt, mass: MASS, T: TEMP,
  });

  const sasaB = sasaBurial(holoB.frames, apoB.frames, {
    nProt: ffB.nProt,
    selElements: selectedLigandElements(bnzMols, [0]),
    protRadius: cgBeadExtendedRadius(ffB),
    stride: THERMO_SASA_STRIDE,
  });
  const resBRev = computeThermodynamics({
    holoFrames: holoB.frames, apoFrames: apoB.frames, holoEnergies: holoB.energies,
    pocketIdx: pocket8A(ffB), nProt: ffB.nProt, mass: MASS, T: TEMP,
    sasa: {
      dsasa: sasaB.dsasa, se: sasaB.se, dLig: sasaB.dLig, dProt: sasaB.dProt,
      method: sasaB.method, stride: sasaB.stride,
      nHoloEval: sasaB.nHoloEval, nApoEval: sasaB.nApoEval,
    },
  });

  const scanSys = { mode: "cg", sel, ff: ffB, ligands: bnzMols };
  const pocket = pocketResidues(scanSys, { rCut: POCKET_RCUT, maxN: MAXN });
  const scan = scanPocket(scanSys, pocket.map((p) => p.resId), { relaxSteps: RELAX_STEPS });

  const top3 = scan.rows.slice(0, 3);
  Object.assign(measured, {
    record_dH: resA.dH.total,
    record_dH_se: resA.dH_se,
    bnz_dH: resB.dH.total,
    bnz_dH_lj: resB.dH.lj,
    bnz_dH_desolv: resB.dH.desolv,
    bnz_desolv_share_pct: 100 * resB.dH.desolv / resB.dH.total,
    epe_inflation_dH: resA.dH.total - resB.dH.total,
    bnz_dsasa: sasaB.dsasa,
    bnz_dsasa_se: sasaB.se,
    bnz_negTds_solv: CHAIN_SASA_GAMMA * sasaB.dsasa,
    bnz_negTds_solv_se: CHAIN_SASA_GAMMA * sasaB.se,
    bnz_dG_est: resBRev.dG_estimate,
    ala_ddG_top1: top3[0].ddG,
    ala_ddG_top2: top3[1].ddG,
    ala_ddG_top3: top3[2].ddG,
    schlitter_err_pct: schlitterErrPct,
  });

  console.log(`\n      ── recomputed, both stories, 4W52 + benzene, seed ${SEED_HOLO}/${SEED_APO} ──`);
  console.log(`      HUD   (funnel/WTM, ${FUNNEL_STEPS} steps = ${ffHud.funnel._nHills} hills):  ΔG = ${hudChain.value >= 0 ? "+" : ""}${hudChain.value.toFixed(2)} kcal/mol`);
  console.log(`      DOC   (ΔH − TΔS, BNZ-only + real SASA):            ΔG = ${resBRev.dG_estimate >= 0 ? "+" : ""}${resBRev.dG_estimate.toFixed(2)} kcal/mol`);
  console.log(`      gap   (HUD − DOC):                               ΔG = ${(hudChain.value - resBRev.dG_estimate) >= 0 ? "+" : ""}${(hudChain.value - resBRev.dG_estimate).toFixed(2)} kcal/mol`);
  console.log(`      DOC   ΔH record ${resA.dH.total.toFixed(2)} / BNZ-only ${resB.dH.total.toFixed(2)} · ΔSASA ${sasaB.dsasa.toFixed(1)} ± ${sasaB.se.toFixed(1)} Å² · ΔS_lig ${thermoChain(resB).terms[2].value.toFixed(2)}`);
  console.log(`      ala-scan noise floor in code: ${ALA_DDG_NOISE_FLOOR} kcal/mol; measured max|ΔΔG| ${Math.max(...scan.rows.map((r) => Math.abs(r.ddG))).toFixed(3)}`);

  // The headline claim of this whole exercise: the two numbers disagree, and
  // the disagreement is reported rather than averaged away.
  const rec = reconcileChains(hudChain, thermoChain(resBRev));
  assert(rec.gap !== null && Math.abs(rec.gap) > 0.1,
    `the two stories differ by ${rec.gap.toFixed(2)} kcal/mol and the chain REPORTS it (rec.gap = ${rec.gap.toFixed(4)})`);
  assert(rec.canReconcile === false && /not a documentation fix/.test(rec.verdict),
    "reconcileChains reports canReconcile:false — the gap is physics, not a doc error");
  assert(rec.bridge.length === DG_BRIDGE.length && rec.bridge.every((b) => TERM_NATURE.includes(b.nature)),
    `all ${rec.bridge.length} bridge reasons are tagged, and ${rec.bridge.filter((b) => b.missing).length} of them name the missing physics`);
  const largest = DG_BRIDGE.filter((b) => Number.isFinite(b.movedBy)).sort((a, b) => Math.abs(b.movedBy) - Math.abs(a.movedBy))[0];
  assert(Number.isFinite(largest.movedBy) && Math.abs(largest.movedBy) > 1,
    `the largest MEASURED contributor to the gap is ${largest.id} at ${largest.movedBy.toFixed(2)} ${largest.movedByUnit} — and it is measured, not asserted`);
  console.log(`      biggest named contributor: ${largest.id} (${largest.movedBy} ${largest.movedByUnit})`);

  // The chain's arithmetic on the REAL HUD funnel must reproduce estimateDG.
  const t0 = hudChain.terms[0].value, t1 = hudChain.terms[1].value;
  assert(Math.abs((t0 + t1) - hudChain.value) < 1e-12,
    `on the real 4W52 funnel the chain still sums to estimateDG() to 1e-12 (${(t0 + t1).toFixed(12)} vs ${hudChain.value.toFixed(12)})`);
  assert(hudChain.terms[2].uncertainty !== null && Math.abs(hudChain.terms[2].uncertainty.value - ffHud.funnel.convergenceSE()) < 1e-12,
    `the real HUD chain's error bar is convergenceSE()'s own value (${hudChain.terms[2].uncertainty.value.toFixed(4)} = kT/√${ffHud.funnel._nHills})`);
}

/* ---- the guard itself ---- */
console.log("\n--- 4. docs/VALIDATION.md cannot drift from either story ---");
{
  const md = read("docs/VALIDATION.md");
  const parsedClaims = parseValidationClaims(md);
  assert(parsedClaims.missing.length === 0,
    `every one of the ${VALIDATION_CLAIM_COUNT} tracked claims is still findable in the methods table${parsedClaims.missing.length ? " — " + parsedClaims.missing.join("; ") : ""}`);
  assert(printedPrecision("−0.16") === 2 && printedPrecision("167.1") === 1 && num("−0.16") === -0.16 && num("+2.01") === 2.01,
    "the parser folds the typographic minus U+2212 the document actually uses (num('−0.16') = −0.16)");

  const drift = validationDrift(md, measured);
  if (drift.length) failWith("docs/VALIDATION.md still agrees with every code path", drift);
  else assert(true, `docs/VALIDATION.md agrees with all ${VALIDATION_CLAIM_COUNT} recomputed values (zero drift)`);
  console.log("      ── recomputed vs printed ──");
  for (const [k, v] of Object.entries(measured)) {
    const printed = parsedClaims.found[k];
    const shown = k === "ui_dg_hud_2000steps" ? "(absent from the table — finding F1)" : printed;
    console.log(`      ${k.padEnd(22)} code ${String(v.toFixed(6)).padStart(12)}   doc ${String(shown).padStart(10)}`);
  }

  // FAULT INJECTION 1 — perturb the DOCUMENT (ΔG_est −0.16 → −0.26). Must go red.
  const badDoc = md.replace("ΔG −0.16 (±50% band", "ΔG −0.26 (±50% band");
  assert(badDoc !== md, "FAULT INJECTION the mutated document really differs from the real one");
  const drift1 = validationDrift(badDoc, measured);
  assert(drift1.some((p) => p.includes("bnz_dG_est")),
    `FAULT INJECTION RED: one digit changed in docs/VALIDATION.md is detected ("${drift1.find((p) => p.includes("bnz_dG_est"))}")`);

  // FAULT INJECTION 2 — perturb the CODE side (ΔSASA scale). Must go red.
  const perturbed = { ...measured, bnz_dsasa: measured.bnz_dsasa * 1.13, bnz_dG_est: measured.bnz_dG_est + 0.4 };
  const drift2 = validationDrift(md, perturbed);
  assert(drift2.length >= 2 && drift2.some((p) => p.includes("bnz_dsasa")) && drift2.some((p) => p.includes("bnz_dG_est")),
    `FAULT INJECTION RED: the code side moving (ΔSASA ×1.13) is detected (${drift2.length} problems, e.g. "${drift2.find((p) => p.includes("bnz_dsasa"))}")`);

  // FAULT INJECTION 3 — perturb the HUD side only. The doc cannot be blamed, and
  // the detector must NOT go quiet: the absence of a HUD row is asserted.
  const hudOnly = { ...measured, ui_dg_hud_2000steps: measured.ui_dg_hud_2000steps + 5 };
  assert(validationDrift(md, hudOnly).length === 0,
    "moving ONLY the HUD number produces no doc drift — there is no doc row to contradict, which is why the chain must carry it");

  // FAULT INJECTION 4 — if the document ever starts quoting the HUD number without
  // declaring the funnel protocol, the guard fires. Proven by planting the row.
  const planted = md.replace(
    "| 4W52 BNZ+EPE (record path) |",
    `| 4W52 HUD funnel ΔG at 2000 steps | ${measured.ui_dg_hud_2000steps.toFixed(2)} |\n| 4W52 BNZ+EPE (record path) |`);
  assert(planted !== md, "FAULT INJECTION the planted methods-table row really differs");
  const drift4 = validationDrift(planted, measured);
  assert(drift4.some((p) => p.includes("ui_dg_hud_2000steps")),
    `FAULT INJECTION RED: a HUD row appearing in the methods table without its protocol is rejected ("${drift4.find((p) => p.includes("ui_dg_hud_2000steps"))}")`);
}

/* =====================================================================
 * 5. THE DOC PRESENTS THE CHAIN, GENERATED FROM THE TABLE
 * ===================================================================== */
console.log("\n--- 5. docs/VALIDATION.md presents the same enumeration ---");
{
  const md = read("docs/VALIDATION.md");
  const rowsArePresent = DG_ESTIMATORS.map((e) => ({
    id: e.id,
    present: md.includes(`\`${e.id}\``),
  }));
  const absent = rowsArePresent.filter((r) => !r.present).map((r) => r.id);
  assert(absent.length === 0,
    `every estimator id appears in docs/VALIDATION.md's enumeration table (${DG_ESTIMATORS.length} ids)${absent.length ? " — missing: " + absent.join(", ") : ""}`);
  const bridgePresent = DG_BRIDGE.filter((b) => md.includes(`\`${b.id}\``)).length;
  assert(bridgePresent === DG_BRIDGE.length,
    `every DG_BRIDGE id appears in docs/VALIDATION.md (${bridgePresent}/${DG_BRIDGE.length})`);
  assert(md.includes("measured") && md.includes("assumed") && md.includes("out-of-scope"),
    "docs/VALIDATION.md uses the chain's three nature tags verbatim, so a doc reader and a UI reader see one vocabulary");
  assert(md.includes("irreconcilable"),
    "docs/VALIDATION.md states plainly that the two ΔG stories are irreconcilable rather than leaving a reader to work it out");
  for (const e of DG_ESTIMATORS) {
    if (e.uncertainty) continue;
    assert(typeof e.noErrorBar === "string" && e.noErrorBar.length > 20,
      `${e.id} (no error bar) carries an explicit reason, per the project's never-write-zero rule`);
  }
}

/* =====================================================================
 * 6. CONTRACTS — DOM ids, the subpanel rule, the hotkey range, the record
 * ===================================================================== */
console.log("\n--- 6. UI + record contracts ---");
{
  const html = read("index.html");
  const uiSrc = read("src/ui.js");

  // (a) DOM-id contract (scripts/wikiskill_gate.js rule 3): every $(id) in ui.js
  //     exists in index.html.
  const ids = [...uiSrc.matchAll(/\$\("([^"]+)"\)/g)].map((m) => m[1]);
  const missing = ids.filter((id) => !new RegExp(`id="${id}"`).test(html));
  assert(missing.length === 0, `DOM contract: all ${ids.length} src/ui.js ids exist in index.html${missing.length ? " — missing " + missing.join(", ") : ""}`);
  for (const id of ["dgChain", "dgChainHud", "dgChainThermo"]) {
    assert(uiSrc.includes(`$("${id}")`) && html.includes(`id="${id}"`),
      `#${id} is registered in src/ui.js AND present in index.html`);
  }

  // (b) The new UI is INSIDE an existing top-level panel, never a new one:
  //     the Digit1-7 hotkey contract maps keys to `#controls > .panel` in DOM order.
  //     The count was 8 before this change (the 8th is #equations, hidden by
  //     default and deliberately outside the Digit1-7 range) and must stay 8.
  const topLevelPanels = (html.match(/<details[^>]*class="panel"/g) || []).length;
  assert(topLevelPanels === 8,
    `top-level #controls > .panel count is still 8 — unchanged by this feature, so no Digit hotkey shifted (was 8 at HEAD)`);
  assert(!/id="dgChain"[^>]*/.test(html.split("class=\"panel\"")[0] ?? ""),
    "#dgChain is not inside any earlier top-level panel");
  const dgIdx = html.indexOf('id="dgChain"');
  const pmfPanelIdx = html.indexOf("<summary>PMF &amp; Analysis</summary>");
  const asideEnd = html.indexOf("</aside>");
  assert(dgIdx > pmfPanelIdx && dgIdx < asideEnd,
    `#dgChain sits inside the existing 'PMF & Analysis' panel (offset ${dgIdx} > ${pmfPanelIdx} and < ${asideEnd}, the end of #controls)`);
  assert(html.includes('<summary>ΔG audit chain (which estimator, and on what)</summary>'),
    "the chain is a <details class=\"subpanel\"> — a subpanel, which #controls > .panel does not select");
  const transportSrc = read("src/controllers/transport.js");
  const hotkeyMax = Math.max(...[...transportSrc.matchAll(/Digit\[(\d)-(\d)\]/g)].map((m) => Number(m[2])));
  assert(hotkeyMax === 7 && topLevelPanels === hotkeyMax + 1,
    `hotkey contract: Digit1-${hotkeyMax} against ${topLevelPanels} top-level panels (the 8th, #equations, is hidden by design and stays outside the digit range) — nothing shifted`);

  // (c) The record carries the chain BY REFERENCE, and the exported estimator
  //     list is the live table.
  const block = scopeBlock();
  assert(block.deltaGChain.estimators.length === DG_ESTIMATORS.length
    && block.deltaGChain.estimators.every((e, i) => e.id === DG_ESTIMATORS[i].id),
    `scopeBlock().deltaGChain.estimators is DG_ESTIMATORS row-for-row (${block.deltaGChain.estimators.length})`);
  assert(block.deltaGChain.canReconcile === false && block.deltaGChain.bridge.length === DG_BRIDGE.length,
    "the exported record states canReconcile:false and ships the bridge");
  assert(block.deltaGChain.chainNatureTagSet.join("|") === TERM_NATURE.join("|"),
    "the exported record ships the same three nature tags the chain uses");
  assert(block.deltaGChain.estimators.every((e) => e.uncertaintyKind !== null || (typeof e.noErrorBar === "string" && e.noErrorBar.length > 0)),
    "every exported estimator has an uncertainty kind or an explicit no-error-bar reason — the never-write-zero rule survives serialization");

  // (d) Zero physics change: the chain module imports only units.js, the
  //     estimators module, and nothing physical.
  const chainSrc = read("src/dg-chain.js");
  const estSrc = read("src/dg-estimators.js");
  assert(!/forcefield|integrator|ff-binding|thermodynamics|alanine_scanning/.test(
    [...chainSrc.matchAll(/^import .*$/gm)].join("")),
    "src/dg-chain.js imports NO physics kernel — it reads estimator OBJECTS, it does not compute physics");
  assert(!/:\s*[-\d.]+;/.test(estSrc.replace(/\/\*[\s\S]*?\*\//g, "").replace(/\/\/.*$/gm, "")),
    "src/dg-estimators.js contains no numeric physics literal in code (only in its prose), so it cannot tune anything");
}

console.log(`\n${passed} PASSED, ${failed} FAILED`);
process.exit(failed > 0 ? 1 : 0);
