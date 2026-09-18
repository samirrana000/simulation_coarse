/**
 * test_rev1_issue4_live_terms.js — Revolution 1 / Issue 4 regression.
 *
 * Scope: src/main.js HUD/bindviz tick wiring + trackTerms live policy ONLY.
 * No FF energy-kernel change (kernels untouched).
 *
 * Problem: per-term contact energies (lj/coul/hb/desolv/pi/cpi/xb) were
 * computed only when trackTerms=true/BindLog enabled; a normal Run showed
 * the scalar U_bind but zeros in the energy channel.
 *
 * Fix verified here (minimal, backward compatible):
 *   - Always-on lightweight accumulators while a ligand is present
 *     (main.js applyLiveTrackTerms: bindLogWanted() || liveTermsWanted) —
 *     U/forces bit-identical ON vs OFF, per-terms live.
 *   - HUD tick shows U_bind + per-term split at 10 Hz via formatLiveTermsHUD
 *     (no BindLog checkbox); bindvizEnergy falls back to the live mirror
 *     (liveTermsBindLogView) at ≤1 Hz when the BindLog has no energy events.
 *   - BindLog event capture itself stays gated on bindLogWanted() (no extra
 *     memory, defaults backward compatible).
 *
 * Coverage (behavioral, no grep):
 *   [0]-[3] CG 4W52: OFF reproduction, ON bit-identity, split sums, L1 live.
 *   [4] main.js live helpers imported behaviorally (liveTermsWanted,
 *       applyLiveTrackTerms, readLiveTerms, formatLiveTermsHUD,
 *       updateLiveTermsMirror, liveTermsBindLogView): policy, HUD fragment,
 *       10 Hz mirror + ≤1 Hz synthetic BindLog view.
 *   [5] heavy weak-on FF (4W52 + benzene): trackTerms false vs true
 *       bit-identical U/forces/bindingU; pi+cpi+xb live when ON.
 *
 * Run: node tests/test_rev1_issue4_live_terms.js (fast, <2s; NOT wired into
 * tests/test_all.js FAST so the 352 gate is untouched).
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { parseCa, selectSystem, parseLigands, parseMol2 } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { parseHeavy, selectHeavy, appendHeavyLigands, HeavyForceField } from "../src/heavy.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

let passed = 0;
let failed = 0;
function assert(condition, message) {
  if (condition) {
    passed++;
    console.log(`  ✓ ${message}`);
  } else {
    failed++;
    console.error(`  ✗ FAIL: ${message}`);
  }
}

function findPdb() {
  const candidates = [
    path.resolve(process.cwd(), "4w52.pdb"),
    path.resolve(__dirname, "..", "4w52.pdb"),
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`4w52.pdb not found (tried ${candidates.join(", ")})`);
}

function findMol2() {
  const candidates = [
    path.resolve(process.cwd(), "benzene.mol2"),
    path.resolve(__dirname, "..", "benzene.mol2"),
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`benzene.mol2 not found (tried ${candidates.join(", ")})`);
}

function ensureDomStubs() {
  if (typeof globalThis.document === "undefined") {
    globalThis.document = {
      querySelectorAll: () => [],
      getElementById: () => null,
      activeElement: null,
    };
  } else {
    if (typeof globalThis.document.querySelectorAll !== "function") {
      globalThis.document.querySelectorAll = () => [];
    }
    if (typeof globalThis.document.getElementById !== "function") {
      globalThis.document.getElementById = () => null;
    }
  }
  if (typeof globalThis.window === "undefined") {
    globalThis.window = { addEventListener: () => {}, devicePixelRatio: 1 };
  } else if (typeof globalThis.window.addEventListener !== "function") {
    globalThis.window.addEventListener = () => {};
  }
}

async function main() {
  console.log("=== Rev1/Issue4: live per-term mirror without BindLog (4W52 CG) ===");
  const pdbText = fs.readFileSync(findPdb(), "utf-8");
  const parsed = parseCa(pdbText);
  const sel = selectSystem(parsed);
  const ligs = parseLigands(pdbText);
  assert(sel.beads.length === 164, `4W52 CG beads = 164 (got ${sel.beads.length})`);
  assert(ligs.length >= 1, `native ligands present (${ligs.length} mols)`);

  // -----------------------------------------------------------------
  // [0] Problem reproduction: trackTerms OFF → scalar live, terms zero
  // -----------------------------------------------------------------
  console.log("\n[0] problem reproduction (BindLog off ≡ trackTerms off)...");
  const ffOff = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligs);
  ffOff.trackTerms = false;
  const Uoff = ffOff.compute(ffOff.ref);
  const Boff = ffOff.bindingU;
  assert(Number.isFinite(Uoff) && Number.isFinite(Boff), `OFF: energy/bindingU finite`);
  assert(Math.abs(Boff) > 1e-9, `OFF: scalar U_bind nonzero (got ${Boff.toFixed(3)})`);
  assert(ffOff.bindLJU === 0 && ffOff.bindCoulU === 0 && ffOff.bindHBU === 0,
    `OFF: per-term accumulators stay zero (the reported bug)`);

  // -----------------------------------------------------------------
  // [1] Always-on accumulators are bit-identical (no kernel change)
  // -----------------------------------------------------------------
  console.log("\n[1] always-on accumulators bit-identical (U/forces)...");
  const ffOn = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligs);
  ffOn.trackTerms = true;
  const Un = ffOn.compute(ffOn.ref);
  const Bon = ffOn.bindingU;
  assert(Un === Uoff, `energy bit-identical OFF/ON (Δ=${Math.abs(Un - Uoff).toExponential(1)})`);
  assert(Bon === Boff, `bindingU bit-identical OFF/ON (Δ=${Math.abs(Bon - Boff).toExponential(1)})`);
  let maxD = 0;
  for (let i = 0; i < ffOff.forces.length; i++) maxD = Math.max(maxD, Math.abs(ffOff.forces[i] - ffOn.forces[i]));
  assert(maxD === 0, `forces bit-identical OFF/ON (maxΔ=${maxD.toExponential(1)})`);

  // -----------------------------------------------------------------
  // [2] Per-terms live when ON and split sums to bindingU
  // -----------------------------------------------------------------
  console.log("\n[2] per-term split live and sums to bindingU...");
  const { bindLJU: lj, bindCoulU: coul, bindHBU: hb, desolvU: desolv } = ffOn;
  assert([lj, coul, hb, desolv].every(Number.isFinite), `ON: all per-terms finite`);
  assert(Math.abs(lj) > 1e-9, `ON: LJ accumulator live (got ${lj.toFixed(3)})`);
  const split = lj + coul + hb + desolv;
  assert(Math.abs(split - Bon) < 1e-9,
    `ON: lj+coul+hb+desolv sums to bindingU (Δ=${Math.abs(split - Bon).toExponential(1)})`);
  assert(Math.abs((ffOn.bindU?.lj ?? 0) - lj) < 1e-12 && Math.abs((ffOn.bindU?.desolv ?? 0) - desolv) < 1e-12,
    `ON: bindU vector mirrors accumulators`);

  // -----------------------------------------------------------------
  // [3] L1 Coulomb path also live (salt-bridge branch, cf. Issue 3)
  // -----------------------------------------------------------------
  console.log("\n[3] L1 Coulomb accumulator live...");
  const ffL1 = new ForceField(sel, { rc: 10, gamma: 1.0, physicsLevel: "L1" }, ligs);
  ffL1.trackTerms = true;
  ffL1.compute(ffL1.ref);
  assert(Number.isFinite(ffL1.bindCoulU) && Math.abs(ffL1.bindCoulU) > 1e-9,
    `L1 bindCoulU nonzero (got ${ffL1.bindCoulU.toFixed(3)})`);

  // -----------------------------------------------------------------
  // [4] Behavioral live-terms helpers from src/main.js (no grep)
  // -----------------------------------------------------------------
  console.log("\n[4] behavioral live-terms helpers (src/main.js imports)...");
  ensureDomStubs();
  const {
    liveTermsWanted,
    applyLiveTrackTerms,
    readLiveTerms,
    formatLiveTermsHUD,
    updateLiveTermsMirror,
    liveTermsBindLogView,
  } = await import("../src/main.js");
  // Same-module state instance main.js mutates (ui.js?v=10 query pin).
  const { state: liveState } = await import("../src/ui.js?v=10");
  assert(typeof liveTermsWanted === "function" &&
    typeof applyLiveTrackTerms === "function" &&
    typeof readLiveTerms === "function" &&
    typeof formatLiveTermsHUD === "function" &&
    typeof updateLiveTermsMirror === "function" &&
    typeof liveTermsBindLogView === "function",
    `all six live-terms helpers imported from src/main.js`);

  // liveTermsWanted: ligand presence only, no checkbox.
  assert(liveTermsWanted(ffOn) === true,
    `liveTermsWanted true with ligand (nLigAtoms=${ffOn.nLigAtoms})`);
  const ffNoLig = new ForceField(sel, { rc: 10, gamma: 1.0 }, []);
  assert(liveTermsWanted(ffNoLig) === false, `liveTermsWanted false without ligand`);
  assert(liveTermsWanted(null) === false, `liveTermsWanted(null) false`);
  assert(liveTermsWanted({ nLigAtoms: 0 }) === false, `liveTermsWanted({nLigAtoms:0}) false`);

  // applyLiveTrackTerms: BindLog off headless ⇒ policy ≡ liveTermsWanted.
  const ffPolicy = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligs);
  ffPolicy.trackTerms = false;
  assert(applyLiveTrackTerms(ffPolicy) === true && ffPolicy.trackTerms === true,
    `applyLiveTrackTerms turns ON with ligand even when BindLog off`);
  const ffPolicyNo = new ForceField(sel, { rc: 10, gamma: 1.0 }, []);
  assert(applyLiveTrackTerms(ffPolicyNo) === false && ffPolicyNo.trackTerms === false,
    `applyLiveTrackTerms stays OFF without ligand`);
  assert(applyLiveTrackTerms(null) === false, `applyLiveTrackTerms(null) false`);

  // readLiveTerms: mirrors last per-term vector, safe on null.
  const snap = readLiveTerms(ffOn);
  assert(snap !== null &&
    Math.abs(snap.lj - ffOn.bindLJU) < 1e-12 &&
    Math.abs(snap.coul - ffOn.bindCoulU) < 1e-12 &&
    Math.abs(snap.hb - ffOn.bindHBU) < 1e-12 &&
    Math.abs(snap.desolv - ffOn.desolvU) < 1e-12 &&
    snap.bindingU === ffOn.bindingU,
    `readLiveTerms mirrors CG accumulators + bindingU (${snap.bindingU.toFixed(3)})`);
  assert(Math.abs((snap.lj + snap.coul + snap.hb + snap.desolv + snap.pi + snap.cpi + snap.xb) - snap.bindingU) < 1e-9,
    `readLiveTerms split sums to bindingU (CG, pi/cpi/xb zero)`);
  assert(readLiveTerms(null) === null, `readLiveTerms(null) null`);

  // formatLiveTermsHUD: U_bind + split fragment, "" when N/A.
  const hud = formatLiveTermsHUD(ffOn);
  assert(hud.includes("U_bind") && hud.includes("LJ") &&
    hud.includes("Coul") && hud.includes("HB") && hud.includes("desolv") &&
    hud.includes("kcal/mol"),
    `formatLiveTermsHUD renders U_bind + per-term split`);
  assert(hud.includes(ffOn.bindingU.toFixed(2)),
    `formatLiveTermsHUD embeds live U_bind (${ffOn.bindingU.toFixed(2)})`);
  assert(formatLiveTermsHUD(ffNoLig) === "" && formatLiveTermsHUD(null) === "",
    `formatLiveTermsHUD "" without ligand / null`);

  // updateLiveTermsMirror + liveTermsBindLogView: 10 Hz mirror, ≤1 Hz view.
  liveState._liveTerms = null;
  liveState._liveTermsHist = [];
  assert(liveTermsBindLogView() === null, `liveTermsBindLogView null with <2 samples`);
  assert(updateLiveTermsMirror(ffNoLig, 1) === null,
    `updateLiveTermsMirror null without ligand (no history write)`);
  assert((liveState._liveTermsHist?.length ?? 0) === 0,
    `mirror writes no history without ligand`);
  const s0 = updateLiveTermsMirror(ffOn, 0);
  assert(s0 !== null && s0.bindingU === ffOn.bindingU,
    `updateLiveTermsMirror returns snap at t=0 (no extra compute)`);
  assert(liveState._liveTermsHist.length === 1,
    `mirror history holds 1 sample after t=0`);
  const s1 = updateLiveTermsMirror(ffOn, 0.2);
  assert(s1 !== null, `second mirror sample at t=0.2 (≥0.09 ps cadence)`);
  assert(liveState._liveTermsHist.length === 2,
    `mirror history holds 2 samples (10 Hz cadence)`);
  const view = liveTermsBindLogView();
  assert(view !== null && view.nEvents === 14,
    `liveTermsBindLogView synthetic 2×7 energy events (got ${view?.nEvents})`);
  assert(view.evTime instanceof Float64Array && view.evType instanceof Uint8Array &&
    view.evA instanceof Int32Array && view.evX instanceof Float32Array,
    `liveTermsBindLogView typed arrays (BindLog-shaped)`);
  let allEnergy = true;
  for (let i = 0; i < view.nEvents; i++) if (view.evType[i] !== 0) { allEnergy = false; break; }
  assert(allEnergy, `synthetic view evType all energy (0)`);
  let chanOk = true;
  for (let i = 0; i < view.nEvents; i++) if (view.evA[i] !== (i % 7)) { chanOk = false; break; }
  assert(chanOk, `synthetic view channels cycle lj/coul/hb/desolv/pi/cpi/xb (0..6)`);
  assert(Math.abs(view.evTime[0] - 0) < 1e-12 && Math.abs(view.evTime[view.nEvents - 1] - 0.2) < 1e-12,
    `synthetic view preserves mirror times (0 → 0.2 ps)`);
  assert(Number.isFinite(view.evX[0]),
    `synthetic view carries live per-term values (evX finite)`);
  // Restore a clean tail for the heavy section (same-instance history).
  liveState._liveTerms = null;
  liveState._liveTermsHist = [];

  // -----------------------------------------------------------------
  // [5] Heavy weak-on FF: trackTerms OFF vs ON bit-identity, pi/cpi/xb live
  // -----------------------------------------------------------------
  console.log("\n[5] heavy weak-on trackTerms OFF vs ON (4W52 + benzene)...");
  const mol2Text = fs.readFileSync(findMol2(), "utf-8");
  const parsedHeavy = parseHeavy(pdbText);
  const benzMols = parseMol2(mol2Text);
  let selH = selectHeavy(parsedHeavy, {
    heteroSelection: { "A|200|BNZ": false, "A|201|EPE": false },
    includePdbLigands: true,
    hasExternalLigand: true,
  });
  selH = appendHeavyLigands(selH, benzMols, { gaff: true });
  assert(selH.atoms.length > 1000, `heavy 4W52+benzene atoms > 1000 (got ${selH.atoms.length})`);

  const ffHOff = new HeavyForceField({ atoms: selH.atoms }, { gamma: 2.0, weak: "on" }, []);
  ffHOff.trackTerms = false;
  const UHOff = ffHOff.compute(ffHOff.ref);
  const BHOff = ffHOff.bindingU;
  assert(Number.isFinite(UHOff) && Number.isFinite(BHOff), `heavy OFF: U/bindingU finite`);
  assert(ffHOff.bindLJU === 0 && ffHOff.bindCoulU === 0 &&
    ffHOff.bindHBU === 0 && ffHOff.desolvU === 0,
    `heavy OFF: lj/coul/hb/desolv accumulators stay zero`);
  assert(ffHOff.bindU.lj === 0 && ffHOff.bindU.coul === 0 && ffHOff.bindU.hb === 0 &&
    ffHOff.bindU.desolv === 0,
    `heavy OFF: bindU vector zero`);

  const ffHOn = new HeavyForceField({ atoms: selH.atoms }, { gamma: 2.0, weak: "on" }, []);
  ffHOn.trackTerms = true;
  const UHOn = ffHOn.compute(ffHOn.ref);
  const BHOn = ffHOn.bindingU;
  assert(UHOn === UHOff,
    `heavy weak-on U bit-identical OFF/ON (Δ=${Math.abs(UHOn - UHOff).toExponential(1)})`);
  assert(BHOn === BHOff,
    `heavy weak-on bindingU bit-identical OFF/ON (Δ=${Math.abs(BHOn - BHOff).toExponential(1)})`);
  let maxDH = 0;
  for (let i = 0; i < ffHOff.forces.length; i++) {
    maxDH = Math.max(maxDH, Math.abs(ffHOff.forces[i] - ffHOn.forces[i]));
  }
  assert(maxDH === 0, `heavy weak-on forces bit-identical OFF/ON (maxΔ=${maxDH.toExponential(1)})`);

  assert(Number.isFinite(ffHOn.bindLJU) && Math.abs(ffHOn.bindLJU) > 1e-9,
    `heavy ON: LJ accumulator live (got ${ffHOn.bindLJU.toFixed(3)})`);
  const hSplit = ffHOn.bindLJU + ffHOn.bindCoulU + ffHOn.bindHBU + ffHOn.desolvU;
  assert(Math.abs(hSplit - BHOn) < 1e-9,
    `heavy ON: lj+coul+hb+desolv sums to bindingU (Δ=${Math.abs(hSplit - BHOn).toExponential(1)})`);
  assert(Math.abs(ffHOn.bindU.lj - ffHOn.bindLJU) < 1e-12 &&
    Math.abs(ffHOn.bindU.desolv - ffHOn.desolvU) < 1e-12,
    `heavy ON: bindU vector mirrors accumulators`);
  assert(Math.abs(ffHOn.weakU - (ffHOn.piU + ffHOn.cpiU + ffHOn.xbU)) < 1e-12,
    `heavy ON: weakU = pi+cpi+xb (${ffHOn.weakU.toFixed(3)})`);
  assert(Math.abs(ffHOn.piU) + Math.abs(ffHOn.cpiU) + Math.abs(ffHOn.xbU) > 1e-9,
    `heavy ON: pi+cpi+xb live (π ${ffHOn.piU.toFixed(2)} · cπ ${ffHOn.cpiU.toFixed(2)} · XB ${ffHOn.xbU.toFixed(2)})`);
  assert(ffHOn.cpiU < 0, `heavy ON: cation-π negative contribution present (${ffHOn.cpiU.toFixed(3)})`);
  const hSnap = readLiveTerms(ffHOn);
  assert(Math.abs(hSnap.pi - ffHOn.piU) < 1e-12 &&
    Math.abs(hSnap.cpi - ffHOn.cpiU) < 1e-12 &&
    Math.abs(hSnap.xb - ffHOn.xbU) < 1e-12,
    `readLiveTerms mirrors heavy pi/cpi/xb`);
  const hHud = formatLiveTermsHUD(ffHOn);
  assert(hHud.includes("U_bind") && hHud.includes("π") && hHud.includes("cπ") && hHud.includes("XB"),
    `heavy HUD fragment shows π/cπ/XB when live`);
  assert(liveTermsWanted(ffHOn) === true, `liveTermsWanted true for heavy+ligand`);
  const ffHPolicy = new HeavyForceField({ atoms: selH.atoms }, { gamma: 2.0, weak: "on" }, []);
  ffHPolicy.trackTerms = false;
  assert(applyLiveTrackTerms(ffHPolicy) === true && ffHPolicy.trackTerms === true,
    `applyLiveTrackTerms turns heavy trackTerms ON with ligand`);
  const hS0 = updateLiveTermsMirror(ffHOn, 10.0);
  const hS1 = updateLiveTermsMirror(ffHOn, 10.2);
  assert(hS0 !== null && hS1 !== null && liveState._liveTermsHist.length === 2,
    `heavy mirror pushes 2 samples (10.0 → 10.2 ps)`);
  const hView = liveTermsBindLogView();
  assert(hView !== null && hView.nEvents === 14,
    `heavy synthetic view 2×7 events (got ${hView?.nEvents})`);
  liveState._liveTerms = null;
  liveState._liveTermsHist = [];

  // -----------------------------------------------------------------
  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log("=================================================");
  if (failed > 0) process.exit(1);
  console.log("PASS: test_rev1_issue4_live_terms.js — Rev1/Issue4 live mirror validated");
}

main().catch((e) => {
  console.error(e?.stack ?? String(e));
  process.exit(1);
});
