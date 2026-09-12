/**
 * test_bindviz.mjs — headless verification of bindviz renderers (R7).
 * Canvas-mock records draw calls; assertions on call counts + data math.
 * Run: node scripts/test_bindviz.mjs
 */
import { BindLog } from "../src/capture/bindlog.js";
import {
  renderInteractionTimeline,
  renderEnergyDecomposition,
  renderPmfFormation,
  renderParetoFrontier,
} from "../src/capture/bindviz.js";

let passed = 0, failed = 0;
const assert = (c, m) => { if (c) { passed++; console.log(`  ✓ ${m}`); } else { failed++; console.error(`  ✗ FAIL: ${m}`); } };

// ---- canvas mock: counts draw ops, captures text ----
function mockCanvas(w = 640, h = 220) {
  const calls = { fillRect: 0, strokeRect: 0, lineTo: 0, moveTo: 0, fillText: 0 };
  const texts = [];
  const ctx = {
    fillStyle: "", strokeStyle: "", font: "", textAlign: "left",
    fillRect: (...a) => calls.fillRect++,
    strokeRect: (...a) => calls.strokeRect++,
    beginPath: () => { calls._paths = (calls._paths || 0) + 1; },
    moveTo: () => calls.moveTo++,
    lineTo: () => calls.lineTo++,
    stroke: () => { calls._stroke = (calls._stroke || 0) + 1; },
    fill: () => { calls._fill = (calls._fill || 0) + 1; },
    closePath: () => {},
    fillText: (t) => { calls.fillText++; texts.push(t); },
  };
  return { width: w, height: h, getContext: () => ctx, _calls: calls, _texts: texts };
}

// ---- synthetic BindLog: known contacts + energies + hills ----
const bl = new BindLog();
// contacts: res 99 continuous 0-40ps (atom 0); res 84 forms 10 breaks 30 (atom 1);
// res 118 forms at 5, stays open (atom 2); res 71 brief 20-22 (atom 3)
bl.pushContact(0, true, 0, 99, 4.0);
bl.pushContact(10, true, 1, 84, 3.2);
bl.pushContact(5, true, 2, 118, 4.4);
bl.pushContact(20, true, 3, 71, 3.0);
bl.pushContact(22, false, 3, 71, 5.6);
bl.pushContact(30, false, 1, 84, 5.4);
bl.pushContact(40, false, 0, 99, 6.0);
// energies: 3 terms × 10 samples
for (let f = 0; f < 10; f++) {
  bl.pushEnergyComponents(f * 4, [-5 - 0.1 * f, -2 + 0.05 * f, -1]);
}
// hills: 5 hills around CV 3-5
bl.pushHill(1, 3.2, 0.10); bl.pushHill(2, 3.4, 0.09); bl.pushHill(3, 3.6, 0.08);
bl.pushHill(4, 4.8, 0.05); bl.pushHill(5, 3.1, 0.07);

// ---- (a) timeline ----
console.log("=== renderInteractionTimeline ===");
const c1 = mockCanvas();
renderInteractionTimeline(c1, bl, { tFrom: 0, tTo: 42 });
assert(c1._calls.fillText > 5, `timeline drew text/labels (${c1._calls.fillText} fillText)`);
assert(c1._calls.fillRect >= 4, `timeline drew ≥4 contact bars (${c1._calls.fillRect} fillRect)`);
assert(c1._texts.some(t => t.includes("99")), "row label includes residue 99");
assert(c1._texts.some(t => t.includes("ps")), "time axis labeled in ps");
// empty state
const c1e = mockCanvas(300, 140);
renderInteractionTimeline(c1e, new BindLog());
assert(c1e._texts.length >= 1 && c1e._texts[0].includes("no binding data"), "empty state text drawn");
// null canvas no-throw
let threw = false;
try { renderInteractionTimeline(null, bl); } catch { threw = true; }
assert(!threw, "null canvas does not throw");

// ---- (b) energy decomposition ----
console.log("=== renderEnergyDecomposition ===");
const c2 = mockCanvas();
renderEnergyDecomposition(c2, bl);
assert(c2._calls.lineTo >= 27, `3 term polylines × 10 pts (${c2._calls.lineTo} lineTo, ${c2._calls._paths || 0} paths)`);
assert(c2._texts.some(t => t.includes("ΔH")), "ΔH total annotated");
assert(c2._texts.some(t => t === "LJ") && c2._texts.some(t => t === "Coul"), "term legend names drawn");
// ΔH value check: means = (−5.45, −1.775, −1) → total −8.225
const dhText = c2._texts.find(t => t.startsWith("ΔH"));
assert(dhText && Math.abs(parseFloat(dhText.replace(/[^\d.\-+]/g, "")) - (-8.225)) < 0.01, `ΔH value ≈ −8.22 (${dhText})`);

// ---- (c) PMF formation ----
console.log("=== renderPmfFormation ===");
const c3 = mockCanvas();
renderPmfFormation(c3, bl, { binCount: 40 });
assert(c3._texts.some(t => t.includes("hills 5")), `hill count annotated (${c3._texts.find(t => t.includes("hills"))})`);
assert(c3._calls._fill >= 1, "PMF area filled");
// max bias > 0 → hills add up
const maxBias = c3._texts.find(t => t.includes("max bias"));
assert(maxBias && parseFloat(maxBias.replace(/[^\d.\-]/g, "")) > 0.1, `max bias positive (${maxBias})`);
// cursor: only first 2 hills
const c3b = mockCanvas();
renderPmfFormation(c3b, bl, { binCount: 40, cursorT: 2.5 });
assert(c3b._texts.some(t => t.includes("hills 2")), "cursorT limits to 2 hills");
// empty hills → no-data
const c3e = mockCanvas(300, 140);
renderPmfFormation(c3e, new BindLog());
assert(c3e._texts.length > 0 && c3e._texts[0].includes("no metadynamics"), "PMF empty state drawn");
// D1 (Loop-1 review): null bindlog → actionable empty state, never throw
const c3n = mockCanvas(300, 140);
let d1threw = false;
try { renderPmfFormation(c3n, null); } catch { d1threw = true; }
assert(!d1threw && c3n._texts.length > 0 && c3n._texts[0].includes("no binding data"), "PMF null-bindlog empty state drawn (D1 fix)");

// ---- (d) pareto frontier ----
console.log("=== renderParetoFrontier ===");
const tiers = [
  { name: "L0 CG", msPerStep: 0.203, top1: 1.0, pareto: true },
  { name: "L2 heavy", msPerStep: 82.5, top1: 0.9, pareto: true },
  { name: "L1 CG+ est", msPerStep: 0.23, top1: 0.95, pareto: false },
  { name: "L4 est", msPerStep: 45, top1: 0.9, pareto: false },
];
const c4 = mockCanvas();
renderParetoFrontier(c4, tiers);
assert(c4._calls.fillRect >= 4, `4 tier points drawn (${c4._calls.fillRect} fillRect)`);
assert(c4._calls._stroke >= 2, "axes + frontier line stroked");
assert(c4._texts.some(t => t.includes("L0 CG")) && c4._texts.some(t => t.includes("L2 heavy")), "tier labels drawn");
// empty tiers
const c4e = mockCanvas(300, 140);
renderParetoFrontier(c4e, []);
assert(c4e._texts.length > 0 && c4e._texts[0].includes("no tier"), "pareto empty state drawn");
// null canvas no-throw
threw = false;
try { renderParetoFrontier(null, tiers); } catch { threw = true; }
assert(!threw, "pareto null canvas does not throw");

// ---- BindLog round-trip through a viz path (integration with R6) ----
console.log("=== R6 integration: blob round-trip → timeline ===");
const blob = bl.toBinaryBlob();
const back = BindLog.fromBinaryBlob(blob);
const c5 = mockCanvas();
renderInteractionTimeline(c5, back, { tFrom: 0, tTo: 42 });
assert(c5._calls.fillRect >= 4, "timeline renders from blob-round-tripped BindLog");

console.log(`\n=== test_bindviz: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed ? 1 : 0);
