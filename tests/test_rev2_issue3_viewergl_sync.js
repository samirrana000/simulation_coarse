/**
 * test_rev2_issue3_viewergl_sync.js — revolution 2 issue 3.
 * ViewerGL.setSystem must mirror ligandStart/colors/pocketCenter/contactIdx/
 * ligandBonds/holoIdx/heavy/segments/secStruct so the GL path API matches Viewer.
 * Run: node tests/test_rev2_issue3_viewergl_sync.js
 */

import { Viewer } from "../src/viewer.js";
import { ViewerGL } from "../src/viewer-gl.js";

function makeMockCtx() {
  return {
    canvas: null, fillStyle: "", strokeStyle: "", lineWidth: 1, font: "",
    fillRect() {}, fillText() {}, beginPath() {}, arc() {},
    fill() {}, stroke() {}, moveTo() {}, lineTo() {}, setLineDash() {},
    save() {}, restore() {}, clearRect() {},
  };
}

function makeMockCanvas(w = 600, h = 450) {
  const ctx = makeMockCtx();
  const canvas = {
    width: w, height: h, clientWidth: w, clientHeight: h, style: {},
    parentElement: null,
    getContext(type) { return type === "2d" ? ctx : null; },
    getBoundingClientRect() { return { left: 0, top: 0, width: w, height: h, right: w, bottom: h }; },
  };
  ctx.canvas = canvas;
  return canvas;
}

let fails = 0;
function assert(cond, msg) {
  if (!cond) { console.error(`  ✗ FAIL: ${msg}`); fails++; process.exitCode = 1; }
  else console.log(`  ✓ ${msg}`);
}

console.log("=== test_rev2_issue3_viewergl_sync — GL/Canvas2D system parity ===");

// Hetero system: 3 protein + 1 hetero (idx 3) + 1 true ligand (idx 4).
// ligandStart (4) !== nProt (3) exercises the hetero-exclusion path.
const sel = {
  heavy: true,
  beads: [
    { element: "C", chain: "A", resName: "GLY", resSeq: 1 },
    { element: "N", chain: "A", resName: "GLY", resSeq: 1 },
    { element: "O", chain: "A", resName: "GLY", resSeq: 1 },
    { element: "O", chain: "A", resName: "HOH", resSeq: 2 },
    { element: "C", chain: "L", resName: "BNZ", resSeq: 1 },
  ],
  atoms: [],
  segments: [[0, 3]],
};
const ff = {
  n: 5, nProt: 3, ligandStart: 4,
  ref: new Float64Array([0, 0, 0, 5, 0, 0, 0, 5, 0, 50, 50, 50, 1, 1, 5]),
  springs: [0, 1, 5.0, 1, 2, 6.0],
  holoSprings: [0, 1, 5.0],
  covalentBonds: [3, 4, 1],
  ligandBonds: null,
  ligandAtoms: [{ element: "C" }],
};

const v = new Viewer(makeMockCanvas());
const g = new ViewerGL(makeMockCanvas());
v.setSystem(sel, ff);
g.setSystem(sel, ff);

assert(g.ligandStart === v.ligandStart, `ligandStart mirrors (${g.ligandStart}===${v.ligandStart})`);
assert(g.ligandStart === 4, "hetero case: ligandStart (4) !== nProt (3)");
assert(JSON.stringify(g.colors) === JSON.stringify(v.colors), "colors equal");
assert(JSON.stringify(g.pocketCenter) === JSON.stringify(v.pocketCenter), `pocketCenter equal (${JSON.stringify(g.pocketCenter)})`);
assert(JSON.stringify(Array.from(g.contactIdx || [])) === JSON.stringify(Array.from(v.contactIdx || [])), "contactIdx equal");
assert(JSON.stringify(Array.from(g.holoIdx || [])) === JSON.stringify(Array.from(v.holoIdx || [])), "holoIdx equal");
assert(JSON.stringify(g.ligandBonds) === JSON.stringify(v.ligandBonds), "ligandBonds equal");
assert(g.heavy === v.heavy, "heavy equal");
assert(JSON.stringify(g.segments) === JSON.stringify(v.segments), "segments equal");
assert(JSON.stringify(g.secStruct) === JSON.stringify(v.secStruct), "secStruct equal");
assert(g.n === v.n && g.nProt === v.nProt, "n/nProt still mirror");
assert(JSON.stringify(g.center) === JSON.stringify(v.center) && g.radius === v.radius, "center/radius still mirror");

if (fails === 0) console.log("\n=== PASS test_rev2_issue3_viewergl_sync ===");
else { console.error(`\n=== FAIL test_rev2_issue3_viewergl_sync (${fails}) ===`); process.exit(1); }
