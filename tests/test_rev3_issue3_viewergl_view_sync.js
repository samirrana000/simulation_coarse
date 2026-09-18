/**
 * test_rev3_issue3_viewergl_view_sync.js — revolution 3 issue 3.
 * ViewerGL view transform (center/radius/rotX/rotY/zoom/panX/panY) must NOT be
 * value copies: fallback interaction (rotate/zoom/pan) and setSystem recompute
 * must be visible via gl immediately. Fix: get/set delegating to fallback.
 * Semantic mirrors from R2 (n/nProt/ligandStart/heavy/colors/segments/secStruct/
 * contactIdx/pocketCenter/ligandBonds/holoIdx/ref) stay as resynced copies.
 * Run: node tests/test_rev3_issue3_viewergl_view_sync.js
 */

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

console.log("=== test_rev3_issue3_viewergl_view_sync — live view-transform delegation ===");

const sel = {
  heavy: true,
  beads: [
    { element: "C", chain: "A", resName: "GLY", resSeq: 1 },
    { element: "N", chain: "A", resName: "GLY", resSeq: 1 },
    { element: "O", chain: "A", resName: "GLY", resSeq: 1 },
    { element: "C", chain: "L", resName: "BNZ", resSeq: 1 },
    { element: "C", chain: "L", resName: "BNZ", resSeq: 1 },
  ],
  atoms: [],
  segments: [[0, 3]],
};
const ff = {
  n: 5, nProt: 3, ligandStart: 4,
  ref: new Float64Array([0, 0, 0, 5, 0, 0, 0, 5, 0, 1, 1, 5, 2, 2, 6]),
  springs: [], holoSprings: [], covalentBonds: null, ligandBonds: null,
  ligandAtoms: [{ element: "C" }],
};

const g = new ViewerGL(makeMockCanvas());
g.setSystem(sel, ff);

const VIEW_FIELDS = ["center", "radius", "rotX", "rotY", "zoom", "panX", "panY"];

// [1] no own value-copy shadows the delegating accessor
for (const f of VIEW_FIELDS) {
  assert(!Object.hasOwn(g, f), `no own '${f}' copy (delegates to fallback)`);
  const d = Object.getOwnPropertyDescriptor(Object.getPrototypeOf(g), f);
  assert(!!(d && typeof d.get === "function" && typeof d.set === "function"),
    `prototype has get/set for '${f}'`);
}

// [2] fallback interaction mutation reflects in gl (the stale-copy bug)
g.fallback.rotX = 1.11; g.fallback.rotY = -2.22;
g.fallback.zoom = 2.5; g.fallback.panX = 7; g.fallback.panY = -9;
g.fallback.radius = 42; g.fallback.center = [11, 22, 33];
assert(g.rotX === 1.11 && g.rotY === -2.22, `rot follows fallback (${g.rotX},${g.rotY})`);
assert(g.zoom === 2.5 && g.panX === 7 && g.panY === -9, `zoom/pan follow fallback (${g.zoom},${g.panX},${g.panY})`);
assert(g.radius === 42, `radius follows fallback (${g.radius})`);
assert(JSON.stringify(g.center) === JSON.stringify([11, 22, 33]), `center follows fallback (${JSON.stringify(g.center)})`);
assert(g.center === g.fallback.center, "center is live reference (=== fallback.center)");

// [3] gl writes through to fallback (two-way delegation)
g.rotX = 0.5; g.rotY = 0.6; g.zoom = 1.7; g.panX = -3; g.panY = 4;
g.radius = 31; g.center = [1, 2, 3];
assert(g.fallback.rotX === 0.5 && g.fallback.rotY === 0.6, "gl.rot writes through");
assert(g.fallback.zoom === 1.7 && g.fallback.panX === -3 && g.fallback.panY === 4, "gl.zoom/pan write through");
assert(g.fallback.radius === 31, "gl.radius writes through");
assert(JSON.stringify(g.fallback.center) === JSON.stringify([1, 2, 3]), "gl.center writes through");

// [4] setSystem recompute visible via gl; rot/zoom/pan preserved live
g.rotX = 1.23; g.zoom = 3.0; g.panX = 5;
const ff2 = {
  ...ff,
  ref: new Float64Array([100, 0, 0, 105, 0, 0, 100, 5, 0, 1, 1, 5, 2, 2, 6]),
};
g.setSystem(sel, ff2);
assert(JSON.stringify(g.center) === JSON.stringify(g.fallback.center), `center resynced after setSystem (${JSON.stringify(g.center)})`);
assert(g.radius === g.fallback.radius, `radius resynced after setSystem (${g.radius})`);
assert(JSON.stringify(g.center) !== JSON.stringify([1, 2, 3]), "center actually recomputed (not stale)");
assert(g.rotX === 1.23 && g.zoom === 3.0 && g.panX === 5, "rot/zoom/pan preserved across setSystem");
assert(g.rotX === g.fallback.rotX && g.zoom === g.fallback.zoom, "rot/zoom still live after setSystem");

// [5] R2 semantic mirrors intact
assert(g.n === 5 && g.nProt === 3 && g.ligandStart === 4, "R2 mirrors: n/nProt/ligandStart");
assert(Array.isArray(g.colors) && g.colors.length === 5, "R2 mirrors: colors");
assert(g.ref === g.fallback.ref, "R2 mirrors: ref");

// [6] render keeps delegation intact
g.render(g.fallback.ref);
assert(g.rotX === g.fallback.rotX && g.zoom === g.fallback.zoom, "live after render");

if (fails === 0) console.log("\n=== PASS test_rev3_issue3_viewergl_view_sync ===");
else { console.error(`\n=== FAIL test_rev3_issue3_viewergl_view_sync (${fails}) ===`); process.exit(1); }
