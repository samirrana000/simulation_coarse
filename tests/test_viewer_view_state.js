/**
 * test_viewer_view_state.js — Viewer system-parity + view-transform freshness.
 *
 * Replaces the two deleted ViewerGL stub suites:
 *   - test_rev2_issue3_viewergl_sync.js      (GL/Canvas2D system parity)
 *   - test_rev3_issue3_viewergl_view_sync.js (live view-transform delegation)
 *
 * Those two files asserted properties of the deleted `ViewerGL` wrapper
 * (prototype accessors, own-property shadowing, `gl === fallback` identity).
 * With the wrapper gone, every one of those assertions was either a tautology
 * or a statement about a class that no longer exists — zero coverage of the
 * real renderer.
 *
 * What they were actually standing in for IS tested here, against the real
 * `Viewer` (src/viewer.js), which is the only renderer the app ships:
 *   [A] setSystem derives ligandStart / colors / pocketCenter / radius on a
 *       system where `ligandStart !== nProt` (a hetero/cofactor atom sits in
 *       the gap). The distant hetero atom at (50,50,50) must be excluded from
 *       BOTH the pocket centre and the view radius.
 *   [B] The view transform (rotX/rotY/zoom/panX/panY) must survive a
 *       setSystem rebuild, while center/radius ARE recomputed from the new
 *       system — i.e. the camera the user set is never silently reset and the
 *       framing is never left stale. render() must not clobber it either.
 *
 * Run: node tests/test_viewer_view_state.js
 */

import { Viewer } from "../src/viewer.js";

function makeMockCtx() {
  return {
    canvas: null, fillStyle: "", strokeStyle: "", lineWidth: 1, font: "",
    fillRect() {}, fillText() {}, strokeText() {}, beginPath() {}, arc() {},
    fill() {}, stroke() {}, moveTo() {}, lineTo() {}, setLineDash() {},
    save() {}, restore() {}, clearRect() {}, closePath() {}, clip() {},
    strokeRect() {}, measureText() { return { width: 0 }; },
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
const close = (a, b, eps = 1e-9) => Math.abs(a - b) < eps;
const arrClose = (a, b, eps = 1e-9) => a.length === b.length && a.every((v, i) => close(v, b[i], eps));

console.log("=== test_viewer_view_state — system parity + view-transform freshness ===");

/* ------------------------------------------------------------------ *
 * [A] Hetero-gap system: ligandStart (4) !== nProt (3).
 *     idx 3 = HOH hetero at (50,50,50); idx 4 = BNZ ligand at (1,1,5).
 * ------------------------------------------------------------------ */
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
v.setSystem(sel, ff);

assert(v.n === 5 && v.nProt === 3, `system: n=5 nProt=3 (got ${v.n}/${v.nProt})`);
assert(v.ligandStart === 4, `hetero gap: ligandStart=4 !== nProt=3 (got ${v.ligandStart})`);
assert(v.heavy === true, "heavy flag carried from sel.heavy");
assert(JSON.stringify(v.segments) === JSON.stringify([[0, 3]]), "segments carried from sel");
assert(v.secStruct.length === 3 && v.secStruct.every((s) => s === "C"), `secStruct sized to nProt and all-coil (${JSON.stringify(v.secStruct)})`);
// holoSprings [i,j,d] -> strip the distance column (k % 3 === 2)
assert(JSON.stringify(Array.from(v.holoIdx)) === JSON.stringify([0, 1]), `holoIdx drops the distance column (${JSON.stringify(Array.from(v.holoIdx))})`);
// springs with d < 9.5 Å become contact pairs; both here are 5.0 and 6.0
assert(JSON.stringify(Array.from(v.contactIdx)) === JSON.stringify([0, 1, 1, 2]), `contactIdx from sub-9.5Å springs (${JSON.stringify(Array.from(v.contactIdx))})`);
assert(JSON.stringify(v.ligandBonds) === JSON.stringify([3, 4, 1]), "ligandBonds prefers covalentBonds over ligandBonds");

// class-before-element colouring across the hetero gap
const PROTEIN_C = [180, 180, 180];   // ELEMENT_COLOR.C
const PROTEIN_O = [235, 70, 70];     // ELEMENT_COLOR.O
const LIGAND_C = [255, 121, 98];     // LIGAND_COLOR.C
const LIGAND_O = [38, 208, 206];     // LIGAND_COLOR.O — must NOT reach idx 3
assert(JSON.stringify(v.colors[0]) === JSON.stringify(PROTEIN_C), "protein C keeps CPK grey");
assert(JSON.stringify(v.colors[3]) === JSON.stringify(PROTEIN_O), "hetero O (idx 3 < ligandStart) uses ELEMENT colour, not ligand palette");
assert(JSON.stringify(v.colors[3]) !== JSON.stringify(LIGAND_O), "hetero O does not pick up LIGAND_COLOR.O (class gate is ligandStart, not element)");
assert(JSON.stringify(v.colors[4]) === JSON.stringify(LIGAND_C), "true ligand C (idx 4 >= ligandStart) uses ligand palette");
assert(JSON.stringify(v.colors[3]) !== JSON.stringify(v.colors[4]), "hetero and ligand never share a colour class");

// The (50,50,50) hetero atom must not drag the view or the pocket with it.
assert(arrClose(v.pocketCenter, [1, 1, 5]), `pocketCenter = true-ligand COM only (${JSON.stringify(v.pocketCenter)})`);
assert(v.radius === 12, `radius ignores the distant hetero atom, stays at the 12Å floor (got ${v.radius})`);
assert(arrClose(v.center, [5 / 3, 5 / 3, 0]), `center = protein COM (${JSON.stringify(v.center.map((v2) => +v2.toFixed(3)))})`);

/* ------------------------------------------------------------------ *
 * [B] View transform survives a system rebuild; framing is refreshed.
 * ------------------------------------------------------------------ */
const sel2 = {
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
const mkFF = (x0, x1, x2) => ({
  n: 5, nProt: 3, ligandStart: 4,
  ref: new Float64Array([x0, 0, 0, x1, 0, 0, x2, 5, 0, 1, 1, 5, 2, 2, 6]),
  springs: [], holoSprings: [], covalentBonds: null, ligandBonds: null,
  ligandAtoms: [{ element: "C" }],
});

const cam = new Viewer(makeMockCanvas());
cam.setSystem(sel2, mkFF(0, 5, 0));
const firstCenter = cam.center.slice();
const firstRadius = cam.radius;

// User rotates / zooms / pans after the first load.
cam.rotX = 1.23; cam.rotY = -0.5; cam.zoom = 3.0; cam.panX = 5; cam.panY = -7;
const setCam = { rotX: cam.rotX, rotY: cam.rotY, zoom: cam.zoom, panX: cam.panX, panY: cam.panY };

// Force-field rebuild / a different system lands.
cam.setSystem(sel2, mkFF(100, 200, 100));
assert(!arrClose(cam.center, firstCenter), `center recomputed for the new system (${JSON.stringify(cam.center.map((x) => +x.toFixed(2)))})`);
assert(cam.radius !== firstRadius, `radius recomputed for the new system (${firstRadius} -> ${cam.radius})`);
assert(cam.rotX === setCam.rotX && cam.rotY === setCam.rotY, "rotX/rotY preserved across setSystem");
assert(cam.zoom === setCam.zoom, "zoom preserved across setSystem");
assert(cam.panX === setCam.panX && cam.panY === setCam.panY, "panX/panY preserved across setSystem");

// render() must not clobber the camera either.
cam.showStates = false; cam.showHBonds = false; cam.showContacts = false;
cam.drawRibbon = false; cam.drawSpheres = false;
cam.render(cam.ref);
assert(cam.rotX === setCam.rotX && cam.rotY === setCam.rotY, "rot preserved across render()");
assert(cam.zoom === setCam.zoom && cam.panX === setCam.panX && cam.panY === setCam.panY, "zoom/pan preserved across render()");

// Picking is invertible against the same live view transform.
const probe = cam.unproject(300, 225, 0);
assert(!!probe && probe.every(Number.isFinite), `unproject returns a finite world point under the live camera (${JSON.stringify(probe && probe.map((x) => +x.toFixed(2)))})`);

if (fails === 0) console.log("\n=== PASS test_viewer_view_state ===");
else { console.error(`\n=== FAIL test_viewer_view_state (${fails}) ===`); process.exit(1); }
