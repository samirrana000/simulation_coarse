/**
 * test_picking.js — H73 Picking invert motionGain
 * Tests Viewer with motionGain 5 picks within 0.5Å (fov unified 800).
 * Run: node tests/test_picking.js
 */

import { Viewer } from "../src/viewer.js";
import fs from "fs";

// mock canvas + ctx for Node (no DOM)
function makeMockCtx() {
  return {
    canvas: null,
    fillStyle: "",
    strokeStyle: "",
    lineWidth: 1,
    font: "",
    // no-op drawing methods used by Viewer.render
    fillRect() {},
    fillText() {},
    strokeText() {},
    beginPath() {},
    arc() {},
    fill() {},
    stroke() {},
    moveTo() {},
    lineTo() {},
    setLineDash() {},
    save() {},
    restore() {},
    clearRect() {},
    closePath() {},
    clip() {},
    strokeRect() {},
    measureText() { return { width: 0 }; },
  };
}

function makeMockCanvas(w=600, h=450) {
  const ctx = makeMockCtx();
  const canvas = {
    width: w,
    height: h,
    clientWidth: w,
    clientHeight: h,
    style: {},
    parentElement: null,
    getContext(type) {
      if (type === "2d") return ctx;
      return null;
    },
    getBoundingClientRect() {
      return { left: 0, top: 0, width: w, height: h, right: w, bottom: h };
    },
  };
  ctx.canvas = canvas;
  return canvas;
}

function assert(cond, msg) {
  if (!cond) {
    console.error(`  ✗ FAIL: ${msg}`);
    process.exitCode = 1;
    throw new Error(msg);
  } else {
    console.log(`  ✓ ${msg}`);
  }
}

function distance3(a, b) {
  const dx = a[0]-b[0], dy=a[1]-b[1], dz=a[2]-b[2];
  return Math.hypot(dx,dy,dz);
}

async function run() {
  console.log("=== test_picking — motionGain invert (H73) ===");

  // verify fov unified in source (grepable)
  const viewerSrc = fs.readFileSync("src/viewer.js", "utf-8");
  const fovHits = (viewerSrc.match(/fov = 800/g) || []).length;
  assert(fovHits >= 2, `fov unified fov=800 appears >=2 times (got ${fovHits}) — src/viewer.js fov unified`);
  assert(viewerSrc.includes("unified with render() fov=800") || viewerSrc.includes("unified with unproject() fov=800"),
    "fov unified comment present for picking invertibility");

  const canvas = makeMockCanvas(600, 450);
  const viewer = new Viewer(canvas);
  // Override to deterministic camera for test
  viewer.rotX = -0.4;
  viewer.rotY = 0.6;
  viewer.zoom = 1.0;
  viewer.panX = 0;
  viewer.panY = 0;

  // synthetic system: 4 beads tetrahedron-like
  const ref = new Float64Array([
    0, 0, 0,
    10, 0, 0,
    0, 10, 0,
    0, 0, 10,
  ]);
  const n = 4;
  const nProt = 4;
  const sel = {
    beads: [
      { element: "C", chain: "A", resName: "GLY", resSeq: 1 },
      { element: "C", chain: "A", resName: "ALA", resSeq: 2 },
      { element: "C", chain: "A", resName: "VAL", resSeq: 3 },
      { element: "C", chain: "A", resName: "LEU", resSeq: 4 },
    ],
    atoms: [
      { element: "C", chain: "A", resName: "GLY", resSeq: 1 },
      { element: "C", chain: "A", resName: "ALA", resSeq: 2 },
      { element: "C", chain: "A", resName: "VAL", resSeq: 3 },
      { element: "C", chain: "A", resName: "LEU", resSeq: 4 },
    ],
    segments: [[0, 4]],
    heavy: false,
  };
  const ff = {
    n,
    nProt,
    ref,
    springs: [],
    holoSprings: [],
    covalentBonds: null,
    ligandBonds: null,
  };

  viewer.setSystem(sel, ff);
  // ensure system built
  assert(viewer.n === n, `Viewer n=${n}`);
  assert(viewer.ref && viewer.ref.length === n*3, "Viewer ref set");
  // deterministic center/radius already computed; keep but also verify
  // force known radius for stable scale
  // viewer.center already computed as avg of ref; keep
  // viewer.radius already computed; keep

  // create displaced pos: each bead offset by small delta so motionGain has effect
  const pos = new Float64Array(ref);
  // bead 0: +1.2, -0.8, +0.5 ; bead1: +0.3 etc ; bead2 bead3 similar
  pos[0] += 1.2; pos[1] -= 0.8; pos[2] += 0.5;
  pos[3] += -0.7; pos[4] += 1.1; pos[5] += -0.4;
  pos[6] += 0.9; pos[7] += 0.2; pos[8] += -1.3;
  pos[9] += -1.0; pos[10] += -0.5; pos[11] += 0.9;

  // test with motionGain 5
  viewer.setMotionGain(5);
  assert(viewer.motionGain === 5, "motionGain 5 set");

  // render to populate _px,_py,_pz and _lastPos
  viewer.render(pos);
  assert(viewer._px && viewer._px.length === n, "_px populated");
  assert(viewer._lastPos && viewer._lastPos.length === n*3, "_lastPos stored");

  const rect = canvas.getBoundingClientRect();
  const W = canvas.width, H = canvas.height;
  const kx = W / Math.max(1, rect.width);
  const ky = H / Math.max(1, rect.height);

  let maxErr = 0;
  for (let i = 0; i < n; i++) {
    const px = viewer._px[i];
    const py = viewer._py[i];
    const pz = viewer._pz[i];
    // convert device px/py back to client coords
    const clientX = rect.left + px / kx;
    const clientY = rect.top + py / ky;
    const world = viewer.unproject(clientX, clientY, pz);
    assert(world !== null, `unproject bead ${i} non-null`);
    // also test screenToWorld alias
    const world2 = viewer.screenToWorld(clientX, clientY, { depth: pz });
    assert(world2 !== null, `screenToWorld bead ${i} non-null`);
    // distance to true pos (unamplified) should be <0.5Å with gain invert
    const truePos = [pos[3*i], pos[3*i+1], pos[3*i+2]];
    const err = distance3(world, truePos);
    const err2 = distance3(world2, truePos);
    console.log(`    bead ${i}: err=${err.toFixed(4)}Å (screenToWorld ${err2.toFixed(4)}Å) pos=[${truePos.map(v=>v.toFixed(2)).join(",")}] world=[${world.map(v=>v.toFixed(2)).join(",")}]`);
    maxErr = Math.max(maxErr, err, err2);
    assert(err < 0.5, `bead ${i} picks within 0.5Å with motionGain 5 (err=${err.toFixed(4)}Å)`);
    assert(err2 < 0.5, `bead ${i} screenToWorld within 0.5Å (err=${err2.toFixed(4)}Å)`);
  }
  console.log(`  max picking error with gain 5: ${maxErr.toFixed(4)}Å < 0.5Å ✓`);

  // also test with gain 1 baseline (should be even smaller)
  viewer.setMotionGain(1);
  viewer.render(pos);
  let maxErr1 = 0;
  for (let i = 0; i < n; i++) {
    const px = viewer._px[i], py = viewer._py[i], pz = viewer._pz[i];
    const clientX = rect.left + px / kx;
    const clientY = rect.top + py / ky;
    const world = viewer.unproject(clientX, clientY, pz);
    const truePos = [pos[3*i], pos[3*i+1], pos[3*i+2]];
    const err = distance3(world, truePos);
    maxErr1 = Math.max(maxErr1, err);
    assert(err < 0.5, `gain 1 bead ${i} within 0.5Å (err=${err.toFixed(4)}Å)`);
  }
  console.log(`  max picking error with gain 1: ${maxErr1.toFixed(4)}Å`);

  // also test that fov unified comment still ensures invertibility when zoom changed
  viewer.setMotionGain(5);
  viewer.zoom = 2.0;
  viewer.render(pos);
  let maxErrZoom = 0;
  for (let i = 0; i < n; i++) {
    const px = viewer._px[i], py = viewer._py[i], pz = viewer._pz[i];
    const clientX = rect.left + px / kx;
    const clientY = rect.top + py / ky;
    const world = viewer.unproject(clientX, clientY, pz);
    const truePos = [pos[3*i], pos[3*i+1], pos[3*i+2]];
    const err = distance3(world, truePos);
    maxErrZoom = Math.max(maxErrZoom, err);
    assert(err < 0.5, `zoom 2 gain5 bead ${i} within 0.5Å (err=${err.toFixed(4)}Å)`);
  }
  console.log(`  max picking error zoom2 gain5: ${maxErrZoom.toFixed(4)}Å`);

  console.log("\n=== PASS test_picking ===");
}

run().catch(e => {
  console.error("Test error:", e);
  process.exit(1);
});
