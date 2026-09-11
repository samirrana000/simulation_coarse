/**
 * test_ligand_colors.js — wiki P1: class-before-element coloring.
 * Heavy mode must assign ligand atoms (i >= nProt) a DISTINCT palette,
 * never the protein CPK map. Verifiable headlessly via Viewer.setSystem.
 * Run: node tests/test_ligand_colors.js
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

const PROTEIN_C = [180, 180, 180];
const LIGAND_C = [255, 121, 98];      // LIGAND_COLOR.C
const LIGAND_H_DEFAULT = [255, 159, 128];

console.log("=== test_ligand_colors — class-before-element (wiki P1) ===");

// Heavy-mode system: 3 protein atoms + 2 ligand atoms (all carbon)
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
  n: 5, nProt: 3,
  ref: new Float64Array([
    0, 0, 0,  5, 0, 0,  0, 5, 0,  0, 0, 5,  1, 1, 5,
  ]),
  springs: [], holoSprings: [], covalentBonds: null, ligandBonds: null,
  ligandAtoms: [
    { element: "C" }, { element: "C" },
  ],
};

const viewer = new Viewer(makeMockCanvas());
viewer.setSystem(sel, ff);

assert(viewer.n === 5 && viewer.nProt === 3, "system: 5 atoms, 3 protein");
assert(
  JSON.stringify(viewer.colors[0]) === JSON.stringify(PROTEIN_C),
  "protein C keeps CPK grey [180,180,180]"
);
assert(
  JSON.stringify(viewer.colors[3]) === JSON.stringify(LIGAND_C),
  "ligand C gets vivid palette [255,121,98], NOT protein CPK"
);
assert(
  JSON.stringify(viewer.colors[4]) === JSON.stringify(LIGAND_C),
  "second ligand C also vivid"
);
assert(
  JSON.stringify(viewer.colors[3]) !== JSON.stringify(viewer.colors[0]),
  "same element (C) differs across classes — visually distinguishable"
);

// CG mode: ligand (via ff.ligandAtoms) also gets LIGAND_COLOR
const selCg = { ...sel, heavy: false, beads: sel.beads.map((b) => ({ ...b })) };
const ffCg = { ...ff, n: 5, nProt: 3 };
const viewer2 = new Viewer(makeMockCanvas());
viewer2.setSystem(selCg, ffCg);
assert(
  JSON.stringify(viewer2.colors[3]) === JSON.stringify(LIGAND_C),
  "CG mode: ligand atom colored with LIGAND_COLOR palette"
);
assert(
  JSON.stringify(viewer2.colors[0]) !== JSON.stringify(LIGAND_C),
  "CG mode: protein bead uses chain palette, not ligand palette"
);

// Unknown element falls to LIGAND_COLOR_DEFAULT for ligands
const selX = { ...sel, beads: sel.beads.map((b, i) => i >= 3 ? { ...b, element: "Xx" } : { ...b }) };
const viewer3 = new Viewer(makeMockCanvas());
viewer3.setSystem(selX, ff);
assert(
  JSON.stringify(viewer3.colors[3]) === JSON.stringify(LIGAND_H_DEFAULT),
  "unknown ligand element falls back to LIGAND_COLOR_DEFAULT"
);

if (fails === 0) console.log("\n=== PASS test_ligand_colors ===");
else { console.error(`\n=== FAIL test_ligand_colors (${fails}) ===`); process.exit(1); }
