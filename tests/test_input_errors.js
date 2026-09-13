/**
 * tests/test_input_errors.js — FP2 input-robustness unit suite (FAST tier).
 *
 * Covers src/input_errors.js (central failure-class → actionable-copy
 * validator): frozen failure classes + thresholds, every class message
 * carries a next-step verb ("Click"), raw failures classify to the right
 * class, formatters keep technical detail secondary, safe parsers never
 * throw uncaught on garbage (structured { ok, data, error }), size/clash
 * guards fire only over threshold, and live 4W52/benzene files stay ok.
 *
 * Runnable: node tests/test_input_errors.js
 * Prints PASS/FAIL per assertion; exits 1 on any failure.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import {
  INPUT_ERROR_CODES,
  INPUT_ERROR_COPY,
  CLASH_RESIDUAL_THRESHOLD,
  MAX_HEAVY_ATOMS,
  MAX_CA_BEADS,
  MAX_SYSTEM_ATOMS,
  createInputError,
  classifyInputError,
  formatInputError,
  validatePdbText,
  guardParse,
  safeParseCa,
  safeParseMol2,
  safeParseLigands,
  checkSystemSize,
  checkPlacement,
} from "../src/input_errors.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

let passed = 0, failed = 0;
function assert(c, m) {
  if (c) { passed++; console.log(`  ✓ ${m}`); }
  else { failed++; console.error(`  ✗ FAIL: ${m}`); }
}

function findFile(name) {
  const candidates = [
    path.resolve(process.cwd(), name),
    path.resolve(__dirname, "..", name),
    path.resolve(__dirname, name),
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`${name} not found (tried ${candidates.join(", ")})`);
}

// ---- 1. frozen contract: codes + thresholds (3) ----
for (const c of ["EMPTY_PDB", "NO_ATOM", "NO_CA", "LIGAND_PARSE_FAIL", "NO_POCKET", "CLASH_HIGH", "SYSTEM_TOO_LARGE"]) {
  assert(INPUT_ERROR_CODES.includes(c), `failure class frozen: ${c}`);
}
assert(CLASH_RESIDUAL_THRESHOLD === 2.0, `clash residual threshold is 2.0 (got ${CLASH_RESIDUAL_THRESHOLD})`);
assert(MAX_HEAVY_ATOMS === 20000 && MAX_CA_BEADS === 5000 && MAX_SYSTEM_ATOMS === 20000,
  `state limits frozen: heavy ${MAX_HEAVY_ATOMS} / Cα ${MAX_CA_BEADS}`);

// ---- 2. every class message carries a next-step verb (9) ----
for (const c of INPUT_ERROR_CODES) {
  const e = createInputError(c, "synthetic technical detail");
  assert(/click/i.test(e.message) && /click/i.test(e.nextStep),
    `${c}: message + nextStep contain a next-step verb ("${e.nextStep.slice(0, 42)}…")`);
  assert(e instanceof Error && e.code === c, `${c}: mapped error is an Error with .code`);
}

// ---- 3. raw failures classify to the right class (8) ----
assert(classifyInputError(new Error("No Cα ATOM records found — is this a valid (protein) PDB file?")).code === "NO_CA",
  "raw parseCa throw → NO_CA");
assert(classifyInputError(new Error("No heavy atoms found in PDB file.")).code === "NO_ATOM",
  "raw parseHeavy throw → NO_ATOM");
assert(classifyInputError(new Error("lig.mol2: no usable molecules in MOL2 file"), { stage: "mol2" }).code === "LIGAND_PARSE_FAIL",
  "empty MOL2 (stage mol2) → LIGAND_PARSE_FAIL");
assert(classifyInputError(new Error("Selection must contain at least 3 Cα beads (need bonded + angle terms).")).code === "SELECTION_EMPTY",
  "raw selectSystem throw → SELECTION_EMPTY");
assert(classifyInputError("").code === "EMPTY_PDB", "empty string → EMPTY_PDB");
assert(classifyInputError(new Error("Could not download PDB XXXX from local presets, RCSB, or PDBe.")).code === "NO_ATOM",
  "fetch failure → NO_ATOM");
assert(classifyInputError(new Error("Placed pose still clashes residual 5.2")).code === "CLASH_HIGH",
  "clash residual prose → CLASH_HIGH");
assert(classifyInputError(new Error("System too large — exceeds the state limit")).code === "SYSTEM_TOO_LARGE",
  "oversize prose → SYSTEM_TOO_LARGE");

// ---- 4. classify totality: never throws, GENERIC fallback, idempotent (4) ----
for (const bad of [null, undefined, 42, {}, "some totally novel failure 12345"]) {
  let code = null;
  try { code = classifyInputError(bad).code; } catch (_) { code = "THREW"; }
  assert(code !== "THREW" && INPUT_ERROR_CODES.includes(code), `classify total on ${JSON.stringify(String(bad)).slice(0, 28)} → ${code}`);
}
const pre = createInputError("NO_CA", "orig detail");
assert(classifyInputError(pre) === pre, "classify passes already-mapped errors through untouched");

// ---- 5. format: actionable primary, detail secondary (3) ----
const fmt = formatInputError(new Error("No Cα ATOM records found — boom"));
assert(fmt.startsWith("⚠ ") && /click/i.test(fmt), "format starts ⚠ + actionable verb");
assert(fmt.includes("(detail:"), "format keeps raw detail as secondary (detail:) suffix");
assert(formatInputError(createInputError("NO_POCKET", INPUT_ERROR_COPY.NO_POCKET.message)) &&
  !formatInputError(createInputError("NO_POCKET", INPUT_ERROR_COPY.NO_POCKET.message)).includes("(detail:)"),
  "format omits redundant detail when technical === message");

// ---- 6. validatePdbText (4) ----
assert(validatePdbText("")?.code === "EMPTY_PDB", "validate: empty string → EMPTY_PDB");
assert(validatePdbText(null)?.code === "EMPTY_PDB", "validate: null → EMPTY_PDB");
assert(validatePdbText("HEADER    TEST\nTITLE     hello world\n")?.code === "NO_ATOM",
  "validate: header-only prose → NO_ATOM");
assert(validatePdbText("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n") === null,
  "validate: minimal CA line passes (null)");

// ---- 7. garbage-input matrix: safe parsers never throw uncaught (12) ----
const GARBAGE = [
  ["empty string", ""],
  ["null", null],
  ["undefined", undefined],
  ["number", 42],
  ["object", {}],
  ["prose garbage", "lorem ipsum dolor sit amet\nnot a pdb file\n"],
  ["header-only", "HEADER    TEST\nTITLE     hello\n"],
];
for (const [label, g] of GARBAGE) {
  let r;
  try { r = safeParseCa(g); } catch (_) { r = { threw: true }; }
  assert(r && r.threw !== true && r.ok === false && INPUT_ERROR_CODES.includes(r.error?.code),
    `safeParseCa never throws on ${label} → structured ${r.error?.code}`);
}
const HET_ONLY = "HETATM    1  C1  BNZ A 200      10.000  10.000  10.000  1.00 20.00           C  \n";
const hetCa = safeParseCa(HET_ONLY);
assert(hetCa.ok === false && hetCa.error?.code === "NO_CA", `HETATM-only text → NO_CA (got ${hetCa.error?.code})`);
assert(safeParseMol2("").error?.code === "LIGAND_PARSE_FAIL", "safeParseMol2: empty → LIGAND_PARSE_FAIL");
assert(safeParseMol2("lorem ipsum\nno sections here\n").error?.code === "LIGAND_PARSE_FAIL",
  "safeParseMol2: section-less prose → LIGAND_PARSE_FAIL");
assert(safeParseMol2("@<TRIPOS>BOND\n1 999 998 1\n").error?.code === "LIGAND_PARSE_FAIL",
  "safeParseMol2: dangling BOND refs → LIGAND_PARSE_FAIL");
const apoLig = safeParseLigands("HEADER    APO\nTITLE     no hetatm\nATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n");
assert(apoLig.ok === true && Array.isArray(apoLig.data) && apoLig.data.length === 0,
  "safeParseLigands: zero HETATM is valid apo (ok:true, []) — ligand errors fire at MOL2/place time");
assert(safeParseLigands(null).ok === false, "safeParseLigands: null → structured error, no throw");
assert(guardParse(() => { throw new Error("kaput"); }, "x", "NO_ATOM").ok === false, "guardParse catches throws → { ok:false }");

// ---- 8. size + placement guards (8) ----
assert(checkSystemSize({ nHeavy: 1308 }) === null, "size: 4W52 heavy 1308 within limit (null)");
assert(checkSystemSize({ nCa: 164 }) === null, "size: 4W52 Cα 164 within limit (null)");
assert(checkSystemSize({ nHeavy: 20000 }) === null, "size: exactly at heavy limit passes");
const big = checkSystemSize({ nHeavy: 20001 });
assert(big?.code === "SYSTEM_TOO_LARGE" && /click/i.test(big.message), "size: 20001 heavy → SYSTEM_TOO_LARGE + Click");
const bigCa = checkSystemSize({ nCa: 5001 });
assert(bigCa?.code === "SYSTEM_TOO_LARGE", "size: 5001 Cα → SYSTEM_TOO_LARGE");
assert(checkSystemSize({}) === null && checkSystemSize() === null, "size: missing counts never throw (null)");
assert(checkPlacement({ converged: true, residualClash: 1.05 }) === null, "placement: converged 1.05 → null (clash-free)");
const clash = checkPlacement({ converged: false, residualClash: 5.0 });
assert(clash?.code === "CLASH_HIGH" && /click/i.test(clash.message), "placement: unconverged 5.0 → CLASH_HIGH + Click");
assert(checkPlacement({ converged: true, residualClash: 3.5 })?.code === "CLASH_HIGH",
  "placement: over-threshold despite converged flag → CLASH_HIGH");
assert(checkPlacement(null)?.code === "CLASH_HIGH" && checkPlacement({ residualClash: NaN })?.code === "CLASH_HIGH",
  "placement: null/NaN residual → CLASH_HIGH (degenerate)");

// ---- 9. live-file sanity: bundled systems stay ok (3) ----
const pdbText = fs.readFileSync(findFile("4w52.pdb"), "utf-8");
const mol2Text = fs.readFileSync(findFile("benzene.mol2"), "utf-8");
const liveCa = safeParseCa(pdbText);
assert(liveCa.ok === true && liveCa.data.beads.length === 164, `live 4W52 safeParseCa ok (164 Cα, got ${liveCa.data?.beads?.length})`);
const liveMol2 = safeParseMol2(mol2Text);
assert(liveMol2.ok === true && liveMol2.data[0].atoms.length === 6 && liveMol2.data[0].bonds.length === 6,
  "live benzene.mol2 safeParseMol2 ok (6 atoms, 6 bonds)");
const liveLig = safeParseLigands(pdbText);
assert(liveLig.ok === true && liveLig.data.length >= 1, `live 4W52 safeParseLigands ok (${liveLig.data?.length} molecules)`);

console.log(`\nTEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
process.exit(failed ? 1 : 0);
