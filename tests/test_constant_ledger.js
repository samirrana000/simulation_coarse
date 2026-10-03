/**
 * test_constant_ledger.js — a SCREAMING_CASE constant has exactly one home.
 *
 * THE PRECEDENT (this file exists because of it)
 * ---------------------------------------------
 * The loop found a real defect: `KB_KCAL` was declared twice with different
 * values —
 *     src/units.js     0.001987204
 *     src/ff-params.js 0.0019872041
 * — a 5.03e-8 relative split. The reported instantaneous temperature depended
 * on which module an engine happened to import from (CG read units.js, heavy
 * read ff-params.js). It survived seven rounds of binding-physics review.
 * tests/test_unit_contract.js locks that ONE pair.
 *
 * A follow-up scan then found the same defect still live elsewhere, under a
 * different name, in a file the first test never looked at:
 *     src/analysis/thermodynamics.js  `const KB = 0.0019872041`
 * and a family of Coulomb-constant splits (332.0 / 332.0637 / 332.06371)
 * across six files, plus a bare `332.0` inside a WGSL template string that no
 * import-based test can see. A per-name test cannot catch that class. This one
 * scans EVERY module, so it does.
 *
 * WHAT IT ASSERTS
 *   1. Single-site: a SCREAMING_CASE name may carry a numeric literal in at
 *      most one file, unless explicitly allowlisted with a reason. Same-value
 *      duplicates are a maintenance hazard (the KB_KCAL pair was one edit away
 *      from becoming a DIFFERENT-value bug).
 *   2. Divergence: two sites for one name with DIFFERENT literals fails
 *      unconditionally — that is a live physics bug, not a style issue.
 *   3. Unit-shadowing: no module may re-derive a units.js contract constant
 *      under a different name (this is what catches `KB = ...`,
 *      `ELC = 332.0637`, `K_ELEC = 332.0`).
 *   4. Runtime: the contract constants are importable and node-safe, and the
 *      re-export sites are bit-identical to units.js.
 *   5. TABLE ledger (added by goal M7): rules 1–4 only see a constant whose
 *      initializer is a bare NUMBER, so every parameter TABLE — an object
 *      literal keyed by chemical element — walked straight through. That is
 *      exactly where the repo's third instance of this defect class lived:
 *          src/ff-params.js:112  LIG_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, … }
 *          src/heavy.js:82       HEAVY_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, … }
 *      Two declarations of one quantity in two modules, one edit from a live
 *      sigma/eps/q split between the CG and heavy paths. Rule 5 makes that
 *      shape un-declarable outside the canonical parameter module unless it
 *      carries a written justification (which then has to still be in use).
 *
 * THE ALLOWLIST IS THE DELIVERABLE, NOT A PERMISSION SLIP
 * Every entry carries a `why` and every entry is asserted to be STILL IN USE,
 * so an entry that stops matching fails the suite. A stale exemption is a
 * silent hole; this makes it impossible to leave one behind.
 *
 * Zero dependencies. Node-importable; no DOM globals.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import * as units from "../src/units.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
let passed = 0, failed = 0;

function assert(cond, msg) {
  if (cond) { passed++; console.log(`  ✓ ${msg}`); }
  else { failed++; console.error(`  ✗ FAIL: ${msg}`); }
}

function srcFiles(dir = path.join(ROOT, "src"), out = []) {
  for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) srcFiles(p, out);
    else if (e.name.endsWith(".js")) out.push(p);
  }
  return out;
}

/**
 * Blank out comments and every string/template literal, preserving offsets
 * and newlines. A constant mentioned in prose, or embedded in a WGSL source
 * string, is not a JavaScript declaration — but a bare literal inside a
 * template IS a real risk, so callers that care about shader text use
 * `shaderLiterals` below rather than this function.
 */
function blankCommentsAndStrings(src) {
  const out = src.split("");
  let i = 0, mode = null;
  while (i < src.length) {
    const c = src[i], d = src[i + 1];
    if (mode === null) {
      if (c === "/" && d === "/") { mode = "line"; out[i] = " "; out[i + 1] = " "; i += 2; continue; }
      if (c === "/" && d === "*") { mode = "block"; out[i] = " "; out[i + 1] = " "; i += 2; continue; }
      if (c === "'") { mode = "sq"; out[i] = " "; i++; continue; }
      if (c === '"') { mode = "dq"; out[i] = " "; i++; continue; }
      if (c === "`") { mode = "tpl"; out[i] = " "; i++; continue; }
      i++; continue;
    }
    if (mode === "line") { if (c === "\n") mode = null; else out[i] = " "; i++; continue; }
    if (mode === "block") {
      if (c === "*" && d === "/") { out[i] = " "; out[i + 1] = " "; i += 2; mode = null; continue; }
      if (c !== "\n") out[i] = " ";
      i++; continue;
    }
    if (c === "\\") { out[i] = " "; if (src[i + 1] !== "\n") out[i + 1] = " "; i += 2; continue; }
    if ((mode === "sq" && c === "'") || (mode === "dq" && c === '"') || (mode === "tpl" && c === "`")) {
      mode = null; out[i] = " "; i++; continue;
    }
    if (c !== "\n") out[i] = " ";
    i++;
  }
  return out.join("");
}

/**
 * Blank COMMENTS only, preserving string/template content. The shader check
 * needs this: it must see literals inside a WGSL template string, but must
 * not be fooled by a doc comment that merely NAMES the old wrong value
 * ("previously hardcoded `332.0`") — which is exactly the kind of note that
 * should be allowed to exist.
 */
function blankCommentsOnly(src) {
  const out = src.split("");
  let i = 0, mode = null;
  while (i < src.length) {
    const c = src[i], d = src[i + 1];
    if (mode === null) {
      if (c === "/" && d === "/") { mode = "line"; out[i] = " "; out[i + 1] = " "; i += 2; continue; }
      if (c === "/" && d === "*") { mode = "block"; out[i] = " "; out[i + 1] = " "; i += 2; continue; }
      if (c === "'") { mode = "sq"; i++; continue; }
      if (c === '"') { mode = "dq"; i++; continue; }
      if (c === "`") { mode = "tpl"; i++; continue; }
      i++; continue;
    }
    if (mode === "line") { if (c === "\n") mode = null; else out[i] = " "; i++; continue; }
    if (mode === "block") {
      if (c === "*" && d === "/") { out[i] = " "; out[i + 1] = " "; i += 2; mode = null; continue; }
      if (c !== "\n") out[i] = " ";
      i++; continue;
    }
    // Inside a string: only the delimiters and escapes are structural.
    if (c === "\\") { i += 2; continue; }
    if ((mode === "sq" && c === "'") || (mode === "dq" && c === '"') || (mode === "tpl" && c === "`")) mode = null;
    i++;
  }
  return out.join("");
}

const NAME = "[A-Z][A-Z0-9]*(?:_[A-Z0-9]+)*";
const NUM = "[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?";
// One declaration statement, possibly a multi-declarator list:
//   const HS = x, RMIN = 2.6, KF = 8.0;
const STMT = new RegExp(
  `\\b(?:export\\s+)?(?:const|let|var)\\s+((?:${NAME}\\s*=[^;,]*)(?:\\s*,\\s*${NAME}\\s*=[^;,]*)*)\\s*[;\\n]`,
  "g");
const ITEM = new RegExp(`(${NAME})\\s*=\\s*(${NUM})(?![\\w.])`, "g");

/** All SCREAMING_CASE numeric declarations under src/, with file:line. */
function scanDeclarations() {
  const rows = [];
  for (const f of srcFiles()) {
    const raw = fs.readFileSync(f, "utf-8");
    const code = blankCommentsAndStrings(raw);
    STMT.lastIndex = 0;
    let s;
    while ((s = STMT.exec(code)) !== null) {
      ITEM.lastIndex = 0;
      let it;
      while ((it = ITEM.exec(s[1])) !== null) {
        if (!/[A-Z]{2}/.test(it[1])) continue;      // drops `let x`, `U`, `M`
        const abs = s.index + it.index;
        rows.push({
          name: it[1], literal: it[2],
          file: path.relative(ROOT, f),
          line: code.slice(0, abs).split("\n").length,
        });
      }
    }
  }
  return rows;
}

// ─────────────────────────────────────────────────────────────────────
// ALLOWLIST — every entry must justify itself AND be still in use.
// ─────────────────────────────────────────────────────────────────────

/**
 * Names allowed to carry a literal in more than one file. Same-value
 * duplicates only; a divergent pair is never allowlistable (see DIVERGENCE).
 */
const ALLOW_MULTI_SITE = [
  {
    name: "R_CUT",
    why:
      "src/heavy/params.js R_CUT=8.5 terminates a smooth C2 switch over " +
      "[6.5, 8.5] plus GB; src/force-worker.js R_CUT=8.5 terminates an " +
      "exponential-screened exp(-r/8)/r^2 kernel with no switch. Different " +
      "potentials that happen to share a truncation radius — binding them would " +
      "assert an identity that does not hold, and the worker value is a " +
      "per-backend tuning constant already parameterised as RCUT_DEFAULT in " +
      "compute/webgpu_backend.js. Both sites carry this justification. (The " +
      "heavy site moved here in the 2026-10 heavy split; this string said " +
      "src/heavy.js for the whole life of the exemption, which after that split " +
      "named a facade that no longer declares anything. A justification that " +
      "points at the wrong file is worse than no justification.)",
  },
];

// ─────────────────────────────────────────────────────────────────────
// RULE 5 — THE TABLE LEDGER (goal M7)
//
// Rules 1–4 match `NAME = <number>`. A force-field PARAMETER TABLE is
// `NAME = { SIGMA: …, EPS: … }`, which none of them can see — and a table is
// where a chemical element's LJ parameters live, so the whole class of
// "two modules declare one physical quantity" defects that this file exists
// for was completely unguarded until M7.
//
// THE SHAPE BEING HUNTED: a SCREAMING_CASE declaration whose initializer is an
// object literal carrying ≥2 of the five canonical element-parameter fields
// (sigma / eps / q / hb / dG). That catches all three cases that matter:
//
//   • the full table          ELEMENT_LJ, METAL_ELEMENT, RES_CLASS
//   • a per-element sub-table  NONBONDED_TABLE (AMBER parm99 — a DIFFERENT
//                              physical model, hence allowlisted below)
//   • a PROJECTION of one     `HEAVY_ELEMENT_DEFAULT = { sigma, eps, q }`
//                              — a single default object carrying element
//                              parameter fields, which is precisely what
//                              src/heavy.js:82 was, and precisely what must
//                              not come back.
//
// `hb`/`dG` are included because they are part of the same row; a table with
// only them is still a per-element parameter table.
// ─────────────────────────────────────────────────────────────────────

/** The canonical home. Element parameters are declared HERE or not at all. */
const PARAM_HOME = path.join("src", "physics", "params.js");

/** The field set that marks an object literal as an element-parameter table. */
const ELEMENT_PARAM_FIELDS = ["sigma", "eps", "q", "hb", "dG"];

/**
 * Element-parameter tables that are legitimately a DIFFERENT physical model
 * and therefore live outside the canonical home. Each needs a substantive
 * reason and must still be present (see the rot check below), so an exemption
 * cannot outlive its justification.
 */
const ALLOW_ELEMENT_TABLE = [
  {
    name: "NONBONDED_TABLE",
    file: path.join("src", "physics", "forcefield", "amber14sb.js"),
    why:
      "AMBER parm99 non-bonded LJ per element, sigma converted from Rmin/2 — a " +
      "PUBLISHED FORCE FIELD, not the coarse bead-scale table in params.js. " +
      "The two genuinely differ and are used for different things: params.js " +
      "ELEMENT_LJ drives both engines' non-bonded LJ grid, this one is the " +
      "AMBER fallback consulted by getNonbondedParams(). Worst-case epsilon gap " +
      "between them is 45% (I: 0.4 vs 0.22), measured 2026-10-03 across every " +
      "element the two tables share; so silently merging or syncing them would " +
      "be a physics change. tests/test_element_params.js " +
      "measures the gap on every run; docs/CHARGES.md scopes the heavy model.",
  },
  {
    name: "NONBONDED_DEFAULT",
    file: path.join("src", "physics", "forcefield", "amber14sb.js"),
    why:
      "The AMBER table's own unknown-element fallback, which must stay with " +
      "that table: getNonbondedParams() returns it for any element parm99 does " +
      "not cover. Its sigma/eps happen to equal params.ELEMENT_LJ_DEFAULT " +
      "(3.4 / 0.12) but they are a different model's fallback, and asserting " +
      "they stay equal would assert an identity that is not a fact about the " +
      "two models — it is a coincidence of two independent choices.",
  },
];

/**
 * All SCREAMING_CASE declarations in src/ whose initializer is an object
 * literal with ≥2 element-parameter fields. Returns [{name, file, line}].
 */
function scanElementParamTables() {
  const rows = [];
  // `(?:export\s+)?(?:const|let|var)\s+NAME\s*=\s*\{` … to the matching `}`.
  // The body scan is brace-counted rather than regexed so a nested object does
  // not truncate the match, and a field is only counted when it is a real key
  // (`sigma:`) rather than a mention in a nested value.
  const DECL = /\b(?:export\s+)?(?:const|let|var)\s+([A-Z][A-Z0-9]*(?:_[A-Z0-9]+)*)\s*=\s*\{/g;
  for (const f of srcFiles()) {
    const raw = fs.readFileSync(f, "utf-8");
    const code = blankCommentsAndStrings(raw);
    DECL.lastIndex = 0;
    let m;
    while ((m = DECL.exec(code)) !== null) {
      const name = m[1];
      if (!/[A-Z]{2}/.test(name)) continue;      // drops `let x`, `U`, `M`
      // Brace-count the literal body from the `{` DECL just matched.
      let depth = 0, i = m.index + m[0].length - 1;
      const start = i;
      for (; i < code.length; i++) {
        if (code[i] === "{") depth++;
        else if (code[i] === "}") { depth--; if (depth === 0) break; }
      }
      const body = code.slice(start, i + 1);
      const hits = ELEMENT_PARAM_FIELDS.filter((fld) =>
        new RegExp(`[{,\\s]${fld}\\s*:`).test(body));
      if (hits.length < 2) continue;
      rows.push({
        name, fields: hits,
        file: path.relative(ROOT, f),
        line: code.slice(0, start).split("\n").length,
      });
    }
  }
  return rows;
}

/**
 * Declared same-value ALIASES inside the contract file itself, where a second
 * name for one quantity is deliberate API surface rather than drift.
 */
const ALLOW_CONTRACT_ALIAS = [
  {
    name: "KCAL_TO_DA_A2_PS2",
    aliases: "KCONV",
    why:
      "units.js exports BOTH names for the same 418.4 Da·A^2/ps^2 factor: " +
      "KCAL_TO_DA_A2_PS2 is the descriptive name, KCONV is the short name " +
      "kept for backward compatibility with integrator.js / funnel.js / the " +
      "test suite. Same value, one file, documented in the file header and " +
      "asserted equal by tests/test_unit_contract.js. Collapsing them would " +
      "be a public-API rename, which is out of scope for a de-duplication.",
  },
];

/**
 * Declared exemptions from the unit-shadow check (rule 3). The value must be
 * a genuinely different physical quantity that merely sits near a contract
 * constant.
 */
const ALLOW_UNIT_SHADOW = [];

/** Contract constants that no module may re-derive under another name. */
const CONTRACT = {
  KB_KCAL: units.KB_KCAL,
  KCONV: units.KCONV,
  STANDARD_VOLUME: units.STANDARD_VOLUME,
  COULOMB_CONST: units.COULOMB_CONST,
};
/**
 * Relative window for rule 3. 1e-3 is far wider than any legitimate
 * coincidence at these magnitudes (332.0 is 1.9e-4 away and MUST be caught)
 * and far narrower than the nearest genuinely-distinct tuning constant in
 * src/ (MEMBRANE_WIDTH=2.0, R_SWITCH_ON=6.5, GB_RADII — none within 1e-3
 * relative of a contract value).
 */
const SHADOW_REL_TOL = 1e-3;

console.log("=== constant ledger: one SCREAMING_CASE name, one home ===");

// ── 0. units.js is the single home and stays importable ──────────────
{
  const raw = fs.readFileSync(path.join(ROOT, "src", "units.js"), "utf-8");
  const code = blankCommentsAndStrings(raw);
  const hasImports = /^\s*(import|export)\s[^;]*\bfrom\b/m.test(code);
  if (hasImports) console.error("      src/units.js has an import — it must stay the leaf");
  assert(!hasImports,
    `units.js is the leaf of the dependency graph (zero static imports), so ` +
    `every re-export edge added here is acyclic by construction`);
  assert(units.KB_KCAL === 0.001987204, `KB_KCAL = ${units.KB_KCAL}`);
  assert(units.KCONV === 418.4, `KCONV = ${units.KCONV}`);
  assert(units.STANDARD_VOLUME === 1660.54, `STANDARD_VOLUME = ${units.STANDARD_VOLUME}`);
  assert(units.COULOMB_CONST === 332.06371, `COULOMB_CONST = ${units.COULOMB_CONST}`);
}

// ── 1 + 2. Single-site, and divergence is never excusable ─────────────
const decls = scanDeclarations();
const byName = new Map();
for (const d of decls) {
  if (!byName.has(d.name)) byName.set(d.name, []);
  byName.get(d.name).push(d);
}

const multiSite = [...byName.entries()].filter(([, list]) =>
  new Set(list.map((d) => d.file)).size >= 2);

let divergent = 0, unallowed = 0;
for (const [name, list] of multiSite) {
  const vals = [...new Set(list.map((d) => d.literal))];
  const sites = list.map((d) => `${d.file}:${d.line}=${d.literal}`).join(", ");
  if (vals.length > 1) {
    // UNCONDITIONAL FAILURE. A divergent pair is a live physics bug.
    divergent++;
    console.error(`      DIVERGENT ${name}: ${sites}`);
    continue;
  }
  const allowed = ALLOW_MULTI_SITE.find((e) => e.name === name);
  if (!allowed) {
    unallowed++;
    console.error(`      DUPLICATE ${name}: ${sites}`);
  } else if (!allowed.why || allowed.why.length < 40) {
    unallowed++;
    console.error(`      ALLOWLIST entry ${name} has no substantive justification`);
  }
}
assert(divergent === 0,
  `no SCREAMING_CASE constant is declared with two different literals; found ${divergent}`);
assert(unallowed === 0,
  `every multi-site SCREAMING_CASE constant is justified in ALLOW_MULTI_SITE; ` +
  `${unallowed} unaccounted for (of ${multiSite.length} multi-site names)`);

// ── 2b. The allowlist must not rot: every entry must still match ──────
{
  let stale = 0;
  for (const e of ALLOW_MULTI_SITE) {
    const list = byName.get(e.name) || [];
    const files = new Set(list.map((d) => d.file));
    if (files.size < 2) {
      stale++;
      console.error(`      STALE ALLOWLIST entry "${e.name}" — now declared in ${files.size} file(s); delete it`);
    }
  }
  assert(stale === 0, `no stale ALLOW_MULTI_SITE entries; ${stale} found`);

  // ALLOW_CONTRACT_ALIAS and ALLOW_UNIT_SHADOW are exemptions too, and the
  // file header promises "every entry is asserted to be STILL IN USE". Until
  // 2026-10-03 only the two lists above were rot-checked, so an alias or a
  // unit-shadow exemption could outlive its justification while the header
  // claimed the opposite. Both lists are checked here for the same reason: an
  // entry for a condition that no longer exists is a permanent hole, and this
  // is where the hole gets closed.
  const aliasOrShadow = [
    ...ALLOW_CONTRACT_ALIAS.map((e) => ({
      key: `${e.name} (alias of ${e.aliases})`,
      live: byName.has(e.name) && byName.get(e.name).some((d) => d.file === path.join("src", "units.js")),
    })),
    ...ALLOW_UNIT_SHADOW.map((e) => ({
      key: `${e.name}@${e.file}`,
      live: decls.some((d) => d.name === e.name && d.file === path.join("src", e.file)),
    })),
  ];
  const deadExemptions = aliasOrShadow.filter((e) => !e.live);
  for (const e of deadExemptions) {
    console.error(`      STALE EXEMPTION "${e.key}" — the declaration it excuses no longer exists; delete it`);
  }
  assert(deadExemptions.length === 0,
    `every ALLOW_CONTRACT_ALIAS / ALLOW_UNIT_SHADOW entry still excuses a live declaration ` +
    `(${aliasOrShadow.length} entries checked); ${deadExemptions.length} stale`);
}

// ── 3. Unit-shadowing: no re-derivation under a different name ────────
{
  const contractEntries = Object.entries(CONTRACT);
  const shadowAllow = new Set(ALLOW_UNIT_SHADOW.map((e) => `${e.name}@${e.file}`));
  const aliasOk = new Set(ALLOW_CONTRACT_ALIAS.map((e) => e.name));
  const UNITS_REL = path.join("src", "units.js");
  let shadows = 0;
  for (const d of decls) {
    for (const [cname, cval] of contractEntries) {
      if (d.name === cname) continue;              // the contract's own name
      const v = Number(d.literal);
      if (cval === 0) continue;
      const rel = Math.abs(v - cval) / Math.abs(cval);
      if (rel > SHADOW_REL_TOL) continue;
      // Inside units.js a second name is a documented alias (allowlisted).
      if (d.file === UNITS_REL && aliasOk.has(d.name)) continue;
      if (shadowAllow.has(`${d.name}@${d.file}`)) continue;
      shadows++;
      console.error(
        `      SHADOW ${d.name}@${d.file}:${d.line} = ${d.literal} re-derives ` +
        `${cname} = ${cval} (relative ${rel.toExponential(3)})`);
    }
  }
  assert(shadows === 0,
    `no module re-derives a units.js contract constant under another name ` +
    `(tolerance ${SHADOW_REL_TOL}); found ${shadows}`);
}

// ── 3b. GLSL/WGSL template text: a bare literal there is invisible to ──
//        the AST scan above but is exactly the 332.0-in-gpu.js case. ───
{
  let shaderHits = 0, scanned = 0;
  for (const f of srcFiles()) {
    const raw = fs.readFileSync(f, "utf-8");
    // A file only counts as shader-bearing if it has an actual shader
    // entry-point attribute in its source. Matching on the word "WGSL" alone
    // would flag prose (units.js documents the whole point of this check).
    if (!/@compute\b|@vertex\b|@fragment\b/.test(raw)) continue;
    // Comments are blanked; string/template bodies are KEPT, because a
    // literal inside a WGSL template is exactly the thing being hunted.
    const text = blankCommentsOnly(raw);
    scanned++;
    const rel = path.relative(ROOT, f);
    // Scan the WHOLE raw text (not the blanked code) for numeric literals
    // that sit within SHADOW_REL_TOL of a contract constant. Matching on the
    // exact contract literal is not enough: the real defect was a shader
    // carrying `332.0`, which is 1.9e-4 away and would slip past an
    // equality test. Interpolated `${COULOMB_CONST}` produces no literal
    // here, so the correct form passes.
    for (const m of text.matchAll(/(?<![\w.$])(-?\d+\.?\d*(?:[eE][-+]?\d+)?)(?![\w.])/g)) {
      const v = Number(m[1]);
      for (const [cname, cval] of Object.entries(CONTRACT)) {
        if (!Number.isFinite(v) || cval === 0) continue;
        const rel2 = Math.abs(v - cval) / Math.abs(cval);
        if (rel2 > SHADOW_REL_TOL) continue;
        shaderHits++;
        console.error(
          `      SHADER LITERAL ${rel}: ${m[1]} re-derives ${cname} = ${cval} ` +
          `(relative ${rel2.toExponential(3)}) at offset ${m.index} — interpolate it`);
      }
    }
  }
  assert(shaderHits === 0,
    `no units.js contract constant is hardcoded in shader source ` +
    `(${scanned} shader-bearing file(s) scanned); found ${shaderHits}`);
}

// ── 5. TABLE ledger: an element parameter table has one home ──────────────
// Added by goal M7. Rules 1–4 cannot see an object-literal table, which is
// precisely where a chemical element's LJ parameters live.
{
  const tables = scanElementParamTables();
  const sites = tables.map((t) => `${t.name}@${t.file}:${t.line}(${t.fields.join("/")})`);
  const allowed = new Map(ALLOW_ELEMENT_TABLE.map((e) => [`${e.name}@${e.file}`, e]));
  const outside = [];

  for (const t of tables) {
    if (t.file === PARAM_HOME) continue;                 // the canonical home
    const key = `${t.name}@${t.file}`;
    const ex = allowed.get(key);
    if (!ex) {
      outside.push(t);
      console.error(
        `      ELEMENT PARAM TABLE ${key}: object literal declares ` +
        `${t.fields.join("/")} outside ${PARAM_HOME} — a per-element LJ/charge ` +
        `table (or a projection of one) may only be declared in the canonical ` +
        `parameter module, or allowlisted as a DIFFERENT physical model`);
    } else if (!ex.why || ex.why.length < 40) {
      outside.push(t);
      console.error(`      ALLOW_ELEMENT_TABLE entry ${key} has no substantive justification`);
    }
  }
  assert(outside.length === 0,
    `every per-element sigma/eps/q table is declared in ${PARAM_HOME} or carries a ` +
    `written justification for being a different model; ${tables.length} table site(s) ` +
    `scanned, ${outside.length} unaccounted for` +
    (tables.length ? ` — ${sites.join(", ")}` : ""));

  // The allowlist must not rot, same as ALLOW_MULTI_SITE: an exemption whose
  // table no longer exists is a silent hole and is itself a failure.
  let stale = 0;
  for (const e of ALLOW_ELEMENT_TABLE) {
    const key = `${e.name}@${e.file}`;
    if (!tables.some((t) => `${t.name}@${t.file}` === key)) {
      stale++;
      console.error(`      STALE ALLOW_ELEMENT_TABLE entry "${key}" — no such table; delete it`);
    }
  }
  assert(stale === 0, `no stale ALLOW_ELEMENT_TABLE entries; ${stale} found`);

  // The canonical module must actually BE canonical: if someone empties it, the
  // single-site rule above becomes vacuous and the duplication guard guards
  // nothing. This is the vacuity check the repo already insists on elsewhere
  // (scripts/check.sh per-extension totals, test_suite_registry's census).
  const params = await import("../src/physics/params.js");
  // Same leaf discipline as units.js (rule 0): the canonical parameter module
  // must import nothing, so ff-params.js, forcefield.js, heavy.js,
  // ligand-panel.js and placement.js can all reach it without a cycle, and so a
  // second literal cannot arrive "via a dependency".
  const paramsRaw = fs.readFileSync(path.join(ROOT, "src", "physics", "params.js"), "utf-8");
  const paramsHasImports = /^\s*(import|export)\s[^;]*\bfrom\b/m.test(blankCommentsAndStrings(paramsRaw));
  assert(!paramsHasImports,
    `${PARAM_HOME} is a leaf too (zero static imports) — every consumer can reach ` +
    `the parameter contract without creating a cycle`);
  assert(Object.keys(params.ELEMENT_LJ).length === 9 &&
    Object.keys(params.METAL_ELEMENT).length === 10 &&
    Object.keys(params.RES_CLASS).length === 5 &&
    Object.keys(params.COVALENT_RADIUS).length === 21,
    `the canonical parameter module is populated, not hollow: ` +
    `ELEMENT_LJ ${Object.keys(params.ELEMENT_LJ).length}, ` +
    `METAL_ELEMENT ${Object.keys(params.METAL_ELEMENT).length}, ` +
    `RES_CLASS ${Object.keys(params.RES_CLASS).length}, ` +
    `COVALENT_RADIUS ${Object.keys(params.COVALENT_RADIUS).length}`);

  // …and the old home must be a re-export, not a second declaration. This is
  // the assertion that fails if HEAVY_ELEMENT_DEFAULT is ever planted again.
  const facade = await import("../src/ff-params.js");
  assert(facade.LIG_ELEMENT === params.ELEMENT_LJ &&
    facade.LIG_ELEMENT_DEFAULT === params.ELEMENT_LJ_DEFAULT &&
    facade.METAL_ELEMENT === params.METAL_ELEMENT &&
    facade.RES_CLASS === params.RES_CLASS &&
    facade.COVALENT_RADIUS === params.COVALENT_RADIUS,
    `src/ff-params.js re-exports the canonical tables by IDENTITY (no copies to drift)`);
}

// ── 4. Runtime: every re-export site is bit-identical to units.js ────
{
  // Each entry lists the contract names that module is expected to
  // re-export. A module that re-exports a name MUST return the units.js
  // value; a module that does not list a name is not asserted on it.
  const SITES = [
    ["../src/physics/gb.js", ["COULOMB_CONST"]],
    ["../src/physics/solvation/gb_obc2.js", ["COULOMB_CONST"]],
    ["../src/physics/solvation/membrane_slab.js", ["COULOMB_CONST"]],
    ["../src/compute/webgpu_backend.js", ["COULOMB_CONST"]],
    ["../src/ff-params.js", ["COULOMB_CONST", "KB_KCAL"]],  // KB_KCAL asserted by
    //   test_unit_contract.js; ff-params re-exports KB_KCAL + KCONV from units.js.
  ];
  for (const [spec, names] of SITES) {
    const m = await import(spec);
    const label = spec.replace("../src/", "");
    const exposed = names.filter((n) => n in m);
    for (const n of names) {
      if (!(n in m)) continue;   // not re-exported here; nothing to compare
      assert(m[n] === units[n], `${label}.${n} re-exports units.js (${m[n]})`);
    }
    if (exposed.length === 0) {
      // ff-params.js deliberately does not re-export COULOMB_CONST; assert the
      // KB_KCAL/KCONV pair it DOES re-export instead.
      const ff = m;
      assert(ff.KB_KCAL === units.KB_KCAL && ff.KCONV === units.KCONV,
        `${label} re-exports units.js KB_KCAL/KCONV (${ff.KB_KCAL}, ${ff.KCONV})`);
    }
  }
  const { PROBE_RADIUS: prA } = await import("../src/physics/sasa.js");
  const { PROBE_RADIUS: prB } = await import("../src/physics/solvation/lcpo_sasa.js");
  assert(prA === prB, `both SASA models share one PROBE_RADIUS (${prA} === ${prB})`);
}

// ── 4c. A module that READS a contract constant must hold a LOCAL binding
//        for it. `export { X } from "./y.js"` publishes X to importers but
//        creates no local binding, so a body that references X throws
//        ReferenceError at CALL time — long after an import-only check
//        passes. Consolidating COULOMB_CONST hit exactly this: the first
//        version of this file asserted only on import, and the suite died
//        later inside GeneralizedBorn.pairInteraction.
//
//        Local bindings come from `import { X } from "…"` (the correct
//        form) or from a `const X = …` declaration (already covered by the
//        shadow rule). A re-export clause is neither.
// ──────────────────────────────────────────────────────────────────────
{
  let missingBinding = 0, checkedPairs = 0;
  for (const f of srcFiles()) {
    const rel = path.relative(ROOT, f);
    if (rel === path.join("src", "units.js")) continue;
    const code = blankCommentsAndStrings(fs.readFileSync(f, "utf-8"));

    for (const cname of Object.keys(CONTRACT)) {
      const NAME_RE = `(?<![\\w.$])${cname}(?![\\w$])`;

      // Blank only the parts that BIND nothing locally:
      //   - import clauses          (`import { X } from "…"`)
      //   - re-export clauses      (`export { X } from "…"`, `export { X };`)
      //   - declaration INITIALISERS (`const y = X` — X on the right-hand
      //     side IS a real read of the name, so only the `= …` tail goes).
      // Statements themselves must survive: `const a = X / b;` is a read.
      const body = code
        .replace(/import\s+[^;]*?from\s*["'][^"']*["']\s*;?/g, " ")
        .replace(/import\s*["'][^"']*["']\s*;?/g, " ")
        .replace(/export\s*\{[^}]*\}\s*(from\s*["'][^"']*["'])?\s*;?/g, " ")
        // A declaration of the name ITSELF binds it; a declaration of some
        // other name does not, so keep those bodies intact.
        .replace(new RegExp(`(export\\s+)?(?:const|let|var)\\s+${cname}\\s*=[^;\\n]*`, "g"), " ");

      const reads = new RegExp(NAME_RE, "g").test(body);
      if (!reads) continue;
      checkedPairs++;

      const hasLocalImport = new RegExp(
        `import\\s*\\{[^}]*\\b${cname}\\b[^}]*\\}\\s*from`).test(code);
      const hasLocalDecl = new RegExp(
        `(?:export\\s+)?(?:const|let|var)\\s+${cname}\\s*=`).test(code);
      const hasNamespaceOrDefault = new RegExp(
        `import\\s+\\*\\s+as\\s+\\w+\\s+from`).test(code);

      if (!hasLocalImport && !hasLocalDecl && !hasNamespaceOrDefault) {
        missingBinding++;
        console.error(
          `      NO LOCAL BINDING ${rel}: ${cname} is read in the module body but ` +
          `only re-exported — \`export { ${cname} } from …\` publishes it without ` +
          `binding it, so this throws ReferenceError at call time`);
      }
    }
  }
  assert(missingBinding === 0,
    `every contract constant a module READS has a local binding ` +
    `(${checkedPairs} module/name pairs checked); ${missingBinding} missing`);
}

// ── 4b. The specific regression: thermodynamics.js entropy must use the
//        contract k_B, not a private copy. This is the live bug the scan
//        found AFTER test_unit_contract.js was written. ───────────────
{
  const { schlitterEntropy, torsionEntropy } = await import("../src/analysis/thermodynamics.js");
  const cov = [[0.30, 0.11, 0.05, 0.02], [0.11, 0.28, 0.09, 0.04],
    [0.05, 0.09, 0.26, 0.10], [0.02, 0.04, 0.10, 0.31]];
  const got = schlitterEntropy(cov, [12, 12, 14, 1], 300);

  // Reference implementation written against units.KB_KCAL directly.
  const KB = units.KB_KCAL, e = Math.E, T = 300;
  const m = [12, 12, 14, 1];
  const A = cov.map((row, i) => row.map((v, j) =>
    (KB * T * e * e / (KB * KB)) * v * Math.sqrt(m[i] * m[j])));
  const M = A.map((row, i) => row.map((v, j) => (i === j ? 1 + v : v)));
  // Cholesky ln-det, computed independently of the module under test.
  const Mc = M.map((r) => r.slice());
  let lnDet = 0;
  for (let i = 0; i < Mc.length; i++) {
    let d = Mc[i][i];
    for (let k = 0; k < i; k++) d -= Mc[i][k] * Mc[i][k];
    lnDet += Math.log(d);
    for (let j = i + 1; j < Mc.length; j++) {
      let t = Mc[i][j];
      for (let k = 0; k < i; k++) t -= Mc[i][k] * Mc[j][k];
      Mc[j][i] = t / d;
    }
  }
  const ref = (KB / 2) * 2 * lnDet;

  assert(got === ref,
    `schlitterEntropy uses units.js KB_KCAL exactly (${got} === ${ref})`);

  // And quantify the defect that was there, so a regression is legible.
  const KB_OLD = 0.0019872041;
  const refOld = (KB_OLD / 2) * 2 * lnDet;
  const relSplit = Math.abs(refOld - ref) / Math.abs(ref);
  assert(relSplit > 1e-8 && relSplit < 1e-6,
    `the reintroduced-split magnitude is ${relSplit.toExponential(3)} ` +
    `(KB 0.0019872041 vs ${units.KB_KCAL}) — nonzero, and the assertion ` +
    `above is exact, so a private copy cannot pass`);

  const frames = [];
  for (let f = 0; f < 24; f++) {
    const th = (f * 2 * Math.PI) / 24;
    frames.push(new Float64Array([0, 0, 0, 1, 0, 0, 1, 1,
      Math.cos(th) * 0.5 + 0.5, Math.sin(th) * 0.5, 0, 1, 0]));
  }
  const dS = torsionEntropy(frames, [[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]], 30);
  assert(Number.isFinite(dS) && dS > 0,
    `torsionEntropy finite and non-degenerate (${dS.toExponential(4)} kcal/mol/K)`);
}

console.log(`=== test_constant_ledger: ${passed} PASSED, ${failed} FAILED ===`);
process.exit(failed > 0 ? 1 : 0);
