// validate_binding_physics_r1.mjs — zero-dep check of every src:line claim in
// docs/BINDING_PHYSICS_R1.md §0. Run: node scripts/validate_binding_physics_r1.mjs
//
// Registered in tests/suites.js (FAST tier) so it actually gates; a validator
// that can fail silently is decoration.
//
// WHY THE INVENTORY BLOCK BELOW IS SPLIT
// ---------------------------------------
// This file originally asserted, as one flat negative list, that heavy.js
// "lacks" halogen / cation-pi / pi-stack / sigma-hole. That was TRUE when
// BINDING_PHYSICS_R1.md was written (step 1 of 7) and became FALSE when phase
// R3 landed `src/physics/weakint.js` and wired π-stack, cation-π and halogen
// σ-hole into the heavy force loop behind `par.weak === "on"` (heavy.js:850,
// 1097-1104, 1391-1459). The validator's premise, not the code, was stale.
//
// Two things follow, and both are enforced below:
//
//   1. The three terms heavy.js DOES implement are now asserted POSITIVELY —
//      and positively means "the kernel is imported, the pair lists are built,
//      the method exists, all three kernels are called inside it, and the
//      result is added to U". A comment or a docstring cannot satisfy that.
//   2. The terms still genuinely absent (chalcogen, PME, explicit water) stay
//      negative, but the absence test now runs against CODE ONLY — comments,
//      docstrings and string literals are blanked first. The old raw-substring
//      test could not tell a prose mention from an implementation; that is
//      exactly how this validator rotted in the first place.
//
// stripJs() over-stripping is caught by a self-check: the three real kernel
// call sites must survive stripping, or the validator fails loudly. Under-
// stripping is caught the other way — a surviving comment makes an absence
// assertion fail loudly. Neither direction can pass silently.
import { readFileSync } from "node:fs";
import { strict as assert } from "node:assert";

const R = (p) => readFileSync(new URL(p, import.meta.url), "utf8");
const has = (s, re, msg) => assert.match(s, re, msg);

let passed = 0;
let failed = 0;
/** Counted assertion that reports and throws on failure (harness counts these). */
function check(cond, msg) {
  if (cond) { passed++; return; }
  failed++;
  assert.fail(msg);
}

// ---------------------------------------------------------------------
// stripJs(src, keepStrings) — blank out comments and (by default) string /
// template literal contents. Newlines are preserved so reported line numbers
// still line up with the original file; every blanked character becomes a
// space otherwise.
//
// Two views are used, and the difference matters:
//   stripJs(s, false) — strings blanked  -> the CODE view. Every
//       "is this implemented?" question is asked here, so a comment or a
//       docstring can neither satisfy a positive assertion nor trip an
//       absence assertion.
//   stripJs(s, true)  — strings kept     -> used only for module-specifier
//       checks (`from "./physics/weakint.js"`), where the specifier is
//       legitimately a string and must stay readable. Comments are still
//       blanked, so a faked import cannot hide in one.
// ---------------------------------------------------------------------
const REGEX_OK_BEFORE = new Set([
  "(", "[", "{", ",", ";", ":", "=", "!", "&", "|", "?", "+", "-", "*",
  "%", "~", "^", "<", ">", "}", ">", "\n",
]);
const REGEX_KEYWORD_BEFORE = new Set([
  "return", "typeof", "instanceof", "in", "of", "new", "delete", "void",
  "do", "else", "case", "yield", "await", "throw",
]);

function stripJs(src, keepStrings = false) {
  let out = "";
  let i = 0;
  const n = src.length;
  let prev = ""; // last non-space code token: a punctuation char or a word
  const isWord = (c) => /[A-Za-z0-9_$]/.test(c);
  const take = (c) => {
    if (!/\s/.test(c)) {
      if (isWord(c)) prev = /[A-Za-z0-9_$]/.test(prev) ? prev + c : c;
      else prev = c;
    }
    out += c;
  };
  while (i < n) {
    const c = src[i], c2 = src[i + 1];
    if (c === "/" && c2 === "/") {
      while (i < n && src[i] !== "\n") { out += " "; i++; }
      continue;
    }
    if (c === "/" && c2 === "*") {
      out += " "; i += 2;
      while (i < n && !(src[i] === "*" && src[i + 1] === "/")) {
        out += src[i] === "\n" ? "\n" : " "; i++;
      }
      out += " "; i += 2;
      continue;
    }
    if (keepStrings) { take(c); i++; continue; }
    if (c === '"' || c === "'" || c === "`") {
      const q = c;
      out += " "; i++;
      while (i < n) {
        if (src[i] === "\\") { out += "  "; i += 2; continue; }
        if (src[i] === q) break;
        if (q === "`" && src[i] === "$" && src[i + 1] === "{") {
          // Template substitution is REAL code — keep it, blank the literal parts.
          out += "${"; i += 2;
          let depth = 1;
          while (i < n && depth > 0) {
            if (src[i] === "{") depth++;
            else if (src[i] === "}") { depth--; if (depth === 0) { out += "}"; i++; break; } }
            take(src[i]); i++;
          }
          continue;
        }
        out += src[i] === "\n" ? "\n" : " "; i++;
      }
      out += " "; i++;
      continue;
    }
    const regexOk =
      c === "/" &&
      (prev === "" || REGEX_OK_BEFORE.has(prev) || REGEX_KEYWORD_BEFORE.has(prev));
    if (regexOk) {
      out += " "; i++;
      let inClass = false;
      while (i < n) {
        if (src[i] === "\\") { out += "  "; i += 2; continue; }
        if (src[i] === "[") inClass = true;
        else if (src[i] === "]") inClass = false;
        else if (src[i] === "/" && !inClass) break;
        else if (src[i] === "\n") break; // not a regex after all
        out += " "; i++;
      }
      out += " "; i++;
      while (i < n && /[a-z]/.test(src[i])) { out += " "; i++; }
      continue;
    }
    take(c);
    i++;
  }
  return out;
}

const ffB = R("../src/ff-binding.js");
const ff = R("../src/forcefield.js");
const ffp = R("../src/ff-params.js");
const heavy = R("../src/heavy.js");
const hb = R("../src/physics/hbond.js");
const gb = R("../src/physics/gb.js");
const sasa = R("../src/physics/sasa.js");
const prot = R("../src/chem/protonation.js");
const r1doc = R("../docs/BINDING_PHYSICS_R1.md");
const weak = R("../src/physics/weakint.js");

// CODE view of heavy.js — comments/docstrings/strings blanked. Every
// "does heavy.js implement X" question below is asked of THIS, not of the raw
// text, so prose can no longer masquerade as implementation (or vice versa).
const heavyCode = stripJs(heavy, false);
// Module-specifier view of heavy.js — comments blanked, strings kept, so
// `from "./physics/weakint.js"` stays readable without letting a comment
// fabricate it.
const heavyImports = stripJs(heavy, true);

// =====================================================================
// CG kernel constants
// =====================================================================
has(ffB, /EPSHB\s*=\s*0\.8/, "EPSHB=0.8");
has(ffB, /HB_R0\s*=\s*3\.2/, "HB_R0=3.2");
has(ffB, /HB_W\s*=\s*0\.6/, "HB_W=0.6");
has(ffB, /4\s*\+\s*76\s*\*\s*Math\.tanh\(r\s*\/\s*8\)/, "eps(r)=4+76 tanh(r/8)");
has(ffB, /R0\s*=\s*4\.5,\s*SIG\s*=\s*1\.8,\s*NS\s*=\s*3\.0/, "R0/SIG/NS");
has(ffB, /_pairKey\(i,\s*la\)[\s\S]{0,120}_excluded\.has\(pk\)/, "holo exclusion");
has(ff, /bindRcut\s*=\s*9\.0/, "bindRcut=9.0");
has(ff, /holoGamma.*0\.5/, "holoGamma=0.5");
has(ff, /_protHB\[i\]\s*=\s*\(cls\s*===\s*"P"/, "protHB flags");
// Protein charges all zero in CG table
has(ffp, /H:\s*\{\s*sigma:\s*4\.0,\s*eps:\s*0\.15,\s*q:\s*0\s*\}/, "RES_CLASS H q=0");
has(ffp, /Cp:\s*\{\s*sigma:\s*3\.6,\s*eps:\s*0\.10,\s*q:\s*0\s*\}/, "RES_CLASS Cp q=0");
// Ligand dG range
for (const v of ["-0.55", "-0.35", "-0.30", "-0.25", "-0.40", "-0.45", "-0.50"])
  check(ffp.includes(v), `lig dG ${v} present`);

// =====================================================================
// Heavy kernel inventory — PRESENT terms, asserted as real code paths
// =====================================================================
// Base inventory: LJ + GB + H-bond. Unchanged.
for (const t of ["_nonBondedGrid", "pairInteraction", "evaluatePair", "GeneralizedBorn", "SasaModel"])
  check(heavyCode.includes(t), `heavy has ${t}`);

// R3 weak interactions ARE implemented in heavy mode. Each link in the chain
// is checked, so removing any single link turns this validator red.
check(
  /import\s*\{[^}]*piStackForces[^}]*cationPiForces[^}]*halogenForces[^}]*\}\s*from\s*"\.\/physics\/weakint\.js/.test(heavyImports),
  "heavy imports piStackForces/cationPiForces/halogenForces from physics/weakint.js",
);
check(
  /this\._weakRings\s*=\s*buildRingFrames\(/.test(heavyCode) &&
  /this\._weakCations\s*=\s*buildCationList\(/.test(heavyCode) &&
  /this\._weakHalogens\s*=\s*buildHalogenList\(/.test(heavyCode),
  "heavy builds ring / cation / halogen pair lists once per topology",
);
check(/_weakInteractions\s*\(\s*pos\s*,\s*f\s*\)\s*\{/.test(heavyCode), "heavy defines _weakInteractions(pos, f)");

// The method body must contain all three kernel calls, and all three must live
// inside the method (not merely somewhere in the file).
const weakBody = heavyCode.match(/_weakInteractions\s*\(\s*pos\s*,\s*f\s*\)\s*\{([\s\S]*?)\n {2}\}/);
check(weakBody !== null, "_weakInteractions body is extractable");
const wb = weakBody ? weakBody[1] : "";
check(/\bpiStackForces\s*\(/.test(wb), "_weakInteractions calls piStackForces (pi-stack)");
check(/\bcationPiForces\s*\(/.test(wb), "_weakInteractions calls cationPiForces (cation-pi)");
check(/\bhalogenForces\s*\(/.test(wb), "_weakInteractions calls halogenForces (halogen sigma-hole)");

// And the result must reach the energy, not be computed and dropped.
check(
  /const\s+w\s*=\s*this\._weakInteractions\(pos,\s*f\)/.test(heavyCode) &&
  /U\s*\+=\s*this\.weakU/.test(heavyCode),
  "weak-interaction energy is added to U in the force loop",
);

// The honest-scope fact that makes the rest of the docs true: opt-in and OFF
// by default. If this ever becomes default-on, the ROADMAP/LIMITATIONS wording
// has to be revisited, and this validator is where that shows up. (Checked on
// the strings-kept view because the sentinel is the string literal "on".)
check(/this\.weakOn\s*=\s*par\.weak\s*===\s*"on"/.test(heavyImports), 'heavy weak terms are opt-in (par.weak === "on")');
check(/if\s*\(\s*this\.weakOn\s*\)/.test(heavyCode), "force loop calls the weak pass only when weakOn");

// The kernels themselves exist in their own module and export forces.
for (const t of ["piStackEnergy", "cationPiEnergy", "halogenEnergy", "HALOGEN_EPS",
                 "piStackForces", "cationPiForces", "halogenForces",
                 "buildRingFrames", "buildCationList", "buildHalogenList"])
  check(new RegExp(`export\\s+(?:const|function)\\s+${t}\\b`).test(weak), `weakint.js exports ${t}`);

// =====================================================================
// Heavy kernel inventory — ABSENT terms, asserted against CODE ONLY
// =====================================================================
const heavyCodeLC = heavyCode.toLowerCase();
for (const t of ["chalcogen", "pme", "explicit water"])
  check(!heavyCodeLC.includes(t), `heavy lacks ${t} (code view — comments/docstrings/strings excluded)`);

// Self-check on stripJs: the real kernel call sites must survive stripping.
// If this ever fails, the stripper is eating code and every absence assertion
// above would have been vacuous.
for (const t of ["piStackForces(pos, f, ra, rb)", "cationPiForces(pos, f, ci, ring)",
                 "halogenForces(pos, f, c, x, d"])
  check(heavyCode.includes(t), `stripJs self-check: ${t} survives comment stripping`);

// =====================================================================
// GB/SASA constants and the rest of the R1 §0 inventory (unchanged)
// =====================================================================
check(heavy.includes("epsIn: 4.0") && heavy.includes("epsOut: 78.5"), "GB dielectrics");
has(sasa, /gamma.*0\.0072/, "SASA gamma");
has(hb, /HBOND_EQ_DIST\s*=\s*2\.9/, "hbond r_eq");
has(hb, /evaluatePair[\s\S]{0,300}180/, "evaluatePair ideal 180");
check(hb.includes("angular gradient not yet propagated") || hb.includes("angular gradient not propagated"), "hbond stub note");
has(heavy, /BOND_SLACK\s*\*\s*\(rA\s*\+\s*rB\)/, "covalent rule");
check(heavy.includes("2.2"), "2.2 cap");
has(prot, /SALTBRIDGE_DIST\s*=\s*4\.5/, "salt bridge cutoff");
has(prot, /HBOND_DIST\s*=\s*3\.5/, "hbond dist cutoff");
check(gb.includes("Hawkins"), "GB HCT provenance");

// =====================================================================
// Doc consistency — the premise this validator certifies
// =====================================================================
// BINDING_PHYSICS_R1.md §0 states, in the present tense, that the heavy kernel
// has "no halogen / cation-π / π-stack / chalcogen terms anywhere". R3 made the
// first three of those false. The doc is a dated R1 snapshot now and must say
// so, or the honest-scope record contradicts the code.
check(
  /as of R1|point-in-time|snapshot/i.test(r1doc) &&
  /BINDING_PHYSICS_R3\.md/.test(r1doc),
  "BINDING_PHYSICS_R1.md is labelled a dated snapshot superseded by BINDING_PHYSICS_R3.md",
);
// The R1 snapshot must defer to the LIVE scope guard, and that guard must
// still be intact. This is the one place the R1 validator touches ROADMAP.md
// §1 — R3's weak terms do not violate it, because they are opt-in and are
// not QM/MM, a lipid bilayer, or PME.
const roadmap = R("../ROADMAP.md");
check(/ROADMAP\.md/.test(r1doc), "R1 snapshot defers to ROADMAP.md as the live scope guard");
check(
  /No QM\/MM, no explicit membrane, no PME/.test(roadmap),
  "ROADMAP.md still carries its live scope guard (No QM/MM, no explicit membrane, no PME)",
);
// docs/LIMITATIONS.md is the second live honest-scope surface.
const limitations = R("../docs/LIMITATIONS.md");
check(
  /no PME/i.test(limitations) && /no membrane/i.test(limitations),
  "LIMITATIONS.md still discloses no PME and no membrane",
);

if (failed > 0) {
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  assert.fail(`${failed} R1 code-citation assertion(s) failed`);
}
console.log(`R1 validation: all code citations OK (${passed} assertions)`);
console.log(`TEST RESULTS: ${passed} PASSED, 0 FAILED`);
