// validate_binding_physics_r1.mjs — zero-dep check of every src:line claim in
// docs/BINDING_PHYSICS_R1.md §0. Run: node scripts/validate_binding_physics_r1.mjs
import { readFileSync } from "node:fs";
import { strict as assert } from "node:assert";

const R = (p) => readFileSync(new URL(p, import.meta.url), "utf8");
const has = (s, re, msg) => assert.match(s, re, msg);

const ffB = R("../src/ff-binding.js");
const ff = R("../src/forcefield.js");
const ffp = R("../src/ff-params.js");
const heavy = R("../src/heavy.js");
const hb = R("../src/physics/hbond.js");
const gb = R("../src/physics/gb.js");
const sasa = R("../src/physics/sasa.js");
const prot = R("../src/chem/protonation.js");

// CG kernel constants
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
  assert.ok(ffp.includes(v), `lig dG ${v} present`);
// Heavy kernel inventory: LJ + GB + H-bond present; halogen/cation-pi/pistack absent
for (const t of ["_nonBondedGrid", "pairInteraction", "evaluatePair", "GeneralizedBorn", "SasaModel"])
  assert.ok(heavy.includes(t), `heavy has ${t}`);
for (const t of ["halogen", "cation-pi", "cationPi", "pi-stack", "pistack", "sigma-hole", "sigmahole", "chalcogen", "PME", "explicit water"])
  assert.ok(!heavy.toLowerCase().includes(t.toLowerCase()), `heavy lacks ${t}`);
// GB/SASA constants
assert.ok(heavy.includes("epsIn: 4.0") && heavy.includes("epsOut: 78.5"), "GB dielectrics");
has(sasa, /gamma.*0\.0072/, "SASA gamma");
has(hb, /HBOND_EQ_DIST\s*=\s*2\.9/, "hbond r_eq");
has(hb, /evaluatePair[\s\S]{0,300}180/, "evaluatePair ideal 180");
assert.ok(hb.includes("angular gradient not yet propagated") || hb.includes("angular gradient not propagated"), "hbond stub note");
has(heavy, /BOND_SLACK\s*\*\s*\(rA\s*\+\s*rB\)/, "covalent rule");
assert.ok(heavy.includes("2.2"), "2.2 cap");
has(prot, /SALTBRIDGE_DIST\s*=\s*4\.5/, "salt bridge cutoff");
has(prot, /HBOND_DIST\s*=\s*3\.5/, "hbond dist cutoff");
assert.ok(gb.includes("Hawkins"), "GB HCT provenance");

console.log("R1 validation: all code citations OK");
