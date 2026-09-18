/**
 * test_l0_default_exposure.js — Revolution 1 / Issue 3 regression.
 *
 * Scope: src/forcefield.js + src/ff-binding.js + physics-level default
 * exposure ONLY. No energy-kernel physics change.
 *
 * Problem: default CG L0 path (charges OFF, hbMode "off", uniform ENM) misses
 * salt-bridge electrostatics + backbone directional H-bonds in bindingU.
 *
 * Fix verified here (minimal, backward compatible):
 *   - Default stays L0 (bit-identical legacy) but is now EXPLICIT via
 *     DEFAULT_PHYSICS_LEVEL / CG_PHYSICS_LEVELS / resolvePhysicsLevel() and
 *     ForceField.describePhysics() + ff-binding bindingTermsActive().
 *   - Opt-in L1 wiring (par.physicsLevel "L1"/"L2" ≡ charges:true +
 *     hbMode:"directional") revives the existing Coulomb + virtual-site HB
 *     paths; explicit par.binding flags always win over the tier.
 *
 * Coverage (headless, 4W52 CG + native charged ligands, ASP+LYS census):
 *   [0] default-exposure constants + resolver (L0 default, unknown→L0,
 *       explicit-flags-win, L2 CG ≡ L1).
 *   [1] charge table: L0 all-zero vs L1 ASP/GLU −1, LYS/ARG +1, HIS 0
 *       (4W52 census 18/26), Lys35 + ASP spot checks.
 *   [2] L0 vs L1 bindingU on native 4W52 CG (BNZ+EPE, charged EPE): both
 *       finite, Δ finite + nonzero + sane (|Δ|<10), |bindingU|<50 sane,
 *       forces all-finite with sane max (<10).
 *   [3] alias wiring: physicsLevel:"L1" ≡ explicit charges+directional
 *       (bit-identical); default ≡ explicit L0 (bit-identical).
 *   [4] term exposure: bindingTermsActive + trackTerms Coulomb nonzero sane.
 *
 * Run: node tests/test_l0_default_exposure.js (fast, <2s; NOT wired into
 * tests/test_all.js FAST so the 352 gate is untouched).
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import {
  ForceField,
  DEFAULT_PHYSICS_LEVEL,
  CG_PHYSICS_LEVELS,
  resolvePhysicsLevel,
} from "../src/forcefield.js";
import { bindingTermsActive, describeBindingTerms } from "../src/ff-binding.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

let passed = 0;
let failed = 0;
function assert(condition, message) {
  if (condition) {
    passed++;
    console.log(`  ✓ ${message}`);
  } else {
    failed++;
    console.error(`  ✗ FAIL: ${message}`);
  }
}

function findPdb() {
  const candidates = [
    path.resolve(process.cwd(), "4w52.pdb"),
    path.resolve(__dirname, "..", "4w52.pdb"),
    path.resolve(__dirname, "4w52.pdb"),
  ];
  for (const p of candidates) if (fs.existsSync(p)) return p;
  throw new Error(`4w52.pdb not found (tried ${candidates.join(", ")})`);
}

function maxAbsFinite(forces) {
  let m = 0;
  for (let i = 0; i < forces.length; i++) {
    const v = forces[i];
    if (!Number.isFinite(v)) return { max: Infinity, allFinite: false };
    const a = Math.abs(v);
    if (a > m) m = a;
  }
  return { max: m, allFinite: true };
}

function main() {
  console.log("=== Rev1/Issue3: L0 default exposure + L1 opt-in wiring (4W52 CG) ===");
  const pdbText = fs.readFileSync(findPdb(), "utf-8");
  const parsed = parseCa(pdbText);
  const sel = selectSystem(parsed);
  const ligs = parseLigands(pdbText);
  assert(sel.beads.length === 164, `4W52 CG beads = 164 (got ${sel.beads.length})`);
  assert(ligs.length >= 1, `native ligands present (${ligs.length} mols)`);

  // -----------------------------------------------------------------
  // [0] Default-exposure constants + resolver
  // -----------------------------------------------------------------
  console.log("\n[0] physics-level default exposure (no kernel change)...");
  assert(DEFAULT_PHYSICS_LEVEL === "L0", `DEFAULT_PHYSICS_LEVEL is "L0"`);
  assert(
    CG_PHYSICS_LEVELS.L0.charges === false && CG_PHYSICS_LEVELS.L0.hbMode === "off"
    && CG_PHYSICS_LEVELS.L0.enm === "uniform",
    `L0 tier = charges OFF / hbMode off / uniform ENM`,
  );
  assert(
    CG_PHYSICS_LEVELS.L1.charges === true && CG_PHYSICS_LEVELS.L1.hbMode === "directional",
    `L1 tier = charges ON / hbMode directional (opt-in)`,
  );
  assert(
    CG_PHYSICS_LEVELS.L2.charges === true && CG_PHYSICS_LEVELS.L2.hbMode === "directional",
    `L2 CG slice ≡ L1 (weak/bindLog live elsewhere)`,
  );
  const rDef = resolvePhysicsLevel({});
  assert(rDef.level === "L0" && rDef.charges === false && rDef.hbMode === "off",
    `resolver default → L0/off (got ${rDef.level}/${rDef.hbMode})`);
  const rBad = resolvePhysicsLevel({ physicsLevel: "L9" });
  assert(rBad.level === "L0", `unknown level falls back to L0 (got ${rBad.level})`);
  const rL1 = resolvePhysicsLevel({ physicsLevel: "L1" });
  assert(rL1.level === "L1" && rL1.charges === true && rL1.hbMode === "directional",
    `resolver L1 → charges ON / directional`);
  const rWin = resolvePhysicsLevel({ physicsLevel: "L1", binding: { charges: false, hbMode: "off" } });
  assert(rWin.charges === false && rWin.hbMode === "off",
    `explicit binding flags win over the tier (L1 + explicit off → off)`);

  // -----------------------------------------------------------------
  // [1] Charge table: L0 silent-zero vs L1 ASP+LYS formal charges
  // -----------------------------------------------------------------
  console.log("\n[1] charge table: 4W52 ASP+LYS census...");
  const ffL0 = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligs);
  const ffL1 = new ForceField(sel, { rc: 10, gamma: 1.0, physicsLevel: "L1" }, ligs);
  assert(ffL0.physicsLevel === "L0" && ffL1.physicsLevel === "L1",
    `requested tiers exposed (L0=${ffL0.physicsLevel}, L1=${ffL1.physicsLevel})`);
  assert(ffL0.chargesOn === false && ffL0.hbMode === "off",
    `L0 flags off (charges=${ffL0.chargesOn}, hb=${ffL0.hbMode})`);
  assert(ffL1.chargesOn === true && ffL1.hbMode === "directional",
    `L1 flags on (charges=${ffL1.chargesOn}, hb=${ffL1.hbMode})`);
  assert(ffL0.enmModel === "uniform" && ffL1.enmModel === "uniform",
    `ENM stays uniform unless separately opted in (both "${ffL0.enmModel}")`);
  let nNeg = 0, nPos = 0, bad = 0, nHis = 0;
  sel.beads.forEach((b, i) => {
    const q = ffL1._protQ[i];
    if (b.resName === "ASP" || b.resName === "GLU") { (q === -1) ? nNeg++ : bad++; }
    else if (b.resName === "LYS" || b.resName === "ARG") { (q === 1) ? nPos++ : bad++; }
    else { if (b.resName === "HIS" && q === 0) nHis++; (q === 0) ? 0 : bad++; }
  });
  assert(bad === 0, `L1 charges match residue identity (${nNeg} −1, ${nPos} +1, 0 mismatches)`);
  assert(nNeg === 18 && nPos === 26, `4W52 census Asp/Glu=18 Lys/Arg=26 (got ${nNeg}/${nPos})`);
  assert(nHis > 0, `HIS stays neutral by default (${nHis} His, q=0)`);
  let allZero = true;
  for (let i = 0; i < ffL0._protQ.length; i++) if (ffL0._protQ[i] !== 0) allZero = false;
  assert(allZero, `L0 _protQ all-zero (legacy silent baseline, now explicit)`);
  const lys35 = sel.beads.findIndex((b) => b.resName === "LYS" && b.resSeq === 35 && b.chain === "A");
  const asp0 = sel.beads.findIndex((b) => b.resName === "ASP");
  assert(lys35 >= 0 && asp0 >= 0, `ASP + LYS35 beads found (asp=${asp0}, lys35=${lys35})`);
  assert(ffL1._protQ[lys35] === 1 && ffL1._protQ[asp0] === -1,
    `spot check: Lys35 +1, ASP ${sel.beads[asp0].resSeq} −1 in L1`);
  assert(ffL0._protQ[lys35] === 0 && ffL0._protQ[asp0] === 0,
    `spot check: same beads 0 in L0 default`);
  const d0 = ffL0.describePhysics();
  const d1 = ffL1.describePhysics();
  assert(d0.isSimplifiedDefault === true && d0.coulombActive === false && d0.directionalHBActive === false,
    `L0 describePhysics flags simplified default (coulomb/directional OFF)`);
  assert(d1.isSimplifiedDefault === false && d1.coulombActive === true && d1.directionalHBActive === true,
    `L1 describePhysics flags opt-in wiring (coulomb/directional ON)`);

  // -----------------------------------------------------------------
  // [2] L0 vs L1 bindingU on native 4W52 CG + charged ligand (EPE)
  // -----------------------------------------------------------------
  console.log("\n[2] L0 vs L1 bindingU (native pocket, charged EPE present)...");
  const U0 = ffL0.compute(ffL0.ref);
  const B0 = ffL0.bindingU;
  const U1 = ffL1.compute(ffL1.ref);
  const B1 = ffL1.bindingU;
  const dB = B1 - B0;
  console.log(`    L0 bindingU=${B0.toFixed(4)} total=${U0.toFixed(4)} :: ${describeBindingTerms(ffL0)}`);
  console.log(`    L1 bindingU=${B1.toFixed(4)} total=${U1.toFixed(4)} :: ${describeBindingTerms(ffL1)}`);
  console.log(`    ΔbindingU(L1−L0)=${dB.toFixed(4)} kcal/mol`);
  assert(Number.isFinite(U0) && Number.isFinite(B0), `L0 energy/bindingU finite`);
  assert(Number.isFinite(U1) && Number.isFinite(B1), `L1 energy/bindingU finite`);
  assert(Number.isFinite(dB) && Math.abs(dB) > 1e-6,
    `ΔbindingU finite + nonzero (wiring revives missing terms, Δ=${dB.toFixed(3)})`);
  assert(Math.abs(dB) < 10, `ΔbindingU sane magnitude |Δ|=${Math.abs(dB).toFixed(3)} < 10 kcal/mol`);
  assert(Math.abs(B0) < 50 && Math.abs(B1) < 50,
    `absolute bindingU sane (L0=${B0.toFixed(2)}, L1=${B1.toFixed(2)}, |U|<50)`);
  const f0 = maxAbsFinite(ffL0.forces);
  // recompute L1 forces (trackTerms below overwrites forces; capture first)
  ffL1.compute(ffL1.ref);
  const f1 = maxAbsFinite(ffL1.forces);
  assert(f0.allFinite && f1.allFinite, `forces all-finite (L0+L1)`);
  assert(f0.max < 10 && f1.max < 10,
    `forces sane max|F| (L0=${f0.max.toFixed(2)}, L1=${f1.max.toFixed(2)} < 10)`);

  // -----------------------------------------------------------------
  // [3] Alias wiring bit-identity (backward compat)
  // -----------------------------------------------------------------
  console.log("\n[3] wiring bit-identity (alias ≡ explicit, default untouched)...");
  const ffL1exp = new ForceField(sel, { rc: 10, gamma: 1.0, binding: { charges: true, hbMode: "directional" } }, ligs);
  const Ue = ffL1exp.compute(ffL1exp.ref);
  const Be = ffL1exp.bindingU;
  assert(Math.abs(Ue - U1) < 1e-12 && Math.abs(Be - B1) < 1e-12,
    `physicsLevel:"L1" ≡ explicit charges+directional (ΔU=${Math.abs(Ue - U1).toExponential(1)})`);
  const ffL0exp = new ForceField(sel, { rc: 10, gamma: 1.0, binding: { charges: false, hbMode: "off" } }, ligs);
  const U0e = ffL0exp.compute(ffL0exp.ref);
  assert(Math.abs(U0e - U0) < 1e-12,
    `default ≡ explicit L0 (bit-identical, ΔU=${Math.abs(U0e - U0).toExponential(1)})`);

  // -----------------------------------------------------------------
  // [4] Term exposure: Coulomb revived, directional wired
  // -----------------------------------------------------------------
  console.log("\n[4] per-term exposure (Coulomb revived in L1)...");
  const t0 = bindingTermsActive(ffL0);
  const t1 = bindingTermsActive(ffL1);
  assert(t0.coulomb === false && t0.directionalHB === false && t0.isotropicHB === true,
    `L0 terms: coulomb OFF / directional OFF / isotropic ON`);
  assert(t0.isotropicHBFallback === false && t0.hasInvalidFallback === false,
    `L0 fallback OFF (primary isotropic covers all pairs)`);
  assert(t1.coulomb === true && t1.directionalHB === true,
    `L1 terms: coulomb ON / directional ON (vSites=${ffL1._vSites ? ffL1._vSites.filter((s) => s.valid).length : 0} valid)`);
  assert(t1.isotropicHB === false && t1.isotropicHBFallback === true && t1.hasInvalidFallback === true,
    `Rev2/Issue2: L1 primary isotropic OFF but invalid-site fallback live (termini invalid)`);
  assert(describeBindingTerms(ffL1).includes("isoHB(fallback)"),
    `Rev2/Issue2: describeBindingTerms reports isoHB(fallback) in L1 (${describeBindingTerms(ffL1)})`);
  ffL1.trackTerms = true;
  ffL1.compute(ffL1.ref);
  assert(Number.isFinite(ffL1.bindCoulU) && Math.abs(ffL1.bindCoulU) > 1e-9,
    `L1 Coulomb accumulator nonzero (bindCoulU=${ffL1.bindCoulU.toFixed(3)} — salt-bridge path live)`);
  assert(Math.abs(ffL1.bindCoulU) < 10,
    `L1 Coulomb sane |coul|=${Math.abs(ffL1.bindCoulU).toFixed(3)} < 10 kcal/mol`);
  assert(Number.isFinite(ffL1.bindHBU), `L1 HB accumulator finite (${ffL1.bindHBU.toFixed(3)})`);

  // -----------------------------------------------------------------
  console.log("\n=================================================");
  console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
  console.log("=================================================");
  if (failed > 0) process.exit(1);
  console.log("PASS: test_l0_default_exposure.js — Rev1/Issue3 L0 exposure + L1 wiring validated");
}

main();
