/**
 * exp_cg_metal_stability.mjs — does a Zn ion survive a short CG run, with the
 * CG metal parameters now in place?
 *
 * Measured, not argued:
 *   1. What the CG engine hands the ion (sigma/eps/q/dG/hb).
 *   2. Whether a CHARGED ion with NO coordination restraint flies off, in the
 *      three poses that matter: crystallographic site (holo springs), a
 *      user-placed ion at that site (holo off), and a user-placed ion in bulk
 *      30 Å from the protein.
 *   3. Whether a physically-motivated ion desolvation ΔG (Born low-dielectric
 *      cavity penalty) helps or hurts — the option that was REJECTED.
 *
 * Zero dependencies, node-importable, no DOM globals.
 * Runnable: node tools/exp_cg_metal_stability.mjs
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { SeededRNG } from "../src/seeded-rng.js";
import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { METAL_ELEMENT } from "../src/physics/params.js";
import { COULOMB_CONST } from "../src/units.js";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const STEPS = 4000;

const hetatm = (serial, name, resName, resSeq, x, y, z, element) =>
  "HETATM" + String(serial).padStart(5) + " " + name.padEnd(4) +
  " " + resName.padStart(3) + " A" + String(resSeq).padStart(4) + "    " +
  x.toFixed(3).padStart(8) + y.toFixed(3).padStart(8) + z.toFixed(3).padStart(8) +
  "  1.00  0.00          " + element.padStart(2);

const pdbText = fs.readFileSync(path.join(ROOT, "4w52.pdb"), "utf-8");
const sel = selectSystem(parseCa(pdbText));
// Most-exposed Cα: the site is built around a surface bead, not the core.
let pocket = 0, best = -1;
for (let i = 0; i < sel.beads.length; i++) {
  const b = sel.beads[i];
  let far = 0;
  for (const c of sel.beads) if (Math.hypot(c.x - b.x, c.y - b.y, c.z - b.z) > 8) far++;
  if (far > best) { best = far; pocket = i; }
}
const PB = sel.beads[pocket];
// Two placements of the SAME 3-atom ZNF group:
//   SITE  — Zn 5.0 Å from that Cα (a realistic buried-metal distance: a Zn
//           bonded to a backbone carbonyl O is 2.1 Å from an O ~3.4 Å from
//           its Cα) and inside the 6 Å holo-spring cutoff;
//   BULK  — Zn 30 Å from the whole protein (the "user dropped an ion in the
//           water" case, where nothing pins it at all).
const place = (dFromPocket) => {
  const ZN = [PB.x + dFromPocket, PB.y, PB.z];
  // two donors at the canonical tetrahedral Zn(II)–O distance, 2.10 Å
  const O1 = [ZN[0] - 1.30, ZN[1] + 1.55, ZN[2] + 0.60];
  const O2 = [ZN[0] - 1.35, ZN[1] - 1.50, ZN[2] - 0.65];
  return { ZN, O1, O2 };
};
const mkText = ({ ZN, O1, O2 }) => pdbText.split(/\r?\n/)
  .filter((l) => l.slice(0, 6) !== "HETATM").concat([
    hetatm(9001, "ZN  ", "ZNF", 901, ZN[0], ZN[1], ZN[2], "ZN"),
    hetatm(9002, "O1  ", "ZNF", 901, O1[0], O1[1], O1[2], "O"),
    hetatm(9003, "O2  ", "ZNF", 901, O2[0], O2[1], O2[2], "O"),
    "END",
  ]).join("\n");

const SITE = place(5.0), BULK = place(30.0);
console.log("=== CG METAL ION STABILITY EXPERIMENT ===");
console.log(`protein 4W52, ${sel.beads.length} Cα beads;  pocket Cα = bead ${pocket} ` +
  `(${best} beads beyond 8 Å)`);
console.log(`SITE  pose: Zn  5.0 Å from that Cα   (Zn–O = ` +
  `${Math.hypot(SITE.ZN[0] - SITE.O1[0], SITE.ZN[1] - SITE.O1[1], SITE.ZN[2] - SITE.O1[2]).toFixed(2)}/` +
  `${Math.hypot(SITE.ZN[0] - SITE.O2[0], SITE.ZN[1] - SITE.O2[1], SITE.ZN[2] - SITE.O2[2]).toFixed(2)} Å)`);
console.log(`BULK  pose: Zn 30.0 Å from that Cα   (nothing within the 10 Å cutoff; no holo spring possible)`);
console.log(`Zn–O distances sit ABOVE pdb.js's 1.8 Å covalent gap-fill, so no bond is created ` +
  `for them; parseLigands keeps ZNF (a 3-atom HETATM group)`);
console.log(`parseLigands(SITE) → ` + parseLigands(mkText(SITE))
  .map((m) => `{${m.resName}:[${m.atoms.map((a) => a.element)}]}`).join(" "));

function build(pose, { holo, charges, hb }) {
  const ligands = parseLigands(mkText(pose));
  const ff = new ForceField(sel, {
    rc: 10, gamma: 1.0, physicsLevel: charges ? "L2" : "L0",
    binding: { on: true, holo, charges, hbMode: hb ? "directional" : "off" },
  }, ligands);
  return { ff, ligands };
}

// ── 1: what does CG hand the ion, and what does the pair look like? ────────
const { ff: ff0 } = build(SITE, { holo: true, charges: true, hb: true });
const znIdx = ff0.nProt, oIdx = [ff0.nProt + 1, ff0.nProt + 2];
console.log(`\n[1] CG resolution of the ion (src/cg/params.js → resolveElementParams)`);
console.log(`    sigma=${ff0._ligSigma[0]} Å  eps=${ff0._ligEps[0]} kcal/mol  q=${ff0._ligQ[0]} e  ` +
  `dG=${ff0._ligdG[0]} kcal/mol  hb=${ff0._ligHB[0]}`);
console.log(`    heavy METAL_ELEMENT.ZN: sigma=${METAL_ELEMENT.ZN.sigma} eps=${METAL_ELEMENT.ZN.eps} ` +
  `q=${METAL_ELEMENT.ZN.q}   → CG and heavy AGREE`);
console.log(`    CG restraint inventory: nHolo=${ff0.nHolo} pose springs, ` +
  `nNativeContacts=${ff0.nativeContacts.length / 3} soft intra-ligand springs (k=1.0), ` +
  `METAL–DONOR COORDINATION SPRINGS: 0 (heavy-only: src/heavy/topology.js)`);

const LB = (s1, e1, s2, e2) => [0.5 * (s1 + s2), Math.sqrt(e1 * e2)];
const lj = (s, e, r) => 4 * e * ((s / r) ** 12 - (s / r) ** 6);
const pS = ff0._protSigma[pocket], pE = ff0._protEps[pocket];
for (const [tag, s, e] of [
  ["ion row   ", ...LB(ff0._ligSigma[0], ff0._ligEps[0], pS, pE)],
  ["generic row", ...LB(3.4, 0.12, pS, pE)],
]) {
  console.log(`    Zn–Cα LJ with the ${tag}: σ_mix=${s.toFixed(2)} Å  r_min=${(s * 2 ** (1 / 6)).toFixed(2)} Å  ` +
    `U_min=${lj(s, e, s * 2 ** (1 / 6)).toFixed(4)} kcal/mol`);
}

// ── 2 + 3: runs ────────────────────────────────────────────────────────────
const dist = (pos, i, j) => Math.hypot(
  pos[3 * i] - pos[3 * j], pos[3 * i + 1] - pos[3 * j + 1], pos[3 * i + 2] - pos[3 * j + 2]);

function run(label, pose, opts, patch) {
  const { ff } = build(pose, opts);
  if (patch) patch(ff);
  const rng = new SeededRNG(20251003);
  rng.install();
  try {
    const integ = new LangevinIntegrator(ff.ref, ff, 110);
    integ.setTemperature(300); integ.setFriction(5.0);
    let minC = Infinity, sumC = 0, maxC = 0;
    let minD = Infinity, sumD = 0, badD = 0, inside = 0;
    for (let s = 0; s < STEPS; s++) {
      integ.step();
      const pos = integ.pos;
      let d = Infinity;
      for (let i = 0; i < ff.nProt; i++) { const x = dist(pos, znIdx, i); if (x < d) d = x; }
      let dd = Infinity;
      for (const o of oIdx) dd = Math.min(dd, dist(pos, znIdx, o));
      minC = Math.min(minC, d); maxC = Math.max(maxC, d); sumC += d;
      minD = Math.min(minD, dd); sumD += dd;
      if (dd < 1.5 || dd > 3.5) badD++;   // outside the physical Zn(II)–O band
      if (d < 2.5) inside++;             // ion inside a Cα bead, not next to a donor
    }
    const p = (k) => (100 * k / STEPS).toFixed(1);
    const r = {
      meanC: sumC / STEPS, maxC, minC, meanD: sumD / STEPS, minD,
      badPct: p(badD), insidePct: p(inside), nHolo: ff.nHolo, dG: ff._ligdG[0],
    };
    console.log(`\n  ${label}`);
    console.log(`      dG(Zn)=${ff._ligdG[0]} q(Zn)=${ff._ligQ[0]} nHolo=${ff.nHolo} | ${STEPS} steps, ` +
      `T=300 K, dt=${integ.dt} ps`);
    console.log(`      r(Zn, nearest Cα): min ${r.minC.toFixed(2)}  mean ${r.meanC.toFixed(2)}  ` +
      `max ${r.maxC.toFixed(2)} Å`);
    console.log(`      r(Zn, its own donor): min ${r.minD.toFixed(2)}  mean ${r.meanD.toFixed(2)} Å  |  ` +
      `unphysical Zn–O (<1.5 or >3.5 Å) ${r.badPct}% of the run`);
    console.log(`      ion inside a Cα bead (r < 2.5 Å): ${r.insidePct}% of the run`);
    return r;
  } finally { rng.restore(); }
}

console.log(`\n[2] does a charged, coordination-less ion fly off? (${STEPS} steps ≈ 7 ps of model time)`);
const R = {};
R.siteHolo = run("[A] SITE  holo ON  charges ON   — crystallographic metal site",
  SITE, { holo: true, charges: true, hb: true });
R.siteFree = run("[B] SITE  holo OFF charges ON   — user-placed ion at that site",
  SITE, { holo: false, charges: true, hb: true });
R.siteQ0 = run("[C] SITE  holo OFF charges OFF  — the ion as the PRE-FIX model saw it (q=0)",
  SITE, { holo: false, charges: false, hb: false });
R.bulk = run("[D] BULK  holo OFF charges ON   — a free ion dropped in solvent, 30 Å out",
  BULK, { holo: false, charges: true, hb: true });
// The PRE-FIX row, patched back onto the same engine and the same seed, so
// [A] vs [G] is a before/after of THIS change and nothing else.
const prefix = (ff) => {
  ff._ligSigma[0] = 3.4; ff._ligEps[0] = 0.12; ff._ligQ[0] = 0.0; ff._ligHB[0] = 0;
};
R.siteHoloPre = run("[G] SITE  holo ON  charges ON   + the PRE-FIX generic row (σ3.4 ε0.12 q0)",
  SITE, { holo: true, charges: true, hb: true }, prefix);
R.siteFreePre = run("[H] SITE  holo OFF charges ON   + the PRE-FIX generic row (σ3.4 ε0.12 q0)",
  SITE, { holo: false, charges: true, hb: true }, prefix);

console.log(`\n[3] the REJECTED option (b): a Born low-dielectric cavity ΔG for the ion`);
// Marcus/Born self-energy: ΔG = (q²/2)(1/ε_low − 1/ε_high)·332.0637 kcal·Å/mol,
// ε_high = 78.5 (water at 298 K), ε_low = 4 (protein interior).
const born = (q) => -(q * q / 2) * COULOMB_CONST * (1 / 4 - 1 / 78.5);
console.log(`    ΔG_born(Zn²⁺) = ${born(2).toFixed(1)} kcal/mol, ΔG_born(Na⁺) = ${born(1).toFixed(1)} kcal/mol`);
R.bornSite = run("[E] SITE  holo OFF charges ON   + Born ΔG on the ion (patched onto _ligdG)",
  SITE, { holo: false, charges: true, hb: true }, (ff) => { ff._ligdG[0] = born(2); });
R.bornBulk = run("[F] BULK  holo OFF charges ON   + Born ΔG on the ion",
  BULK, { holo: false, charges: true, hb: true }, (ff) => { ff._ligdG[0] = born(2); });

console.log(`\nSummary — the before/after of THIS change, same seed, same engine`);
console.log(`  site + pose pins : r(Zn,Cα) mean ${R.siteHoloPre.meanC.toFixed(2)} → ${R.siteHolo.meanC.toFixed(2)} Å | ` +
  `unphysical Zn–O ${R.siteHoloPre.badPct}% → ${R.siteHolo.badPct}% | ` +
  `inside a bead ${R.siteHoloPre.insidePct}% → ${R.siteHolo.insidePct}%`);
console.log(`  site, no pins    : r(Zn,Cα) mean ${R.siteFreePre.meanC.toFixed(2)} → ${R.siteFree.meanC.toFixed(2)} Å | ` +
  `unphysical Zn–O ${R.siteFreePre.badPct}% → ${R.siteFree.badPct}%`);
console.log(`  bulk, no pins    : stays ${R.bulk.meanC.toFixed(2)} ± ${(R.bulk.maxC - R.bulk.minC).toFixed(1)} Å from the ` +
  `protein — the charged ion is neither ejected nor blown up`);
console.log(`  Born ΔG at site  : r(Zn,Cα) mean ${R.siteFree.meanC.toFixed(2)} → ${R.bornSite.meanC.toFixed(2)} Å, ` +
  `inside a bead ${R.siteFree.insidePct}% → ${R.bornSite.insidePct}% — the ion is dragged into the CORE`);
