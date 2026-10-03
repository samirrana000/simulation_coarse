/**
 * repro_cg_metal_params.mjs — REPRODUCTION for the CG metal-parameter gap.
 *
 * Builds a minimal PDB with a Zn inside a multi-atom HETATM group (plus a Na),
 * runs it through pdb.js parseLigands → cg/system.js → cg/params.js, and prints
 * the resolved sigma / eps / q / dG for every ligand atom on BOTH engines.
 *
 * Zero dependencies. Node-importable, no DOM globals.
 * Runnable: node tools/repro_cg_metal_params.mjs
 */

import { parseCa, selectSystem, parseLigands } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { HeavyForceField } from "../src/heavy.js";
import { METAL_ELEMENT, ELEMENT_LJ_DEFAULT } from "../src/physics/params.js";

/** A 3-residue poly-ALA Cα backbone + two HETATM groups: ZNF (Zn + 2 O) and NAI (Na + Cl). */
function hetatm(serial, name, resName, chain, resSeq, x, y, z, element) {
  return "HETATM" + String(serial).padStart(5) + " " + name.padEnd(4) +
    " " + resName.padStart(3) + " " + chain + String(resSeq).padStart(4) + "    " +
    x.toFixed(3).padStart(8) + y.toFixed(3).padStart(8) + z.toFixed(3).padStart(8) +
    "  1.00  0.00          " + element.padStart(2);
}
function atom(serial, name, resName, chain, resSeq, x, y, z, element) {
  return "ATOM  " + String(serial).padStart(5) + " " + name.padEnd(4) +
    " " + resName.padStart(3) + " A" + String(resSeq).padStart(4) + "    " +
    x.toFixed(3).padStart(8) + y.toFixed(3).padStart(8) + z.toFixed(3).padStart(8) +
    "  1.00  0.00           " + element.padStart(2);
}

// Three Cα beads, 3.8 Å apart, one per ALA (only the Cα is read by parseCa).
const CA = [
  [10.000, 10.000, 10.000],
  [13.800, 10.000, 10.000],
  [17.600, 10.000, 10.000],
];
const LINES = [
  "HEADER    SYNTHETIC ZINC-SITE REPRO",
  ...CA.map(([x, y, z], i) => atom(i + 1, " CA ", "ALA", "A", i + 1, x, y, z, "C")),
  // ZNF: a zinc coordinated by two oxygens (a zinc site as one HETATM group)
  hetatm(101, "ZN  ", "ZNF", "A", 101, 12.100, 10.500, 10.300, "ZN"),
  hetatm(102, "O1  ", "ZNF", "A", 101, 10.400, 10.900, 10.500, "O"),
  hetatm(103, "O2  ", "ZNF", "A", 101, 13.600, 10.900, 10.500, "O"),
  // NAI: a sodium + chloride pair
  hetatm(201, "NA  ", "NAI", "A", 201, 20.000, 12.000, 10.000, "NA"),
  hetatm(202, "CL  ", "NAI", "A", 201, 22.400, 12.000, 10.000, "CL"),
  "END",
];
const pdbText = LINES.join("\n");

const parsed = parseCa(pdbText);
const sel = selectSystem(parsed);
const ligands = parseLigands(pdbText);

console.log("=== REPRODUCTION: a Zn/Na in a hetero group, CG engine parameters ===");
console.log(`\nparseCa → ${sel.beads.length} Cα beads;  parseLigands → ${ligands.length} group(s) kept ` +
  `(rule: HETATM group with >= 2 heavy atoms):`);
for (const m of ligands) {
  console.log(`  {res:"${m.resName}", els:[${m.atoms.map((a) => `"${a.element}"`).join(",")}]}`);
}

const ff = new ForceField(sel, { rc: 10, gamma: 1.0 }, ligands);
console.log(`\nCG system: nProt=${ff.nProt} nLigAtoms=${ff.nLigAtoms} n=${ff.n}`);
console.log("\n metal | CG sigma / eps / q      / dG     | heavy sigma / eps / q   | sigma | eps    | dq");
const heavyAtoms = [];
for (const m of ligands) for (const a of m.atoms) heavyAtoms.push({
  element: a.element, atomName: "X1", resName: m.resName, x: a.x, y: a.y, z: a.z,
  isProtein: false, isLigand: true, isMetal: !!METAL_ELEMENT[a.element], isHetero: true,
});
const hff = new HeavyForceField({ atoms: heavyAtoms }, { gamma: 1.0 }, []);

let a = 0, metalRows = [];
for (const m of ligands) for (const at of m.atoms) {
  const el = at.element;
  const c = { sigma: ff._ligSigma[a], eps: ff._ligEps[a], q: ff._ligQ[a], dG: ff._ligdG[a] };
  const h = hff._elem[a];
  a++;
  if (!METAL_ELEMENT[el]) {
    console.log(` ${el.padEnd(5)} | ${`${c.sigma} / ${c.eps} / q=${c.q} / dG=${c.dG}`.padEnd(26)} | ` +
      `${`${h.sigma} / ${h.eps} / q=${h.q}`.padEnd(21)} | ${(c.sigma / h.sigma).toFixed(2)}x  | ` +
      `${(c.eps / h.eps).toFixed(2)}x  | ${(c.q - h.q).toFixed(1)}`);
    continue;
  }
  metalRows.push({ el, c, h });
  console.log(` ${el.padEnd(5)} | ${`${c.sigma} / ${c.eps} / q=${c.q} / dG=${c.dG}`.padEnd(26)} | ` +
    `${`${h.sigma} / ${h.eps} / q=${h.q}`.padEnd(21)} | ${(c.sigma / h.sigma).toFixed(2)}x  | ` +
    `${(c.eps / h.eps).toFixed(2)}x  | ${(c.q - h.q).toFixed(1)}  <== METAL`);
}

const neutralIon = metalRows.filter((r) => r.c.q === 0);
console.log(`\nVERDICT on the CG side: ${metalRows.length} metal atom(s) reached as ligand atoms.`);
console.log(`  • exactly ELEMENT_LJ_DEFAULT (σ ${ELEMENT_LJ_DEFAULT.sigma} / ε ${ELEMENT_LJ_DEFAULT.eps} ` +
  `/ q ${ELEMENT_LJ_DEFAULT.q} / dG ${ELEMENT_LJ_DEFAULT.dG}): ${neutralIon.length}/${metalRows.length}`);
console.log(`  • carrying their formal ion charge: ${metalRows.length - neutralIon.length}/${metalRows.length}`);
console.log(`  • coordination restraint terms in the CG system: ` +
  `${ff.holoSprings.length / 3} holo contact spring(s), 0 metal–donor springs ` +
  `(src/cg/* builds no metal coordination; only src/heavy/topology.js does)`);
console.log(`\nWorst sigma inflation: ` +
  `${Math.max(...metalRows.map((r) => r.c.sigma / r.h.sigma)).toFixed(2)}x`);
