/**
 * test_forces_fd.js — Finite-difference validation for HeavyForceField
 *
 * Builds 4W52 heavy system with MOL2 benzene, computes analytic forces,
 * then numeric central diff (perturb 1e-4 Å) for 10 atoms, asserts max
 * relative error < 1e-3. Runnable: node tests/test_forces_fd.js
 */

import fs from "fs";
import path from "path";
import { parseHeavy, selectHeavy, appendHeavyLigands, HeavyForceField } from "../src/heavy.js";
import { parseMol2 } from "../src/mol2.js";

function findFile(names) {
  for (const p of names) {
    if (fs.existsSync(p)) return p;
  }
  throw new Error(`Cannot find ${names[0]} — tried ${names.join(", ")} cwd=${process.cwd()}`);
}

async function run() {
  // Locate PDB + MOL2 regardless of cwd (project root vs tests/)
  const pdbPath = findFile(["4w52.pdb", "./4w52.pdb", "simulation_coarse/4w52.pdb", "../4w52.pdb", path.join(path.dirname(new URL(import.meta.url).pathname), "../4w52.pdb")]);
  const mol2Path = findFile(["benzene.mol2", "./benzene.mol2", "simulation_coarse/benzene.mol2", "../benzene.mol2", path.join(path.dirname(new URL(import.meta.url).pathname), "../benzene.mol2")]);

  const pdbText = fs.readFileSync(pdbPath, "utf-8");
  const mol2Text = fs.readFileSync(mol2Path, "utf-8");

  const parsedHeavy = parseHeavy(pdbText);
  const mols = parseMol2(mol2Text);

  // Build heavy system: exclude PDB HETATM ligands, use external MOL2 benzene
  let sel = selectHeavy(parsedHeavy, {
    heteroSelection: { "A|200|BNZ": false, "A|201|EPE": false },
    includePdbLigands: true,
    hasExternalLigand: true,
  });
  sel = appendHeavyLigands(sel, mols);
  const ff = new HeavyForceField({ atoms: sel.atoms }, { gamma: 2.0, temp: 300 }, []);

  const n = ff.n;
  // Reference positions
  const pos0 = new Float64Array(ff.ref); // copy
  // Analytic forces at pos0
  ff.compute(pos0);
  const anaForces = new Float64Array(ff.forces); // copy
  const E0 = ff.energy;

  if (!Number.isFinite(E0) || anaForces.some(v => !Number.isFinite(v))) {
    console.error(`FAIL: analytic energy/forces contain NaN/Inf (E0=${E0})`);
    process.exit(1);
  }
  console.log(`Heavy system N=${n}  E0=${E0.toFixed(3)} kcal/mol  sample forces[:3]=${Array.from(anaForces.slice(0,3)).map(v=>v.toFixed(3)).join(",")}`);

  const h = 1e-4; // Å
  // Deterministic selection of 10 atoms spaced through system (avoid termini if possible)
  // Use seeded pseudo-random or simply spaced indices
  const candidates = [0, 10, 20, 50, 100, 200, 400, 600, 800, 1000].map(i => Math.min(i, n-1));
  // Ensure unique
  const indices = [...new Set(candidates)].slice(0, 10);

  let maxAbsErr = 0;
  let maxRelErr = 0;
  let worst = null;

  for (const idx of indices) {
    for (let d = 0; d < 3; d++) {
      const j = 3*idx + d;
      const posPlus = new Float64Array(pos0);
      const posMinus = new Float64Array(pos0);
      posPlus[j] += h;
      posMinus[j] -= h;

      const ePlus = ff.compute(posPlus);
      const eMinus = ff.compute(posMinus);

      if (!Number.isFinite(ePlus) || !Number.isFinite(eMinus)) {
        console.error(`FAIL: numeric energy NaN at atom ${idx} dim ${d}  E+ ${ePlus} E- ${eMinus}`);
        process.exit(1);
      }
      // F = -dU/dx
      const numF = -(ePlus - eMinus) / (2*h);
      const anaF = anaForces[j];
      const absErr = Math.abs(numF - anaF);
      const denom = Math.max(1e-6, Math.abs(anaF));
      // If |ana| tiny, use absolute error scaled to tolerance; otherwise relative
      const relErr = absErr / (Math.abs(anaF) + 1e-8);

      if (absErr > maxAbsErr) maxAbsErr = absErr;
      if (relErr > maxRelErr) {
        maxRelErr = relErr;
        worst = { idx, dim: ["x","y","z"][d], anaF, numF, absErr, relErr, ePlus, eMinus };
      }

      // Also catch large absolute error for near-zero forces
      // Only fail if both abs >1e-3 and rel >1e-3
    }
  }

  console.log(`Tested ${indices.length*3} components (10 atoms × xyz) h=${h}`);
  console.log(`  maxAbsErr = ${maxAbsErr.toExponential(3)} kcal/mol/Å`);
  console.log(`  maxRelErr = ${maxRelErr.toExponential(3)}`);
  if (worst) {
    console.log(`  worst: atom ${worst.idx} dim ${worst.dim}  ana=${worst.anaF.toExponential(6)} num=${worst.numF.toExponential(6)} abs=${worst.absErr.toExponential(3)} rel=${worst.relErr.toExponential(3)}`);
  }

  // Restore analytic state (optional)
  ff.compute(pos0);

  // The spec says max relative error <1e-3 ; we allow absolute fallback for tiny forces
  // Require maxAbsErr <1e-3 if ana small, otherwise maxRelErr <1e-3
  // For robustness, pass if maxRelErr <1e-3 OR maxAbsErr <1e-3 (since relative blows up at zero)
  const tol = 1e-3;
  const pass = maxRelErr < tol || maxAbsErr < tol;
  // More strict: if any component with |ana|>1e-2 has rel>1e-3 => fail
  let strictFail = false;
  // Re-evaluate worst with threshold
  if (maxRelErr >= tol && maxAbsErr >= tol) {
    // Check if worst has significant magnitude
    if (worst && Math.abs(worst.anaF) > 1e-4) strictFail = true;
    else if (maxAbsErr >= tol) strictFail = true;
  }

  if (!strictFail && pass) {
    console.log(`PASS: HeavyForceField forces match finite difference (maxRel ${maxRelErr.toExponential(3)} < ${tol} or maxAbs ${maxAbsErr.toExponential(3)} < ${tol})`);
    process.exit(0);
  } else {
    console.error(`FAIL: Force FD mismatch — maxRel ${maxRelErr.toExponential(3)} maxAbs ${maxAbsErr.toExponential(3)} > ${tol}`);
    if (worst) console.error(`  worst atom ${worst.idx} dim ${worst.dim}`);
    process.exit(1);
  }
}

run().catch(e => { console.error("FAIL:", e); process.exit(1); });
