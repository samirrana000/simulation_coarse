/**
 * test_placement_hetero.js — rev1-issue2: placement must see hetero atoms.
 *
 * Regression: the heavy placement collision set was protein-only
 * (getProteinCoordsAndSigma used subarray 0..3*nProt; clashGrad /
 * relaxClash / placeLigand truncated to nProt), so hetero atoms, PDB
 * ligands/cofactors/metals and existing-ligand atoms were invisible and a
 * second placement could land inside them with false converged:true.
 *
 * This test places the same ligand at the same target with the same seed
 * against (a) a protein-only collision set and (b) a protein+hetero set
 * (hetero cage around the target). (a) must converge; (b) must NOT falsely
 * converge and must report a sharply higher residual. It also covers the
 * backward-compatible relaxClash(pos, mol, protein) signature and the
 * opts/protein .excludeFrom slot-exclusion policy.
 *
 * Run: node tests/test_placement_hetero.js
 */

import { placeLigand, relaxClash, findPocketCenter } from "../src/placement.js";
import { getProteinCoordsAndSigma } from "../src/ligand-panel.js";
// NOTE: ligand-panel.js binds `state` from "./ui.js?v=10"; the query string
// makes that a distinct module instance from "../src/ui.js", so the test
// must inject through the same ?v=10 instance the panel reads.
import { state } from "../src/ui.js?v=10";

let fails = 0;
let passes = 0;
function assert(cond, msg) {
  if (!cond) { console.error(`  ✗ FAIL: ${msg}`); fails++; process.exitCode = 1; }
  else { console.log(`  ✓ ${msg}`); passes++; }
}

function centroid(pos) {
  const n = pos.length / 3;
  let x = 0, y = 0, z = 0;
  for (let i = 0; i < n; i++) { x += pos[3 * i]; y += pos[3 * i + 1]; z += pos[3 * i + 2]; }
  return [x / n, y / n, z / n];
}

function dist(a, b) {
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

console.log("=== test_placement_hetero — hetero-aware collision set (rev1-issue2) ===");

// 3-atom ligand fragment (elements resolve via LIG_ELEMENT).
const mol = {
  atoms: [
    { element: "C", x: 0, y: 0, z: 0 },
    { element: "C", x: 1.4, y: 0, z: 0 },
    { element: "O", x: -1.2, y: 0.5, z: 0 },
  ],
};
const target = [0, 0, 0];
const SEED = 7;

// Protein beads far from the target: placement there is clash-free.
const protPos = new Float64Array([40, 0, 0, 40, 5, 0, 45, 0, 0, 40, 0, 5]);
const protSig = new Float64Array(4).fill(3.8);

// Hetero cage (e.g. cofactor/metal shell) around the target: every escape
// direction stays inside clash distance, so a hetero-aware relaxation
// cannot converge at the target.
const cageXYZ = [2.5, 0, 0, -2.5, 0, 0, 0, 2.5, 0, 0, -2.5, 0, 0, 0, 2.5, 0, 0, -2.5];
const hetPos = new Float64Array(cageXYZ);
const hetSig = new Float64Array(6).fill(3.4);

// (a) protein-only collision set — legacy behavior, must converge.
const only = placeLigand(mol, target, {
  protein: { pos: protPos, sigma: protSig },
  seed: SEED,
});
assert(only.converged === true, `protein-only placement converges (got converged=${only.converged})`);
assert(Number.isFinite(only.residualClash) && only.residualClash < 1.0,
  `protein-only residual clash-free (${only.residualClash.toFixed(3)} < 1.0)`);
assert(dist(centroid(only.pos), target) < 0.5,
  `protein-only pose lands on target (drift ${dist(centroid(only.pos), target).toFixed(3)} Å)`);

// (b) protein+hetero collision set — must NOT falsely converge.
const both = placeLigand(mol, target, {
  protein: {
    pos: Float64Array.from([...protPos, ...hetPos]),
    sigma: Float64Array.from([...protSig, ...hetSig]),
  },
  seed: SEED,
});
assert(both.converged === false,
  `protein+hetero placement refuses occupied site (converged=${both.converged})`);
assert(both.residualClash > only.residualClash + 1.0,
  `hetero residual sharply higher (${both.residualClash.toFixed(3)} vs ${only.residualClash.toFixed(3)})`);
assert(both.iterations > 0,
  `hetero-aware relaxation attempted escape (${both.iterations} iters)`);

// Backward-compatible relaxClash(pos, mol, protein) with protein-only arrays.
const nativePos = Float64Array.from([0, 0, 0, 1.4, 0, 0, -1.2, 0.5, 0]);
const legacy = relaxClash(nativePos, mol, { pos: protPos, sigma: protSig });
assert(legacy.converged === true && legacy.iterations === 0,
  `legacy protein-only relaxClash signature converges immediately (${legacy.residualClash.toFixed(3)})`);

// excludeFrom policy: a stale incoming-slot atom sitting exactly on the
// target distorts the pose unless excluded (ligand-panel truncates to
// ligandStart; relaxClash also honors opts/protein .excludeFrom).
const stalePos = Float64Array.from([...protPos, 0, 0, 0]);
const staleSig = Float64Array.from([...protSig, 3.4]);
const stale = { pos: stalePos, sigma: staleSig };
const noEx = relaxClash(Float64Array.from([0, 0, 0, 1.4, 0, 0, -1.2, 0.5, 0]), mol, stale);
const withEx = relaxClash(Float64Array.from([0, 0, 0, 1.4, 0, 0, -1.2, 0.5, 0]), mol, stale, { excludeFrom: 4 });
assert(withEx.converged === true && withEx.iterations === 0,
  `excludeFrom skips stale incoming slot (converged, 0 iters, residual ${withEx.residualClash.toFixed(3)})`);
assert(noEx.iterations > 0,
  `without exclusion the stale slot forces flee (${noEx.iterations} iters)`);
assert(noEx.residualClash > withEx.residualClash,
  `stale slot distorts pose (${noEx.residualClash.toFixed(3)} vs ${withEx.residualClash.toFixed(3)})`);

// findPocketCenter: default count covers the full array (explicit nProt unchanged).
const pocketExplicit = findPocketCenter(protPos, 4);
const pocketDefault = findPocketCenter(protPos);
assert(pocketExplicit.every(Number.isFinite) && pocketDefault.every(Number.isFinite) &&
  dist(pocketExplicit, pocketDefault) === 0,
  `findPocketCenter default count matches explicit nProt (${pocketDefault.map((v) => v.toFixed(2)).join(", ")})`);

// -----------------------------------------------------------------
// ligand-panel getProteinCoordsAndSigma: heavy-mode collision set must
// be protein+hetero (incoming slot excluded), sigma from ff._elem.
// Pre-fix this returned protein-only, so a placement at an occupied
// hetero site falsely reported converged:true.
// -----------------------------------------------------------------
{
  // Fake heavy system: 4 protein beads (far) + 6 hetero cage atoms at the
  // target (first hetero is an FE metal) + 3-atom stale incoming slot.
  const fullPos = Float64Array.from([
    ...protPos, ...hetPos, 100, 100, 100, 101, 100, 100, 100, 101, 100,
  ]);
  const elemSigma = [3.8, 3.8, 3.8, 3.8, 1.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4];
  state.ff = {
    nProt: 4, n: 13, ligandStart: 10,
    _elem: elemSigma.map((sigma) => ({ sigma })),
  };
  state.integ = { pos: fullPos };
  state.sel = { atoms: elemSigma.map(() => ({ element: "C" })) };

  const coll = getProteinCoordsAndSigma();
  assert(coll.pos.length === 30 && coll.sigma.length === 10,
    `heavy collision set is protein+hetero, incoming slot excluded (got ${coll.sigma.length} atoms)`);
  assert(coll.sigma[4] === 1.4 && coll.sigma[5] === 3.4,
    `hetero sigma sourced per-atom from ff._elem incl. metal (got [${coll.sigma[4]}, ${coll.sigma[5]}])`);

  const placedHet = placeLigand(mol, target, { protein: coll, seed: SEED });
  assert(placedHet.converged === false,
    `panel-built hetero set refuses occupied site (converged=${placedHet.converged})`);

  // Old caller semantics (protein-only slice) on the same system falsely converge.
  const legacySlice = {
    pos: fullPos.subarray(0, 3 * 4),
    sigma: Float64Array.from(elemSigma.slice(0, 4)),
  };
  const placedLegacy = placeLigand(mol, target, { protein: legacySlice, seed: SEED });
  assert(placedLegacy.converged === true,
    `protein-only slice falsely converges inside hetero cage (reproduces pre-fix bug)`);

  // CG mode (no ligandStart/_elem): unchanged protein-only legacy behavior.
  state.ff = { nProt: 4, n: 7, _protSigma: Float64Array.from([4.0, 4.1, 3.8, 3.6]) };
  state.integ = { pos: Float64Array.from([...protPos, 0, 0, 0, 1, 1, 1, 2, 2, 2]) };
  state.sel = null;
  const collCg = getProteinCoordsAndSigma();
  assert(collCg.sigma.length === 4 && collCg.sigma[1] === 4.1,
    `CG collision set stays protein-only via _protSigma (got ${collCg.sigma.length} atoms)`);

  state.ff = null; state.integ = null; state.sel = null;
}

if (fails === 0) console.log(`\n=== PASS test_placement_hetero (${passes} assertions) ===`);
else { console.error(`\n=== FAIL test_placement_hetero (${fails} failed, ${passes} passed) ===`); process.exit(1); }
