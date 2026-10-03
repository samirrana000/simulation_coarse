/**
 * heavy.js — All-atom "heavy mode" force field with biophysics & high
 * performance. FACADE.
 *
 * Capabilities:
 *  - Explicit heavy atoms for protein, ligands, cofactors, and metal ions
 *  - Geometric covalent bond detection + metal coordination springs (N/O/S donors)
 *  - Fast O(N) spatial hashing grid for non-bonded interactions (LJ + GB + electrostatics)
 *  - Fast analytic proper and improper dihedral gradients (Blondel-Karplus / Bekker vector formulas)
 *  - Generalized Born + Debye-Hückel implicit solvent & AMBER partial charges
 *  - Hydrophobic SASA burial & directional H-bonding
 *  - Full compatibility with MOL2 ligand files, library placement, and PDB HETATM
 *
 * Handshake with CG mode: see docs/CG_HEAVY.md — CG (Cα ENM) and heavy
 * share the same reference coordinates (Å) and KCONV/KB_KCAL unit contract
 * (src/units.js:418.4). Switching is `selectSystem` vs `selectHeavy` +
 * `ForceField` vs `HeavyForceField`; viewer center/radius is protein-only
 * in both modes (src/viewer.js:163-182) so a distant ligand does not
 * inflate the camera (docs/CG_HEAVY.md §Handshake center/radius fix).
 *
 * THE SPLIT (2026-10) — WHY THIS FILE IS NOW A TABLE OF CONTENTS
 * ------------------------------------------------------------------
 * This file was 1674 LOC, the largest module in the repo and the only live
 * gradient in the automated gate. It held nine unrelated responsibilities, and
 * because it is the root of the heavy engine's import graph, one unparseable
 * line in any of them — a parsing loop, a ring DFS, a force kernel — took down
 * every heavy consumer with nothing but a SyntaxError in a browser console.
 * It is now eleven modules under src/heavy/, one per responsibility:
 *
 *   src/heavy/parse.js        PDB text -> heavy atoms + hetero groups    (171)
 *   src/heavy/select.js      chain/residue selection + ligand append     (90)
 *   src/heavy/topology.js    bond rows, rings, angles, propers,
 *                            metal coordination, flat typed lists       (308)
 *   src/heavy/params.js      cutoffs, C2 switch, K_ELEC, METAL_K,
 *                            heavyMass                                    (107)
 *   src/heavy/level.js       L0/L1/L2 resolver + describePhysics         (76)
 *   src/heavy/kernels.js     harmonicFlat / angleFlat (+ per-entry k)    (134)
 *   src/heavy/nonbonded.js   the three spatial-grid kernels              (240)
 *   src/heavy/weak.js        pi-stack / cation-pi / halogen + metal upgrade(162)
 *   src/heavy/energy.js      compute(): the term order                   (202)
 *   src/heavy/observables.js kineticTemp / rmsd / rmsdLig / rmsdAll      (64)
 *   src/heavy/forcefield.js  the HeavyForceField class + prototype wiring (385)
 *   src/heavy/pairs.js       the exclusion key                             (13)
 *
 * `src/heavy/` rather than `src/physics/heavy/`, because everything under
 * src/physics/ is a primitive BOTH engines share (observables, gb, sasa,
 * hbond, weakint, solvation, forcefield tables) and nothing here is shared
 * with CG mode — parseHeavy has no CG counterpart. Putting heavy-only input
 * parsing and a heavy-only force field inside the shared-physics directory
 * would break that distinction. src/heavy/ also mirrors the src/controllers/
 * split main.js got, so the repo now has one convention: a top-level engine
 * file plus a directory of its responsibilities.
 *
 * COMPATIBILITY
 * -------------
 * This module re-exports the same 21 names it always exported, so every
 * existing `import { … } from "./heavy.js"` keeps working untouched — no
 * consumer was repointed. Four of the validators in this repo read src/heavy.js
 * as TEXT (tests/test_rev3_issue1_heavy_physics.js, tests/test_element_params.js,
 * scripts/validate_binding_physics_r1.mjs, tests/test_doc_citations.js); those
 * were repointed at the heavy module family and the facts they pin are
 * unchanged. tests/test_module_size.js fails if any src/ module regrows past
 * 600 LOC, so this cannot silently happen again.
 *
 * Zero dependencies. No build step. Zero DOM globals.
 */

export {
  K_ELEC, SCREEN_LEN, R_CUT, R_SWITCH_ON, switchFunc, switchDeriv, METAL_K, heavyMass,
} from "./heavy/params.js";

export { countWarnings, parseHeavy } from "./heavy/parse.js";

export { selectHeavy, appendHeavyLigands } from "./heavy/select.js";

export {
  buildTopology, buildTopologyChunked, buildMetalCoordination, buildLists, angleAt,
} from "./heavy/topology.js";

export {
  DEFAULT_PHYSICS_LEVEL, HEAVY_PHYSICS_LEVELS, resolveHeavyPhysicsLevel,
} from "./heavy/level.js";

export { HeavyForceField } from "./heavy/forcefield.js";