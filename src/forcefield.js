/**
 * forcefield.js — Cα coarse-grained potential with anisotropic-network-style
 * native contacts and implicit-solvent excluded volume.
 *
 * THIS IS A FACADE. The implementation lives in src/cg/, split by
 * responsibility. This file keeps the historical import contract intact
 * (`import { ForceField, resolvePhysicsLevel, KB_KCAL, KCONV } from
 * "./forcefield.js"`) for every existing consumer — src/analysis-panel.js,
 * src/controllers/*, 25 tests, 4 bench scripts, cli.js and 7 scripts/*.mjs —
 * none of which had to change. Adding a new CG responsibility therefore
 * touches one file under src/cg/ instead of the force field's single point of
 * syntactic failure; a parse error in one term's kernel can no longer take
 * down every importer of the whole CG engine.
 *
 * The layout under src/cg/ mirrors the src/heavy/ convention established by
 * the heavy-atom split — a directory named for the engine, one module per
 * responsibility, and this file as its table of contents:
 *
 *   src/cg/level.js      L0/L1/L2 tier resolution + describePhysics
 *   src/cg/params.js     constructor scalars; per-residue and per-element tables
 *   src/cg/system.js     particle set: ref coords, appended ligand atoms, masses
 *   src/cg/topology.js   connectivity, ENM contacts, exclusion ledger, native contacts
 *   src/cg/springs.js    spring-network mutators (holo, ML prior, seq weights, Tirion)
 *   src/cg/grid.js       spatial-hash codec + every preallocated compute buffer
 *   src/cg/compute.js    the compute() term assembly (the summation order lives here)
 *   src/cg/forcefield.js the ForceField class: build order + thin wrappers
 *
 * Kernels that were already separate are NOT here and were not moved:
 * ff-harmonic.js (pairs/angles/springs), ff-repulsion.js (grid excluded
 * volume), ff-binding.js (cross LJ + electrostatics + H-bond + EEF1-lite),
 * ff-params.js (parameter tables), physics/observables.js (kineticTemp/rmsd),
 * physics/forcefield/tirion_anm.js (Tirion γ + SS basins), physics/virtual-sites.js
 * (directional-HB O sites). src/funnel.js is attached by the caller via
 * setFunnel() and is never imported here.
 *
 * Energy model (all standard CG/ENM physics; units Å, ps, kcal/mol, Da):
 *
 *   1. PEPTIDE BONDS (backbone connectivity, per contiguous segment):
 *        U_bond = Σ ½ k_b (r − r0)² ,       k_b ≈ 100 kcal/mol/Å²
 *      with r0 = equilibrium Cα–Cα distance (≈ 3.81 Å) from the input structure.
 *
 *   2. BACKBONE ANGLES (pseudo-angle i−1, i, i+1 within a segment):
 *        U_ang = Σ ½ k_θ (θ − θ0)² ,        k_θ ≈ 20 kcal/mol/rad²
 *      θ0 taken from the input coordinates (native-like local stiffness).
 *
 *   3. NATIVE CONTACTS — Elastic Network Model (K Tirion-style / ANM-uniform
 *      spring): every bead pair whose NATIVE distance r0 ≤ Rc and that is not
 *      already a 1-2 or 1-3 sequence neighbour gets a Hookean spring:
 *        U_ENM = Σ ½ γ (r − r0)² ,          H(Rc − r0) implicit by construction
 *      This preserves the native fold while allowing thermal fluctuations that
 *      reproduce experimental B-factor statistics at T = 300 K for γ ≈ 1.
 *
 *   4. NON-NATIVE EXCLUDED VOLUME (implicit-solvent, repulsive 12-10-free
 *      Lennard-Jones form, attractive tail switched off so beads never
 *      collapse spuriously):
 *        U_rep = ε [ (r_e/r)¹² − 2 (r_e/r)⁶ + 1 ]   for r < r_e,
 *              = 0                                      otherwise,
 *      with  r_e = 2^(1/6) σ ,  σ = 4.0 Å (Cα bead diameter), ε = 0.3 kcal/mol.
 *      This term is smooth at the cutoff and prevents chain crossing.
 *
 *   5. PROTEIN–LIGAND BINDING (pair pass): protein beads 0..nProt−1 vs. ligand
 *      atoms nProt..n−1 are excluded from the grid repulsion above and handled
 *      by a dedicated binding term — cross 12-6 LJ (attractive well, Lorentz–
 *      Berthelot combining rules), screened electrostatics (εr = 4+76·tanh(r/8))
 *      and an H-bond term, each smoothstep-switched off at bindRcut ≈ 9 Å.
 *
 *   5b. HOLO CONTACT SPRINGS (native protein–ligand pose): every protein
 *      bead–ligand atom pair whose NATIVE distance r0 ≤ 6 Å gets a Hookean
 *      spring U = ½ γ_lig (r − r0)² with γ_lig = 0.5 kcal/mol/Å², encoding the
 *      observed bound pose of a holo complex (e.g. T4 lysozyme L99A + benzene,
 *      PDB 4W52). The pair is excluded from the binding pair pass above so the
 *      spring alone governs the native-contact distance. Gated by
 *      par.binding.holo (default true); no-op without ligands.
 *
 *   6. EEF1-LITE DESOLVATION/BURIAL (apolar binding): while the same 27-cell
 *      scan drives the pair terms above, each ligand atom also accumulates a
 *      "soft protein-occupancy count" n_a = Σ_j g(r_aj) over the protein beads
 *      it sees, with g(r) = exp(−(r−r0)²/(2σ²)), r0 = 4.5 Å, σ = 1.8 Å. g has
 *      decayed to ≈ 0.04 by r ≈ 9 Å = bindRcut, so the existing grid pass fully
 *      covers the density field — no second grid, no second cutoff. The burial
 *      fraction B_a = 1 − exp(−n_a/3) ∈ [0,1] (0 exposed, 1 buried) converts the
 *      per-atom hydrophobic transfer energy ΔG_a (≈ −0.55 kcal/mol for C) into
 *      a burial cost U_desolv = Σ_a ΔG_a·B_a (negative ⇒ burial favorable).
 *      Forces follow the chain rule dU/dr = ΔG_a·(dB/dn)·(dg/dr) and are
 *      applied in a second pass over the recorded (a, j, r) pairs (see
 *      `_binding`).
 *
 *   7. OPTIONAL FUNNEL BIAS (off by default; par.funnel.on): external
 *      collective-variable bias on the ligand↔pocket COM distance — a binding
 *      funnel that is flat inside the bound state plus well-tempered
 *      metadynamics for PMF reconstruction (see src/funnel.js). Attached by
 *      the caller via setFunnel(); contributes nothing while funnelOn is false.
 *
 * The system is a unified multi-particle set: the Cα beads (global indices
 * 0..nProt−1) followed by any ligand heavy atoms (indices nProt..n−1). When no
 * ligands are given the protein-only model is reproduced exactly. Ligands get
 * a united-atom internal force field (bonds/angles/impropers from ligand.js)
 * plus the same excluded-volume grid; covalently linked ligand pairs are
 * excluded from the repulsion so the rigid molecule never repels itself.
 *
 * Performance: bonds/angles/ENM springs are O(N); the repulsive term uses a
 * uniform grid (cell list) in 3D rebuilt each force call — O(N) average.
 * All arithmetic is done on flat Float64Array buffers for speed; forces are
 * accumulated in a preallocated flat array — zero allocations per step.
 */

// Canonical units live in units.js; re-export keeps backward compat for
// integrator.js / funnel.js / tests that historically imported from here. The
// import below now exists ONLY to feed this re-export — the physics moved to
// physics/observables.js.
export { KB_KCAL, KCONV } from "./units.js";

// The class. src/cg/forcefield.js owns it; the tier resolver
// (DEFAULT_PHYSICS_LEVEL / CG_PHYSICS_LEVELS / resolvePhysicsLevel) lives in
// src/cg/level.js and is re-exported here unchanged — its export shape is part
// of the public contract tests/test_l0_default_exposure.js and
// tests/test_rev2_issue1_physics_level.js read. This export list is exactly
// the pre-split one; describePhysics() is a class method and stays one.
export { ForceField } from "./cg/forcefield.js";
export {
  DEFAULT_PHYSICS_LEVEL, CG_PHYSICS_LEVELS, resolvePhysicsLevel,
} from "./cg/level.js";