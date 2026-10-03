/**
 * ff-params.js — the STABLE FACADE for model constants and per-element /
 * per-residue-class parameter tables, plus the backbone + ENM constants that
 * only the coarse-grained engine uses.
 *
 * Split out of forcefield.js (item 5, "move never rewrite": every value below
 * is byte-identical to the original tables). forcefield.js re-exports
 * KB_KCAL and KCONV so downstream modules (integrator.js, funnel.js) keep
 * importing them from "./forcefield.js" unchanged.
 *
 * Since goal M7 this file OWNS almost no parameter literals: the per-element
 * LJ/charge tables, the metal-ion table, the covalent radii and the CG
 * residue-class table live once in src/physics/params.js and are re-exported
 * here under their historical names, so no downstream import moves. What is
 * still declared below is the CG-only model: the backbone Boltzmann-inversion
 * constants, the Bahar sequence weights, and the holo-pair compression floor.
 *
 * ── Backbone bond/angle constants (Boltzmann inversion) ──────────────────
 *   k_b = 100 kcal/mol/Å²  — Cα–Cα peptide bond stiffness
 *   k_θ = 20  kcal/mol/rad² — Cα pseudo-angle stiffness
 *
 * Source & rationale (cite Tirion 1996, PRL 77:1905; AMBER ff14SB Cα):
 *   • Tirion proposed ENM uniform harmonic springs to reproduce crystallographic
 *     B-factors: ½γ(r−r0)² with γ~1 kcal/mol/Å² for non-bonded contacts. For
 *     covalent peptide geometry (r0≈3.81 Å, θ0≈90–150°) the local curvature must
 *     be ~2 orders stronger to keep bond-length RMS <0.05 Å at 300 K. From
 *     Boltzmann inversion on AMBER all-atom Cα–Cα distances:
 *         k = k_B T / σ²  where σ² = ⟨(r−r0)²⟩ from explicit-solvent MD.
 *     AMBER ff14SB equilibrium Cα variance σ_b≈0.045 Å ⇒ k_b = RT/σ² ≈
 *     0.6/0.002 ≈ 120, rounded to 100 for Langevin stability (dt=4 fs, see
 *     integrator.js _pickDt). Similarly σ_θ≈0.14 rad (≈8°) ⇒
 *     k_θ = RT/σ_θ² ≈ 0.6/0.02 = 30, softened to 20 to allow loop flexibility
 *     while preserving helicity. Values validated in tests/test_bond_dist.js:
 *     200-step Langevin on 1crn keeps ⟨b⟩ = 3.81±0.05 Å and ⟨θ⟩ within 2° of
 *     native (see docs/CG_HEAVY.md and forcefield.js:107 comment).
 *
 *   References:
 *     [1] Tirion MM, PRL 77, 1905 (1996) — ENM single-parameter harmonic model.
 *     [2] Atilgan et al., Biophys J 80, 505 (2001) — ANM & B-factor correlation.
 *     [3] Bahar et al., Folding & Design 2, 173 (1997) — GNM sequence weighting.
 *     [4] Case et al., AMBER 2020 Manual, ff14SB Cα parameters (r0=3.81 Å).
 *
 * ── Sequence-dependent ENM weights (Bahar-style) ────────────────────────
 *   SEQ_WEIGHT[resClass] = dimensionless stiffness multiplier w_i per bead
 *   used to build springK_seq (see forcefield.js:234). The stub formula
 *       K_seq(i,j) = gamma * (1 + 0.2*(w_i + w_j)/2)
 *   makes hydrophobic contacts (H) reference 1.0, aromatic (A) slightly softer
 *   (0.9, delocalized π cloud), polar/charged (P/Cp/Cn) slightly stiffer
 *   (1.1/1.05 via H-bond capacity) — consistent with Miyazawa–Jernigan
 *   contact energies rescaled to ±10% (Bahar 1997). w_i retrieved as
 *   SEQ_WEIGHT[RES_CLASS_OF[bead.resName]]; w_j likewise. Uniform model is
 *   default; sequence-weighted is opt-in via ForceField.applySeqWeights().
 */

// Boltzmann constant and the kcal→mechanical conversion are RE-EXPORTED from
// src/units.js, which is the declared single unit contract (docs/UNITS.md,
// TRANSFORMATION_PLAN_100 A05). This file previously carried its own copies:
// KB_KCAL = 0.0019872041 here vs 0.001987204 in units.js — a 5.03e-8 relative
// split that made the reported temperature depend on which module a value was
// imported from (CG via units.js, heavy via ff-params.js). Re-exporting makes
// divergence structurally impossible rather than merely fixed once.
export { KB_KCAL, KCONV } from "./units.js";

// ── The PARAMETER contract is RE-EXPORTED from src/physics/params.js ─────
// This file used to OWN the per-element LJ / partial-charge tables, the metal
// table, the covalent radii and the CG residue-class table. Two of its
// consumers then grew private copies: src/heavy.js declared its own
// `HEAVY_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0 }` (dead, but one
// edit away from being a live sigma/eps/q split between the CG and heavy
// paths). Moving the literals into src/physics/params.js and re-exporting them
// here makes a second declaration a structural impossibility rather than a
// thing to remember — the same remedy, and the same facade shape, as the
// units.js re-exports above.
//
// The PUBLIC NAMES ARE UNCHANGED. Every downstream import keeps working:
//     LIG_ELEMENT        ← ELEMENT_LJ          (the old name was misleading:
//                           heavy.js uses this table for protein backbone and
//                           side-chain atoms too, not only for ligands)
//     LIG_ELEMENT_DEFAULT← ELEMENT_LJ_DEFAULT
//     METAL_ELEMENT, METAL_ELEMENT_DEFAULT, COVALENT_RADIUS, BOND_SLACK,
//     RES_CLASS, RES_CLASS_OF, CG_FORMAL_CHARGES  ← same name
//     resolveElementParams / resolveHeavyElementParams  (new: the ONE
//                           element-resolution rule, called by both engines)
//     elementCoverage    (new: which table owns which element)
export {
  ELEMENT_LJ as LIG_ELEMENT,
  ELEMENT_LJ_DEFAULT as LIG_ELEMENT_DEFAULT,
  METAL_ELEMENT,
  METAL_ELEMENT_DEFAULT,
  COVALENT_RADIUS,
  BOND_SLACK,
  RES_CLASS,
  RES_CLASS_OF,
  CG_FORMAL_CHARGES,
  resolveElementParams,
  resolveHeavyElementParams,
  elementCoverage,
} from "./physics/params.js";

// ── Sequence-dependent ENM (Bahar-style) ───────────────────────────────
// Per-residue-class stiffness multiplier w_i (dimensionless). Used by
// forcefield.js springK_seq stub: K = gamma * (1 + 0.2*(w_i+w_j)/2).
// w_i = SEQ_WEIGHT[RES_CLASS_OF[resName]] (w_j likewise). Uniform is
// default; sequence-weighted ENM is opt-in via ForceField.applySeqWeights().
// Scale chosen so (w_i+w_j)/2 ∈ [0.9,1.1] → K ∈ [0.98,1.02]*gamma for
// typical pairs (±2% modulation), up to ±10% for P–P vs A–A extremes,
// preserving mean γ and ENM stability while encoding chemistry. See
// docs/CG_HEAVY.md §CG ENM for validation (1ubq B-factor Pearson).
export const SEQ_WEIGHT = {
  H:  1.0,   // hydrophobic (ALA, VAL, LEU, ILE, PRO, MET, GLY, CYS) — reference
  A:  0.9,   // aromatic (PHE, TRP, TYR, HIS) — softer, π stacking
  P:  1.1,   // polar (SER, THR, ASN, GLN) — stiffer, H-bond network
  Cp: 1.05,  // positively charged (LYS, ARG) — slight stiffening
  Cn: 1.05,  // negatively charged (ASP, GLU) — slight stiffening
};

// Backbone constants — re-exported for integrator / docs provenance
// Derived via Boltzmann inversion (see header). k_B T(300K)=0.6 kcal/mol.
export const KBOND_DEFAULT = 100.0;  // kcal/mol/Å²  — Cα–Cα bond (Tirion 1996; AMBER ff14SB)
export const KANGLE_DEFAULT = 20.0;  // kcal/mol/rad² — Cα pseudo-angle (Tirion 1996)

// ── Holo-pair one-sided compression floor ─────────────────────────────
// A flat-bottom harmonic wall on holo-pinned pairs: no energy while r >=
// HOLO_FLOOR_RMIN, then ½·HOLO_FLOOR_K·(r − r_min)² below it. Stops the
// funnel or the desolvation term from squeezing the ligand through a pocket
// wall (holo pairs are excluded from both grid repulsion and the binding
// pass, so without this term nothing resists r < r_min).
//
// These two numbers were typed into BOTH forcefield.js (the _holoFloor
// kernel) and physics/integrators/respa.js (splitCoarsegrained's bondedFn
// replica) — the duplication is load-bearing for parity, not accidental, but
// a silent edit to one copy is exactly the divergence class this repo has
// been bitten by, so the values live here once. respa.js's docstring
// ("keep in sync with src/forcefield.js if that kernel changes") is now
// satisfied structurally: there is only one copy to keep in sync.
export const HOLO_FLOOR_RMIN = 2.6;   // Å
export const HOLO_FLOOR_K = 8.0;     // kcal/mol/Å²
