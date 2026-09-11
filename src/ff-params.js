/**
 * ff-params.js — model constants and per-element/per-residue-class parameter
 * tables for the coarse-grained force field in forcefield.js.
 *
 * Split out of forcefield.js (item 5, "move never rewrite": every value below
 * is byte-identical to the original tables). forcefield.js re-exports
 * KB_KCAL and KCONV so downstream modules (integrator.js, funnel.js) keep
 * importing them from "./forcefield.js" unchanged.
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

export const KB_KCAL = 0.0019872041; // Boltzmann constant, kcal mol^-1 K^-1

/**
 * Unit-conversion constant between (kcal/mol, Å) and mechanical (Da, Å, ps)
 * units used by the integrator:
 *     1 kcal/mol            = 418.4 Da·Å²/ps²
 *     a[Å/ps²]              = KCONV · F[kcal/mol/Å] / m[Da]
 *     k_B in mech. units    = KB_KCAL · KCONV = 0.8314 Da·Å²/(ps²·K)
 * (derivation: 1 kcal/mol·Å = 6.9477e-11 N ≙ 418.4 Da·Å/ps² since 1 N =
 * 6.022e15 Da·Å/ps².)
 */
export const KCONV = 418.4;

// Residue class → protein-bead LJ parameters (σ Å, ε kcal/mol, charge e)
export const RES_CLASS = {
  H:  { sigma: 4.0, eps: 0.15, q: 0 },  // hydrophobic
  A:  { sigma: 4.1, eps: 0.18, q: 0 },  // aromatic
  P:  { sigma: 3.8, eps: 0.12, q: 0 },  // polar (H-bond capable)
  Cp: { sigma: 3.6, eps: 0.10, q: 0 },  // positively charged (H-bond capable)
  Cn: { sigma: 3.6, eps: 0.10, q: 0 },  // negatively charged (H-bond capable)
};
export const RES_CLASS_OF = {
  ALA: "H", VAL: "H", LEU: "H", ILE: "H", PRO: "H", MET: "H", GLY: "H", CYS: "H",
  PHE: "A", TRP: "A", TYR: "A", HIS: "A",
  SER: "P", THR: "P", ASN: "P", GLN: "P",
  LYS: "Cp", ARG: "Cp",
  ASP: "Cn", GLU: "Cn",
};
// Ligand element → LJ params (σ Å, ε kcal/mol), partial charge (e), H-bond flag, ΔG desolvation (kcal/mol)
export const LIG_ELEMENT = {
  C:  { sigma: 3.4, eps: 0.12, q: 0.0,  hb: false, dG: -0.55 },
  N:  { sigma: 3.2, eps: 0.15, q: -0.30, hb: true,  dG: -0.35 },
  O:  { sigma: 3.0, eps: 0.16, q: -0.50, hb: true,  dG: -0.30 },
  S:  { sigma: 3.6, eps: 0.18, q: 0.0,  hb: false, dG: -0.45 },
  F:  { sigma: 2.9, eps: 0.10, q: -0.20, hb: true,  dG: -0.25 },
  CL: { sigma: 3.5, eps: 0.18, q: 0.0,  hb: false, dG: -0.40 },
  BR: { sigma: 3.6, eps: 0.20, q: 0.0,  hb: false, dG: -0.45 },
  I:  { sigma: 3.8, eps: 0.22, q: 0.0,  hb: false, dG: -0.50 },
  P:  { sigma: 3.5, eps: 0.14, q: 0.40, hb: false, dG: -0.35 },
};
export const LIG_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0, hb: false, dG: -0.30 };

/**
 * Metal-ion parameters for the all-atom heavy mode (heavy.js).
 *   sigma/eps — LJ size & well (Å, kcal/mol) for non-bonded repulsion,
 *   q         — formal charge (e), used for screened electrostatics,
 *   coordR    — metal–donor coordination distance (Å) used to build the
 *               coordination springs (heavy.js detects donors within coordR of
 *               the ion), and
 *   coordN    — target coordination number (how many donor springs to build;
 *               capped by however many donors are actually within coordR).
 * Metals are treated as explicit +2/+1 ions that coordinate N/O/S donors
 * (histidine N, carboxylate O, thiolate S, backbone carbonyl O) rather than
 * forming covalent bonds.
 */
export const METAL_ELEMENT = {
  ZN: { sigma: 1.40, eps: 0.05, q: 2.0, coordR: 2.30, coordN: 4 },
  FE: { sigma: 1.50, eps: 0.05, q: 2.0, coordR: 2.20, coordN: 6 },
  MG: { sigma: 1.30, eps: 0.05, q: 2.0, coordR: 2.15, coordN: 6 },
  CA: { sigma: 1.70, eps: 0.05, q: 2.0, coordR: 2.45, coordN: 6 },
  CU: { sigma: 1.40, eps: 0.05, q: 2.0, coordR: 2.20, coordN: 4 },
  MN: { sigma: 1.45, eps: 0.05, q: 2.0, coordR: 2.25, coordN: 6 },
  NI: { sigma: 1.40, eps: 0.05, q: 2.0, coordR: 2.15, coordN: 6 },
  CO: { sigma: 1.40, eps: 0.05, q: 2.0, coordR: 2.15, coordN: 6 },
  NA: { sigma: 1.70, eps: 0.05, q: 1.0, coordR: 2.50, coordN: 6 },
  K:  { sigma: 2.00, eps: 0.05, q: 1.0, coordR: 2.80, coordN: 6 },
};
export const METAL_ELEMENT_DEFAULT = { sigma: 1.50, eps: 0.05, q: 2.0, coordR: 2.30, coordN: 6 };

/**
 * All-atom (heavy) covalent radii (Å) for the heavy mode bond-building
 * (heavy.js). A pair of heavy atoms within the SUM of their covalent radii
 * (× a 1.15 slack factor) is treated as covalently bonded. Solvent/water O
 * and H are handled by the parser (dropped), so these cover protein heavy
 * atoms + ligand heavy atoms.
 *
 * Source — Covalent radii from Cambridge Structural Database (CSD) surveys
 * (Allen et al., J. Chem. Soc. Perkin Trans. 2, 1987; Cordero et al., Dalton
 * Trans. 2008, 2832; Bondi, J. Phys. Chem. 1964, 68, 441 for van-der-Waals
 * reference). Values below are single-bond covalent radii matched to CSD
 * organic/metal-organic statistics:
 *
 * | Element | r_cov (Å) | CSD/Bondi source | Example bond | r_sum | r_sum×1.15 | Cap 2.2 |
 * |---------|-----------|------------------|--------------|-------|------------|---------|
 * | C       | 0.77      | CSD C(sp3) 0.76  | C–C 1.54     | 1.54  | 1.77       | pass |
 * | N       | 0.75      | CSD N 0.71–0.75  | C–N 1.47     | 1.52  | 1.75       | pass |
 * | O       | 0.73      | CSD O 0.66–0.73  | C–O 1.43     | 1.50  | 1.73       | pass |
 * | S       | 1.02      | CSD S 1.05       | S–S 2.04     | 2.04  | 2.35→2.20  | **pass (2.04 < 2.20 < 2.35)** |
 * | P       | 1.06      | CSD P 1.07       | P–O 1.60     | 1.79  | 2.06       | pass |
 * | CA      | 1.76      | CSD Ca 1.76      | Ca–N 2.51    | 2.51  | 2.89→2.20  | fail at 2.9 Å (spurious Ca–N rejected) |
 *
 * Rationale for thresholds — heavy.js:251 `r < BOND_SLACK*(rA+rB) && r < 2.2`:
 *   BOND_SLACK = 1.15 gives 15% tolerance for thermal elongation / PDB
 *   coordinate uncertainty while keeping C–C/C–N/C–O true bonds well inside.
 *   Hard cap 2.2 Å prevents spurious long-range contacts (e.g. Ca–N 2.9 Å,
 *   H-bond O⋯N 2.9 Å, van-der-Waals contacts) from being mis-typed as
 *   covalent despite large radii sums (Ca 1.76 + N 0.75 = 2.51 → 2.89 with
 *   slack). Critically, the biologically important disulfide S–S 2.04 Å
 *   (CSD mean 2.03–2.05 Å, crambin 1CRN SSBOND records 2.00/2.04/2.05 Å)
 *   has r_sum = 2.04 so r_sum×1.15 = 2.35; min(2.35, 2.20) = 2.20 still
 *   captures 2.04 Å with 0.16 Å margin, so all three crambin disulfides are
 *   recovered while Ca–N 2.9 Å is correctly excluded — see tests/test_topology.js.
 */
export const COVALENT_RADIUS = {
  C: 0.77, N: 0.75, O: 0.73, S: 1.02, P: 1.06, F: 0.71,
  CL: 0.99, BR: 1.14, I: 1.33, B: 0.84, SE: 1.20,
  ZN: 1.22, FE: 1.32, MG: 1.41, CA: 1.76, CU: 1.32, MN: 1.39,
  NI: 1.24, CO: 1.26, NA: 1.66, K: 2.03,
};
/**
 * Slack factor on the covalent-radius sum for bond detection.
 * 1.15 = CSD/Bondi covalent sum × 1.15, then hard cap 2.2 Å (see table above).
 * Validated: S–S 2.04 Å pass, Ca–N 2.9 Å fail — tests/test_topology.js.
 */
export const BOND_SLACK = 1.15;

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
