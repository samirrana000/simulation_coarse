/**
 * ff-params.js — model constants and per-element/per-residue-class parameter
 * tables for the coarse-grained force field in forcefield.js.
 *
 * Split out of forcefield.js (item 5, "move never rewrite": every value below
 * is byte-identical to the original tables). forcefield.js re-exports
 * KB_KCAL and KCONV so downstream modules (integrator.js, funnel.js) keep
 * importing them from "./forcefield.js" unchanged.
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
 */
export const COVALENT_RADIUS = {
  C: 0.77, N: 0.75, O: 0.73, S: 1.02, P: 1.06, F: 0.71,
  CL: 0.99, BR: 1.14, I: 1.33, B: 0.84, SE: 1.20,
  ZN: 1.22, FE: 1.32, MG: 1.41, CA: 1.76, CU: 1.32, MN: 1.39,
  NI: 1.24, CO: 1.26, NA: 1.66, K: 2.03,
};
/** Slack factor on the covalent-radius sum for bond detection. */
export const BOND_SLACK = 1.15;
