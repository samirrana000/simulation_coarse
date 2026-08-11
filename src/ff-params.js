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
