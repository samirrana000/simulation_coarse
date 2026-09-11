/**
 * charges.js — Partial charge assignments, Generalized Born radii, and SASA parameters.
 *
 * AMBER ff14SB reference:
 *   Maier, J.A. et al. "ff14SB: Improving the Accuracy of Protein Side Chain
 *   and Backbone Parameters from ff99SB." J. Chem. Theory Comput. 2015, 11,
 *   3696-3713. DOI:10.1021/acs.jctc.5b00255
 *   Original ff94 charge derivation: Cornell et al. JACS 1995, 117, 5179.
 *   This module uses an approximate, united-atom mapping of ff14SB partial
 *   charges (heavy-atom only; hydrogens summed into heavy). Values are
 *   heuristic/rounded for coarse-grained visualization and teaching — see
 *   docs/CHARGES.md. NOT a full ff14SB force-field implementation.
 *
 * Provides standard AMBER ff14SB / united-atom partial charges for standard
 * amino acid residues, terminal caps, ions/metals, and fallback rules for ligands.
 *
 * Residue net-charge validation (internal, without OXT):
 *   ASP  -1.00 (calc -1.00) | GLU  -1.00 (calc -1.00)
 *   ARG  +1.00 (calc +1.00) | LYS  +1.00 (calc +1.00)
 *   — validated by validateCharges() below; OXT (-0.80) is a C-terminal cap
 *   and excluded from the internal-residue sum.
 */

// United-atom AMBER/GROMOS partial charges (neutral residues = 0.0, ARG/LYS = +1.0, ASP/GLU = -1.0)
// Approximated from ff14SB (heavy-atom sums). See docs/CHARGES.md for table and caveats.
// Formal net charges: ASP/GLU = -1, ARG/LYS = +1, others = 0 (internal, no OXT/termini).
export const AA_CHARGES = {
  ALA: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, OXT: -0.8 },
  ARG: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, CD: 0.10, NE: -0.40, CZ: 0.64, NH1: 0.33, NH2: 0.33, OXT: -0.8 },
  ASN: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.55, OD1: -0.55, ND2: 0.0, OXT: -0.8 },
  ASP: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.60, OD1: -0.80, OD2: -0.80, OXT: -0.8 },
  CYS: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.10, SG: -0.10, OXT: -0.8 },
  GLN: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, CD: 0.55, OE1: -0.55, NE2: 0.0, OXT: -0.8 },
  GLU: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, CD: 0.60, OE1: -0.80, OE2: -0.80, OXT: -0.8 },
  GLY: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, OXT: -0.8 },
  HIS: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, ND1: -0.36, CD2: 0.18, CE1: 0.18, NE2: 0.0, OXT: -0.8 },
  ILE: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG1: 0.0, CG2: 0.0, CD1: 0.0, OXT: -0.8 },
  LEU: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, CD1: 0.0, CD2: 0.0, OXT: -0.8 },
  LYS: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, CD: 0.0, CE: 0.25, NZ: 0.75, OXT: -0.8 },
  MET: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.06, SD: -0.12, CE: 0.06, OXT: -0.8 },
  PHE: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, CD1: 0.0, CD2: 0.0, CE1: 0.0, CE2: 0.0, CZ: 0.0, OXT: -0.8 },
  PRO: { N: -0.20, CA: 0.20, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, CD: 0.0, OXT: -0.8 },
  SER: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.20, OG: -0.20, OXT: -0.8 },
  THR: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.20, OG1: -0.20, CG2: 0.0, OXT: -0.8 },
  TRP: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, CD1: 0.10, CD2: 0.0, NE1: -0.10, CE2: 0.0, CE3: 0.0, CZ2: 0.0, CZ3: 0.0, CH2: 0.0, OXT: -0.8 },
  TYR: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG: 0.0, CD1: 0.0, CD2: 0.0, CE1: 0.0, CE2: 0.0, CZ: 0.20, OH: -0.20, OXT: -0.8 },
  VAL: { N: -0.28, CA: 0.28, C: 0.55, O: -0.55, CB: 0.0, CG1: 0.0, CG2: 0.0, OXT: -0.8 },
};

// Generalized Born intrinsic atomic radii (Å)
export const GB_RADII = {
  H: 1.20,
  C: 1.70,
  N: 1.55,
  O: 1.50,
  F: 1.47,
  P: 1.85,
  S: 1.80,
  CL: 1.75,
  BR: 1.85,
  I: 1.98,
  ZN: 1.20,
  FE: 1.30,
  MG: 1.18,
  CA: 1.52,
  CU: 1.20,
  MN: 1.25,
  NI: 1.20,
  CO: 1.20,
  NA: 1.40,
  K: 1.70,
  DEFAULT: 1.60,
};

// Atomic surface area probe parameters (Å)
export const SASA_RADII = {
  C: 1.70,
  N: 1.55,
  O: 1.52,
  S: 1.80,
  P: 1.80,
  F: 1.47,
  CL: 1.75,
  BR: 1.85,
  I: 1.98,
  DEFAULT: 1.60,
};

/**
 * Assign physically realistic partial charges to heavy atoms.
 * Protein atoms receive standard AMBER charges; ligands retain parsed charges or element defaults.
 *
 * @param {Array<{atomName: string, resName: string, element: string, isProtein: boolean, charge?: number}>} atoms
 * @returns {Float64Array} array of partial charges in elementary charge (e)
 */
export function assignCharges(atoms) {
  const n = atoms.length;
  const q = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    const a = atoms[i];
    if (a.isProtein) {
      const res = AA_CHARGES[a.resName];
      if (res && res[a.atomName] !== undefined) {
        q[i] = res[a.atomName];
      } else {
        // Fallback by element
        if (a.element === "O") q[i] = -0.5;
        else if (a.element === "N") q[i] = -0.3;
        else if (a.element === "S") q[i] = -0.2;
        else q[i] = 0.0;
      }
    } else if (a.isMetal) {
      if (a.element === "NA" || a.element === "K") q[i] = 1.0;
      else q[i] = 2.0;
    } else {
      // Ligand or cofactor
      if (a.charge !== undefined && !Number.isNaN(a.charge) && a.charge !== 0) {
        q[i] = a.charge;
      } else {
        if (a.element === "O") q[i] = -0.45;
        else if (a.element === "N") q[i] = -0.30;
        else if (a.element === "F") q[i] = -0.20;
        else if (a.element === "CL" || a.element === "BR") q[i] = -0.10;
        else if (a.element === "P") q[i] = 0.40;
        else q[i] = 0.0;
      }
    }
  }

  return q;
}

/**
 * Validation: check that 4 charged residues reproduce formal integer net charges.
 *
 * Internal residues (no termini) should sum to formal charge when OXT is excluded.
 * Formal: ASP -1, GLU -1, ARG +1, LYS +1. Use literature ff14SB rounded values
 * (e.g., ASP OD1 -0.8014, OD2 -0.8014, CG 0.7172) approximated here as OD -0.80,
 * CG +0.60 — see docs/CHARGES.md for note "approximate, not full ff14SB".
 *
 * @param {number} [tol=0.05] tolerance in e for sum vs. formal
 * @returns {{residue:string, calc:number, formal:number, ok:boolean}[]}
 */
export function validateCharges(tol = 0.05) {
  const formal = { ASP: -1.0, GLU: -1.0, ARG: +1.0, LYS: +1.0 };
  const results = [];
  for (const [res, qexp] of Object.entries(formal)) {
    const table = AA_CHARGES[res];
    let sum = 0;
    for (const [atom, q] of Object.entries(table)) {
      if (atom === "OXT") continue; // terminal cap excluded
      sum += q;
    }
    const ok = Math.abs(sum - qexp) < tol;
    results.push({ residue: res, calc: sum, formal: qexp, ok });
  }
  return results;
}

/**
 * Table comment: residue-level net charges (heavy-atom sums, no OXT, no hydrogens):
 *
 * | Residue | Atoms summed                          | Net calc | Formal |
 * |---------|---------------------------------------|----------|--------|
 * | ASP     | N CA C O CB CG OD1 OD2                | -1.00    |  -1    |
 * | GLU     | N CA C O CB CG CD OE1 OE2             | -1.00    |  -1    |
 * | ARG     | N CA C O CB CG CD NE CZ NH1 NH2      | +1.00    |  +1    |
 * | LYS     | N CA C O CB CG CD CE NZ               | +1.00    |  +1    |
 * Literature ff14SB ASP example: OD1 -0.8014 OD2 -0.8014 CG 0.7172 CB -0.2826
 * → rounded here to OD -0.80 CG 0.60 CB 0.0 for united-atom teaching model.
 * Statement: "approximate, not full ff14SB" — full ff14SB would include
 * hydrogens, separate Nter/Cter patches, and exact Cornell et al. charges.
 */
