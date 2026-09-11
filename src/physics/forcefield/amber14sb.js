/**
 * amber14sb.js — AMBER ff14SB-style bonded + non-bonded parameter tables.
 *
 * Pure data + getters, no Node dependencies, ES-module safe (`?v=` compatible).
 * Units: Length Å, Energy kcal/mol, Angle rad internally (tables store deg,
 * getters convert), Mass Da. See src/units.js and docs/UNITS.md.
 *
 * Scope: standard amino-acid heavy-atom topology. Hydrogens are united into
 * heavy atoms (this engine drops H); listed bonds/angles/dihedrals involve
 * heavy atoms only (N, CA, C, O, CB, CG, ...). Values are ff14SB/parm99/GAFF
 * consensus rounded for a geometry-driven engine — suitable as physically
 * grounded defaults for src/heavy.js covalent terms, NOT a byte-exact AMBER
 * port (full ff14SB needs atom types CT, C2, N-star and H plus 1-4 scalings).
 *
 * References:
 *   [1] Maier et al., JCTC 11, 3696 (2015) — ff14SB.
 *   [2] Cornell et al., JACS 117, 5179 (1995) — ff94/parm99 base.
 *   [3] Wang et al., JCC 25, 1157 (2004) — GAFF bond/angle complements.
 *   [4] Case et al., AMBER 2020 Manual — bond/angle/dihedral functional forms.
 *
 * Functional forms used by callers:
 *   Bond:     U = k (r - r0)^2            (AMBER uses k with 1/2 folded in;
 *             here k is the full harmonic constant as in heavy.js harmonicFlat,
 *             i.e. U = 1/2 K (r-r0)^2 with K = 2k_AMBER. Tables below store K
 *             directly (kcal/mol/A^2) matching heavy.js kBond semantics.)
 *   Angle:    U = 1/2 k (th - th0)^2      (k in kcal/mol/rad^2)
 *   Proper:   U = sum_m (V_m/2) [1 + cos(n_m φ - γ_m)]   (Fourier, multi-term)
 *   Improper: U = 1/2 k (φ - φ0)^2        (harmonic planar restraint)
 *
 * Key conventions:
 *   - Atom-name keys use PDB heavy-atom names (N, CA, C, O, CB, ...).
 *   - Element-fallback keys use uppercase elements (C, N, O, S, ...).
 *   - Bond keys are order-independent ("CA-C" == "C-CA").
 *   - Angle keys are reversal-symmetric ("N-CA-C" == "C-CA-N").
 *   - Dihedral keys support `X` wildcards ("X-CA-CB-X"); matching tries exact,
 *     then single/double wildcard, then element-fallback.
 */

/** AMBER ff14SB heavy-atom masses (Da, IUPAC). */
export const AMBER_MASS = {
  C: 12.011, N: 14.007, O: 15.999, S: 32.06, P: 30.974,
  F: 18.998, CL: 35.45, BR: 79.904, I: 126.904,
  ZN: 65.38, FE: 55.845, MG: 24.305, CA: 40.078, CU: 63.546,
  MN: 54.938, NI: 58.693, CO: 58.933, NA: 22.99, K: 39.098,
  SE: 78.971,
};

/**
 * Bond table: key "A-B" (sorted) -> { r0: Å, k: kcal/mol/Å² }.
 * Specific atom-name entries first; ELEMENT_BONDS provides element fallback.
 * Sources: ff14SB parm99 bond list (CT-CT 310/1.526, C-N 490/1.335,
 * C-O 570/1.229, N-CT 337/1.449, S-S 166/2.038, C-S 237/1.810).
 */
export const BOND_TABLE = {
  // --- peptide backbone (atom-name specific) ---
  "C-N":   { r0: 1.335, k: 490.0 },
  "C-O":   { r0: 1.229, k: 570.0 },
  "C-CA":  { r0: 1.522, k: 317.0 },
  "CA-C":  { r0: 1.522, k: 317.0 },
  "CA-N":  { r0: 1.449, k: 337.0 },
  "N-CA":  { r0: 1.449, k: 337.0 },
  "CA-CB": { r0: 1.526, k: 310.0 },
  "CB-CA": { r0: 1.526, k: 310.0 },
  // --- side-chain / generic heavy ---
  "CB-CG":   { r0: 1.526, k: 310.0 },
  "CG-CD":   { r0: 1.526, k: 310.0 },
  "CD-CE":   { r0: 1.526, k: 310.0 },
  "CE-NZ":   { r0: 1.471, k: 337.0 },
  "CG-OD1":  { r0: 1.229, k: 570.0 },
  "CG-OD2":  { r0: 1.250, k: 570.0 },
  "CD-OE1":  { r0: 1.229, k: 570.0 },
  "CD-OE2":  { r0: 1.250, k: 570.0 },
  "CG-ND2":  { r0: 1.335, k: 490.0 },
  "CD-NE2":  { r0: 1.335, k: 490.0 },
  "CB-SG":   { r0: 1.810, k: 237.0 },
  "SG-SG":   { r0: 2.038, k: 166.0 },
  "CB-OG":   { r0: 1.430, k: 320.0 },
  "CB-OG1":  { r0: 1.430, k: 320.0 },
  "CG-SD":   { r0: 1.810, k: 237.0 },
  "SD-CE":   { r0: 1.810, k: 237.0 },
  "CG-ND1":  { r0: 1.380, k: 400.0 },
  "CE1-NE2": { r0: 1.380, k: 400.0 },
  "CD2-NE2": { r0: 1.380, k: 400.0 },
  "CZ-NH1":  { r0: 1.335, k: 490.0 },
  "CZ-NH2":  { r0: 1.335, k: 490.0 },
  "NE-CZ":   { r0: 1.335, k: 490.0 },
  "CD-NE":   { r0: 1.471, k: 337.0 },
  "CZ-OH":   { r0: 1.410, k: 320.0 },
  // aromatic C-C (PHE/TYR/TRP/HIS ring)
  "CG-CD1": { r0: 1.400, k: 469.0 },
  "CG-CD2": { r0: 1.400, k: 469.0 },
  "CD1-CE1": { r0: 1.400, k: 469.0 },
  "CD2-CE2": { r0: 1.400, k: 469.0 },
  "CE1-CZ":  { r0: 1.400, k: 469.0 },
  "CE2-CZ":  { r0: 1.400, k: 469.0 },
};

/** Element-fallback bonds: key sorted "C-N" etc. -> { r0, k }. */
export const ELEMENT_BONDS = {
  "C-C": { r0: 1.510, k: 310.0 },
  "C-N": { r0: 1.440, k: 350.0 },
  "C-O": { r0: 1.300, k: 450.0 },
  "C-S": { r0: 1.810, k: 237.0 },
  "C-P": { r0: 1.800, k: 230.0 },
  "N-N": { r0: 1.450, k: 300.0 },
  "N-O": { r0: 1.420, k: 300.0 },
  "O-O": { r0: 1.480, k: 300.0 },
  "O-P": { r0: 1.600, k: 300.0 },
  "S-S": { r0: 2.038, k: 166.0 },
  "C-F": { r0: 1.350, k: 400.0 },
  "C-CL": { r0: 1.760, k: 250.0 },
  "C-BR": { r0: 1.940, k: 220.0 },
  "C-I": { r0: 2.100, k: 200.0 },
};

/** Default bond when nothing matches (heavy.js legacy k=200-300 range). */
export const BOND_DEFAULT = { r0: 1.500, k: 300.0 };

/**
 * Angle table: key "A-B-C" with B apex, reversal-symmetric.
 * { th0: degrees, k: kcal/mol/rad² }. ff14SB angle K ~ 40-80.
 */
export const ANGLE_TABLE = {
  "N-CA-C":   { th0: 111.2, k: 80.0 },
  "CA-C-N":   { th0: 116.6, k: 70.0 },
  "CA-C-O":   { th0: 120.4, k: 80.0 },
  "O-C-N":    { th0: 122.9, k: 80.0 },
  "C-N-CA":   { th0: 121.9, k: 50.0 },
  "N-CA-CB":  { th0: 110.1, k: 80.0 },
  "C-CA-CB":  { th0: 110.1, k: 80.0 },
  "CA-CB-CG": { th0: 113.6, k: 58.0 },
  "CB-CG-CD": { th0: 113.6, k: 58.0 },
  "CA-CB-OG": { th0: 109.5, k: 50.0 },
  "CA-CB-SG": { th0: 114.0, k: 50.0 },
  "CB-SG-SG": { th0: 104.0, k: 60.0 },
  "CB-CG-OD1": { th0: 118.0, k: 80.0 },
  "CB-CG-OD2": { th0: 118.0, k: 80.0 },
  "OD1-CG-OD2": { th0: 124.0, k: 80.0 },
  "CG-CD-NE":  { th0: 112.0, k: 50.0 },
  "CD-NE-CZ":  { th0: 123.0, k: 50.0 },
  "NE-CZ-NH1": { th0: 120.0, k: 80.0 },
  "NE-CZ-NH2": { th0: 120.0, k: 80.0 },
  "NH1-CZ-NH2": { th0: 120.0, k: 80.0 },
  "CB-CG-ND1": { th0: 122.0, k: 60.0 },
  "CB-CG-CD1": { th0: 120.0, k: 60.0 },
  "CB-CG-CD2": { th0: 120.0, k: 60.0 },
  "CG-CD1-CE1": { th0: 120.0, k: 60.0 },
  "CD1-CE1-CZ":  { th0: 120.0, k: 60.0 },
  "CE1-CZ-OH":   { th0: 120.0, k: 60.0 },
};

/** Element-fallback angles (apex B = middle element). */
export const ELEMENT_ANGLES = {
  "C-C-C": { th0: 112.0, k: 58.0 },
  "C-C-N": { th0: 111.0, k: 60.0 },
  "C-C-O": { th0: 115.0, k: 60.0 },
  "N-C-N": { th0: 120.0, k: 60.0 },
  "N-C-O": { th0: 122.0, k: 70.0 },
  "C-N-C": { th0: 122.0, k: 50.0 },
  "C-O-P": { th0: 120.0, k: 50.0 },
  "C-S-C": { th0: 100.0, k: 50.0 },
  "C-S-S": { th0: 104.0, k: 60.0 },
  "C-C-S": { th0: 114.0, k: 50.0 },
  "C-C-F": { th0: 109.5, k: 50.0 },
  "O-C-O": { th0: 126.0, k: 80.0 },
  "O-P-O": { th0: 109.5, k: 50.0 },
};

/** Default angle (heavy.js legacy kAngle=40). */
export const ANGLE_DEFAULT = { th0: 112.0, k: 40.0 };

/**
 * Proper dihedral table: key "A-B-C-D" (`X` = wildcard).
 * Value: array of Fourier terms [{ Vn: kcal/mol, n: int, gamma: degrees }].
 * Energy: U = sum_m (V_m/2)[1 + cos(n_m φ - γ_m)].
 *
 * Backbone (ff99SB/ff14SB consensus, heavy-atom):
 *   omega C-CA-N-C / CA-C-N-CA: strong trans barrier V2 ~ 10-15, γ=180.
 *   phi   C-N-CA-C: multi-term (V1/V2/V3 small, backbone-corrected).
 *   psi   N-CA-C-N: multi-term.
 * Side-chain chi (CT-CT-CT-X): V3 ~ 0.5-2.0, n=3, γ=0.
 */
export const DIHEDRAL_TABLE = {
  // peptide omega — trans planar, large V2 barrier
  "X-C-N-X": [
    { Vn: 12.0, n: 2, gamma: 180.0 },
  ],
  "CA-C-N-CA": [
    { Vn: 12.0, n: 2, gamma: 180.0 },
  ],
  "C-CA-N-C": [
    { Vn: 12.0, n: 2, gamma: 180.0 },
  ],
  // backbone phi C-N-CA-C (ff99SB-style multi-term)
  "C-N-CA-C": [
    { Vn: 0.60, n: 1, gamma: 0.0 },
    { Vn: 1.20, n: 2, gamma: 180.0 },
    { Vn: 0.20, n: 3, gamma: 0.0 },
  ],
  // backbone psi N-CA-C-N (multi-term)
  "N-CA-C-N": [
    { Vn: 0.50, n: 1, gamma: 0.0 },
    { Vn: 1.00, n: 2, gamma: 180.0 },
    { Vn: 0.25, n: 3, gamma: 0.0 },
  ],
  // side-chain chi1 N-CA-CB-CG
  "N-CA-CB-CG": [
    { Vn: 1.80, n: 3, gamma: 0.0 },
  ],
  "N-CA-CB-OG": [
    { Vn: 1.20, n: 3, gamma: 0.0 },
  ],
  "N-CA-CB-SG": [
    { Vn: 1.60, n: 3, gamma: 0.0 },
  ],
  // side-chain chi2 CA-CB-CG-CD
  "CA-CB-CG-CD": [
    { Vn: 1.80, n: 3, gamma: 0.0 },
  ],
  "CA-CB-CG-OD1": [
    { Vn: 1.00, n: 3, gamma: 0.0 },
  ],
  "CA-CB-CG-ND1": [
    { Vn: 1.00, n: 3, gamma: 0.0 },
  ],
  // generic sp3 chain
  "X-CT-CT-X": [
    { Vn: 1.40, n: 3, gamma: 0.0 },
  ],
  "X-CA-CB-X": [
    { Vn: 1.40, n: 3, gamma: 0.0 },
  ],
  "X-CB-CG-X": [
    { Vn: 1.40, n: 3, gamma: 0.0 },
  ],
  // aromatic ring torsions kept near-planar via impropers; small proper kept
  "X-CG-CD1-X": [
    { Vn: 2.00, n: 2, gamma: 180.0 },
  ],
};

/** Default proper (alkane-like) when nothing matches. */
export const DIHEDRAL_DEFAULT = [{ Vn: 1.40, n: 3, gamma: 0.0 }];

/**
 * Improper table: harmonic planar restraints U = 1/2 k (φ-φ0)².
 * Key "A-B-C-D" with C the central atom (AMBER convention i-j-k-l, k central).
 * φ0 in degrees; k in kcal/mol/rad².
 */
export const IMPROPER_TABLE = {
  // peptide plane: O-C-N-CA and C-N-CA-C kept planar (trans 180)
  "O-C-N-CA":  { k: 20.0, phase: 180.0 },
  "C-N-CA-C":  { k: 20.0, phase: 180.0 },
  "CA-C-N-CA": { k: 20.0, phase: 180.0 },
  // carboxylate / guanidinium / aromatic planes (φ0 = 0)
  "OD1-CG-OD2-CB": { k: 20.0, phase: 0.0 },
  "NH1-CZ-NH2-NE": { k: 20.0, phase: 0.0 },
  "CD1-CG-CD2-CB": { k: 15.0, phase: 0.0 },
  "CE1-CG-CE2-CZ": { k: 15.0, phase: 0.0 },
  "X-CG-CD1-X": { k: 15.0, phase: 0.0 },
  "X-CA-N-X":   { k: 10.0, phase: 180.0 },
};

/** Default improper (weak planar). */
export const IMPROPER_DEFAULT = { k: 10.0, phase: 0.0 };

/**
 * Non-bonded LJ per element (AMBER parm99 Rmin/2 + eps converted to σ).
 * σ = Rmin/2 / 2^(1/6)... stored directly as σ (Å) for 4ε[(σ/r)^12-(σ/r)^6].
 * Rmin/2 values: C 1.908/0.1094, N 1.824/0.1700, O 1.661/0.2100,
 * S 2.000/0.2500, P 2.100/0.2000, F 1.750/0.0610, CL 1.948/0.2650.
 */
export const NONBONDED_TABLE = {
  C:  { sigma: 3.400, eps: 0.1094 },
  N:  { sigma: 3.250, eps: 0.1700 },
  O:  { sigma: 2.960, eps: 0.2100 },
  S:  { sigma: 3.563, eps: 0.2500 },
  P:  { sigma: 3.742, eps: 0.2000 },
  F:  { sigma: 3.119, eps: 0.0610 },
  CL: { sigma: 3.470, eps: 0.2650 },
  BR: { sigma: 3.600, eps: 0.3200 },
  I:  { sigma: 3.800, eps: 0.4000 },
  SE: { sigma: 3.600, eps: 0.2500 },
};
export const NONBONDED_DEFAULT = { sigma: 3.400, eps: 0.1200 };

/* ------------------------------------------------------------------ */
/* Normalization helpers                                              */
/* ------------------------------------------------------------------ */

/**
 * Normalize a bond/angle/dihedral atom spec to an uppercase name string.
 * @param {string|object} a  element string, atom-name string, or {element, atomName}
 * @returns {string}
 */
export function normName(a) {
  if (typeof a === "string") return a.trim().toUpperCase();
  if (a && typeof a === "object") {
    const nm = (a.atomName ?? a.name ?? a.type ?? a.element ?? "").toString().trim().toUpperCase();
    return nm || "C";
  }
  return "C";
}

/**
 * Normalize to element symbol (for ELEMENT_* fallback).
 * @param {string|object} a
 * @returns {string}
 */
export function normElement(a) {
  if (typeof a === "string") {
    const s = a.trim().toUpperCase();
    if (s.length <= 2 && /^[A-Z]{1,2}$/.test(s)) {
      // Could be atom name (CA) or element (CA=calcium). Single letters are
      // elements; two-letter atom names starting with C (CA/CB/CG...) are carbon.
      if (s.length === 1) return s;
      if (s === "CA" || s === "CB" || s === "CG" || s === "CD" || s === "CE" || s === "CZ" || s === "CT") return "C";
      if (s === "CL" || s === "BR") return s;
      if (s === "NA" || s === "ZN" || s === "FE" || s === "MG" || s === "CU" || s === "MN" || s === "NI" || s === "CO") return s;
      return s[0];
    }
    return s;
  }
  if (a && typeof a === "object") {
    const el = (a.element ?? "").toString().trim().toUpperCase();
    if (el) return el;
    return normElement(normName(a));
  }
  return "C";
}

function bondKey(a, b) {
  const A = normName(a), B = normName(b);
  return A < B ? `${A}-${B}` : `${B}-${A}`;
}

function elementBondKey(a, b) {
  const A = normElement(a), B = normElement(b);
  return A < B ? `${A}-${B}` : `${B}-${A}`;
}

function angleKey(a, b, c) {
  const A = normName(a), B = normName(b), C = normName(c);
  const fwd = `${A}-${B}-${C}`, rev = `${C}-${B}-${A}`;
  return { fwd, rev };
}

/* ------------------------------------------------------------------ */
/* Lookup API                                                         */
/* ------------------------------------------------------------------ */

/**
 * Bond parameters for atom pair (a, b).
 * Tries atom-name table, then element fallback, then default.
 * @param {string|object} a
 * @param {string|object} b
 * @returns {{r0:number,k:number}}  r0 in Å, k in kcal/mol/Å²
 */
export function getBondParams(a, b) {
  const k = bondKey(a, b);
  if (BOND_TABLE[k]) return { ...BOND_TABLE[k] };
  const ek = elementBondKey(a, b);
  if (ELEMENT_BONDS[ek]) return { ...ELEMENT_BONDS[ek] };
  // Try element-of-name cross (e.g. CA is carbon): map names to elements first
  const eA = normElement(a), eB = normElement(b);
  const ek2 = eA < eB ? `${eA}-${eB}` : `${eB}-${eA}`;
  if (ELEMENT_BONDS[ek2]) return { ...ELEMENT_BONDS[ek2] };
  return { ...BOND_DEFAULT };
}

/**
 * Angle parameters for triple (a, b, c), b = apex.
 * @param {string|object} a
 * @param {string|object} b
 * @param {string|object} c
 * @returns {{th0:number,k:number}}  th0 in radians, k in kcal/mol/rad²
 */
export function getAngleParams(a, b, c) {
  const D2R = Math.PI / 180;
  const { fwd, rev } = angleKey(a, b, c);
  const hit = ANGLE_TABLE[fwd] ?? ANGLE_TABLE[rev];
  if (hit) return { th0: hit.th0 * D2R, k: hit.k };
  const eA = normElement(a), eB = normElement(b), eC = normElement(c);
  const ef = `${eA}-${eB}-${eC}`, er = `${eC}-${eB}-${eA}`;
  const ehit = ELEMENT_ANGLES[ef] ?? ELEMENT_ANGLES[er];
  if (ehit) return { th0: ehit.th0 * D2R, k: ehit.k };
  return { th0: ANGLE_DEFAULT.th0 * D2R, k: ANGLE_DEFAULT.k };
}

/**
 * Proper Fourier dihedral terms for quadruple (a, b, c, d).
 * Always returns an array (multi-term); empty never (falls back to default).
 * @param {string|object} a
 * @param {string|object} b
 * @param {string|object} c
 * @param {string|object} d
 * @returns {Array<{Vn:number,n:number,gamma:number}>}  gamma in radians
 */
export function getDihedralParams(a, b, c, d) {
  const D2R = Math.PI / 180;
  const A = normName(a), B = normName(b), C = normName(c), D = normName(d);
  const exact = `${A}-${B}-${C}-${D}`;
  if (DIHEDRAL_TABLE[exact]) {
    return DIHEDRAL_TABLE[exact].map((t) => ({ Vn: t.Vn, n: t.n, gamma: t.gamma * D2R }));
  }
  // Wildcard sweep: X-B-C-X, A-B-C-X, X-B-C-D, then element-level
  const pats = [
    `X-${B}-${C}-X`,
    `${A}-${B}-${C}-X`,
    `X-${B}-${C}-${D}`,
    `${A}-${B}-${C}-${D}`,
  ];
  for (const p of pats) {
    if (DIHEDRAL_TABLE[p]) {
      return DIHEDRAL_TABLE[p].map((t) => ({ Vn: t.Vn, n: t.n, gamma: t.gamma * D2R }));
    }
  }
  const eA = normElement(a), eB = normElement(b), eC = normElement(c), eD = normElement(d);
  const epats = [
    `X-${eB}-${eC}-X`,
    `${eA}-${eB}-${eC}-${eD}`,
  ];
  for (const p of epats) {
    if (DIHEDRAL_TABLE[p]) {
      return DIHEDRAL_TABLE[p].map((t) => ({ Vn: t.Vn, n: t.n, gamma: t.gamma * D2R }));
    }
  }
  return DIHEDRAL_DEFAULT.map((t) => ({ Vn: t.Vn, n: t.n, gamma: t.gamma * D2R }));
}

/**
 * Harmonic improper planar restraint for quadruple (a, b, c, d), c central.
 * @param {string|object} a
 * @param {string|object} b
 * @param {string|object} c
 * @param {string|object} d
 * @returns {{k:number,phase:number}}  k kcal/mol/rad², phase radians
 */
export function getImproperParams(a, b, c, d) {
  const D2R = Math.PI / 180;
  const A = normName(a), B = normName(b), C = normName(c), D = normName(d);
  const exact = `${A}-${B}-${C}-${D}`;
  if (IMPROPER_TABLE[exact]) {
    const t = IMPROPER_TABLE[exact];
    return { k: t.k, phase: t.phase * D2R };
  }
  for (const key of Object.keys(IMPROPER_TABLE)) {
    // wildcard match where stored key has X entries
    const parts = key.split("-");
    const q = [A, B, C, D];
    let ok = true;
    for (let i = 0; i < 4; i++) {
      if (parts[i] !== "X" && parts[i] !== q[i]) { ok = false; break; }
    }
    if (ok) {
      const t = IMPROPER_TABLE[key];
      return { k: t.k, phase: t.phase * D2R };
    }
  }
  return { k: IMPROPER_DEFAULT.k, phase: IMPROPER_DEFAULT.phase * D2R };
}

/**
 * Non-bonded LJ parameters for an element.
 * @param {string|object} a
 * @returns {{sigma:number,eps:number}}
 */
export function getNonbondedParams(a) {
  const e = normElement(a);
  if (NONBONDED_TABLE[e]) return { ...NONBONDED_TABLE[e] };
  return { ...NONBONDED_DEFAULT };
}

/**
 * Heavy-atom mass (Da).
 * @param {string|object} a
 * @returns {number}
 */
export function getMass(a) {
  const e = normElement(a);
  return AMBER_MASS[e] ?? 14.0;
}

/**
 * Proper-dihedral Fourier energy for one torsion (no forces).
 * U = Σ_m (V_m/2)[1 + cos(n_m φ − γ_m)].
 * @param {number} phi  dihedral angle in radians
 * @param {Array<{Vn:number,n:number,gamma:number}>} terms
 * @returns {number}  energy in kcal/mol
 */
export function dihedralFourierEnergy(phi, terms) {
  let U = 0;
  for (const t of terms) U += 0.5 * t.Vn * (1 + Math.cos(t.n * phi - t.gamma));
  return U;
}

/**
 * Derivative dU/dφ for Fourier dihedral (for force projection).
 * dU/dφ = Σ_m −(V_m n_m /2) sin(n_m φ − γ_m).
 * @param {number} phi
 * @param {Array<{Vn:number,n:number,gamma:number}>} terms
 * @returns {number}
 */
export function dihedralFourierTorque(phi, terms) {
  let t = 0;
  for (const m of terms) t += -0.5 * m.Vn * m.n * Math.sin(m.n * phi - m.gamma);
  return t;
}
