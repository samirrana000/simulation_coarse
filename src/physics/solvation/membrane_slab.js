/**
 * membrane_slab.js — IMM1/HDGB-style implicit bilayer slab.
 *
 * Flat infinite bilayer with hydrophobic core |z| ≤ T (T = 15 Å half-thickness,
 * total 30 Å) centered at z = 0, smooth water↔membrane transitions of width w.
 *
 * Profiles (all C∞ via tanh/sech², so forces are continuous):
 *   S(z)  = 1/2·[1 + tanh((|z| − T)/w)]      water fraction (0 core → 1 water)
 *   M(z)  = 1 − S(z)                          membrane burial fraction
 *   ε(z)  = ε_mem + (ε_water − ε_mem)·S(z)    depth-dependent dielectric
 *   dS/dz = sech²((|z|−T)/w)/(2w)·sign(z),    dε/dz = (ε_w − ε_m)·dS/dz,
 *                                          dM/dz = −dS/dz.
 *
 * Per-atom energy (kcal/mol):
 *   Non-polar transfer:  U_np,i = ΔG_tr,i · M(z_i),
 *     ΔG_tr,i < 0 hydrophobic (favors core), > 0 polar (penalizes core).
 *     Default ΔG_tr by element: C −0.50, S −0.40 (hydrophobic); N +0.60,
 *     O +0.60, P +1.00 (polar); halogens −0.30; metals +1.50. Matches the
 *     IMM1 water→cyclohexane transfer scale order (Wimley–White-like).
 *   Polar Born cost:     U_pol,i = 166.0317·q_i²/R_i·(1/ε(z_i) − 1/ε_water),
 *     zero in bulk water, positive in the low-dielectric core (Born penalty
 *     for burying charge). R_i = Born/inherent radius (fallback 2.0 Å).
 *   U_mem = Σ_i (U_np,i + U_pol,i).
 *
 * Forces are z-only by slab symmetry:
 *   F_z,i = −ΔG_tr,i·dM/dz − dU_pol/dz,
 *   dU_pol/dz = −166.0317·q_i²/R_i·(1/ε²)·dε/dz.
 * Accumulated into the caller's flat 3n force array (x/y untouched).
 *
 * Units: Å, kcal/mol, e. Float32Array/Float64Array safe. No Node deps.
 *
 * References:
 *   [1] Lazaridis, Proteins 52, 176 (2003) — IMM1 slab + transfer energies.
 *   [2] Im–Feig–Brooks, Biophys J 85, 2900 (2003) — HDGB/GB with ε(z).
 *   [3] Wimley–White interfacial hydrophobicity scale (transfer ordering).
 */

export const MEMBRANE_HALF_THICKNESS = 15.0;
export const MEMBRANE_WIDTH = 2.0;
export const COULOMB_CONST = 332.06371;

/**
 * Default water→membrane transfer free energies (kcal/mol) per element.
 * Negative = hydrophobic (partition into core), positive = polar cost.
 * @type {Record<string, number>}
 */
export const TRANSFER_DG = {
  C: -0.50, S: -0.40, SE: -0.40,
  N: 0.60, O: 0.60, P: 1.00, F: 0.20,
  CL: -0.30, BR: -0.35, I: -0.40,
  ZN: 1.50, FE: 1.50, MG: 1.50, CA: 1.50, CU: 1.50,
  MN: 1.50, NI: 1.50, CO: 1.50, NA: 1.50, K: 1.50,
  DEFAULT: 0.0,
};

/**
 * Transfer ΔG lookup.
 * @param {string|object} a  element string or {element}
 * @returns {number} kcal/mol
 */
export function transferDgFor(a) {
  let e = "DEFAULT";
  if (typeof a === "string") e = a.trim().toUpperCase() || "DEFAULT";
  else if (a && typeof a === "object") e = ((a.element ?? "DEFAULT").toString().trim().toUpperCase() || "DEFAULT");
  return TRANSFER_DG[e] ?? TRANSFER_DG.DEFAULT;
}

/**
 * Water fraction profile S(z) ∈ [0,1].
 * @param {number} z
 * @param {number} [T=15]  half-thickness (Å)
 * @param {number} [w=2]   transition width (Å)
 * @returns {number}
 */
export function waterFraction(z, T = MEMBRANE_HALF_THICKNESS, w = MEMBRANE_WIDTH) {
  const u = (Math.abs(z) - T) / w;
  // tanh clamp avoids overflow for |u| large
  const t = Math.tanh(Math.max(-20, Math.min(20, u)));
  return 0.5 * (1 + t);
}

/**
 * Membrane burial fraction M(z) = 1 − S(z).
 * @param {number} z
 * @param {number} [T]
 * @param {number} [w]
 * @returns {number}
 */
export function membraneFraction(z, T = MEMBRANE_HALF_THICKNESS, w = MEMBRANE_WIDTH) {
  return 1 - waterFraction(z, T, w);
}

/**
 * Depth-dependent dielectric ε(z).
 * @param {number} z
 * @param {number} [epsWater=78.5]
 * @param {number} [epsMem=2.0]
 * @param {number} [T]
 * @param {number} [w]
 * @returns {number}
 */
export function epsilonOfZ(z, epsWater = 78.5, epsMem = 2.0, T = MEMBRANE_HALF_THICKNESS, w = MEMBRANE_WIDTH) {
  return epsMem + (epsWater - epsMem) * waterFraction(z, T, w);
}

/**
 * dS/dz analytic derivative (1/Å).
 * @param {number} z
 * @param {number} [T]
 * @param {number} [w]
 * @returns {number}
 */
export function dWaterFractionDz(z, T = MEMBRANE_HALF_THICKNESS, w = MEMBRANE_WIDTH) {
  const u = (Math.abs(z) - T) / w;
  const uc = Math.max(-20, Math.min(20, u));
  const sech2 = 1 - Math.tanh(uc) * Math.tanh(uc);
  const sgn = z >= 0 ? 1 : -1;
  // At z=0 the profile is flat (|z| cusp smoothed: limit → 0 since sech² small deep in core)
  if (z === 0) return 0;
  return (sech2 / (2 * w)) * sgn;
}

/**
 * Implicit-membrane slab energy + forces.
 *
 * @param {ArrayLike<number>} positions  flat 3n
 * @param {ArrayLike<number>} charges    length n (e)
 * @param {ArrayLike<number>} radii      length n (Å; Born/vdW, fallback 2.0)
 * @param {ArrayLike<number>|Array<string|object>|null} transferOrElements  length n transfer ΔG (kcal/mol) OR element specs (converted via TRANSFER_DG); null → all 0
 * @param {ArrayLike<number>} forces     flat 3n accumulator (added to)
 * @param {object} [opts]
 * @param {number} [opts.thickness=15]  half-thickness T (Å)
 * @param {number} [opts.width=2]  transition width w (Å)
 * @param {number} [opts.epsWater=78.5]
 * @param {number} [opts.epsMem=2.0]
 * @param {number} [opts.zCenter=0]  bilayer midplane z
 * @returns {{energy:number, npEnergy:number, polEnergy:number}}
 */
export function membraneEnergyForces(positions, charges, radii, transferOrElements, forces, opts = {}) {
  const T = opts.thickness ?? MEMBRANE_HALF_THICKNESS;
  const w = opts.width ?? MEMBRANE_WIDTH;
  const epsW = opts.epsWater ?? 78.5;
  const epsM = opts.epsMem ?? 2.0;
  const z0 = opts.zCenter ?? 0;
  const n = charges.length;
  // Resolve per-atom transfer energies
  let dG = null;
  if (transferOrElements && transferOrElements.length >= n) {
    const first = transferOrElements[0];
    if (typeof first === "number") {
      dG = transferOrElements;
    } else {
      dG = new Float64Array(n);
      for (let i = 0; i < n; i++) dG[i] = transferDgFor(/** @type {any} */ (transferOrElements[i]));
    }
  } else {
    dG = new Float64Array(n);
  }
  let npE = 0, polE = 0;
  for (let i = 0; i < n; i++) {
    const z = positions[3 * i + 2] - z0;
    const S = waterFraction(z, T, w);
    const M = 1 - S;
    const dS = dWaterFractionDz(z, T, w);
    const dM = -dS;
    const gtr = dG[i] ?? 0;
    npE += gtr * M;
    let fz = -gtr * dM;
    const q = charges[i];
    if (q !== 0) {
      const R = radii && radii[i] > 0.3 ? radii[i] : 2.0;
      const eps = epsM + (epsW - epsM) * S;
      const Upol = ((COULOMB_CONST / 2) * q * q) / R * (1 / eps - 1 / epsW);
      polE += Upol;
      const deps = (epsW - epsM) * dS;
      const dUpol = (-(COULOMB_CONST / 2) * q * q) / R / (eps * eps) * deps;
      fz += -dUpol;
    }
    forces[3 * i + 2] += fz;
  }
  return { energy: npE + polE, npEnergy: npE, polEnergy: polE };
}
