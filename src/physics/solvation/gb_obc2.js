/**
 * gb_obc2.js — Generalized Born OBC-II (Onufriev–Bashford–Case) implicit solvent.
 *
 * Implements:
 *   HCT pairwise descreening (Hawkins–Cramer–Truhlar) +
 *   OBC-II rescaling (α=1.0, β=0.8, γ=4.85, offset 0.09 Å) +
 *   Still equation  f_GB = sqrt(r² + Ri·Rj·exp(−r²/4RiRj)) +
 *   salt-screened reaction field  ΔG_GB (Debye–Hückel κ from ionic strength).
 *
 * Energy (kcal/mol, COULOMB_CONST = 332.06371 kcal·Å/mol/e²):
 *   Self:   G_self,i = −166.0317·q_i²/R_i·(1/ε_in − e^(−κR_i)/ε_out)
 *   Pair:   G_pair   = −332.0637·(q_i q_j/f_GB)·(1/ε_in − e^(−κf)/ε_out)
 *           summed once over i<j (the 1/2 double-sum folded in).
 *   Total ΔG_GB = Σ_i G_self,i + Σ_{i<j} G_pair,ij.
 * Coulomb in solute dielectric is NOT included here (caller adds it or uses
 * pairCoulomb); this module is the pure GB reaction field + optional screened
 * Coulomb helper. gbEnergyForces returns gbEnergy (+ optionally coulomb).
 *
 * Forces: analytic, fixed-Radii explicit term plus per-atom dG/dR outputs for
 * the outer Born-radius chain rule:
 *   dU/dr  via df_GB/dr = r(1 − e/4)/f,  e = exp(−r²/4RiRj),
 *           dP/df with P(f) = 1/ε_in − e^(−κf)/ε_out, dP/df = κe^(−κf)/ε_out,
 *           dg/df with g = P/f, dg/df = (f·dP/df − P)/f²,
 *           dU/dr = −C·qq·dg/df·df/dr   (C = 332.0637 per pair).
 *   dU/dR_i via df/dR_i = Rj·e·(1 + r²/4a)/(2f), a = Ri·Rj, plus self-term
 *           dG_self/dR analytic. Returned as dGdR[i] so callers doing
 *           dR/dx propagation (full OBC forces) can contract it.
 *
 * Units: Å, kcal/mol, e. Works on Float32Array or Float64Array positions.
 * No Node deps, `?v=` compatible.
 *
 * References:
 *   [1] Still et al., JACS 112, 6127 (1990) — f_GB form.
 *   [2] Hawkins–Cramer–Truhlar, JPCB 100, 19824 (1996) — pairwise descreening.
 *   [3] Onufriev–Bashford–Case, Proteins 55, 383 (2004) — OBC α/β/γ + offset.
 *   [4] Onufriev et al., JPCB 104, 3712 (2000) — GB salt screening κ.
 */

export const OBC_ALPHA = 1.0;
export const OBC_BETA = 0.8;
export const OBC_GAMMA = 4.85;
export const OBC_OFFSET = 0.09;
export const HCT_SCALE = 0.8;
export const COULOMB_CONST = 332.06371;
export const R_MAX = 30.0;

/**
 * Inverse Debye length κ (Å⁻¹) from ionic strength (Debye–Hückel).
 * κ ≈ 0.329·sqrt(I·(298.15/T)·(78.5/ε_out)).
 * @param {number} saltM  ionic strength in mol/L (M)
 * @param {number} tempK  temperature in K
 * @param {number} epsOut solvent dielectric
 * @returns {number} κ in Å⁻¹
 */
export function debyeKappa(saltM = 0.15, tempK = 300, epsOut = 78.5) {
  if (!(saltM > 0)) return 0;
  const T = Math.max(10, tempK);
  const eps = Math.max(1, epsOut);
  return 0.329 * Math.sqrt(saltM * (298.15 / T) * (78.5 / eps));
}

/**
 * Still GB function f_GB(r, Ri, Rj).
 * @param {number} r
 * @param {number} Ri
 * @param {number} Rj
 * @returns {number}
 */
export function stillF(r, Ri, Rj) {
  const a = Math.max(1e-6, Ri * Rj);
  const r2 = r * r;
  return Math.sqrt(r2 + a * Math.exp(-r2 / (4 * a)));
}

/**
 * Derivative df_GB/dr at (r, Ri, Rj).
 * @param {number} r
 * @param {number} Ri
 * @param {number} Rj
 * @param {number} [f]  precomputed f_GB (recomputed if omitted)
 * @returns {number}
 */
export function dStillFdr(r, Ri, Rj, f) {
  const a = Math.max(1e-6, Ri * Rj);
  const r2 = r * r;
  const e = Math.exp(-r2 / (4 * a));
  const ff = f ?? Math.sqrt(r2 + a * e);
  if (ff < 1e-12) return 0;
  return (r * (1 - 0.25 * e)) / ff;
}

/**
 * Derivative df_GB/dRi at (r, Ri, Rj).
 * df/da = e(1 + r²/4a)/(2f), df/dRi = df/da·Rj.
 * @param {number} r
 * @param {number} Ri
 * @param {number} Rj
 * @param {number} [f] precomputed f_GB
 * @returns {number}
 */
export function dStillFdRi(r, Ri, Rj, f) {
  const a = Math.max(1e-6, Ri * Rj);
  const r2 = r * r;
  const e = Math.exp(-r2 / (4 * a));
  const ff = f ?? Math.sqrt(r2 + a * e);
  if (ff < 1e-12) return 0;
  const dfda = (e * (1 + r2 / (4 * a))) / (2 * ff);
  return dfda * Rj;
}

/**
 * HCT descreening sum S_i = Σ_{j≠i} I_ij for one atom (internal helper).
 * I_ij: far  r > ρ_i+ρ_j : V_j/(3 r^4)... near: overlap integral branch.
 * Mirrors src/physics/gb.js HCT for consistency (cutoff 12 Å).
 * @param {number} i
 * @param {ArrayLike<number>} pos  flat 3n positions
 * @param {ArrayLike<number>} rho  intrinsic radii (offset-subtracted outside)
 * @param {number} n
 * @param {number} cutoff
 * @returns {number}
 */
function hctSum(i, pos, rho, n, cutoff = 12) {
  const xi = 3 * i;
  const rhoI = rho[i];
  const cut2 = cutoff * cutoff;
  let sum = 0;
  for (let j = 0; j < n; j++) {
    if (j === i) continue;
    const xj = 3 * j;
    const rhoJ = rho[j] * HCT_SCALE;
    const dx = pos[xj] - pos[xi];
    const dy = pos[xj + 1] - pos[xi + 1];
    const dz = pos[xj + 2] - pos[xi + 2];
    const r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < 1e-4 || r2 >= cut2) continue;
    const r = Math.sqrt(r2);
    if (r > rhoI + rhoJ) {
      sum += (rhoJ * rhoJ * rhoJ) / (3 * r2 * r2);
    } else if (r > Math.abs(rhoI - rhoJ) && r > 1e-6) {
      const lij = Math.max(rhoI, r - rhoJ);
      const uij = r + rhoJ;
      if (lij > 1e-9 && uij > lij) {
        sum += (1 / (4 * r)) * (Math.log(uij / lij) + (rhoJ * rhoJ - (r - lij) * (r - lij)) / (2 * lij * lij));
      }
    }
  }
  return sum;
}

/**
 * OBC-II effective Born radii from positions + intrinsic radii.
 *
 * Steps per atom i:
 *   ρ̃_i = max(0.5, ρ_i − 0.09)            (OBC offset radii)
 *   ψ_i  = ρ̃_i · Σ_j I_ij                  (dimensionless burial, HCT sums on ρ̃)
 *   t    = αψ − βψ² + γψ³
 *   R_i⁻¹ = ρ̃_i⁻¹ − tanh(t)/ρ_i  →  R_i = 1/R_i⁻¹, clamped [ρ_i, 30 Å].
 * For ψ→0, tanh(t)≈ψ so R⁻¹ ≈ 1/ρ̃ − ψ/ρ ≈ 1/ρ̃ − ΣI (HCT limit). For deep
 * burial ψ≫1, tanh→1, R⁻¹→(1/ρ̃ − 1/ρ) small → R large (correct screening).
 * Isolated atom (ΣI=0) returns ρ̃ (≈ ρ−0.09) — intrinsic limit.
 *
 * @param {ArrayLike<number>} positions  flat 3n (Float32Array|Float64Array)
 * @param {ArrayLike<number>} intrinsicRadii  length n (Å, e.g. GB_RADII)
 * @param {object} [opts]
 * @param {number} [opts.cutoff=12]  HCT pair cutoff (Å)
 * @param {number} [opts.offset=0.09] OBC radius offset (Å)
 * @param {number} [opts.alpha=1.0]
 * @param {number} [opts.beta=0.8]
 * @param {number} [opts.gamma=4.85]
 * @param {number} [opts.rMax=30]  cap for buried radii
 * @returns {Float64Array} effective Born radii R_i (Å)
 */
export function computeBornRadii(positions, intrinsicRadii, opts = {}) {
  const cutoff = opts.cutoff ?? 12;
  const offset = opts.offset ?? OBC_OFFSET;
  const alpha = opts.alpha ?? OBC_ALPHA;
  const beta = opts.beta ?? OBC_BETA;
  const gamma = opts.gamma ?? OBC_GAMMA;
  const rMax = opts.rMax ?? R_MAX;
  const n = intrinsicRadii.length;
  const rho = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const r0 = intrinsicRadii[i] > 0 ? intrinsicRadii[i] : 1.6;
    rho[i] = Math.max(0.5, r0 - offset);
  }
  const out = new Float64Array(n);
  if (!positions || positions.length < 3 * n || n === 0) {
    for (let i = 0; i < n; i++) out[i] = rho[i];
    return out;
  }
  for (let i = 0; i < n; i++) {
    const sumI = hctSum(i, positions, rho, n, cutoff);
    const psi = rho[i] * sumI;
    const t = alpha * psi - beta * psi * psi + gamma * psi * psi * psi;
    const th = Math.tanh(Math.max(-10, Math.min(10, t)));
    const rhoOrig = intrinsicRadii[i] > 0 ? intrinsicRadii[i] : 1.6;
    const invR = 1 / rho[i] - th / rhoOrig;
    // Floor is the offset radius rho_tilde (OBC convention); isolated atoms
    // return rho_tilde = rho_orig − 0.09. Cap buried radii at rMax.
    let R = invR > 1e-6 ? 1 / invR : rMax;
    if (R < rho[i]) R = rho[i];
    if (R > rMax) R = rMax;
    if (!Number.isFinite(R)) R = rMax;
    out[i] = R;
  }
  return out;
}

/**
 * GB-OBC2 reaction-field energy + analytic forces (fixed-Radii explicit part).
 *
 * Accumulates into `forces` (flat 3n, same layout as positions; added to, not
 * zeroed). Optionally fills `dGdR` (length n) with ∂G/∂R_i (self + pair parts)
 * for callers propagating Born-radius derivatives.
 *
 * @param {ArrayLike<number>} positions  flat 3n
 * @param {ArrayLike<number>} charges    length n (e)
 * @param {ArrayLike<number>} bornRadii  length n (Å, from computeBornRadii)
 * @param {ArrayLike<number>} forces     flat 3n accumulator (Float64/Float32)
 * @param {object} [opts]
 * @param {number} [opts.epsIn=4.0]
 * @param {number} [opts.epsOut=78.5]
 * @param {number} [opts.saltM=0.15]  used when opts.kappa omitted
 * @param {number} [opts.tempK=300]
 * @param {number} [opts.kappa]  direct inverse Debye length (Å⁻¹, overrides saltM)
 * @param {number} [opts.cutoff=1e9]  pair cutoff (Å); default no cutoff
 * @param {Set<number>|null} [opts.excluded]  numeric pair keys i*1e6+j to skip
 * @param {Map<number,number>|null} [opts.scale14]  pair key -> Coulomb/GB scale
 * @param {Float64Array|null} [opts.dGdR]  length-n output for ∂G/∂R_i
 * @param {boolean} [opts.includeSelf=true]  add Born self terms
 * @param {boolean} [opts.includeCoulomb=false] add 332·qq/(ε_in·r) + forces
 * @returns {{energy:number, gbEnergy:number, coulEnergy:number, selfEnergy:number, pairEnergy:number, dGdR:Float64Array|null}}
 */
export function gbEnergyForces(positions, charges, bornRadii, forces, opts = {}) {
  const n = charges.length;
  const epsIn = opts.epsIn ?? 4.0;
  const epsOut = opts.epsOut ?? 78.5;
  const kappa = opts.kappa ?? debyeKappa(opts.saltM ?? 0.15, opts.tempK ?? 300, epsOut);
  const cutoff = opts.cutoff ?? 1e9;
  const cut2 = cutoff * cutoff;
  const excluded = opts.excluded ?? null;
  const scale14 = opts.scale14 ?? null;
  const includeSelf = opts.includeSelf ?? true;
  const includeCoulomb = opts.includeCoulomb ?? false;
  /** @type {Float64Array|null} */
  let dGdR = opts.dGdR ?? null;
  const ownDgdR = !dGdR;
  if (ownDgdR) dGdR = new Float64Array(n);

  const C = COULOMB_CONST;
  const invIn = 1 / epsIn, invOut = 1 / epsOut;
  let selfE = 0, pairE = 0, coulE = 0;

  if (includeSelf) {
    for (let i = 0; i < n; i++) {
      const q = charges[i];
      if (q === 0) continue;
      const R = Math.max(0.5, bornRadii[i]);
      const eK = Math.exp(-kappa * R);
      const P = invIn - eK * invOut;
      const g = P / R;
      selfE += -0.5 * C * q * q * g;
      // dG_self/dR = −0.5 C q² [(κe^(−κR)/ε_out·R − P)/R²]
      const dP = (kappa * eK * invOut);
      const dg = (dP * R - P) / (R * R);
      dGdR[i] += -0.5 * C * q * q * dg;
    }
  }

  for (let i = 0; i < n; i++) {
    const qi = charges[i];
    if (qi === 0) continue;
    const xi = 3 * i;
    const Ri = Math.max(0.5, bornRadii[i]);
    for (let j = i + 1; j < n; j++) {
      const qj = charges[j];
      if (qj === 0) continue;
      if (excluded) {
        const key = i < j ? i * 1e6 + j : j * 1e6 + i;
        if (excluded.has(key)) continue;
      }
      let s14 = 1;
      if (scale14) {
        const key = i < j ? i * 1e6 + j : j * 1e6 + i;
        const s = scale14.get(key);
        if (s !== undefined) s14 = s;
      }
      const xj = 3 * j;
      const dx = positions[xj] - positions[xi];
      const dy = positions[xj + 1] - positions[xi + 1];
      const dz = positions[xj + 2] - positions[xi + 2];
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 >= cut2 || r2 < 1e-12) continue;
      const r = Math.sqrt(r2);
      const Rj = Math.max(0.5, bornRadii[j]);
      const qq = qi * qj * s14;
      const a = Math.max(1e-6, Ri * Rj);
      const e = Math.exp(-r2 / (4 * a));
      const f2 = r2 + a * e;
      const f = Math.sqrt(f2);
      const eK = Math.exp(-kappa * f);
      const P = invIn - eK * invOut;
      const g = P / f;
      pairE += -C * qq * g;

      // dU/dr chain through f
      const dfdr = (r * (1 - 0.25 * e)) / f;
      const dP = kappa * eK * invOut;
      const dgdf = (dP * f - P) / f2;
      const dudr = -C * qq * dgdf * dfdr;
      // optional screened Coulomb handled by caller; here GB-only unless asked
      let dudrTot = dudr;
      if (includeCoulomb) {
        const uC = ((C / epsIn) * qq) / r;
        coulE += uC;
        dudrTot += -uC / r;
      }
      const fmag = dudrTot / r;
      const fx = fmag * dx, fy = fmag * dy, fz = fmag * dz;
      forces[xi] += fx; forces[xi + 1] += fy; forces[xi + 2] += fz;
      forces[xj] -= fx; forces[xj + 1] -= fy; forces[xj + 2] -= fz;

      // ∂G/∂R pieces (for outer dR/dx chain rule)
      const dfdRi = (e * (1 + r2 / (4 * a))) / (2 * f) * Rj;
      const dfdRj = (e * (1 + r2 / (4 * a))) / (2 * f) * Ri;
      const dgdR = dgdf; // dg/df shared
      dGdR[i] += -C * qq * dgdR * dfdRi;
      dGdR[j] += -C * qq * dgdR * dfdRj;
    }
  }

  const gbEnergy = selfE + pairE;
  return {
    energy: gbEnergy + (includeCoulomb ? coulE : 0),
    gbEnergy,
    coulEnergy: coulE,
    selfEnergy: selfE,
    pairEnergy: pairE,
    dGdR: ownDgdR ? dGdR : dGdR,
  };
}
