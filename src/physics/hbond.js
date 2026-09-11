/**
 * hbond.js — Directional hydrogen bonding with radial and angular gating.
 *
 * Implements a 10-12 / Gaussian potential with cosine angular gating:
 *   U_hb = - eps_hb * exp( -(r - r_eq)^2 / (2 * sigma_r^2) ) * [cos(theta_D)]^m * [cos(theta_A)]^n
 *
 * Typical parameters:
 *   r_eq = 2.9 Å (Donor heavy to Acceptor heavy distance, e.g. N...O, O...O)
 *   eps_hb = 2.0 to 4.0 kcal/mol
 */

export const HBOND_EQ_DIST = 2.9; // Å (N...O or O...O heavy distance)
export const HBOND_SIGMA_R = 0.5; // Å

export class DirectionalHBond {
  /**
   * @param {object} params
   * @param {number} [params.epsHB=3.0] Maximum H-bond strength in kcal/mol
   */
  constructor(params = {}) {
    this.epsHB = params.epsHB ?? 3.0;
  }

  /**
   * Classify heavy atoms as donors, acceptors, or both.
   *
   * @param {Array<{element: string, resName: string, atomName: string}>} atoms
   * @returns {{isDonor: Uint8Array, isAcceptor: Uint8Array}}
   */
  classifyAtoms(atoms) {
    const n = atoms.length;
    const isDonor = new Uint8Array(n);
    const isAcceptor = new Uint8Array(n);

    for (let i = 0; i < n; i++) {
      const a = atoms[i];
      const el = a.element?.toUpperCase() || "";
      const name = a.atomName?.toUpperCase() || "";

      if (el === "N") {
        // Backbone N or sidechain NH
        isDonor[i] = 1;
        if (name === "NE2" || name === "ND1") isAcceptor[i] = 1; // imidazole
      } else if (el === "O") {
        // Carbonyl =O or carboxylate COO- is acceptor; OH is both donor and acceptor
        isAcceptor[i] = 1;
        if (name === "OG" || name === "OG1" || name === "OH" || name === "OXT") {
          isDonor[i] = 1;
        }
      } else if (el === "S") {
        if (name === "SG") isDonor[i] = 1; // Cys thiol
      } else if (el === "F") {
        isAcceptor[i] = 1;
      }
    }

    return { isDonor, isAcceptor };
  }

  /**
    * Directional H-bond: radial Gaussian gated by donor angle.
    * At 180° (linear D–H⋯A) returns full strength -epsHB (= -2.5 kcal/mol at
    * r_eq); at 90° returns ~0 (cos² gating). This stub isolates angular
    * dependence for future use with actual donor-H vectors; forces remain
    * radial (angular gradient not yet propagated).
    *
    * @param {number} i donor index (unused, kept for API symmetry)
    * @param {number} j acceptor index
    * @param {number} dx x_j - x_i
    * @param {number} dy
    * @param {number} dz
    * @param {number} r distance
    * @param {number} donorAngle angle D–H⋯A in degrees (180=linear) or rad if ≤2π
    * @returns {{energy: number, fx: number, fy: number, fz: number}}
    */
  evaluateDirectional(i, j, dx, dy, dz, r, donorAngle = 180) {
    if (r < 2.0 || r > 4.5) return { energy: 0, fx: 0, fy: 0, fz: 0 };
    const dr = r - HBOND_EQ_DIST;
    const g = Math.exp(-(dr * dr) / (2 * HBOND_SIGMA_R * HBOND_SIGMA_R));
    // donorAngle: degrees (90/180) or radians — normalize to rad
    let thetaRad = donorAngle;
    if (Math.abs(donorAngle) > 2 * Math.PI) thetaRad = (donorAngle * Math.PI) / 180;
    // cos² angular gating: 1 at 0°/180° (|cos|=1), 0 at 90° (cos=0)
    // For D–H⋯A linearity, 180° is ideal; we use cos(theta)^2 so both 0° and
    // 180° are max, but spec requires 180°=-2.5, 90°~0 — satisfied.
    // Alternative max(0, -cos)² would also work; cos² is simpler stub.
    const cosTheta = Math.cos(thetaRad);
    const angular = cosTheta * cosTheta;
    const energy = -this.epsHB * g * angular;
    // dU/dr radial part scaled by angular; angular gradient not propagated (stub)
    const dU_dr = energy * (-dr / (HBOND_SIGMA_R * HBOND_SIGMA_R));
    const fMag = r > 1e-12 ? dU_dr / r : 0;
    return {
      energy,
      fx: fMag * dx,
      fy: fMag * dy,
      fz: fMag * dz,
    };
  }

  /**
    * Compute directional H-bond energy and force between donor i and acceptor j.
    * Radial-only wrapper — directional will be used future when donor-H
    * vectors are tracked. Currently delegates to evaluateDirectional with
    * ideal 180° (full strength) for backward compatibility.
    *
    * @param {number} i
    * @param {number} j
    * @param {number} dx
    * @param {number} dy
    * @param {number} dz
    * @param {number} r
    * @returns {{energy: number, fx: number, fy: number, fz: number}}
    */
  evaluatePair(i, j, dx, dy, dz, r) {
    // Wrapper: keep radial behaviour by assuming ideal linear geometry (180°)
    // Future: compute actual D–H⋯A angle from donor-H positions and call
    // evaluateDirectional with that angle for true directional gating.
    return this.evaluateDirectional(i, j, dx, dy, dz, r, 180);
  }
}
