/**
 * sasa.js — Non-polar hydrophobic solvation energy via LCPO / neighbor burial.
 *
 * Implements the LCPO (Linear Combinations of Pairwise Overlaps) model:
 *   SASA_i = P1 * 4*pi*R_i^2 - P2 * sum_j A_ij + P3 * sum_j sum_k A_jk ...
 *
 * Non-polar solvation energy:
 *   U_np = gamma_surf * sum_i SASA_i + sum_i b_i
 *
 * Typically gamma_surf = 0.005 to 0.015 kcal/(mol·Å²).
 * Hydrophobic burial drives non-polar ligand binding into protein hydrophobic pockets.
 */

import { SASA_RADII } from "./charges.js?v=10";

// Probe radius of water (Å)
export const PROBE_RADIUS = 1.4;

export class SasaModel {
  /**
   * @param {object} params
   * @param {number} [params.gamma=0.0072] Surface tension coefficient in kcal/(mol·Å²)
   * @param {number} [params.probeRadius=1.4] Solvent probe radius in Å
   */
  constructor(params = {}) {
    this.gamma = params.gamma ?? 0.0072;
    this.probeRadius = params.probeRadius ?? PROBE_RADIUS;
  }

  /**
   * Compute pairwise burial / overlap non-polar forces and energy.
   *
   * @param {Float64Array} pos
   * @param {Float64Array} forces
   * @param {Array<{element: string}>} atoms
   * @param {number} n
   * @param {number} [nProt]
   * @param {number} [ligStart]
   * @returns {{energy: number, sasa: number, bindSasaE: number}}
   */
  compute(pos, forces, atoms, n, nProt = 0, ligStart = 0) {
    let totalEnergy = 0;
    let bindEnergy = 0;
    let totalSasa = 0;

    // 1. Initial isolated surface areas
    const s0 = new Float64Array(n);
    const radii = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const el_i = atoms[i]?.element || "C";
      radii[i] = (SASA_RADII[el_i] ?? SASA_RADII.DEFAULT) + this.probeRadius;
      s0[i] = 4.0 * Math.PI * radii[i] * radii[i];
    }
    const burial = new Float64Array(n);

    // 2. Symmetric pairwise overlap evaluation (i < j)
    for (let i = 0; i < n; i++) {
      const xi = 3 * i;
      const Ri = radii[i];
      for (let j = i + 1; j < n; j++) {
        const xj = 3 * j;
        const Rj = radii[j];
        const cut = Ri + Rj;
        const dx = pos[xj] - pos[xi];
        const dy = pos[xj + 1] - pos[xi + 1];
        const dz = pos[xj + 2] - pos[xi + 2];
        const r2 = dx * dx + dy * dy + dz * dz;

        if (r2 < cut * cut && r2 > 1e-4) {
          const r = Math.sqrt(r2);
          // Overlap sphere cap area for pair (i, j)
          const overlap_ij = Math.max(0, Math.PI * Ri * (cut - r) * (1.0 + (Rj - Ri) / r));
          const overlap_ji = Math.max(0, Math.PI * Rj * (cut - r) * (1.0 + (Ri - Rj) / r));
          
          burial[i] += overlap_ij * 0.25;
          burial[j] += overlap_ji * 0.25;

          // Derivative of total burial (overlap_ij+overlap_ji) w.r.t. r:
          // S = PI*(cut - r)*(C - D/r) where C=Ri+Rj, D=(Ri-Rj)^2
          // dS/dr = -PI*C + PI*C*D/r^2 = -PI*C*(1 - D/r^2)
          const C = Ri + Rj;
          const D = (Ri - Rj) * (Ri - Rj);
          const dOverlap_dr = -Math.PI * C * (1.0 - D / (r * r));
          const fMag = -this.gamma * 0.25 * dOverlap_dr / r;

          // Newton's third law: F_i = -F_j
          forces[xi] += fMag * dx;
          forces[xi + 1] += fMag * dy;
          forces[xi + 2] += fMag * dz;

          forces[xj] -= fMag * dx;
          forces[xj + 1] -= fMag * dy;
          forces[xj + 2] -= fMag * dz;

          if (i < nProt && j >= ligStart && ligStart > 0) {
            bindEnergy += -this.gamma * 0.25 * (overlap_ij + overlap_ji);
          }
        }
      }
    }

    for (let i = 0; i < n; i++) {
      // Unclamped for force-energy consistency (analytic forces = -dE/dr).
      // Clamping to max(0, s0 - burial) would make E non-differentiable and
      // break finite-difference validation (SASA would be 0 while forces non-zero).
      // Keep unclamped here; totalSasa reported as physical (clamped) for display.
      const sasa_i = s0[i] - burial[i];
      const sasa_phys = Math.max(0, sasa_i);
      totalSasa += sasa_phys;
      totalEnergy += this.gamma * sasa_i;
    }

    return { energy: totalEnergy, sasa: totalSasa, bindSasaE: bindEnergy };
  }
}
