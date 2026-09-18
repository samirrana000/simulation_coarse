/**
 * placement.js — rigid-body clash-relaxation placement for library and MOL2 ligands.
 *
 * Given a parsed ligand (mol2.js shape) and a target point (or automatic pocket/surface sampling),
 * this module produces a sterically relaxed rigid pose:
 *
 *   1. The molecule is re-centered on its own centroid, given a random rotation,
 *      and translated so the centroid lands on the target point.
 *   2. Steepest descent over the 6 rigid degrees of freedom only (3 translations +
 *      3 exact Rodrigues rotations about the centroid).
 *   3. Intramolecular distances and bond geometries are preserved to float64 roundoff.
 */

import { parseMol2 } from "./mol2.js?v=10";
import { LIG_ELEMENT, LIG_ELEMENT_DEFAULT } from "./ff-params.js?v=10";

// ============================================================================
// Deterministic RNG (mulberry32)
// ============================================================================

function mulberry32(seed) {
  let a = (seed >>> 0) || 1;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ============================================================================
// Library parsing (cached per entry id)
// ============================================================================

const _libCache = new Map();

export function parseLibraryLigand(entry) {
  let mol = _libCache.get(entry.id);
  if (!mol) {
    mol = parseMol2(entry.mol2)[0];
    _libCache.set(entry.id, mol);
  }
  return mol;
}

// ============================================================================
// Small 3×3 rotation helpers (Rodrigues exponential map)
// ============================================================================

function rotExp(ox, oy, oz, R, out) {
  const th = Math.sqrt(ox * ox + oy * oy + oz * oz);
  let kx, ky, kz, A, B;
  if (th < 1e-12) {
    kx = ox; ky = oy; kz = oz; A = 1; B = 0;
  } else {
    kx = ox / th; ky = oy / th; kz = oz / th;
    A = Math.sin(th); B = 1 - Math.cos(th);
  }
  const c = th < 1e-12 ? 1 : Math.cos(th);
  const E0 = c + B * kx * kx, E1 = B * kx * ky - A * kz, E2 = B * kx * kz + A * ky;
  const E3 = B * ky * kx + A * kz, E4 = c + B * ky * ky, E5 = B * ky * kz - A * kx;
  const E6 = B * kz * kx - A * ky, E7 = B * kz * ky + A * kx, E8 = c + B * kz * kz;
  const r0 = R[0], r1 = R[1], r2 = R[2], r3 = R[3], r4 = R[4], r5 = R[5];
  const r6 = R[6], r7 = R[7], r8 = R[8];
  out[0] = E0 * r0 + E1 * r3 + E2 * r6;
  out[1] = E0 * r1 + E1 * r4 + E2 * r7;
  out[2] = E0 * r2 + E1 * r5 + E2 * r8;
  out[3] = E3 * r0 + E4 * r3 + E5 * r6;
  out[4] = E3 * r1 + E4 * r4 + E5 * r7;
  out[5] = E3 * r2 + E4 * r5 + E5 * r8;
  out[6] = E6 * r0 + E7 * r3 + E8 * r6;
  out[7] = E6 * r1 + E7 * r4 + E8 * r7;
  out[8] = E6 * r2 + E7 * r5 + E8 * r8;
}

function clashGrad(rel, ligSig, pp, ps, nColl, cx, cy, cz, margin, out) {
  const m = ligSig.length;
  let U = 0, Fx = 0, Fy = 0, Fz = 0, Tx = 0, Ty = 0, Tz = 0, minRatio = Infinity;
  for (let a = 0; a < m; a++) {
    const a3 = 3 * a;
    const rx = rel[a3], ry = rel[a3 + 1], rz = rel[a3 + 2];
    const ax = rx + cx, ay = ry + cy, az = rz + cz;
    const sL = ligSig[a];
    for (let i = 0; i < nColl; i++) {
      const i3 = 3 * i;
      const dx = ax - pp[i3], dy = ay - pp[i3 + 1], dz = az - pp[i3 + 2];
      const r2 = dx * dx + dy * dy + dz * dz;
      const rE = 0.5 * (ps[i] + sL);
      if (r2 < 1e-6) {
        // rev2-issue5: deterministic index-hashed escape (replaces the fixed
        // +x+y+z diagonal, brittle when +x+y+z is walled in a dense hetero
        // cage). Direction varies per (a,i) but is fixed given indices, so the
        // full placement stays deterministic given seed. The minRatio floor
        // keeps the residual (1/minRatio) finite instead of Infinity.
        if (1e-3 / rE < minRatio) minRatio = 1e-3 / rE;
        const h = ((a + 1) * 73856093 ^ (i + 1) * 19349663 ^ 0x9e3779b9) >>> 0;
        const jx = ((h % 2001) - 1000) / 1000 || 0.317;
        const jy = ((((h / 2001) | 0) % 2001) - 1000) / 1000 || -0.413;
        const jz = ((((h / 4004001) | 0) % 2001) - 1000) / 1000 || 0.521;
        const jl = Math.hypot(jx, jy, jz) || 1;
        const ex = jx / jl, ey = jy / jl, ez = jz / jl;
        Fx += ex; Fy += ey; Fz += ez;
        Tx += ry * ez - rz * ey; Ty += rz * ex - rx * ez; Tz += rx * ey - ry * ex;
        U += 50.0;
        continue;
      }
      const ratio = Math.sqrt(r2) / rE;
      if (ratio < minRatio) minRatio = ratio;
      if (r2 >= margin * margin * rE * rE) continue;
      const r = Math.sqrt(r2);
      const dr = rE - r;
      U += 0.5 * dr * dr;
      const s = dr / r;
      const fx = s * dx, fy = s * dy, fz = s * dz;
      Fx += fx; Fy += fy; Fz += fz;
      Tx += ry * fz - rz * fy;
      Ty += rz * fx - rx * fz;
      Tz += rx * fy - ry * fx;
    }
  }
  out[0] = U; out[1] = Fx; out[2] = Fy; out[3] = Fz;
  out[4] = Tx; out[5] = Ty; out[6] = Tz; out[7] = minRatio;
}

function applyRot(R, ref, rel, m) {
  const r0 = R[0], r1 = R[1], r2 = R[2], r3 = R[3], r4 = R[4], r5 = R[5];
  const r6 = R[6], r7 = R[7], r8 = R[8];
  for (let a = 0; a < m; a++) {
    const a3 = 3 * a;
    const x = ref[a3], y = ref[a3 + 1], z = ref[a3 + 2];
    rel[a3]     = r0 * x + r1 * y + r2 * z;
    rel[a3 + 1] = r3 * x + r4 * y + r5 * z;
    rel[a3 + 2] = r6 * x + r7 * y + r8 * z;
  }
}

/**
 * Relax a rigid ligand pose against the collision set.
 *
 * Collision policy (rev1-issue2): `protein` is the FULL clash set, not
 * necessarily protein-only. Legacy callers pass protein-only arrays
 * ({ pos: 3*nProt, sigma: nProt }) — still supported. Heavy-mode callers
 * pass the extended set (protein + hetero/PDB-ligand/cofactor/metal +
 * existing-ligand atoms, sigma per atom from ff._elem). The incoming
 * ligand's own slot must be excluded by the caller — either by truncating
 * the arrays or via `opts.excludeFrom` / `protein.excludeFrom` (indices
 * [excludeFrom, nTotal) are skipped; defaults to the full arrays).
 */
export function relaxClash(pos, mol, protein, opts = {}) {
  const maxIters = opts.maxIters ?? 350;
  const step0 = opts.step0 ?? 0.06;
  const clashMargin = opts.clashMargin ?? 0.85;

  const m = mol.atoms.length;
  const pp = protein.pos, ps = protein.sigma;
  const nTotal = Math.min(pp.length / 3, ps.length);
  const nColl = Math.max(0, Math.min(opts.excludeFrom ?? protein.excludeFrom ?? nTotal, nTotal));

  const ref = new Float64Array(3 * m);
  const rel = new Float64Array(3 * m);
  const R = new Float64Array(9);
  const Rn = new Float64Array(9);
  const grad = new Float64Array(8);

  const ligSig = new Float64Array(m);
  for (let a = 0; a < m; a++) {
    ligSig[a] = (LIG_ELEMENT[mol.atoms[a].element] || LIG_ELEMENT_DEFAULT).sigma;
  }

  let cx = 0, cy = 0, cz = 0;
  for (let a = 0; a < m; a++) {
    cx += pos[3 * a]; cy += pos[3 * a + 1]; cz += pos[3 * a + 2];
  }
  cx /= m; cy /= m; cz /= m;
  for (let a = 0; a < m; a++) {
    ref[3 * a] = pos[3 * a] - cx;
    ref[3 * a + 1] = pos[3 * a + 1] - cy;
    ref[3 * a + 2] = pos[3 * a + 2] - cz;
  }
  R[0] = R[4] = R[8] = 1;
  rel.set(ref);

  clashGrad(rel, ligSig, pp, ps, nColl, cx, cy, cz, clashMargin, grad);
  let E = grad[0], minRatio = grad[7];
  let iter = 0, converged = minRatio >= clashMargin;

  if (!converged) {
    let step = step0;
    for (iter = 1; iter <= maxIters; iter++) {
      const Fx = grad[1], Fy = grad[2], Fz = grad[3];
      const Tx = grad[4], Ty = grad[5], Tz = grad[6];
      if (Fx === 0 && Fy === 0 && Fz === 0 && Tx === 0 && Ty === 0 && Tz === 0) {
        iter--; break;
      }

      rotExp(step * Tx, step * Ty, step * Tz, R, Rn);
      applyRot(Rn, ref, rel, m);
      const ncx = cx + step * Fx, ncy = cy + step * Fy, ncz = cz + step * Fz;
      clashGrad(rel, ligSig, pp, ps, nColl, ncx, ncy, ncz, clashMargin, grad);

      if (grad[0] < E) {
        E = grad[0]; minRatio = grad[7];
        R.set(Rn); cx = ncx; cy = ncy; cz = ncz;
        if (minRatio >= clashMargin) { converged = true; break; }
      } else {
        applyRot(R, ref, rel, m);
        clashGrad(rel, ligSig, pp, ps, nColl, cx, cy, cz, clashMargin, grad);
        step *= 0.5;
        if (step < 1e-10) break;
      }
    }
  }

  for (let a = 0; a < m; a++) {
    pos[3 * a] = rel[3 * a] + cx;
    pos[3 * a + 1] = rel[3 * a + 1] + cy;
    pos[3 * a + 2] = rel[3 * a + 2] + cz;
  }

  const residualClash = minRatio > 0 ? 1 / minRatio : Infinity;
  return {
    pos,
    iterations: Math.min(iter, maxIters),
    residualClash,
    converged,
  };
}

/**
 * Place a ligand at a world-space target point with clash avoidance.
 *
 * `opts.protein` is the full collision set (see relaxClash): protein-only
 * for legacy callers, or the extended protein+hetero/existing-ligand set
 * built by ligand-panel getProteinCoordsAndSigma in heavy mode. The
 * incoming ligand's own slot must already be excluded (or passed via
 * `opts.excludeFrom`, forwarded to relaxClash).
 */
export function placeLigand(mol, target, opts = {}) {
  const m = mol.atoms.length;
  const seed = opts.seed ?? 1;

  let cx = 0, cy = 0, cz = 0;
  for (const at of mol.atoms) { cx += at.x; cy += at.y; cz += at.z; }
  cx /= m; cy /= m; cz /= m;

  const rand = mulberry32(seed);
  // Isotropic uniform random rotation on SO(3) via Shoemake quaternion subgroup algorithm
  const u1 = rand(), u2 = rand(), u3 = rand();
  const q0 = Math.sqrt(1 - u1) * Math.sin(2 * Math.PI * u2);
  const q1 = Math.sqrt(1 - u1) * Math.cos(2 * Math.PI * u2);
  const q2 = Math.sqrt(u1) * Math.sin(2 * Math.PI * u3);
  const q3 = Math.sqrt(u1) * Math.cos(2 * Math.PI * u3);

  const r00 = 1 - 2 * (q2 * q2 + q3 * q3);
  const r01 = 2 * (q1 * q2 - q0 * q3);
  const r02 = 2 * (q1 * q3 + q0 * q2);
  const r10 = 2 * (q1 * q2 + q0 * q3);
  const r11 = 1 - 2 * (q1 * q1 + q3 * q3);
  const r12 = 2 * (q2 * q3 - q0 * q1);
  const r20 = 2 * (q1 * q3 - q0 * q2);
  const r21 = 2 * (q2 * q3 + q0 * q1);
  const r22 = 1 - 2 * (q1 * q1 + q2 * q2);

  const pos = new Float64Array(3 * m);
  for (let a = 0; a < m; a++) {
    const x = mol.atoms[a].x - cx, y = mol.atoms[a].y - cy, z = mol.atoms[a].z - cz;
    pos[3 * a]     = r00 * x + r01 * y + r02 * z + target[0];
    pos[3 * a + 1] = r10 * x + r11 * y + r12 * z + target[1];
    pos[3 * a + 2] = r20 * x + r21 * y + r22 * z + target[2];
  }

  if (!opts.protein) {
    return { pos, iterations: 0, residualClash: 1.0, converged: true };
  }

  return relaxClash(pos, mol, opts.protein, opts);
}

/**
 * Automatically detect a prominent concave pocket / binding cavity center.
 *
 * @param {Float64Array} protPos coordinates (3 * nProt or longer)
 * @param {number} [nProt] atoms to scan. Defaults to the full array
 *   (protPos.length / 3). Pass the protein-only count to restrict the
 *   cavity search to protein beads while clashing against a longer
 *   protein+hetero collision array (ligand-panel placeInPocket policy).
 * @returns {[number, number, number]} pocket center coordinates
 */
export function findPocketCenter(protPos, nProt = protPos.length / 3) {
  if (nProt === 0) return [0, 0, 0];

  // Calculate geometric center of mass (COM) and maximum radius
  let comX = 0, comY = 0, comZ = 0;
  for (let i = 0; i < nProt; i++) {
    comX += protPos[3 * i];
    comY += protPos[3 * i + 1];
    comZ += protPos[3 * i + 2];
  }
  comX /= nProt; comY /= nProt; comZ /= nProt;

  let maxR2 = 0;
  for (let i = 0; i < nProt; i++) {
    const dx = protPos[3 * i] - comX, dy = protPos[3 * i + 1] - comY, dz = protPos[3 * i + 2] - comZ;
    maxR2 = Math.max(maxR2, dx * dx + dy * dy + dz * dz);
  }
  const maxR = Math.sqrt(maxR2) || 1.0;

  // Search for surface/pocket residues (exclude buried core R < 0.45 * maxR)
  let bestScore = -Infinity;
  let bestPoint = [comX + 0.6 * maxR, comY, comZ];

  for (let i = 0; i < nProt; i += 2) {
    const xi = protPos[3 * i], yi = protPos[3 * i + 1], zi = protPos[3 * i + 2];
    const rCOM = Math.hypot(xi - comX, yi - comY, zi - comZ);

    // Skip deeply buried interior atoms
    if (rCOM < 0.45 * maxR) continue;

    // Measure local pocket concavity: neighbor shell density at 5-11 Å
    let nClose = 0;
    let nDirectClash = 0;
    let localComX = 0, localComY = 0, localComZ = 0;

    for (let j = 0; j < nProt; j += 2) {
      if (i === j) continue;
      const dx = protPos[3 * j] - xi, dy = protPos[3 * j + 1] - yi, dz = protPos[3 * j + 2] - zi;
      const d2 = dx * dx + dy * dy + dz * dz;
      if (d2 < 12.25) nDirectClash++; // d < 3.5 Å
      if (d2 < 121.0) { // d < 11.0 Å
        nClose++;
        localComX += protPos[3 * j];
        localComY += protPos[3 * j + 1];
        localComZ += protPos[3 * j + 2];
      }
    }

    if (nClose >= 6) {
      localComX /= nClose;
      localComY /= nClose;
      localComZ /= nClose;

      // Pocket score: high surrounding wall density (nClose) with open cavity access
      const cavityScore = nClose / (1.0 + 0.5 * nDirectClash);
      if (cavityScore > bestScore) {
        bestScore = cavityScore;
        // Offset candidate pocket center into cavity space away from direct atom collision
        const outDx = xi - localComX, outDy = yi - localComY, outDz = zi - localComZ;
        const outLen = Math.hypot(outDx, outDy, outDz) || 1.0;
        const offsetDist = 4.0; // Å outward probe offset
        bestPoint = [
          xi + (outDx / outLen) * offsetDist,
          yi + (outDy / outLen) * offsetDist,
          zi + (outDz / outLen) * offsetDist,
        ];
      }
    }
  }

  return bestPoint;
}
