/**
 * tirion_anm.js — Tirion distance-weighted Cα elastic network + secondary-
 * structure-dependent backbone pseudo-dihedral basins.
 *
 * ENM springs (Tirion 1996 single-parameter form, distance-weighted variant):
 *   U_ENM = Σ_{i<j, r0_ij ≤ Rc, |i−j|>2} 1/2·γ_ij·(r_ij − r0_ij)²,
 *   γ_ij  = γ0·(R0 / r0_ij)^6.
 * The (R0/r0)^6 falloff (≈ Lennard-Jones attractive-tail curvature) keeps
 * sequential neighbors stiff while softening long 8–10 Å contacts, giving a
 * broader, more physical normal-mode spectrum than uniform-γ ANM while
 * preserving the native fold. R0 ≈ 3.81 Å (Cα–Cα virtual-bond length) is the
 * reference; γ0 sets the overall scale (kcal/mol/Å²).
 *
 * Backbone pseudo-dihedrals (CG dihedral basins from local geometry):
 * for each i..i+3 quadruple inside one contiguous segment, the native
 * pseudo-dihedral φ0 (dihedral of four successive Cα) and pseudo-angle θ0
 * classify local secondary structure:
 *   helix: r_{i,i+3} < 6.5 Å and |φ0| < 80°      → stiff basin k=2.0
 *   sheet: r_{i,i+3} > 8.0 Å and θ_mid > 110°   → medium basin k=1.5
 *   loop/coil: otherwise                        → soft basin k=0.5
 * Thresholds from Cα geometry statistics (helix i..i+3 ≈ 5.5 Å, sheet ≈ 9.5 Å;
 * helix φ ≈ +50°, sheet φ ≈ ±180°). Basins are harmonic in wrapped angle:
 *   U_dih = Σ 1/2·k_ss·wrap(φ − φ0)².
 * Callers without a CG dihedral kernel can consume pseudoDihedrals as restraints
 * via ff-harmonic dihedralForcesAnalytic (stride-5 [i,j,k,l,φ0] + per-entry k).
 *
 * Output is compatible with ForceField: { springs: Float64Array [i,j,r0…],
 * springK: Float64Array per-spring γ_ij, pseudoDihedrals: Float64Array
 * [i,j,k,l,φ0…], pseudoDihedralK: Float64Array, ss: string[] per-quadruple }.
 *
 * Units: Å, kcal/mol(/rad² for dihedrals). No Node deps, `?v=` compatible.
 *
 * References:
 *   [1] Tirion, PRL 77, 1905 (1996).  [2] Atilgan et al., Biophys J 80, 505
 *   (2001) — ANM.  [3] Bahar et al., Folding & Design 2, 173 (1997).
 */

function dist3(pos, i, j) {
  const dx = pos[3 * j] - pos[3 * i];
  const dy = pos[3 * j + 1] - pos[3 * i + 1];
  const dz = pos[3 * j + 2] - pos[3 * i + 2];
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

function angle3(pos, i, j, k) {
  const ax = pos[3 * i] - pos[3 * j], ay = pos[3 * i + 1] - pos[3 * j + 1], az = pos[3 * i + 2] - pos[3 * j + 2];
  const bx = pos[3 * k] - pos[3 * j], by = pos[3 * k + 1] - pos[3 * j + 1], bz = pos[3 * k + 2] - pos[3 * j + 2];
  const la = Math.hypot(ax, ay, az) || 1e-12;
  const lb = Math.hypot(bx, by, bz) || 1e-12;
  let c = (ax * bx + ay * by + az * bz) / (la * lb);
  c = Math.min(1, Math.max(-1, c));
  return Math.acos(c);
}

function dihedral4(pos, i, j, k, l) {
  const ijx = pos[3 * j] - pos[3 * i], ijy = pos[3 * j + 1] - pos[3 * i + 1], ijz = pos[3 * j + 2] - pos[3 * i + 2];
  const jkx = pos[3 * k] - pos[3 * j], jky = pos[3 * k + 1] - pos[3 * j + 1], jkz = pos[3 * k + 2] - pos[3 * j + 2];
  const klx = pos[3 * l] - pos[3 * k], kly = pos[3 * l + 1] - pos[3 * k + 1], klz = pos[3 * l + 2] - pos[3 * k + 2];
  const mx = ijy * jkz - ijz * jky, my = ijz * jkx - ijx * jkz, mz = ijx * jky - ijy * jkx;
  const nx = jky * klz - jkz * kly, ny = jkz * klx - jkx * klz, nz = jkx * kly - jky * klx;
  const m2 = mx * mx + my * my + mz * mz, n2 = nx * nx + ny * ny + nz * nz;
  const jk2 = jkx * jkx + jky * jky + jkz * jkz;
  if (m2 < 1e-12 || n2 < 1e-12 || jk2 < 1e-12) return 0;
  const im = 1 / Math.sqrt(m2), inn = 1 / Math.sqrt(n2);
  const jk = Math.sqrt(jk2);
  const cosPhi = Math.max(-1, Math.min(1, (mx * nx + my * ny + mz * nz) * im * inn));
  const mxnX = my * nz - mz * ny, mxnY = mz * nx - mx * nz, mxnZ = mx * ny - my * nx;
  const sinPhi = ((mxnX * jkx + mxnY * jky + mxnZ * jkz) * im * inn) / jk;
  return Math.atan2(sinPhi, cosPhi);
}

/**
 * Classify local secondary structure from Cα geometry.
 * @param {number} r03  distance Cα_i..Cα_{i+3} (Å)
 * @param {number} phi  pseudo-dihedral (rad)
 * @param {number} thetaMid  middle pseudo-angle (rad)
 * @returns {"helix"|"sheet"|"loop"}
 */
export function classifySS(r03, phi, thetaMid) {
  const adeg = Math.abs(phi) * 180 / Math.PI;
  const tdeg = thetaMid * 180 / Math.PI;
  if (r03 < 6.5 && adeg < 80) return "helix";
  if (r03 > 8.0 && tdeg > 110) return "sheet";
  return "loop";
}

/** Per-SS pseudo-dihedral stiffness (kcal/mol/rad²). */
export const SS_DIHEDRAL_K = { helix: 2.0, sheet: 1.5, loop: 0.5 };

/**
 * Build a Tirion distance-weighted ENM + SS-dependent dihedral basins.
 *
 * @param {ArrayLike<number>} cAlphaPositions  flat 3n Cα array (native)
 * @param {object} [opts]
 * @param {number} [opts.gamma0=1.0]  overall stiffness (kcal/mol/Å²)
 * @param {number} [opts.R0=3.81]  reference distance (Å)
 * @param {number} [opts.cutoff=10]  ENM cutoff Rc (Å)
 * @param {number} [opts.seqSep=3]  min sequence separation (|i−j| ≥ seqSep)
 * @param {Array<Array<number>>|null} [opts.segments=null]  contiguous [s,e) index ranges; cross-segment pairs use |i−j| rule only when segments omitted
 * @param {boolean} [opts.includeDihedrals=true]  build pseudo-dihedral basins
 * @returns {{springs:Float64Array, springK:Float64Array, pseudoDihedrals:Float64Array, pseudoDihedralK:Float64Array, ss:string[]}}
 */
export function buildTirionNetwork(cAlphaPositions, opts = {}) {
  const gamma0 = opts.gamma0 ?? 1.0;
  const R0 = opts.R0 ?? 3.81;
  const cutoff = opts.cutoff ?? 10;
  const seqSep = opts.seqSep ?? 3;
  const segments = opts.segments ?? null;
  const includeDihedrals = opts.includeDihedrals ?? true;
  const n = Math.floor(cAlphaPositions.length / 3);

  const inSameSegment = (i, j) => {
    if (!segments) return true;
    for (const [s, e] of segments) {
      if (i >= s && i < e && j >= s && j < e) return true;
    }
    return false;
  };
  const seqExcluded = (i, j) => {
    if (segments) {
      // Same rule as ForceField: 1-2/1-3 within a segment excluded; cross-segment kept
      for (const [s, e] of segments) {
        if (i >= s && j < e && j - i <= 2 && j > i) return true;
        if (j >= s && i < e && i - j <= 2 && i > j) return true;
      }
      return false;
    }
    return Math.abs(i - j) < seqSep;
  };

  const S = [];
  const K = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (seqExcluded(i, j)) continue;
      if (segments && !inSameSegment(i, j)) {
        // cross-segment (multi-chain complex): keep — interfacial contacts matter
      }
      const r0 = dist3(cAlphaPositions, i, j);
      if (r0 > cutoff || r0 < 1e-6) continue;
      const gij = gamma0 * Math.pow(R0 / r0, 6);
      S.push(i, j, r0);
      K.push(gij);
    }
  }

  const D = [];
  const DK = [];
  const ss = [];
  if (includeDihedrals) {
    const quads = [];
    if (segments) {
      for (const [s, e] of segments) {
        for (let i = s; i + 3 < e; i++) quads.push(i);
      }
    } else {
      for (let i = 0; i + 3 < n; i++) quads.push(i);
    }
    for (const i of quads) {
      const phi0 = dihedral4(cAlphaPositions, i, i + 1, i + 2, i + 3);
      const thMid = angle3(cAlphaPositions, i + 1, i + 2, i + 3);
      const r03 = dist3(cAlphaPositions, i, i + 3);
      const cls = classifySS(r03, phi0, thMid);
      D.push(i, i + 1, i + 2, i + 3, phi0);
      DK.push(SS_DIHEDRAL_K[cls]);
      ss.push(cls);
    }
  }
  return {
    springs: new Float64Array(S),
    springK: new Float64Array(K),
    pseudoDihedrals: new Float64Array(D),
    pseudoDihedralK: new Float64Array(DK),
    ss,
  };
}

/**
 * Tirion pair stiffness for one native distance (helper for tests/adapters).
 * @param {number} r0  native distance (Å)
 * @param {number} [gamma0=1.0]
 * @param {number} [R0=3.81]
 * @returns {number} γ_ij
 */
export function tirionGamma(r0, gamma0 = 1.0, R0 = 3.81) {
  if (!(r0 > 1e-9)) return gamma0;
  return gamma0 * Math.pow(R0 / r0, 6);
}

/**
 * Apply a Tirion network to a ForceField-compatible object in place.
 * Sets ff.springs / ff.springK and stores CG dihedral basins on
 * ff.tirionDihedrals / ff.tirionDihedralK / ff.tirionSS (non-breaking extras;
 * existing compute() ignores them unless the caller wires a dihedral kernel).
 * Also refreshes ff._excluded for the new spring set when present.
 * @param {object} ff  ForceField instance (or compatible {springs, springK})
 * @param {{springs:Float64Array, springK:Float64Array, pseudoDihedrals:Float64Array, pseudoDihedralK:Float64Array, ss:string[]}} net
 */
export function applyTirionToForceField(ff, net) {
  ff.springs = net.springs;
  ff.springK = net.springK;
  ff.tirionDihedrals = net.pseudoDihedrals;
  ff.tirionDihedralK = net.pseudoDihedralK;
  ff.tirionSS = net.ss;
  ff.springScaleActive = true;
  if (ff._excluded && typeof ff._pairKey === "function") {
    for (let k = 0; k < net.springs.length; k += 3) {
      ff._excluded.add(ff._pairKey(net.springs[k], net.springs[k + 1]));
    }
  }
}
