/**
 * thermodynamics.js — ΔH/ΔS binding decomposition (Loop-2 S5, R4 §2 pipeline).
 *
 * ΔG = ΔH − TΔS with:
 *   ΔH       = ⟨U_bind⟩_holo − ⟨U_bind⟩_apo  (apo binding ≡ 0 by construction;
 *              component split: LJ / Coul / HB / desolv / π / cation-π / halogen
 *              from the S4 per-term accumulators)
 *   ΔS_pocket = Schlitter(holo pocket) − Schlitter(apo pocket)   [holo − apo: negative]
 *   ΔS_lig    = torsion-Shannon(bound) − torsion-Shannon(bulk); rigid ⇒ 0
 *   ΔS_solv   = −ΔSASA × 0.012 kcal/mol/Å²  (R4 scale, ±50% band)
 *
 * Schlitter: S ≤ (kB/2)·ln det(I + kBT·e²·σ̄/kB²), σ̄ mass-weighted covariance,
 * Cholesky ln-det (verified vs analytic Gaussian: err 0.20%, r4 prototype).
 */

const KB = 0.0019872041; // kcal/mol/K

import { findRotatableBonds, autoTorsions } from "./rotbonds.js";
export { findRotatableBonds, autoTorsions };

/** ln det of symmetric PD matrix via Cholesky. −Infinity if not PD. */
function lnDetCholesky(M) {
  const n = M.length;
  let lnD = 0;
  for (let i = 0; i < n; i++) {
    let s = M[i][i];
    for (let k = 0; k < i; k++) s -= M[i][k] * M[i][k];
    if (s <= 0) return -Infinity;
    lnD += Math.log(s);
    for (let j = i + 1; j < n; j++) {
      let t = M[j][i];
      for (let k = 0; k < i; k++) t -= M[j][k] * M[i][k];
      M[j][i] = t / s;
    }
    M[i][i] = s;
  }
  return 2 * lnD;
}

/**
 * Schlitter quasi-harmonic entropy (upper bound).
 * @param {number[][]} cov covariance matrix (DOF × DOF), Å²
 * @param {number[]|number} masses per-DOF mass (Da) or single value
 * @param {number} T Kelvin
 * @returns {number} kcal/mol/K
 */
export function schlitterEntropy(cov, masses, T = 300) {
  const n = cov.length;
  const m = typeof masses === "number" ? new Array(n).fill(masses) : masses;
  const e = Math.E;
  const A = cov.map((row, i) =>
    row.map((v, j) => (KB * T * e * e / (KB * KB)) * v * Math.sqrt(m[i] * m[j]))
  );
  const M = A.map((row, i) => row.map((v, j) => (i === j ? 1 + v : v)));
  const lnDet = lnDetCholesky(M.map((r) => [...r]));
  return (KB / 2) * lnDet;
}

/** Kabsch-align all frames onto frame 0 (in-place copy, drift removal) and
 *  return per-DOF-aligned frames. */
function kabschAlignFrames(frames) {
  const n3 = frames[0].length;
  // centroid of frame 0
  const c0 = centroid(frames[0]);
  // rotation of each frame onto frame 0 via Kabsch (SVD-free 3×3 via quaternion
  // — simplified Jacobi eigen for the 3×3 cross-covariance)
  const out = frames.map((f) => Float32Array.from(f));
  const ref = new Float32Array(n3);
  for (let k = 0; k < n3; k += 3) {
    ref[k] = frames[0][k] - c0[0]; ref[k + 1] = frames[0][k + 1] - c0[1]; ref[k + 2] = frames[0][k + 2] - c0[2];
  }
  for (let fi = 0; fi < out.length; fi++) {
    const f = out[fi];
    const cf = centroid(f);
    // center
    for (let k = 0; k < n3; k += 3) {
      f[k] -= cf[0]; f[k + 1] -= cf[1]; f[k + 2] -= cf[2];
    }
    // covariance H = Σ ref ⊗ f  (3×3), find rotation R: f' = R f minimizing RMSD
    const H = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let k = 0; k < n3; k += 3) {
      for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) {
        H[a][b] += ref[k + a] * f[k + b];
      }
    }
    // quaternion-trace method (Horn): build 4×4 K matrix, power-iterate largest eigvec
    const R = hornRotation(H);
    for (let k = 0; k < n3; k += 3) {
      const x = f[k], y = f[k + 1], z = f[k + 2];
      f[k] = R[0][0] * x + R[0][1] * y + R[0][2] * z;
      f[k + 1] = R[1][0] * x + R[1][1] * y + R[1][2] * z;
      f[k + 2] = R[2][0] * x + R[2][1] * y + R[2][2] * z;
    }
  }
  return out;
}
function centroid(f) {
  const n = f.length / 3;
  let x = 0, y = 0, z = 0;
  for (let k = 0; k < f.length; k += 3) { x += f[k]; y += f[k + 1]; z += f[k + 2]; }
  return [x / n, y / n, z / n];
}
/** Horn quaternion rotation from cross-covariance H (maps f → ref). */
function hornRotation(H) {
  // K matrix (Horn 1987 closed form via largest eigenvalue power iteration)
  const Sxx = H[0][0], Sxy = H[0][1], Sxz = H[0][2];
  const Syx = H[1][0], Syy = H[1][1], Syz = H[1][2];
  const Szx = H[2][0], Szy = H[2][1], Szz = H[2][2];
  const K = [
    [Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
    [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
    [Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy],
    [Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz],
  ];
  // power iteration for dominant eigenvector
  let v = [1, 0, 0, 0];
  for (let it = 0; it < 64; it++) {
    const nv = [
      K[0][0] * v[0] + K[0][1] * v[1] + K[0][2] * v[2] + K[0][3] * v[3],
      K[1][0] * v[0] + K[1][1] * v[1] + K[1][2] * v[2] + K[1][3] * v[3],
      K[2][0] * v[0] + K[2][1] * v[1] + K[2][2] * v[2] + K[2][3] * v[3],
      K[3][0] * v[0] + K[3][1] * v[1] + K[3][2] * v[2] + K[3][3] * v[3],
    ];
    const l = Math.hypot(...nv) || 1;
    v = nv.map((x) => x / l);
  }
  const [w, x, y, z] = v;
  return [
    [w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
    [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
    [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z],
  ];
}

/** Covariance (unbiased) of aligned frames restricted to DOF indices. */
function covariance(frames, dofIdx) {
  const nF = frames.length, n = dofIdx.length;
  const mean = new Float64Array(n);
  for (const f of frames) for (let i = 0; i < n; i++) mean[i] += f[dofIdx[i]] / nF;
  const cov = Array.from({ length: n }, () => new Float64Array(n));
  for (let a = 0; a < n; a++)
    for (let b = a; b < n; b++) {
      let s = 0;
      for (const f of frames) s += (f[dofIdx[a]] - mean[a]) * (f[dofIdx[b]] - mean[b]);
      cov[a][b] = cov[b][a] = s / (nF - 1);
    }
  return cov.map((r) => Array.from(r));
}

/** Torsion Shannon entropy: S = −kB Σ p ln p over 30° bins. */
export function torsionEntropy(frames, torsions, binDeg = 30) {
  if (!frames.length || !torsions.length) return 0;
  const nbins = Math.round(360 / binDeg);
  const total = new Array(torsions.length).fill(0).map(() => new Float64Array(nbins));
  for (const f of frames) {
    for (let t = 0; t < torsions.length; t++) {
      const [a, b, c, d] = torsions[t];
      const phi = dihedral(f, a, b, c, d);
      let bin = Math.floor(((phi + 180) / 360) * nbins);
      if (bin >= nbins) bin = nbins - 1;
      total[t][bin]++;
    }
  }
  let S = 0;
  const n = frames.length;
  for (const hist of total) {
    for (const c of hist) {
      if (c > 0) { const p = c / n; S += -KB * p * Math.log(p); }
    }
  }
  return S;
}

function dihedral(f, a, b, c, d) {
  // signed dihedral (IUPAC convention, −180..180)
  const r = (i) => [f[3 * i], f[3 * i + 1], f[3 * i + 2]];
  const sub = (u, v) => [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
  const cr = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const dt = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const b1 = r(a), b2 = r(b), b3 = r(c), b4 = r(d);
  const n1 = cr(sub(b2, b1), sub(b3, b2));
  const n2 = cr(sub(b3, b2), sub(b4, b3));
  const m1 = Math.hypot(...n1) || 1, m2 = Math.hypot(...n2) || 1;
  const cphi = dt(n1, n2) / (m1 * m2);
  const ang = Math.acos(Math.max(-1, Math.min(1, cphi)));
  const sign = dt(cr(n1, sub(b3, b2)), n2) >= 0 ? 1 : -1;
  return sign * (ang * 180) / Math.PI;
}

/**
 * Full ΔH/ΔS pipeline (R4 §2).
 * @param {object} p
 * @param {number[][]|Float32Array[]} p.holoFrames
 * @param {number[][]|Float32Array[]} p.apoFrames
 * @param {number[][]} [p.holoEnergies] per-frame 7-vector [lj,coul,hb,desolv,pi,cpi,xb]
 * @param {number[]} p.pocketIdx residue indices (bead space)
 * @param {number} [p.nProt] protein bead count (pocket indices absolute)
 * @param {number} [p.mass] bead mass Da (110; used when p.masses absent)
 * @param {number[]|Float64Array} [p.masses] per-atom masses Da in atom space
 *   (Stage-1 heavy path: pocket DOFs inherit their atom's mass; falls back to
 *   p.mass when absent or too short — CG callers unaffected)
 * @param {number} [p.T] 300
 * @param {number[]} [p.torsions] ligand rotatable-bond quadruples (absolute atom idx)
 * @param {object} [p.ligand] ligand graph for auto-detection fallback
 *   { atoms, bonds, offset?, orders?, aromatic? } — when p.torsions is empty,
 *   rotatable torsions are auto-detected via rotbonds.autoTorsions (Stage-4).
 *   Shorthand: p.ligandAtoms + p.ligandBonds (+ p.ligandStart | p.nProt as
 *   offset) is also accepted. Additive only: explicit p.torsions always wins;
 *   no ligand graph ⇒ legacy rigid-ligand path (ΔS_lig = 0), bit-identical.
 * @param {number} [p.dsasa] precomputed ΔSASA Å² (else contact-count × 10 Å² proxy)
 */
export function computeThermodynamics(p) {
  const T = p.T ?? 300, mass = p.mass ?? 110;
  const TERMS = ["lj", "coul", "hb", "desolv", "pi", "cpi", "xb"];
  // ---- ΔH ----
  const dH = { total: 0 };
  TERMS.forEach((t) => (dH[t] = 0));
  let dH_se = 0;
  if (p.holoEnergies && p.holoEnergies.length) {
    const nF = p.holoEnergies.length;
    const sums = new Array(7).fill(0);
    for (const row of p.holoEnergies) for (let i = 0; i < 7; i++) sums[i] += row[i] / nF;
    for (let i = 0; i < 7; i++) dH[TERMS[i]] = sums[i];
    dH.total = sums.reduce((a, b) => a + b, 0);
    // block bootstrap SE on total (20 blocks)
    const B = 20, bs = Math.max(1, Math.floor(nF / B));
    const blockMeans = [];
    for (let bI = 0; bI < B; bI++) {
      let s = 0, c = 0;
      for (let k = bI * bs; k < Math.min(nF, (bI + 1) * bs); k++) {
        s += p.holoEnergies[k].reduce((a, x) => a + x, 0); c++;
      }
      if (c) blockMeans.push(s / c);
    }
    if (blockMeans.length > 1) {
      const mb = blockMeans.reduce((a, b) => a + b, 0) / blockMeans.length;
      dH_se = Math.sqrt(blockMeans.reduce((s, m) => s + (m - mb) ** 2, 0) / (blockMeans.length - 1) / blockMeans.length) * blockMeans.length ** 0.25;
      // (SE of the mean over blocks, scaled back to per-frame mean)
      dH_se = Math.sqrt(blockMeans.reduce((s, m) => s + (m - mb) ** 2, 0) / (blockMeans.length - 1)) / Math.sqrt(blockMeans.length);
    }
  }
  // ---- ΔS pocket ----
  let dS_pocket = 0, framesPerDof = 0, massModel = `uniform ${mass} Da`;
  let S_holo_pocket = 0, S_apo_pocket = 0;
  if (p.holoFrames?.length >= 10 && p.apoFrames?.length >= 10 && p.pocketIdx?.length) {
    const dofIdx = [];
    for (const r of p.pocketIdx) for (let k = 0; k < 3; k++) dofIdx.push(3 * r + k);
    // Stage-1: per-atom masses when the caller supplies an atom-space table
    // (heavy FF ff.masses); else the legacy uniform bead mass. Additive only.
    let dofMasses = null;
    if (p.masses && p.masses.length > Math.max(...p.pocketIdx)) {
      dofMasses = [];
      for (const r of p.pocketIdx) for (let k = 0; k < 3; k++) dofMasses.push(p.masses[r]);
      massModel = "per-atom";
    } else {
      dofMasses = new Array(dofIdx.length).fill(mass);
    }
    const holoA = kabschAlignFrames(p.holoFrames);
    const apoA = kabschAlignFrames(p.apoFrames);
    const S_holo = schlitterEntropy(covariance(holoA, dofIdx), dofMasses, T);
    const S_apo = schlitterEntropy(covariance(apoA, dofIdx), dofMasses, T);
    S_holo_pocket = S_holo; S_apo_pocket = S_apo;
    dS_pocket = S_holo - S_apo; // holo − apo (sign is model-dependent: see R4 §5)
    framesPerDof = Math.min(holoA.length, apoA.length) / dofIdx.length;
  }
  // ---- ΔS ligand (torsion) ----
  // Stage-4: autoTorsions fallback — when no explicit torsion list is given
  // but a ligand graph is supplied, detect rotatable bonds automatically.
  let dS_lig = 0, ligNote = "rigid ligand (no rotatable bonds) ⇒ 0";
  let rotatableBonds = p.torsions?.length ?? 0;
  let torsions = p.torsions;
  let ligAuto = false;
  if ((!torsions || !torsions.length) && p.holoFrames?.length) {
    const lig = p.ligand
      ?? ((p.ligandAtoms && p.ligandBonds)
        ? { atoms: p.ligandAtoms, bonds: p.ligandBonds, offset: p.ligandStart ?? p.nProt ?? 0 }
        : null);
    if (lig?.atoms?.length && lig?.bonds?.length) {
      try {
        torsions = autoTorsions(lig.atoms, lig.bonds, lig.offset ?? 0, { orders: lig.orders, aromatic: lig.aromatic });
        rotatableBonds = torsions.length;
        ligAuto = true;
      } catch { torsions = []; rotatableBonds = 0; }
    }
  }
  if (torsions?.length && p.holoFrames?.length) {
    dS_lig = torsionEntropy(p.holoFrames, torsions) - torsionEntropy(p.apoFrames ?? [], torsions);
    ligNote = ligAuto
      ? `${torsions.length} auto rotatable bonds (${torsions.length} torsions)`
      : `${torsions.length} torsions`;
  } else if (ligAuto) {
    ligNote = "rigid ligand (0 rotatable bonds) ⇒ 0";
  }
  // ---- ΔS solvent (SASA proxy) ----
  let dS_solv = 0, dsasa = p.dsasa ?? null;
  if (dsasa === null) {
    // contact-count proxy: ΔSASA ≈ nContacts(holo) × 10 Å² (R1 §2 hydrophobic burial scale)
    // derive from holo frames if a contact count series exists; else 0 with note
    dsasa = p.contactCount ? p.contactCount * 10 : 0;
  }
  dS_solv = -dsasa * 0.012 / T; // kcal/mol/K — favorable release (negative ΔG contribution)
  const dS_total = dS_pocket + dS_lig + dS_solv;
  const dG = dH.total - T * dS_total;
  return {
    dH, dH_se,
    dS: { pocket: dS_pocket, ligand: dS_lig, solvent: dS_solv, total: dS_total },
    // Stage-1: absolute pocket entropies (diagnose holo-tightening vs apo-loosening).
    S_pocket: { holo: S_holo_pocket, apo: S_apo_pocket },
    dG_estimate: dG,
    meta: {
      T, framesPerDof, ligNote,
      rotatableBonds,
      dsasa, sasaScale: "0.012 kcal/mol/Å² (±50%)",
      pocketResidues: p.pocketIdx?.length ?? 0,
      holoFrames: p.holoFrames?.length ?? 0, apoFrames: p.apoFrames?.length ?? 0,
      massModel,
      warning: framesPerDof && framesPerDof < 10 ? "UNDER-SAMPLED (frames/DOF < 10) — Schlitter unstable" : null,
    },
  };
}

/** Compact mono table for the analysis <pre>. */
export function formatThermoTable(r) {
  const L = [];
  L.push("── Thermodynamic decomposition (ΔG = ΔH − TΔS) ──");
  L.push(`ΔH  total ${r.dH.total >= 0 ? "+" : ""}${r.dH.total.toFixed(2)} ± ${r.dH_se.toFixed(2)} kcal/mol`);
  for (const [k, v] of Object.entries(r.dH)) {
    if (k === "total") continue;
    L.push(`     ${k.padEnd(7)} ${v >= 0 ? "+" : ""}${v.toFixed(2)}`);
  }
  L.push(`ΔS  pocket  ${r.dS.pocket.toFixed(4)} kcal/mol/K  (−TΔS = ${(-r.meta.T * r.dS.pocket).toFixed(2)})`);
  // Stage-4: ligand line always carries the rotatable-bond count — via ligNote
  // ("N torsions" / "N auto rotatable bonds") plus the explicit suffix below
  // when meta.rotatableBonds is known.
  const rotSuffix = Number.isFinite(r.meta?.rotatableBonds) ? ` [${r.meta.rotatableBonds} rotatable]` : "";
  L.push(`     ligand  ${r.dS.ligand.toFixed(4)}  (${r.meta.ligNote}${rotSuffix})`);
  L.push(`     solvent ${r.dS.solvent.toFixed(4)}  (ΔSASA ${r.meta.dsasa} Å² × 0.012)`);
  L.push(`ΔS total    ${r.dS.total.toFixed(4)} kcal/mol/K  (−TΔS = ${(-r.meta.T * r.dS.total).toFixed(2)})`);
  L.push(`ΔG estimate ${r.dG_estimate.toFixed(2)} kcal/mol  (±bootstrap SE on ΔH; ±50% on SASA term)`);
  if (r.meta.warning) L.push(`⚠ ${r.meta.warning}`);
  if (r.meta.framesPerDof) L.push(`frames/DOF ${r.meta.framesPerDof.toFixed(1)} (need ≥10 for Schlitter)`);
  return L.join("\n");
}
