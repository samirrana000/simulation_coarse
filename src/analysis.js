/**
 * analysis.js — post-run analysis of a recorded trajectory.
 *
 * Turns a Recorder's captured frames into research-grade numbers:
 *
 *  1. Fluctuations → simulated B-factors
 *     Per-bead mean-square fluctuation about the native structure,
 *     B = (8π²/3)·⟨Δr²⟩, Pearson-correlated against the experimental
 *     temperature factors carried by the loaded PDB (beads[i].bfac).
 *
 *  2. Essential dynamics + RMSIP
 *     Each frame is least-squares superposed onto the native structure
 *     (Kabsch), the trajectory covariance of the protein Cα beads is
 *     eigen-decomposed (subspace iteration → top-k principal components),
 *     and the root-mean-square inner product (RMSIP) against the k softest
 *     ENM normal modes (shifted-inverse subspace iteration on the elastic
 *     Hessian) measures how much of the observed dynamics the elastic
 *     network explains:
 *         RMSIP = sqrt( (1/k) Σ_{i,j≤k} (v_i·w_j)² )
 *
 *  3. Ligand occupancy (benzene)
 *     A 1 Å atom-density grid over the native ligand volume: the peak
 *     occupancy cell and its distance to the crystal ligand COM (pose
 *     recovery), plus the fraction of frames the ligand COM spends inside
 *     the pocket ("bound fraction").
 *
 *  4. Contact lifetimes
 *     Protein–ligand pairs within 6 Å are tracked frame-to-frame; the mean
 *     and longest uninterrupted residence streaks (ps) quantify how
 *     persistent each native/non-native contact is.
 *
 *  5. Binding free energy
 *     ΔG_bind from the well-tempered metadynamics bias (funnel.estimateDG)
 *     and a CSV of the reconstructed PMF for export.
 *
 * The module is DOM-free (main.js handles downloads) so it runs under plain
 * Node for unit testing. Units: Å, ps, kcal/mol.
 */

/** Pearson correlation coefficient between two arrays of equal length. */
export function pearson(a, b) {
  const n = a.length;
  if (n < 2) return { r: NaN, n };
  let ma = 0, mb = 0;
  for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  let cov = 0, va = 0, vb = 0;
  for (let i = 0; i < n; i++) {
    const da = a[i] - ma, db = b[i] - mb;
    cov += da * db; va += da * da; vb += db * db;
  }
  if (va <= 0 || vb <= 0) return { r: NaN, n };
  return { r: cov / Math.sqrt(va * vb), n };
}

/** dot(a, b) for equal-length Float64Array/Array. */
function dot(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}

/* ------------------------------------------------------------------ */
/*  Kabsch superposition (SVD-based)                                   */
/* ------------------------------------------------------------------ */

/** Jacobi eigensolver for a symmetric 3×3 matrix (row-major, 9 entries). */
function eigenSymmetric3(M) {
  const A = M.slice(), V = [1, 0, 0, 0, 1, 0, 0, 0, 1];
  for (let sweep = 0; sweep < 60; sweep++) {
    let p = 0, q = 1, mx = 0;
    for (let r = 0; r < 3; r++)
      for (let c = r + 1; c < 3; c++) {
        const v = Math.abs(A[r * 3 + c]);
        if (v > mx) { mx = v; p = r; q = c; }
      }
    if (mx < 1e-14) break;
    const App = A[p * 3 + p], Aqq = A[q * 3 + q], Apq = A[p * 3 + q];
    const theta = (Aqq - App) / (2 * Apq);
    const t = Math.sign(theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
    const c = 1 / Math.sqrt(t * t + 1), s = t * c;
    // classic full update (avoids double-applying to the rotated off-diagonal)
    for (let k = 0; k < 3; k++) {
      if (k === p || k === q) continue;
      const akp = A[k * 3 + p], akq = A[k * 3 + q];
      A[k * 3 + p] = c * akp - s * akq;
      A[k * 3 + q] = s * akp + c * akq;
      A[p * 3 + k] = A[k * 3 + p];
      A[q * 3 + k] = A[k * 3 + q];
    }
    const dApp = c * c * App - 2 * s * c * Apq + s * s * Aqq;
    const dAqq = s * s * App + 2 * s * c * Apq + c * c * Aqq;
    A[p * 3 + p] = dApp;
    A[q * 3 + q] = dAqq;
    A[p * 3 + q] = 0;
    A[q * 3 + p] = 0;
    // accumulate rotation into eigenvectors: V ← V·J
    for (let k = 0; k < 3; k++) {
      const vkp = V[k * 3 + p], vkq = V[k * 3 + q];
      V[k * 3 + p] = c * vkp - s * vkq;
      V[k * 3 + q] = s * vkp + c * vkq;
    }
  }
  const vals = [A[0], A[4], A[8]];
  const order = [0, 1, 2].sort((a, b) => vals[b] - vals[a]);
  return { vals: order.map((i) => vals[i]), vecs: order.map((i) => [V[0 * 3 + i], V[1 * 3 + i], V[2 * 3 + i]]) };
}

/**
 * Best-fit rotation + translation mapping point set P onto Q (Kabsch).
 * @param {ArrayLike} P  coordinates, length 3n (may be Float32Array)
 * @param {ArrayLike} Q  coordinates, length 3n
 * @param {number} n     point count
 * @returns {{R: number[9], t: number[3]}}  such that R·p + t ≈ q
 */
export function kabsch(P, Q, n) {
  let px = 0, py = 0, pz = 0, qx = 0, qy = 0, qz = 0;
  for (let i = 0; i < n; i++) {
    px += P[3 * i]; py += P[3 * i + 1]; pz += P[3 * i + 2];
    qx += Q[3 * i]; qy += Q[3 * i + 1]; qz += Q[3 * i + 2];
  }
  px /= n; py /= n; pz /= n; qx /= n; qy /= n; qz /= n;
  const H = [0, 0, 0, 0, 0, 0, 0, 0, 0];
  for (let i = 0; i < n; i++) {
    const ax = P[3 * i] - px, ay = P[3 * i + 1] - py, az = P[3 * i + 2] - pz;
    const bx = Q[3 * i] - qx, by = Q[3 * i + 1] - qy, bz = Q[3 * i + 2] - qz;
    H[0] += ax * bx; H[1] += ax * by; H[2] += ax * bz;
    H[3] += ay * bx; H[4] += ay * by; H[5] += ay * bz;
    H[6] += az * bx; H[7] += az * by; H[8] += az * bz;
  }
  // H = U Σ V^T  via eigen of H^T H (V, Σ²) then U = H V Σ⁻¹
  const HtH = [];
  for (let r = 0; r < 3; r++)
    for (let c = 0; c < 3; c++) {
      let s = 0;
      for (let k = 0; k < 3; k++) s += H[k * 3 + r] * H[k * 3 + c];
      HtH.push(s);
    }
  const { vals, vecs } = eigenSymmetric3(HtH);
  const U = [0, 0, 0, 0, 0, 0, 0, 0, 0];
  for (let c = 0; c < 3; c++) {
    const s = Math.sqrt(Math.max(1e-30, vals[c]));
    const v = vecs[c];
    for (let r = 0; r < 3; r++) {
      const Hr = [H[r * 3 + 0], H[r * 3 + 1], H[r * 3 + 2]];
      U[r * 3 + c] = (Hr[0] * v[0] + Hr[1] * v[1] + Hr[2] * v[2]) / s;
    }
  }
  // proper rotation: reflect columns so det(V U^T) = +1
  let det = 0;
  {
    const M = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    // vecs[c] is COLUMN c of V ⇒ V_{r,k} = vecs[k][r], V·U^T uses V_{r,k}·U_{c,k}
    for (let r = 0; r < 3; r++)
      for (let c = 0; c < 3; c++)
        for (let k = 0; k < 3; k++) M[r * 3 + c] += vecs[k][r] * U[c * 3 + k];
    det = M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) + M[2] * (M[3] * M[7] - M[4] * M[6]);
  }
  const d = det < 0 ? -1 : 1;
  // R = V · diag(1,1,d) · U^T
  const R = [0, 0, 0, 0, 0, 0, 0, 0, 0];
  for (let r = 0; r < 3; r++)
    for (let c = 0; c < 3; c++) {
      let s = 0;
      for (let k = 0; k < 3; k++) s += vecs[k][r] * (k === 2 ? d : 1) * U[c * 3 + k];
      R[r * 3 + c] = s;
    }
  const t = [
    qx - (R[0] * px + R[1] * py + R[2] * pz),
    qy - (R[3] * px + R[4] * py + R[5] * pz),
    qz - (R[6] * px + R[7] * py + R[8] * pz),
  ];
  return { R, t };
}

/** Apply rotation + translation to n points. */
function applyRt(R, t, P, out, n) {
  for (let i = 0; i < n; i++) {
    const x = P[3 * i], y = P[3 * i + 1], z = P[3 * i + 2];
    out[3 * i] = R[0] * x + R[1] * y + R[2] * z + t[0];
    out[3 * i + 1] = R[3] * x + R[4] * y + R[5] * z + t[1];
    out[3 * i + 2] = R[6] * x + R[7] * y + R[8] * z + t[2];
  }
}

/* ------------------------------------------------------------------ */
/*  Subspace iteration (top-k eigenvectors of a symmetric operator)    */
/* ------------------------------------------------------------------ */

/**
 * Block power iteration returning the k dominant eigenvectors of a
 * symmetric linear operator given as a matvec closure.
 * @param {(v: Float64Array) => Float64Array} matvec
 * @param {number} dim
 * @param {number} k
 * @param {number} iters
 * @returns {Array<{vec: Float64Array, eig: number}>} sorted by eig desc
 */
export function subspaceIteration(matvec, dim, k, iters = 40) {
  let V = [];
  for (let c = 0; c < k; c++) {
    const v = new Float64Array(dim);
    for (let i = 0; i < dim; i++) v[i] = Math.random() * 2 - 1;
    V.push(v);
  }
  let W = V.map((v) => v.slice());
  for (let it = 0; it < iters; it++) {
    for (let c = 0; c < k; c++) W[c] = matvec(V[c]);
    // full Gram–Schmidt orthonormalization (in place)
    for (let c = 0; c < k; c++) {
      for (let r = 0; r < c; r++) {
        const d = dot(W[r], W[c]);
        for (let i = 0; i < dim; i++) W[c][i] -= d * W[r][i];
      }
      let nr = 0;
      for (let i = 0; i < dim; i++) nr += W[c][i] * W[c][i];
      nr = Math.sqrt(nr);
      if (nr > 1e-12) for (let i = 0; i < dim; i++) W[c][i] /= nr;
    }
    V = W.map((v) => v.slice());
  }
  const out = [];
  for (let c = 0; c < k; c++) {
    const w = matvec(V[c]);
    out.push({ vec: V[c], eig: dot(V[c], w) });
  }
  out.sort((a, b) => b.eig - a.eig);
  return out;
}

/* ------------------------------------------------------------------ */
/*  ENM Hessian matvec (finite difference of the spring forces)        */
/* ------------------------------------------------------------------ */

/** Pure-ENM harmonic force at configuration x (protein beads only). */
function enmForce(ff, x, nProt) {
  const f = new Float64Array(3 * nProt);
  const S = ff.springs, K = ff.springK;
  for (let s = 0, q = 0; s < S.length; s += 3, q++) {
    const i = 3 * S[s], j = 3 * S[s + 1], r0 = S[s + 2], kk = K[q];
    const dx = x[j] - x[i], dy = x[j + 1] - x[i + 1], dz = x[j + 2] - x[i + 2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
    const sc = (kk * (r - r0)) / r;
    const fx = sc * dx, fy = sc * dy, fz = sc * dz;
    f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
    f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
  }
  return f;
}

/**
 * Hessian·v of the ENM at the native reference via central difference:
 *   H·v = (F(x − εv) − F(x + εv)) / 2ε.
 */
function hessianTimesV(ff, ref, nProt, v) {
  const dim = 3 * nProt;
  const eps = 1e-4;
  const xp = ref.slice(), xm = ref.slice();
  for (let i = 0; i < dim; i++) {
    xp[i] += eps * v[i];
    xm[i] -= eps * v[i];
  }
  const fp = enmForce(ff, xp, nProt);
  const fm = enmForce(ff, xm, nProt);
  const out = new Float64Array(dim);
  for (let i = 0; i < dim; i++) out[i] = (fm[i] - fp[i]) / (2 * eps);
  return out;
}

/** Conjugate-gradient solve of (H + δI) z = b with a matvec closure. */
function cgSolve(Ax, b, dim, maxIt = 80, tol = 1e-8) {
  const x = new Float64Array(dim);
  const r = new Float64Array(b);
  const p = new Float64Array(b);
  let rs = dot(r, r);
  const rs0 = rs;
  if (rs0 < 1e-30) return x;
  for (let it = 0; it < maxIt; it++) {
    const Ap = Ax(p);
    const alpha = rs / Math.max(1e-30, dot(p, Ap));
    for (let i = 0; i < dim; i++) x[i] += alpha * p[i];
    for (let i = 0; i < dim; i++) r[i] -= alpha * Ap[i];
    const rsNew = dot(r, r);
    if (rsNew < tol * rs0) break;
    const beta = rsNew / rs;
    for (let i = 0; i < dim; i++) p[i] = r[i] + beta * p[i];
    rs = rsNew;
  }
  return x;
}

/** Rigid-body (translation + rotation) basis of the protein, orthonormalized. */
function rigidBasis(ref, nProt, dim) {
  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < nProt; i++) {
    cx += ref[3 * i]; cy += ref[3 * i + 1]; cz += ref[3 * i + 2];
  }
  cx /= nProt; cy /= nProt; cz /= nProt;
  const modes = [];
  const push = (arr) => {
    for (const m of modes) {
      const d = dot(m, arr);
      for (let i = 0; i < dim; i++) arr[i] -= d * m[i];
    }
    let nr = 0;
    for (let i = 0; i < dim; i++) nr += arr[i] * arr[i];
    nr = Math.sqrt(nr);
    if (nr > 1e-12) {
      for (let i = 0; i < dim; i++) arr[i] /= nr;
      modes.push(arr);
    }
  };
  // translations
  const tx = new Float64Array(dim), ty = new Float64Array(dim), tz = new Float64Array(dim);
  for (let i = 0; i < nProt; i++) { tx[3 * i] = 1; ty[3 * i + 1] = 1; tz[3 * i + 2] = 1; }
  push(tx); push(ty); push(tz);
  // rotations about the centroid
  for (const [ux, uy, uz] of [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) {
    const w = new Float64Array(dim);
    for (let i = 0; i < nProt; i++) {
      const rx = ref[3 * i] - cx, ry = ref[3 * i + 1] - cy, rz = ref[3 * i + 2] - cz;
      w[3 * i] = uy * rz - uz * ry;
      w[3 * i + 1] = uz * rx - ux * rz;
      w[3 * i + 2] = ux * ry - uy * rx;
    }
    push(w);
  }
  return modes;
}

/**
 * The k softest internal ENM normal modes via shifted-inverse subspace
 * iteration: largest eigenvectors of (H + δI)⁻¹, with the 6 rigid-body
 * modes deflated at every step.
 */
function softEnmModes(ff, ref, nProt, k, dim) {
  const delta = 0.5;
  const rigid = rigidBasis(ref, nProt, dim);
  const inv = (v) => cgSolve((w) => {
    const hw = hessianTimesV(ff, ref, nProt, w);
    for (let i = 0; i < dim; i++) hw[i] += delta * w[i];
    return hw;
  }, v, dim, 80, 1e-7);
  const deflate = (arr) => {
    for (const m of rigid) {
      const d = dot(m, arr);
      for (let i = 0; i < dim; i++) arr[i] -= d * m[i];
    }
    return arr;
  };
  const modes = subspaceIteration((v) => deflate(inv(v)), dim, k, 12);
  // eigenvalues of (H+δI)⁻¹ → λ_H = 1/μ − δ
  return modes.map(({ vec, eig }) => ({ vec, eig: 1 / Math.max(1e-9, eig) - delta }));
}

/** RMSIP between two sets of k unit modes: sqrt((1/k)ΣΣ (v_i·w_j)²). */
function rmsip(a, b, k) {
  let s = 0;
  for (let i = 0; i < k; i++)
    for (let j = 0; j < k; j++) {
      const d = dot(a[i].vec, b[j].vec);
      s += d * d;
    }
  return Math.sqrt(s / k);
}

/* ------------------------------------------------------------------ */
/*  Main analysis driver                                               */
/* ------------------------------------------------------------------ */

const RMSF_TO_B = (8 * Math.PI * Math.PI) / 3; // B = 8π²/3 ⟨Δr²⟩ ≈ 26.32

/**
 * Analyze a recorded trajectory.
 * @param {object} o
 * @param {Float32Array[]} o.frames   captured snapshots (length 3n each)
 * @param {number[]} o.times          simulation time (ps) of each frame
 * @param {Float64Array} o.ref        native coordinates, length 3n
 * @param {number} o.nProt            protein bead count (indices 0..nProt−1)
 * @param {number} o.n                total particle count
 * @param {Array} o.beads             selected protein beads (length nProt) with
 *   {resSeq, resName, chain, bfac}
 * @param {object} o.ff               ForceField (springs / springK / holo etc.)
 * @param {object|null} o.funnel      Funnel (pocket COM, PMF, ΔG) or null
 * @param {number} [o.k=10]           number of modes for PCA/RMSIP
 * @returns {object} report (numbers + human-readable `lines`)
 */
export function analyzeTrajectory({ frames, times, ref, nProt, n, beads, ff, funnel, k = 10 }) {
  const nF = frames.length;
  if (nF < 2) throw new Error("Analysis needs ≥ 2 recorded frames — press ● Rec, run, then Stop.");
  const dim = 3 * nProt;
  const lines = [];
  const rep = { nFrames: nF, spanPs: (times[nF - 1] - times[0]) || 0 };
  const dtPs = nF > 1 ? rep.spanPs / (nF - 1) : 0;
  lines.push(`Trajectory: ${nF} frames · ${rep.spanPs.toFixed(1)} ps (${(rep.spanPs / 1000).toFixed(3)} ns) · dt ≈ ${dtPs.toFixed(2)} ps`);

  /* ---- 1. Fluctuations → simulated B-factors vs experimental --------- */
  const msf = new Float64Array(nProt);
  for (const fr of frames) {
    for (let i = 0; i < nProt; i++) {
      const dx = fr[3 * i] - ref[3 * i], dy = fr[3 * i + 1] - ref[3 * i + 1], dz = fr[3 * i + 2] - ref[3 * i + 2];
      msf[i] += dx * dx + dy * dy + dz * dz;
    }
  }
  const Bsim = new Float64Array(nProt);
  let msfSum = 0;
  for (let i = 0; i < nProt; i++) { Bsim[i] = RMSF_TO_B * (msf[i] / nF); msfSum += msf[i] / nF; }
  const rmsfMean = Math.sqrt(msfSum / nProt);
  const Bexp = beads.map((b) => (Number.isFinite(b.bfac) ? b.bfac : 0));
  const bs = [], be = [];
  for (let i = 0; i < nProt; i++) if (Bexp[i] > 0 && Number.isFinite(Bsim[i])) { bs.push(Bsim[i]); be.push(Bexp[i]); }
  let bfRes = null;
  if (bs.length >= 2) {
    const pc = pearson(bs, be);
    let mb = 0; for (const v of be) mb += v; mb /= be.length;
    bfRes = { r: pc.r, n: bs.length, meanBsim: meanOf(Bsim), meanBexp: mb };
    lines.push(`[1] B-factors: <RMSF> = ${rmsfMean.toFixed(2)} Å · <B_sim> = ${bfRes.meanBsim.toFixed(1)} Å² · <B_exp> = ${mb.toFixed(1)} Å²`);
    lines.push(`    Pearson r(B_sim, B_exp) = ${pc.r.toFixed(3)} over ${bs.length} residues`);
  } else {
    lines.push("[1] B-factors: no experimental B-values in the structure (skipped)");
  }

  /* ---- 2. Essential dynamics + RMSIP --------------------------------- */
  let rmsipRes = null;
  if (dim <= 720 && nF >= 3) {
    // superpose every frame onto the native structure (Kabsch, protein Cα)
    const sup = [];
    for (const fr of frames) {
      const { R, t } = kabsch(fr, ref, nProt);
      const out = new Float64Array(3 * n);
      for (let i = 0; i < nProt; i++) applyRt(R, t, fr, out, nProt);
      sup.push(out);
    }
    // mean structure of superposed frames
    const mean = new Float64Array(dim);
    for (const fr of sup)
      for (let i = 0; i < dim; i++) mean[i] += fr[i] / nF;
    // covariance matvec C·v = (1/nF) Σ (x_f − x̄)((x_f − x̄)·v)
    const covMatvec = (v) => {
      const out = new Float64Array(dim);
      for (const fr of sup) {
        let d = 0;
        for (let i = 0; i < dim; i++) d += (fr[i] - mean[i]) * v[i];
        for (let i = 0; i < dim; i++) out[i] += (fr[i] - mean[i]) * d / nF;
      }
      return out;
    };
    const pcModes = subspaceIteration(covMatvec, dim, k, 40);
    // total variance = trace(C) = (1/nF) Σ_f |x_f − x̄|² (cheap, exact)
    let trace = 0;
    for (const fr of sup) {
      let s = 0;
      for (let i = 0; i < dim; i++) { const d = fr[i] - mean[i]; s += d * d; }
      trace += s / nF;
    }
    const totVar = pcModes.reduce((a, m) => a + m.eig, 0) / Math.max(1e-30, trace);
    // softest ENM modes (deflated rigid-body subspace) — the elastic network
    // exists only in Cα mode, so skip the comparison when there are no springs
    let rip = null;
    if (ff.springs.length > 0) {
      const enmModes = softEnmModes(ff, ref, nProt, k, dim);
      rip = rmsip(pcModes, enmModes, k);
    }
    const var1 = pcModes.length ? pcModes[0].eig / Math.max(1e-30, trace) : 0;
    rmsipRes = {
      k,
      rmsip: rip,
      varTopK: Math.min(1, totVar),
      varMode1: Math.min(1, var1),
      pcaEigs: pcModes.slice(0, 5).map((m) => m.eig),
    };
    lines.push(`[2] Essential dynamics (k=${k}, protein Cα): mode 1 variance = ${(var1 * 100).toFixed(1)}% · top-${k} = ${(Math.min(1, totVar) * 100).toFixed(1)}%`);
    if (rip !== null) {
      lines.push(`    RMSIP(PCA, ENM soft modes) = ${rip.toFixed(3)}  (1 = fully captured, 0 = orthogonal)`);
    } else {
      lines.push("    RMSIP vs ENM skipped (elastic network exists only in Cα mode)");
    }
  } else {
    lines.push(dim > 720 ? `[2] RMSIP skipped: ${nProt} residues (${dim} dof) exceeds the 720-dim cap` : "[2] RMSIP skipped: need ≥ 3 frames");
  }

  /* ---- 3. Ligand occupancy ------------------------------------------- */
  const nLig = n - nProt;
  let occ = null;
  if (nLig > 0 && nF >= 1) {
    const lig0 = nProt;
    let lx = 0, ly = 0, lz = 0;
    for (let a = 0; a < nLig; a++) { lx += ref[3 * (lig0 + a)]; ly += ref[3 * (lig0 + a) + 1]; lz += ref[3 * (lig0 + a) + 2]; }
    const ligCOM0 = [lx / nLig, ly / nLig, lz / nLig];
    // pocket COM + radius (funnel, or a native-ligand fallback)
    let pcom = ligCOM0, rPocket = 8;
    if (funnel && funnel.active) { pcom = funnel.pocketCOM; rPocket = funnel.rPocket; }
    // 3D occupancy grid (1 Å) over the native ligand volume, padded 4 Å
    let minX = Infinity, minY = Infinity, minZ = Infinity, maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
    for (let a = 0; a < nLig; a++) {
      const c = 3 * (lig0 + a);
      minX = Math.min(minX, ref[c]); minY = Math.min(minY, ref[c + 1]); minZ = Math.min(minZ, ref[c + 2]);
      maxX = Math.max(maxX, ref[c]); maxY = Math.max(maxY, ref[c + 1]); maxZ = Math.max(maxZ, ref[c + 2]);
    }
    minX -= 4; minY -= 4; minZ -= 4;
    const nx = Math.max(1, Math.ceil(maxX - minX + 8) + 1), ny = Math.max(1, Math.ceil(maxY - minY + 8) + 1), nz = Math.max(1, Math.ceil(maxZ - minZ + 8) + 1);
    const grid = new Float64Array(nx * ny * nz);
    let boundCount = 0;
    for (const fr of frames) {
      // ligand COM distance → bound fraction
      let cx = 0, cy = 0, cz = 0;
      for (let a = 0; a < nLig; a++) { const c = 3 * (lig0 + a); cx += fr[c]; cy += fr[c + 1]; cz += fr[c + 2]; }
      cx /= nLig; cy /= nLig; cz /= nLig;
      const drx = cx - pcom[0], dry = cy - pcom[1], drz = cz - pcom[2];
      if (Math.sqrt(drx * drx + dry * dry + drz * drz) <= rPocket) boundCount++;
      // per-atom density
      for (let a = 0; a < nLig; a++) {
        const c = 3 * (lig0 + a);
        const gx = Math.min(nx - 1, Math.max(0, Math.floor((fr[c] - minX))));
        const gy = Math.min(ny - 1, Math.max(0, Math.floor((fr[c + 1] - minY))));
        const gz = Math.min(nz - 1, Math.max(0, Math.floor((fr[c + 2] - minZ))));
        grid[gx * ny * nz + gy * nz + gz]++;
      }
    }
    let peak = -1, pk = 0;
    for (let i = 0; i < grid.length; i++) if (grid[i] > peak) { peak = grid[i]; pk = i; }
    const pxx = pk / (ny * nz), rem = pk % (ny * nz), pyy = Math.floor(rem / nz), pzz = rem % nz;
    const peakCenter = [minX + pxx + 0.5, minY + pyy + 0.5, minZ + pzz + 0.5];
    const dPeak = Math.sqrt(
      (peakCenter[0] - ligCOM0[0]) ** 2 + (peakCenter[1] - ligCOM0[1]) ** 2 + (peakCenter[2] - ligCOM0[2]) ** 2
    );
    const totalAtomFrames = nF * nLig;
    occ = {
      boundFrac: boundCount / nF,
      rPocket,
      peakCenter,
      peakCount: peak,
      peakFrac: peak / totalAtomFrames,
      peakToNative: dPeak,
    };
    lines.push(`[3] Ligand occupancy (${nLig} atoms): bound fraction = ${(occ.boundFrac * 100).toFixed(1)}% (COM within ${rPocket} Å of pocket COM)`);
    lines.push(`    occupancy peak at (${peakCenter.map((v) => v.toFixed(1)).join(", ")}) · ${occ.peakCount} of ${totalAtomFrames} atom-frames (${(occ.peakFrac * 100).toFixed(2)}%)`);
    lines.push(`    peak → crystal ligand COM: ${dPeak.toFixed(2)} Å`);
  } else {
    lines.push("[3] Ligand occupancy: no ligand in the system (skipped)");
  }

  /* ---- 4. Contact lifetimes ------------------------------------------ */
  let ct = null;
  if (nLig > 0 && nF >= 1) {
    const cutoff = 6.0, cut2 = cutoff * cutoff;
    const streaks = new Map(); // key -> {cur, max, pres}
    for (const fr of frames) {
      const seen = new Set();
      for (let i = 0; i < nProt; i++) {
        const ix = 3 * i;
        for (let a = 0; a < nLig; a++) {
          const c = 3 * (nProt + a);
          const dx = fr[c] - fr[ix], dy = fr[c + 1] - fr[ix + 1], dz = fr[c + 2] - fr[ix + 2];
          if (dx * dx + dy * dy + dz * dz < cut2) seen.add(i * 100000 + a);
        }
      }
      for (const key of streaks.keys()) if (!seen.has(key)) streaks.get(key).cur = 0;
      for (const key of seen) {
        let s = streaks.get(key);
        if (!s) { s = { cur: 0, max: 0, pres: 0 }; streaks.set(key, s); }
        s.cur++; s.pres++;
        if (s.cur > s.max) s.max = s.cur;
      }
    }
    let sumMax = 0, maxMax = 0;
    const top = [...streaks.entries()].sort((a, b) => b[1].pres - a[1].pres).slice(0, 5);
    for (const [, s] of streaks) { sumMax += s.max; maxMax = Math.max(maxMax, s.max); }
    const meanLife = streaks.size ? sumMax / streaks.size * dtPs : 0;
    ct = {
      nPairs: streaks.size,
      meanLifetimePs: meanLife,
      maxLifetimePs: maxMax * dtPs,
      top: top.map(([key, s]) => {
        const i = Math.floor(key / 100000), a = key % 100000;
        return { bead: i, lig: a, res: beads[i] ? `${beads[i].resName}${beads[i].resSeq}:${beads[i].chain}` : String(i), residence: s.pres / nF, maxStreakPs: s.max * dtPs };
      }),
    };
    lines.push(`[4] Contact lifetimes (P–L pairs < ${cutoff} Å): mean streak = ${meanLife.toFixed(1)} ps · max = ${(maxMax * dtPs).toFixed(1)} ps · ${streaks.size} pairs`);
    lines.push(`    most persistent: ${ct.top.map((t) => `${t.res}–lig${t.lig + 1} (res ${(t.residence * 100).toFixed(0)}%, max ${t.maxStreakPs.toFixed(1)} ps)`).join(" · ")}`);
  } else {
    lines.push("[4] Contact lifetimes: no ligand in the system (skipped)");
  }

  /* ---- 5. Binding free energy ---------------------------------------- */
  let dg = null;
  if (funnel && funnel.active) {
    const val = funnel.estimateDG();
    dg = { value: Number.isNaN(val) ? null : val, nHills: funnel._nHills };
    if (dg.value !== null) lines.push(`[5] Binding ΔG (WTM bias, bound vs +6 Å) ≈ ${dg.value.toFixed(2)} kcal/mol (${funnel._nHills} hills)`);
    else lines.push("[5] Binding ΔG: not enough metadynamics hills deposited yet");
  } else {
    lines.push("[5] Binding ΔG: funnel not active (no ligand / funnel off)");
  }

  rep.lines = lines;
  rep.bf = bfRes;
  rep.rmsip = rmsipRes;
  rep.occupancy = occ;
  rep.contacts = ct;
  rep.dg = dg;
  return rep;
}

function meanOf(arr) {
  let s = 0;
  for (const v of arr) s += v;
  return arr.length ? s / arr.length : NaN;
}

/**
 * Serialize the reconstructed PMF as CSV (r, pmf) plus a ΔG footer line.
 * @param {object} funnel  active Funnel
 * @returns {string} CSV text
 */
export function pmfCsv(funnel) {
  if (!funnel || !funnel.active) throw new Error("Funnel is not active — no PMF to export.");
  const { r, pmf } = funnel.getPMF();
  const dg = funnel.estimateDG();
  const out = ["r_Ang,pMF_kcal_per_mol"];
  for (let k = 0; k < r.length; k++) out.push(`${r[k].toFixed(3)},${pmf[k].toFixed(4)}`);
  out.push(`# dG_bind_kcal_per_mol,${Number.isNaN(dg) ? "nan" : dg.toFixed(4)}`);
  out.push(`# hills,${funnel._nHills}`);
  return out.join("\n") + "\n";
}
