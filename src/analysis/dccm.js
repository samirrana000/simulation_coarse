/**
 * dccm.js — dynamical cross-correlation matrix (DCCM) from a trajectory.
 *
 * Normalized Cα displacement covariance (Å units cancel):
 *
 *   Cij = ⟨Δri · Δrj⟩ / sqrt(⟨Δri²⟩ ⟨Δrj²⟩)   ∈ [−1, 1],  Cii = 1,
 *
 * where Δri(t) = ri(t) − ⟨ri⟩ after least-squares superposition of every
 * frame onto the reference (Kabsch, imported from src/analysis.js). Cij > 0
 * marks residues moving together (e.g. a domain hinge closing as a rigid
 * unit); Cij < 0 marks anticorrelated motion (e.g. two lobes breathing
 * against each other).
 *
 * Exports: computeDCCM(frames, opts), renderDCCMHeatmap(canvas, matrix, opts),
 * highlightCorrelatedPair(viewer, i, j, opts), clearCorrelatedPair(viewer).
 * DOM-free except the canvas renderer (null-safe under Node) and the viewer
 * hook (duck-typed, never throws).
 */

import { kabsch } from "../analysis.js?v=10";

/**
 * @typedef {object} DCCMResult
 * @property {number} n               residue count (Cα subset size)
 * @property {Float64Array} matrix    row-major n×n correlation matrix
 * @property {Float64Array} variances per-residue ⟨Δri²⟩ (Å²)
 * @property {number} nFrames         frames consumed
 */

/**
 * Compute the normalized cross-correlation matrix.
 *
 * @param {Array<Float32Array|Float64Array>} frames  snapshots, length 3n each
 * @param {object} opts
 * @param {number} opts.nProt   protein bead count (Cα subset = indices 0..nProt−1)
 * @param {Float64Array} [opts.ref] reference coords, length ≥ 3·nProt
 *   (default: Cα part of frames[0])
 * @param {boolean} [opts.align=true]  Kabsch-superpose frames onto ref first
 *   (removes rigid-body translation/rotation that would otherwise correlate
 *   every pair spuriously toward +1)
 * @returns {DCCMResult}
 */
export function computeDCCM(frames, opts = {}) {
  const nProt = opts.nProt;
  if (!Number.isInteger(nProt) || nProt < 2) {
    throw new Error("computeDCCM: opts.nProt (≥ 2) is required.");
  }
  if (!frames || frames.length < 2) {
    throw new Error(`computeDCCM: need ≥ 2 frames (got ${frames ? frames.length : 0}).`);
  }
  const n = nProt;
  const dim = 3 * n;
  for (let f = 0; f < frames.length; f++) {
    if (!frames[f] || frames[f].length < dim) {
      throw new Error(`computeDCCM: frame ${f} too short (${frames[f] ? frames[f].length : 0} < ${dim}).`);
    }
  }
  const refFull = opts.ref ?? frames[0];
  const ref = Float64Array.from(refFull.subarray ? refFull.subarray(0, dim) : refFull.slice(0, dim));

  // Superpose + accumulate mean in one pass (store aligned Cα frames).
  const aligned = new Array(frames.length);
  const mean = new Float64Array(dim);
  const tmp = new Float64Array(dim);
  for (let f = 0; f < frames.length; f++) {
    const fr = frames[f];
    let src = fr;
    if (opts.align !== false) {
      const { R, t } = kabsch(fr, ref, n);
      for (let i = 0; i < n; i++) {
        const x = fr[3 * i], y = fr[3 * i + 1], z = fr[3 * i + 2];
        tmp[3 * i] = R[0] * x + R[1] * y + R[2] * z + t[0];
        tmp[3 * i + 1] = R[3] * x + R[4] * y + R[5] * z + t[1];
        tmp[3 * i + 2] = R[6] * x + R[7] * y + R[8] * z + t[2];
      }
      src = tmp;
    }
    const cp = Float64Array.from(src.subarray ? src.subarray(0, dim) : src.slice(0, dim));
    aligned[f] = cp;
    for (let i = 0; i < dim; i++) mean[i] += cp[i] / frames.length;
  }

  // Covariance numerator Sij = Σ_f Δri·Δrj and variances Vi = Σ_f |Δri|².
  const S = new Float64Array(n * n);
  const V = new Float64Array(n);
  const dx = new Float64Array(dim);
  for (let f = 0; f < frames.length; f++) {
    const cp = aligned[f];
    for (let i = 0; i < dim; i++) dx[i] = cp[i] - mean[i];
    for (let i = 0; i < n; i++) {
      const ax = dx[3 * i], ay = dx[3 * i + 1], az = dx[3 * i + 2];
      V[i] += ax * ax + ay * ay + az * az;
      for (let j = i; j < n; j++) {
        S[i * n + j] += ax * dx[3 * j] + ay * dx[3 * j + 1] + az * dx[3 * j + 2];
      }
    }
  }
  const matrix = new Float64Array(n * n);
  const variances = new Float64Array(n);
  for (let i = 0; i < n; i++) variances[i] = V[i] / frames.length;
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      let c;
      if (i === j) {
        c = 1; // exact by definition, even for a frozen residue
      } else {
        const den = Math.sqrt(Math.max(0, V[i]) * Math.max(0, V[j]));
        c = den > 0 ? (S[i * n + j] / den) : 0;
        if (c > 1) c = 1;
        else if (c < -1) c = -1;
      }
      matrix[i * n + j] = c;
      matrix[j * n + i] = c;
    }
  }
  return { n, matrix, variances, nFrames: frames.length };
}

/**
 * Strongest off-diagonal correlations, sorted by |C|.
 * @param {DCCMResult|{n:number, matrix:Float64Array}} dcc
 * @param {number} [top=8]
 * @returns {Array<{i:number, j:number, c:number}>}
 */
export function topCorrelations(dcc, top = 8) {
  const { n, matrix } = dcc;
  const out = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) out.push({ i, j, c: matrix[i * n + j] });
  }
  out.sort((a, b) => Math.abs(b.c) - Math.abs(a.c));
  return out.slice(0, top);
}

/** Diverging blue → slate → red colormap for C ∈ [vmin, vmax]. */
function dccColor(v, vmin, vmax) {
  const t = Math.min(1, Math.max(0, (v - vmin) / Math.max(1e-12, vmax - vmin)));
  // 0 → deep blue, 0.5 → slate background, 1 → warm red
  const lo = [37, 99, 235], mid = [30, 41, 59], hi = [239, 68, 68];
  let r, g, b;
  if (t < 0.5) {
    const u = t / 0.5;
    r = lo[0] + (mid[0] - lo[0]) * u;
    g = lo[1] + (mid[1] - lo[1]) * u;
    b = lo[2] + (mid[2] - lo[2]) * u;
  } else {
    const u = (t - 0.5) / 0.5;
    r = mid[0] + (hi[0] - mid[0]) * u;
    g = mid[1] + (hi[1] - mid[1]) * u;
    b = mid[2] + (hi[2] - mid[2]) * u;
  }
  return `rgb(${Math.round(r)},${Math.round(g)},${Math.round(b)})`;
}

/**
 * Paint the DCCM onto a compact canvas (square, cell = size/n). Stores the
 * matrix on canvas._dccm for click→(i,j) picking by the panel.
 *
 * @param {HTMLCanvasElement|null} canvas
 * @param {Float64Array} matrix  row-major n×n
 * @param {object} [opts]
 * @param {number} [opts.n]      inferred from sqrt(matrix.length) if omitted
 * @param {number} [opts.size=220] CSS px (backing store honors devicePixelRatio)
 * @param {number} [opts.vmin=-1]
 * @param {number} [opts.vmax=1]
 * @returns {HTMLCanvasElement|null} the canvas, or null when unavailable (Node)
 */
export function renderDCCMHeatmap(canvas, matrix, opts = {}) {
  if (!canvas || typeof canvas.getContext !== "function" || !matrix) return null;
  const n = opts.n ?? Math.round(Math.sqrt(matrix.length));
  if (!Number.isFinite(n) || n < 2 || matrix.length < n * n) return null;
  const size = opts.size ?? 220;
  const vmin = opts.vmin ?? -1, vmax = opts.vmax ?? 1;
  const dpr = Math.min(2, (typeof window !== "undefined" && window.devicePixelRatio) || 1);
  canvas.width = Math.round(size * dpr);
  canvas.height = Math.round(size * dpr);
  if (canvas.style) { canvas.style.width = `${size}px`; canvas.style.height = `${size}px`; }
  const ctx = canvas.getContext("2d");
  if (!ctx) return null;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const cell = size / n;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      ctx.fillStyle = dccColor(matrix[i * n + j], vmin, vmax);
      ctx.fillRect(j * cell, i * cell, Math.ceil(cell * 10) / 10, Math.ceil(cell * 10) / 10);
    }
  }
  canvas._dccm = { n, matrix, size };
  return canvas;
}

/**
 * Map a click on a heatmap canvas to a residue pair.
 * @param {HTMLCanvasElement} canvas  previously painted by renderDCCMHeatmap
 * @param {number} clientX
 * @param {number} clientY
 * @returns {{i:number, j:number}|null}
 */
export function dccmPick(canvas, clientX, clientY) {
  try {
    const meta = canvas && canvas._dccm;
    if (!meta || typeof canvas.getBoundingClientRect !== "function") return null;
    const rect = canvas.getBoundingClientRect();
    const x = (clientX - rect.left) / Math.max(1e-9, rect.width);
    const y = (clientY - rect.top) / Math.max(1e-9, rect.height);
    if (x < 0 || x > 1 || y < 0 || y > 1) return null;
    const j = Math.min(meta.n - 1, Math.floor(x * meta.n));
    const i = Math.min(meta.n - 1, Math.floor(y * meta.n));
    return { i, j };
  } catch (_) {
    return null;
  }
}

/**
 * Viewer hook: mark residues (i, j) as the correlated pair of interest.
 * Idempotently wraps viewer.render() so two highlight rings are drawn at the
 * residues' projected positions after the normal scene. Duck-typed and
 * exception-safe — a viewer without canvas/projection support just records
 * the pair and returns it.
 *
 * @param {object} viewer  Viewer instance (or any {render, _px, _py} object)
 * @param {number} i
 * @param {number} j
 * @param {object} [opts]
 * @param {string} [opts.color="#fbbf24"]
 * @returns {{i:number, j:number}}
 */
export function highlightCorrelatedPair(viewer, i, j, opts = {}) {
  const color = opts.color ?? "#fbbf24";
  if (viewer && typeof viewer === "object") {
    try {
      viewer._dccmHighlight = { i, j, color };
      if (typeof viewer.render === "function" && !viewer._dccmWrapped) {
        const base = viewer.render.bind(viewer);
        viewer.render = function dccmHighlightRender(pos) {
          base(pos);
          try {
            const hl = viewer._dccmHighlight;
            const ctx = viewer.ctx;
            if (!hl || !ctx || !viewer._px || !viewer._py) return;
            const dpr = viewer.dpr || 1;
            ctx.save();
            ctx.strokeStyle = hl.color;
            ctx.lineWidth = 2 * dpr;
            for (const k of [hl.i, hl.j]) {
              if (k < 0 || k >= viewer.n) continue;
              ctx.beginPath();
              ctx.arc(viewer._px[k], viewer._py[k], 9 * dpr, 0, 2 * Math.PI);
              ctx.stroke();
            }
            ctx.restore();
          } catch (_) { /* overlay must never break rendering */ }
        };
        viewer._dccmWrapped = true;
      }
      if (typeof viewer.render === "function" && viewer._lastPos && viewer.canvas) {
        viewer.render(viewer._lastPos);
      }
    } catch (_) { /* record-only fallback */ }
  }
  return { i, j };
}

/**
 * Clear the viewer correlation highlight (keeps the idempotent wrapper).
 * @param {object} viewer
 */
export function clearCorrelatedPair(viewer) {
  try {
    if (viewer && typeof viewer === "object") {
      viewer._dccmHighlight = null;
      if (typeof viewer.render === "function" && viewer._lastPos && viewer.canvas) {
        viewer.render(viewer._lastPos);
      }
    }
  } catch (_) {}
}
