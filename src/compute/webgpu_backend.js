/**
 * webgpu_backend.js — WebGPU WGSL orchestrator for GPU non-bonded forces.
 *
 * Owns the full GPU dispatch sequence for Phase 3:
 *   1. upload unified buffers (positions + LJ/charge/Born params),
 *   2. cell_list.wgsl `bin_count` -> read back counts -> CPU prefix-sum ->
 *      cell_list.wgsl `bin_fill` (counting-sort neighbor list),
 *   3. nonbonded_forces.wgsl tiled LJ + screened Coulomb + Still-GB kernel,
 *   4. zero-copy-ish readback of forces+per-atom energies via mapped buffers.
 *
 * Unified ArrayBuffer layout (all Float32, 16 B stride per atom):
 *   positions : [x, y, z, w] x N   (w unused padding, Angstrom)
 *   params    : [sigma, eps, q, bornR] x N  (A, kcal/mol, e, A)
 *   forces    : [fx, fy, fz, eAtom] x N    (kcal/mol/A, kcal/mol)
 * Cell arrays are Uint32: atomCell[N], cellStart[numCells+1], sortedIdx[N].
 *
 * Fallback chain (all inside try/catch, never throws for missing hardware):
 *   WebGPU -> opts.workerPool.computeParallel -> cpuNonbondedReference.
 * Headless Node (no navigator.gpu) lands on the CPU reference, which mirrors
 * the WGSL kernel term-by-term (same switch polynomial, Still-GB chain rule,
 * sr clamp, exclusion keys) so parity can be validated without a GPU.
 *
 * Units: A, kcal/mol, e. Zero new npm deps. Vanilla ES module.
 *
 * @module compute/webgpu_backend
 */

import { SpatialGrid } from "../spatial-grid.js?v=10";

/** Default non-bonded cutoff == GPU cell edge (A). Matches gb_obc2 HCT cutoff. */
export const RCUT_DEFAULT = 12.0;
/** Default LJ switch start (A). */
export const SWITCH_ON_DEFAULT = 10.0;
/** Coulomb constant (kcal·A/mol/e^2). */
export const COULOMB_CONST = 332.06371;
/** Threads per workgroup (must match both .wgsl files). */
export const WORKGROUP = 64;
/** Max dense-grid cells per axis (caps VRAM for pathological boxes). */
export const MAX_GRID_DIM = 64;

// ---------------------------------------------------------------------------
// Module-singleton GPU state
// ---------------------------------------------------------------------------

const _gpu = {
  adapter: null,
  device: null,
  pipelines: null, // { count, fill, nonbonded }
  wgsl: null, // { cellList, nonbonded }
  ready: false,
  reason: "not initialized",
  bufs: null, // persistent buffers, see _ensureBuffers
  bufN: 0,
  bufCells: 0,
  lastExclSig: "",
  lastScaleSig: "",
};

/**
 * Synchronous capability probe. True only where navigator.gpu exists.
 * @returns {boolean}
 */
export function isSupported() {
  try {
    return typeof navigator !== "undefined" && !!navigator.gpu;
  } catch (_) {
    return false;
  }
}

/**
 * Load a WGSL source file sibling to this module. Browser: fetch();
 * Node: node:fs fallback. Never throws — returns null on failure.
 * @param {string} name file name under ./wgsl/
 * @returns {Promise<string|null>}
 */
async function _loadWgslText(name) {
  const url = new URL(`./wgsl/${name}`, import.meta.url);
  try {
    if (typeof fetch !== "undefined") {
      const res = await fetch(url);
      if (res && res.ok) return await res.text();
    }
  } catch (_) { /* fall through to fs */ }
  try {
    const fs = await import("node:fs/promises");
    const u = await import("node:url");
    return await fs.readFile(u.fileURLToPath(url), "utf-8");
  } catch (_) {
    return null;
  }
}

/**
 * Initialize the WebGPU device and compile the Phase 3 pipelines.
 * Safe to call anywhere: returns false (with reason) when WebGPU is
 * unavailable instead of throwing.
 * @param {object} [opts]
 * @param {boolean} [opts.force=false] re-initialize even if ready
 * @returns {Promise<boolean>} true when the GPU path is usable
 */
export async function initWebGPU(opts = {}) {
  if (_gpu.ready && !opts.force) return true;
  if (!isSupported()) {
    _gpu.ready = false;
    _gpu.reason = "navigator.gpu unavailable (headless or unsupported browser)";
    return false;
  }
  try {
    const [cellList, nonbonded] = await Promise.all([
      _loadWgslText("cell_list.wgsl"),
      _loadWgslText("nonbonded_forces.wgsl"),
    ]);
    if (!cellList || !nonbonded) {
      _gpu.reason = "WGSL sources failed to load";
      _gpu.ready = false;
      return false;
    }
    // Balanced-braces sanity before handing to the driver (cheap early error).
    for (const [nm, src] of [["cell_list", cellList], ["nonbonded", nonbonded]]) {
      const open = (src.match(/{/g) || []).length;
      const close = (src.match(/}/g) || []).length;
      if (open !== close || !/entry|@compute/.test(src)) {
        _gpu.reason = `WGSL sanity failed for ${nm} ({${open}} vs }{${close}})`;
        _gpu.ready = false;
        return false;
      }
    }
    _gpu.wgsl = { cellList, nonbonded };
    _gpu.adapter = await navigator.gpu.requestAdapter();
    if (!_gpu.adapter) {
      _gpu.reason = "no GPU adapter found";
      _gpu.ready = false;
      return false;
    }
    _gpu.device = await _gpu.adapter.requestDevice();
    const mk = (code, entry) => {
      const module = _gpu.device.createShaderModule({ code });
      return _gpu.device.createComputePipeline({
        layout: "auto",
        compute: { module, entryPoint: entry },
      });
    };
    _gpu.pipelines = {
      count: mk(cellList, "bin_count"),
      fill: mk(cellList, "bin_fill"),
      nonbonded: mk(nonbonded, "main"),
    };
    _gpu.bufs = null;
    _gpu.bufN = 0;
    _gpu.bufCells = 0;
    _gpu.ready = true;
    _gpu.reason = "ready";
    return true;
  } catch (e) {
    _gpu.ready = false;
    _gpu.reason = `init failed: ${e?.message ?? e}`;
    return false;
  }
}

/**
 * Why the GPU path is (un)available. Useful for the settings status line.
 * @returns {string}
 */
export function webgpuStatus() {
  return _gpu.reason;
}

// ---------------------------------------------------------------------------
// Grid geometry (CPU; shared by GPU + CPU-reference paths)
// ---------------------------------------------------------------------------

/**
 * Dense-grid geometry covering the positions with cells >= cellSize.
 * @param {Float32Array|Float64Array} pos4 flat 4n [x,y,z,w] or 3n
 * @param {number} n atom count
 * @param {number} stride 4 or 3
 * @param {number} cellSize cell edge (A)
 * @returns {{min: number[], dim: number[], numCells: number}}
 */
export function planGrid(pos4, n, stride, cellSize) {
  let mnx = Infinity, mny = Infinity, mnz = Infinity;
  let mxx = -Infinity, mxy = -Infinity, mxz = -Infinity;
  for (let i = 0; i < n; i++) {
    const x = pos4[stride * i], y = pos4[stride * i + 1], z = pos4[stride * i + 2];
    if (x < mnx) mnx = x; if (y < mny) mny = y; if (z < mnz) mnz = z;
    if (x > mxx) mxx = x; if (y > mxy) mxy = y; if (z > mxz) mxz = z;
  }
  const pad = cellSize + 1e-3;
  const dim = [0, 1, 2].map((a) => {
    const lo = [mnx, mny, mnz][a] - pad;
    const hi = [mxx, mxy, mxz][a] + pad;
    return Math.max(1, Math.min(MAX_GRID_DIM, Math.ceil((hi - lo) / cellSize)));
  });
  // Re-derive min so the box is centered on the span (matches clamp rule).
  const min = [mnx - pad, mny - pad, mnz - pad];
  return { min, dim, numCells: dim[0] * dim[1] * dim[2] };
}

// ---------------------------------------------------------------------------
// Persistent GPU buffers
// ---------------------------------------------------------------------------

function _createBuf(size, usage) {
  return _gpu.device.createBuffer({ size: Math.max(16, size), usage });
}

function _ensureBuffers(n, numCells, exclCount, scaleCount) {
  const B = _gpu.bufs;
  const S = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST;
  const need = !B || _gpu.bufN < n || _gpu.bufCells < numCells;
  if (need) {
    if (B) for (const k of Object.keys(B)) { try { B[k].destroy(); } catch (_) {} }
    const R = GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST;
    _gpu.bufs = {
      pos: _createBuf(n * 16, S),
      prm: _createBuf(n * 16, S),
      atomCell: _createBuf(n * 4, S | GPUBufferUsage.COPY_SRC),
      counts: _createBuf(numCells * 4, S | GPUBufferUsage.COPY_SRC),
      countsRead: _createBuf(numCells * 4, R),
      start: _createBuf((numCells + 1) * 4, S),
      cursor: _createBuf(numCells * 4, S),
      sorted: _createBuf(n * 4, S),
      frc: _createBuf(n * 16, S | GPUBufferUsage.COPY_SRC),
      frcRead: _createBuf(n * 16, R),
      uniCell: _createBuf(48, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST),
      uniNB: _createBuf(80, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST),
      excl: _createBuf(Math.max(1, exclCount) * 4, S),
      scaled: _createBuf(Math.max(1, scaleCount) * 4, S),
    };
    _gpu.bufN = n;
    _gpu.bufCells = numCells;
    _gpu.lastExclSig = "";
    _gpu.lastScaleSig = "";
  } else if (exclCount > 0 && _gpu.bufs.excl.size < exclCount * 4) {
    try { _gpu.bufs.excl.destroy(); } catch (_) {}
    _gpu.bufs.excl = _createBuf(exclCount * 4, S);
    _gpu.lastExclSig = "";
  }
}

function _bind(pipeline, entries) {
  return _gpu.device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries,
  });
}

async function _readU32(srcBuf, readBuf, count) {
  const enc = _gpu.device.createCommandEncoder();
  enc.copyBufferToBuffer(srcBuf, 0, readBuf, 0, count * 4);
  _gpu.device.queue.submit([enc.finish()]);
  await readBuf.mapAsync(GPUMapMode.READ);
  const out = new Uint32Array(readBuf.getMappedRange().slice(0));
  readBuf.unmap();
  return out;
}

async function _readF32(srcBuf, readBuf, count) {
  const enc = _gpu.device.createCommandEncoder();
  enc.copyBufferToBuffer(srcBuf, 0, readBuf, 0, count * 4);
  _gpu.device.queue.submit([enc.finish()]);
  await readBuf.mapAsync(GPUMapMode.READ);
  const out = new Float32Array(readBuf.getMappedRange().slice(0));
  readBuf.unmap();
  return out;
}

// ---------------------------------------------------------------------------
// GPU dispatch
// ---------------------------------------------------------------------------

/**
 * Normalize caller positions to the unified Float32 4n layout.
 * @param {Float32Array|Float64Array} positions flat 4n [x,y,z,w] or 3n
 * @param {number} n
 * @returns {{pos4: Float32Array, stride: number}}
 */
export function toUnifiedPositions(positions, n) {
  if (positions.length === 4 * n) {
    return {
      pos4: positions instanceof Float32Array ? positions : Float32Array.from(positions),
      stride: 4,
    };
  }
  if (positions.length === 3 * n) {
    const pos4 = new Float32Array(4 * n);
    for (let i = 0; i < n; i++) {
      pos4[4 * i] = positions[3 * i];
      pos4[4 * i + 1] = positions[3 * i + 1];
      pos4[4 * i + 2] = positions[3 * i + 2];
    }
    return { pos4, stride: 4 };
  }
  throw new Error(`webgpu_backend: positions length ${positions.length} != 4n or 3n (n=${n})`);
}

/**
 * Build the unified Float32 4n params table from per-atom records.
 * @param {Array<{sigma: number, eps: number, q: number}>} elem length n
 * @param {ArrayLike<number>|null} [bornRadii] length n (0/omitted -> 1.7 A default, clamped >= 0.5 in kernel)
 * @returns {Float32Array} [sigma, eps, q, bornR] x N
 */
export function buildParamsTable(elem, bornRadii = null) {
  const n = elem.length;
  const out = new Float32Array(4 * n);
  for (let i = 0; i < n; i++) {
    const e = elem[i];
    out[4 * i] = e.sigma ?? 3.4;
    out[4 * i + 1] = e.eps ?? 0.12;
    out[4 * i + 2] = e.q ?? 0;
    out[4 * i + 3] = bornRadii ? (bornRadii[i] || 1.7) : 1.7;
  }
  return out;
}

/**
 * GPU pass: bin_count -> CPU prefix-sum -> bin_fill -> nonbonded kernel.
 * Caller must have awaited initWebGPU() === true.
 * @param {Float32Array} pos4 unified positions (4n)
 * @param {Float32Array} prm4 unified params (4n)
 * @param {object} o physics + exclusion options (see dispatchNonbonded)
 * @returns {Promise<{forces: Float32Array, numCells: number, gridDim: number[]}>}
 */
async function _dispatchGPU(pos4, prm4, o) {
  const { device, pipelines, bufs } = _gpu;
  const { n, cutOff, switchOn, epsIn, epsOut, kappa, useGB, excl, scaled } = o;
  const grid = planGrid(pos4, n, 4, cutOff);
  const { min, dim, numCells } = grid;

  _ensureBuffers(n, numCells, excl.length, scaled.length);
  const B = _gpu.bufs;
  const q = device.queue;
  const groups = Math.ceil(n / WORKGROUP);

  q.writeBuffer(B.pos, 0, pos4);
  q.writeBuffer(B.prm, 0, prm4);

  // Uniforms: CellParams (48 B).
  const uc = new ArrayBuffer(48);
  const ucf = new Float32Array(uc);
  const ucu = new Uint32Array(uc);
  ucf[0] = min[0]; ucf[1] = min[1]; ucf[2] = min[2]; ucf[3] = cutOff;
  ucu[4] = dim[0]; ucu[5] = dim[1]; ucu[6] = dim[2]; ucu[7] = n;
  ucu[8] = numCells; ucu[9] = 0; ucu[10] = 0; ucu[11] = 0;
  q.writeBuffer(B.uniCell, 0, uc);

  // Pass 1: zero counts, dispatch bin_count.
  let enc = device.createCommandEncoder();
  enc.clearBuffer(B.counts, 0, numCells * 4);
  {
    const pass = enc.beginComputePass();
    pass.setPipeline(pipelines.count);
    pass.setBindGroup(0, _bind(pipelines.count, [
      { binding: 0, resource: { buffer: B.pos } },
      { binding: 1, resource: { buffer: B.uniCell } },
      { binding: 2, resource: { buffer: B.atomCell } },
      { binding: 3, resource: { buffer: B.counts } },
    ]));
    pass.dispatchWorkgroups(groups);
    pass.end();
  }
  q.submit([enc.finish()]);
  const counts = await _readU32(B.counts, B.countsRead, numCells);

  // CPU prefix-sum -> cellStart[numCells+1] (sentinel = n).
  const start = new Uint32Array(numCells + 1);
  let acc = 0;
  for (let c = 0; c < numCells; c++) { start[c] = acc; acc += counts[c]; }
  start[numCells] = n;
  q.writeBuffer(B.start, 0, start);

  // Pass 2: zero cursor, dispatch bin_fill.
  enc = device.createCommandEncoder();
  enc.clearBuffer(B.cursor, 0, numCells * 4);
  {
    const pass = enc.beginComputePass();
    pass.setPipeline(pipelines.fill);
    // NOTE: bin_fill only reads bindings 0,1,2 + 4,5,6; binding 3 is
    // inactive for this entry point. With layout:"auto" the bind group must
    // contain exactly the *active* bindings, so counts is omitted here.
    pass.setBindGroup(0, _bind(pipelines.fill, [
      { binding: 0, resource: { buffer: B.pos } },
      { binding: 1, resource: { buffer: B.uniCell } },
      { binding: 2, resource: { buffer: B.atomCell } },
      { binding: 4, resource: { buffer: B.start } },
      { binding: 5, resource: { buffer: B.cursor } },
      { binding: 6, resource: { buffer: B.sorted } },
    ]));
    pass.dispatchWorkgroups(groups);
    pass.end();
  }

  // Exclusion tables (re-upload only when contents change).
  const sigOf = (a) => `${a.length}:${a.length ? a[0] : 0}:${a.length ? a[a.length - 1] : 0}`;
  if (sigOf(excl) !== _gpu.lastExclSig && excl.length) {
    q.writeBuffer(B.excl, 0, excl);
    _gpu.lastExclSig = sigOf(excl);
  }
  if (sigOf(scaled) !== _gpu.lastScaleSig && scaled.length) {
    q.writeBuffer(B.scaled, 0, scaled);
    _gpu.lastScaleSig = sigOf(scaled);
  }

  // Uniforms: NBParams (80 B).
  const un = new ArrayBuffer(80);
  const unf = new Float32Array(un);
  const unu = new Uint32Array(un);
  unu[0] = dim[0]; unu[1] = dim[1]; unu[2] = dim[2]; unu[3] = n;
  unf[4] = min[0]; unf[5] = min[1]; unf[6] = min[2]; unf[7] = cutOff;
  unf[8] = cutOff; unf[9] = switchOn; unf[10] = COULOMB_CONST; unf[11] = epsIn;
  unf[12] = epsOut; unf[13] = kappa;
  unu[14] = useGB ? 1 : 0; unu[15] = excl.length;
  unu[16] = scaled.length; unu[17] = 0; unu[18] = 0; unu[19] = 0;
  q.writeBuffer(B.uniNB, 0, un);

  // Pass 3: nonbonded forces (appended to the same encoder as pass 2).
  {
    const pass = enc.beginComputePass();
    pass.setPipeline(pipelines.nonbonded);
    pass.setBindGroup(0, _bind(pipelines.nonbonded, [
      { binding: 0, resource: { buffer: B.pos } },
      { binding: 1, resource: { buffer: B.prm } },
      { binding: 2, resource: { buffer: B.atomCell } },
      { binding: 3, resource: { buffer: B.start } },
      { binding: 4, resource: { buffer: B.sorted } },
      { binding: 5, resource: { buffer: B.frc } },
      { binding: 6, resource: { buffer: B.uniNB } },
      { binding: 7, resource: { buffer: B.excl } },
      { binding: 8, resource: { buffer: B.scaled } },
    ]));
    pass.dispatchWorkgroups(groups);
    pass.end();
  }
  q.submit([enc.finish()]);
  const forces = await _readF32(B.frc, B.frcRead, 4 * n);
  return { forces, numCells, gridDim: dim };
}

// ---------------------------------------------------------------------------
// CPU reference (term-by-term mirror of nonbonded_forces.wgsl)
// ---------------------------------------------------------------------------

function _bsearch(arr, key) {
  let lo = 0, hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    const v = arr[mid];
    if (v === key) return true;
    if (v < key) lo = mid + 1; else hi = mid;
  }
  return false;
}

function _switchS(r, on, off) {
  if (r <= on) return 1;
  if (r >= off) return 0;
  const rsq = r * r, on2 = on * on, cut2 = off * off;
  const denom = (cut2 - on2) ** 3;
  return ((cut2 - rsq) ** 2) * (cut2 + 2 * rsq - 3 * on2) / denom;
}

function _switchDS(r, on, off) {
  if (r <= on || r >= off) return 0;
  const rsq = r * r, on2 = on * on, cut2 = off * off;
  const denom = (cut2 - on2) ** 3;
  const a = (cut2 - rsq) * (cut2 - rsq);
  const b = cut2 + 2 * rsq - 3 * on2;
  return ((-4 * r * (cut2 - rsq)) * b + a * (4 * r)) / denom;
}

/**
 * Single-core CPU mirror of the WGSL non-bonded kernel (fallback + parity
 * reference). Same formulas, same exclusion keys, same per-atom energy in
 * w. Uses SpatialGrid so the neighbor definition matches the GPU cell list.
 *
 * @param {Float32Array} pos4 unified positions (4n)
 * @param {Float32Array} prm4 unified params (4n)
 * @param {object} o { n, cutOff, switchOn, epsIn, epsOut, kappa, useGB, excl, scaled }
 * @returns {{forces: Float32Array, backend: string}}
 */
export function cpuNonbondedReference(pos4, prm4, o) {
  const { n, cutOff, switchOn, epsIn, epsOut, kappa, useGB, excl, scaled } = o;
  // SpatialGrid works on flat 3n Float64; repack once (fallback path).
  const pos3 = new Float64Array(3 * n);
  for (let i = 0; i < n; i++) {
    pos3[3 * i] = pos4[4 * i];
    pos3[3 * i + 1] = pos4[4 * i + 1];
    pos3[3 * i + 2] = pos4[4 * i + 2];
  }
  const grid = new SpatialGrid(cutOff, n);
  grid.build(pos3, n);
  const F = new Float32Array(4 * n);
  const invIn = 1 / epsIn, invOut = 1 / epsOut;
  grid.forEachPair(pos3, n, cutOff, (i, j, dx, dy, dz, r2, r) => {
    if (r2 < 1e-6) return;
    const key = i < j ? i * 1e6 + j : j * 1e6 + i;
    if (_bsearch(excl, key)) return;
    const s14 = _bsearch(scaled, key) ? 0.5 : 1.0;
    const sgi = prm4[4 * i], epi = prm4[4 * i + 1], qi = prm4[4 * i + 2];
    const Ri = Math.max(prm4[4 * i + 3] || 1.7, 0.5);
    const sgj = prm4[4 * j], epj = prm4[4 * j + 1], qj = prm4[4 * j + 2];
    const Rj = Math.max(prm4[4 * j + 3] || 1.7, 0.5);
    const s = 0.5 * (sgi + sgj);
    const eps = Math.sqrt(Math.max(epi, 0) * Math.max(epj, 0));
    const sr = Math.min(s / r, 5);
    const sr6 = sr ** 6;
    const ljE = 4 * eps * (sr6 * sr6 - sr6);
    const ljF = (4 * eps * (12 * sr6 * sr6 - 6 * sr6)) / r;
    const S = _switchS(r, switchOn, cutOff);
    const dS = _switchDS(r, switchOn, cutOff);
    let dudr = s14 * (S * -ljF + dS * ljE);
    let ePair = s14 * S * ljE;
    if (qi !== 0 && qj !== 0) {
      const qq = qi * qj * s14;
      const uC = ((COULOMB_CONST / epsIn) * qq) / r;
      let dudrE = -uC / r;
      ePair += uC;
      if (useGB) {
        const a = Math.max(Ri * Rj, 1e-6);
        const e = Math.exp(-r2 / (4 * a));
        const f2 = r2 + a * e;
        const f = Math.sqrt(Math.max(f2, 1e-12));
        const eK = Math.exp(-kappa * f);
        const P = invIn - eK * invOut;
        ePair += -COULOMB_CONST * qq * (P / f);
        const dfdr = (r * (1 - 0.25 * e)) / f;
        const dgdf = ((kappa * eK * invOut) * f - P) / f2;
        dudrE += -COULOMB_CONST * qq * dgdf * dfdr;
      }
      dudr += dudrE;
    }
    const fmag = dudr / r;
    const fx = fmag * dx, fy = fmag * dy, fz = fmag * dz;
    // dx points i -> j (SpatialGrid callback convention matches kernel d).
    F[4 * i] += fx; F[4 * i + 1] += fy; F[4 * i + 2] += fz; F[4 * i + 3] += 0.5 * ePair;
    F[4 * j] -= fx; F[4 * j + 1] -= fy; F[4 * j + 2] -= fz; F[4 * j + 3] += 0.5 * ePair;
  });
  if (useGB) {
    for (let i = 0; i < n; i++) {
      const q = prm4[4 * i + 2];
      if (q === 0) continue;
      const R = Math.max(prm4[4 * i + 3] || 1.7, 0.5);
      const eK = Math.exp(-kappa * R);
      F[4 * i + 3] += -0.5 * COULOMB_CONST * q * q * (invIn - eK * invOut) / R;
    }
  }
  return { forces: F, backend: "cpu" };
}

// ---------------------------------------------------------------------------
// Public dispatch with fallback chain
// ---------------------------------------------------------------------------

/**
 * Evaluate non-bonded forces, preferring WebGPU and falling back through
 * worker-pool to the single-core reference. Never throws for missing
 * hardware — the returned `backend` field reports which path ran.
 *
 * @param {Float32Array|Float64Array} positions flat 4n [x,y,z,w] or 3n
 * @param {Float32Array|Array<{sigma:number,eps:number,q:number}>} params
 *   unified 4n table (use buildParamsTable) or per-atom records
 * @param {object} [opts]
 * @param {number} opts.n atom count (required)
 * @param {number} [opts.cutOff=12] pair cutoff / cell edge (A)
 * @param {number} [opts.switchOn=10] LJ switch start (A)
 * @param {number} [opts.epsIn=4] solute dielectric
 * @param {number} [opts.epsOut=78.5] solvent dielectric
 * @param {number} [opts.kappa=0.127] inverse Debye length (A^-1, ~150 mM)
 * @param {boolean} [opts.useGB=true] Still-GB reaction field on/off
 * @param {ArrayLike<number>|null} [opts.bornRadii] used when params are records
 * @param {Uint32Array} [opts.excluded] sorted skip keys (i*1e6+j)
 * @param {Uint32Array} [opts.scaled14] sorted 1-4 keys (s14 = 0.5)
 * @param {object|null} [opts.workerPool] optional WorkerPool (throughput fallback)
 * @returns {Promise<{forces: Float32Array, backend: string, numCells: number}>}
 */
export async function dispatchNonbonded(positions, params, opts = {}) {
  const n = opts.n;
  if (!Number.isInteger(n) || n <= 0) throw new Error("dispatchNonbonded: opts.n required");
  const cutOff = opts.cutOff ?? RCUT_DEFAULT;
  const switchOn = opts.switchOn ?? SWITCH_ON_DEFAULT;
  const o = {
    n,
    cutOff,
    switchOn,
    epsIn: opts.epsIn ?? 4.0,
    epsOut: opts.epsOut ?? 78.5,
    kappa: opts.kappa ?? 0.127,
    useGB: opts.useGB ?? true,
    excl: opts.excluded ?? new Uint32Array(0),
    scaled: opts.scaled14 ?? new Uint32Array(0),
  };
  const { pos4 } = toUnifiedPositions(positions, n);
  const prm4 = (params instanceof Float32Array && params.length === 4 * n)
    ? params
    : buildParamsTable(params, opts.bornRadii ?? null);

  // 1. GPU path.
  if (_gpu.ready) {
    try {
      const r = await _dispatchGPU(pos4, prm4, o);
      return { forces: r.forces, backend: "webgpu", numCells: r.numCells };
    } catch (e) {
      console.warn(`[webgpu_backend] GPU dispatch failed (${e?.message ?? e}) — CPU fallback`);
    }
  }
  // 2. Worker-pool throughput fallback (caller-owned pool, guarded n > 800
  //    like src/main.js useWorkerPoolIfNeeded; physics = worker kernel).
  if (opts.workerPool) {
    try {
      const pos64 = new Float64Array(3 * n);
      for (let i = 0; i < n; i++) {
        pos64[3 * i] = pos4[4 * i];
        pos64[3 * i + 1] = pos4[4 * i + 1];
        pos64[3 * i + 2] = pos4[4 * i + 2];
      }
      const par = await opts.workerPool.computeParallel(pos64, n);
      if (par && par.forces) {
        const F = new Float32Array(4 * n);
        for (let i = 0; i < n; i++) {
          F[4 * i] = par.forces[3 * i];
          F[4 * i + 1] = par.forces[3 * i + 1];
          F[4 * i + 2] = par.forces[3 * i + 2];
        }
        return { forces: F, backend: "worker", numCells: 0 };
      }
    } catch (e) {
      console.warn(`[webgpu_backend] worker fallback failed (${e?.message ?? e}) — single-core fallback`);
    }
  }
  // 3. Single-core reference (always available, kernel-parity physics).
  try {
    const r = cpuNonbondedReference(pos4, prm4, o);
    return { forces: r.forces, backend: "cpu", numCells: 0 };
  } catch (e) {
    console.warn(`[webgpu_backend] CPU reference failed (${e?.message ?? e}) — zero forces`);
    return { forces: new Float32Array(4 * n), backend: "none", numCells: 0 };
  }
}

/**
 * Release GPU buffers/device (page hide / tests). Safe when uninitialized.
 */
export function disposeWebGPU() {
  try {
    if (_gpu.bufs) for (const k of Object.keys(_gpu.bufs)) {
      try { _gpu.bufs[k].destroy(); } catch (_) {}
    }
    if (_gpu.device) { try { _gpu.device.destroy(); } catch (_) {} }
  } catch (_) { /* never throw from dispose */ }
  _gpu.bufs = null;
  _gpu.device = null;
  _gpu.adapter = null;
  _gpu.pipelines = null;
  _gpu.ready = false;
  _gpu.reason = "disposed";
}
