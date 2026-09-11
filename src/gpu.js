/**
 * gpu.js — WebGPU compute shader pipeline for GPU-accelerated non-bonded forces.
 *
 * Implements GPU dispatch using WebGPU WGSL compute shaders when supported by the browser,
 * with automatic fallback detection.
 * G64 — GPU vs CPU correlation R>0.999 is aspirational, not yet benchmarked (placeholder target)
 */

const WGSL_NONBONDED = `
struct AtomParams {
  sigma: f32,
  eps: f32,
  q: f32,
  padding: f32,
};

@group(0) @binding(0) var<storage, read> pos : array<vec4<f32>>;
@group(0) @binding(1) var<storage, read> params : array<AtomParams>;
@group(0) @binding(2) var<storage, read_write> forces : array<vec4<f32>>;
@group(0) @binding(3) var<uniform> numAtoms : u32;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) global_id : vec3<u32>) {
  let i = global_id.x;
  let n = numAtoms;
  if (i >= n) { return; }

  let p_i = pos[i].xyz;
  let param_i = params[i];
  var f_i = vec3<f32>(0.0, 0.0, 0.0);

  let R_CUT = 8.5;
  let cut2 = R_CUT * R_CUT;

  for (var j = 0u; j < n; j = j + 1u) {
    if (i == j) { continue; }

    let p_j = pos[j].xyz;
    let dx = p_j.x - p_i.x;
    let dy = p_j.y - p_i.y;
    let dz = p_j.z - p_i.z;
    let r2 = dx * dx + dy * dy + dz * dz;

    if (r2 < cut2 && r2 > 0.0001) {
      let r = sqrt(r2);
      let param_j = params[j];

      // LJ — G64 clamp sr<5 to avoid overflow at near-zero r
      let sigma = 0.5 * (param_i.sigma + param_j.sigma);
      let eps = sqrt(param_i.eps * param_j.eps);
      let sr = min(sigma / r, 5.0); // clamp sigma/r to 5
      let sr6 = sr * sr * sr * sr * sr * sr;
      let ljF = 4.0 * eps * (12.0 * sr6 * sr6 - 6.0 * sr6) / r;

      // Screened Coulomb
      var ef = 0.0;
      if (param_i.q != 0.0 && param_j.q != 0.0) {
        let qq = param_i.q * param_j.q;
        let fac = exp(-r / 8.0) / (r * r);
        let ee = 332.0 * qq * fac;
        ef = ee * (-0.125 - 2.0 / r);
      }

      let totF = -ljF + ef;
      f_i.x = f_i.x + totF * (dx / r);
      f_i.y = f_i.y + totF * (dy / r);
      f_i.z = f_i.z + totF * (dz / r);
    }
  }

  forces[i] = vec4<f32>(f_i, 0.0);
}
`;

export class GpuAccelerator {
  constructor() {
    this.adapter = null;
    this.device = null;
    this.pipeline = null;
    this.isSupported = false;
    this.deviceInfo = "Not initialized";
  }

  async init() {
    if (typeof navigator === "undefined" || !navigator.gpu) {
      this.isSupported = false;
      this.deviceInfo = "WebGPU not supported in this environment";
      return false;
    }

    try {
      this.adapter = await navigator.gpu.requestAdapter();
      if (!this.adapter) {
        this.isSupported = false;
        this.deviceInfo = "No GPU adapter found";
        return false;
      }

      this.device = await this.adapter.requestDevice();
      const info = await this.adapter.requestAdapterInfo?.() || {};
      this.deviceInfo = `${info.vendor || "GPU"} ${info.architecture || ""} (${info.description || "WebGPU"})`;

      const module = this.device.createShaderModule({ code: WGSL_NONBONDED });
      this.pipeline = this.device.createComputePipeline({
        layout: "auto",
        compute: { module, entryPoint: "main" },
      });

      this.isSupported = true;
      return true;
    } catch (e) {
      this.isSupported = false;
      this.deviceInfo = `GPU Init error: ${e.message}`;
      return false;
    }
  }

  /**
   * Compute non-bonded forces on the GPU.
   *
   * @param {Float64Array} pos
   * @param {Array<{sigma: number, eps: number, q: number}>} elemParams
   * @returns {Promise<Float32Array|null>}
   */
  async computeForces(pos, elemParams) {
    if (!this.isSupported || !this.device || !this.pipeline) return null;

    const n = elemParams.length;
    const posF32 = new Float32Array(n * 4);
    const paramsF32 = new Float32Array(n * 4);

    for (let i = 0; i < n; i++) {
      posF32[4 * i] = pos[3 * i];
      posF32[4 * i + 1] = pos[3 * i + 1];
      posF32[4 * i + 2] = pos[3 * i + 2];
      posF32[4 * i + 3] = 0;

      const p = elemParams[i];
      paramsF32[4 * i] = p.sigma;
      paramsF32[4 * i + 1] = p.eps;
      paramsF32[4 * i + 2] = p.q;
      paramsF32[4 * i + 3] = 0;
    }

    // Allocate or reuse GPU buffers
    const numAtomsArray = new Uint32Array([n]);
    const numAtomsBuffer = this.device.createBuffer({
      size: 16,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(numAtomsBuffer, 0, numAtomsArray);

    const posBuffer = this.device.createBuffer({
      size: Math.max(16, posF32.byteLength),
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(posBuffer, 0, posF32);

    const paramsBuffer = this.device.createBuffer({
      size: Math.max(16, paramsF32.byteLength),
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(paramsBuffer, 0, paramsF32);

    const forcesBuffer = this.device.createBuffer({
      size: Math.max(16, n * 4 * 4),
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
    });

    const readbackBuffer = this.device.createBuffer({
      size: Math.max(16, n * 4 * 4),
      usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
    });

    const bindGroup = this.device.createBindGroup({
      layout: this.pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: posBuffer } },
        { binding: 1, resource: { buffer: paramsBuffer } },
        { binding: 2, resource: { buffer: forcesBuffer } },
        { binding: 3, resource: { buffer: numAtomsBuffer } },
      ],
    });

    const encoder = this.device.createCommandEncoder();
    const pass = encoder.beginComputePass();
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(Math.ceil(n / 64));
    pass.end();

    encoder.copyBufferToBuffer(forcesBuffer, 0, readbackBuffer, 0, n * 4 * 4);
    this.device.queue.submit([encoder.finish()]);

    await readbackBuffer.mapAsync(GPUMapMode.READ);
    const result = new Float32Array(readbackBuffer.getMappedRange().slice(0));
    readbackBuffer.unmap();

    // Clean up GPU buffers
    numAtomsBuffer.destroy();
    posBuffer.destroy();
    paramsBuffer.destroy();
    forcesBuffer.destroy();
    readbackBuffer.destroy();

    return result;
  }
}

// G64 — CPU fallback clamp helper for testing: sr = sigma/r clamped to 5
// G64 clamp sr<5 — aspirational GPU vs CPU R>0.999 validation target
export function gpuClampSr(sigma, r) {
  return Math.min(sigma / r, 5); // clamp sigma/r to 5
}
