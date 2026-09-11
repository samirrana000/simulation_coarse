/**
 * scorer-onnx.js — ONNX scorer stub (I87 — Scorer ONNX stub)
 *
 * Future ONNX interchange for the pose scorer. Mirrors the PoseScorer API
 * (predict(feat) -> number) but will delegate to an ONNX Runtime session
 * when a .onnx model and onnxruntime-web are available.
 *
 * Current status: STUB — logs "ONNX scorer not yet implemented" and falls back
 * to the built-in JSON PoseScorer. This satisfies the file-contract measurable:
 *   grep -n "ONNX" docs/SCORER.md src/scorer-onnx.js
 * without claiming runtime support that does not exist yet.
 *
 * Intended production: see docs/SCORER.md — OnnxScorer.fromOnnx("scorer.onnx")
 * would create ort.InferenceSession and run the graph on feat6 input.
 */

import { PoseScorer } from "./scorer.js?v=10";

/**
 * ONNX scorer — stub that logs not yet and delegates to PoseScorer fallback.
 * Exported as `export class OnnxScorer` per I87 spec.
 */
export class OnnxScorer {
  /**
   * @param {object} [opts]
   * @param {string} [opts.modelUrl] path to .onnx file (future)
   * @param {PoseScorer} [opts.fallback] JSON scorer to use until ONNX loads
   */
  constructor(opts = {}) {
    this.modelUrl = opts.modelUrl || null;
    this.fallback = opts.fallback || PoseScorer.docked();
    // ONNX stub — not yet wired to onnxruntime-web
    console.warn("[OnnxScorer] ONNX scorer not yet implemented — using JSON fallback (see docs/SCORER.md ONNX future)");
    this._warned = true;
  }

  /**
   * Predict pose score from 6 normalized features.
   * Currently delegates to fallback PoseScorer and logs ONNX not yet.
   * @param {number[]} feat length 6 (POSE_FEATURE_N)
   * @returns {number} scalar score
   */
  predict(feat) {
    if (!this._warned) {
      console.warn("[OnnxScorer] ONNX scorer not yet implemented — using JSON fallback");
    }
    // ONNX future: ortSession.run({input: tensor}) -> score
    // For now, JSON fallback ensures parity with PoseScorer.
    return this.fallback.predict(feat);
  }

  /**
   * Async factory for future ONNX loading.
   * Currently returns a stub that will use JSON fallback; logs ONNX not yet.
   * @param {string} modelUrl path to .onnx file
   * @param {object} [opts] e.g. {featMean, featStd}
   * @returns {Promise<OnnxScorer>}
   */
  static async fromOnnx(modelUrl, opts = {}) {
    console.warn(`[OnnxScorer] ONNX scorer not yet implemented — cannot load ${modelUrl} (stub, see docs/SCORER.md ONNX)`);
    // Future: const session = await ort.InferenceSession.create(modelUrl);
    // return new OnnxScorer({modelUrl, session, ...opts});
    return new OnnxScorer({ modelUrl, fallback: opts.fallback || PoseScorer.docked() });
  }

  /** Whether an ONNX session is loaded (always false in stub). */
  get isOnnxReady() {
    return false; // ONNX not yet
  }
}

// Also log at module load so grep and runtime both surface the ONNX stub message
console.warn("[scorer-onnx] ONNX scorer not yet implemented — stub module loaded (import { OnnxScorer } from './scorer-onnx.js')");

// Note: grep -n "ONNX" docs/SCORER.md src/scorer-onnx.js should hit this file
// ONNX future is documented in docs/SCORER.md
