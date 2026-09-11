/**
 * recorder.js — Trajectory capture and export.
 *
 * Frames are stored as compact Float32Array snapshots of the bead positions,
 * plus the simulation time of each snapshot, so PDB CONECT/REMARK info and
 * XYZ headers carry real elapsed picoseconds. Recording is controlled by a
 * stride in picoseconds of *simulation time* (not wall time), so replay is
 * deterministic regardless of frame rate.
 *
 * Export formats:
 *   - XYZ : standard multi-frame (comment line holds time in ps).
 *   - PDB : MODEL/ENDMDL blocks with one ATOM record per Cα bead using the
 *           original residue names/numbers/chain IDs (re-importable into
 *           PyMOL / VMD / MDAnalysis).
 *
 * Provenance (A09): every exported file starts with
 *   `REMARK simulation_coarse vX, T, gamma, seed, date` from src/version.js
 *   so trajectories are traceable to code version, temperature, friction,
 *   random seed and build date. Call `buildFile(fmt, beads, {T, gamma, seed})`
 *   to embed run-time values; defaults are 300 K, gamma 2.0, seed 0.
 */

import { VERSION, BUILD_DATE } from "./version.js?v=10";

export class Recorder {
  // G69 — Memory leak guard: maxFrames cap 500 documented; auto-stop when reached to prevent unbounded heap growth
  constructor() {
    this.frames = [];        // Float32Array(3n) clones
    this.times = [];         // ps at each frame
    this.recording = false;
    this.stridePs = 2.0;
    this.maxFrames = 500;    // 0 = unlimited — G69 cap documented (500 frames ~ few MB; prevents leak)
    this._nextAt = 0;        // simulation time of next capture
  }

  /** Begin recording from `nowPs` (sim time). */
  start(nowPs, stridePs, maxFrames) {
    this.frames.length = 0;
    this.times.length = 0;
    this.stridePs = Math.max(0.001, stridePs);
    this.maxFrames = Math.max(0, maxFrames);
    this._nextAt = nowPs;    // capture immediately
    this.recording = true;
  }

  stop() { this.recording = false; }
  clear() {
    this.frames.length = 0; this.times.length = 0; this.recording = false;
  }

  get count() { return this.frames.length; }
  /** Simulated nanoseconds spanned by the recording. */
  get spanNs() {
    return this.times.length > 1 ? (this.times[this.times.length - 1] - this.times[0]) / 1000 : 0;
  }

  // H78 — Trajectory scrubbing: random-access frame retrieval for <input type=range scrub>
  // UI placeholder: <input type="range" id="scrub" min="0" max="count-1" value="0">
  // oninput="const f = recorder.getFrame(Number(scrub.value)); viewer.render(f.pos)"
  /**
   * Retrieve a single frame by index for scrubbing / seeking.
   * @param {number} i frame index (0 .. count-1)
   * @returns {{pos: Float32Array, time: number}|null}
   */
  getFrame(i) {
    if (!Number.isInteger(i) || i < 0 || i >= this.frames.length) return null;
    return { pos: this.frames[i], time: this.times[i] };
  }

  /** Alias for scrub UI: total frame count */
  get numFrames() { return this.frames.length; }

  /**
   * Maybe capture a frame.
   * @param {Float64Array} pos   live positions (copied defensively)
   * @param {number} timePs      current simulation time
   */
  /**
   * Maybe capture a frame — G69 maxFrames cap guard.
   * @param {Float64Array} pos   live positions (copied defensively)
   * @param {number} timePs      current simulation time
   */
  maybeCapture(pos, timePs) {
    if (!this.recording) return false;
    if (timePs + 1e-9 < this._nextAt) return false;
    if (this.maxFrames > 0 && this.frames.length >= this.maxFrames) {
      this.recording = false;      // G69 auto-stop: maxFrames cap reached — prevents memory leak
      return false;
    }
    this.frames.push(Float32Array.from(pos));
    this.times.push(timePs);
    this._nextAt = timePs + this.stridePs;
    return true;
  }

  /* ----------------------------------------------------------------- */
  /*  Export                                                           */
  /* ----------------------------------------------------------------- */

  /**
   * @param {"xyz"|"pdb"} fmt
   * @param {Array} beads   bead metadata from pdb.selectSystem()
   * @param {object} [opts] provenance: {T, gamma, seed, date}
   * @returns {string} file contents
   */
  buildFile(fmt, beads, opts = {}) {
    if (this.frames.length === 0) throw new Error("No frames recorded yet.");
    return fmt === "pdb" ? this._toPdb(beads, opts) : this._toXyz(beads, opts);
  }

  /** Provenance header for exported files (A09). */
  _provenanceLine(opts = {}) {
    const T = opts.T ?? opts.temp ?? 300;
    const gamma = opts.gamma ?? 2;
    const seed = opts.seed ?? opts.rngSeed ?? 0;
    const date = opts.date ?? BUILD_DATE;
    return `REMARK simulation_coarse v${VERSION} T=${T}K gamma=${gamma} seed=${seed} date=${date}`;
  }

  _toXyz(beads, opts = {}) {
    const n = beads.length;
    const out = [this._provenanceLine(opts)];
    this.frames.forEach((fr, k) => {
      out.push(String(n));
      out.push(`frame ${k}  time = ${this.times[k].toFixed(3)} ps  (CG Cα model)`);
      for (let i = 0; i < n; i++) {
        out.push(
          `C  ${fr[3 * i].toFixed(3)}  ${fr[3 * i + 1].toFixed(3)}  ${fr[3 * i + 2].toFixed(3)}`
        );
      }
    });
    return out.join("\n") + "\n";
  }

  _toPdb(beads, opts = {}) {
    const n = beads.length;
    const out = [this._provenanceLine(opts), "REMARK  CG Cα Langevin trajectory (BAOAB integrator)"];
    this.frames.forEach((fr, k) => {
      out.push(`MODEL     ${String(k + 1).padStart(4)}`);
      out.push(`REMARK  time = ${this.times[k].toFixed(3)} ps`);
      for (let i = 0; i < n; i++) {
        const b = b0(beads[i]);
        out.push(
          "ATOM  " +
          String(i + 1).padStart(5) + " " +
          " CA ".padEnd(4) + " " +
          b.resName.padStart(3) + " " +
          (b.chain === "_" ? " " : b.chain) +
          String(b.resSeq).padStart(4) + "    " +
          fr[3 * i].toFixed(3).padStart(8) +
          fr[3 * i + 1].toFixed(3).padStart(8) +
          fr[3 * i + 2].toFixed(3).padStart(8) +
          "  1.00" + "  0.00".padStart(6) + "           C"
        );
      }
      out.push("ENDMDL");
    });
    out.push("END");
    return out.join("\n") + "\n";
  }
}

function b0(b) { // fallback metadata if fields missing
  return { resName: b.resName || "GLY", chain: b.chain || "A", resSeq: b.resSeq ?? 1 };
}

/** Trigger a browser download of `text` as `filename`. */
export function downloadText(text, filename) {
  const blob = new Blob([text], { type: "chemical/x-pdb;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
}
