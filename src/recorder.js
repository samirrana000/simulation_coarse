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
 */

export class Recorder {
  constructor() {
    this.frames = [];        // Float32Array(3n) clones
    this.times = [];         // ps at each frame
    this.recording = false;
    this.stridePs = 2.0;
    this.maxFrames = 500;    // 0 = unlimited
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

  /**
   * Maybe capture a frame.
   * @param {Float64Array} pos   live positions (copied defensively)
   * @param {number} timePs      current simulation time
   */
  maybeCapture(pos, timePs) {
    if (!this.recording) return false;
    if (timePs + 1e-9 < this._nextAt) return false;
    if (this.maxFrames > 0 && this.frames.length >= this.maxFrames) {
      this.recording = false;      // auto-stop: reached requested length
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
   * @returns {string} file contents
   */
  buildFile(fmt, beads) {
    if (this.frames.length === 0) throw new Error("No frames recorded yet.");
    return fmt === "pdb" ? this._toPdb(beads) : this._toXyz(beads);
  }

  _toXyz(beads) {
    const n = beads.length;
    const out = [];
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

  _toPdb(beads) {
    const n = beads.length;
    const out = ["REMARK  CG Cα Langevin trajectory (BAOAB integrator)"];
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
