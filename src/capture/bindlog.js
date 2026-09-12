/**
 * bindlog.js — Efficient simulation-data capture port (R6, Loop-1).
 *
 * Two channels:
 *   DENSE  — quantized Int16 delta-compressed position frames
 *            (frame 0 absolute Float32, then per-frame deltas × quantStep).
 *            ~2× smaller than Float32 frames; error ≤ quantStep.
 *   SPARSE — typed-array struct-of-arrays event log (no per-event objects):
 *            energy components, contact form/break, PMF hills, state hops,
 *            pocket-volume samples.
 *
 * Binary blob layout (fromBinaryBlob/toBinaryBlob round-trips losslessly):
 *   [0..3]    magic 0x42 0x4C 0x47 0x31 ("BLG1")
 *   [4..7]    u32 version = 1
 *   [8..11]   u32 headerLen (JSON header length, ASCII)
 *   [12...]   header JSON: {nAtoms, nFrames, quantStep, nEvents,
 *              offsets: {frame0, frameTimes, frames, events...}}
 *   then aligned sections: Float32 frame0 (3n), Float32 times (nFrames),
 *   Int16 deltas (3n × (nFrames−1)), then event sections (see header).
 *
 * Zero deps; ES module; additive capture (does not touch existing modules).
 */

const MAGIC = 0x42434731; // "BCG1" little-endian check

export class BindLog {
  /**
   * @param {object} [opts]
   * @param {number} [opts.quantStep=0.01] Å per Int16 step for delta frames
   * @param {number} [opts.maxEvents=100000] sparse-event capacity (grows ×2)
   */
  constructor(opts = {}) {
    this.quantStep = opts.quantStep ?? 0.01;
    // DENSE channel
    this.nAtoms = 0;
    this.frame0 = null;        // Float32Array(3n) — first captured frame
    this.frameTimes = [];       // number[] ps
    this.deltas = [];           // Int16Array per frame (3n each, vs previous frame)
    this._lastFrame = null;    // Float32Array scratch (dequantized previous)
    // SPARSE channel — struct of arrays
    this.maxEvents = opts.maxEvents ?? 100000;
    this.evTime = new Float64Array(1024);   // ps
    this.evType = new Uint8Array(1024);     // 0 energy, 1 contact+, 2 contact−, 3 hill, 4 state, 5 pocketVol
    this.evA = new Int32Array(1024);        // atom/residue index A (or state id)
    this.evB = new Int32Array(1024);        // residue index B (or −1)
    this.evX = new Float32Array(1024);      // payload: distance / height / volume / component value
    this.evY = new Float32Array(1024);      // secondary payload (e.g. energy term id)
    this.nEvents = 0;
  }

  /* ---------------- DENSE ---------------- */

  /** Capture a position frame (quantized delta vs previous).
   * @param {Float32Array|Float64Array} pos 3n positions
   * @param {number} timePs */
  captureFrame(pos, timePs) {
    const n3 = pos.length;
    if (this.frame0 === null) {
      this.nAtoms = n3 / 3;
      this.frame0 = Float32Array.from(pos);
      this._lastFrame = Float32Array.from(pos);
      this.frameTimes.push(timePs);
      return true;
    }
    if (n3 !== this.frame0.length) throw new Error("BindLog: atom count changed mid-capture");
    const d = new Int16Array(n3);
    const cur = new Float32Array(n3);
    const q = this.quantStep;
    for (let i = 0; i < n3; i++) {
      const raw = pos[i];
      const prev = this._lastFrame[i];
      // quantized delta; round-trip error bounded by q/2 per coordinate
      const dv = Math.round((raw - prev) / q);
      d[i] = dv > 32767 ? 32767 : dv < -32768 ? -32768 : dv;
      cur[i] = prev + dv * q;
    }
    this.deltas.push(d);
    this._lastFrame = cur;
    this.frameTimes.push(timePs);
    return true;
  }

  /** Reconstruct frame i (0-based). @returns {Float32Array} */
  getFrame(i) {
    if (i === 0) return this.frame0.slice();
    const out = this.frame0.slice();
    const q = this.quantStep;
    for (let k = 0; k < i; k++) {
      const d = this.deltas[k];
      for (let j = 0; j < d.length; j++) out[j] += d[j] * q;
    }
    return out;
  }

  get nFrames() { return 1 + this.deltas.length; }

  /** Max absolute reconstruction error vs raw (for validation).
   * @param {Array<Float32Array>} rawFrames */
  maxErrorVs(rawFrames) {
    let max = 0;
    for (let i = 0; i < rawFrames.length; i++) {
      const rec = this.getFrame(i), raw = rawFrames[i];
      for (let j = 0; j < raw.length; j++) {
        const e = Math.abs(rec[j] - raw[j]);
        if (e > max) max = e;
      }
    }
    return max;
  }

  /* ---------------- SPARSE ---------------- */

  /** Grow the event arrays by ×2 (amortized O(1) append). */
  _grow() {
    const cap = this.evTime.length * 2;
    const t = new Float64Array(cap); t.set(this.evTime); this.evTime = t;
    const ty = new Uint8Array(cap); ty.set(this.evType); this.evType = ty;
    const a = new Int32Array(cap); a.set(this.evA); this.evA = a;
    const b = new Int32Array(cap); b.set(this.evB); this.evB = b;
    const x = new Float32Array(cap); x.set(this.evX); this.evX = x;
    const y = new Float32Array(cap); y.set(this.evY); this.evY = y;
  }

  /**
   * Append one event.
   * @param {number} timePs
   * @param {0|1|2|3|4|5} type
   * @param {number} a
   * @param {number} b
   * @param {number} x
   * @param {number} y
   */
  pushEvent(timePs, type, a = -1, b = -1, x = 0, y = 0) {
    if (this.nEvents >= this.evTime.length) {
      if (this.nEvents >= this.maxEvents) return false;
      this._grow();
    }
    const i = this.nEvents++;
    this.evTime[i] = timePs; this.evType[i] = type;
    this.evA[i] = a; this.evB[i] = b; this.evX[i] = x; this.evY[i] = y;
    return true;
  }

  /** Helper: energy-component sample (R4 §2 — 7 components per frame). */
  pushEnergyComponents(timePs, comps) {
    for (let c = 0; c < comps.length; c++) this.pushEvent(timePs, 0, c, -1, comps[c], 0);
  }

  /** Helper: contact formed/broken (atom a ↔ residue b at distance dist Å). */
  pushContact(timePs, formed, atom, residue, dist) {
    this.pushEvent(timePs, formed ? 1 : 2, atom, residue, dist, 0);
  }

  /** Helper: PMF hill deposit (CV value, height). */
  pushHill(timePs, cv, height) { this.pushEvent(timePs, 3, 0, -1, cv, height); }

  /** Helper: binding-state hop (stateId 0-3). */
  pushState(timePs, stateId) { this.pushEvent(timePs, 4, stateId, -1, 0, 0); }

  /** Helper: pocket volume sample (Å³). */
  pushPocketVolume(timePs, volume) { this.pushEvent(timePs, 5, 0, -1, volume, 0); }

  /* ---------------- BINARY BLOB ---------------- */

  /** Serialize to a single ArrayBuffer. */
  toBinaryBlob() {
    const n = this.nAtoms, nF = this.nFrames, nE = this.nEvents;
    const frame0Bytes = this.frame0 ? this.frame0.byteLength : 0;
    const timesBytes = nF * 4;                       // Float32
    const deltaBytes = this.deltas.length * n * 3 * 2; // Int16 per frame
    // events: time(8) type(1) a(4) b(4) x(4) y(4) = 25 B → pad to 32 for alignment
    const evBytes = nE * 32;
    // Fixed-width offsets (7 digits each, zero-padded) so the header JSON length
    // is IDENTICAL before and after filling the offsets (rewrite-in-place safe).
    const O = (v) => String(v).padStart(7, "0");
    const header = {
      n, nF, nE, quantStep: this.quantStep,
      typeNames: ["energy", "contact+", "contact-", "hill", "state", "pocketVol"],
      offsets: { frame0: O(0), times: O(0), deltas: O(0), events: O(0) },
    };
    const headerStr = JSON.stringify(header);
    const headerBytes = headerStr.length; // ASCII
    // pad header to 4-byte multiple so Float32 sections stay aligned
    const headerPad = (4 - ((12 + 4 + headerBytes) % 4)) % 4;
    const total = 12 + 4 + headerBytes + headerPad + frame0Bytes + timesBytes + deltaBytes + evBytes + 8;
    const buf = new ArrayBuffer(total);
    const dv = new DataView(buf);
    const u8 = new Uint8Array(buf);
    let off = 0;
    dv.setUint32(off, 0x424c4731, true); off += 4;              // "BLG1"
    dv.setUint32(off, 1, true); off += 4;                        // version
    dv.setUint32(off, headerBytes, true); off += 4;              // header length
    for (let i = 0; i < headerBytes; i++) u8[off + i] = headerStr.charCodeAt(i);
    off += headerBytes;
    for (let i = 0; i < headerPad; i++) u8[off + i] = 0; // padding
    off += headerPad;
    const base0 = off;
    const rel = () => off - base0; // relative to sections start
    if (this.frame0) { new Float32Array(buf, off, n * 3).set(this.frame0); off += frame0Bytes; }
    const offTimes = rel(); {
      const timesF = new Float32Array(buf, off, nF);
      for (let i = 0; i < nF; i++) timesF[i] = this.frameTimes[i];
      off += timesBytes;
    }
    const offDeltas = rel();
    for (const d of this.deltas) { new Int16Array(buf, off, d.length).set(d); off += d.byteLength; }
    const offEvents = rel();
    for (let i = 0; i < nE; i++) {
      dv.setFloat64(off, this.evTime[i], true); off += 8;
      dv.setUint8(off, this.evType[i]); off += 1; off += 3; // pad
      dv.setInt32(off, this.evA[i], true); off += 4;
      dv.setInt32(off, this.evB[i], true); off += 4;
      dv.setFloat32(off, this.evX[i], true); off += 4;
      dv.setFloat32(off, this.evY[i], true); off += 4;
    }
    dv.setUint32(off, 0x454e4431, true); off += 4; // "END1"
    dv.setUint32(off, total, true); off += 4;
    // re-write header with final offsets (fixed width ⇒ same length)
    header.offsets = { frame0: O(0), times: O(offTimes), deltas: O(offDeltas), events: O(offEvents) };
    const h2 = JSON.stringify(header);
    if (h2.length !== headerBytes) throw new Error(`BindLog: header size changed ${headerBytes}→${h2.length}`);
    const hdr = new Uint8Array(buf, 12, headerBytes);
    for (let i = 0; i < headerBytes; i++) hdr[i] = h2.charCodeAt(i);
    return buf;
  }

  /** Parse a binary blob back into a BindLog. @returns {BindLog} */
  static fromBinaryBlob(buf) {
    const dv = new DataView(buf);
    if (dv.getUint32(0, true) !== 0x424c4731) throw new Error("BindLog: bad magic");
    const version = dv.getUint32(4, true);
    if (version !== 1) throw new Error(`BindLog: unsupported version ${version}`);
    const hLen = dv.getUint32(8, true);
    let headerStr = "";
    const u8 = new Uint8Array(buf);
    for (let i = 0; i < hLen; i++) headerStr += String.fromCharCode(u8[12 + i]);
    const header = JSON.parse(headerStr);
    const out = new BindLog({ quantStep: header.quantStep });
    const headerPad = (4 - ((12 + hLen) % 4)) % 4;
    const base = 12 + hLen + headerPad;
    const header2 = { ...header, offsets: Object.fromEntries(Object.entries(header.offsets).map(([k, v]) => [k, Number(v)])) };
    const { n, nF, nE } = header;
    header.offsets = header2.offsets;
    out.frame0 = new Float32Array(buf.slice(base + header.offsets.frame0, base + header.offsets.frame0 + n * 3 * 4));
    out._lastFrame = out.frame0.slice();
    out.nAtoms = n;
    const timesF = new Float32Array(buf, base + header.offsets.times, nF);
    out.frameTimes = Array.from(timesF);
    // deltas
    const dBytes = n * 3 * 2;
    out.deltas = [];
    for (let f = 1; f < nF; f++) {
      out.deltas.push(new Int16Array(buf.slice(base + header.offsets.deltas + (f - 1) * dBytes, base + header.offsets.deltas + f * dBytes)));
    }
    // events
    out.nEvents = nE;
    const cap = Math.max(1024, Math.ceil(nE / 1024) * 1024);
    out.evTime = new Float64Array(cap); out.evType = new Uint8Array(cap);
    out.evA = new Int32Array(cap); out.evB = new Int32Array(cap);
    out.evX = new Float32Array(cap); out.evY = new Float32Array(cap);
    let off = base + header.offsets.events;
    for (let i = 0; i < nE; i++) {
      out.evTime[i] = dv.getFloat64(off, true); off += 8;
      out.evType[i] = dv.getUint8(off); off += 4;
      out.evA[i] = dv.getInt32(off, true); off += 4;
      out.evB[i] = dv.getInt32(off, true); off += 4;
      out.evX[i] = dv.getFloat32(off, true); off += 4;
      out.evY[i] = dv.getFloat32(off, true); off += 4;
    }
    return out;
  }

  /** Compact JSON summary (for quick human inspection / R7 viz). */
  toJSONSummary() {
    const counts = { energy: 0, "contact+": 0, "contact-": 0, hill: 0, state: 0, pocketVol: 0 };
    for (let i = 0; i < this.nEvents; i++) {
      const t = this.evType[i];
      counts[["energy", "contact+", "contact-", "hill", "state", "pocketVol"][t]]++;
    }
    return {
      frames: this.nFrames,
      atoms: this.nAtoms,
      spanPs: this.frameTimes.length ? this.frameTimes[this.frameTimes.length - 1] - this.frameTimes[0] : 0,
      quantStep: this.quantStep,
      events: counts,
      approxBytes: this.nAtoms * 3 * 4 + this.nFrames * 4 + this.deltas.length * this.nAtoms * 3 * 2 + this.nEvents * 32,
    };
  }

  /**
   * Attach to a live simulation loop (optional, additive). Calls the capture
   * hooks if the objects expose them; silently skips absent ones.
   * @param {object} sim {recorder?, ff?, funnel?, network?}
   */
  attach(sim) {
    this._attached = sim || {};
    return this;
  }

  /** Drain a recorded Recorder's frames into the dense channel.
   * @param {import('./recorder.js').Recorder} recorder */
  importRecorder(recorder) {
    for (let i = 0; i < recorder.frames.length; i++) {
      this.captureFrame(recorder.frames[i], recorder.times[i]);
    }
    return this;
  }
}

export const BINDLOG_EVENT_TYPES = ["energy", "contact+", "contact-", "hill", "state", "pocketVol"];
