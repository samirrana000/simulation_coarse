/**
 * viewer.js — Dependency-free real-time 3D Canvas viewer with multi-scale biophysical rendering.
 *
 * Features:
 *  - DSSP-inspired secondary structure cartoon ribbon rendering (Helices, Strands, Coils)
 *  - 2.5D depth-shaded CPK atom spheres and covalent wireframe bonds
 *  - Dynamic directional hydrogen bond (H-bond) visualization in real-time
 *  - 3D Chemical Network macrostate spatial overlays (Bulk, Encounter, Intermediate, Bound)
 *  - Dynamic Pocket Cavity Isosurface Contours & Vector Force Overlays
 *  - Multi-touch and mouse trackball controls with perspective picking
 */

import { VERSION } from "./version.js?v=10"; // A01 provenance — viewer knows build version (unused but validates deterministic import)

const CHAIN_PALETTE = [
  [86, 156, 214], // blue
  [197, 134, 192], // purple
  [106, 203, 166], // teal
  [220, 150, 86], // orange
  [181, 206, 168], // sage
  [240, 113, 120], // coral
  [255, 214, 102], // gold
  [156, 220, 254], // cyan
];

const ELEMENT_COLOR = {
  C: [180, 180, 180], N: [90, 130, 235], O: [235, 70, 70], S: [200, 180, 60],
  P: [200, 120, 40], F: [120, 210, 120], CL: [60, 200, 120], BR: [160, 80, 40],
  I: [140, 60, 160],
  ZN: [125, 128, 176], FE: [200, 120, 60], MG: [138, 255, 138], CA: [61, 89, 171],
  CU: [200, 130, 80], MN: [156, 122, 199], NI: [80, 160, 180], CO: [240, 140, 160],
  NA: [102, 102, 255], K: [143, 64, 212],
};
const ELEMENT_COLOR_DEFAULT = [230, 160, 200];

// Ligand palette — saturated, colorblind-distinct from protein CPK greys/blues.
// Class (protein vs ligand) is encoded BEFORE element (wiki pattern P1):
// protein keeps classic muted CPK; ligand gets vivid hues so a benzene C
// can never read as a protein C.
const LIGAND_COLOR = {
  C: [255, 121, 98],   N: [232, 121, 249], O: [38, 208, 206], S: [253, 224, 71],
  P: [168, 249, 122], F: [94, 246, 146], CL: [126, 252, 220], BR: [255, 160, 111],
  I: [249, 121, 205],
  ZN: [80, 250, 216], FE: [255, 140, 105], MG: [151, 255, 185], CA: [120, 190, 255],
  CU: [255, 170, 130], MN: [214, 143, 255], NI: [130, 230, 240], CO: [255, 160, 180],
  NA: [140, 180, 255], K: [200, 130, 255],
};
const LIGAND_COLOR_DEFAULT = [255, 159, 128];

const STATE_COLORS = [
  { stroke: "rgba(56, 189, 248, 0.55)", fill: "rgba(56, 189, 248, 0.08)", name: "Bulk Solvated (S0)" },
  { stroke: "rgba(251, 191, 36, 0.75)", fill: "rgba(251, 191, 36, 0.12)", name: "Encounter Complex (S1)" },
  { stroke: "rgba(192, 132, 252, 0.8)", fill: "rgba(192, 132, 252, 0.14)", name: "Intermediate (S2)" },
  { stroke: "rgba(52, 211, 153, 0.9)", fill: "rgba(52, 211, 153, 0.18)", name: "Native Bound (S3)" },
];

export class Viewer {
  /** @param {HTMLCanvasElement} canvas */
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas ? canvas.getContext("2d") : null;
    this.rotX = -0.4;
    this.rotY = 0.6;
    this.zoom = 1.0;
    this.panX = 0;
    this.panY = 0;

    this.showContacts = true;
    this.drawSpheres = true;
    this.drawRibbon = true;
    this.showHBonds = true;
    this.showStates = true;
    this.motionGain = 1.0;
    this._amp = null;
    this._mean = null;

    this.n = 0;
    this.nProt = 0;
    this.ligandBonds = null;
    this.holoIdx = null;
    this.colors = null;
    this.segments = [];
    this.secStruct = null; // "H" (helix), "E" (sheet), "C" (coil)
    this.contactIdx = null;
    this.center = [0, 0, 0];
    this.radius = 30;

    this.activeState = 0;
    this.pocketCenter = null;
    this._order = null;
    this._dragging = false; // touch: drag state shared with touchstart/touchmove (see _bindMouse)
    this._lastEmptyDraw = 0; // three-state: throttle no-data redraw to ≤1 Hz
    this._emptyWarned = false;

    if (typeof window !== "undefined" && this.canvas) {
      this._bindMouse();
      this._resize();
      window.addEventListener("resize", () => this._resize());
      // Observe flex layout changes (details toggle, sidebar collapse, DPR switch)
      if (typeof ResizeObserver !== "undefined") {
        try { new ResizeObserver(() => this._resize()).observe(this.canvas); } catch (_) {}
        try { new ResizeObserver(() => this._resize()).observe(this.canvas.parentElement); } catch (_) {}
      }
      if (typeof window.matchMedia === "function") {
        try {
          window.matchMedia(`(resolution: ${window.devicePixelRatio}dppx)`)
            .addEventListener("change", () => this._resize());
        } catch (_) {}
      }
    }
  }

  /* ------------------------------ setup ------------------------------- */

  setSystem(sel, ff) {
    if (!sel || !ff) return;
    const beads = sel.beads || sel.atoms || [];
    const segments = sel.segments || [];
    this.n = ff.n || beads.length;
    this.nProt = ff.nProt ?? this.n;
    this.ligandStart = ff.ligandStart ?? this.nProt;
    this.heavy = !!sel.heavy;
    this.ligandBonds = ff.covalentBonds && ff.covalentBonds.length ? ff.covalentBonds
      : (ff.ligandBonds && ff.ligandBonds.length ? ff.ligandBonds : null);
    this.holoIdx = ff.holoSprings && ff.holoSprings.length
      ? Uint32Array.from(ff.holoSprings.filter((_, k) => k % 3 !== 2))
      : null; // heuristic, not DSSP — ribbon assignment is Ca-distance heuristic, see docs/VIEWER.md
    this.segments = segments;
    this._mean = null;

    // DSSP Secondary structure assignment — heuristic, not DSSP (Ca d3 distance, not Kabsch-Sander H-bonds)
    this.secStruct = new Array(this.nProt).fill("C");
    const r = ff.ref || new Float64Array(this.n * 3);
    this.ref = r;
    for (const [s, e] of this.segments) {
      for (let i = s + 1; i < e - 2 && i < this.nProt - 2; i++) {
        const dx1 = r[3*(i+1)] - r[3*i], dy1 = r[3*(i+1)+1] - r[3*i+1], dz1 = r[3*(i+1)+2] - r[3*i+2];
        const dx3 = r[3*(i+3)] - r[3*i], dy3 = r[3*(i+3)+1] - r[3*i+1], dz3 = r[3*(i+3)+2] - r[3*i+2];
        const d3 = Math.hypot(dx3, dy3, dz3);
        if (d3 < 5.8) {
          this.secStruct[i] = "H";
          this.secStruct[i+1] = "H";
          this.secStruct[i+2] = "H";
        } else if (d3 > 8.5) {
          this.secStruct[i] = "E";
          this.secStruct[i+1] = "E";
        }
      }
    }

    this.colors = new Array(this.n);
    if (this.heavy) {
      for (let i = 0; i < this.n; i++) {
        const el = (beads[i] && beads[i].element) ? beads[i].element.toUpperCase() : "C";
        const isLigand = i >= this.ligandStart;
        this.colors[i] = isLigand
          ? (LIGAND_COLOR[el] || LIGAND_COLOR_DEFAULT)
          : (ELEMENT_COLOR[el] || ELEMENT_COLOR_DEFAULT);
      }
    } else {
      const chainIdx = new Map();
      let ci = 0;
      for (let i = 0; i < Math.min(this.nProt, beads.length); i++) {
        const b = beads[i];
        const chain = (b && b.chain) || "A";
        if (!chainIdx.has(chain)) chainIdx.set(chain, ci++ % CHAIN_PALETTE.length);
        this.colors[i] = CHAIN_PALETTE[chainIdx.get(chain)];
      }

      if (ff.ligandAtoms && ff.ligandAtoms.length) {
        const ligStart = this.ligandStart ?? this.nProt;
        for (let i = ligStart; i < this.n; i++) {
          const la = ff.ligandAtoms[i - ligStart];
          const el = (la && la.element) ? la.element.toUpperCase() : "C";
          this.colors[i] = LIGAND_COLOR[el] || LIGAND_COLOR_DEFAULT;
        }
      }
    }

    for (let i = 0; i < this.n; i++) {
      if (!this.colors[i]) this.colors[i] = ELEMENT_COLOR_DEFAULT;
    }

    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < this.nProt; i++) {
      cx += r[3 * i] || 0;
      cy += r[3 * i + 1] || 0;
      cz += r[3 * i + 2] || 0;
    }
    const safeN = Math.max(1, this.nProt);
    cx /= safeN; cy /= safeN; cz /= safeN;
    this.center = [cx, cy, cz];

    let maxR2 = 1;
    // radius from protein only — distant ligand must not inflate view
    for (let i = 0; i < this.nProt; i++) {
      const dx = (r[3 * i] || 0) - cx;
      const dy = (r[3 * i + 1] || 0) - cy;
      const dz = (r[3 * i + 2] || 0) - cz;
      maxR2 = Math.max(maxR2, dx * dx + dy * dy + dz * dz);
    }
    this.radius = Math.max(12, Math.sqrt(maxR2) * 1.2);

    const s = ff.springs || [];
    const pairs = [];
    for (let k = 0; k < s.length; k += 3) {
      if (s[k + 2] < 9.5 && pairs.length / 2 < 800) pairs.push(s[k], s[k + 1]);
    }
    this.contactIdx = Uint32Array.from(pairs);

    // Compute pocket center — external-ligand COM only (hetero excluded via ligandStart)
    const ligStartPc = this.ligandStart ?? this.nProt;
    if (this.n > ligStartPc) {
      let px = 0, py = 0, pz = 0, count = 0;
      for (let i = ligStartPc; i < this.n; i++) {
        px += r[3 * i]; py += r[3 * i + 1]; pz += r[3 * i + 2]; count++;
      }
      if (count > 0) this.pocketCenter = [px / count, py / count, pz / count];
    } else {
      this.pocketCenter = [cx, cy, cz];
    }

    this._order = new Int32Array(this.n);
    this._resize();
  }

  clear() {
    this.n = 0;
    this._requestClear = true;
  }

  setMotionGain(g) {
    this.motionGain = Math.min(6, Math.max(1, g));
  }

  setActiveState(stateIdx) {
    this.activeState = Math.max(0, Math.min(3, stateIdx));
  }

  /* ------------------------------ input ------------------------------- */

  _bindMouse() {
    const c = this.canvas;
    if (!c) return;
    let lx = 0, ly = 0;
    c.addEventListener("mousedown", (e) => { this._dragging = true; lx = e.clientX; ly = e.clientY; });
    window.addEventListener("mouseup", () => (this._dragging = false));
    window.addEventListener("mousemove", (e) => {
      if (!this._dragging) return;
      const dx = e.clientX - lx, dy = e.clientY - ly;
      lx = e.clientX; ly = e.clientY;
      this.rotY += dx * 0.008;
      this.rotX += dy * 0.008;
    });
    c.addEventListener("wheel", (e) => {
      e.preventDefault();
      this.zoom *= e.deltaY > 0 ? 0.92 : 1.087;
      this.zoom = Math.min(8, Math.max(0.1, this.zoom));
    }, { passive: false });
    // touch support for mobile/trackball (H80): single-finger rotate, pinch zoom
    c.addEventListener("touchstart", (e) => {
      if (e.touches.length === 1) { this._dragging = true; lx = e.touches[0].clientX; ly = e.touches[0].clientY; }
      else if (e.touches.length === 2) {
        this._touchDist = Math.hypot(e.touches[0].clientX - e.touches[1].clientX, e.touches[0].clientY - e.touches[1].clientY);
      }
      e.preventDefault();
    }, { passive: false });
    c.addEventListener("touchmove", (e) => {
      if (e.touches.length === 1 && this._dragging) {
        const dx = e.touches[0].clientX - lx, dy = e.touches[0].clientY - ly;
        lx = e.touches[0].clientX; ly = e.touches[0].clientY;
        this.rotY += dx * 0.008;
        this.rotX += dy * 0.008;
      } else if (e.touches.length === 2 && this._touchDist) {
        const d = Math.hypot(e.touches[0].clientX - e.touches[1].clientX, e.touches[0].clientY - e.touches[1].clientY);
        const factor = d / this._touchDist;
        this.zoom = Math.min(8, Math.max(0.1, this.zoom * factor));
        this._touchDist = d;
      }
      e.preventDefault();
    }, { passive: false });
    c.addEventListener("touchend", (e) => {
      if (e.touches.length === 0) { this._dragging = false; this._touchDist = null; }
      else if (e.touches.length === 1) { lx = e.touches[0].clientX; ly = e.touches[0].clientY; }
    }, { passive: false });
  }

  _resize() {
    if (!this.canvas) return;
    const dpr = Math.min(2, (typeof window !== "undefined" && window.devicePixelRatio) || 1);
    // Use getBoundingClientRect for sub-pixel accuracy; handle true 0x0 flex collapse
    const rect = (typeof this.canvas.getBoundingClientRect === "function")
      ? this.canvas.getBoundingClientRect()
      : { width: this.canvas.clientWidth, height: this.canvas.clientHeight };
    let w = Math.round(rect.width);
    let h = Math.round(rect.height);
    // Fallback only when rect is unreliable (e.g. display:none); otherwise
    // a true 0 means the flex container collapsed — retry next frame.
    if (w === 0 || h === 0) {
      if (this.canvas.clientWidth === 0 && this.canvas.clientHeight === 0) {
        // Genuine layout collapse — schedule retry, keep last good backing
        if (typeof requestAnimationFrame !== "undefined") {
          requestAnimationFrame(() => this._resize());
        }
        if (this.canvas.width === 0 || this.canvas.height === 0) {
          w = 600; h = 450;
        } else {
          return;
        }
      } else {
        w = w || 600;
        h = h || 450;
      }
    }
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.dpr = dpr;
  }

  /* ------------------------------ picking ------------------------------ */

  unproject(clientX, clientY, depth = 0) {
    if (this.n === 0 || !this.canvas) return null;
    const rect = this.canvas.getBoundingClientRect();
    const W = this.canvas.width, H = this.canvas.height;
    const kx = W / Math.max(1, rect.width), ky = H / Math.max(1, rect.height);
    const pxdev = (clientX - rect.left) * kx;
    const pydev = (clientY - rect.top) * ky;

    const scale = (Math.min(W, H) / (2.6 * Math.max(1, this.radius))) * this.zoom;
    const fov = 800; // unified with render() fov=800 for picking invertibility
    const z2 = depth;
    const persp = fov / Math.max(1e-3, fov + z2);
    const x1 = (pxdev - W / 2 - this.panX) / (scale * persp);
    const y1 = -(pydev - H / 2 - this.panY) / (scale * persp);

    const sx = Math.sin(this.rotX), cxx = Math.cos(this.rotX);
    const sy = Math.sin(this.rotY), cyy = Math.cos(this.rotY);
    const y0 = cxx * y1 + sx * z2;
    const z1 = -sx * y1 + cxx * z2;
    const x0 = cyy * x1 - sy * z1;
    const z0 = sy * x1 + cyy * z1;

    const [cx, cy, cz] = this.center;
    let wx = x0 + cx, wy = y0 + cy, wz = z0 + cz;
    // H73 — picking must invert motionGain amplification (gain 5 within 0.5Å)
    // Render amplifies positions via ref + (pos-ref)*motionGain; unproject inverts by finding
    // nearest amplified bead and de-amplifying: true = ref + (amplified - ref)/gain
    if (this.motionGain > 1.0 && this.motionGain !== 1 && this.ref && this.ref.length >= this.n * 3 && this.n > 0) {
      const gain = this.motionGain;
      const inv = 1.0 / gain;
      const pos = this._lastPos;
      const usePos = pos && pos.length === this.ref.length;
      let bestIdx = 0;
      let bestD2 = Infinity;
      for (let i = 0; i < this.n; i++) {
        const rx = this.ref[3 * i], ry = this.ref[3 * i + 1], rz = this.ref[3 * i + 2];
        if (!Number.isFinite(rx) || !Number.isFinite(ry) || !Number.isFinite(rz)) continue;
        let ax = rx, ay = ry, az = rz;
        if (usePos) {
          const dx = pos[3 * i] - rx, dy = pos[3 * i + 1] - ry, dz = pos[3 * i + 2] - rz;
          if (Number.isFinite(dx) && Number.isFinite(dy) && Number.isFinite(dz)) {
            ax = rx + dx * gain;
            ay = ry + dy * gain;
            az = rz + dz * gain;
          }
        }
        const ddx = wx - ax, ddy = wy - ay, ddz = wz - az;
        const d2 = ddx * ddx + ddy * ddy + ddz * ddz;
        if (d2 < bestD2) { bestD2 = d2; bestIdx = i; }
      }
      const rx = this.ref[3 * bestIdx], ry = this.ref[3 * bestIdx + 1], rz = this.ref[3 * bestIdx + 2];
      if (Number.isFinite(rx) && Number.isFinite(ry) && Number.isFinite(rz)) {
        wx = rx + (wx - rx) * inv;
        wy = ry + (wy - ry) * inv;
        wz = rz + (wz - rz) * inv;
      }
    }
    return [wx, wy, wz];
  }

  screenToWorld(clientX, clientY, opts = {}) {
    const depth = opts.depth ?? 0;
    return this.unproject(clientX, clientY, depth);
  }

  /* ------------------------------ render ------------------------------ */

  render(pos) {
    if (!this.canvas || !this.ctx) return;
    if (pos) this._lastPos = pos;
    const ctx = this.ctx;

    if (!ctx) return;
    const W = this.canvas.width;
    const H = this.canvas.height;

    ctx.fillStyle = "#0f172a";
    ctx.fillRect(0, 0, W, H);

    // Three-state (P2): no-data renders an actionable empty state at ≤1 Hz
    // (not every frame) so an empty canvas never reads as dead or spammy.
    if (!pos || this.n === 0) {
      const now = (typeof performance !== "undefined" && performance.now()) || Date.now();
      if (now - this._lastEmptyDraw < 1000 && this._emptyDrawn) return;
      this._lastEmptyDraw = now;
      this._emptyDrawn = true;
      if (!this._emptyWarned) {
        this._emptyWarned = true;
        if (typeof console !== "undefined") console.info("[viewer] No system — load a PDB and click Build System");
      }
      ctx.fillStyle = "#94a3b8";
      ctx.font = "12px sans-serif";
      ctx.textAlign = "left";
      ctx.fillText("No system — 1-click 4W52 sample (Structure) or load a PDB, then Build", 12, 24);
      ctx.fillStyle = "#475569";
      ctx.font = "10.5px sans-serif";
      ctx.fillText("Ligand atoms render vivid + white halo ring (see legend)", 12, 40);
      return;
    }
    this._emptyDrawn = false;

    if (!this._px || this._px.length !== this.n) {
      this._px = new Float64Array(this.n);
      this._py = new Float64Array(this.n);
      this._pz = new Float64Array(this.n);
    }
    const px = this._px, py = this._py, pz = this._pz;
    const [cx, cy, cz] = this.center;

    const scale = (Math.min(W, H) / (2.6 * this.radius)) * this.zoom;
    const cyy = Math.cos(this.rotY), sy = Math.sin(this.rotY);
    const cxx = Math.cos(this.rotX), sx = Math.sin(this.rotX);
    const fov = 800; // unified with unproject() fov=800 for picking invertibility

    const useGain = this.motionGain > 1.0 && this.ref && this.ref.length === pos.length;
    for (let i = 0; i < this.n; i++) {
      // handle Infinity (not just ||0): fall back to ref if pos is non-finite
      let rx = Number.isFinite(pos[3 * i]) ? pos[3 * i] : this.ref[3 * i];
      let ry = Number.isFinite(pos[3 * i + 1]) ? pos[3 * i + 1] : this.ref[3 * i + 1];
      let rz = Number.isFinite(pos[3 * i + 2]) ? pos[3 * i + 2] : this.ref[3 * i + 2];
      if (useGain) {
        rx = this.ref[3 * i] + (rx - this.ref[3 * i]) * this.motionGain;
        ry = this.ref[3 * i + 1] + (ry - this.ref[3 * i + 1]) * this.motionGain;
        rz = this.ref[3 * i + 2] + (rz - this.ref[3 * i + 2]) * this.motionGain;
      }
      const x0 = rx - cx;
      const y0 = ry - cy;
      const z0 = rz - cz;
      const x1 = cyy * x0 + sy * z0;
      const z1 = -sy * x0 + cyy * z0;
      const y1 = cxx * y0 - sx * z1;
      const z2 = sx * y0 + cxx * z1;
      const persp = fov / Math.max(1e-3, fov + z2);
      px[i] = W / 2 + x1 * scale * persp + this.panX;
      py[i] = H / 2 - y1 * scale * persp + this.panY;
      pz[i] = z2;
    }

    // 1. 3D Chemical Network Overlays
    if (this.showStates && this.pocketCenter) {
      const pk = this.pocketCenter;
      const pk_x0 = pk[0] - cx, pk_y0 = pk[1] - cy, pk_z0 = pk[2] - cz;
      const pk_x1 = cyy * pk_x0 + sy * pk_z0;
      const pk_z1 = -sy * pk_x0 + cyy * pk_z0;
      const pk_y1 = cxx * pk_y0 - sx * pk_z1;
      const pk_z2 = sx * pk_y0 + cxx * pk_z1;
      const pk_persp = fov / Math.max(1e-3, fov + pk_z2);
      const pk_px = W / 2 + pk_x1 * scale * pk_persp + this.panX;
      const pk_py = H / 2 - pk_y1 * scale * pk_persp + this.panY;

      const stateInfo = STATE_COLORS[this.activeState] || STATE_COLORS[0];
      const radMap = [40, 26, 16, 9];
      const rPix = radMap[this.activeState] * scale * pk_persp;

      ctx.save();
      ctx.beginPath();
      ctx.arc(pk_px, pk_py, Math.max(6, rPix), 0, 2 * Math.PI);
      ctx.fillStyle = stateInfo.fill;
      ctx.fill();
      ctx.strokeStyle = stateInfo.stroke;
      ctx.lineWidth = 1.5 * (this.dpr || 1);
      ctx.setLineDash([4, 4]);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.restore();
    }

    // 2. ENM contact sticks
    if (this.showContacts && this.contactIdx) {
      ctx.lineWidth = (this.dpr || 1) * 0.5;
      ctx.strokeStyle = "rgba(100, 120, 150, 0.14)";
      ctx.beginPath();
      for (let k = 0; k < this.contactIdx.length; k += 2) {
        const i = this.contactIdx[k], j = this.contactIdx[k + 1];
        if (i < this.n && j < this.n) {
          ctx.moveTo(px[i], py[i]);
          ctx.lineTo(px[j], py[j]);
        }
      }
      ctx.stroke();
    }

    // 3. Dynamic Protein-Ligand Hydrogen Bonds (j-set = true ligand only, hetero excluded)
    const ligStartHb = this.ligandStart ?? this.nProt;
    if (this.showHBonds && this.n > ligStartHb) {
      ctx.lineWidth = (this.dpr || 1) * 1.5;
      ctx.strokeStyle = "#38bdf8";
      ctx.setLineDash([4, 3]);
      ctx.beginPath();
      for (let i = 0; i < this.nProt; i += 2) {
        const xi = pos[3 * i], yi = pos[3 * i + 1], zi = pos[3 * i + 2];
        for (let j = ligStartHb; j < this.n; j++) {
          const dx = pos[3 * j] - xi, dy = pos[3 * j + 1] - yi, dz = pos[3 * j + 2] - zi;
          const d2 = dx * dx + dy * dy + dz * dz;
          if (d2 > 6.0 && d2 < 11.5) {
            ctx.moveTo(px[i], py[i]);
            ctx.lineTo(px[j], py[j]);
          }
        }
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // 4. Covalent bonds
    if (this.ligandBonds) {
      ctx.lineWidth = (this.dpr || 1) * 2.0;
      ctx.strokeStyle = "rgba(226, 232, 240, 0.85)";
      ctx.beginPath();
      for (let k = 0; k < this.ligandBonds.length; k += 3) {
        const i = this.ligandBonds[k], j = this.ligandBonds[k + 1];
        if (i < this.n && j < this.n) {
          ctx.moveTo(px[i], py[i]);
          ctx.lineTo(px[j], py[j]);
        }
      }
      ctx.stroke();
    }
    // painter sort fallback: Canvas2D uses order.sort by pz (painter sort); future WebGL uses DEPTH_TEST depth buffer — GL depth buffer future
    // 5. DSSP Secondary Structure Smooth Ribbon Backbone (heuristic, not DSSP — see docs/VIEWER.md)
    if (this.drawRibbon && this.segments && this.segments.length) {
      for (const [s, e] of this.segments) {
        if (e - s < 2) continue;
        for (let i = s; i < e - 1 && i < this.nProt - 1; i++) {
          const ss = this.secStruct ? this.secStruct[i] : "C";
          const p0x = px[Math.max(s, i - 1)], p0y = py[Math.max(s, i - 1)];
          const p1x = px[i], p1y = py[i];
          const p2x = px[i + 1], p2y = py[i + 1];
          const p3x = px[Math.min(e - 1, i + 2)], p3y = py[Math.min(e - 1, i + 2)];

          ctx.beginPath();
          ctx.moveTo(p1x, p1y);
          // 4-step Catmull-Rom spline interpolation
          for (let st = 1; st <= 4; st++) {
            const t = st / 4.0;
            const t2 = t * t, t3 = t2 * t;
            const qx = 0.5 * ((2 * p1x) + (-p0x + p2x) * t + (2 * p0x - 5 * p1x + 4 * p2x - p3x) * t2 + (-p0x + 3 * p1x - 3 * p2x + p3x) * t3);
            const qy = 0.5 * ((2 * p1y) + (-p0y + p2y) * t + (2 * p0y - 5 * p1y + 4 * p2y - p3y) * t2 + (-p0y + 3 * p1y - 3 * p2y + p3y) * t3);
            ctx.lineTo(qx, qy);
          }

          if (ss === "H") {
            ctx.lineWidth = (this.dpr || 1) * 3.8; // Thick Helical Ribbon
            ctx.strokeStyle = "rgba(251, 191, 36, 0.85)"; // Amber gold
          } else if (ss === "E") {
            ctx.lineWidth = (this.dpr || 1) * 3.0; // Flat Beta-Strand Ribbon
            ctx.strokeStyle = "rgba(56, 189, 248, 0.85)"; // Cyan blue
          } else {
            ctx.lineWidth = (this.dpr || 1) * 1.5; // Slender Loop
            ctx.strokeStyle = "rgba(148, 163, 184, 0.65)"; // Slate grey
          }
          ctx.stroke();
        }
      }
    }

    // 6. Beads / Spheres (Painter's sort fallback: far -> near — painter sort; GL depth buffer future via DEPTH_TEST)
    if (!this._order || this._order.length !== this.n) {
      this._order = new Int32Array(this.n);
    }
    const order = this._order;
    for (let i = 0; i < this.n; i++) order[i] = i;
    order.sort((a, b) => pz[b] - pz[a]);

    const baseR = this.drawSpheres ? 4.2 * (this.dpr || 1) * Math.sqrt(this.zoom) : 1.8 * (this.dpr || 1);
    for (const i of order) {
      const depth = Math.max(0.25, Math.min(1, 1 - pz[i] / (this.radius * 2.2)));
      const isLig = i >= this.ligandStart;
      const rScale = (this.heavy || isLig) ? (isLig ? 0.9 : 0.65) : 1.0;
      const persp = fov / Math.max(1e-3, fov + pz[i]);
      const r = Math.max(0.6, rScale * baseR * persp);

      const col = (this.colors && this.colors[i]) || ELEMENT_COLOR_DEFAULT;
      const [cr, cg, cb] = col;
      const rr = Math.max(0, Math.min(255, Math.round(cr * depth + (isLig ? 35 : 15))));
      const gg = Math.max(0, Math.min(255, Math.round(cg * depth + (isLig ? 35 : 15))));
      const bb = Math.max(0, Math.min(255, Math.round(cb * depth + (isLig ? 35 : 15))));

      ctx.beginPath();
      ctx.arc(px[i], py[i], r, 0, 2 * Math.PI);
      ctx.fillStyle = `rgb(${rr},${gg},${bb})`;
      ctx.fill();

      // Ligand halo ring — second identity channel (wiki P1): a white ring
      // marks covalently-bonded ligand atoms even where hues look close.
      if (isLig && this.drawSpheres) {
        ctx.beginPath();
        ctx.arc(px[i], py[i], Math.max(1.2, r + 1.6 * (this.dpr || 1)), 0, 2 * Math.PI);
        ctx.strokeStyle = "rgba(255, 255, 255, 0.95)";
        ctx.lineWidth = (this.dpr || 1) * 0.9;
        ctx.stroke();
      }

      if (this.drawSpheres && r > 2.2) {
        ctx.beginPath();
        ctx.arc(px[i] - r * 0.35, py[i] - r * 0.35, Math.max(0.3, r * 0.38), 0, 2 * Math.PI);
        ctx.fillStyle = "rgba(255,255,255,0.28)";
        ctx.fill();
      }
    }
  }
}
