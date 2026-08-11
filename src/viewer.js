/**
 * viewer.js — Dependency-free real-time 3D Canvas viewer.
 *
 * Renders Cα beads as depth-shaded spheres (or points) connected by peptide
 * "cartoon" lines, plus optional elastic-network contact sticks. Uses a simple
 * trackball rotation (quaternion-free YXZ Euler) + perspective projection —
 * sufficient for a few-thousand-bead CG model at 60 fps.
 *
 * Render cost is O(N + nContacts) per frame with precomputed native colors
 * (chain-based palette). Canvas 2D `arc()` with a radial gradient fake-shading
 * keeps it fast; bead sort every frame gives correct painter's-algorithm
 * depth ordering.
 *
 * Picking: unproject()/pickDepth()/screenToWorld() invert the geometric
 * projection (pan → perspective divide → rotateX/rotateY → center) so a
 * cursor position maps back to a world point (Å) for ligand placement.
 * The render-only motion magnification is NOT inverted — picking always
 * targets the true, un-amplified world frame.
 */

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
};
const ELEMENT_COLOR_DEFAULT = [230, 160, 200];

export class Viewer {
  /** @param {HTMLCanvasElement} canvas */
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext("2d");
    this.rotX = -0.4;                // trackball Euler angles (rad)
    this.rotY = 0.6;
    this.zoom = 1.0;                 // multiplicative
    this.panX = 0; this.panY = 0;

    this.showContacts = true;
    this.drawSpheres = true;

    this.motionGain = 1.3;         // render-only × amplification of the thermal
                                   // jitter around the running mean (HUD/recorder
                                   // keep the exact, un-amplified coordinates)
    this._amp = null;              // scratch: amplified coordinates (3n)
    this._mean = null;             // running mean of pos (slow EMA)

    this.n = 0;
    this.nProt = 0;
    this.ligandBonds = null;
    this.holoIdx = null;
    this.colors = null;              // per-bead [r,g,b]
    this.segments = [];              // backbone polylines
    this.contactIdx = null;          // Uint32Array [i0,j0,i1,j1,...]
    this.center = [0, 0, 0];         // model centroid (for rotation pivot)
    this.radius = 30;                // bounding radius (Å) → base scale

    this._order = null;              // painter sort
    this._dragging = false;
    this._bindMouse();
    this._resize();
    window.addEventListener("resize", () => this._resize());
  }

  /* ------------------------------ setup ------------------------------- */

  /** (Re)bind a CG system: beads metadata + topology. */
  setSystem(sel, ff) {
    const { beads, segments } = sel;
    this.n = ff.n;
    this.nProt = ff.nProt;
    this.ligandBonds = ff.ligandBonds && ff.ligandBonds.length ? ff.ligandBonds : null;
    this.holoIdx = ff.holoSprings && ff.holoSprings.length ? Uint32Array.from(ff.holoSprings.filter((_, k) => k % 3 !== 2)) : null;
    this.segments = segments;
    this._mean = null; // new structure ⇒ reset the motion-magnification mean

    // per-bead colors by chain (protein, indices 0..nProt-1)
    const chainIdx = new Map();
    let ci = 0;
    this.colors = new Array(this.n);
    beads.forEach((b, i) => {
      if (!chainIdx.has(b.chain)) chainIdx.set(b.chain, ci++ % CHAIN_PALETTE.length);
      this.colors[i] = CHAIN_PALETTE[chainIdx.get(b.chain)];
    });

    // ligand atoms colored by element (indices nProt..n-1)
    if (ff.ligandAtoms) {
      for (let i = this.nProt; i < this.n; i++) {
        const la = ff.ligandAtoms[i - this.nProt];
        this.colors[i] = ELEMENT_COLOR[la.element] || ELEMENT_COLOR_DEFAULT;
      }
    }

    // centroid + bounding radius from native reference
    const r = ff.ref;
    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < this.n; i++) { cx += r[3 * i]; cy += r[3 * i + 1]; cz += r[3 * i + 2]; }
    cx /= this.n; cy /= this.n; cz /= this.n;
    this.center = [cx, cy, cz];
    let maxR2 = 1;
    for (let i = 0; i < this.n; i++) {
      const dx = r[3 * i] - cx, dy = r[3 * i + 1] - cy, dz = r[3 * i + 2] - cz;
      maxR2 = Math.max(maxR2, dx * dx + dy * dy + dz * dz);
    }
    this.radius = Math.sqrt(maxR2) * 1.25;

    // strip contacts for drawing: every 3rd ENM spring to avoid clutter
    const s = ff.springs, pairs = [];
    for (let k = 0; k < s.length; k += 3) {
      if (s[k + 2] < 9.5 && pairs.length / 2 < 800) pairs.push(s[k], s[k + 1]);
    }
    this.contactIdx = Uint32Array.from(pairs);

    this._order = new Int32Array(this.n);
    this._resize();
  }

  clear() { this.n = 0; this._requestClear = true; }

  /** View-only amplification of thermal jitter; clamped to [1, 6]. */
  setMotionGain(g) {
    this.motionGain = Math.min(6, Math.max(1, g));
  }

  /* ------------------------------ input ------------------------------- */

  _bindMouse() {
    const c = this.canvas;
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
  }

  _resize() {
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    if (w === 0 || h === 0) return;
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.dpr = dpr;
  }

  /* ------------------------------ picking ------------------------------ */

  /**
   * Invert the geometric projection: screen (CSS px) → world (Å).
   * `depth` is the camera-space z2 of the plane to intersect
   * (default 0 = plane through the model center). The render-only motion
   * magnification is NOT inverted — picking targets the true world frame.
   * @returns {[number,number,number]|null} world point, or null if no system.
   */
  unproject(clientX, clientY, depth = 0) {
    if (this.n === 0) return null;
    const rect = this.canvas.getBoundingClientRect();
    const W = this.canvas.width, H = this.canvas.height;
    // CSS px → device px via the rect/buffer ratio (robust to CSS transforms;
    // equals this.dpr for the standard 1:1 layout).
    const kx = W / rect.width, ky = H / rect.height;
    const pxdev = (clientX - rect.left) * kx;
    const pydev = (clientY - rect.top) * ky;

    const scale = (Math.min(W, H) * 0.45 * this.zoom) / this.radius;
    const fov = this.radius * 4;
    const z2 = depth;
    const persp = fov / (fov + z2);
    const x1 = (pxdev - W / 2 - this.panX) / (scale * persp);
    const y1 = -(pydev - H / 2 - this.panY) / (scale * persp);

    // forward: rotateY(rotY) then rotateX(rotX); inverse: rotateX(-rotX) then rotateY(-rotY)
    const sx = Math.sin(this.rotX), cxx = Math.cos(this.rotX);
    const sy = Math.sin(this.rotY), cyy = Math.cos(this.rotY);
    const y0 = cxx * y1 + sx * z2;
    const z1 = -sx * y1 + cxx * z2;
    const x0 = cyy * x1 - sy * z1;
    const z0 = sy * x1 + cyy * z1;

    const [cx, cy, cz] = this.center;
    return [x0 + cx, y0 + cy, z0 + cz];
  }

  /**
   * Camera-space depth (z2) of the bead nearest the cursor within 12 device
   * px, preferring the nearest bead along z among candidates. Returns 0 when
   * nothing is under the cursor. `pos` defaults to the last rendered
   * (un-amplified) position buffer.
   */
  pickDepth(clientX, clientY, pos = this._lastPos) {
    if (this.n === 0 || !pos) return 0;
    const rect = this.canvas.getBoundingClientRect();
    const W = this.canvas.width, H = this.canvas.height;
    const kx = W / rect.width, ky = H / rect.height;
    const mx = (clientX - rect.left) * kx;
    const my = (clientY - rect.top) * ky;

    const scale = (Math.min(W, H) * 0.45 * this.zoom) / this.radius;
    const fov = this.radius * 4;
    const [cx, cy, cz] = this.center;
    const sx = Math.sin(this.rotX), cxx = Math.cos(this.rotX);
    const sy = Math.sin(this.rotY), cyy = Math.cos(this.rotY);

    const tol = 12, tol2 = tol * tol; // device px
    let best = null, bestD2 = Infinity;
    for (let i = 0; i < this.n; i++) {
      const x0 = pos[3 * i] - cx, y0 = pos[3 * i + 1] - cy, z0 = pos[3 * i + 2] - cz;
      const x1 = cyy * x0 + sy * z0;
      const z1 = -sy * x0 + cyy * z0;
      const y1 = cxx * y0 - sx * z1;
      const z2 = sx * y0 + cxx * z1;
      const persp = fov / (fov + z2);
      const px = W / 2 + x1 * scale * persp + this.panX;
      const py = H / 2 - y1 * scale * persp + this.panY;
      const dx = px - mx, dy = py - my;
      const d2 = dx * dx + dy * dy;
      if (d2 > tol2) continue;
      if (!best || z2 < best.z2 || (z2 === best.z2 && d2 < bestD2)) {
        best = { z2 }; bestD2 = d2;
      }
    }
    return best ? best.z2 : 0;
  }

  /**
   * screen → world convenience: unproject() at the bead-picked depth when
   * opts.snapToBead (and positions are available), else at opts.depth ?? 0.
   */
  screenToWorld(clientX, clientY, opts = {}) {
    const depth = opts.snapToBead ? this.pickDepth(clientX, clientY, opts.pos)
                                  : (opts.depth ?? 0);
    return this.unproject(clientX, clientY, depth);
  }

  /* ------------------------------ render ------------------------------ */

  /** Draw one frame; pos = Float64Array(3n) current bead coordinates. */
  render(pos) {
    if (pos) this._lastPos = pos;   // kept for pickDepth()/screenToWorld()
    const ctx = this.ctx;
    // Keep the backing buffer in sync with the CSS size every frame (covers
    // late layout, panel toggles and pinch-zoom without a resize event).
    if (this.canvas.clientWidth > 0 && this.canvas.clientHeight > 0) {
      const wantW = Math.round(this.canvas.clientWidth * (this.dpr || 1));
      const wantH = Math.round(this.canvas.clientHeight * (this.dpr || 1));
      if (wantW !== this.canvas.width || wantH !== this.canvas.height) this._resize();
    }
    const W = this.canvas.width, H = this.canvas.height;
    ctx.fillStyle = "#0a0c10";
    ctx.fillRect(0, 0, W, H);
    if (!pos || this.n === 0) return;

    // Scratch arrays are sized to the current bead count — reallocate if a
    // larger (or differently sized) system was built since the last render.
    if (!this._px || this._px.length !== this.n) {
      this._px = new Float64Array(this.n);
      this._py = new Float64Array(this.n);
      this._pz = new Float64Array(this.n);
    }

    const scale = (Math.min(W, H) * 0.45 * this.zoom) / this.radius;
    const [cx, cy, cz] = this.center;
    const sx = Math.sin(this.rotX), cxx = Math.cos(this.rotX);
    const sy = Math.sin(this.rotY), cyy = Math.cos(this.rotY);
    const fov = this.radius * 4; // perspective denominator

    // ---- render-only motion magnification ----------------------------------
    // A folded Cα ENM only exhibits sub-Å thermal jitter — invisible at 60 fps.
    // Amplify each particle's fluctuation around a slowly-relaxing running mean
    // so the dynamics are actually visible. Only the *displayed* coords are
    // amplified; physics (HUD, recorder) always uses the exact `pos`.
    let amp = pos;
    if (this.motionGain > 1) {
      if (!this._amp || this._amp.length !== this.n * 3) this._amp = new Float64Array(this.n * 3);
      if (!this._mean) this._mean = Float64Array.from(pos);
      const m = this._mean, a = this._amp, g = this.motionGain - 1, f = 0.05;
      const n3 = this.n * 3;
      for (let i = 0; i < n3; i++) {
        m[i] += f * (pos[i] - m[i]);          // EMA mean (≈ native state)
        a[i] = pos[i] + g * (pos[i] - m[i]);  // amplified fluctuation
      }
      amp = a;
    }

    // ---- project all beads (store in scratch arrays) -------------------
    const px = this._px, py = this._py, pz = this._pz;
    for (let i = 0; i < this.n; i++) {
      const x0 = amp[3 * i] - cx, y0 = amp[3 * i + 1] - cy, z0 = amp[3 * i + 2] - cz;
      // rotate Y then X
      const x1 = cyy * x0 + sy * z0;
      const z1 = -sy * x0 + cyy * z0;
      const y1 = cxx * y0 - sx * z1;
      const z2 = sx * y0 + cxx * z1;
      const persp = fov / (fov + z2);
      px[i] = W / 2 + x1 * scale * persp + this.panX;
      py[i] = H / 2 - y1 * scale * persp + this.panY;
      pz[i] = z2;
    }

    // ---- ENM contact sticks (behind everything) -------------------------
    if (this.showContacts && this.contactIdx) {
      ctx.lineWidth = this.dpr * 0.6;
      ctx.strokeStyle = "rgba(120,140,170,0.16)";
      ctx.beginPath();
      for (let k = 0; k < this.contactIdx.length; k += 2) {
        const i = this.contactIdx[k], j = this.contactIdx[k + 1];
        ctx.moveTo(px[i], py[i]);
        ctx.lineTo(px[j], py[j]);
      }
      ctx.stroke();
    }

    // holo contact springs — amber sticks showing the native protein–ligand pose
    if (this.holoIdx && this.showContacts) {
      ctx.lineWidth = this.dpr * 1.2;
      ctx.strokeStyle = "rgba(255,176,80,0.75)";
      ctx.beginPath();
      for (let k = 0; k < this.holoIdx.length; k += 2) {
        const i = this.holoIdx[k], j = this.holoIdx[k + 1];
        ctx.moveTo(px[i], py[i]);
        ctx.lineTo(px[j], py[j]);
      }
      ctx.stroke();
    }

    // ---- ligand covalent bonds (behind protein, above contacts) -------------
    if (this.ligandBonds) {
      ctx.lineWidth = this.dpr * 1.9;
      ctx.strokeStyle = "rgba(200,205,215,0.9)";
      ctx.beginPath();
      for (let k = 0; k < this.ligandBonds.length; k += 3) {
        const i = this.ligandBonds[k], j = this.ligandBonds[k + 1];
        ctx.moveTo(px[i], py[i]);
        ctx.lineTo(px[j], py[j]);
      }
      ctx.stroke();
    }

    // ---- backbone trace --------------------------------------------------
    ctx.lineWidth = this.dpr * 1.6;
    ctx.strokeStyle = "rgba(230,235,245,0.85)";
    ctx.beginPath();
    for (const [s, e] of this.segments) {
      ctx.moveTo(px[s], py[s]);
      for (let i = s + 1; i < e; i++) ctx.lineTo(px[i], py[i]);
    }
    ctx.stroke();

    // ---- beads (painter's algorithm, far → near) -------------------------
    const order = this._order;
    for (let i = 0; i < this.n; i++) order[i] = i;
    // insertion sort on depth keys (fast for mostly-sorted trajectories)
    order.sort((a, b) => pz[b] - pz[a]);

    const baseR = this.drawSpheres ? 4.0 * this.dpr * Math.sqrt(this.zoom) : 1.8 * this.dpr;
    for (const i of order) {
      const depth = Math.max(0.25, Math.min(1, 1 - pz[i] / (this.radius * 2.2)));
      const rScale = i >= this.nProt ? 0.62 : 1.0;
      const r = rScale * baseR * Math.max(0.5, fov / (fov + pz[i]));
      const [cr, cg, cb] = this.colors[i];
      const rr = Math.round(cr * depth + 20), gg = Math.round(cg * depth + 20), bb = Math.round(cb * depth + 20);

      ctx.beginPath();
      ctx.arc(px[i], py[i], r, 0, 6.2832);
      ctx.fillStyle = `rgb(${rr},${gg},${bb})`;
      ctx.fill();
      if (this.drawSpheres) {
        // cheap specular highlight
        ctx.beginPath();
        ctx.arc(px[i] - r * 0.35, py[i] - r * 0.35, r * 0.4, 0, 6.2832);
        ctx.fillStyle = "rgba(255,255,255,0.22)";
        ctx.fill();
      }
    }
  }
}
