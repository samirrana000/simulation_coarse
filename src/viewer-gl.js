/**
 * viewer-gl.js — WebGL prototype stub (H71).
 *
 * This is an honest placeholder: the app currently renders with Canvas2D
 * (`src/viewer.js:41` `Viewer` — painter's sort, no depth buffer).
 * ViewerGL wraps Viewer and logs a single warning so callers can
 * feature-flag WebGL without breaking:
 *   "WebGL not yet, fallback to Canvas2D"
 *
 * Future implementation will use WebGL2 / Three.js / regl with:
 *   - depth buffer (correct hidden-surface removal vs painter's sort)
 *   - instanced sphere geometry for 1308 heavy atoms
 *   - GPU-based ribbon generation
 * Until then, all rendering delegates to Canvas2D.
 */

import { Viewer } from "./viewer.js?v=10";

export class ViewerGL {
  /**
   * @param {HTMLCanvasElement} canvas
   * @param {object} [opts]
   */
  constructor(canvas, opts = {}) {
    // Intentionally log exactly this string — measurable criterion for H71
    console.warn("[ViewerGL] WebGL not yet, fallback to Canvas2D");
    this.canvas = canvas || null;
    this.opts = opts;
    // Delegate to the proven Canvas2D viewer
    this.fallback = new Viewer(canvas);
    // Mirror the Viewer public API so drop-in replacement is safe
    this.n = this.fallback.n;
    this.nProt = this.fallback.nProt;
    this.center = this.fallback.center;
    this.radius = this.fallback.radius;
    this.rotX = this.fallback.rotX;
    this.rotY = this.fallback.rotY;
    this.zoom = this.fallback.zoom;
    this.panX = this.fallback.panX;
    this.panY = this.fallback.panY;
  }

  // ---- delegated API (mirrors Viewer) ----

  setSystem(sel, ff) {
    const r = this.fallback.setSystem(sel, ff);
    this.n = this.fallback.n;
    this.nProt = this.fallback.nProt;
    this.center = this.fallback.center;
    this.radius = this.fallback.radius;
    return r;
  }

  render(pos) {
    // Future: bind WebGL framebuffer, depth test, instanced draw.
    // Today: fallback to Canvas2D painter's sort (correct for demo, no depth buffer).
    return this.fallback.render(pos);
  }

  setActiveState(s) { return this.fallback.setActiveState(s); }
  setMotionGain(g) { return this.fallback.setMotionGain(g); }
  clear() { return this.fallback.clear(); }
  unproject(x, y, d) { return this.fallback.unproject(x, y, d); }
  screenToWorld(x, y, opts) { return this.fallback.screenToWorld(x, y, opts); }

  // Passthrough for viewer flags so UI code works unchanged
  get showContacts() { return this.fallback.showContacts; }
  set showContacts(v) { this.fallback.showContacts = v; }
  get drawSpheres() { return this.fallback.drawSpheres; }
  set drawSpheres(v) { this.fallback.drawSpheres = v; }
  get drawRibbon() { return this.fallback.drawRibbon; }
  set drawRibbon(v) { this.fallback.drawRibbon = v; }
  get showHBonds() { return this.fallback.showHBonds; }
  set showHBonds(v) { this.fallback.showHBonds = v; }
  get showStates() { return this.fallback.showStates; }
  set showStates(v) { this.fallback.showStates = v; }

  // Resize / input delegated via fallback's constructor bindings
  _resize() { return this.fallback._resize(); }
}
