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
    // Semantic mirrors (R2): system data re-synced in setSystem below.
    // View transform (center/radius/rotX/rotY/zoom/panX/panY) is NOT copied
    // here — live get/set accessors below delegate to fallback so interaction
    // (rotate/zoom/pan) and setSystem recompute never go stale (R3 issue 3).
    this.n = this.fallback.n;
    this.nProt = this.fallback.nProt;
    this.ligandStart = this.fallback.ligandStart;
    this.heavy = this.fallback.heavy;
    this.colors = this.fallback.colors;
    this.segments = this.fallback.segments;
    this.secStruct = this.fallback.secStruct;
    this.contactIdx = this.fallback.contactIdx;
    this.pocketCenter = this.fallback.pocketCenter;
    this.ligandBonds = this.fallback.ligandBonds;
    this.holoIdx = this.fallback.holoIdx;
    this.ref = this.fallback.ref;
  }

  // ---- live view-transform delegation (R3 issue 3) ----
  // Value copies of center/radius/rotX/rotY/zoom/panX/panY go stale: primitives
  // never update after fallback interaction, and center is reassigned by
  // Viewer.setSystem. Delegating accessors keep gl === fallback always.
  get center() { return this.fallback.center; }
  set center(v) { this.fallback.center = v; }
  get radius() { return this.fallback.radius; }
  set radius(v) { this.fallback.radius = v; }
  get rotX() { return this.fallback.rotX; }
  set rotX(v) { this.fallback.rotX = v; }
  get rotY() { return this.fallback.rotY; }
  set rotY(v) { this.fallback.rotY = v; }
  get zoom() { return this.fallback.zoom; }
  set zoom(v) { this.fallback.zoom = v; }
  get panX() { return this.fallback.panX; }
  set panX(v) { this.fallback.panX = v; }
  get panY() { return this.fallback.panY; }
  set panY(v) { this.fallback.panY = v; }

  // ---- delegated API (mirrors Viewer) ----

  setSystem(sel, ff) {
    const r = this.fallback.setSystem(sel, ff);
    // Semantic mirrors (R2) — re-sync references reassigned by Viewer.setSystem.
    // View transform needs no copy here: accessors above delegate live to
    // fallback, so setSystem recompute (center/radius) and preserved
    // rot/zoom/pan are visible via gl immediately (R3 issue 3).
    this.n = this.fallback.n;
    this.nProt = this.fallback.nProt;
    this.ligandStart = this.fallback.ligandStart;
    this.heavy = this.fallback.heavy;
    this.colors = this.fallback.colors;
    this.segments = this.fallback.segments;
    this.secStruct = this.fallback.secStruct;
    this.contactIdx = this.fallback.contactIdx;
    this.pocketCenter = this.fallback.pocketCenter;
    this.ligandBonds = this.fallback.ligandBonds;
    this.holoIdx = this.fallback.holoIdx;
    this.ref = this.fallback.ref;
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
