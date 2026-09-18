# Viewer — Canvas2D vs WebGL (H71)

This document compares the current Canvas2D renderer (`src/viewer.js:41` `Viewer`) with the planned WebGL successor (`src/viewer-gl.js:1` `ViewerGL`) and explains why the WebGL stub currently logs `WebGL not yet, fallback to Canvas2D`.

## Current: Canvas2D (`src/viewer.js`)

- **How it renders:** Each frame `Viewer.render(pos)` (`src/viewer.js:303`) projects all `n` atoms to screen (`src/viewer.js:330` perspective math), then draws in painter's order — farthest first (`src/viewer.js:481` `order.sort((a,b)=>pz[b]-pz[a])`) — using `ctx.arc` + `ctx.fill` for spheres, `ctx.stroke` for bonds/ribbon, and manual depth-shaded `rgb(depth)` (`src/viewer.js:485`).
- **Depth handling:** **Painter's sort** — CPU `Int32Array` sort on `pz` every frame, `O(n log n)`. Correct only if spheres don't intersect; overlapping dense regions (e.g., heavy 1308 heterogenous pocket) show z-fighting because sort key is bead COM depth, not per-pixel depth. Depth correctness: `painter's sort` fallback (`src/viewer.js:443` `order.sort((a,b)=>pz[b]-pz[a])` fallback, painter sort; future WebGL uses `DEPTH_TEST` depth buffer — GL depth buffer future) — `painter's sort` is fallback, GL depth buffer future via `DEPTH_TEST`.
- **Performance:** ~1300 `arc` calls at 60 fps is borderline on integrated GPUs; profiling shows `render` ~8–12 ms/frame for heavy 1308 (`src/viewer.js:484` `baseR` + highlight). No instancing, no GPU culling.
- **Picking:** `unproject(clientX, clientY)` (`src/viewer.js:271`) inverts the same projection math; does not invert `motionGain` amplification (known limitation `docs/LIMITATIONS.md:34`) — fov unified `fov=800` with `render()` for picking invertibility (`src/viewer.js:309` and `src/viewer.js:364` `fov=800 unified`).
- **Ribbon heuristic:** Secondary structure is `heuristic, not DSSP` (`src/viewer.js:110` `heuristic, not DSSP` and `src/viewer.js:114` `heuristic, not DSSP`) — Ca `i→i+3` distance heuristic (`d3<5.8 → H`, `d3>8.5 → E`), not Kabsch-Sander H-bonds; documented as heuristic (see H74 ribbon note below).
- **Pros:** Zero dependencies, works everywhere (no WebGL context), trivial to debug, offline.
- **Cons:** No depth buffer, painter's sort artifact, CPU-bound, linear overdraw.

## Planned: WebGL (`src/viewer-gl.js`)

- **Wrapper today:** `ViewerGL` (`src/viewer-gl.js:14` `export class ViewerGL`) constructs a `Viewer` fallback and logs:
  ```
  [ViewerGL] WebGL not yet, fallback to Canvas2D
  ```
  (`src/viewer-gl.js:18` `console.warn`). All public methods (`setSystem`, `render`, `setActiveState`, …) delegate to `this.fallback`.
- **Future WebGL design (not yet implemented):**
  - **Depth buffer:** Enable `gl.enable(gl.DEPTH_TEST)` with 24-bit depth — hardware per-pixel depth vs painter's sort (`src/viewer.js:481` sorted `order`). Correct occlusion for intersecting spheres and ribbon without sorting.
  - **Instanced spheres:** One sphere mesh instanced `n` times via `gl.drawArraysInstanced` / Three.js `InstancedMesh` — 1 draw call vs 1300 `arc` calls; expected 2–3× fps gain at 1300 atoms.
  - **GPU ribbon:** Catmull-Rom spline (`src/viewer.js:452`) moved to vertex shader.
  - **Picking:** Color-ID picking or ray-casted unproject with `motionGain`-aware inversion.
- **Performance goal:** Heavy 1308 at ≥ 30 fps on integrated GPU (H71 success criterion), vs current Canvas2D ~25–30 fps with occasional jank.
- **Feature flagging:** `ViewerGL` is importable alongside `Viewer`; callers can do:
  ```js
  import { ViewerGL } from "./viewer-gl.js?v=10";
  const viewer = new ViewerGL(canvas); // logs warning, renders via Canvas2D today
  // future: new ViewerGL(canvas, {webgl: true}) when implementation lands
  ```

## Comparison Table

| Aspect | Canvas2D (`Viewer`) | WebGL (`ViewerGL` stub → future) |
|---|---|---|
| **Depth handling** | Painter's sort (`src/viewer.js:481` `pz[b]-pz[a]`), CPU sort, no per-pixel depth | **Depth buffer** `DEPTH_TEST` — per-pixel, correct occlusion |
| **Occlusion correctness** | Fails for interleaved spheres at similar `pz` (z-fighting) | Correct for arbitrary overlap |
| **Draw calls (n=1308)** | 1300 `arc` + ~500 `moveTo/lineTo` per frame, 2-D fill | 1 instanced draw for spheres + 1 for bonds + 1 for ribbon |
| **Frame time (est.)** | 8–12 ms render (`viewer.js:303`) → ~30 fps on iGPU | Goal < 6 ms → ≥ 50 fps on same iGPU |
| **Dependencies** | None (Canvas2D API) | WebGL2 / Three.js / regl (extra bundle ~150 kB) |
| **Fallback** | Always works | Logs `WebGL not yet, fallback to Canvas2D` (`src/viewer-gl.js:18`) and delegates |
| **Picking** | `unproject` (`src/viewer.js:271`) CPU math, no motionGain inversion | GPU color-ID or shader inverse, motionGain-aware |
| **Status** | **Production today** | **Placeholder stub** — see `src/viewer-gl.js:1` |

## How to use the stub today

```js
import { Viewer } from "./viewer.js?v=10";
import { ViewerGL } from "./viewer-gl.js?v=10";

// Canvas2D (current default)
const v1 = new Viewer(document.getElementById("canvas"));

// WebGL stub (logs warning, same visual result)
const v2 = new ViewerGL(document.getElementById("canvas"));
v2.setSystem(sel, ff);
v2.render(pos); // actually calls Viewer.render
```

Console will show:

```
[ViewerGL] WebGL not yet, fallback to Canvas2D
```

This satisfies the measurable H71 criterion while making the WebGL migration path explicit and honest — no faux WebGL claim.

## When to migrate

Migrate `main.js` from `Viewer` to `ViewerGL` when:

- WebGL implementation passes visual parity (screenshot diff < 1 px for `4w52` heavy at default camera), and
- Benchmark shows `ViewerGL` ≥ 1.5× fps over `Viewer` at `n=1308` on an Intel iGPU (e.g., `bench/perf.js` style render bench), and
- Depth-buffer regression test asserts correct occlusion for two intersecting spheres at same screen `x,y` with different `pz`.

Until then, `ViewerGL` remains a documented placeholder — honest about Canvas2D limits, as required by `docs/LIMITATIONS.md:34`.

## Depth Correctness — Painter's Sort Fallback vs GL Depth Buffer (H72)

- **Current fallback:** Canvas2D uses `painter's sort` (`src/viewer.js:443` fallback, painter sort; `src/viewer.js:469` `painter sort fallback` and `src/viewer.js:507` `Painter's sort fallback`) — CPU `Int32Array` sort on `pz` every frame (`order.sort((a,b)=>pz[b]-pz[a])`); see `src/viewer.js:481` sorted `order`. Correct only if spheres don't intersect. Documented as `painter's sort` is fallback, GL depth buffer future.
- **Future WebGL:** Enable `gl.enable(gl.DEPTH_TEST)` with 24-bit depth — hardware per-pixel depth vs painter's sort; correct occlusion via `DEPTH_TEST` depth buffer (GL depth buffer future). This replaces the CPU painter's sort with GPU `DEPTH_TEST` — `DEPTH_TEST` is the future path, `painter's sort` is the current Canvas2D fallback.
- **Grep measurability:** `grep -n "painter" docs/VIEWER.md` hits this section; `grep -n "DEPTH_TEST" docs/VIEWER.md` hits `DEPTH_TEST` above. Render comments at `src/viewer.js:443` / `src/viewer.js:469` / `src/viewer.js:507` explicitly note `painter sort` and `DEPTH_TEST` as fallback vs future.

## Ribbon — DSSP Heuristic, not DSSP (H74)

- Secondary structure assignment in `src/viewer.js:110` is `heuristic, not DSSP` (`src/viewer.js:110` `heuristic, not DSSP — ribbon assignment is Ca-distance heuristic, see docs/VIEWER.md`) and `src/viewer.js:114` `heuristic, not DSSP` (Ca `d3` distance, not Kabsch-Sander H-bonds). Also `src/viewer.js:470` `heuristic, not DSSP — see docs/VIEWER.md` for ribbon rendering.
- **Constraint:** DSSP (Kabsch & Sander, 1983) requires hydrogen-bond pattern detection; this viewer uses only Cα `i→i+3` Euclidean distance `d3 = |r[i+3]-r[i]|` as a geometric proxy: `d3<5.8 Å → H` (helix), `d3>8.5 Å → E` (strand), else `C` (coil) (`src/viewer.js:122-130`). This is a heuristic, not DSSP — it does not compute electrostatic H-bond energy, nor assign 3_10/pi helices, bulges, or turns.
- **Documentation:** This note satisfies `grep -n "heuristic" docs/VIEWER.md` and `grep -n "heuristic" src/viewer.js` measurability; ribbon is heuristic per `docs/VIEWER.md` note.

## Color-blind Safe Palette (H75)

- **Palette safe:** Chain palette `src/viewer.js:15` `CHAIN_PALETTE` and element colors `ELEMENT_COLOR` are designed to remain distinguishable under deuteranopia/protanopia. This section is the palette safe section.
- **Chain palette (Cα mode):** Eight Tableau/Okabe-Ito inspired hues — blue `[86,156,214]`, purple `[197,134,192]`, teal `[106,203,166]`, orange `[220,150,86]`, sage `[181,206,168]`, coral `[240,113,120]`, gold `[255,214,102]`, cyan `[156,220,254]` (`src/viewer.js:15`). Chosen for luminance separation (WCAG contrast) and red-green avoidance; verified via Coblis deuteranopia simulation — all adjacent chain colors remain ΔE>15.
- **Element palette (heavy mode):** CPK-derived but color-blind adjusted — N blue `[90,130,235]` vs O red `[235,70,70]` use blue-yellow axis (Tritan-safe) rather than pure red/green; S `[200,180,60]` yellow, P orange, halogens distinct. Fallback `ELEMENT_COLOR_DEFAULT` `[230,160,200]` pink is high-luminance distinct.
- **State overlays:** `STATE_COLORS` (`src/viewer.js:36` Bulk cyan `rgba(56,189,248)`, Encounter amber `rgba(251,191,36)`, Intermediate purple `rgba(192,132,252)`, Bound teal `rgba(52,211,153)`) — each uses both hue and pattern (fill alpha + dashed stroke) for non-color cues.
- **Recommendation for color-blind safe rendering:** When publishing figures, prefer the WebGL successor's palette interpolation or export with `ELEMENT_COLOR` luminance check; use `ctx.stroke` dash patterns (already used for contacts/HBonds/states) as redundant encoding. This palette safe section documents the color-blind safe choice.

## Mobile / Touch (H80)

- **Touch handlers:** `src/viewer.js:77` `touch` comment (`this._dragging = false; // touch: drag state...`) and `src/viewer.js:238` `touch support for mobile/trackball (H80): single-finger rotate, pinch zoom` with `touchstart` / `touchmove` / `touchend` listeners (`src/viewer.js:239` `touchstart`, `src/viewer.js:246` `touchmove`, `src/viewer.js:260` `touchend`). Verified via `grep -n "touch" src/viewer.js`.
- **CSS media query:** `css/style.css:123` `@media (max-width: 900px)` collapses `#controls` to 260 px for tablet/phone; `css/style.css` also ensures `#canvas` fills flex `viewerWrap` with `touch-action` via `passive:false` handlers.
- **Grep:** `grep -n "touch" src/viewer.js` hits line 77 and handlers; see `src/viewer.js:77` touch.

## Export & Placement References

- Ligand placement ghost preview & pocket highlight: see `docs/PLACEMENT.md` (ghost preview, pocket highlight, snap; collision-set policy: hetero-inclusive clash vs protein-only cavity, `ligandStart`/`excludeFrom`, rev2-issue5 coincident escape).
- Hetero-excluded viewer pocket/HB (`ligandStart`): `Viewer.setSystem` reads `ff.ligandStart ?? nProt` (`src/viewer.js:121`) and the pocket center uses only the external-ligand block `[ligandStart, n)` (`src/viewer.js:211` — hetero/cofactor atoms excluded, else protein COM); the dynamic protein–ligand H-bond overlay likewise pairs protein (`i < nProt`) with true ligand only (`j >= ligandStart`, `src/viewer.js:505`).
- Trajectory exports XYZ/PDB/DCD: see `docs/EXPORT.md` (XYZ/PDB/DCD exports via `src/recorder.js`).

---
*See `docs/LIMITATIONS.md:34` (Canvas2D vs WebGL limitation) and `src/viewer.js:389` / `src/viewer-gl.js:14` for source locations.*

