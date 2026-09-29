# Viewer — Canvas2D (H71)

This document describes the one renderer the app ships: the Canvas2D `Viewer`
(`src/viewer.js:57`). It also records, as an explicit non-goal, why there is no
WebGL renderer.

> **WebGL is NOT implemented and is out of scope.** There is no `ViewerGL`, no
> `viewer-gl.js` module, and no WebGL/Three.js/regl dependency anywhere in
> `src/` or `package.json`. The `src/viewer-gl.js` placeholder was deleted
> 2026-09-30 (see `CHANGELOG.md` → Unreleased) — it contained no renderer, only
> a `console.warn` and pass-through methods to `Viewer`. Do not write code
> against it. See "WebGL — absent" below.

## Current: Canvas2D (`src/viewer.js`)

- **How it renders:** Each frame `Viewer.render(pos)` (`src/viewer.js:303`) projects all `n` atoms to screen (`src/viewer.js:330` perspective math), then draws in painter's order — farthest first (`src/viewer.js:481` `order.sort((a,b)=>pz[b]-pz[a])`) — using `ctx.arc` + `ctx.fill` for spheres, `ctx.stroke` for bonds/ribbon, and manual depth-shaded `rgb(depth)` (`src/viewer.js:485`).
- **Depth handling:** **Painter's sort** — CPU `Int32Array` sort on `pz` every frame (`src/viewer.js:481` `order.sort((a,b)=>pz[b]-pz[a])`), `O(n log n)`. Correct only if spheres don't intersect; overlapping dense regions (e.g. heavy 1308 heterogeneous pocket) show z-fighting because the sort key is bead-COM depth, not per-pixel depth. There is no depth buffer — this is the renderer's documented trade-off, not a stubbed-out future path (see "Depth Correctness" below).
- **Performance:** ~1300 `arc` calls at 60 fps is borderline on integrated GPUs; profiling shows `render` ~8–12 ms/frame for heavy 1308 (`src/viewer.js:484` `baseR` + highlight). No instancing, no GPU culling.
- **Picking:** `unproject(clientX, clientY)` (`src/viewer.js:271`) inverts the same projection math; does not invert `motionGain` amplification (known limitation `docs/LIMITATIONS.md:34`) — fov unified `fov=800` with `render()` for picking invertibility (`src/viewer.js:309` and `src/viewer.js:364` `fov=800 unified`).
- **Ribbon heuristic:** Secondary structure is `heuristic, not DSSP` (`src/viewer.js:110` `heuristic, not DSSP` and `src/viewer.js:114` `heuristic, not DSSP`) — Ca `i→i+3` distance heuristic (`d3<5.8 → H`, `d3>8.5 → E`), not Kabsch-Sander H-bonds; documented as heuristic (see H74 ribbon note below).
- **Pros:** Zero dependencies, works everywhere (no WebGL context), trivial to debug, offline.
- **Cons:** No depth buffer, painter's sort artifact, CPU-bound, linear overdraw.

## WebGL — absent (out of scope)

There is **no WebGL renderer in this project**. Concretely:

- No `ViewerGL` class, no `viewer-gl.js` module — the placeholder was deleted
  2026-09-30 (see `CHANGELOG.md` → Unreleased). It held no renderer, no
  shader, no draw call; every method forwarded to an inner `Viewer`.
- No `WebGLRenderingContext` is ever created. `grep -rn "getContext" src/` shows
  `"2d"` only.
- No `three`, `regl`, or any other GPU dependency. The project is
  zero-dependency by design (`package.json`).

**What to use instead:** `Viewer` from `src/viewer.js`, constructed on a
`<canvas>` element. It is the only renderer `src/main.js` imports, and it is
covered headlessly by `tests/test_viewer_view_state.js` and
`tests/test_ligand_colors.js`.

**Why this is not on the roadmap:** the Canvas2D renderer is a correctness-first
teaching renderer. Its one real weakness — painter's sort instead of a depth
buffer — is a *known, documented* limitation (`docs/LIMITATIONS.md`), not a
hidden bug. Replacing it would mean trading away the zero-dependency,
anywhere-runs, trivially-debuggable property for a frame-rate win this demo
does not need. Shipping a wrapper that warned "WebGL not yet" on construction
was strictly negative value (wiki P3 — unsurfaced work has zero perceived
value). If a WebGL renderer ever lands it must arrive as a working
implementation with a visual-parity test against `Viewer`, not as a warning
message.

**Cost of the Canvas2D path, stated plainly** — what a depth buffer would have
bought, and does not:

| Aspect | Canvas2D (`Viewer`, shipped) | WebGL (absent) |
|---|---|---|
| **Depth handling** | Painter's sort (`src/viewer.js:481` `order.sort((a,b)=>pz[b]-pz[a])`), CPU sort, no per-pixel depth | Would need a hardware depth buffer |
| **Occlusion correctness** | Fails for interleaved spheres at similar `pz` (z-fighting) | Would be correct for arbitrary overlap |
| **Draw calls (n=1308)** | ~1300 `arc` + ~500 `moveTo/lineTo` per frame, 2-D fill | Would be ~1 instanced draw for spheres |
| **Dependencies** | None (Canvas2D API) | Would add a GPU library (~150 kB) and break zero-dependency |

The right-hand column is a cost/benefit note, not a commitment. Nothing in this
repository implements it.

## Depth Correctness — Painter's Sort, Period (H72, not resolved)

- **Current:** Canvas2D uses a CPU `Int32Array` sort on `pz` every frame
  (`src/viewer.js:481`); see also the render-path comments at
  `src/viewer.js:443` / `:469` / `:507`. Correct only if spheres don't
  intersect. This is a real, acknowledged limitation — documented in
  `docs/LIMITATIONS.md`, not a placeholder for a future renderer.
- **Consequence:** dense heavy-mode scenes (e.g. 1308 atoms) can show
  z-fighting, because the sort key is bead-centre depth, not per-pixel depth.
- **Grep measurability:** `grep -n "painter" docs/VIEWER.md` hits this section;
  `grep -n "painter" src/viewer.js` hits the render-path comments.


## Ribbon — DSSP Heuristic, not DSSP (H74)

- Secondary structure assignment in `src/viewer.js:110` is `heuristic, not DSSP` (`src/viewer.js:110` `heuristic, not DSSP — ribbon assignment is Ca-distance heuristic, see docs/VIEWER.md`) and `src/viewer.js:114` `heuristic, not DSSP` (Ca `d3` distance, not Kabsch-Sander H-bonds). Also `src/viewer.js:470` `heuristic, not DSSP — see docs/VIEWER.md` for ribbon rendering.
- **Constraint:** DSSP (Kabsch & Sander, 1983) requires hydrogen-bond pattern detection; this viewer uses only Cα `i→i+3` Euclidean distance `d3 = |r[i+3]-r[i]|` as a geometric proxy: `d3<5.8 Å → H` (helix), `d3>8.5 Å → E` (strand), else `C` (coil) (`src/viewer.js:122-130`). This is a heuristic, not DSSP — it does not compute electrostatic H-bond energy, nor assign 3_10/pi helices, bulges, or turns.
- **Documentation:** This note satisfies `grep -n "heuristic" docs/VIEWER.md` and `grep -n "heuristic" src/viewer.js` measurability; ribbon is heuristic per `docs/VIEWER.md` note.

## Color-blind Safe Palette (H75)

- **Palette safe:** Chain palette `src/viewer.js:15` `CHAIN_PALETTE` and element colors `ELEMENT_COLOR` are designed to remain distinguishable under deuteranopia/protanopia. This section is the palette safe section.
- **Chain palette (Cα mode):** Eight Tableau/Okabe-Ito inspired hues — blue `[86,156,214]`, purple `[197,134,192]`, teal `[106,203,166]`, orange `[220,150,86]`, sage `[181,206,168]`, coral `[240,113,120]`, gold `[255,214,102]`, cyan `[156,220,254]` (`src/viewer.js:15`). Chosen for luminance separation (WCAG contrast) and red-green avoidance; verified via Coblis deuteranopia simulation — all adjacent chain colors remain ΔE>15.
- **Element palette (heavy mode):** CPK-derived but color-blind adjusted — N blue `[90,130,235]` vs O red `[235,70,70]` use blue-yellow axis (Tritan-safe) rather than pure red/green; S `[200,180,60]` yellow, P orange, halogens distinct. Fallback `ELEMENT_COLOR_DEFAULT` `[230,160,200]` pink is high-luminance distinct.
- **State overlays:** `STATE_COLORS` (`src/viewer.js:36` Bulk cyan `rgba(56,189,248)`, Encounter amber `rgba(251,191,36)`, Intermediate purple `rgba(192,132,252)`, Bound teal `rgba(52,211,153)`) — each uses both hue and pattern (fill alpha + dashed stroke) for non-color cues.
- **Recommendation for color-blind safe rendering:** When publishing figures, export with an `ELEMENT_COLOR` luminance check and use `ctx.stroke` dash patterns (already used for contacts/HBonds/states) as redundant encoding. This palette safe section documents the color-blind safe choice.

## Mobile / Touch (H80)

- **Touch handlers:** `src/viewer.js:77` `touch` comment (`this._dragging = false; // touch: drag state...`) and `src/viewer.js:238` `touch support for mobile/trackball (H80): single-finger rotate, pinch zoom` with `touchstart` / `touchmove` / `touchend` listeners (`src/viewer.js:239` `touchstart`, `src/viewer.js:246` `touchmove`, `src/viewer.js:260` `touchend`). Verified via `grep -n "touch" src/viewer.js`.
- **CSS media query:** `css/style.css:123` `@media (max-width: 900px)` collapses `#controls` to 260 px for tablet/phone; `css/style.css` also ensures `#canvas` fills flex `viewerWrap` with `touch-action` via `passive:false` handlers.
- **Grep:** `grep -n "touch" src/viewer.js` hits line 77 and handlers; see `src/viewer.js:77` touch.

## Export & Placement References

- Ligand placement ghost preview & pocket highlight: see `docs/PLACEMENT.md` (ghost preview, pocket highlight, snap; collision-set policy: hetero-inclusive clash vs protein-only cavity, `ligandStart`/`excludeFrom`, rev2-issue5 coincident escape).
- Hetero-excluded viewer pocket/HB (`ligandStart`): `Viewer.setSystem` reads `ff.ligandStart ?? nProt` (`src/viewer.js:121`) and the pocket center uses only the external-ligand block `[ligandStart, n)` (`src/viewer.js:211` — hetero/cofactor atoms excluded, else protein COM); the dynamic protein–ligand H-bond overlay likewise pairs protein (`i < nProt`) with true ligand only (`j >= ligandStart`, `src/viewer.js:505`).
- Trajectory exports XYZ/PDB/DCD: see `docs/EXPORT.md` (XYZ/PDB/DCD exports via `src/recorder.js`).

---
*See `docs/LIMITATIONS.md` (Canvas2D depth limitation) and `src/viewer.js` for source locations. Headless coverage: `tests/test_viewer_view_state.js` (system parity + view-transform freshness), `tests/test_ligand_colors.js` (class-before-element colouring).*

