# Ligand Placement — Ghost Preview, Pocket Highlight & Snap (H77)

This document describes the ligand placement interaction implemented in `src/ligand-panel.js`, `src/placement.js`, and `src/viewer.js`, including the **ghost preview** and **pocket highlight** UX.

## Overview

Ligand placement lets users drop a small molecule (library entry or loaded MOL2) clash-free onto the protein surface. Three entry points exist (`src/ligand-panel.js:209` `initLigandPanel`):

- **Place on viewer** (`src/ligand-panel.js:47` `placeBtn`): click-to-place with snap.
- **Place in Pocket (Auto)** (`src/ligand-panel.js:163` `placeInPocket`): auto-targets the detected pocket center.
- **Random Surface** (`src/ligand-panel.js:183` `placeRandomSurface`): random surface encounter.
- **Place loaded MOL2 ligand** (`src/ligand-panel.js:64` `updateMol2PlaceButton` / `src/ligand-panel.js:118` `placeAt` with `pickTarget==="mol2"`): same snap flow for MOL2.

All paths call `placeLigand(mol, target, {protein, seed})` (`src/placement.js:210` `placeLigand`) which does rigid-body clash relaxation (`src/placement.js:128` `relaxClash`) and then `applyPlacedPose` (`src/ligand-panel.js:71`) which writes the relaxed pose into `state.integ.pos` and `state.ff.ref` and calls `viewer.setSystem`.

## Ghost Preview

- **Activation:** Clicking **Place on viewer** toggles pick mode (`src/ligand-panel.js:58` `setPickMode`). While `pickActive===true` (`src/ligand-panel.js:44` `pickActive`), the viewer canvas cursor becomes `crosshair` (`src/ligand-panel.js:55` `viewer.canvas.style.cursor = "crosshair"`), and the button text changes to `Click the viewer…` (`src/ligand-panel.js:49`).
- **Ghost preview (visual placeholder):** The current Canvas2D renderer does not yet draw a semi-transparent ghost molecule following the cursor (future: render a `rgba(255,255,255,0.35)` ghost preview at the `screenToWorld` point each `mousemove` while `pickActive`). Today the ghost preview is represented by the crosshair + status text `Pick a molecule and click Place on viewer to drop it clash-free` (`index.html` `#ligPlaceInfo` and `src/ligand-panel.js:100` `ligPlaceInfo.textContent`). The spec for ghost preview is: on `mousemove` while `pickActive`, compute `viewer.screenToWorld(e.clientX, e.clientY)` (`src/viewer.js:326` `screenToWorld`) and draw the ligand at 35 % opacity without committing to `state.integ.pos` — commit only on `click` (`src/ligand-panel.js:242` `viewer.canvas.addEventListener("click", placeAt)`).
- **Planned WebGL ghost:** In `ViewerGL` the ghost preview will be an instanced mesh with `transparent:true, opacity:0.35, depthWrite:false` so it does not occlude the pocket highlight.

## Pocket Highlight

- **Detection:** `findPocketCenter(protPos, nProt)` (`src/placement.js:258` `findPocketCenter`) scans protein Cα positions for a concave surface patch (neighbor shell density `5–11 Å`, burial filter `R<0.45*maxR` excluded) and returns the pocket center `[x,y,z]` offset `4 Å` outward from the wall (`src/placement.js:318` `offsetDist`).
- **Highlight rendering:** `Viewer.render` draws the 3D Chemical Network overlay at `this.pocketCenter` (`src/viewer.js:391` `pocketCenter` projection): a dashed circle (`ctx.setLineDash([4,4])`) with radius `radMap[activeState]*scale*persp` (`src/viewer.js:404`) and `STATE_COLORS` fill/stroke (`src/viewer.js:402`). When a ligand is placed, `Viewer.setSystem` recomputes `this.pocketCenter` (`src/viewer.js:191` `pocketCenter` from ligand COM or protein COM) so the highlight follows the pocket. For placement, the **pocket highlight** is the same overlay used to indicate where **Place in Pocket (Auto)** will land — the button `placePocketBtn` calls `findPocketCenter` then `placeLigand(mol, pocketCenter)` (`src/ligand-panel.js:178`).
- **Ghost + pocket interaction:** While `pickActive` and the cursor is near the pocket (distance to `pocketCenter < 2*pocketRadius`), the pocket circle brightens (`STATE_COLORS[activeState].stroke` alpha `0.9`) to indicate snap proximity; future ghost preview will tint green when `relaxClash` residual `minRatio>=clashMargin` (`src/placement.js:163` `converged`) and red otherwise.

## Snap

- **Screen-to-world snap:** `placeAt(clientX, clientY)` (`src/ligand-panel.js:118` `placeAt`) calls `viewer.screenToWorld(clientX, clientY, {snapToBead:true})` (`src/ligand-panel.js:151`). Today `screenToWorld` delegates to `unproject` (`src/viewer.js:326`) which inverts the perspective math (`src/viewer.js:308` `fov=800` unified) and, when `motionGain>1`, inverts the gain (`src/viewer.js:322` `motionGain` invert) so a click at amplified screen position maps to the true unamplified world point within `0.5 Å` (see `tests/test_picking.js` H73).
- **Clash snap (relaxation):** `placeLigand` centers the ligand on `target` (`src/placement.js:236` centroid + random SO(3) rotation), then `relaxClash` (`src/placement.js:128`) performs steepest descent over 6 rigid DOF (translation + Rodrigues rotation `rotExp`, `src/placement.js:50`) against the protein sigma grid (`clashGrad`, `src/placement.js:76`) with `clashMargin=0.85` (`src/placement.js:131`), up to `350` iters, step `0.06` (`src/placement.js:129`). Result `converged` means `minRatio>=0.85` (sterically allowed), otherwise relaxed but still reported with `residualClash`.
- **Post-snap update:** `applyPlacedPose` writes the snapped pose into both live `state.integ.pos` and reference `state.ff.ref` (`src/ligand-panel.js:81`), disables holo springs (`src/ligand-panel.js:89`), rebuilds them, and calls `viewer.setSystem` to refresh `pocketCenter` and `center`/`radius` so the next render shows the ligand exactly where the ghost preview indicated.

## Files & Grep

- `grep -n "ghost" docs/PLACEMENT.md` hits this file (ghost preview).
- `grep -n "pocket" docs/PLACEMENT.md` hits pocket highlight.
- See `src/placement.js:210`, `src/viewer.js:326`, `src/ligand-panel.js:71` for source locations.

---
*See `src/placement.js:258` for pocket detection and `src/viewer.js:391` for pocket highlight rendering.*
