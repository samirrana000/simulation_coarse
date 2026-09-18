# Ligand Placement — Ghost Preview, Pocket Highlight & Snap (H77)

This document describes the ligand placement interaction implemented in `src/ligand-panel.js`, `src/placement.js`, and `src/viewer.js`, including the **ghost preview** and **pocket highlight** UX.

## Overview

Ligand placement lets users drop a small molecule (library entry or loaded MOL2) clash-free onto the protein surface. Three entry points exist (`src/ligand-panel.js:305` `initLigandPanel`):

- **Place on viewer** (`src/ligand-panel.js:312` `placeBtn`): click-to-place with snap.
- **Place in Pocket (Auto)** (`src/ligand-panel.js:223` `placeInPocket`): auto-targets the detected pocket center.
- **Random Surface** (`src/ligand-panel.js:267` `placeRandomSurface`): random surface encounter.
- **Place loaded MOL2 ligand** (`src/ligand-panel.js:66` `updateMol2PlaceButton` / `src/ligand-panel.js:150` `placeAt` with `pickTarget==="mol2"`): same snap flow for MOL2.

All paths call `placeLigand(mol, target, {protein, seed})` (`src/placement.js:238` `placeLigand`) which does rigid-body clash relaxation (`src/placement.js:149` `relaxClash`) and then `applyPlacedPose` (`src/ligand-panel.js:73`) which writes the relaxed pose into `state.integ.pos` and `state.ff.ref` and calls `viewer.setSystem`.

## Collision set — hetero-inclusive clash vs protein-only cavity (rev1-issue2)

- **Clash set:** `getProteinCoordsAndSigma` (`src/ligand-panel.js:116`) builds indices `[0, collEnd)` with `collEnd = ligandStart ?? nProt` (`src/ligand-panel.js:127`). In heavy mode this is protein + hetero (PDB ligands/cofactors/metals, per-atom sigma from `ff._elem`), EXCLUDING the incoming ligand slot `[ligandStart, n)`; in CG mode (`ForceField` has no `ligandStart`) it stays protein-only via `_protSigma` (unchanged legacy behavior).
- **Slot exclusion:** `relaxClash` additionally honors `opts.excludeFrom` / `protein.excludeFrom` (`src/placement.js:157` — indices `[excludeFrom, nTotal)` are skipped), so a stale incoming-slot atom sitting on the target never distorts the pose.
- **Cavity vs clash split:** `placeInPocket` searches the cavity protein-only — `findPocketCenter(protein.pos, state.ff.nProt)` (`src/ligand-panel.js:251`) — while the relaxation clashes against the extended hetero-inclusive `protein` set (`src/ligand-panel.js:243`). Hetero atoms shape the escape, never the cavity vote.
- **Coincident escape (rev2-issue5):** an exactly coincident ligand/collider pair (`r2 < 1e-6`, `src/placement.js:89`) escapes along a deterministic index-hashed direction per `(a,i)` with its torque arm (`src/placement.js:89-104`), NOT the old fixed `+x+y+z` diagonal (brittle when `+x+y+z` is walled in a dense hetero shell). Fixed given indices, so placement stays deterministic given `seed` (seed owns the initial SO(3) rotation, `src/placement.js:246`); the `minRatio` floor (`1e-3/rE`) keeps `residualClash` (`src/placement.js:220`) finite instead of `Infinity`.

## Ghost Preview

- **Activation:** Clicking **Place on viewer** toggles pick mode (`src/ligand-panel.js:60` `setPickMode`). While `pickActive===true` (`src/ligand-panel.js:46` `pickActive`), the viewer canvas cursor becomes `crosshair` (`src/ligand-panel.js:57` `viewer.canvas.style.cursor = "crosshair"`), and the button text changes to `Click the viewer…` (`src/ligand-panel.js:50`).
- **Ghost preview (visual placeholder):** The current Canvas2D renderer does not yet draw a semi-transparent ghost molecule following the cursor (future: render a `rgba(255,255,255,0.35)` ghost preview at the `screenToWorld` point each `mousemove` while `pickActive`). Today the ghost preview is represented by the crosshair + status text `Pick a molecule and click Place on viewer to drop it clash-free` (`index.html` `#ligPlaceInfo` and `src/ligand-panel.js:103` `ligPlaceInfo.textContent`). The spec for ghost preview is: on `mousemove` while `pickActive`, compute `viewer.screenToWorld(e.clientX, e.clientY)` (`src/viewer.js:382` `screenToWorld`) and draw the ligand at 35 % opacity without committing to `state.integ.pos` — commit only on `click` (`src/ligand-panel.js:339` `viewer.canvas.addEventListener("click", placeAt)`).
- **Planned WebGL ghost:** In `ViewerGL` the ghost preview will be an instanced mesh with `transparent:true, opacity:0.35, depthWrite:false` so it does not occlude the pocket highlight.

## Pocket Highlight

- **Detection:** `findPocketCenter(protPos, nProt)` (`src/placement.js:289` `findPocketCenter`) scans protein Cα positions for a concave surface patch (neighbor shell density `5–11 Å`, burial filter `R<0.45*maxR` excluded) and returns the pocket center `[x,y,z]` offset `4 Å` outward from the wall (`src/placement.js:349` `offsetDist`).
- **Highlight rendering:** `Viewer.render` draws the 3D Chemical Network overlay at `this.pocketCenter` (`src/viewer.js:461` `pocketCenter` projection): a dashed circle (`ctx.setLineDash([4,4])`) with radius `radMap[activeState]*scale*persp` (`src/viewer.js:473`) and `STATE_COLORS` fill/stroke (`src/viewer.js:472`). When a ligand is placed, `Viewer.setSystem` recomputes `this.pocketCenter` (`src/viewer.js:211` `pocketCenter` from external-ligand COM, hetero excluded via `ligandStart`, else protein COM) so the highlight follows the pocket. For placement, the **pocket highlight** is the same overlay used to indicate where **Place in Pocket (Auto)** will land — the button `placePocketBtn` calls `findPocketCenter` then `placeLigand(mol, pocketCenter)` (`src/ligand-panel.js:259`).
- **Ghost + pocket interaction:** While `pickActive` and the cursor is near the pocket (distance to `pocketCenter < 2*pocketRadius`), the pocket circle brightens (`STATE_COLORS[activeState].stroke` alpha `0.9`) to indicate snap proximity; future ghost preview will tint green when `relaxClash` residual `minRatio>=clashMargin` (`src/placement.js:185` `converged`) and red otherwise.

## Snap

- **Screen-to-world snap:** `placeAt(clientX, clientY)` (`src/ligand-panel.js:150` `placeAt`) calls `viewer.screenToWorld(clientX, clientY, {snapToBead:true})` (`src/ligand-panel.js:197`). Today `screenToWorld` delegates to `unproject` (`src/viewer.js:382`) which inverts the perspective math (`src/viewer.js:434` `fov=800` unified) and, when `motionGain>1`, inverts the gain (`src/viewer.js:349` `motionGain` invert) so a click at amplified screen position maps to the true unamplified world point within `0.5 Å` (see `tests/test_picking.js` H73).
- **Clash snap (relaxation):** `placeLigand` centers the ligand on `target` (`src/placement.js:242` centroid + random SO(3) rotation), then `relaxClash` (`src/placement.js:149`) performs steepest descent over 6 rigid DOF (translation + Rodrigues rotation `rotExp`, `src/placement.js:50`) against the collision sigma set (`clashGrad`, `src/placement.js:76`) with `clashMargin=0.85` (`src/placement.js:152`), up to `350` iters, step `0.06` (`src/placement.js:150`). Result `converged` means `minRatio>=0.85` (sterically allowed), otherwise relaxed but still reported with `residualClash`.
- **Post-snap update:** `applyPlacedPose` writes the snapped pose into both live `state.integ.pos` and reference `state.ff.ref` (`src/ligand-panel.js:81`), disables holo springs (`src/ligand-panel.js:92`), rebuilds them, and calls `viewer.setSystem` to refresh `pocketCenter` and `center`/`radius` so the next render shows the ligand exactly where the ghost preview indicated.

## Cross-module policies (line refs)

- **RMSD split:** `HeavyForceField.rmsd` is protein-only (`src/heavy.js:1395`), `rmsdLig` covers the external-ligand block `[ligandStart, n)` (`src/heavy.js:1413`; hetero excluded — same slice the HUD `ligRMSD` uses), `rmsdAll` preserves the legacy all-atom value.
- **Live per-term mirror:** per-term energies stay live without BindLog via `liveTermsWanted` / `applyLiveTrackTerms` / `readLiveTerms` / `formatLiveTermsHUD` (`src/main.js:571`; mirror `updateLiveTermsMirror`, synthetic view `liveTermsBindLogView`); the wiring descriptor is `bindingTermsActive` (`src/ff-binding.js:51`).

## Files & Grep

- `grep -n "ghost" docs/PLACEMENT.md` hits this file (ghost preview).
- `grep -n "pocket" docs/PLACEMENT.md` hits pocket highlight.
- `grep -n "ligandStart" docs/PLACEMENT.md` hits the collision-set policy; `grep -n "excludeFrom" docs/PLACEMENT.md` hits slot exclusion.
- See `src/placement.js:238`, `src/viewer.js:382`, `src/ligand-panel.js:73` for source locations.

---
*See `src/placement.js:289` for pocket detection and `src/viewer.js:461` for pocket highlight rendering.*
