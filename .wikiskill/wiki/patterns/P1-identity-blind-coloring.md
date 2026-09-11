# Pattern P1 — Identity-blind color assignment

**Root cause:** Color semantics encode only element/chain; molecule class
(protein vs ligand) is never encoded, so chemically-similar atoms are
visually identical across classes.

**Evidence:** src/viewer.js:135-139 — heavy mode loops `i < this.n` and
assigns `ELEMENT_COLOR[el]` to every atom; ligand indices (i >= nProt) share
the same map as protein. CG mode partially avoids this (chain palette for
protein, element color for ligand, src/viewer.js:140-156) — proving the data
boundary was known but not applied in heavy mode.

**Fix that works (validated):** Keep protein CPK; shift ligand to a distinct
saturated palette + white halo ring (i >= nProt test already exists in the
render loop at src/viewer.js:552-561 as `isLig`). Measure: visual
differentiation is testable headlessly — assert `colors[i]` for ligand
indices differs from ELEMENT_COLOR of same element.

**Generalization:** Whenever two semantic classes share one visual channel,
the larger class wins and the minority becomes invisible. Encode class in
the channel BEFORE element (or via a second channel: halo/outline).
