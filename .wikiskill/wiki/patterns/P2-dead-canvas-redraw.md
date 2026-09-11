# Pattern P2 — Event-driven-only redraw reads as "nothing renders"

**Root cause:** A canvas visualization updated only on state *change* is
indistinguishable from a dead one in the steady state, and shows no
empty-state or legend.

**Evidence:** src/network-panel.js:11 updateNetworkPlot() called (a) once at
init (src/network-panel.js:186-191) and (b) only inside
`if (curMacroState !== networkModel.currentState)` (src/main.js:652-659).
With no ligand placed or stable state, zero redraws after paint;
index.html:174 canvas has no empty-state text, no axes, no legend.

**Fix that works:** (1) redraw on a low cadence (e.g. 1 Hz) not only on
change; (2) explicit empty-state message when no ligand/system;
(3) one-line caption + legend so the canvas self-explains;
(4) explanatory text visible by default, not behind a nested <details>.

**Generalization:** Every live canvas needs three states drawn explicitly:
no-data, steady-state, and active-change. Missing any one reads as "broken"
to a user who didn't build it.
