# Phase 5 — Anti-Slop Cockpit

Per https://github.com/miqdadbadjuber/anti-slop and P1–P5 decluttering protocol.

- Rank, never remove: one visible flow Load → Build → Run → Record (2 `<details open>`); 5 collapsed; 0 numbered panels;
  Phase 1–4 controls preserved inside Advanced chemistry / Advanced sampling / ML Tier / Mutagenesis / Correlations / Kinetics.
- Top status bar (`#sysState #topPdb #topEngine #topStep`), Metrics HUD (T/Etot/PMF/DCCM), bottom dock (`#scrub #cvStrip #hudSpark` + hints).
- Three-state canvases: viewer/pmF/DCCM/network each render no-data text + ≤1 Hz steady + active caption + legend.
- Color class before element: ligand vivid palette + white halo + legend chip; protein C [180,180,180] vs ligand C [255,121,98] (asserted).
- Typography: ui-monospace scalars/HUD 10–12px; flat `#0B0D0E/#121517/#1B1F22` + 1px `#2A3036`; accents only
  #38BDF8/#F87171/#FBBF24/#34D399/#E879F9.
- Sliders: 6 slider+numeric pairs, `Math.fround` Float32, hot-rebuild. Hotkeys: Space run/pause, R reset, C record, M mutagenesis, 1–7 panels, Esc close.
- DOM contract: `node scripts/wikiskill_gate.js` OPEN (90 ui ids ⊆ index.html; Digit1-7 vs 7 panels). Full suite 29/29 files exit 0; `cli --pdb 4w52 --steps 10` E=12.13, RMSD 0.044 Å.
