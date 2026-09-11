# Pattern P5 — UI-element contract fragility under refactor

**Root cause:** DOM IDs are a distributed contract: ui.js registry
(src/ui.js:10-33), main.js hotkey Digit1-7 panel indices (src/main.js:507-514),
and per-panel modules all bind by getElementById; any index.html restructure
can silently null a handle (guarded) or misindex a hotkey (unguarded).

**Evidence:** audit history — blank-render was caused by a syntax error, but
the same fragility class appears as `if (ui.x)` guards everywhere; hotkeys
1-7 assume panel order `#controls > .panel` (src/main.js:509-513).

**Fix that works:** After any index.html restructure: (1) run a headless
contract test that every ui.* id exists in the served HTML; (2) keep hotkey
indices in lockstep with panel count; (3) `node --check` all src.

**Generalization:** When refactoring a DOM-heavy app, write the contract
test BEFORE moving elements; IDs are a public API.
