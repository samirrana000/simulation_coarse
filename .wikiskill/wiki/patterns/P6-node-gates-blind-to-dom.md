# Pattern P6 — Node-only gates cannot catch DOM-structure breakage

**Root cause:** ES-module syntax checks and JS regression suites validate
code, but HTML structure errors (orphaned closing tags, mis-nested panels)
silently produce broken layouts that still "pass" — the modules load, the
IDs exist, and node tests are green while the sidebar collapses to 24px.

**Evidence:** 2026-09-02 session — panel restructure left duplicate
`</div></details>` after the Chemical Network panel; all 62 DOM-contract IDs
still resolved (gate check #3 passed), 29/29 node tests passed, yet the
browser test could not click a preset link because the panels rendered
outside the collapsed `#controls` (Playwright: "subtree intercepts pointer
events"; debug: `#controls` height 24px, only 5 of 7 panels inside).

**Fix that works:** (1) Browser-level test (Playwright) is a REQUIRED gate
for any index.html change — not optional. (2) Debug recipe that found it in
one pass: `document.elementFromPoint()` on the click target + walk
`parentElement` chain printing bounding rects. (3) Structural sanity:
assert `#controls` scrollHeight > 500px and panel count inside `#controls`
equals expected (now part of visual verify script).

**Generalization:** The closer a test sits to the real runtime (real layout
engine), the more classes of bug it catches. For UI work, gate on the
browser test — node gates catch logic, browsers catch layout.
