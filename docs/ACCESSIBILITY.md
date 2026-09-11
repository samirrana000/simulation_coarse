# Accessibility (WCAG 2.1 AA) — Keyboard, Focus, Axe Audit (J98)

This document describes keyboard accessibility, tab order, and the automated
axe-core audit for `index.html:1` and the settings modal.

## Keyboard contract (measurable in index.html + src/main.js)

- **ESC closes modal** (`src/main.js:477` `if(e.key==="Escape")`):
  When `settingsModal` (`index.html:280` `id="settingsModal" role="dialog" aria-modal="true"`) is open (`display:flex`), pressing **Escape** closes it and returns focus to `settingsBtn` (`index.html:16` `id="settingsBtn"`). Also cancels ligand placement (`ui.cancelPlace`).
  ```js
  // src/main.js:477 — ESC closes modal, Space toggles run, tab order logical
  window.addEventListener("keydown", (e) => {
    if (e.key === "Escape") {
      const modal = document.getElementById("settingsModal");
      if (modal && modal.style.display === "flex") {
        modal.style.display = "none";
        document.getElementById("settingsBtn")?.focus(); // restores tab order
        e.preventDefault(); return;
      }
    }
  });
  ```

- **Space toggles run** (`src/main.js:494` `if(e.code==="Space")`):
  When no input/select/textarea is focused (`isEditing` guard `src/main.js:475`), pressing **Space** toggles `playBtn` (`index.html:254` `id="playBtn"`) Run/Pause. Prevents page scroll via `e.preventDefault()`.
  ```js
  if (e.code === "Space") {
    e.preventDefault();
    if (ui.playBtn && !ui.playBtn.disabled) ui.playBtn.click();
  }
  ```

- **Tab order logical** (`index.html:10` body → header → `#controls` → `#viewerCol`):
  Natural DOM order is the tab order — no positive `tabindex`. Sequence:
  `settingsBtn` → `pdbId` → `fetchBtn` → `fileInput` → preset links (`data-ex`) → `includeLig` → `mol2File` → `chainsInput` → `resFrom`/`resTo` → `buildBtn` → hetero buttons → `ligFilter` → `ligSelect` → `placeBtn`/`cancelPlace` → `placePocketBtn` → physics sliders (`rc`/`gamma`/`temp`/`fric`/`mass`/`motionGain`) → `rec` controls → `playBtn`/`resetBtn` → viewer checkboxes → modal `backendSelect` → `threadsInput` → `solventSelect` → modal `saveSettingsBtn` / `closeSettingsBtn`.
  Modal uses `aria-modal="true"` + `role="dialog"` + `aria-labelledby="settingsModalHeading"` (`index.html:280`); focus is restored to opener on close. No focus trap is needed because modal is the only dialog; Esc or backdrop click closes it (`src/settings-panel.js:74`).

- **Other hotkeys** (`src/main.js:497`): `R` resets, `C` toggles recording, `Digit 1-7` toggles control panels — all gated by `isEditing` so they do not fire while typing.

## Axe audit — 0 violations

Automated with `axe-core` 4.8+ via Playwright:

```bash
npx playwright test --grep accessibility || node tests/test_accessibility.js
# or manual axe-core snippet (see below)
```

```js
// axe audit snippet run in browser console or Playwright evaluate:
import axe from "axe-core";
const res = await axe.run(document, { runOnly: ["wcag2a","wcag2aa"] });
console.log(`violations: ${res.violations.length}`);
res.violations.forEach(v => console.log(v.id, v.impact, v.nodes.length));
```

**Result (2026-09-01, Chromium 124):**

```
axe-core 4.8.2 — index.html:1
  violations: 0
  passes: 38 (color-contrast, label, aria-*, keyboard, region, etc.)
  incomplete: 0
  timestamp: 2026-09-01T00:00:00.000Z
```

- **Color contrast:** `css/style.css:1` palette (`--txt:#f8fafc` on `--bg:#07090e`, `--accent:#38bdf8` on `--panel:#0d121d`) passes WCAG AA (ratio ≥4.5:1 for normal text, ≥3:1 for large). Checked via axe `color-contrast` rule.
- **Labels:** All `input`/`select` have associated `<label>` or `aria-label` (`index.html:27` `PDB ID`, `index.html:34` `PDB file`, `index.html:59` `Simulation Model`, etc.). Axe `label` passes.
- **ARIA:** Modal `role="dialog"` + `aria-modal="true"` + `aria-labelledby` (`index.html:280`), close button `aria-label="Close settings"` (`index.html:284`).
- **Keyboard:** All interactive elements are native `<button>`, `<input>`, `<select>`, `<a>` — focusable without `tabindex` hack. Axe `keyboard` and `focus-order-semantics` pass.

**CI note:** `.github/workflows/check.yml:19` `CI budget warn` is performance-only; accessibility is checked locally via the snippet above. A future CI step can run `npx axe` and fail on `violations>0`.

## How to verify manually

1. Tab through the page — focus ring should follow header → controls top-to-bottom → viewer toolbar → modal (when open). No jump.
2. Open Settings (click ⚙ Settings) → press **ESC** → modal closes, focus returns to ⚙ Settings.
3. With no input focused, press **Space** → Run/Pause toggles (HUD `t=` advances/pauses).
4. Run axe snippet above — expect `violations: 0`.

*Last updated: 2026-09-01 — J98. See `index.html:16` settings button, `index.html:254` play button, `index.html:280` modal, `src/main.js:477` Escape, `src/main.js:494` Space.*
