/**
 * controllers/transport.js — Run/Pause, Reset, viewer display toggles and the
 * global scientific hotkeys.
 *
 * PUBLIC API
 *   initTransport()  Wire #playBtn, #resetBtn, the five viewer display
 *                    checkboxes, and the window keydown handler.
 *
 * NEEDS (imports): ui + state + viewer + recorder from ../ui.js; updateGuide
 * from ./guide.js.
 *
 * Hotkey contract (docs/ACCESSIBILITY.md, FP3): every binding is gated on
 * `isEditing` (INPUT/SELECT/TEXTAREA focused) and Space additionally on
 * `activatesNatively` (BUTTON/A/SUMMARY focused), so a keyboard user tabbed
 * onto a natively-activatable control never has Space stolen from them.
 * Digit1-7 toggles `#controls > .panel` by INDEX, so the panel order in
 * index.html and the regex range below are a contract — see
 * scripts/wikiskill_gate.js check 4 and docs/ACCESSIBILITY.md.
 */

import { ui, state, viewer, recorder } from "../ui.js";
import { updateGuide } from "./guide.js";
import { ignore } from "../errors.js";

/** Wire the transport controls, the viewer display toggles and the hotkeys. */
export function initTransport() {
  if (ui.playBtn) {
    ui.playBtn.addEventListener("click", () => {
      if (!state.integ) return;
      state.running = !state.running;
      ui.playBtn.textContent = state.running ? "⏸ Pause" : "▶ Run";
      try { updateGuide(); } catch (e) { ignore(e, "updateGuide@Run toggle", "checklist DOM absent headless; the run state itself is already latched"); } // FP1: Run ✓ (latches via steps/time)
    });
  }

  if (ui.resetBtn) {
    ui.resetBtn.addEventListener("click", () => {
      if (!state.integ) return;
      state.integ.reset();
      state.nanWarning = false;
    });
  }

  if (ui.showContacts) {
    ui.showContacts.addEventListener("change", () => { if (viewer) viewer.showContacts = ui.showContacts.checked; else console.warn("[viewer] not ready"); });
  }
  if (ui.spheres) {
    ui.spheres.addEventListener("change", () => { if (viewer) viewer.drawSpheres = ui.spheres.checked; else console.warn("[viewer] not ready"); });
  }
  if (ui.showRibbon) {
    ui.showRibbon.addEventListener("change", () => { if (viewer) viewer.drawRibbon = ui.showRibbon.checked; else console.warn("[viewer] not ready"); });
  }
  if (ui.showHBonds) {
    ui.showHBonds.addEventListener("change", () => { if (viewer) viewer.showHBonds = ui.showHBonds.checked; else console.warn("[viewer] not ready"); });
  }
  if (ui.showStates) {
    ui.showStates.addEventListener("change", () => { if (viewer) viewer.showStates = ui.showStates.checked; else console.warn("[viewer] not ready"); });
  }

  /* ------------------------------------------------------------------ */
  /*  Global Scientific Hotkeys (Space, R, M, C, 1-7, Escape)            */
  /* ------------------------------------------------------------------ */
  if (typeof window !== "undefined") {
    window.addEventListener("keydown", (e) => {
      const tag = document.activeElement ? document.activeElement.tagName : "";
      const isEditing = tag === "INPUT" || tag === "SELECT" || tag === "TEXTAREA";

      if (e.key === "Escape") {
        const modal = document.getElementById("settingsModal");
        if (modal && modal.style.display === "flex") {
          modal.style.display = "none";
          document.getElementById("settingsBtn")?.focus();
          e.preventDefault();
          return;
        }
        if (ui.cancelPlace && !ui.cancelPlace.disabled) {
          ui.cancelPlace.click();
          e.preventDefault();
          return;
        }
      }

      if (isEditing) return;

      // FP3 a11y: Space must not hijack natively-activatable elements. When a
      // <button>/<a>/<summary> has focus, Space already does the right thing
      // (activate the button, follow/toggle natively) — stealing it for Run
      // would break keyboard users tabbed onto any other control. (Space on a
      // focused Run button still toggles via its native click.)
      const activatesNatively = tag === "BUTTON" || tag === "A" || tag === "SUMMARY";
      if (e.code === "Space") {
        if (activatesNatively) return;
        e.preventDefault();
        if (ui.playBtn && !ui.playBtn.disabled) ui.playBtn.click();
      } else if (e.code === "KeyR" && !e.ctrlKey && !e.metaKey) {
        e.preventDefault();
        if (ui.resetBtn && !ui.resetBtn.disabled) ui.resetBtn.click();
      } else if (e.code === "KeyM" && !e.ctrlKey && !e.metaKey) {
        // [M] mutagenesis panel: open PMF & Analysis + its Mutagenesis disclosure
        e.preventDefault();
        const panels = document.querySelectorAll("#controls > .panel");
        const pmf = [...panels].find((p) => p.querySelector("#alaScanBtn")) || panels[6];
        if (pmf) {
          pmf.open = true;
          pmf.scrollIntoView({ behavior: "auto", block: "nearest" });
          const sub = pmf.querySelector(".subpanel summary");
          if (sub && ui.alaScanBtn) {
            const det = ui.alaScanBtn.closest("details");
            if (det) det.open = true;
          }
          ui.alaRes?.focus?.();
        }
      } else if (e.code === "KeyC" && !e.ctrlKey && !e.metaKey) {
        e.preventDefault();
        if (recorder.recording) {
          ui.recStopBtn?.click();
        } else {
          ui.recBtn?.click();
        }
      } else if (/^Digit[1-7]$/.test(e.code)) {
        // FP3: index from e.code (layout-independent) — e.key yields symbols
        // with Shift held (e.g. "!" for Digit1) and would misindex.
        const idx = Number(e.code.slice(5)) - 1;
        const panels = document.querySelectorAll("#controls > .panel");
        if (panels[idx]) {
          panels[idx].open = !panels[idx].open;
          panels[idx].scrollIntoView({ behavior: "auto", block: "nearest" });
        }
      }
    });
  }
}