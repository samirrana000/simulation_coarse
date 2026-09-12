/**
 * settings-panel.js — System Settings modal dialog.
 *
 * Controls:
 *   - Compute Acceleration (WebGPU, Web Workers, CPU, thread count)
 *   - Biophysics & Solvent parameters (GB/SA, Debye-Hückel salt, dielectrics)
 *   - Chemical Network Model & Multi-scale dynamics
 *   - Integrator performance & sub-stepping
 */

import { ui, state } from "./ui.js?v=10";
import { GpuAccelerator } from "./gpu.js?v=10";
import { WorkerPool } from "./worker-pool.js?v=10";
import { initWebGPU, isSupported as isWebgpuSupported, webgpuStatus } from "./compute/webgpu_backend.js?v=10";

export const gpuAccelerator = new GpuAccelerator();
export const workerPool = new WorkerPool();

export const settingsState = {
  backend: "auto", // "auto" | "gpu" | "workers" | "cpu"
  numThreads: typeof navigator !== "undefined" && navigator.hardwareConcurrency ? Math.max(1, navigator.hardwareConcurrency - 1) : 2,
  solventModel: "gb", // "gb" | "dist_dep" | "eef1"
  saltM: 0.15, // 150 mM
  epsIn: 4.0,
  epsOut: 78.5,
  sasaGamma: 0.0072,
  subSteps: 10,
  chemicalNetworkOn: true,
  // Phase 3: WebGPU WGSL backend auto-detect (opt-in via backendSelect;
  // default path stays CPU/worker — see webgpu backend notes in Phase 4).
  webgpuSupported: isWebgpuSupported(),
  webgpuReady: false,
  webgpuNote: "probing…",
  // Phase 3: r-RESPA multiple-time-stepping (opt-in, default OFF).
  respaOn: false,
  respaOuterFs: 4, // outer non-bonded step in fs (inner bonded fixed 1 fs)
  // Loop-2 S7: physics fidelity level (opt-in, default "L0" = baseline).
  // In-memory persist for the Dynamics-panel selector; "L1" adds CG charges
  // + directional HB, "L2" adds heavy weakint + BindLog accumulators.
  physicsLevel: "L0",
};

/**
 * Stage-5: localStorage key for the physics-level selector.
 * No other setting persisted before (all in-memory); this follows the
 * specs/PLAN.md "persist to localStorage" convention. Default stays L0.
 */
export const PHYSICS_LEVEL_KEY = "sim.physicsLevel";
const _PHYSICS_LEVELS = ["L0", "L1", "L2"];

/** Guarded localStorage read → valid tier or null (fresh default L0). */
export function readStoredPhysicsLevel() {
  try {
    if (typeof localStorage === "undefined") return null;
    const v = localStorage.getItem(PHYSICS_LEVEL_KEY);
    return _PHYSICS_LEVELS.includes(v) ? v : null;
  } catch (_) { return null; }
}

/**
 * Guarded persist of the physics level (validates, never throws).
 * @param {string} lvl "L0"|"L1"|"L2"
 */
export function persistPhysicsLevel(lvl) {
  try {
    if (!_PHYSICS_LEVELS.includes(lvl)) return;
    settingsState.physicsLevel = lvl;
    if (typeof localStorage !== "undefined") localStorage.setItem(PHYSICS_LEVEL_KEY, lvl);
  } catch (_) { /* storage unavailable (private mode/headless) */ }
}

/** Restore the select element from stored state (guarded, default L0). */
export function restorePhysicsLevelSelect() {
  try {
    const stored = readStoredPhysicsLevel();
    if (stored) settingsState.physicsLevel = stored;
    const sel = (ui && ui.physicsLevel) ||
      (typeof document !== "undefined" ? document.getElementById("physicsLevel") : null);
    if (sel && _PHYSICS_LEVELS.includes(settingsState.physicsLevel)) sel.value = settingsState.physicsLevel;
  } catch (_) {}
}

// Stage-5: restore persisted tier at startup (default L0 when absent/invalid).
try { restorePhysicsLevelSelect(); } catch (_) {}

// Phase 3 auto-detect: probe the WGSL backend once at startup. Never throws;
// failure simply leaves the CPU/worker default in place.
try {
  initWebGPU().then((ok) => {
    settingsState.webgpuReady = !!ok;
    settingsState.webgpuNote = webgpuStatus();
    updateSettingsUI();
  }).catch((e) => {
    settingsState.webgpuReady = false;
    settingsState.webgpuNote = `probe failed: ${e?.message ?? e}`;
    try { updateSettingsUI(); } catch (_) {}
  });
} catch (e) {
  settingsState.webgpuReady = false;
  settingsState.webgpuNote = `probe failed: ${e?.message ?? e}`;
}

// Initialize GPU detection asynchronously
gpuAccelerator.init().then(() => {
  updateSettingsUI();
});

export function updateSettingsUI() {
  if (typeof document === "undefined") return; // headless/Node: no DOM to update
  const modal = document.getElementById("settingsModal");
  if (!modal) return;

  const gpuStatus = document.getElementById("gpuStatus");
  if (gpuStatus) {
    // Phase 3: prefer the WGSL backend status; fall back to legacy probe.
    const wgsl = settingsState.webgpuSupported
      ? `WGSL backend ${settingsState.webgpuReady ? "ready" : "not ready"} (${settingsState.webgpuNote})`
      : "WGSL backend n/a (no navigator.gpu)";
    const legacy = gpuAccelerator.isSupported
      ? `legacy GPU: ${gpuAccelerator.deviceInfo}`
      : `legacy GPU: ${gpuAccelerator.deviceInfo}`;
    gpuStatus.textContent = `${wgsl} · ${legacy}`;
    gpuStatus.style.color = (settingsState.webgpuReady || gpuAccelerator.isSupported) ? "#34d399" : "#94a3b8";
  }

  const cpuStatus = document.getElementById("cpuStatus");
  if (cpuStatus) {
    const cores = typeof navigator !== "undefined" ? (navigator.hardwareConcurrency || "unknown") : "unknown";
    cpuStatus.textContent = `${cores} logical cores detected (${workerPool.workers.length} active workers)`;
  }
}

export function initSettingsModal() {
  if (typeof document === "undefined") return; // headless/Node: no DOM to wire
  try { restorePhysicsLevelSelect(); } catch (_) {}
  const settingsBtn = document.getElementById("settingsBtn");
  const settingsModal = document.getElementById("settingsModal");
  const closeSettingsBtn = document.getElementById("closeSettingsBtn");
  const saveSettingsBtn = document.getElementById("saveSettingsBtn");

  if (settingsBtn && settingsModal) {
    settingsBtn.addEventListener("click", () => {
      updateSettingsUI();
      settingsModal.style.display = "flex";
    });
  }

  if (closeSettingsBtn && settingsModal) {
    closeSettingsBtn.addEventListener("click", () => {
      settingsModal.style.display = "none";
    });
  }

  if (settingsModal) {
    settingsModal.addEventListener("click", (e) => {
      if (e.target === settingsModal) {
        settingsModal.style.display = "none";
      }
    });
  }
  const threadsInput = document.getElementById("threadsInput");
  const threadsVal = document.getElementById("threadsVal");
  if (threadsInput && threadsVal) {
    threadsInput.addEventListener("input", () => {
      threadsVal.textContent = threadsInput.value;
    });
  }

  if (saveSettingsBtn) {
    saveSettingsBtn.addEventListener("click", () => {
      // Apply backend settings
      const backendSel = document.getElementById("backendSelect");
      if (backendSel) settingsState.backend = backendSel.value;

      const threadsInput = document.getElementById("threadsInput");
      if (threadsInput) {
        settingsState.numThreads = parseInt(threadsInput.value, 10);
        workerPool.setNumWorkers(settingsState.numThreads);
      }

      // Apply biophysics parameters
      const solventSel = document.getElementById("solventSelect");
      if (solventSel) settingsState.solventModel = solventSel.value;

      const saltInput = document.getElementById("saltInput");
      if (saltInput) settingsState.saltM = parseFloat(saltInput.value) / 1000.0;

      const epsInInput = document.getElementById("epsInInput");
      if (epsInInput) settingsState.epsIn = parseFloat(epsInInput.value);

      const epsOutInput = document.getElementById("epsOutInput");
      if (epsOutInput) settingsState.epsOut = parseFloat(epsOutInput.value);

      const sasaInput = document.getElementById("sasaInput");
      if (sasaInput) settingsState.sasaGamma = parseFloat(sasaInput.value);

      const netToggle = document.getElementById("netToggleModal");
      if (netToggle) settingsState.chemicalNetworkOn = netToggle.checked;

      // Phase 3: r-RESPA opt-in (default OFF). Read live from the Dynamics
      // panel checkbox when present so no save round-trip is required.
      try {
        const respaT = document.getElementById("respaToggle");
        if (respaT) settingsState.respaOn = !!respaT.checked;
        const respaO = document.getElementById("respaOuter");
        if (respaO) settingsState.respaOuterFs = Number(respaO.value) || 4;
      } catch (_) { /* respa stays default OFF */ }

      // Update active forcefield if present
      if (state.ff) {
        if (state.ff.gb) {
          state.ff.gb.epsIn = settingsState.epsIn;
          state.ff.gb.epsOut = settingsState.epsOut;
          state.ff.gb.saltM = settingsState.saltM;
          state.ff.gb.updateKappa();
        }
        if (state.ff.sasa) {
          state.ff.sasa.gamma = settingsState.sasaGamma;
        }
      }

      settingsModal.style.display = "none";
      if (ui.hud) {
        ui.hud.textContent = `Settings applied: Backend=${settingsState.backend}, Solvation=${settingsState.solventModel}, Salt=${(settingsState.saltM * 1000).toFixed(0)}mM`;
      }
    });
  }
}
