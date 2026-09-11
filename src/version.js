/**
 * version.js — Deterministic build metadata for simulation_coarse.
 *
 * Licensed under the MIT License — see LICENSE at the repository root.
 *
 * VERSION is the canonical, human-readable release string used in all
 * provenance headers (trajectory files, HUD, logs). BUILD_DATE is the
 * ISO date of this transform (deterministic, not generated at runtime).
 * Both are imported by src/main.js (HUD), src/recorder.js (file headers)
 * and optionally src/units.js / src/viewer.js for display.
 */

export const VERSION = "1.0.0-transform";
export const BUILD_DATE = "2026-08-29";
export const BUILD_TIME = "2026-08-29T00:00:00Z";
