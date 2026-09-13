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
 *
 * Release history: 1.0.0-transform (2026-08-29, transformation era) →
 * 1.1.0-fp7 (2026-09-13, FP7 release closeout). Minor bump: Loop-2 + FP1–FP7
 * shipped backward-compatible additive features only (no breaking changes,
 * no parameter retuning, defaults unchanged); patch would understate seven
 * finish-product increments. No package.json in this repo — this constant
 * IS the version flow (HUD prefix + trajectory REMARK provenance).
 */

export const VERSION = "1.1.0-fp7";
export const BUILD_DATE = "2026-09-13";
export const BUILD_TIME = "2026-09-13T00:00:00Z";
