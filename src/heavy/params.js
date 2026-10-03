/**
 * heavy/params.js — every number heavy mode owns locally, plus the cutoff switch.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10). src/heavy.js was 1674 LOC and the
 * largest module in the repo; this file is one of the modules that replaced it.
 * Responsibilities here, and nothing else:
 *
 *   - the non-bonded cutoff contract: R_CUT / R_SWITCH_ON / SCREEN_LEN and the
 *     C2 switch + its derivative (switchFunc / switchDeriv). Every non-bonded
 *     kernel (see heavy/nonbonded.js) switches through THESE two functions, so
 *     the smooth cut is defined in exactly one place.
 *   - K_ELEC, the public Coulomb constant, single-sourced from physics/gb.js.
 *   - METAL_K, the metal-coordination spring stiffness.
 *   - heavyMass, the per-element mass table.
 *
 * Per-element LJ / charge / Born-radius parameters are NOT here: they live in
 * the canonical table in src/physics/params.js and are read through
 * resolveHeavyElementParams (see heavy/forcefield.js).
 *
 * Zero DOM globals. Node-importable.
 */
import { COULOMB_CONST } from "../physics/gb.js";

/**
 * Coulomb constant, kcal·Å/(mol·e²). Public name kept (docs/UNITS.md lists
 * it), value now single-sourced: this used to read `332.0` — 1.9e-4 below the
 * CODATA value every other non-bonded kernel used, and it was dead (nothing
 * in the repo ever imported it), so unifying it costs no numeric result and
 * removes a trap for the next caller.
 */
export const K_ELEC = COULOMB_CONST;
export const SCREEN_LEN = 8.0; // Å
export const R_CUT = 8.5;      // Å
export const R_SWITCH_ON = 6.5;// Å

export function switchFunc(r) {
  if (r <= R_SWITCH_ON) return 1;
  if (r >= R_CUT) return 0;
  const rsq = r * r, rOn = R_SWITCH_ON, rCut = R_CUT;
  const on2 = rOn * rOn, cut2 = rCut * rCut;
  const denom = (cut2 - on2) ** 3;
  const num = (cut2 - rsq) * (cut2 - rsq) * (cut2 + 2 * rsq - 3 * on2);
  return num / denom;
}

export function switchDeriv(r) {
  if (r <= R_SWITCH_ON || r >= R_CUT) return 0;
  const rsq = r * r, rOn = R_SWITCH_ON, rCut = R_CUT;
  const on2 = rOn * rOn, cut2 = rCut * rCut;
  const denom = (cut2 - on2) ** 3;
  const a = (cut2 - rsq) * (cut2 - rsq);
  const b = cut2 + 2 * rsq - 3 * on2;
  const da = -4 * r * (cut2 - rsq);
  const db = 4 * r;
  return (da * b + a * db) / denom;
}

export const METAL_K = 40.0;
// M7: `HEAVY_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0 }` used to be
// declared HERE, a second table for a quantity that also lived in
// ff-params.js LIG_ELEMENT_DEFAULT. Measured: sigma/eps/q agreed to the last bit
// (3.4 / 0.12 / 0.0), so it was never a live divergence — but it was DEAD
// (nothing read it; the path below resolved through LIG_ELEMENT_DEFAULT) and
// one edit away from being a silent CG-vs-heavy LJ split, the same defect class
// as the KB_KCAL 5.03e-8 split and the 332.0-vs-332.06371 Coulomb split. The
// literal is gone; `resolveHeavyElementParams` (src/physics/params.js) is the
// one resolution rule both engines now call. Guard:
// tests/test_element_params.js + tests/test_constant_ledger.js rule 5.

/**
 * Physical atomic masses (Da) for heavy atoms — AMBER ff14SB / IUPAC.
 * Exhaustive over METAL_ELEMENT (ZN,FE,MG,CA,CU,MN,NI,CO,NA,K) and
 * LIG_ELEMENT (C,N,O,S,F,CL,BR,I,P) plus common extras B, SE, SI, AL.
 * Fallback for unknown elements warns via console.warn (not silent 14.0)
 * and returns 14.0 (approx N mass) so dynamics remain stable while the
 * caller is alerted — see tests for coverage.
 */
export function heavyMass(element) {
  switch (element) {
    case "C": return 12.011;
    case "N": return 14.007;
    case "O": return 15.999;
    case "S": return 32.06;
    case "P": return 30.974;
    case "F": return 18.998;
    case "CL": return 35.45;
    case "BR": return 79.904;
    case "I": return 126.904;
    case "B": return 10.81;
    case "SE": return 78.971;
    case "SI": return 28.085;
    case "AL": return 26.982;
    case "ZN": return 65.38;
    case "FE": return 55.845;
    case "MG": return 24.305;
    case "CA": return 40.078;
    case "CU": return 63.546;
    case "MN": return 54.938;
    case "NI": return 58.693;
    case "CO": return 58.933;
    case "NA": return 22.99;
    case "K": return 39.098;
    default:
      console.warn(`[heavyMass] unknown element "${element}" — fallback 14.0 Da (N mass); add to METAL_ELEMENT/LIG_ELEMENT if needed`);
      return 14.0;
  }
}
