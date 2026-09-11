/**
 * units.js — Single unit contract for simulation_coarse.
 *
 * Canonical base units (all physics, integrator, and I/O use these):
 *   Length : Ångström (Å)
 *   Time   : picosecond (ps)
 *   Energy : kcal/mol
 *   Mass   : Dalton (Da ≈ g/mol)
 *   Temp   : Kelvin (K)
 *
 * Conversion notes:
 *   1 kcal/mol = 418.4 Da·Å²/ps²  (mechanical units)
 *   a[Å/ps²]   = KCONV · F[kcal/mol/Å] / m[Da]
 *   k_B        = 0.0019872041 kcal/mol/K  (KB_KCAL)
 *   V°         = 1660.54 Å³  (≈ 1 M standard-state volume per molecule)
 *
 * No imports (except version.js which is also leaf) — this is the leaf of the
 * dependency graph so every other module can safely import from here without
 * creating a cycle. Version provenance: see src/version.js (VERSION, BUILD_DATE).
 */

/** kcal/mol → Da·Å²/ps²  (also exported as KCONV for backward compat). */
export const KCAL_TO_DA_A2_PS2 = 418.4;

/** Boltzmann constant in kcal/mol/K. */
export const KB_KCAL = 0.001987204;

/** Force → acceleration conversion: KCONV = 418.4 Da·Å²/ps² per kcal/mol. */
export const KCONV = 418.4;

/** Standard-state volume for 1 M concentration: 1660.54 Å³. */
export const STANDARD_VOLUME = 1660.54;
