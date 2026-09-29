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
 *   k_B        = 0.001987204 kcal/mol/K  (KB_KCAL)
 *   V°         = 1660.54 Å³  (≈ 1 M standard-state volume per molecule)
 *   C = e²/4πε₀ = 332.06371 kcal·Å/(mol·e²)  (COULOMB_CONST)
 *
 * EVERY physical constant in this file is a UNIT CONVERSION FACTOR, and each
 * one used to be re-typed as a literal somewhere else in src/ with a different
 * number of significant figures. Those are not style nits: k_B split by
 * 5.03e-8 across two modules made the reported temperature depend on which
 * file an engine imported, and the Coulomb constant was live at 332.0 in two
 * force kernels while the main one used 332.06371 (1.9e-4). This file is the
 * single home for all of them; every other module re-exports (see
 * tests/test_constant_ledger.js, which fails if a literal reappears).
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

/**
 * Coulomb constant C = e²/(4πε₀), in kcal·Å/(mol·e²).
 *
 * Derived from CODATA 2018 (e = 1.602176634e-19 C, exact by definition;
 * ε₀ = 8.8541878128e-12 F/m):
 *   e²/4πε₀ = 2.30707755e-28 J·m = 1.38935458e6 J·Å/mol
 *           / 4184 J/kcal      = 332.06371 kcal·Å/(mol·e²)
 * which is the value AMBER ff14SB / parm99 and CHARMM use. 332.0 is the
 * 4-significant-figure MD convention and is 1.9e-4 low — enough to make the
 * main-thread, worker and GPU non-bonded kernels disagree.
 *
 * Was declared with a literal in four files (physics/gb.js,
 * physics/solvation/gb_obc2.js, physics/solvation/membrane_slab.js,
 * compute/webgpu_backend.js), plus 332.0637 in ff-binding.js, 332.0 in
 * heavy.js, and a bare 332.0 inside the WGSL string in gpu.js and inside
 * force-worker.js. All of those now read this one.
 */
export const COULOMB_CONST = 332.06371;
