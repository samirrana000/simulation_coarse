/**
 * heavy/pairs.js — the non-bonded exclusion key.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10). This existed as a private function
 * at the bottom of the monolith; it is shared by the exclusion bookkeeping in
 * heavy/forcefield.js, both grid kernels in heavy/nonbonded.js and the weak
 * pass in heavy/weak.js, so it is a real module edge rather than a copy.
 *
 * Capacity contract: the key is min*1e6 + max, collision-free for any atom
 * count below 1e6 (tests/test_exclusions.js pins the 1e6-factor argument). Do
 * not change the factor without re-pinning that test.
 */
export function pairKey(i, j) { return i < j ? i * 1e6 + j : j * 1e6 + i; }
