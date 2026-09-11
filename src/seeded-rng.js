/**
 * seeded-rng.js — Deterministic mulberry32 PRNG for reproducible simulations.
 *
 * Uses the mulberry32 algorithm (https://gist.github.com/tommyettinger/46a874533244883189143505d203312c):
 *   a fast 32-bit PRNG with period 2^32, suitable for deterministic Langevin
 *   noise and regression golden-file generation. The same integer seed always
 *   yields the same sequence of doubles in [0,1).
 *
 * Approach: the constructor seeds an internal 32-bit state. Each call to
 * rand() advances the state by 0x6D2B79F5 (a Weyl increment) then applies the
 * mulberry32 mixing:
 *     t = state += 0x6D2B79F5
 *     t = imul(t ^ t>>>15, t|1)
 *     t ^= t + imul(t ^ t>>>7, t|61)
 *     return ((t ^ t>>>14) >>>0) / 4294967296
 *
 * Determinism: Math.random is NOT patched by default; callers that need a
 * globally seeded run (e.g. tests/generate_golden.js) explicitly save and
 * restore Math.random via install()/restore() or by constructing a SeededRNG
 * and passing its rand() as the noise source. This avoids implicit global
 * side-effects in library code.
 *
 * Test hook: when executed directly (`node src/seeded-rng.js`) prints a
 * short self-test sequence so `node --check` and manual runs both succeed.
 */

export class SeededRNG {
  /**
   * @param {number} seed 32-bit integer seed (any integer; floats are truncated).
   */
  constructor(seed = 0) {
    this.seed = seed >>> 0;
    this.state = this.seed;
    this._originalRandom = null;
  }

  /**
   * Next uniform double in [0,1).
   * @returns {number}
   */
  rand() {
    let t = (this.state += 0x6D2B79F5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  /**
   * Uniform integer in [min, max) or [0, max) if only one arg.
   * @param {number} [minOrMax=0x7fffffff] lower bound or upper bound
   * @param {number} [max] upper bound (exclusive)
   * @returns {number}
   */
  randInt(minOrMax, max) {
    if (max === undefined) {
      if (minOrMax === undefined) return Math.floor(this.rand() * 0x100000000) >>> 0;
      const hi = minOrMax | 0;
      return Math.floor(this.rand() * hi);
    }
    const lo = minOrMax | 0;
    const hi = max | 0;
    if (hi <= lo) return lo;
    return lo + Math.floor(this.rand() * (hi - lo));
  }

  /**
   * Floating-point uniform in [min, max).
   * @param {number} [min=0]
   * @param {number} [max=1]
   */
  randFloat(min = 0, max = 1) {
    return min + this.rand() * (max - min);
  }

  /**
   * Gaussian N(0,1) via Box-Muller (consumes two uniform draws, caches spare).
   * @returns {number}
   */
  randn() {
    if (this._hasSpare) {
      this._hasSpare = false;
      return this._spare;
    }
    let u, v, s;
    do {
      u = this.rand() * 2 - 1;
      v = this.rand() * 2 - 1;
      s = u * u + v * v;
    } while (s >= 1 || s === 0);
    const mul = Math.sqrt((-2 * Math.log(s)) / s);
    this._spare = v * mul;
    this._hasSpare = true;
    return u * mul;
  }

  /** Install this RNG as global Math.random (save original for restore). */
  install() {
    if (this._originalRandom === null) this._originalRandom = Math.random;
    const self = this;
    Math.random = () => self.rand();
  }

  /** Restore original Math.random if install() was called. */
  restore() {
    if (this._originalRandom !== null) {
      Math.random = this._originalRandom;
      this._originalRandom = null;
    }
  }

  /** Reset to initial seed state. */
  reset() {
    this.state = this.seed;
    this._hasSpare = false;
    this._spare = 0;
  }
}

// Small test hook — runs when file is executed directly: `node src/seeded-rng.js`
// This is intentionally import-free so `node --check` (syntax only) still passes.
if (typeof process !== "undefined" && process.argv && process.argv[1] && process.argv[1].endsWith("seeded-rng.js")) {
  const rng = new SeededRNG(42);
  const seq = Array.from({ length: 5 }, () => rng.rand());
  // Expected first 5 values for seed 42 (mulberry32) — handy for regression
  // seq ≈ [0.601, 0.448, 0.852, 0.669, 0.174] — printed for manual inspection
  console.log("SeededRNG self-test seed=42:", seq.map((v) => v.toFixed(6)).join(", "));
  // Light assert
  const rng2 = new SeededRNG(42);
  const a = rng2.rand();
  const b = rng2.rand();
  if (Math.abs(a - 0.6011037519201636) > 1e-12) {
    console.error("SeededRNG self-test FAILED: first value mismatch", a);
    process.exit(1);
  }
  console.log("SeededRNG self-test PASSED");
}
