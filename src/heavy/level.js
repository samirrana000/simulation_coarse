/**
 * heavy/level.js — the heavy physics-level (L0/L1/L2) resolver and its
 * descriptor.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10).
 *
 * Same tier semantics as ForceField/resolvePhysicsLevel in src/forcefield.js:
 * L0 (default) = charges OFF / hbMode off; L1/L2 = charges ON / hbMode
 * directional. The heavy kernels already run full physics unconditionally
 * (assignCharges + GB + DirectionalHBond always on), so the resolver output is
 * stored queryably and never gates energy/force math. Explicit
 * par.binding.charges/hbMode win over the tier; unknown physicsLevel falls back
 * to L0. Pure + headless-safe; never throws.
 *
 * describePhysics lives here rather than on the class because it is this
 * resolver's own read-only projection — it reads only physicsLevel / chargesOn /
 * hbMode, which resolveHeavyPhysicsLevel produced.
 */
// ── Revolution 3 / Issue 1: heavy physics-level mirror (queryable only) ────
// Same tier semantics as ForceField/resolvePhysicsLevel (src/forcefield.js):
// L0 (default) = charges OFF / hbMode "off"; L1/L2 = charges ON /
// hbMode "directional". Heavy kernels already run full physics unconditionally
// (assignCharges + GB + DirectionalHBond always on), so the resolver output is
// stored queryably (physicsLevel/chargesOn/hbMode/describePhysics) and never
// gates energy/force math. Explicit par.binding.charges/hbMode win over the
// tier; unknown physicsLevel → L0. Pure + headless-safe; never throws.
export const DEFAULT_PHYSICS_LEVEL = "L0";
export const HEAVY_PHYSICS_LEVELS = {
  L0: { charges: false, hbMode: "off" },
  L1: { charges: true, hbMode: "directional" },
  // Heavy slice of UI L2 == L1 for charges/hbMode (L2 adds weakint via
  // par.weak + BindLog accumulators, consumed elsewhere).
  L2: { charges: true, hbMode: "directional" },
};
/**
 * Resolve the effective heavy binding flags from a HeavyForceField `par`.
 * @param {object} [par] constructor params ({physicsLevel, binding:{charges,hbMode}})
 * @returns {{level:string, charges:boolean, hbMode:string}}
 */
export function resolveHeavyPhysicsLevel(par = {}) {
  const raw = par?.physicsLevel;
  const level = (raw === "L1" || raw === "L2" || raw === "L0") ? raw : DEFAULT_PHYSICS_LEVEL;
  const tier = HEAVY_PHYSICS_LEVELS[level] ?? HEAVY_PHYSICS_LEVELS.L0;
  const b = par?.binding ?? {};
  const charges = (b.charges === undefined) ? tier.charges : (b.charges === true);
  const hbMode = b.hbMode ?? tier.hbMode;
  return { level, charges, hbMode };
}

/**
 * Read-only, headless-safe projection of the mirrored tier. Attached to
 * HeavyForceField.prototype by heavy/forcefield.js.
 *
 * Reports the mirrored tier flags, NOT a kernel gate: the heavy kernels
 * (charges + GB + directional H-bonds) run unconditionally, so
 * coulombActive/directionalHBActive are the tier flags rather than evidence
 * that a kernel is switched on.
 * @returns {{level:string, charges:boolean, hbMode:string,
 *   coulombActive:boolean, directionalHBActive:boolean,
 *   isSimplifiedDefault:boolean}}
 */
export const physicsLevelMethods = {
  describePhysics() {
    return {
      level: this.physicsLevel ?? DEFAULT_PHYSICS_LEVEL,
      charges: this.chargesOn === true,
      hbMode: this.hbMode,
      coulombActive: this.chargesOn === true,
      directionalHBActive: this.hbMode === "directional",
      isSimplifiedDefault:
        (this.physicsLevel ?? DEFAULT_PHYSICS_LEVEL) === "L0"
        && !(this.chargesOn === true)
        && this.hbMode === "off",
    };
  },
};
