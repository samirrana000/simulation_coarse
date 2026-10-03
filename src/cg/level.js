/**
 * cg/level.js — the CG physics-level (tier) resolver and its read-only
 * descriptor. Split out of src/forcefield.js; moved verbatim, no behaviour
 * change.
 *
 * WHAT THIS FILE OWNS
 * -------------------
 * Exactly one decision: which binding flags a ForceField `par` object means.
 * Everything downstream (charges, directional H-bonds, the ENM model label)
 * reads its answer from here and nothing else computes a tier.
 *
 * WHY IT IS ITS OWN MODULE
 * ------------------------
 * The tier answer is a *policy*, not physics: no arithmetic, no allocation,
 * no DOM, no dependency beyond nothing at all. It was previously the first
 * 40 lines of a 947-line file, so the one function a UI panel or a test
 * needed to reason about had to be found inside the force field that
 * implements it. Isolated here it is directly testable, and it is a leaf —
 * zero imports — so it can never introduce a cycle.
 */

// ── Revolution 1 / Issue 3: explicit CG physics-level default ─────────────
// Default CG L0 path is intentionally simplified (fast baseline) and therefore
// misses salt-bridge electrostatics + backbone directional H-bonds unless the
// caller opts in. This block makes that default EXPLICIT (not silent) without
// changing any energy kernel:
//
//   L0 (default, bit-identical legacy): charges OFF, hbMode "off", uniform ENM.
//     bindingU = cross LJ + isotropic bead-flag HB (P/Cp/Cn only) + desolv;
//     screened-Coulomb path stays dead (q=0), directional virtual-site HB off.
//   L1 / L2 (opt-in CG part of the UI tiers in src/main.js): charges ON +
//     hbMode "directional" (uniform ENM unless par.enmModel/par.seqWeight opt
//     in separately). Revives CG_FORMAL_CHARGES Coulomb + virtual-site HB.
//
// Wiring: `par.physicsLevel` ("L0"|"L1"|"L2", default "L0") selects the tier
// defaults below. Explicit `par.binding.charges` / `par.binding.hbMode` ALWAYS
// win over the tier, so existing callers passing raw flags are unaffected and
// the default path (no physicsLevel, no binding flags) is bit-identical to
// pre-change. Unknown levels fall back to "L0" (same as settings-panel).
export const DEFAULT_PHYSICS_LEVEL = "L0";
export const CG_PHYSICS_LEVELS = {
  L0: { charges: false, hbMode: "off", enm: "uniform" },
  L1: { charges: true, hbMode: "directional", enm: "uniform" },
  // CG slice of UI L2 == L1 (L2 adds heavy weakint + BindLog, consumed elsewhere).
  L2: { charges: true, hbMode: "directional", enm: "uniform" },
};
/**
 * Resolve the effective CG binding flags from a ForceField `par` object.
 * Pure + headless-safe; never throws; unknown physicsLevel → L0.
 * @param {object} [par] constructor params ({physicsLevel, binding:{charges,hbMode}})
 * @returns {{level:string, charges:boolean, hbMode:string}}
 */
export function resolvePhysicsLevel(par = {}) {
  const raw = par?.physicsLevel;
  const level = (raw === "L1" || raw === "L2" || raw === "L0") ? raw : DEFAULT_PHYSICS_LEVEL;
  const tier = CG_PHYSICS_LEVELS[level] ?? CG_PHYSICS_LEVELS.L0;
  const b = par?.binding ?? {};
  const charges = (b.charges === undefined) ? tier.charges : (b.charges === true);
  const hbMode = b.hbMode ?? tier.hbMode;
  return { level, charges, hbMode };
}

/**
 * Revolution 1 / Issue 3: explicit physics-fidelity descriptor (read-only,
 * headless-safe). Makes the L0 simplification queryable instead of silent:
 * L0 default reports charges:false / hbMode:"off" / uniform ENM with
 * isSimplifiedDefault:true; L1/L2 opt-in reports the revived Coulomb +
 * directional-HB wiring. No energy effect — pure introspection for UI/spec
 * wiring and regression tests. Lives beside the resolver because it is the
 * same question asked from the other side.
 * @param {object} ff  a ForceField instance
 * @returns {{level:string, charges:boolean, hbMode:string, enm:string,
 *   coulombActive:boolean, directionalHBActive:boolean,
 *   isSimplifiedDefault:boolean}}
 */
export function describePhysics(ff) {
  // Revolution 3 / Issue 2: directional live only with >=1 valid vSite
  // (zero-valid ⇒ no dir; presence of the _vSites array alone is not enough).
  let nValidDir = 0;
  const vsDir = ff._vSites;
  if (vsDir) {
    const nDir = Number.isFinite(ff.nProt) ? ff.nProt : vsDir.length;
    for (let i = 0; i < nDir; i++) if (vsDir[i] && vsDir[i].valid) nValidDir++;
  }
  const directional = ff.hbMode === "directional" && nValidDir >= 1;
  return {
    level: ff.physicsLevel ?? DEFAULT_PHYSICS_LEVEL,
    charges: ff.chargesOn === true,
    hbMode: ff.hbMode,
    enm: ff.enmModel ?? "uniform",
    coulombActive: ff.chargesOn === true,
    directionalHBActive: directional,
    isSimplifiedDefault:
      (ff.physicsLevel ?? DEFAULT_PHYSICS_LEVEL) === "L0"
      && !(ff.chargesOn === true)
      && ff.hbMode === "off",
  };
}