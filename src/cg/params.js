/**
 * cg/params.js — CG parameter resolution: constructor scalars and the
 * per-particle type tables. Split out of src/forcefield.js; moved verbatim.
 *
 * WHAT THIS FILE OWNS
 * -------------------
 * Turning a `par` object plus a bead list into the numbers the kernels read:
 *   - the scalar block (rc, gamma, kBond, kAngle, epsRep, sigmaRep, the
 *     binding/funnel flags, and the tier-resolved chargesOn/hbMode);
 *   - per-protein-bead residue class + LJ sigma/eps + formal charge + HB flag
 *     (RES_CLASS / CG_FORMAL_CHARGES);
 *   - the Cα-triplet virtual O-sites built when hbMode === "directional";
 *   - the per-ligand-atom element parameters (σ, ε, q, hb flag, ΔG), INCLUDING
 *     metals, which resolve from the shared METAL_ELEMENT ion table (a
 *     one-line user warning is emitted because CG builds no coordination
 *     restraint — see the block at the foot of assignParticleParams).
 *
 * WHY IT IS ITS OWN MODULE
 * ------------------------
 * This is the ONE place the CG engine decides what a residue or an element
 * IS. The energy kernels in ff-repulsion.js / ff-binding.js are handed
 * pre-built tables and never re-derive them, so every divergence risk
 * ("which σ does this bead get?") lives in a single readable file. It is
 * also the table that must stay in lockstep with src/heavy/params.js — both
 * call the SAME resolveElementParams from src/physics/params.js, and
 * tests/test_element_params.js reads THIS file to prove it.
 *
 * ZERO-ARITHMETM-ORDER NOTE: every loop below is a per-particle fill. No
 * cross-particle accumulation happens here, so the order in which the tables
 * are built cannot change any floating-point result.
 */

import { resolvePhysicsLevel } from "./level.js";
import {
  RES_CLASS, RES_CLASS_OF, CG_FORMAL_CHARGES, resolveElementParams,
  KBOND_DEFAULT, KANGLE_DEFAULT,
} from "../ff-params.js";
import { buildVirtualSites, coneAxisOf } from "../physics/virtual-sites.js";

/**
 * Constructor scalar block. Must run FIRST: the per-particle tables read
 * `chargesOn` / `hbMode`, and the topology builders read rc / gamma /
 * sigmaRep, so this cannot be reordered against the steps that follow.
 * @param {object} ff     the ForceField under construction
 * @param {object} par    constructor params
 */
export function initScalars(ff, par) {
  ff.nLigAtoms = 0;
  ff.rc = par.rc ?? 10.0;
  ff.gamma = par.gamma ?? 1.0;
  // Backbone stiffness from Boltzmann inversion (Tirion 1996, AMBER ff14SB Cα):
  //   k_b ≈ 100 kcal/mol/Å², k_θ ≈ 20 kcal/mol/rad² — see src/ff-params.js header
  //   and docs/CG_HEAVY.md §Backbone. Defaults equal KBOND_DEFAULT/KANGLE_DEFAULT
  //   (re-exported from ff-params.js) and match test_bond_dist.js: 200-step
  //   Langevin on 1crn keeps ⟨r⟩=3.81±0.05 Å.
  ff.kBond = par.kBond ?? KBOND_DEFAULT ?? 100.0;   // kcal/mol/Å²  (Cα–Cα peptide bond, Tirion 1996)
  ff.kAngle = par.kAngle ?? KANGLE_DEFAULT ?? 20.0;  // kcal/mol/rad² (backbone pseudo-angle, Tirion 1996)
  ff.epsRep = par.epsRep ?? 0.3;   // kcal/mol      (excluded-volume depth)
  ff.sigmaRep = par.sigmaRep ?? 4.0; // Å           (Cα bead diameter ≈ 2·2.0 Å)
  ff.bindOn = par.binding?.on ?? true; // protein–ligand binding potentials
  // Revolution 1 / Issue 3: tier-aware defaults (explicit, backward compatible).
  // par.physicsLevel "L0" (default) → charges OFF / hbMode "off" (legacy,
  // bit-identical). "L1"/"L2" → charges ON / "directional" (opt-in CG path).
  // Explicit par.binding.charges / par.binding.hbMode always win over the
  // tier (see resolvePhysicsLevel); ff.physicsLevel records the REQUESTED
  // tier so the simplification is queryable via describePhysics().
  // Loop-2 S2 (R2 §2a): H-bond mode. "off" (default) = legacy isotropic
  // bead-flag term, bit-identical to pre-S2. "directional" = Cα-triplet
  // virtual O-sites + angular gates. "all-flags" = legacy term with the
  // 52% class-A flag bug fixed via _protHB=1 everywhere (A/B stopgap).
  const _phys = resolvePhysicsLevel(par);
  ff.physicsLevel = _phys.level;
  ff.hbMode = _phys.hbMode;
  // CG salt-bridge charges (Loop-2 S1, R2 term b): when true, _protQ is
  // filled from CG_FORMAL_CHARGES (ASP/GLU −1, LYS/ARG +1; HIS 0 — neutral
  // default, HIP hookup deferred), reviving the screened-Coulomb path in
  // ff-binding.js. DEFAULT OFF: preserves the pre-S1 q=0 behavior exactly
  // (RES_CLASS.q stays 0), so existing trajectories are bit-identical
  // unless the caller opts in (par.physicsLevel "L1"/"L2" or
  // par.binding.charges === true). See src/ff-params.js:CG_FORMAL_CHARGES.
  ff.chargesOn = _phys.charges;
  ff.bindRcut = 9.0;                    // Å, cross-term cutoff
  ff.holoGamma = par.binding?.gammaLig ?? 0.5; // kcal/mol/Å²  (holo-pose spring)
  ff.holoOn = par.binding?.holo ?? true;
  ff.nHolo = 0;  // number of holo contact springs
  ff.funnel = null;                        // optional Funnel instance (set via setFunnel)
  ff.funnelOn = par.funnel?.on ?? false;   // bias enabled (default OFF — unbiased binding is the default)
}

/**
 * Per-protein-bead residue class + LJ params, the directional-HB virtual
 * sites, and the per-ligand-atom element tables. Runs after initScalars()
 * (needs chargesOn / hbMode) and before buildTopology() (needs _protQ for
 * nothing, but _ligSigma is read by the excluded-volume tables).
 * @param {object} ff     the ForceField under construction
 * @param {Array}  beads  pdb.selectSystem() beads
 */
export function assignParticleParams(ff, beads) {
  const nProt = ff.nProt;
  // Per-protein-bead residue class + LJ params
  ff.resClass = new Uint8Array(nProt);
  ff._protSigma = new Float64Array(nProt);
  ff._protEps = new Float64Array(nProt);
  ff._protQ = new Float64Array(nProt);
  ff._protHB = new Uint8Array(nProt);
  beads.forEach((b, i) => {
    const cls = RES_CLASS_OF[b.resName] ?? "H";
    ff.resClass[i] = ["H", "A", "P", "Cp", "Cn"].indexOf(cls);
    const p = RES_CLASS[cls];
    ff._protSigma[i] = p.sigma; ff._protEps[i] = p.eps; ff._protQ[i] = p.q;
    // Loop-2 S1 (R2 term b): opt-in formal charges override the class-level
    // q=0 with the residue-identity value (ASP/GLU −1, LYS/ARG +1, HIS 0).
    if (ff.chargesOn) ff._protQ[i] = CG_FORMAL_CHARGES[b.resName] ?? 0;
    ff._protHB[i] = (cls === "P" || cls === "Cp" || cls === "Cn") ? 1 : 0;
  });
  // Loop-2 S2: virtual interaction sites (backbone O acceptors for ALL
  // residues — fixes the 52% class-A gap where only P/Cp/Cn beads could
  // ever accept; hydrophobic backbones accept too in reality).
  ff._vSites = null;
  if (ff.hbMode === "directional") {
    ff._vSites = buildVirtualSites(beads);
    beads.forEach((b, i) => {
      const s = ff._vSites[i];
      if (s.valid) coneAxisOf(s, [b.x, b.y, b.z]);
    });
    ff._vSiteValence = new Uint8Array(nProt); // per-frame acceptor valence caps
  }
  // Per-ligand-atom element params
  ff._ligSigma = new Float64Array(ff.nLigAtoms);
  ff._ligEps = new Float64Array(ff.nLigAtoms);
  ff._ligQ = new Float64Array(ff.nLigAtoms);
  ff._ligHB = new Uint8Array(ff.nLigAtoms);
  ff._ligdG = new Float64Array(ff.nLigAtoms);
  const metalsHere = [];
  for (let a = 0; a < ff.nLigAtoms; a++) {
    // M7: `LIG_ELEMENT[el] ?? LIG_ELEMENT_DEFAULT` is now the shared
    // `resolveElementParams` from src/physics/params.js — the SAME function
    // heavy.js calls, so an element cannot resolve to a different sigma/eps
    // in the two engines. Identical values (ELEMENT_LJ is LIG_ELEMENT).
    // 2026-10: the resolver is METAL-aware (ELEMENT_LJ → METAL_ELEMENT →
    // ELEMENT_LJ_DEFAULT), so a metal ligand atom gets its ion row (sigma/eps/
    // q/coordR/coordN) on BOTH engines. This is a sound lookup in CG because
    // CG ligand atoms are explicit united-atom particles, not beads
    // (cg/system.js: one particle per ligand heavy atom) — only the protein is
    // coarsened to Cα. The `hb` and `dG` fields the CG binding pass reads
    // (_ligHB / _ligdG) are carried on the metal row too: hb=false (an ion
    // coordinates, it is not a donor/acceptor) and dG=-0.30, the EEF1-lite
    // placeholder. The placeholder is a declared LIMITATION, not a claim about
    // ion hydration — see the dG block in src/physics/params.js and
    // docs/LIMITATIONS.md. CG still builds NO metal-coordination restraint
    // (heavy/topology.js buildMetalCoordination is heavy-only) — hence the
    // one-line warning below.
    const e = resolveElementParams(ff.ligandAtoms[a].element);
    if (e.coordR !== undefined) metalsHere.push(ff.ligandAtoms[a].element);
    ff._ligSigma[a] = e.sigma; ff._ligEps[a] = e.eps; ff._ligQ[a] = e.q;
    ff._ligHB[a] = e.hb ? 1 : 0; ff._ligdG[a] = e.dG;
  }
  // User-facing honesty about what CG cannot do with a metal. Measured before
  // this warning existed (tools/exp_cg_metal_stability.mjs, 4 000 steps): the
  // ion is NOT ejected — a +2 charge here does not blow up, because the CG
  // Coulomb term is screened by ε(r)=4+76·tanh(r/8) and there is no
  // ion-atmosphere monopole to diverge. What CG cannot do is hold the ion's
  // COORDINATION SPHERE: with no metal–donor springs a Zn spends ~13 % of a
  // short run outside the physical 1.5–3.5 Å Zn–O band, and its formal charge
  // attracts it to the whole-residue Cα bead of a charged residue, so it can
  // sit ~2 Å from a bead instead of 2 Å from a donor ATOM. Silent metal-site
  // results are the failure mode this warns about; heavy mode is the path that
  // builds the coordination polyhedron. Once per ForceField, never per atom.
  if (metalsHere.length) {
    const kinds = [...new Set(metalsHere)].join(", ");
    console.warn(`[cg] metal ion(s) ${kinds} present: CG gives them real ion parameters ` +
      `(METAL_ELEMENT: σ/ε/formal charge) but builds NO metal–donor coordination ` +
      `restraint, so the coordination sphere is not held. Expect a ~2 Å metal–bead ` +
      `approach and occasional unphysical metal–donor distances. Use heavy mode for ` +
      `metal sites (see docs/LIMITATIONS.md).`);
  }
}