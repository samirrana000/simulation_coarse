/**
 * params.js — THE canonical home for the coarse-grained / heavy-atom
 * PARAMETER contract: per-element Lennard-Jones + partial charge, per-element
 * covalent radii, and the CG residue-class bead table.
 *
 * ── WHY THIS FILE EXISTS (goal M7, P1) ─────────────────────────────────────
 * A LJ sigma, a LJ epsilon and a partial charge are properties of a chemical
 * element. Declaring them in two places means the two places can disagree, and
 * when they do, the coarse-grained and heavy-atom engines silently compute two
 * different potentials for the same molecule. That is not a style nit; it is
 * the defect class this repo has already paid for twice:
 *
 *   • KB_KCAL was 0.001987204 in src/units.js and 0.0019872041 in
 *     src/ff-params.js — a 5.03e-8 split that made the reported temperature
 *     depend on which module an engine imported from.
 *   • The Coulomb constant was live at 332.0 in two force kernels (and inside
 *     a WGSL template string) while the main kernel used 332.06371 — 1.9e-4,
 *     so the worker and GPU kernels computed electrostatics 0.02% off the CPU.
 *
 * The third instance of the same class was hiding in plain sight:
 *
 *     src/ff-params.js:112  LIG_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0, hb: false, dG: -0.30 }
 *     src/heavy.js:82       HEAVY_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0 }
 *
 * Two declarations of the same quantity in two modules, one edit apart from a
 * live physics bug. (Measured: sigma/eps/q agreed to the last bit — 3.4/0.12/0.0
 * — so this was a maintenance hazard, not yet a divergence. The heavy copy was
 * also DEAD: nothing ever read it, the heavy path reads LIG_ELEMENT_DEFAULT.
 * See tests/test_element_params.js and tests/test_constant_ledger.js rule 5.)
 *
 * The fix follows the pattern this repo already proved with src/units.js: the
 * contract lives in ONE leaf module, every other module re-exports from it, and
 * a guard test fails if a second literal definition reappears. src/ff-params.js
 * is now that re-export facade — the public names (LIG_ELEMENT,
 * LIG_ELEMENT_DEFAULT, METAL_ELEMENT, COVALENT_RADIUS, RES_CLASS, …) are
 * unchanged, so no downstream import moves.
 *
 * ── WHY HERE, AND WHY A NEW FILE RATHER THAN src/units.js ─────────────────
 * 1. src/units.js is asserted to be the LEAF of the dependency graph: zero
 *    static imports (tests/test_constant_ledger.js rule 0). A module that every
 *    engine imports must be able to sit at or near that leaf, and this one is
 *    (it imports nothing at all), so it can.
 * 2. src/units.js holds UNIT CONVERSION FACTORS. A Lennard-Jones sigma is not
 *    a conversion factor — it is a force-field parameter with a literature
 *    provenance (CSD covalent radii, Miyazawa–Jernigan contact energies,
 *    AMBER parm99). Mixing the two would make "the single unit contract" mean
 *    two different things and would put 21 element charges next to k_B in the
 *    shadow-distance scan.
 * 3. src/physics/ is already where the parameter tables this one must not be
 *    confused with live: physics/forcefield/amber14sb.js (AMBER ff14SB
 *    bond/angle/nonbonded), physics/charges.js (AMBER partial charges,
 *    GB_RADII), physics/gb.js. Putting the shared element contract there makes
 *    "these are DIFFERENT models, do not merge" a directory fact.
 * 4. Name: `params.js`, not `element-params.js`, because it also owns the
 *    metal-ion and covalent-radius tables that both engines consume — all of
 *    it is "the parameter contract", and only the first block is per-element LJ.
 *
 * ── WHAT IS DELIBERATELY *NOT* MERGED HERE (two different models) ─────────
 * Keeping these apart is a physics decision, not an omission. Each is listed
 * with the test that holds the line.
 *
 *  (a) ELEMENT_LJ (per element, coarse, both engines)
 *      vs RES_CLASS (per CG residue-class bead, forcefield.js only).
 *      A Cα bead is not an atom: its sigma is a 3.8–4.1 Å effective diameter
 *      for a whole residue, not a van-der-Waals radius, and its weights come
 *      from Miyazawa–Jernigan contact energies, not from an element. Forcing
 *      them into one table would be a physics change wearing a refactor's
 *      clothes. tests/test_element_params.js asserts their key sets are
 *      DISJOINT, so a future merge cannot happen by accident.
 *
 *  (b) ELEMENT_LJ (coarse, per element)
 *      vs physics/forcefield/amber14sb.js NONBONDED_TABLE (AMBER parm99 σ
 *      converted from Rmin/2, per element, heavy mode only).
 *      Both are per-element LJ, and both are legitimately different models:
 *      the coarse table is a hand-rounded bead-scale approximation, the AMBER
 *      one is a published force field. E.g. for O the two give σ 3.0 vs 2.960
 *      (1.35% apart) but ε 0.16 vs 0.2100 (31% apart). tests/test_element_
 *      params.js measures that gap for all 9 shared elements on every run, and
 *      tests/test_constant_ledger.js rule 5 requires the second table to carry
 *      a written justification rather than existing silently.
 *
 *  (c) The partial CHARGE is per-element in ELEMENT_LJ but per-ATOM in the
 *      heavy path: heavy.js overwrites `q` with physics/charges.js
 *      assignCharges(), the approximate AMBER ff14SB mapping documented in
 *      docs/CHARGES.md. That is the heavy engine's defining feature (a united-
 *      atom ff14SB model), so it stays — but tests/test_element_params.js
 *      asserts both q columns separately and pins each, so neither can drift
 *      into the other unnoticed.
 *
 * Zero imports (like src/units.js) so it can sit anywhere in the graph.
 * Node-importable with no DOM globals.
 */

// ═══════════════════════════════════════════════════════════════════════════
// 1. PER-ELEMENT LJ + PARTIAL CHARGE  (σ Å, ε kcal/mol, q e, hb flag,
//    ΔG desolvation kcal/mol). The CG engine uses this for every ligand atom
//    (forcefield.js); the heavy engine uses it for EVERY non-metal atom —
//    protein backbone and side-chain heavy atoms included, which is why it is
//    named ELEMENT_LJ and not LIG_ELEMENT (the old name is kept as a re-export
//    alias from src/ff-params.js and is actively misleading).
// ═══════════════════════════════════════════════════════════════════════════
export const ELEMENT_LJ = {
  C:  { sigma: 3.4, eps: 0.12, q: 0.0,  hb: false, dG: -0.55 },
  N:  { sigma: 3.2, eps: 0.15, q: -0.30, hb: true,  dG: -0.35 },
  O:  { sigma: 3.0, eps: 0.16, q: -0.50, hb: true,  dG: -0.30 },
  S:  { sigma: 3.6, eps: 0.18, q: 0.0,  hb: false, dG: -0.45 },
  F:  { sigma: 2.9, eps: 0.10, q: -0.20, hb: true,  dG: -0.25 },
  CL: { sigma: 3.5, eps: 0.18, q: 0.0,  hb: false, dG: -0.40 },
  BR: { sigma: 3.6, eps: 0.20, q: 0.0,  hb: false, dG: -0.45 },
  I:  { sigma: 3.8, eps: 0.22, q: 0.0,  hb: false, dG: -0.50 },
  P:  { sigma: 3.5, eps: 0.14, q: 0.40, hb: false, dG: -0.35 },
};
/**
 * Fallback for any element absent from ELEMENT_LJ (B, SE, SI, AL, H, and
 * anything a PDB brings that nobody typed a row for).
 *
 * NOTE the two paths resolve this object IDENTICALLY — heavy.js reads this
 * same reference, not a projection of it. That was the whole point: it used to
 * also carry a private `HEAVY_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0 }`
 * in src/heavy.js, one edit away from disagreeing on the LJ well depth.
 *
 * PROVENANCE / OPEN PHYSICS QUESTION (not a refactor — reported, not changed):
 * this fallback is what the CG engine hands a METAL. ELEMENT_LJ has no metal
 * rows, so `resolveElementParams("ZN")` in CG mode returns these generic values
 * — σ 3.4 Å, ε 0.12 kcal/mol, q 0 — while the heavy engine hands the same atom
 * METAL_ELEMENT.ZN — σ 1.40 Å, ε 0.05 kcal/mol, q +2.0 e — plus coordination
 * springs. Measured gap: σ up to 2.62× too large (MG), ε uniformly 2.40× too
 * deep, and the ion's entire formal charge (Δq = −2.0 e for the divalents,
 * −1.0 e for Na⁺/K⁺) missing. It is REACHABLE: pdb.js parseLigands keeps any
 * HETATM group with ≥ 2 atoms, so a Zn coordinated inside a multi-atom
 * hetero group arrives in CG mode as a neutral, oversized, carbon-sized bead.
 * tests/test_element_params.js pins the current numbers so the gap cannot
 * silently widen, and prints the full table on every run. Closing it is a
 * physics decision (does the CG model want explicit ions? a metal-aware
 * desolvation term? coordination?), so M7 does not make it.
 */
export const ELEMENT_LJ_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0, hb: false, dG: -0.30 };

/**
 * THE resolution rule, in one place, so the two engines cannot resolve an
 * element by two different code paths ever again.
 *
 * Both engines call this for the non-metal case (forcefield.js for ligand
 * atoms, heavy.js via resolveHeavyElementParams for everything that is not a
 * metal). It is a pure function of the element symbol: `ELEMENT_LJ[el] ??
 * ELEMENT_LJ_DEFAULT`, i.e. exactly the expression both modules used before.
 *
 * @param {string} el uppercase element symbol
 * @returns {{sigma:number, eps:number, q:number, hb:boolean, dG:number}}
 *   the live table object (NOT a copy) — callers that mutate must copy first
 */
export function resolveElementParams(el) {
  return ELEMENT_LJ[el] ?? ELEMENT_LJ_DEFAULT;
}

// ═══════════════════════════════════════════════════════════════════════════
// 2. METAL IONS — a DIFFERENT model (explicit ions with coordination
//    geometry), consumed ONLY by the heavy engine. Kept out of ELEMENT_LJ on
//    purpose: see the note on ELEMENT_LJ_DEFAULT for the CG-side gap.
// ═══════════════════════════════════════════════════════════════════════════
/**
 * Metal-ion parameters for the all-atom heavy mode (heavy.js).
 *   sigma/eps — LJ size & well (Å, kcal/mol) for non-bonded repulsion,
 *   q         — formal charge (e), used for screened electrostatics,
 *   coordR    — metal–donor coordination distance (Å) used to build the
 *               coordination springs (heavy.js detects donors within coordR of
 *               the ion), and
 *   coordN    — target coordination number (how many donor springs to build;
 *               capped by however many donors are actually within coordR).
 * Metals are treated as explicit +2/+1 ions that coordinate N/O/S donors
 * (histidine N, carboxylate O, thiolate S, backbone carbonyl O) rather than
 * forming covalent bonds.
 */
export const METAL_ELEMENT = {
  ZN: { sigma: 1.40, eps: 0.05, q: 2.0, coordR: 2.30, coordN: 4 },
  FE: { sigma: 1.50, eps: 0.05, q: 2.0, coordR: 2.20, coordN: 6 },
  MG: { sigma: 1.30, eps: 0.05, q: 2.0, coordR: 2.15, coordN: 6 },
  CA: { sigma: 1.70, eps: 0.05, q: 2.0, coordR: 2.45, coordN: 6 },
  CU: { sigma: 1.40, eps: 0.05, q: 2.0, coordR: 2.20, coordN: 4 },
  MN: { sigma: 1.45, eps: 0.05, q: 2.0, coordR: 2.25, coordN: 6 },
  NI: { sigma: 1.40, eps: 0.05, q: 2.0, coordR: 2.15, coordN: 6 },
  CO: { sigma: 1.40, eps: 0.05, q: 2.0, coordR: 2.15, coordN: 6 },
  NA: { sigma: 1.70, eps: 0.05, q: 1.0, coordR: 2.50, coordN: 6 },
  K:  { sigma: 2.00, eps: 0.05, q: 1.0, coordR: 2.80, coordN: 6 },
};
export const METAL_ELEMENT_DEFAULT = { sigma: 1.50, eps: 0.05, q: 2.0, coordR: 2.30, coordN: 6 };

/**
 * THE heavy-engine element rule: metals first, then the shared coarse table.
 *
 * This is the ONE asymmetry between the engines, and it is now a named function
 * in the parameter contract instead of a bare ternary that lived 640 lines away
 * from the table it consults. It is exactly the expression heavy.js used before
 * (METAL_ELEMENT wins, otherwise resolveElementParams), so every number is
 * unchanged; what changed is that the rule is now visible, reviewable and
 * asserted in one place.
 *
 * The asymmetry itself is NOT removed here. Making CG metals match heavy metals
 * is a physics change (see the ELEMENT_LJ_DEFAULT note above): it needs a
 * decision about the CG model's ion model first.
 *
 * @param {string} el uppercase element symbol
 * @returns {object} metal row (with coordR/coordN) or the ELEMENT_LJ row
 */
export function resolveHeavyElementParams(el) {
  return METAL_ELEMENT[el] ?? resolveElementParams(el);
}

// ═══════════════════════════════════════════════════════════════════════════
// 3. COVALENT RADII (heavy-mode bond detection)
// ═══════════════════════════════════════════════════════════════════════════
/**
 * All-atom (heavy) covalent radii (Å) for the heavy mode bond-building
 * (heavy.js). A pair of heavy atoms within the SUM of their covalent radii
 * (× a 1.15 slack factor) is treated as covalently bonded. Solvent/water O
 * and H are handled by the parser (dropped), so these cover protein heavy
 * atoms + ligand heavy atoms.
 *
 * Source — Covalent radii from Cambridge Structural Database (CSD) surveys
 * (Allen et al., J. Chem. Soc. Perkin Trans. 2, 1987; Cordero et al., Dalton
 * Trans. 2008, 2832; Bondi, J. Phys. Chem. 1964, 68, 441 for van-der-Waals
 * reference). Values below are single-bond covalent radii matched to CSD
 * organic/metal-organic statistics:
 *
 * | Element | r_cov (Å) | CSD/Bondi source | Example bond | r_sum | r_sum×1.15 | Cap 2.2 |
 * |---------|-----------|------------------|--------------|-------|------------|---------|
 * | C       | 0.77      | CSD C(sp3) 0.76  | C–C 1.54     | 1.54  | 1.77       | pass |
 * | N       | 0.75      | CSD N 0.71–0.75  | C–N 1.47     | 1.52  | 1.75       | pass |
 * | O       | 0.73      | CSD O 0.66–0.73  | C–O 1.43     | 1.50  | 1.73       | pass |
 * | S       | 1.02      | CSD S 1.05       | S–S 2.04     | 2.04  | 2.35→2.20  | **pass (2.04 < 2.20 < 2.35)** |
 * | P       | 1.06      | CSD P 1.07       | P–O 1.60     | 1.79  | 2.06       | pass |
 * | CA      | 1.76      | CSD Ca 1.76      | Ca–N 2.51    | 2.51  | 2.89→2.20  | fail at 2.9 Å (spurious Ca–N rejected) |
 *
 * Rationale for thresholds — heavy.js:251 `r < BOND_SLACK*(rA+rB) && r < 2.2`:
 *   BOND_SLACK = 1.15 gives 15% tolerance for thermal elongation / PDB
 *   coordinate uncertainty while keeping C–C/C–N/C–O true bonds well inside.
 *   Hard cap 2.2 Å prevents spurious long-range contacts (e.g. Ca–N 2.9 Å,
 *   H-bond O⋯N 2.9 Å, van-der-Waals contacts) from being mis-typed as
 *   covalent despite large radii sums (Ca 1.76 + N 0.75 = 2.51 → 2.89 with
 *   slack). Critically, the biologically important disulfide S–S 2.04 Å
 *   (CSD mean 2.03–2.05 Å, crambin 1CRN SSBOND records 2.00/2.04/2.05 Å)
 *   has r_sum = 2.04 so r_sum×1.15 = 2.35; min(2.35, 2.20) = 2.20 still
 *   captures 2.04 Å with 0.16 Å margin, so all three crambin disulfides are
 *   recovered while Ca–N 2.9 Å is correctly excluded — see tests/test_topology.js.
 */
export const COVALENT_RADIUS = {
  C: 0.77, N: 0.75, O: 0.73, S: 1.02, P: 1.06, F: 0.71,
  CL: 0.99, BR: 1.14, I: 1.33, B: 0.84, SE: 1.20,
  ZN: 1.22, FE: 1.32, MG: 1.41, CA: 1.76, CU: 1.32, MN: 1.39,
  NI: 1.24, CO: 1.26, NA: 1.66, K: 2.03,
};
/**
 * Slack factor on the covalent-radius sum for bond detection.
 * 1.15 = CSD/Bondi covalent sum × 1.15, then hard cap 2.2 Å (see table above).
 * Validated: S–S 2.04 Å pass, Ca–N 2.9 Å fail — tests/test_topology.js.
 */
export const BOND_SLACK = 1.15;

// ═══════════════════════════════════════════════════════════════════════════
// 4. CG RESIDUE-CLASS BEADS — a DIFFERENT model from §1 (see the header's
//    "deliberately not merged" note (a)). Consumed only by forcefield.js.
// ═══════════════════════════════════════════════════════════════════════════
/** Residue class → protein-bead LJ parameters (σ Å, ε kcal/mol, charge e) */
export const RES_CLASS = {
  H:  { sigma: 4.0, eps: 0.15, q: 0 },  // hydrophobic
  A:  { sigma: 4.1, eps: 0.18, q: 0 },  // aromatic
  P:  { sigma: 3.8, eps: 0.12, q: 0 },  // polar (H-bond capable)
  Cp: { sigma: 3.6, eps: 0.10, q: 0 },  // positively charged (H-bond capable)
  Cn: { sigma: 3.6, eps: 0.10, q: 0 },  // negatively charged (H-bond capable)
};
export const RES_CLASS_OF = {
  ALA: "H", VAL: "H", LEU: "H", ILE: "H", PRO: "H", MET: "H", GLY: "H", CYS: "H",
  PHE: "A", TRP: "A", TYR: "A", HIS: "A",
  SER: "P", THR: "P", ASN: "P", GLN: "P",
  LYS: "Cp", ARG: "Cp",
  ASP: "Cn", GLU: "Cn",
};

/**
 * Formal bead charges (e) by residue identity — CG salt-bridge term
 * (R2 term b, docs/BINDING_PHYSICS_R2.md §2b; Loop-2 S1).
 *
 * Revives the dead screened-Coulomb path in ff-binding.js
 * (`E_coul = COULOMB_CONST·q_i·q_a/(ε(r)·r)·sw(r)`, ε(r) = 4+76·tanh(r/8)) by
 * giving the Cα bead of each charged residue its formal charge:
 *   ASP/GLU −1 (deprotonated carboxylate), LYS/ARG +1 (ammonium/guanidinium).
 *
 * Design decisions (R2 §2b + Loop-1 review S1):
 *   • HIS is deliberately absent ⇒ q = 0 (neutral default). A HIP +1 charge
 *     is assigned only when a protonation heuristic says so; the hookup to
 *     chem/protonation.js is explicitly deferred (documented, not guessed).
 *   • No mean-neutralization of net charge: the distance-dependent dielectric
 *     ε(r) together with the EEF1-lite burial/desolvation counterweight
 *     (Hendsch & Tidor 1994 lesson — bare Coulomb over-praises salt bridges,
 *     burial penalty rescues it) already temper net-charge effects; the
 *     formal ±1 values are kept as-is.
 *   • DEFAULT OFF: ForceField fills `_protQ` from this map only when
 *     `par.binding.charges === true` (opt-in flag). RES_CLASS.q stays 0, so
 *     the default path is bit-identical to the pre-S1 force field and the
 *     rollback is a single flag flip.
 */
export const CG_FORMAL_CHARGES = {
  ASP: -1, GLU: -1,
  LYS: 1, ARG: 1,
};

/**
 * Every element the project parameterises, and which table owns it. Exported
 * so a caller (or a test) can ask "is this element covered, and by what?" without
 * re-deriving the precedence rules. `owner` is one of "element" | "metal".
 * @returns {Record<string, {owner:"element"|"metal", hasRow:boolean}>}
 */
export function elementCoverage() {
  const out = {};
  for (const el of Object.keys(ELEMENT_LJ)) out[el] = { owner: "element", hasRow: true };
  for (const el of Object.keys(METAL_ELEMENT)) out[el] = { owner: "metal", hasRow: true };
  return out;
}