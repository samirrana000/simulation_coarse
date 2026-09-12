/**
 * forcefield.js — Cα coarse-grained potential with anisotropic-network-style
 * native contacts and implicit-solvent excluded volume.
 *
 * Energy model (all standard CG/ENM physics; units Å, ps, kcal/mol, Da):
 *
 *   1. PEPTIDE BONDS (backbone connectivity, per contiguous segment):
 *        U_bond = Σ ½ k_b (r − r0)² ,       k_b ≈ 100 kcal/mol/Å²
 *      with r0 = equilibrium Cα–Cα distance (≈ 3.81 Å) from the input structure.
 *
 *   2. BACKBONE ANGLES (pseudo-angle i−1, i, i+1 within a segment):
 *        U_ang = Σ ½ k_θ (θ − θ0)² ,        k_θ ≈ 20 kcal/mol/rad²
 *      θ0 taken from the input coordinates (native-like local stiffness).
 *
 *   3. NATIVE CONTACTS — Elastic Network Model (K Tirion-style / ANM-uniform
 *      spring): every bead pair whose NATIVE distance r0 ≤ Rc and that is not
 *      already a 1-2 or 1-3 sequence neighbour gets a Hookean spring:
 *        U_ENM = Σ ½ γ (r − r0)² ,          H(Rc − r0) implicit by construction
 *      This preserves the native fold while allowing thermal fluctuations that
 *      reproduce experimental B-factor statistics at T = 300 K for γ ≈ 1.
 *
 *   4. NON-NATIVE EXCLUDED VOLUME (implicit-solvent, repulsive 12-10-free
 *      Lennard-Jones form, attractive tail switched off so beads never
 *      collapse spuriously):
 *        U_rep = ε [ (r_e/r)¹² − 2 (r_e/r)⁶ + 1 ]   for r < r_e,
 *              = 0                                      otherwise,
 *      with  r_e = 2^(1/6) σ ,  σ = 4.0 Å (Cα bead diameter), ε = 0.3 kcal/mol.
 *      This term is smooth at the cutoff and prevents chain crossing.
 *
 *   5. PROTEIN–LIGAND BINDING (pair pass): protein beads 0..nProt−1 vs. ligand
 *      atoms nProt..n−1 are excluded from the grid repulsion above and handled
 *      by a dedicated binding term — cross 12-6 LJ (attractive well, Lorentz–
 *      Berthelot combining rules), screened electrostatics (εr = 4+76·tanh(r/8))
 *      and an H-bond term, each smoothstep-switched off at bindRcut ≈ 9 Å.
 *
 *   5b. HOLO CONTACT SPRINGS (native protein–ligand pose): every protein
 *      bead–ligand atom pair whose NATIVE distance r0 ≤ 6 Å gets a Hookean
 *      spring U = ½ γ_lig (r − r0)² with γ_lig = 0.5 kcal/mol/Å², encoding the
 *      observed bound pose of a holo complex (e.g. T4 lysozyme L99A + benzene,
 *      PDB 4W52). The pair is excluded from the binding pair pass above so the
 *      spring alone governs the native-contact distance. Gated by
 *      par.binding.holo (default true); no-op without ligands.
 *
 *   6. EEF1-LITE DESOLVATION/BURIAL (apolar binding): while the same 27-cell
 *      scan drives the pair terms above, each ligand atom also accumulates a
 *      "soft protein-occupancy count" n_a = Σ_j g(r_aj) over the protein beads
 *      it sees, with g(r) = exp(−(r−r0)²/(2σ²)), r0 = 4.5 Å, σ = 1.8 Å. g has
 *      decayed to ≈ 0.04 by r ≈ 9 Å = bindRcut, so the existing grid pass fully
 *      covers the density field — no second grid, no second cutoff. The burial
 *      fraction B_a = 1 − exp(−n_a/3) ∈ [0,1] (0 exposed, 1 buried) converts the
 *      per-atom hydrophobic transfer energy ΔG_a (≈ −0.55 kcal/mol for C) into
 *      a burial cost U_desolv = Σ_a ΔG_a·B_a (negative ⇒ burial favorable).
 *      Forces follow the chain rule dU/dr = ΔG_a·(dB/dn)·(dg/dr) and are
 *      applied in a second pass over the recorded (a, j, r) pairs (see
 *      `_binding`).
 *
 *   7. OPTIONAL FUNNEL BIAS (off by default; par.funnel.on): external
 *      collective-variable bias on the ligand↔pocket COM distance — a binding
 *      funnel that is flat inside the bound state plus well-tempered
 *      metadynamics for PMF reconstruction (see src/funnel.js). Attached by
 *      the caller via setFunnel(); contributes nothing while funnelOn is false.
 *
 * The system is a unified multi-particle set: the Cα beads (global indices
 * 0..nProt−1) followed by any ligand heavy atoms (indices nProt..n−1). When no
 * ligands are given the protein-only model is reproduced exactly. Ligands get
 * a united-atom internal force field (bonds/angles/impropers from ligand.js)
 * plus the same excluded-volume grid; covalently linked ligand pairs are
 * excluded from the repulsion so the rigid molecule never repels itself.
 *
 * Performance: bonds/angles/ENM springs are O(N); the repulsive term uses a
 * uniform grid (cell list) in 3D rebuilt each force call — O(N) average.
 * All arithmetic is done on flat Float64Array buffers for speed; forces are
 * accumulated in a preallocated flat array — zero allocations per step.
 */

import { buildLigandInternalFF, improperAngle } from "./ligand.js?v=10";
import {
  harmonicPairs, springForces, angleForces, ligandBondForces, improperForces,
} from "./ff-harmonic.js?v=10";
import { repulsion } from "./ff-repulsion.js?v=10";
import { binding } from "./ff-binding.js?v=10";
import { KB_KCAL, KCONV } from "./units.js?v=10";
import {
  RES_CLASS, RES_CLASS_OF, CG_FORMAL_CHARGES, LIG_ELEMENT, LIG_ELEMENT_DEFAULT,
  SEQ_WEIGHT, KBOND_DEFAULT, KANGLE_DEFAULT,
} from "./ff-params.js?v=10";
import { buildTirionNetwork, applyTirionToForceField } from "./physics/forcefield/tirion_anm.js?v=10";
import { buildVirtualSites, coneAxisOf } from "./physics/virtual-sites.js";

// canonical units live in units.js; re-export keeps backward compat for
// integrator.js / funnel.js / tests that historically imported from here.
export { KB_KCAL, KCONV };

export class ForceField {
  /**
   * @param {object} sel  output of pdb.selectSystem(): {beads, segments}
   * @param {object} par  { rc, gamma }  — ENM cutoff (Å) and spring constant
   * @param {Array}  ligands  output of pdb.parseLigands() (empty ⇒ protein-only)
   */
  constructor(sel, par = {}, ligands = []) {
    const { beads, segments } = sel;
    const nProt = this.nProt = beads.length;
    this.nLigAtoms = 0;
    this.rc = par.rc ?? 10.0;
    this.gamma = par.gamma ?? 1.0;
    // Backbone stiffness from Boltzmann inversion (Tirion 1996, AMBER ff14SB Cα):
    //   k_b ≈ 100 kcal/mol/Å², k_θ ≈ 20 kcal/mol/rad² — see src/ff-params.js header
    //   and docs/CG_HEAVY.md §Backbone. Defaults equal KBOND_DEFAULT/KANGLE_DEFAULT
    //   (re-exported from ff-params.js) and match test_bond_dist.js: 200-step
    //   Langevin on 1crn keeps ⟨r⟩=3.81±0.05 Å.
    this.kBond = par.kBond ?? KBOND_DEFAULT ?? 100.0;   // kcal/mol/Å²  (Cα–Cα peptide bond, Tirion 1996)
    this.kAngle = par.kAngle ?? KANGLE_DEFAULT ?? 20.0;  // kcal/mol/rad² (backbone pseudo-angle, Tirion 1996)
    this.epsRep = par.epsRep ?? 0.3;   // kcal/mol      (excluded-volume depth)
    this.sigmaRep = par.sigmaRep ?? 4.0; // Å           (Cα bead diameter ≈ 2·2.0 Å)
    this.bindOn = par.binding?.on ?? true; // protein–ligand binding potentials
    // Loop-2 S2 (R2 §2a): H-bond mode. "off" (default) = legacy isotropic
    // bead-flag term, bit-identical to pre-S2. "directional" = Cα-triplet
    // virtual O-sites + angular gates. "all-flags" = legacy term with the
    // 52% class-A flag bug fixed via _protHB=1 everywhere (A/B stopgap).
    this.hbMode = par.binding?.hbMode ?? "off";
    // CG salt-bridge charges (Loop-2 S1, R2 term b): when true, _protQ is
    // filled from CG_FORMAL_CHARGES (ASP/GLU −1, LYS/ARG +1; HIS 0 — neutral
    // default, HIP hookup deferred), reviving the screened-Coulomb path in
    // ff-binding.js. DEFAULT OFF: preserves the pre-S1 q=0 behavior exactly
    // (RES_CLASS.q stays 0), so existing trajectories are bit-identical
    // unless the caller opts in. See src/ff-params.js:CG_FORMAL_CHARGES.
    this.chargesOn = par.binding?.charges === true;
    this.bindRcut = 9.0;                    // Å, cross-term cutoff
    this.holoGamma = par.binding?.gammaLig ?? 0.5; // kcal/mol/Å²  (holo-pose spring)
    this.holoOn = par.binding?.holo ?? true;
    this.nHolo = 0;  // number of holo contact springs
    this.funnel = null;                        // optional Funnel instance (set via setFunnel)
    this.funnelOn = par.funnel?.on ?? false;   // bias enabled (default OFF — unbiased binding is the default)

    // Reference coordinates (native state): copied, never mutated.
    this.ref = new Float64Array(nProt * 3);
    beads.forEach((b, i) => {
      this.ref[3 * i] = b.x; this.ref[3 * i + 1] = b.y; this.ref[3 * i + 2] = b.z;
    });

    // ---- united-atom ligand force field (appended particles) --------------
    // Ligand heavy atoms become explicit particles with global indices
    // nProt..n−1. buildLigandInternalFF already emits ABSOLUTE indices and
    // united-atom masses, so its flat arrays are consumed directly. isLigand
    // flags the appended particles for the analysis helpers.
    this.isLigand = new Uint8Array(nProt);
    this.ligandAtoms = [];
    this.ligandMasses = new Float64Array(0);
    this.ligandBonds = new Float64Array(0);
    this.ligandAngles = new Float64Array(0);
    this.ligandImpropers = new Float64Array(0);
    if (ligands.length) {
      const lig = buildLigandInternalFF(ligands, nProt);
      this.ligandAtoms = lig.atoms;
      this.ligandMasses = lig.masses;
      this.ligandBonds = lig.bonds;
      this.ligandAngles = lig.angles;
      this.ligandImpropers = lig.impropers;

      // Append ligand coordinates in molecule/atom concatenation order — the
      // atom table from ligand.js has no x/y/z, so read them off the input.
      const ref = new Float64Array((nProt + lig.nLigAtoms) * 3);
      ref.set(this.ref.subarray(0, nProt * 3));
      let a = 0;
      for (const mol of ligands) {
        for (const at of mol.atoms) {
          ref[3 * (nProt + a)] = at.x;
          ref[3 * (nProt + a) + 1] = at.y;
          ref[3 * (nProt + a) + 2] = at.z;
          a++;
        }
      }
      this.ref = ref;
      this.n = nProt + lig.nLigAtoms;
      this.nLigAtoms = lig.nLigAtoms;
      this.isLigand = new Uint8Array(this.n);
      this.isLigand.fill(1, nProt);
    } else {
      this.n = nProt;
    }

    // Masses: Cα beads default to 110 Da (united ALA), ligand atoms to their
    // united-atom mass (H implicitly folded in).
    this.masses = new Float64Array(this.n).fill(par.mass ?? 110);
    if (ligands.length) this.masses.set(this.ligandMasses, nProt);

    // Per-protein-bead residue class + LJ params
    this.resClass = new Uint8Array(nProt);
    this._protSigma = new Float64Array(nProt);
    this._protEps = new Float64Array(nProt);
    this._protQ = new Float64Array(nProt);
    this._protHB = new Uint8Array(nProt);
    beads.forEach((b, i) => {
      const cls = RES_CLASS_OF[b.resName] ?? "H";
      this.resClass[i] = ["H", "A", "P", "Cp", "Cn"].indexOf(cls);
      const p = RES_CLASS[cls];
      this._protSigma[i] = p.sigma; this._protEps[i] = p.eps; this._protQ[i] = p.q;
      // Loop-2 S1 (R2 term b): opt-in formal charges override the class-level
      // q=0 with the residue-identity value (ASP/GLU −1, LYS/ARG +1, HIS 0).
      if (this.chargesOn) this._protQ[i] = CG_FORMAL_CHARGES[b.resName] ?? 0;
      this._protHB[i] = (cls === "P" || cls === "Cp" || cls === "Cn") ? 1 : 0;
    });
    // Loop-2 S2: virtual interaction sites (backbone O acceptors for ALL
    // residues — fixes the 52% class-A gap where only P/Cp/Cn beads could
    // ever accept; hydrophobic backbones accept too in reality).
    this._vSites = null;
    if (this.hbMode === "directional") {
      this._vSites = buildVirtualSites(beads);
      beads.forEach((b, i) => {
        const s = this._vSites[i];
        if (s.valid) coneAxisOf(s, [b.x, b.y, b.z]);
      });
      this._vSiteValence = new Uint8Array(nProt); // per-frame acceptor valence caps
    }
    // Per-ligand-atom element params
    this._ligSigma = new Float64Array(this.nLigAtoms);
    this._ligEps = new Float64Array(this.nLigAtoms);
    this._ligQ = new Float64Array(this.nLigAtoms);
    this._ligHB = new Uint8Array(this.nLigAtoms);
    this._ligdG = new Float64Array(this.nLigAtoms);
    for (let a = 0; a < this.nLigAtoms; a++) {
      const e = LIG_ELEMENT[this.ligandAtoms[a].element] ?? LIG_ELEMENT_DEFAULT;
      this._ligSigma[a] = e.sigma; this._ligEps[a] = e.eps; this._ligQ[a] = e.q;
      this._ligHB[a] = e.hb ? 1 : 0; this._ligdG[a] = e.dG;
    }

    // ---- 1-2 bonds & 1-2-3 angles from segments --------------------------
    /** bonds[idx] = [i, j, r0] packed flat as [i0,j0,r0, i1,j1,r1, ...] */
    const bonds = [];
    const angles = []; // [i, j, k, theta0]
    for (const [s, e] of segments) {
      for (let i = s; i + 1 < e; i++) {
        bonds.push(i, i + 1, this._dist(this.ref, i, i + 1));
      }
      for (let i = s; i + 2 < e; i++) {
        angles.push(i, i + 1, i + 2, this._angle(this.ref, i, i + 1, i + 2));
      }
    }
    this.bonds = new Float64Array(bonds);
    this.angles = new Float64Array(angles);

    // ---- ENM contact set (also used for drawing) --------------------------
    // Pair (i,j) gets a spring if r0(i,j) <= rc and |i-j| > 2 within the same
    // segment-family OR i,j belong to different segments/chains (interfacial
    // contacts in complexes matter and must be elastic too).
    // Protein pairs only — ligand geometry is frozen by its own internal FF,
    // and protein–ligand elastic contacts are left to the excluded volume.
    const springs = [];      // [i, j, r0]
    const isBound13 = (i, j) => {
      for (const [s, e] of segments) {
        if (i >= s && j < e) return j - i <= 2; // 1-2 and 1-3 excluded
        if (j >= s && i < e) return i - j <= 2;
      }
      return false;
    };
    for (let i = 0; i < nProt; i++) {
      for (let j = i + 1; j < nProt; j++) {
        const r0 = this._dist(this.ref, i, j);
        if (r0 <= this.rc && !isBound13(i, j)) {
          springs.push(i, j, r0);
        }
      }
    }
    this.springs = new Float64Array(springs);
    // Per-spring force constants (kcal/mol/Å²). Starts uniform at this.gamma;
    // an ML contact prior (contacts.json from ml/export_esm_contacts.py) can
    // rescale individual pairs via setSpringScale() — e.g. stiffen residues
    // the model predicts to be in contact.
    this.springK = new Float64Array(this.springs.length / 3).fill(this.gamma);
    this.springScaleActive = false;
    // ── Sequence-dependent ENM (Bahar-style, stub) ──────────────────────
    // Alternative springK_seq: K = gamma * (1 + 0.2*(w_i + w_j)/2)
    // where w_i = SEQ_WEIGHT[resClass(i)], w_j = SEQ_WEIGHT[resClass(j)].
    // H=1.0, A=0.9, P=1.1, Cp/Cn=1.05 (see src/ff-params.js:SEQ_WEIGHT).
    // Mean (w_i+w_j)/2 ≈1.0 so ⟨K_seq⟩≈gamma; modulation ±2% typical, ±4%
    // extremes. Preserves ENM average flexibility while encoding chemistry
    // (hydrophobic vs polar contacts). Disabled by default for backward
    // compat; enable via par.seqWeight or ForceField.applySeqWeights().
    // Example:
    //   const w_i = SEQ_WEIGHT[RES_CLASS_OF[beads[i].resName]] ?? 1.0;
    //   const w_j = SEQ_WEIGHT[RES_CLASS_OF[beads[j].resName]] ?? 1.0;
    //   const K_seq = this.gamma * (1 + 0.2 * (w_i + w_j) / 2);
    // See docs/CG_HEAVY.md and tests/test_enm_seq.js for validation.
    // If par.seqWeight is true, build springK_seq instead of uniform:
    if (par.seqWeight) this.applySeqWeights();
    // ---- Tirion distance-weighted ENM (opt-in, Phase 1) --------------------
    // par.enmModel === "tirion" (or par.tirion === true) replaces the uniform
    // springK with γ_ij = γ0·(R0/r0_ij)^6 plus SS-dependent backbone basins.
    // Default ("uniform") preserves legacy behavior so existing tests pass.
    // See src/physics/forcefield/tirion_anm.js.
    this.enmModel = par.enmModel ?? (par.tirion ? "tirion" : "uniform");
    this.tirionDihedrals = new Float64Array(0);
    this.tirionDihedralK = new Float64Array(0);
    this.tirionSS = [];
    if (this.enmModel === "tirion") {
      try {
        this.useTirionNetwork({ gamma0: this.gamma, cutoff: this.rc, segments });
      } catch (e) {
        console.warn(`[ForceField] Tirion network failed (${e.message}) — falling back to uniform ENM`);
        this.enmModel = "uniform";
      }
    }

    // ---- holo contact springs ----------------------------------------------
    // Harmonic protein–ligand springs for every pair whose NATIVE distance
    // r0 <= 6.0 Å, encoding the observed bound pose of a holo complex (e.g. T4
    // lysozyme L99A + benzene, PDB 4W52): pins the ligand into its binding
    // pocket while still allowing fluctuation. Gated by par.binding.holo
    // (default true); no-op when no ligands are present.
    this.holoSprings = new Float64Array(0);
    this.nHolo = 0;
    this.rebuildHoloSprings();

    // ---- repulsive-pair bookkeeping ---------------------------------------
    // Non-repulsive pairs (bonded 1-2, 1-3 and ENM springs) are excluded to
    // avoid double counting; the grid only tests the remainder.
    this._excluded = new Set();
    for (let k = 0; k < this.bonds.length; k += 3) {
      this._excluded.add(this._pairKey(this.bonds[k], this.bonds[k + 1]));
    }
    for (let k = 0; k < this.angles.length; k += 4) {
      this._excluded.add(this._pairKey(this.angles[k], this.angles[k + 2]));
    }
    for (let k = 0; k < this.springs.length; k += 3) {
      this._excluded.add(this._pairKey(this.springs[k], this.springs[k + 1]));
    }
    // Holo contact springs govern native protein–ligand contact distances, so
    // those pairs are excluded from the binding pair pass to avoid double
    // counting (the spring alone holds them at r0).
    for (let k = 0; k < this.holoSprings.length; k += 3) {
      this._excluded.add(this._pairKey(this.holoSprings[k], this.holoSprings[k + 1]));
    }
    // Covalently linked ligand pairs (1-2 bonds, 1-3 angles, 1-4 improper
    // partners) must never be repelled by the grid: the internal FF already
    // governs those distances and they sit far inside the repulsion range
    // (e.g. para carbons of an aromatic ring at ~2.8 Å ≪ r_e ≈ 5.0 Å).
    for (let a = 0; a < this.ligandBonds.length; a += 3) {
      this._excluded.add(this._pairKey(this.ligandBonds[a], this.ligandBonds[a + 1]));
    }
    for (let a = 0; a < this.ligandAngles.length; a += 4) {
      this._excluded.add(this._pairKey(this.ligandAngles[a], this.ligandAngles[a + 2]));
    }
    for (let a = 0; a < this.ligandImpropers.length; a += 5) {
      const i = this.ligandImpropers[a], j = this.ligandImpropers[a + 1];
      const l = this.ligandImpropers[a + 3];
      this._excluded.add(this._pairKey(i, l));
      this._excluded.add(this._pairKey(j, l));
      this._excluded.add(this._pairKey(this.ligandImpropers[a + 2], l));
    }
    // Ligand 1-4 pairs (three bonds apart, e.g. the N…N of a N–C–C–N fragment)
    // are NOT repelled by the grid: gauche 1-4 conformers sit at ~2.9 Å — far
    // inside the Cα-sized repulsion range — and with no explicit non-ring
    // torsion term they are legitimately free. Only ring 1-4 pairs were
    // already covered by the improper exclusions above; this walk covers every
    // 1-4 pair in any molecule (incl. poly-atomic buffers like HEPES).
    if (this.ligandBonds.length) {      const adjLig = new Map();
      for (let a = 0; a < this.ligandBonds.length; a += 3) {
        const i = this.ligandBonds[a], j = this.ligandBonds[a + 1];
        if (!adjLig.has(i)) adjLig.set(i, []);
        if (!adjLig.has(j)) adjLig.set(j, []);
        adjLig.get(i).push(j); adjLig.get(j).push(i);
      }
      for (const [x, nbrs] of adjLig) {
        for (const a of nbrs) {
          const nb = adjLig.get(a);
          if (!nb) continue;
          for (const b of nb) {
            const nb2 = adjLig.get(b);
            if (!nb2) continue;
            for (const c of nb2) {
              if (c !== x && c !== a) this._excluded.add(this._pairKey(x, c));
            }
          }
        }
      }
    }

    // Per-particle repulsive σ: protein beads use the Cα size; ligand atoms use
    // their element's LJ σ (binding tables) so intra-ligand excluded volume is
    // not over-estimated — a folded 15-atom co-solute must not sit at 10 kcal.
    // Allocated here (before the native-contact scan, which reads the table).
    this._repSigma = new Float64Array(this.n);
    for (let i = 0; i < this.nProt; i++) this._repSigma[i] = this.sigmaRep;
    for (let a = 0; a < this.nLigAtoms; a++) this._repSigma[this.nProt + a] = this._ligSigma[a];

    // ---- Gö-like native contacts (ligand 1-5+ pairs) -----------------------
    // Folded ligands (sugars, inhibitors, buffers) place non-bonded 1-5+
    // pairs at 2.9–4.0 Å — inside their own LJ repulsion range r_e. Without
    // treatment the grid sees a huge native strain (PDBBind probes: U(native)
    // up to ~745 kcal/mol, max|F| ~140 kcal/mol/Å) that launches the light
    // united atoms as projectiles and blasts the protein. Standard CG fix:
    // every non-excluded pair whose NATIVE distance lies inside its own
    // repulsion range is excluded from the grid and gets a soft Gö-like
    // harmonic contact spring (k = 1.0 kcal/mol/Å², [i,j,r0] packing like
    // holoSprings) that preserves the native distance while still resisting
    // compression. Applied to ANY particle type (L–L, P–P, P–L) — a genuine
    // native clash is a modeling problem regardless of chemistry.
    //
    // Double-coverage diagram (see docs/CG_HEAVY.md: nativeContacts):
    //   native r0 ──►  [excluded from grid repulsion?] ──► [already ENM spring?]
    //         │                     │                               │
    //         │  r0 < r_e && not in _excluded                yes → keep ENM only
    //         │         │                               no → create nativeContacts spring
    //         │         └─► yes → push (i,j,r0) to nativeContacts, add to _excluded
    //         │                     (prevents 12-6 repulsion + harmonic double count)
    //         └─► r0 ≥ r_e → leave to grid repulsion (smooth WCA, zero at r_e)
    //   Intra-molecule only (P–P, L–L): cross P–L pairs are holo-pinned or
    //   binding-evaluated, never grid-repelled, so no Gö spring there (would
    //   double-count attractive binding). The scan is O(N²) at construction
    //   only; compute() then sees each pair exactly once (CG_HEAVY.md Fig.1).
    this.nativeContacts = new Float64Array(0); // docs/CG_HEAVY.md: nativeContacts double-coverage
    {
      const nx = [];
      for (let i = 0; i < this.n; i++) {
        for (let j = i + 1; j < this.n; j++) {
          // Intra-molecule scan only: cross protein–ligand pairs never reach
          // the grid repulsion (they are either holo-spring pinned — excluded
          // from binding AND grid — or evaluated by the binding pass with
          // attractive wells), so adding a Gö-like spring there would double-
          // count the interaction.
          if ((i < this.nProt) !== (j < this.nProt)) continue;
          if (this._excluded.has(this._pairKey(i, j))) continue;
          const r0 = this._dist(this.ref, i, j);
          const s = 0.5 * (this._repSigma[i] + this._repSigma[j]);
          const re = (2 ** (1 / 6)) * s;
          if (r0 < re) nx.push(i, j, r0);
        }
      }
      if (nx.length) {
        this.nativeContacts = new Float64Array(nx);
        for (let k = 0; k < nx.length; k += 3) {
          this._excluded.add(this._pairKey(nx[k], nx[k + 1]));
        }
      }
    }

    // Scratch / grid buffers (spatial hash with numeric keys — GC-free)
    this.forces = new Float64Array(this.n * 3);
    this.rcRep = (2 ** (1 / 6)) * this.sigmaRep; // r_e = 2^(1/6) σ ≈ 5.61 Å
    this._cell = this.rcRep;                    // cell size ≥ repulsive range
    this._grid = new Map();                     // hashCellKey -> bead index list
    this._gridB = new Map();                    // protein-only grid for binding pass
    this.energy = 0; // last computed potential energy (kcal/mol)
    this.bindingU = 0;   // protein–ligand nonbonded energy (last compute)
    this.desolvU = 0;    // EEF1-lite burial/desolvation energy (last compute)
    // Loop-2 S4 (R4 §5 item 1, R6 §5): per-term binding-accumulator opt-in.
    // trackTerms = true fills bindLJU/bindCoulU/bindHBU/desolvU + the
    // bindU vector {lj, coul, hb, desolv} each compute() — consumed by the
    // BindLog energy channel (main.js tick). DEFAULT false: kernel skips
    // all tracker branches → bit-identical to pre-S4, zero overhead.
    this.trackTerms = false;
    this.bindLJU = 0; this.bindCoulU = 0; this.bindHBU = 0;
    this.bindU = { lj: 0, coul: 0, hb: 0, desolv: 0 };

    // EEF1-lite desolvation scratch (reused every compute call — GC-free).
    // dens[a]  = soft protein-occupancy count n_a (Pass 1 accumulation)
    // dBdn[a]  = dB/dn = exp(−n_a/3)/3 at the current density (Pass 2 lookup)
    // bpA/bpJ/bpR = flat pair records (ligand table index, protein bead index,
    //   distance r) written in Pass 1 so Pass 2 adds forces without a re-scan.
    //   Parallel Float64Arrays are chosen over object arrays to stay allocation-
    //   free. The worst case is every ligand atom paired with every protein bead,
    //   so nProt·nLigAtoms is a hard size cap and no growth is ever needed.
    this._dens = new Float64Array(Math.max(1, this.nLigAtoms));
    this._dBdn = new Float64Array(Math.max(1, this.nLigAtoms));
    const maxBP = this.nLigAtoms ? nProt * this.nLigAtoms : 1;
    this._bpA = new Float64Array(maxBP);
    this._bpJ = new Float64Array(maxBP);
    this._bpR = new Float64Array(maxBP);
  }

  _pairKey(i, j) { return i < j ? i * 1e6 + j : j * 1e6 + i; }

  _dist(p, i, j) {
    const dx = p[3 * j] - p[3 * i], dy = p[3 * j + 1] - p[3 * i + 1], dz = p[3 * j + 2] - p[3 * i + 2];
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  _angle(p, i, j, k) {
    const ax = p[3 * i] - p[3 * j], ay = p[3 * i + 1] - p[3 * j + 1], az = p[3 * i + 2] - p[3 * j + 2];
    const bx = p[3 * k] - p[3 * j], by = p[3 * k + 1] - p[3 * j + 1], bz = p[3 * k + 2] - p[3 * j + 2];
    const la = Math.hypot(ax, ay, az), lb = Math.hypot(bx, by, bz);
    let c = (ax * bx + ay * by + az * bz) / (la * lb);
    c = Math.min(1, Math.max(-1, c));
    return Math.acos(c);
  }

  /**
   * compute(pos) → fills this.forces and returns total potential energy.
   * pos: Float64Array(3n). Forces are −∇U in kcal/mol/Å.
   *
   * G65 ZERO-ALLOC AUDIT — hot loop (called every integration step):
   *   Expected heap allocs per compute(): ~0 in steady state.
   *   - f.fill(0) reuses preallocated Float64Array(this.forces) — no alloc.
   *   - harmonicPairs / springForces / angleForces / ligandInternal: pure loops on flat
   *     Float64Array views; only scalar locals (no `new` inside loops).
   *   - _repulsion / _binding: uniform-grid cell list uses preallocated Maps
   *     (`this._grid`, `this._gridB`) and reused scratch buffers (`_dens`, `_dBdn`,
   *     `_bpA/J/R`, `_repSigma`); `arr.length=0` clears in place, `arr.push(i)` reuses
   *     existing bucket arrays (capacity grows once, then stable). The `for (const [k,arr] of grid)`
   *     iterator does allocate a lightweight iterator object (~O(cells)), but no per-pair arrays.
   *   - No `new Float64Array` or `new Array` inside the hot path after construction.
   *   Measured via bench/alloc.js (heap delta over 1000 steps). If a future change
   *   adds `new Float64Array` inside compute(), hoist it to a scratch buffer (see HeavyForceField/SasaModel fix below).
   */
  compute(pos) {
    const f = this.forces;
    f.fill(0);
    let U = 0;

    // --- 1/3. Harmonic two-body terms (bonds + ENM springs share the kernel)
    U += this._harmonicPairs(pos, f, this.bonds, 3, this.kBond);
    U += this._springForces(pos, f);
    if (this.holoSprings.length) U += this._harmonicPairs(pos, f, this.holoSprings, 3, this.holoGamma);
    if (this.nativeContacts.length) U += this._harmonicPairs(pos, f, this.nativeContacts, 3, 1.0);
    // One-sided compression floor on holo-pinned pairs — keeps the funnel or
    // desolvation from ever squeezing the ligand through a pocket wall
    // (holo pairs are excluded from both grid repulsion and the binding pass,
    // so without this term nothing resists r < r_min).
    if (this.holoSprings.length) {
      const HS = this.holoSprings, RMIN = 2.6, KF = 8.0;
      for (let a = 0; a < HS.length; a += 3) {
        const i = 3 * HS[a], j = 3 * HS[a + 1];
        const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
        if (r >= RMIN) continue;
        const dr = r - RMIN;                      // < 0
        U += 0.5 * KF * dr * dr;
        const s = (KF * dr) / r;                  // dU/dr ÷ r (dr<0 ⇒ repulsive)
        const fx = s * dx, fy = s * dy, fz = s * dz;
        f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
        f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
      }
    }

    // --- 2. Angle bending ---------------------------------------------------
    U += this._angleForces(pos, f);

    // --- 3b. Ligand united-atom internal terms (no-op without ligands) ------
    U += this._ligandInternal(pos, f);

    // --- 4. Excluded volume via grid ----------------------------------------
    U += this._repulsion(pos, f);

    // --- 5. Protein–ligand binding (cross LJ, electrostatics, H-bonds) -------
    // (physical documentation lives on the _binding wrapper below; the
    //  kernel itself is in ff-binding.js)
    U += this._binding(pos, f);

    // --- 5c. Optional funnel + well-tempered metadynamics bias (off by default)
    if (this.funnel && this.funnelOn) U += this.funnel.addForces(pos, f);

    // --- 6. NaN/∞ guard ------------------------------------------------------
    // A non-finite coordinate or force kills every downstream pair test and
    // (worse) makes the spatial-hash Map grow unboundedly via garbage cell
    // keys. Detect here (zero cost in the normal case: a single isFinite on
    // the accumulated energy plus a fused isfinite-OR loop on the force
    // buffer), sanitize the force buffer, and leave this.energy = NaN so the
    // app can auto-pause on the next frame.
    if (!Number.isFinite(U)) {
      let bad = 0;
      for (let i = 0; i < f.length; i++) {
        if (!Number.isFinite(f[i])) { f[i] = 0; bad++; }
      }
      this.nanStrikes = (this.nanStrikes ?? 0) + 1;
      if (bad === 0) {
        // energy NaN but forces finite — scan positions for the culprit
        for (let i = 0; i < pos.length && bad < 4; i++) {
          if (!Number.isFinite(pos[i])) bad++;
        }
      }
      U = NaN;
    }
    this.energy = U;
    return U;
  }

  /** Attach an optional Funnel instance (constructed by main.js, not imported here). */
  setFunnel(fn) { this.funnel = fn; }

  /**
   * Rebuild the holo (native-pose) protein–ligand springs from the CURRENT
   * reference coordinates. In the constructor this derives the initial holo
   * springs from the native reference; after a placement it is called with
   * holoOn=false to clear them (a placed library ligand is a hypothesis, not
   * a known crystallographic pose, so it is not pinned). Also keeps the
   * non-repulsive exclusion set in sync with whatever holo pairs exist.
   * No-op when holo springs are off or there are no ligand atoms.
   */
  rebuildHoloSprings() {
    // remove the old holo pairs from the non-repulsive exclusion set
    // (guarded: the constructor calls this before _excluded is created)
    if (this._excluded) {
      for (let k = 0; k < this.holoSprings.length; k += 3) {
        this._excluded.delete(this._pairKey(this.holoSprings[k], this.holoSprings[k + 1]));
      }
    }
    this.holoSprings = new Float64Array(0);
    this.nHolo = 0;
    if (!this.holoOn || this.nLigAtoms <= 0) return;
    const hs = [];
    for (let i = 0; i < this.nProt; i++) {
      for (let la = this.nProt; la < this.n; la++) {
        const r0 = this._dist(this.ref, i, la);
        if (r0 <= 6.0) hs.push(i, la, r0);
      }
    }
    this.holoSprings = new Float64Array(hs);
    this.nHolo = hs.length / 3;
    if (this._excluded) {
      for (let k = 0; k < this.holoSprings.length; k += 3) {
        this._excluded.add(this._pairKey(this.holoSprings[k], this.holoSprings[k + 1]));
      }
    }
  }

  /** Σ 1⁄2 k (r−r0)2 over a flat pair list; kernel lives in ff-harmonic.js. */
  _harmonicPairs(pos, f, list, stride, k) { return harmonicPairs(pos, f, list, stride, k); }

  /** ENM spring forces with per-spring k — kernel in ff-harmonic.js. */
  _springForces(pos, f) { return springForces(this, pos, f); }

  /**
   * Apply an ML contact prior to the ENM spring constants.
   * @param {Array} contacts  [[i, j, p], ...] residue pairs (0-indexed) with p ∈ [0,1].
   * @param {number} alpha    global scale (k_per_pair = gamma · (1 + alpha·p)).
   * Marks springScaleActive so the UI knows a custom map is live.
   */
  setSpringScale(contacts, alpha = 1) {
    const S = this.springs, K = this.springK, gamma = this.gamma;
    // reset to baseline gamma, then rescale by the prior
    K.fill(gamma);
    const pairScale = new Map();
    for (const [i, j, p] of contacts) {
      if (p <= 0) continue;
      pairScale.set(this._pairKey(i, j), 1 + alpha * p);
    }
    for (let k = 0, s = 0; k < S.length; k += 3, s++) {
      const sc = pairScale.get(this._pairKey(S[k], S[k + 1]));
      if (sc !== undefined) K[s] = gamma * sc;
    }
    this.springScaleActive = true;
  }

  /** Restore uniform ENM stiffness (undo an ML contact-prior scaling). */
  clearSpringScale() {
    this.springK.fill(this.gamma);
    this.springScaleActive = false;
  }

  /**
   * Sequence-dependent ENM (Bahar-style) — opt-in stub.
   * Rebuilds springK as K_seq(i,j) = gamma * (1 + 0.2*(w_i + w_j)/2)
   * where w_i = SEQ_WEIGHT[RES_CLASS_OF[beads[i].resName]] (w_j likewise).
   * H=1.0, A=0.9, P=1.1, Cp/Cn=1.05; mean modulation ~0 so ⟨K⟩≈gamma.
   * Call after construction or pass par.seqWeight=true. Preserves the
   * uniform topology (same springs list, same H(Rc−r0) cutoff) — only the
   * stiffness is chemistry-weighted. See src/ff-params.js:SEQ_WEIGHT,
   * docs/CG_HEAVY.md and tests/test_enm_seq.js. Also available as
   * par.seqWeight flag in constructor.
   * @param {Array} beads  optional override (defaults to constructor beads if stored)
   */
  applySeqWeights(beads = null) {
    // beads not stored on ForceField; caller may pass original bead list,
    // otherwise we fall back to resClass table which was built from beads.
    // resClass is Uint8Array over 5 classes; invert to weight via SEQ_WEIGHT.
    const weightByClass = [SEQ_WEIGHT.H, SEQ_WEIGHT.A, SEQ_WEIGHT.P, SEQ_WEIGHT.Cp, SEQ_WEIGHT.Cn];
    // w_i per bead from this.resClass (built in constructor)
    const nProt = this.nProt;
    const wPerBead = new Float64Array(nProt);
    for (let i = 0; i < nProt; i++) {
      // resClass[i] is 0:H,1:A,2:P,3:Cp,4:Cn (see constructor)
      wPerBead[i] = weightByClass[this.resClass[i]] ?? 1.0;
    }
    // If caller supplied explicit beads (e.g. for test), prefer SEQ_WEIGHT via resName
    if (beads && beads.length === nProt) {
      for (let i = 0; i < nProt; i++) {
        const cls = RES_CLASS_OF[beads[i].resName] ?? "H";
        const w = SEQ_WEIGHT[cls];
        if (w !== undefined) wPerBead[i] = w;
      }
    }
    const S = this.springs, K = this.springK, gamma = this.gamma;
    for (let k = 0, s = 0; k < S.length; k += 3, s++) {
      const i = S[k], j = S[k + 1];
      const w_i = wPerBead[i], w_j = wPerBead[j];
      // K_seq = gamma * (1 + 0.2*(w_i + w_j)/2)  — Bahar-style ±2% modulation
      K[s] = gamma * (1 + 0.2 * (w_i + w_j) * 0.5);
    }
    this.springScaleActive = true; // marks non-uniform (reuse flag for UI)
    // For debugging: this._seqW = wPerBead; // keep per-bead w_i if needed
  }

  /**
   * Opt-in Tirion distance-weighted ENM + SS dihedral basins (Phase 1).
   * Rebuilds springs/springK as γ_ij = γ0·(R0/r0_ij)^6 over the same
   * H(Rc−r0) topology and stores backbone pseudo-dihedral basins on
   * ff.tirionDihedrals / ff.tirionDihedralK / ff.tirionSS. Falls back to the
   * existing uniform network when the module is unavailable.
   * @param {object} [opts]  { gamma0, R0, cutoff, segments }
   */
  useTirionNetwork(opts = {}) {
    const ref = this.ref.subarray(0, this.nProt * 3);
    // Reconstruct segments if the caller did not pass them: derive from the
    // current bond list (bonds only join contiguous pairs).
    let segments = opts.segments ?? null;
    if (!segments) {
      segments = [];
      if (this.nProt > 0) {
        let s = 0;
        const bondedNext = new Set();
        for (let k = 0; k < this.bonds.length; k += 3) {
          const i = this.bonds[k], j = this.bonds[k + 1];
          if (j === i + 1) bondedNext.add(i);
        }
        for (let i = 0; i < this.nProt - 1; i++) {
          if (!bondedNext.has(i)) { segments.push([s, i + 1]); s = i + 1; }
        }
        segments.push([s, this.nProt]);
      }
    }
    const net = buildTirionNetwork(ref, {
      gamma0: opts.gamma0 ?? this.gamma,
      R0: opts.R0 ?? 3.81,
      cutoff: opts.cutoff ?? this.rc,
      segments,
    });
    // Replace uniform springs; keep exclusion set consistent.
    if (this._excluded) {
      for (let k = 0; k < this.springs.length; k += 3) {
        this._excluded.delete(this._pairKey(this.springs[k], this.springs[k + 1]));
      }
    }
    applyTirionToForceField(this, net);
    this.enmModel = "tirion";
  }

  /** U_θ = ½ kθ (θ−θ0)² — angle-bending kernel lives in ff-harmonic.js. */
  _angleForces(pos, f, list = this.angles, k = this.kAngle) { return angleForces(this, pos, f, list, k); }

  /**
   * Ligand united-atom internal energy (bonds, angles, improper planarity).
   * No-op (returns 0) when no ligands are present.
   */
  _ligandInternal(pos, f) {
    if (this.ligandAtoms.length === 0) return 0;
    let U = 0;
    U += this._ligandBondForces(pos, f);
    U += this._angleForces(pos, f, this.ligandAngles, 40);
    U += this._improperForces(pos, f, this.ligandImpropers);
    return U;
  }

  /** Ligand bonds U = Σ ½ k (r−r0)² — kernel in ff-harmonic.js. */
  _ligandBondForces(pos, f) { return ligandBondForces(this, pos, f); }

  /** Improper (out-of-plane) term — FD kernel lives in ff-harmonic.js. */
  _improperForces(pos, f, list) { return improperForces(pos, f, list); }

  /**
   * Protein–ligand binding potential — pair pass over every protein bead
   * (0..nProt−1) vs every ligand atom (global la = nProt + a, table a).
   * Three short-range cross terms, each cut off smoothly at bindRcut:
   *
   *   U_LJ = 4ε[(σ/r)¹² − (σ/r)⁶]·sw            cross 12-6 (attractive well)
   *   U_el = 332.0637·q1·q2/(εr·r)·sw            screened electrostatics,
   *                                              εr(r) = 4 + 76·tanh(r/8)
   *   U_HB = −epsHB·exp(−(r−r0)²/2w²)·sw          H-bond well, r0=3.2, w=0.6
   *
   * with Lorentz–Berthelot combining rules σ = (σi+σa)/2, ε = √(εi·εa).
   *
   * SWITCHING: the smoothstep is identically 1 for r ≤ rsw and falls to 0 at
   * rc with a vanishing derivative, so both energy and force are continuous:
   *     t   = (r − rsw)/(rc − rsw),   Δ = rc − rsw
   *     s   = 1 − t³(10 − 15t + 6t²)          [cubic smoothstep, s(0)=1, s(1)=0]
   *     s′  = ds/dr = −(30/Δ)·t²(1−t)²
   * (ds/dt = −30t²(1−t)² by direct differentiation; dt/dr = 1/Δ by chain rule.)
   *
   * SIGN CONVENTION: the force on protein bead i is −∂U/∂pos_i. With
   * dx = pos_i − pos_la pointing from the ligand atom to the bead,
   * ∂r/∂pos_i = dx/r, so f_i += −(dU/dr)·(dx/r) and f_la −= that (Newton's
   * 3rd law). Each term accumulates its dU/dr into one scalar and a single
   * vector update is applied — no sign ambiguity.
   *
   * EEF1-LITE BURIAL (Pass 2): the same grid scan additionally accumulates a
   * soft protein-occupancy count for each ligand atom a,
   *     n_a = Σ_j g(r_aj),   g(r) = exp(−(r−r0)²/(2σ²)),   r0=4.5 Å, σ=1.8 Å
   * (g has decayed to ≈ 0.04 at r = 9 Å = bindRcut, so no second cutoff is
   * needed). The burial fraction
   *     B_a = 1 − exp(−n_a/3) ∈ [0,1]     (0 exposed, 1 fully buried)
   * turns the per-atom hydrophobic transfer free energy ΔG_a into the burial
   * cost U_desolv = Σ_a ΔG_a·B_a (ΔG_a < 0 ⇒ burial lowers the binding energy).
   * The recorded (a, j, r) pairs from Pass 1 feed the force without re-scanning:
   *     dU_desolv/dr = ΔG_a·(dB/dn)·(dg/dr)
   *     dB/dn = exp(−n_a/3)/3 ,   dg/dr = −g(r)·(r−r0)/σ²
   * applied with the same sign convention as the pair terms above. Passing
   * through a burial maximum at r0 gives the correct physics: a fully exposed
   * atom (no stored pairs) feels no desolvation force, an approaching atom is
   * pulled into the protein (favorable ΔG·B), and one pushed inside below r0
   * is pushed back out. The ± sign is verified numerically against total-U
   * finite differences in the test suite.
   *
   * Two-pass structure:
   *   Pass 1 — pair terms (LJ / electrostatics / H-bond) computed exactly as
   *            before, plus dens[a] accumulation and (a, j, r) recording.
   *   Pass 2 — burial fraction/energy per ligand atom, then desolvation forces
   *            over the recorded pairs.
   */
  _binding(pos, f) { return binding(this, pos, f); }

  /** Repulsive-only 12-6 LJ excluded-volume pass — kernel in ff-repulsion.js. */
  _repulsion(pos, f) { return repulsion(this, pos, f); }

  /* --- numeric spatial-hash cell keys (no string allocs) -------------- */
  _encodeCell(cx, cy, cz) { return ((cx + 2048) * 4096 + (cy + 2048)) * 4096 + (cz + 2048); }
  _cellKey(x, y, z, cell) {
    // Non-finite coordinates must never reach the grid: Math.floor(NaN) is
    // NaN and the Map would grow unboundedly with a fresh garbage key on
    // every compute call (the 1HVR hang). Route everything non-finite to a
    // single sentinel cell at the grid corner — the pair distance tests will
    // then reject those particles (r2 is NaN, comparisons false).
    if (!Number.isFinite(x + y + z)) return 0;
    return this._encodeCell(Math.floor(x / cell), Math.floor(y / cell), Math.floor(z / cell));
  }
  _decodeX(k) { return Math.floor(k / (4096 * 4096)) - 2048; }
  _decodeY(k) { return Math.floor(k / 4096) % 4096 - 2048; }
  _decodeZ(k) { return k % 4096 - 2048; }

  /* ------------------------------------------------------------------ */
  /*  Analysis helpers (HUD)                                             */
  /* ------------------------------------------------------------------ */

  /**
   * Instantaneous kinetic temperature from equipartition:
   *   ke_mech = ½ Σ m v²   (Da·Å²/ps²)  →  /KCONV → kcal/mol
   *   T = ke_kcal / (3/2 N k_B)
   * @param {Float64Array} vel   velocities, one component per coordinate (3n)
   * @param {Float64Array} mass  mass per coordinate (3n) — ligand atoms carry
   *   their united-atom mass, protein beads the Cα mass. Includes all particles.
   */
  kineticTemp(vel, mass) {
    let ke = 0;
    for (let i = 0; i < this.n * 3; i++) ke += mass[i] * vel[i] * vel[i];
    ke *= 0.5 / KCONV;                   // → kcal/mol
    return ke / (1.5 * this.n * KB_KCAL);
  }

  /**
   * RMSD to native, restricted to the protein Cα beads (ligand coords are
   * rigid internal DOF and excluded from the fold-space metric). No alignment;
   * ENM keeps the COM/orientation nearly fixed.
   */
  rmsd(pos) {
    let s = 0;
    const r = this.ref;
    for (let i = 0; i < this.nProt * 3; i++) {
      const d = pos[i] - r[i];
      s += d * d;
    }
    return Math.sqrt(s / this.nProt);
  }
}
