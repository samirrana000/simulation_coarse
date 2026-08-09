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

import { buildLigandInternalFF, improperAngle } from "./ligand.js?v=8";

export const KB_KCAL = 0.0019872041; // Boltzmann constant, kcal mol^-1 K^-1

/**
 * Unit-conversion constant between (kcal/mol, Å) and mechanical (Da, Å, ps)
 * units used by the integrator:
 *     1 kcal/mol            = 418.4 Da·Å²/ps²
 *     a[Å/ps²]              = KCONV · F[kcal/mol/Å] / m[Da]
 *     k_B in mech. units    = KB_KCAL · KCONV = 0.8314 Da·Å²/(ps²·K)
 * (derivation: 1 kcal/mol·Å = 6.9477e-11 N ≙ 418.4 Da·Å/ps² since 1 N =
 * 6.022e15 Da·Å/ps².)
 */
export const KCONV = 418.4;

// Residue class → protein-bead LJ parameters (σ Å, ε kcal/mol, charge e)
const RES_CLASS = {
  H:  { sigma: 4.0, eps: 0.15, q: 0 },  // hydrophobic
  A:  { sigma: 4.1, eps: 0.18, q: 0 },  // aromatic
  P:  { sigma: 3.8, eps: 0.12, q: 0 },  // polar (H-bond capable)
  Cp: { sigma: 3.6, eps: 0.10, q: 0 },  // positively charged (H-bond capable)
  Cn: { sigma: 3.6, eps: 0.10, q: 0 },  // negatively charged (H-bond capable)
};
const RES_CLASS_OF = {
  ALA: "H", VAL: "H", LEU: "H", ILE: "H", PRO: "H", MET: "H", GLY: "H", CYS: "H",
  PHE: "A", TRP: "A", TYR: "A", HIS: "A",
  SER: "P", THR: "P", ASN: "P", GLN: "P",
  LYS: "Cp", ARG: "Cp",
  ASP: "Cn", GLU: "Cn",
};
// Ligand element → LJ params (σ Å, ε kcal/mol), partial charge (e), H-bond flag, ΔG desolvation (kcal/mol)
const LIG_ELEMENT = {
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
const LIG_ELEMENT_DEFAULT = { sigma: 3.4, eps: 0.12, q: 0.0, hb: false, dG: -0.30 };

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
    this.kBond = par.kBond ?? 100.0;   // kcal/mol/Å²  (Cα–Cα peptide bond)
    this.kAngle = par.kAngle ?? 20.0;  // kcal/mol/rad² (backbone pseudo-angle)
    this.epsRep = par.epsRep ?? 0.3;   // kcal/mol      (excluded-volume depth)
    this.sigmaRep = par.sigmaRep ?? 4.0; // Å           (Cα bead diameter ≈ 2·2.0 Å)
    this.bindOn = par.binding?.on ?? true; // protein–ligand binding potentials
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
      this._protHB[i] = (cls === "P" || cls === "Cp" || cls === "Cn") ? 1 : 0;
    });
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

    // ---- holo contact springs ----------------------------------------------
    // Harmonic protein–ligand springs for every pair whose NATIVE distance
    // r0 <= 6.0 Å, encoding the observed bound pose of a holo complex (e.g. T4
    // lysozyme L99A + benzene, PDB 4W52): pins the ligand into its binding
    // pocket while still allowing fluctuation. Gated by par.binding.holo
    // (default true); no-op when no ligands are present.
    this.holoSprings = new Float64Array(0);
    if (this.holoOn && this.nLigAtoms > 0) {
      const hs = [];
      for (let i = 0; i < this.nProt; i++) {
        for (let la = this.nProt; la < this.n; la++) {
          const r0 = this._dist(this.ref, i, la);
          if (r0 <= 6.0) hs.push(i, la, r0);
        }
      }
      this.holoSprings = new Float64Array(hs);
      this.nHolo = this.nLigAtoms > 0 ? hs.length / 3 : 0;
    }

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
    if (this.ligandBonds.length) {
      const adjLig = new Map();
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

    // Scratch / grid buffers (spatial hash with numeric keys — GC-free)
    this.forces = new Float64Array(this.n * 3);
    this.rcRep = Math.cbrt(2) * this.sigmaRep; // r_e = 2^(1/6) σ ≈ 5.6 Å
    this._cell = this.rcRep;                    // cell size ≥ repulsive range
    this._grid = new Map();                     // hashCellKey -> bead index list
    this._gridB = new Map();                    // protein-only grid for binding pass
    // Per-particle repulsive σ: protein beads use the Cα size; ligand atoms use
    // their element's LJ σ (binding tables) so intra-ligand excluded volume is
    // not over-estimated — a folded 15-atom co-solute must not sit at 10 kcal.
    this._repSigma = new Float64Array(this.n);
    for (let i = 0; i < this.nProt; i++) this._repSigma[i] = this.sigmaRep;
    for (let a = 0; a < this.nLigAtoms; a++) this._repSigma[this.nProt + a] = this._ligSigma[a];
    this.energy = 0; // last computed potential energy (kcal/mol)
    this.bindingU = 0;   // protein–ligand nonbonded energy (last compute)
    this.desolvU = 0;    // EEF1-lite burial/desolvation energy (last compute)

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
   */
  compute(pos) {
    const f = this.forces;
    f.fill(0);
    let U = 0;

    // --- 1/3. Harmonic two-body terms (bonds + ENM springs share the kernel)
    U += this._harmonicPairs(pos, f, this.bonds, 3, this.kBond);
    U += this._springForces(pos, f);
    if (this.holoSprings.length) U += this._harmonicPairs(pos, f, this.holoSprings, 3, this.holoGamma);

    // --- 2. Angle bending ---------------------------------------------------
    U += this._angleForces(pos, f);

    // --- 3b. Ligand united-atom internal terms (no-op without ligands) ------
    U += this._ligandInternal(pos, f);

    // --- 4. Excluded volume via grid ----------------------------------------
    U += this._repulsion(pos, f);

    // --- 5. Protein–ligand binding (cross LJ, electrostatics, H-bonds) -------
    U += this._binding(pos, f);

    // --- 5c. Optional funnel + well-tempered metadynamics bias (off by default)
    if (this.funnel && this.funnelOn) U += this.funnel.addForces(pos, f);

    this.energy = U;
    return U;
  }

  /** Attach an optional Funnel instance (constructed by main.js, not imported here). */
  setFunnel(fn) { this.funnel = fn; }

  /** Σ ½ k (r−r0)² over a flat pair list; returns energy, accumulates forces. */
  _harmonicPairs(pos, f, list, stride, k) {
    let U = 0;
    for (let a = 0; a < list.length; a += stride) {
      const i = 3 * list[a], j = 3 * list[a + 1], r0 = list[a + 2];
      const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
      const dr = r - r0;
      U += 0.5 * k * dr * dr;
      // Force on j: −∂U/∂r_j = −k·dr·(r̂);  on i the opposite
      const s = (k * dr) / r;
      const fx = s * dx, fy = s * dy, fz = s * dz;
      f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
      f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
    }
    return U;
  }

  /**
   * Σ ½ k_s (r−r0)² for the ENM springs, with a per-spring k from this.springK.
   * Lets an ML contact prior stiffen/weaken individual residue pairs without
   * rebuilding the whole pair list (e.g. map a model's contact probability
   * onto the ENM stiffness of that pair).
   */
  _springForces(pos, f) {
    const S = this.springs, K = this.springK;
    let U = 0;
    for (let k = 0, s = 0; k < S.length; k += 3, s++) {
      const i = 3 * S[k], j = 3 * S[k + 1], r0 = S[k + 2], kk = K[s];
      const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
      const dr = r - r0;
      U += 0.5 * kk * dr * dr;
      const sc = (kk * dr) / r;
      const fx = sc * dx, fy = sc * dy, fz = sc * dz;
      f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
      f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
    }
    return U;
  }

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
   * U_θ = ½ kθ (θ−θ0)²; forces via finite-chain-rule on cosθ.
   * Shares one kernel between the protein backbone angles and the ligand
   * angle set (the latter is stiffer: k = 40 kcal/mol/rad², see ligand.js).
   */
  _angleForces(pos, f, list = this.angles, k = this.kAngle) {
    let U = 0;
    const A = list;
    for (let a = 0; a < A.length; a += 4) {
      const i = 3 * A[a], j = 3 * A[a + 1], kk = 3 * A[a + 2], th0 = A[a + 3];
      // vectors from j
      const ax = pos[i] - pos[j], ay = pos[i + 1] - pos[j + 1], az = pos[i + 2] - pos[j + 2];
      const bx = pos[kk] - pos[j], by = pos[kk + 1] - pos[j + 1], bz = pos[kk + 2] - pos[j + 2];
      const la = Math.hypot(ax, ay, az) || 1e-12;
      const lb = Math.hypot(bx, by, bz) || 1e-12;
      let c = (ax * bx + ay * by + az * bz) / (la * lb);
      c = Math.min(1, Math.max(-1, c));
      const th = Math.acos(c);
      const dth = th - th0;
      U += 0.5 * k * dth * dth;

      // dθ/dc = −1/sinθ; guard against linear geometry
      const sinTh = Math.sqrt(Math.max(1e-12, 1 - c * c));
      const pref = (k * dth) / sinTh; // −dU/dθ · dθ/dc → pref = kθ·Δθ / sinθ

      // ∂c/∂r_i etc. (standard 3-body angle gradient; Allen & Tildesley App.)
      const ga = 1 / la, gb = 1 / lb;
      const axy = ax * ga, ayy = ay * ga, azy = az * ga; // unit a
      const bxy = bx * gb, byy = by * gb, bzy = bz * gb; // unit b
      // force on i:  pref * ∂c/∂r_i = pref * (b̂ − c·â)/la
      let fix = pref * (bxy - c * axy) * ga;
      let fiy = pref * (byy - c * ayy) * ga;
      let fiz = pref * (bzy - c * azy) * ga;
      // force on k:  pref * (â − c·b̂)/lb
      let fkx = pref * (axy - c * bxy) * gb;
      let fky = pref * (ayy - c * byy) * gb;
      let fkz = pref * (azy - c * bzy) * gb;

      f[i] += fix; f[i + 1] += fiy; f[i + 2] += fiz;
      f[kk] += fkx; f[kk + 1] += fky; f[kk + 2] += fkz;
      f[j] -= fix + fkx; f[j + 1] -= fiy + fky; f[j + 2] -= fiz + fkz; // Newton's 3rd law
    }
    return U;
  }

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

  /**
   * Ligand bonds U = Σ ½ k (r−r0)². The list stores [i,j,r0] only — the force
   * constant is recovered from r0: aromatic ring bonds are clamped into the
   * resonance window r0 ≤ 1.44 Å ⇒ k = 200 kcal/mol/Å², all others (measured
   * r0 ≈ 1.5 Å) ⇒ k = 300 (see constants in ligand.js). Same harmonic kernel
   * as the protein bonds, but with per-bond k.
   */
  _ligandBondForces(pos, f) {
    let U = 0;
    const B = this.ligandBonds;
    for (let a = 0; a < B.length; a += 3) {
      const i = 3 * B[a], j = 3 * B[a + 1], r0 = B[a + 2];
      const k = r0 <= 1.44 ? 200 : 300;
      const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
      const dr = r - r0;
      U += 0.5 * k * dr * dr;
      // Force on j: −∂U/∂r_j = −k·dr·(r̂);  on i the opposite
      const s = (k * dr) / r;
      const fx = s * dx, fy = s * dy, fz = s * dz;
      f[i] += fx; f[i + 1] += fy; f[i + 2] += fz;
      f[j] -= fx; f[j + 1] -= fy; f[j + 2] -= fz;
    }
    return U;
  }

  /**
   * Improper (out-of-plane) term U = Σ ½ k (φ−φ0)² with k = 20 kcal/mol/rad²
   * for [i, j, k, l, φ0] (central atom j, planar φ0 = 0). Forces are taken by
   * central finite differences of φ (step h = 1e-5 Å) because φ is defined
   * through an absolute atan2 in ligand.improperAngle, whose analytic sign is
   * easy to get wrong; FD on the same function is self-consistent by
   * construction (verified against total-U finite differences in the tests).
   */
  _improperForces(pos, f, list) {
    const k = 20;
    const h = 1e-5;
    let U = 0;
    for (let a = 0; a < list.length; a += 5) {
      const i = list[a], j = list[a + 1], kk = list[a + 2], l = list[a + 3], phi0 = list[a + 4];
      const phi = improperAngle(pos, i, j, kk, l);
      U += 0.5 * k * (phi - phi0) * (phi - phi0);
      // −dU/dx_m = −k·(φ−φ0)·dφ/dx_m ; pos is restored after each sweep
      const pref = -k * (phi - phi0);
      for (const m of [i, j, kk, l]) {
        const c = 3 * m;
        for (let ax = 0; ax < 3; ax++) {
          const ci = c + ax, save = pos[ci];
          pos[ci] = save + h;
          const phiP = improperAngle(pos, i, j, kk, l);
          pos[ci] = save - h;
          const phiM = improperAngle(pos, i, j, kk, l);
          pos[ci] = save;
          f[ci] += pref * (phiP - phiM) / (2 * h);
        }
      }
    }
    return U;
  }

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
  _binding(pos, f) {
    this.bindingU = 0;
    this.desolvU = 0;
    if (!this.bindOn || this.nLigAtoms === 0) return 0;

    const nProt = this.nProt;
    const rc = this.bindRcut, rsw = 0.85 * rc;
    const rc2 = rc * rc, invDelta = 1 / (rc - rsw);
    const cell = rc;                  // cell = cutoff ⇒ 27-cell scan is complete
    const grid = this._gridB;

    // EEF1-lite constants: contact distance r0, Gaussian width σ, and the
    // density→burial scale nScale (n_a ≈ nScale ≈ 3 ⇒ ~63% buried).
    const R0 = 4.5, SIG = 1.8, NS = 3.0;
    const inv2sig2 = 1 / (2 * SIG * SIG); // 1/(2σ²) for g(r)
    const invSig2 = 1 / (SIG * SIG);      // 1/σ²   for dg/dr

    // Per-atom occupancy density + burial-derivative scratch (GC-free)
    const dens = this._dens, dBdn = this._dBdn;
    dens.fill(0);

    // Reused flat pair records for Pass 2 (a, j, r); capped at nProt·nLigAtoms.
    let bpN = 0;
    const bpA = this._bpA, bpJ = this._bpJ, bpR = this._bpR;

    // rebuild the protein-only grid in place (GC-free, like _repulsion)
    for (const arr of grid.values()) arr.length = 0;
    for (let i = 0; i < nProt; i++) {
      const key = this._cellKey(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2], cell);
      let arr = grid.get(key);
      if (!arr) grid.set(key, (arr = []));
      arr.push(i);
    }

    // =============================================================
    // PASS 1 — pair terms (LJ / electrostatics / H-bond) + density
    // =============================================================
    let U = 0;
    const EPSHB = 0.8, HB_R0 = 3.2, HB_W = 0.6, ELC = 332.0637;
    const w2 = HB_W * HB_W;

    for (let a = 0; a < this.nLigAtoms; a++) {
      const la = nProt + a, lx = 3 * la;
      const laX = pos[lx], laY = pos[lx + 1], laZ = pos[lx + 2];
      const ligSig = this._ligSigma[a], ligEps = this._ligEps[a];
      const ligQ = this._ligQ[a], ligHB = this._ligHB[a];
      const ckey = this._cellKey(laX, laY, laZ, cell);
      const cx = this._decodeX(ckey), cy = this._decodeY(ckey), cz = this._decodeZ(ckey);

      // 27 neighbouring cells (directed protein→ligand, no double counting)
      for (let ox = -1; ox <= 1; ox++)
        for (let oy = -1; oy <= 1; oy++)
          for (let oz = -1; oz <= 1; oz++) {
            const arr = grid.get(this._encodeCell(cx + ox, cy + oy, cz + oz));
            if (!arr) continue;
            for (let bi = 0; bi < arr.length; bi++) {
              const i = arr[bi], ix = 3 * i;
              const dx = pos[ix] - laX, dy = pos[ix + 1] - laY, dz = pos[ix + 2] - laZ;
              const r2 = dx * dx + dy * dy + dz * dz;
              if (r2 >= rc2 || r2 < 1e-10) continue;
              const r = Math.sqrt(r2);

              // Holo contact springs already govern this native pair — skip it
              // so the binding pair pass does not double-count the contact.
              const pk = this._pairKey(i, la);
              if (this._excluded.has(pk)) continue;

              // --- EEF1-lite occupancy density + pair record for Pass 2 ------
              // g(r) feeds the soft count dens[a]; the (a, j, r) triple lets
              // Pass 2 add the desolvation force without re-scanning the grid.
              const gr = Math.exp(-((r - R0) * (r - R0)) * inv2sig2);
              dens[a] += gr;
              bpA[bpN] = a; bpJ[bpN] = i; bpR[bpN] = r; bpN++;

              // smooth switch + its r-derivative
              let sw, dsw;
              if (r <= rsw) { sw = 1; dsw = 0; }
              else {
                const t = (r - rsw) * invDelta, t2 = t * t;
                sw = 1 - t2 * t * (10 - 15 * t + 6 * t2);
                dsw = -(30 * invDelta) * t2 * (1 - t) * (1 - t);
              }

              let dUdr = 0;   // total dU/dr of this pair (all cross terms)

              // --- cross LJ 12-6 (attractive well) ---------------------------
              const s = 0.5 * (this._protSigma[i] + ligSig);
              const e = Math.sqrt(this._protEps[i] * ligEps);
              const s6 = Math.pow(s / r, 6), s12 = s6 * s6;
              const phi = s12 - s6;                       // φ = (σ/r)¹² − (σ/r)⁶
              const dphi = (6 * s6 - 12 * s12) / r;       // dφ/dr
              U += 4 * e * phi * sw;
              dUdr += 4 * e * (dphi * sw + phi * dsw);

              // --- screened electrostatics ----------------------------------
              const q1 = this._protQ[i];
              if (q1 !== 0 && ligQ !== 0) {
                const ch = Math.cosh(r / 8);
                const epsr = 4 + 76 * Math.tanh(r / 8);
                const epsrP = 9.5 / (ch * ch);            // dεr/dr = 9.5·sech²(r/8)
                const g = 1 / (epsr * r);                 // g = 1/(εr·r)
                const gp = -(epsr + r * epsrP) / (epsr * epsr * r * r); // dg/dr
                const A = ELC * q1 * ligQ;
                U += A * g * sw;
                dUdr += A * (gp * sw + g * dsw);
              }

              // --- H-bond (donor/acceptor flags both set) -------------------
              if (this._protHB[i] && ligHB) {
                const g = Math.exp(-((r - HB_R0) * (r - HB_R0)) / (2 * w2));
                const gp = -g * (r - HB_R0) / w2;         // dg/dr
                const B = -EPSHB;                          // U = −epsHB·g·sw
                U += B * g * sw;
                dUdr += B * (gp * sw + g * dsw);
              }

              // force on protein i = −(dU/dr)·(dx/r); on the ligand, opposite
              const F = -dUdr / r;
              const fX = F * dx, fY = F * dy, fZ = F * dz;
              f[ix] += fX; f[ix + 1] += fY; f[ix + 2] += fZ;
              f[lx] -= fX; f[lx + 1] -= fY; f[lx + 2] -= fZ;
            }
          }
    }

    // =============================================================
    // PASS 2 — EEF1-lite burial: burial fraction, energy, forces
    // =============================================================
    // B_a = 1 − exp(−n_a/3); U_desolv = Σ_a ΔG_a·B_a. dB/dn is cached per atom
    // so every pair of atom a reuses the same factor (B depends on total n_a).
    let Udesolv = 0;
    for (let a = 0; a < this.nLigAtoms; a++) {
      const n = dens[a];
      const e = Math.exp(-n / NS);           // exp(−n_a/3)
      dBdn[a] = e / NS;                      // dB/dn = exp(−n_a/3)/3
      Udesolv += this._ligdG[a] * (1 - e);   // ΔG_a·B_a
    }

    // Chain rule dU_desolv/dr = ΔG_a·(dB/dn)·(dg/dr); force on protein j is
    // −(dU/dr)·(dx/r) with dx = pos_j − pos_la — same convention as Pass 1.
    for (let k = 0; k < bpN; k++) {
      const a = bpA[k], j = bpJ[k], r = bpR[k];
      const g = Math.exp(-((r - R0) * (r - R0)) * inv2sig2);
      const dgdr = -g * (r - R0) * invSig2;              // dg/dr
      const dUdr = this._ligdG[a] * dBdn[a] * dgdr;      // dU_desolv/dr
      const Fs = -dUdr / r;                              // −(dU/dr)/r, vector factor
      const la = nProt + a, lx = 3 * la, jx = 3 * j;
      const dx = pos[jx] - pos[lx], dy = pos[jx + 1] - pos[lx + 1], dz = pos[jx + 2] - pos[lx + 2];
      f[jx] += Fs * dx; f[jx + 1] += Fs * dy; f[jx + 2] += Fs * dz;
      f[lx] -= Fs * dx; f[lx + 1] -= Fs * dy; f[lx + 2] -= Fs * dz;
    }

    this.bindingU = U + Udesolv;
    this.desolvU = Udesolv;
    return U + Udesolv;
  }

  /**
   * Repulsive-only 12-6 LJ for all non-bonded/non-contact pairs (implicit
   * solvent: beads cannot overlap). Spatial hash grid with numeric keys.
   *   U_rep(r) = ε[(r_e/r)¹² − 2(r_e/r)⁶ + 1],   r < r_e = 2^{1/6}σ
   *   F(r)     = 24ε[2(σ/r)¹²·σ? − (σ/r)⁶]/r² — see algebra in code body.
   */
  _repulsion(pos, f) {
    const n = this.n, cell = this._cell, rc = this.rcRep, eps = this.epsRep;
    const grid = this._grid, sig = this._repSigma;
    // clear lists in place (avoid GC churn)
    for (const arr of grid.values()) arr.length = 0;

    // insert beads into grid cells
    for (let i = 0; i < n; i++) {
      const key = this._cellKey(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2], cell);
      let arr = grid.get(key);
      if (!arr) grid.set(key, (arr = []));
      arr.push(i);
    }

    let U = 0;
    // iterate over cells; for each, self-pairs handled separately
    for (const [key, arr] of grid) {
      const cx = this._decodeX(key), cy = this._decodeY(key), cz = this._decodeZ(key);
      // visit 14 of 27 neighbours (half-shell) to count each pair once
      for (let ox = -1; ox <= 1; ox++)
        for (let oy = -1; oy <= 1; oy++)
          for (let oz = -1; oz <= 1; oz++) {
            if (ox < 0 || (ox === 0 && oy < 0) || (ox === 0 && oy === 0 && oz < 0)) continue;
            const same = ox === 0 && oy === 0 && oz === 0;
            const nbr = same ? arr : grid.get(this._encodeCell(cx + ox, cy + oy, cz + oz));
            if (!nbr) continue;
            for (let ai = 0; ai < arr.length; ai++) {
              const i = arr[ai], ix = 3 * i;
              for (let bi = same ? ai + 1 : 0; bi < nbr.length; bi++) {
                const j = nbr[bi];
                const pk = this._pairKey(i, j);
                if ((i < this.nProt) !== (j < this.nProt)) continue;
                if (this._excluded.has(pk)) continue;
                const jx = 3 * j;
                const dx = pos[jx] - pos[ix], dy = pos[jx + 1] - pos[ix + 1], dz = pos[jx + 2] - pos[ix + 2];
                const r2 = dx * dx + dy * dy + dz * dz;
                // per-pair σ = arithmetic mean of the two bead sizes; its own
                // r_e = 2^(1/6)σ per pair (keeps U(contact edge) = 0 exactly)
                const s = 0.5 * (sig[i] + sig[j]);
                const re = rc * (s / this.sigmaRep);
                const re2 = re * re;
                if (r2 >= re2 || r2 < 1e-10) continue;
                const r = Math.sqrt(r2);
                const sr = s / r;
                const sr2 = sr * sr;
                const sr6 = sr2 * sr2 * sr2;              // (σ/r)^6
                const uLJ = sr6 * sr6 - sr6;              // [(σ/r)¹² − (σ/r)⁶]
                // energy: ε·(4·uLJ + 1)   (shifted so U(rcRep) = 0 exactly)
                U += eps * (4 * uLJ + 1);
                // −dU/dr = 24ε[2(σ/r)¹² − (σ/r)⁶]/r ; vector form divides by r²
                const fm = (24 * eps * (2 * sr6 * sr6 - sr6)) / r2;
                // i feels −∇_i U → away from j (repulsive out of overlap)
                f[ix] -= fm * dx; f[ix + 1] -= fm * dy; f[ix + 2] -= fm * dz;
                f[jx] += fm * dx; f[jx + 1] += fm * dy; f[jx + 2] += fm * dz;
              }
            }
          }
    }
    return U;
  }

  /* --- numeric spatial-hash cell keys (no string allocs) -------------- */
  _encodeCell(cx, cy, cz) { return ((cx + 2048) * 4096 + (cy + 2048)) * 4096 + (cz + 2048); }
  _cellKey(x, y, z, cell) {
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
