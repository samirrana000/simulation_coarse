/**
 * cg/topology.js — CG connectivity and native-contact topology.
 * Split out of src/forcefield.js; moved verbatim.
 *
 * WHAT THIS FILE OWNS
 * -------------------
 * Everything about WHICH PAIRS ARE BONDED, built once at construction and
 * never touched by the integrator:
 *   1. 1-2 peptide bonds and 1-2-3 pseudo-angles from the segment table,
 *      with r0/θ0 measured off the reference coordinates;
 *   2. the uniform-γ ENM contact set (H(Rc − r0), 1-2/1-3 excluded);
 *   3. the non-repulsive `_excluded` key set (bonds, angles, ENM, holo and
 *      covalent ligand pairs incl. the full 1-4 walk);
 *   4. per-particle repulsive σ;
 *   5. the Gö-like `nativeContacts` springs for native pairs that sit inside
 *      their own repulsion range.
 *
 * WHY IT IS ITS OWN MODULE
 * ------------------------
 * This is the double-counting ledger. docs/EXCLUSIONS.md documents the rule
 * that every pair is seen by exactly ONE term, and this file is the whole of
 * that rule; the exclusion set, the native-contact scan and the ENM cutoff
 * are three views of one decision. Splitting it out of the constructor means
 * a future term has to be added in ONE file next to the ledger that must be
 * updated with it, instead of inside a 390-line constructor.
 *
 * ORDERING (load-bearing): this must run LAST in the constructor. It reads
 * ff.ref / ff.rc / ff._dist / ff._pairKey for every equilibrium distance, it
 * calls ff.rebuildHoloSprings() before _excluded exists (that method is
 * guarded), and allocBuffers() at the very end must see the final n and the
 * final ligand count. The steps are called in the original order and each
 * loop is a straight per-pair scan, so no floating-point accumulation order
 * is disturbed by the split.
 */

/**
 * Build the whole CG topology on `ff`, in the original construction order.
 * @param {object} ff        the ForceField under construction
 * @param {object} par       constructor params (rc, gamma, seqWeight, enmModel, tirion)
 * @param {Array}  segments  pdb.selectSystem() segments [[start, end), ...]
 */
export function buildTopology(ff, par, segments) {
  buildConnectivity(ff, segments);
  buildEnmContacts(ff, par, segments);

  // ---- holo contact springs ----------------------------------------------
  // Harmonic protein–ligand springs for every pair whose NATIVE distance
  // r0 <= 6.0 Å, encoding the observed bound pose of a holo complex (e.g. T4
  // lysozyme L99A + benzene, PDB 4W52): pins the ligand into its binding
  // pocket while still allowing fluctuation. Gated by par.binding.holo
  // (default true); no-op when no ligands are present.
  ff.holoSprings = new Float64Array(0);
  ff.nHolo = 0;
  ff.rebuildHoloSprings();

  buildExcludedPairs(ff);
  buildRepulsionSigma(ff);
  buildNativeContacts(ff);
}

/**
 * ---- 1-2 bonds & 1-2-3 angles from segments --------------------------
 * @param {object} ff        the ForceField under construction
 * @param {Array}  segments  [[start, end), ...]
 */
function buildConnectivity(ff, segments) {
  /** bonds[idx] = [i, j, r0] packed flat as [i0,j0,r0, i1,j1,r1, ...] */
  const bonds = [];
  const angles = []; // [i, j, k, theta0]
  for (const [s, e] of segments) {
    for (let i = s; i + 1 < e; i++) {
      bonds.push(i, i + 1, ff._dist(ff.ref, i, i + 1));
    }
    for (let i = s; i + 2 < e; i++) {
      angles.push(i, i + 1, i + 2, ff._angle(ff.ref, i, i + 1, i + 2));
    }
  }
  ff.bonds = new Float64Array(bonds);
  ff.angles = new Float64Array(angles);
}

/**
 * ---- ENM contact set (also used for drawing) --------------------------
 * Pair (i,j) gets a spring if r0(i,j) <= rc and |i-j| > 2 within the same
 * segment-family OR i,j belong to different segments/chains (interfacial
 * contacts in complexes matter and must be elastic too).
 * Protein pairs only — ligand geometry is frozen by its own internal FF,
 * and protein–ligand elastic contacts are left to the excluded volume.
 * @param {object} ff        the ForceField under construction
 * @param {object} par       constructor params
 * @param {Array}  segments  [[start, end), ...]
 */
function buildEnmContacts(ff, par, segments) {
  const springs = [];      // [i, j, r0]
  const isBound13 = (i, j) => {
    for (const [s, e] of segments) {
      if (i >= s && j < e) return j - i <= 2; // 1-2 and 1-3 excluded
      if (j >= s && i < e) return i - j <= 2;
    }
    return false;
  };
  for (let i = 0; i < ff.nProt; i++) {
    for (let j = i + 1; j < ff.nProt; j++) {
      const r0 = ff._dist(ff.ref, i, j);
      if (r0 <= ff.rc && !isBound13(i, j)) {
        springs.push(i, j, r0);
      }
    }
  }
  ff.springs = new Float64Array(springs);
  // Per-spring force constants (kcal/mol/Å²). Starts uniform at ff.gamma;
  // an ML contact prior (contacts.json from ml/export_esm_contacts.py) can
  // rescale individual pairs via setSpringScale() — e.g. stiffen residues
  // the model predicts to be in contact.
  ff.springK = new Float64Array(ff.springs.length / 3).fill(ff.gamma);
  ff.springScaleActive = false;
  // ── Sequence-dependent ENM (Bahar-style; opt-in, applied by springs.js) ──
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
  if (par.seqWeight) ff.applySeqWeights();
  // ---- Tirion distance-weighted ENM (opt-in, Phase 1) --------------------
  // par.enmModel === "tirion" (or par.tirion === true) replaces the uniform
  // springK with γ_ij = γ0·(R0/r0_ij)^6 plus SS-dependent backbone basins.
  // Default ("uniform") preserves legacy behavior so existing tests pass.
  // See src/physics/forcefield/tirion_anm.js.
  ff.enmModel = par.enmModel ?? (par.tirion ? "tirion" : "uniform");
  ff.tirionDihedrals = new Float64Array(0);
  ff.tirionDihedralK = new Float64Array(0);
  ff.tirionSS = [];
  if (ff.enmModel === "tirion") {
    try {
      ff.useTirionNetwork({ gamma0: ff.gamma, cutoff: ff.rc, segments });
    } catch (e) {
      console.warn(`[ForceField] Tirion network failed (${e.message}) — falling back to uniform ENM`);
      ff.enmModel = "uniform";
    }
  }
}

/**
 * ---- repulsive-pair bookkeeping ---------------------------------------
 * Non-repulsive pairs (bonded 1-2, 1-3 and ENM springs) are excluded to
 * avoid double counting; the grid only tests the remainder.
 * @param {object} ff  the ForceField under construction
 */
function buildExcludedPairs(ff) {
  ff._excluded = new Set();
  for (let k = 0; k < ff.bonds.length; k += 3) {
    ff._excluded.add(ff._pairKey(ff.bonds[k], ff.bonds[k + 1]));
  }
  for (let k = 0; k < ff.angles.length; k += 4) {
    ff._excluded.add(ff._pairKey(ff.angles[k], ff.angles[k + 2]));
  }
  for (let k = 0; k < ff.springs.length; k += 3) {
    ff._excluded.add(ff._pairKey(ff.springs[k], ff.springs[k + 1]));
  }
  // Holo contact springs govern native protein–ligand contact distances, so
  // those pairs are excluded from the binding pair pass to avoid double
  // counting (the spring alone holds them at r0).
  for (let k = 0; k < ff.holoSprings.length; k += 3) {
    ff._excluded.add(ff._pairKey(ff.holoSprings[k], ff.holoSprings[k + 1]));
  }
  // Covalently linked ligand pairs (1-2 bonds, 1-3 angles, 1-4 improper
  // partners) must never be repelled by the grid: the internal FF already
  // governs those distances and they sit far inside the repulsion range
  // (e.g. para carbons of an aromatic ring at ~2.8 Å ≪ r_e ≈ 5.0 Å).
  for (let a = 0; a < ff.ligandBonds.length; a += 3) {
    ff._excluded.add(ff._pairKey(ff.ligandBonds[a], ff.ligandBonds[a + 1]));
  }
  for (let a = 0; a < ff.ligandAngles.length; a += 4) {
    ff._excluded.add(ff._pairKey(ff.ligandAngles[a], ff.ligandAngles[a + 2]));
  }
  for (let a = 0; a < ff.ligandImpropers.length; a += 5) {
    const i = ff.ligandImpropers[a], j = ff.ligandImpropers[a + 1];
    const l = ff.ligandImpropers[a + 3];
    ff._excluded.add(ff._pairKey(i, l));
    ff._excluded.add(ff._pairKey(j, l));
    ff._excluded.add(ff._pairKey(ff.ligandImpropers[a + 2], l));
  }
  // Ligand 1-4 pairs (three bonds apart, e.g. the N…N of a N–C–C–N fragment)
  // are NOT repelled by the grid: gauche 1-4 conformers sit at ~2.9 Å — far
  // inside the Cα-sized repulsion range — and with no explicit non-ring
  // torsion term they are legitimately free. Only ring 1-4 pairs were
  // already covered by the improper exclusions above; this walk covers every
  // 1-4 pair in any molecule (incl. poly-atomic buffers like HEPES).
  if (ff.ligandBonds.length) {      const adjLig = new Map();
    for (let a = 0; a < ff.ligandBonds.length; a += 3) {
      const i = ff.ligandBonds[a], j = ff.ligandBonds[a + 1];
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
            if (c !== x && c !== a) ff._excluded.add(ff._pairKey(x, c));
          }
        }
      }
    }
  }
}

/**
 * Per-particle repulsive σ: protein beads use the Cα size; ligand atoms use
 * their element's LJ σ (binding tables) so intra-ligand excluded volume is
 * not over-estimated — a folded 15-atom co-solute must not sit at 10 kcal.
 * Allocated here (before the native-contact scan, which reads the table).
 * @param {object} ff  the ForceField under construction
 */
function buildRepulsionSigma(ff) {
  ff._repSigma = new Float64Array(ff.n);
  for (let i = 0; i < ff.nProt; i++) ff._repSigma[i] = ff.sigmaRep;
  for (let a = 0; a < ff.nLigAtoms; a++) ff._repSigma[ff.nProt + a] = ff._ligSigma[a];
}

/**
 * ---- Gö-like native contacts (ligand 1-5+ pairs) -----------------------
 * Folded ligands (sugars, inhibitors, buffers) place non-bonded 1-5+
 * pairs at 2.9–4.0 Å — inside their own LJ repulsion range r_e. Without
 * treatment the grid sees a huge native strain (PDBBind probes: U(native)
 * up to ~745 kcal/mol, max|F| ~140 kcal/mol/Å) that launches the light
 * united atoms as projectiles and blasts the protein. Standard CG fix:
 * every non-excluded pair whose NATIVE distance lies inside its own
 * repulsion range is excluded from the grid and gets a soft Gö-like
 * harmonic contact spring (k = 1.0 kcal/mol/Å², [i,j,r0] packing like
 * holoSprings) that preserves the native distance while still resisting
 * compression. Applied to ANY particle type (L–L, P–P, P–L) — a genuine
 * native clash is a modeling problem regardless of chemistry.
 *
 * Double-coverage diagram (see docs/CG_HEAVY.md: nativeContacts):
 *   native r0 ──►  [excluded from grid repulsion?] ──► [already ENM spring?]
 *         │                     │                               │
 *         │  r0 < r_e && not in _excluded                yes → keep ENM only
 *         │         │                               no → create nativeContacts spring
 *         │         └─► yes → push (i,j,r0) to nativeContacts, add to _excluded
 *         │                     (prevents 12-6 repulsion + harmonic double count)
 *         └─► r0 ≥ r_e → leave to grid repulsion (smooth WCA, zero at r_e)
 *   Intra-molecule only (P–P, L–L): cross P–L pairs are holo-pinned or
 *   binding-evaluated, never grid-repelled, so no Gö spring there (would
 *   double-count attractive binding). The scan is O(N²) at construction
 *   only; compute() then sees each pair exactly once (CG_HEAVY.md Fig.1).
 * @param {object} ff  the ForceField under construction
 */
function buildNativeContacts(ff) {
  ff.nativeContacts = new Float64Array(0); // docs/CG_HEAVY.md: nativeContacts double-coverage
  {
    const nx = [];
    for (let i = 0; i < ff.n; i++) {
      for (let j = i + 1; j < ff.n; j++) {
        // Intra-molecule scan only: cross protein–ligand pairs never reach
        // the grid repulsion (they are either holo-spring pinned — excluded
        // from binding AND grid — or evaluated by the binding pass with
        // attractive wells), so adding a Gö-like spring there would double-
        // count the interaction.
        if ((i < ff.nProt) !== (j < ff.nProt)) continue;
        if (ff._excluded.has(ff._pairKey(i, j))) continue;
        const r0 = ff._dist(ff.ref, i, j);
        const s = 0.5 * (ff._repSigma[i] + ff._repSigma[j]);
        const re = (2 ** (1 / 6)) * s;
        if (r0 < re) nx.push(i, j, r0);
      }
    }
    if (nx.length) {
      ff.nativeContacts = new Float64Array(nx);
      for (let k = 0; k < nx.length; k += 3) {
        ff._excluded.add(ff._pairKey(nx[k], nx[k + 1]));
      }
    }
  }
}