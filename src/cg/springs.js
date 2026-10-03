/**
 * cg/springs.js — mutation kernels for the CG harmonic spring networks.
 * Split out of src/forcefield.js; moved verbatim.
 *
 * WHAT THIS FILE OWNS
 * -------------------
 * Every operation that REWRITES a flat spring list or its stiffness table
 * after construction:
 *   - rebuildHoloSprings  — native protein–ligand pose pins (6 Å cutoff)
 *   - setSpringScale      — ML contact prior, K = γ·(1 + α·p)
 *   - clearSpringScale    — back to uniform γ
 *   - applySeqWeights     — Bahar-style sequence-dependent K_seq
 *   - useTirionNetwork    — distance-weighted γ_ij = γ0·(R0/r0)^6 + SS basins
 *
 * WHY THEY SHARE A MODULE
 * -----------------------
 * All five do the same two things: rebuild a spring list (or K table), then
 * keep `ff._excluded` consistent with it. That second half is the part that is
 * easy to get wrong — an ENM rebuild that forgets to drop the old keys leaves
 * the repulsion grid double-counting against a spring that no longer exists.
 * Sharing a module makes the invariant "whatever you rebuild, re-sync
 * _excluded" a property of the file rather than a discipline of the author.
 *
 * These are all MUTATORS on a constructed ff. None of them allocates in the
 * integrator hot path; they are called from the constructor, from placement,
 * and from the ML/UI layers.
 */

import { RES_CLASS_OF, SEQ_WEIGHT } from "../ff-params.js";
import { buildTirionNetwork, applyTirionToForceField } from "../physics/forcefield/tirion_anm.js";

/**
 * Rebuild the holo (native-pose) protein–ligand springs from the CURRENT
 * reference coordinates. In the constructor this derives the initial holo
 * springs from the native reference; after a placement it is called with
 * holoOn=false to clear them (a placed library ligand is a hypothesis, not
 * a known crystallographic pose, so it is not pinned). Also keeps the
 * non-repulsive exclusion set in sync with whatever holo pairs exist.
 * No-op when holo springs are off or there are no ligand atoms.
 * @param {object} ff  the ForceField
 */
export function rebuildHoloSprings(ff) {
  // remove the old holo pairs from the non-repulsive exclusion set
  // (guarded: the constructor calls this before _excluded is created)
  if (ff._excluded) {
    for (let k = 0; k < ff.holoSprings.length; k += 3) {
      ff._excluded.delete(ff._pairKey(ff.holoSprings[k], ff.holoSprings[k + 1]));
    }
  }
  ff.holoSprings = new Float64Array(0);
  ff.nHolo = 0;
  if (!ff.holoOn || ff.nLigAtoms <= 0) return;
  const hs = [];
  for (let i = 0; i < ff.nProt; i++) {
    for (let la = ff.nProt; la < ff.n; la++) {
      const r0 = ff._dist(ff.ref, i, la);
      if (r0 <= 6.0) hs.push(i, la, r0);
    }
  }
  ff.holoSprings = new Float64Array(hs);
  ff.nHolo = hs.length / 3;
  if (ff._excluded) {
    for (let k = 0; k < ff.holoSprings.length; k += 3) {
      ff._excluded.add(ff._pairKey(ff.holoSprings[k], ff.holoSprings[k + 1]));
    }
  }
}

/**
 * Apply an ML contact prior to the ENM spring constants.
 * @param {object} ff         the ForceField
 * @param {Array}  contacts   [[i, j, p], ...] residue pairs (0-indexed) with p ∈ [0,1].
 * @param {number} alpha      global scale (k_per_pair = gamma · (1 + alpha·p)).
 * Marks springScaleActive so the UI knows a custom map is live.
 */
export function setSpringScale(ff, contacts, alpha = 1) {
  const S = ff.springs, K = ff.springK, gamma = ff.gamma;
  // reset to baseline gamma, then rescale by the prior
  K.fill(gamma);
  const pairScale = new Map();
  for (const [i, j, p] of contacts) {
    if (p <= 0) continue;
    pairScale.set(ff._pairKey(i, j), 1 + alpha * p);
  }
  for (let k = 0, s = 0; k < S.length; k += 3, s++) {
    const sc = pairScale.get(ff._pairKey(S[k], S[k + 1]));
    if (sc !== undefined) K[s] = gamma * sc;
  }
  ff.springScaleActive = true;
}

/**
 * Restore uniform ENM stiffness (undo an ML contact-prior scaling).
 * @param {object} ff  the ForceField
 */
export function clearSpringScale(ff) {
  ff.springK.fill(ff.gamma);
  ff.springScaleActive = false;
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
 * @param {object} ff      the ForceField
 * @param {Array}  beads   optional override (defaults to constructor beads if stored)
 */
export function applySeqWeights(ff, beads = null) {
  // beads not stored on ForceField; caller may pass original bead list,
  // otherwise we fall back to resClass table which was built from beads.
  // resClass is Uint8Array over 5 classes; invert to weight via SEQ_WEIGHT.
  const weightByClass = [SEQ_WEIGHT.H, SEQ_WEIGHT.A, SEQ_WEIGHT.P, SEQ_WEIGHT.Cp, SEQ_WEIGHT.Cn];
  // w_i per bead from ff.resClass (built in constructor)
  const nProt = ff.nProt;
  const wPerBead = new Float64Array(nProt);
  for (let i = 0; i < nProt; i++) {
    // resClass[i] is 0:H,1:A,2:P,3:Cp,4:Cn (see constructor)
    wPerBead[i] = weightByClass[ff.resClass[i]] ?? 1.0;
  }
  // If caller supplied explicit beads (e.g. for test), prefer SEQ_WEIGHT via resName
  if (beads && beads.length === nProt) {
    for (let i = 0; i < nProt; i++) {
      const cls = RES_CLASS_OF[beads[i].resName] ?? "H";
      const w = SEQ_WEIGHT[cls];
      if (w !== undefined) wPerBead[i] = w;
    }
  }
  const S = ff.springs, K = ff.springK, gamma = ff.gamma;
  for (let k = 0, s = 0; k < S.length; k += 3, s++) {
    const i = S[k], j = S[k + 1];
    const w_i = wPerBead[i], w_j = wPerBead[j];
    // K_seq = gamma * (1 + 0.2 * (w_i + w_j)/2)  — Bahar-style ±2% modulation
    K[s] = gamma * (1 + 0.2 * (w_i + w_j) * 0.5);
  }
  ff.springScaleActive = true; // marks non-uniform (reuse flag for UI)
  // For debugging: ff._seqW = wPerBead; // keep per-bead w_i if needed
}

/**
 * Opt-in Tirion distance-weighted ENM + SS dihedral basins (Phase 1).
 * Rebuilds springs/springK as γ_ij = γ0·(R0/r0_ij)^6 over the same
 * H(Rc−r0) topology and stores backbone pseudo-dihedral basins on
 * ff.tirionDihedrals / ff.tirionDihedralK / ff.tirionSS. Falls back to the
 * existing uniform network when the module is unavailable.
 * @param {object} ff    the ForceField
 * @param {object} [opts] { gamma0, R0, cutoff, segments }
 */
export function useTirionNetwork(ff, opts = {}) {
  const ref = ff.ref.subarray(0, ff.nProt * 3);
  // Reconstruct segments if the caller did not pass them: derive from the
  // current bond list (bonds only join contiguous pairs).
  let segments = opts.segments ?? null;
  if (!segments) {
    segments = [];
    if (ff.nProt > 0) {
      let s = 0;
      const bondedNext = new Set();
      for (let k = 0; k < ff.bonds.length; k += 3) {
        const i = ff.bonds[k], j = ff.bonds[k + 1];
        if (j === i + 1) bondedNext.add(i);
      }
      for (let i = 0; i < ff.nProt - 1; i++) {
        if (!bondedNext.has(i)) { segments.push([s, i + 1]); s = i + 1; }
      }
      segments.push([s, ff.nProt]);
    }
  }
  const net = buildTirionNetwork(ref, {
    gamma0: opts.gamma0 ?? ff.gamma,
    R0: opts.R0 ?? 3.81,
    cutoff: opts.cutoff ?? ff.rc,
    segments,
  });
  // Replace uniform springs; keep exclusion set consistent.
  if (ff._excluded) {
    for (let k = 0; k < ff.springs.length; k += 3) {
      ff._excluded.delete(ff._pairKey(ff.springs[k], ff.springs[k + 1]));
    }
  }
  applyTirionToForceField(ff, net);
  ff.enmModel = "tirion";
}