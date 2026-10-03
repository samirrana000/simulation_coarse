/**
 * heavy/weak.js — the opt-in weak-interaction layer: p-stacking, cation-pi and
 * halogen sigma-hole, plus the once-per-topology pair lists and the metal
 * coordination upgrade that feeds them.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10). Attached to
 * HeavyForceField.prototype by heavy/forcefield.js.
 *
 *   _buildWeakAndMetalLists  once per topology: ring frames, cations, halogens,
 *                            acceptors; and the chem/metals.js coordination
 *                            detection that upgrades metal distance springs into
 *                            radial+angular restraints (par.metalAngles)
 *   _weakInteractions        the per-step pass. All three kernels must be called
 *                            INSIDE this method and the result must reach U —
 *                            scripts/validate_binding_physics_r1.mjs asserts
 *                            exactly that, on the comment-stripped source, so a
 *                            kernel cannot be "implemented" in prose only.
 *
 * Zero DOM globals. Node-importable.
 */
import {
  piStackForces, cationPiForces, halogenForces,
  buildRingFrames, buildCationList, buildHalogenList, HALOGEN_EPS,
} from "../physics/weakint.js";
import { detectCoordination } from "../chem/metals.js";
import { pairKey } from "./pairs.js";

/** Attached to HeavyForceField.prototype by heavy/forcefield.js. */
export const weakKernels = {
  /**
   * Build the once-per-topology weak-interaction lists + metal coordination
   * upgrade state (Loop-2 S3).
   * Ring frames from buildRingFrames (protein name sets + GAFF2 types +
   * geometric fallback); cation list from Lys NZ / Arg CZ / HIP / charged
   * ligand N; halogen list from Cl/Br/I with a bonded C; metal upgrade from
   * detectCoordination + classifyGeometry at the construction pose.
   * @param {Array} atoms
   * @param {Array} bonds  flat [i,j,r0,...] list from buildLists
   */
  _buildWeakAndMetalLists(atoms, bonds) {
    const bondPairs = [];
    for (let a = 0; a < bonds.length; a += 3) bondPairs.push([bonds[a], bonds[a + 1]]);
    if (this.weakOn) {
      this._weakRings = buildRingFrames(atoms, bondPairs);
      this._weakRingAtoms = new Set();
      for (const r of this._weakRings) for (const i of r.atomIdx) this._weakRingAtoms.add(i);
      this._weakCations = buildCationList(atoms);
      this._weakHalogens = buildHalogenList(atoms, bondPairs);
      // Acceptor list for halogen bonds: hbond acceptor classification
      // (backbone/sidechain O, S, aromatic N) — same set as the H-bond pass.
      const acc = this._hbClassification.isAcceptor;
      this._weakAcceptors = [];
      for (let i = 0; i < this.n; i++) if (acc[i]) this._weakAcceptors.push(i);
    }
    if (this.metalAngles) {
      // R3 §1e metal upgrade: detect coordination at the native pose, classify
      // the polyhedron, and pin ideal angles. Metals without a detected
      // geometry (coordinationNumber ≤ 1) keep the legacy k=40 springs.
      const elements = atoms.map((a) => a.element ?? "C");
      const metals = [];
      const hasGeometry = new Set();
      for (let i = 0; i < this.n; i++) {
        if (!atoms[i].isMetal) continue;
        const det = detectCoordination(this.ref, i, elements);
        metals.push({ index: i, element: atoms[i].element, donors: det.donorIndices });
        if (det.coordinationNumber >= 2) hasGeometry.add(i);
      }
      this._metalEnforce = { metals, elements, hasGeometry };
      // Remove those metals' radial springs from the legacy coord list so
      // enforceCoordination is the ONLY radial term for them (no double radial).
      if (hasGeometry.size > 0 && this.coord.length > 0) {
        const keep = [];
        for (let a = 0; a < this.coord.length; a += 3) {
          if (hasGeometry.has(this.coord[a])) continue; // metal side of [metal, donor]
          keep.push(this.coord[a], this.coord[a + 1], this.coord[a + 2]);
        }
        this.coord = new Float64Array(keep);
      }
    }
  },

  /**
   * Weak-interaction pass (π-stack, cation-π, halogen σ-hole) — opt-in,
   * runs after the LJ/GB grid kernel. Cutoffs: ring-ring 5.5 Å (centroid),
   * cation-ring 6 Å, X···D 4 Å (all inside the kernels' Gaussian tails;
   * pre-screened by centroid distance to skip far pairs cheaply).
   * Exclusions: pairs present in _excluded (bonded/1-3/intra-ligand) and
   * same-ring pairs are skipped. Energies returned per term.
   * @param {ArrayLike<number>} pos
   * @param {Float64Array} f
   * @returns {{pi:number, cpi:number, xb:number}}
   */
  _weakInteractions(pos, f) {
    const rings = this._weakRings;
    const excl = this._excluded;
    let pi = 0, cpi = 0, xb = 0;

    // --- π-stack: ring-ring pairs within 5.5 Å centroid cutoff ---
    for (let a = 0; a < rings.length; a++) {
      const ra = rings[a];
      const a0 = 3 * ra.atomIdx[0];
      const ax = pos[a0], ay = pos[a0 + 1], az = pos[a0 + 2];
      for (let b = a + 1; b < rings.length; b++) {
        const rb = rings[b];
        // cheap prescreen on first atom (within ring diameter of centroid)
        const b0 = 3 * rb.atomIdx[0];
        const dx = pos[b0] - ax, dy = pos[b0 + 1] - ay, dz = pos[b0 + 2] - az;
        if (dx * dx + dy * dy + dz * dz > 121) continue; // 11 Å atom prescreen ≫ 5.5 + 2×2.8 ring radius
        // fused rings share a bond (e.g. Trp 5+6 rings) — those pairs sit in
        // _excluded and must not double-stack; same-molecule NON-bonded rings
        // (Phe–Phe′ stacking) keep the term (R3 §4). Intra-ligand pairs are
        // all-excluded → ligand-internal stacking off (geometry already fixed).
        let skip = false;
        for (const ia of ra.atomIdx) {
          for (const ib of rb.atomIdx) {
            if (ia === ib || excl.has(pairKey(ia, ib))) { skip = true; break; }
          }
          if (skip) break;
        }
        if (skip) continue;
        pi += piStackForces(pos, f, ra, rb);
      }
    }

    // --- cation-π: cation-ring pairs within 6 Å ---
    for (const ci of this._weakCations) {
      const cx = pos[3 * ci], cy = pos[3 * ci + 1], cz = pos[3 * ci + 2];
      for (const ring of rings) {
        // cations inside their own ring (pyridinium N) skip
        if (ring.atomIdx.includes(ci)) continue;
        const i0 = 3 * ring.atomIdx[0];
        const dx = pos[i0] - cx, dy = pos[i0 + 1] - cy, dz = pos[i0 + 2] - cz;
        if (dx * dx + dy * dy + dz * dz > 100) continue; // 10 Å prescreen
        // cation covalently tied to the ring (aniline-type N) — excluded pair
        let skip = false;
        for (const ia of ring.atomIdx) {
          if (excl.has(pairKey(ci, ia))) { skip = true; break; }
        }
        if (skip) continue;
        cpi += cationPiForces(pos, f, ci, ring);
      }
    }

    // --- halogen σ-hole: C–X···D triples, X···D within 4 Å ---
    for (const { x, c } of this._weakHalogens) {
      const xx = pos[3 * x], xy = pos[3 * x + 1], xz = pos[3 * x + 2];
      const el = String(this.atoms[x]?.element ?? "").toUpperCase();
      const epsX = HALOGEN_EPS[el] ?? 1.2;
      for (const d of this._weakAcceptors) {
        if (d === x || d === c) continue;
        const dx = pos[3 * d] - xx, dy = pos[3 * d + 1] - xy, dz = pos[3 * d + 2] - xz;
        const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > 16) continue; // 4 Å
        // bonded X···D / C···D pairs skipped (intra-ligand included)
        if (excl.has(pairKey(x, d)) || excl.has(pairKey(c, d))) continue;
        xb += halogenForces(pos, f, c, x, d, { eps: epsX });
      }
    }

    return { pi, cpi, xb };
  },
};
