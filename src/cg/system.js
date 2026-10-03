/**
 * cg/system.js — CG particle-set assembly: reference coordinates, the
 * appended united-atom ligand force field, masses, and the rmsd fold mask.
 * Split out of src/forcefield.js; moved verbatim.
 *
 * WHAT THIS FILE OWNS
 * -------------------
 * The one decision that defines the CG model's index space: the Cα beads take
 * global indices 0..nProt−1 and every ligand heavy atom is APPENDED at
 * nProt..n−1. `this.n`, `this.ref`, `this.masses`, `this.isLigand` and
 * `this._maskProt` all follow from that choice, and so does every kernel's
 * assumption that `a < nProt` means protein. Nothing else in the CG engine
 * re-derives the layout.
 *
 * WHY IT IS ITS OWN MODULE
 * ------------------------
 * This is the only part of construction that can grow the particle set, and
 * it is the only part that reads ligand.js. Keeping it apart means the
 * "protein-only ⇒ exactly the legacy model" property is readable in one
 * screen: with no ligands, this function sets n = nProt and returns.
 *
 * ORDERING: must run after initScalars() and before assignParticleParams()
 * (which sizes the ligand tables from nLigAtoms) and before
 * buildTopology() (which reads ref for every equilibrium distance).
 */

import { buildLigandInternalFF } from "../ligand.js";

/**
 * Build the unified multi-particle set on `ff`.
 * @param {object} ff       the ForceField under construction
 * @param {Array}  beads    pdb.selectSystem() beads (protein Cα only)
 * @param {object} par      constructor params (mass)
 * @param {Array}  ligands  pdb.parseLigands() output; empty ⇒ protein-only
 */
export function buildParticles(ff, beads, par, ligands) {
  const nProt = ff.nProt;

  // Reference coordinates (native state): copied, never mutated.
  ff.ref = new Float64Array(nProt * 3);
  beads.forEach((b, i) => {
    ff.ref[3 * i] = b.x; ff.ref[3 * i + 1] = b.y; ff.ref[3 * i + 2] = b.z;
  });

  // ---- united-atom ligand force field (appended particles) --------------
  // Ligand heavy atoms become explicit particles with global indices
  // nProt..n−1. buildLigandInternalFF already emits ABSOLUTE indices and
  // united-atom masses, so its flat arrays are consumed directly. isLigand
  // flags the appended particles for the analysis helpers.
  ff.isLigand = new Uint8Array(nProt);
  ff.ligandAtoms = [];
  ff.ligandMasses = new Float64Array(0);
  ff.ligandBonds = new Float64Array(0);
  ff.ligandAngles = new Float64Array(0);
  ff.ligandImpropers = new Float64Array(0);
  if (ligands.length) {
    const lig = buildLigandInternalFF(ligands, nProt);
    ff.ligandAtoms = lig.atoms;
    ff.ligandMasses = lig.masses;
    ff.ligandBonds = lig.bonds;
    ff.ligandAngles = lig.angles;
    ff.ligandImpropers = lig.impropers;

    // Append ligand coordinates in molecule/atom concatenation order — the
    // atom table from ligand.js has no x/y/z, so read them off the input.
    const ref = new Float64Array((nProt + lig.nLigAtoms) * 3);
    ref.set(ff.ref.subarray(0, nProt * 3));
    let a = 0;
    for (const mol of ligands) {
      for (const at of mol.atoms) {
        ref[3 * (nProt + a)] = at.x;
        ref[3 * (nProt + a) + 1] = at.y;
        ref[3 * (nProt + a) + 2] = at.z;
        a++;
      }
    }
    ff.ref = ref;
    ff.n = nProt + lig.nLigAtoms;
    ff.nLigAtoms = lig.nLigAtoms;
    ff.isLigand = new Uint8Array(ff.n);
    ff.isLigand.fill(1, nProt);
  } else {
    ff.n = nProt;
  }

  // Masses: Cα beads default to 110 Da (united ALA), ligand atoms to their
  // united-atom mass (H implicitly folded in).
  ff.masses = new Float64Array(ff.n).fill(par.mass ?? 110);
  if (ligands.length) ff.masses.set(ff.ligandMasses, nProt);

  // Fold-metric selection for rmsd(): protein Cα beads [0, nProt) — ligand
  // coords are rigid internal DOF. Built once so the per-frame call allocates
  // nothing. Formula: physics/observables.js.
  ff._maskProt = new Uint8Array(ff.n);
  ff._maskProt.fill(1, 0, nProt);
}