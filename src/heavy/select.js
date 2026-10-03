/**
 * heavy/select.js — turn a parsed heavy system into the system the force field
 * is built from: chain/residue/hetero filtering, then external-ligand append.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10).
 *
 *   selectHeavy         chain + residue-range + hetero-group selection, and the
 *                       ordering that fixes the force field index blocks:
 *                       protein, then metals/cofactors, then PDB ligands
 *   appendHeavyLigands  append resolved MOL2/library ligand molecules, with the
 *                       opt-in GAFF2 typing + charge pass (opts.gaff)
 */
import { typeMolecule, assignCharges as gaffAssignCharges } from "../chem/gaff2_mapper.js";

/**
 * Filter parseHeavy() output by chain / residue range and hetero-group selection.
 */
export function selectHeavy(parsedHeavy, { chains = null, resFrom = null, resTo = null, heteroSelection = null, includePdbLigands = true, hasExternalLigand = false } = {}) {
  const inRange = (a) =>
    (!chains || chains.includes(a.chain)) &&
    (resFrom === null || a.resSeq >= resFrom) &&
    (resTo === null || a.resSeq <= resTo);
  const proteinAtoms = parsedHeavy.atoms.filter((a) => a.isProtein && inRange(a));
  const heteroAtoms = parsedHeavy.atoms.filter((a) =>
    !a.isProtein &&
    inRange(a) &&
    (heteroSelection == null || heteroSelection[a.heteroKey] === true)
  );

  const metalsAndCofactors = [];
  const pdbLigandAtoms = [];

  for (const a of heteroAtoms) {
    if (a.isMetal || hasExternalLigand || !includePdbLigands) {
      metalsAndCofactors.push({ ...a, isLigand: false });
    } else {
      pdbLigandAtoms.push({ ...a, isLigand: true });
    }
  }

  const atoms = proteinAtoms.concat(metalsAndCofactors).concat(pdbLigandAtoms);
  const beads = atoms.map((a) => ({ ...a, x: a.x, y: a.y, z: a.z }));
  return { atoms, beads, segments: [], heavy: true, pdbLigandAtoms };
}

/**
 * Append resolved external ligand molecules to selection.
 *
 * Phase 2 opt-in: opts.gaff === true runs the GAFF2-lite typer
 * (src/chem/gaff2_mapper.js typeMolecule + Gasteiger/AM1-BCC-lite charges)
 * per molecule and stores the result on the appended atoms
 * (atom.charge, atom.gaffType). Any per-molecule failure falls back to the
 * parsed charges so a bad ligand can never break system construction.
 * Default (opts.gaff falsy) preserves legacy behavior exactly.
 */
export function appendHeavyLigands(sel, molecules, opts = {}) {
  if (!molecules || molecules.length === 0) return sel;
  const atoms = sel.atoms.slice();
  const beads = sel.beads.slice();
  let resSeq = 1;
  for (const mol of molecules) {
    const resName = (mol.resName || "LIG").slice(0, 3).toUpperCase();
    const chain = mol.chain || "L";
    if (opts.gaff) {
      try {
        typeMolecule(mol.atoms, mol.bonds ?? []);
        gaffAssignCharges(mol.atoms, mol.bonds ?? [], { writeBack: true });
      } catch (e) {
        console.warn(`[appendHeavyLigands] GAFF2 fallback failed for ${resName} (${e?.message}) — parsed charges kept`);
      }
    }
    for (const at of mol.atoms) {
      const atom = {
        x: at.x, y: at.y, z: at.z,
        element: at.element,
        resName, chain, resSeq,
        atomName: at.atomName || at.element,
        serial: at.serial ?? 0,
        // GAFF2 fields exist only on the opt-in path so legacy charge
        // assignment (element defaults via assignCharges) is untouched.
        ...(opts.gaff ? { charge: at.charge ?? 0, gaffType: at.gaffType ?? null } : {}),
        isProtein: false, isMetal: false, isLigand: true, isHetero: false, heteroKey: null,
      };
      atoms.push(atom);
      beads.push({ ...atom });
    }
    resSeq++;
  }
  return { atoms, beads, segments: sel.segments, heavy: true };
}
