/**
 * thermo_ligand.js — Stage-2 ligand picker for the ΔH/ΔS thermo path.
 *
 * Pain (Stage-6 finding, preserved): the Loop-2 thermo record path used ALL
 * HETATM hetero groups (BNZ+EPE) for the pocket COM + binding-energy
 * attribution. The surface EPE buffer pulls the all-mol COM off the true
 * benzene cavity (record pocket ALA74/MET102…GLY110 surface loop vs the BNZ
 * cavity ILE78/LEU84/…/LEU118) and inflates |ΔH| by −3.18 kcal/mol
 * (calibration row A −6.82 vs row B BNZ-only −3.64; see
 * scripts/calibration_4w52.mjs and docs/BINDING_LOOP2_DONE.md §13).
 * History is NOT rewritten — the record stands; this module only adds an
 * honest default (auto → BNZ cavity) for new runs.
 *
 * Rule (stated, additive, zero deps):
 *   - "auto" (default) = BNZ-first fallback: when any molecule has
 *     resName "BNZ", the subset is all BNZ molecules (the 4W52 benzene
 *     cavity). Otherwise the first molecule. Otherwise none.
 *     This is a documented fallback, NOT a burial/SASA computation — the
 *     task allows "BNZ-first fallback documented; state rule".
 *   - "all" = every hetero group (the Loop-2 record path, bit-identical).
 *   - explicit index ("0", "1", …) or resName ("BNZ", "EPE") = that subset.
 * Pocket rule unchanged: residues within 8 Å of the SELECTED-ligand COM.
 * Binding-energy attribution likewise comes from the SELECTED subset only
 * (caller slices recorded frames / rebuilds the apo leg with the subset).
 *
 * Units: Å for coordinates/distances; energies kcal/mol (caller-side).
 */

export const THERMO_LIG_AUTO = "auto";
export const THERMO_LIG_ALL = "all";
/** Pocket rule (Å) shared with the thermo handler + calibration. */
export const THERMO_POCKET_RCUT = 8.0;

/**
 * Human label for one ligand molecule in the picker.
 * @param {object} mol parseLigands()/parseMol2() molecule { resName, chain, atoms }
 * @param {number} idx molecule index in state.ligands order
 * @returns {string}
 */
export function ligandLabel(mol, idx) {
  const name = String(mol?.resName ?? "LIG").toUpperCase() || "LIG";
  const n = mol?.atoms?.length ?? 0;
  return `${name} #${idx + 1} (${n} atoms)`;
}

/**
 * Resolve the picker value to a ligand subset.
 * @param {Array} ligands state.ligands (parseLigands/parseMol2 molecules)
 * @param {string} [sel] picker value ("auto" default, "all", index, or resName)
 * @returns {{mode:string, subset:Array, molIdx:number[], label:string}}
 */
export function resolveThermoLigand(ligands, sel) {
  const mols = Array.isArray(ligands) ? ligands : [];
  const v = String(sel ?? THERMO_LIG_AUTO).trim() || THERMO_LIG_AUTO;
  if (!mols.length) return { mode: "none", subset: [], molIdx: [], label: "no ligand" };
  const vl = v.toLowerCase();
  if (vl === THERMO_LIG_ALL || vl === "all-ligands") {
    return { mode: "all", subset: mols.slice(), molIdx: mols.map((_, i) => i), label: "all ligands (record path)" };
  }
  // Explicit molecule index ("0", "1", …).
  if (/^\d+$/.test(v)) {
    const i = Number(v);
    if (i >= 0 && i < mols.length) {
      return { mode: "explicit", subset: [mols[i]], molIdx: [i], label: ligandLabel(mols[i], i) };
    }
  }
  // Explicit resName ("BNZ", "EPE", …; case-insensitive, all matches).
  const byName = mols.map((m, i) => ({ m, i }))
    .filter(({ m }) => String(m?.resName ?? "").toUpperCase() === v.toUpperCase());
  if (byName.length) {
    return {
      mode: "explicit",
      subset: byName.map(({ m }) => m),
      molIdx: byName.map(({ i }) => i),
      label: byName.map(({ m, i }) => ligandLabel(m, i)).join(" + "),
    };
  }
  // Default "auto": BNZ-first fallback (documented; 4W52 cavity = benzene).
  const bnz = mols.map((m, i) => ({ m, i }))
    .filter(({ m }) => String(m?.resName ?? "").toUpperCase() === "BNZ");
  if (bnz.length) {
    return {
      mode: "auto-bnz",
      subset: bnz.map(({ m }) => m),
      molIdx: bnz.map(({ i }) => i),
      label: `auto → ${bnz.map(({ m, i }) => ligandLabel(m, i)).join(" + ")} (BNZ cavity default)`,
    };
  }
  return {
    mode: "auto-first",
    subset: [mols[0]],
    molIdx: [0],
    label: `auto → ${ligandLabel(mols[0], 0)} (first group fallback)`,
  };
}

/**
 * Starting ligand-atom offset per molecule in concatenation order
 * (matches ForceField ref-suffix layout: molecule/atom concatenation).
 * @param {Array} ligands
 * @returns {number[]} offsets[i] = first ligand-table index of molecule i
 */
export function ligandAtomOffsets(ligands) {
  const offs = [];
  let k = 0;
  for (const mol of ligands ?? []) {
    offs.push(k);
    k += mol?.atoms?.length ?? 0;
  }
  return offs;
}

/**
 * Global ligand-table atom indices (0..nLigAtoms-1) for a resolved subset.
 * @param {Array} ligands full list (for offsets)
 * @param {number[]} molIdx resolved molecule indices
 * @returns {number[]}
 */
export function selectedLigandAtomIndices(ligands, molIdx) {
  const offs = ligandAtomOffsets(ligands);
  const out = [];
  for (const i of molIdx ?? []) {
    const mol = ligands?.[i];
    if (!mol) continue;
    for (let a = 0; a < (mol.atoms?.length ?? 0); a++) out.push(offs[i] + a);
  }
  return out;
}

/**
 * COM of the selected ligand atoms in a reference/position buffer.
 * @param {ArrayLike} ref flat 3n coords (protein nProt + ligand suffix)
 * @param {number} nProt protein particle count
 * @param {number[]} selAtomIdx ligand-table indices for the selection
 * @returns {number[]} [x, y, z]
 */
export function selectedLigandCom(ref, nProt, selAtomIdx) {
  const com = [0, 0, 0];
  if (!selAtomIdx?.length) return com;
  for (const a of selAtomIdx) {
    com[0] += ref[3 * (nProt + a)] / selAtomIdx.length;
    com[1] += ref[3 * (nProt + a) + 1] / selAtomIdx.length;
    com[2] += ref[3 * (nProt + a) + 2] / selAtomIdx.length;
  }
  return com;
}

/**
 * Pocket residue indices within rCut of a COM (8 Å rule unchanged).
 * @param {ArrayLike} ref flat 3n reference coords
 * @param {number} nProt protein count
 * @param {number[]} com [x, y, z]
 * @param {number} [rCut] Å (default 8.0)
 * @returns {number[]}
 */
export function pocketFromCom(ref, nProt, com, rCut = THERMO_POCKET_RCUT) {
  const idx = [];
  for (let i = 0; i < nProt; i++) {
    const d = Math.hypot(ref[3 * i] - com[0], ref[3 * i + 1] - com[1], ref[3 * i + 2] - com[2]);
    if (d < rCut) idx.push(i);
  }
  return idx;
}

/**
 * Slice one full-system frame (protein + ALL ligand atoms) down to
 * protein + SELECTED ligand atoms (for selected-only energy recomputation).
 * @param {ArrayLike} fullPos flat 3*(nProt+nLigTotal) coords
 * @param {number} nProt protein count
 * @param {number[]} selAtomIdx selected ligand-table indices
 * @returns {Float64Array} flat 3*(nProt+sel.length) coords
 */
export function sliceSelectedPositions(fullPos, nProt, selAtomIdx) {
  const out = new Float64Array((nProt + selAtomIdx.length) * 3);
  out.set(fullPos.subarray(0, nProt * 3), 0);
  for (let k = 0; k < selAtomIdx.length; k++) {
    const a = selAtomIdx[k];
    out[3 * (nProt + k)] = fullPos[3 * (nProt + a)];
    out[3 * (nProt + k) + 1] = fullPos[3 * (nProt + a) + 1];
    out[3 * (nProt + k) + 2] = fullPos[3 * (nProt + a) + 2];
  }
  return out;
}

/**
 * (Re)build the <select id="thermoLig"> options from the live ligand list.
 * Additive: preserves the current value when still valid, else resets to auto.
 * No-op headless (no DOM) or when the element is absent.
 * @param {HTMLSelectElement|null} selectEl
 * @param {Array} ligands
 * @param {string} [current] current value to preserve
 * @returns {string} effective value after refresh ("auto" default)
 */
export function refreshThermoLigOptions(selectEl, ligands, current) {
  if (!selectEl || typeof document === "undefined") {
    return String(current ?? THERMO_LIG_AUTO) || THERMO_LIG_AUTO;
  }
  const mols = Array.isArray(ligands) ? ligands : [];
  const prev = String(current ?? selectEl.value ?? THERMO_LIG_AUTO) || THERMO_LIG_AUTO;
  selectEl.textContent = "";
  const add = (value, text) => {
    const o = document.createElement("option");
    o.value = value;
    o.textContent = text;
    selectEl.appendChild(o);
  };
  add(THERMO_LIG_AUTO, mols.some((m) => String(m?.resName ?? "").toUpperCase() === "BNZ")
    ? "auto/cavity (BNZ-only default)"
    : "auto (first group)");
  mols.forEach((mol, i) => add(String(i), ligandLabel(mol, i)));
  add(THERMO_LIG_ALL, `all (${mols.length} group(s), record path)`);
  const valid = [THERMO_LIG_AUTO, THERMO_LIG_ALL, ...mols.map((_, i) => String(i)),
    ...mols.map((m) => String(m?.resName ?? "").toUpperCase())];
  selectEl.value = valid.includes(prev) ? prev
    : valid.includes(String(selectEl.value)) ? selectEl.value : THERMO_LIG_AUTO;
  return selectEl.value;
}
