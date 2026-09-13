/**
 * alanine_scanning.js — computational alanine scanning via a thermodynamic cycle.
 *
 * Click (or script) a residue → its sidechain is truncated to CB (alanine),
 * the topology is updated, and both the mutant apo and holo states are
 * short-relaxed. The binding-ΔΔG follows from the standard cycle
 * (ligand self-energy cancels):
 *
 *   ΔG_bind(X) = G_holo(X) − G_apo(X) − G_lig
 *   ΔΔG_bind   = ΔG_bind(mut) − ΔG_bind(WT)
 *              = [G_holo(mut) − G_holo(WT)] − [G_apo(mut) − G_apo(WT)]
 *              = ΔG_mut,holo − ΔG_mut,apo                         (kcal/mol)
 *
 * Each G is the steepest-descent relaxed potential energy from relaxEnergy().
 * Positive ΔΔG ⇒ the WT sidechain stabilizes binding (hotspot); negative ⇒
 * the truncation improves binding.
 *
 * Mode handling (units Å / ps / kcal/mol / Da, vanilla ES modules):
 *   CG Cα  — no sidechain atoms exist, so the mutation is an ENM perturbation:
 *            resName → ALA (resets the residue-class LJ/Coulomb/H-bond tables)
 *            and every spring touching the bead is rescaled,
 *              k_ij(mut) = alaScale · k_ij(WT),   alaScale = 0.7 default,
 *            encoding the smaller alanine contact surface. WT spring weights
 *            (uniform / Tirion / sequence / ML prior) are copied verbatim so
 *            the WT↔mut comparison is apples-to-apples.
 *   heavy  — explicit topology surgery: protein atoms of the residue whose
 *            atomName ∉ {N, CA, C, O, OXT, CB} are deleted, resName → ALA,
 *            and the HeavyForceField (bonds/angles/GB/SASA tables) is rebuilt
 *            from the reduced atom set.
 *
 * Glycine (no CB) and alanine (identity) mutants are defined as ΔΔG = 0 with
 * a note — finite by construction so pocket tables never contain NaN.
 *
 * DOM-free: runs under plain Node for unit testing. UI wiring lives in
 * src/analysis-panel.js (opt-in buttons, graceful degradation).
 */

import { ForceField } from "../forcefield.js?v=10";
import { HeavyForceField } from "../heavy.js?v=10";

/** Heavy-atom names kept by an alanine truncation (backbone + CB + terminal O). */
export const ALA_KEPT_ATOMS = new Set(["N", "CA", "C", "O", "OXT", "CB"]);

/** Default ENM spring rescale for a CG alanine mutation (smaller contact surface). */
export const ALA_SPRING_SCALE = 0.7;

/**
 * Default rest-length contraction (Å) applied to a CG alanine mutant's springs.
 * Alanine is smaller than any other sidechain, so neighbouring Cα beads pack
 * closer; without an r0 shift a pure stiffness change relaxes back to zero
 * energy at the native minimum. Shift is shared-contact physics (≈ one
 * methyl radius split between partners), floored at R0_MIN.
 */
export const ALA_R0_SHIFT = 0.4;
/** Minimum ENM rest length after the alanine shift (Cα–Cα ≈ 3.8 Å native). */
export const ALA_R0_MIN = 3.2;

/**
 * Display noise floor for CG alanine-scan ΔΔG (kcal/mol).
 *
 * Stage-4 honest-display resolution (measure/honest-display only, no
 * force-field retuning): the CG Cα-ENM perturbation (ALA_SPRING_SCALE +
 * ALA_R0_SHIFT) largely relaxes back and cancels in the holo−apo cycle, so
 * live 4W52 BNZ-cavity magnitudes are all |ΔΔG| < 0.02 — below any plausible
 * relaxation/convergence precision. Rows with |ΔΔG| below this floor are
 * flagged "~noise" by formatMutationTable: ranking-only, no hotspot
 * resolution at CG. Heavy-mode explicit sidechain surgery is the honest
 * resolution path (open bench item). Deterministic pure constant.
 */
export const ALA_DDG_NOISE_FLOOR = 0.05;

/**
 * Below-floor test for a scan ΔΔG (deterministic, no RNG).
 * Non-finite inputs count as noise (never claim signal on NaN).
 * @param {number} ddG  binding ΔΔG (kcal/mol)
 * @param {number} [floor=ALA_DDG_NOISE_FLOOR]  noise floor (kcal/mol)
 * @returns {boolean} true when |ΔΔG| < floor (ranking-only)
 */
export function isNoiseDdG(ddG, floor = ALA_DDG_NOISE_FLOOR) {
  if (!Number.isFinite(ddG) || !Number.isFinite(floor)) return true;
  return Math.abs(ddG) < floor;
}

/**
 * Tag scan rows with honest-display noise fields (additive, in place).
 * Sets r.noise (bool) + r.noiseFloor (floor used). Deterministic.
 * @param {Array<object>} rows  scanPocket().rows
 * @param {number} [floor=ALA_DDG_NOISE_FLOOR]
 * @returns {Array<object>} same array (chainable)
 */
export function annotateScanNoise(rows, floor = ALA_DDG_NOISE_FLOOR) {
  for (const r of rows) {
    r.noise = isNoiseDdG(r.ddG, floor);
    r.noiseFloor = floor;
  }
  return rows;
}

/**
 * @typedef {object} ScanSystem
 * @property {"cg"|"heavy"} mode
 * @property {object} sel      selectSystem() output (cg) or selectHeavy() output (heavy)
 * @property {object} ff       ForceField (cg) or HeavyForceField (heavy)
 * @property {object} [par]    force-field parameters used at build (optional)
 * @property {Array}  [ligands] ligand molecules (cg holo only; heavy carries atoms in sel)
 */

/**
 * Canonical residue key.
 * @param {string} chain
 * @param {number} resSeq
 * @returns {string} e.g. "A|101"
 */
export function resKey(chain, resSeq) {
  return `${chain}|${resSeq}`;
}

/**
 * Human label for a residue.
 * @param {object} r {resName, chain, resSeq}
 * @returns {string} e.g. "LEU101:A"
 */
export function resLabel(r) {
  return `${r.resName}${r.resSeq}:${r.chain}`;
}

/**
 * Parse a user/click residue identifier against a system.
 *
 * Accepted resId forms: CG bead index (number); "CHAIN|SEQ" / "CHAIN:SEQ" /
 * "SEQ" string; {chain, resSeq} object.
 *
 * @param {ScanSystem} system
 * @param {number|string|object} resId
 * @returns {{beadIndex: number, chain: string, resSeq: number, resName: string,
 *   label: string, key: string, atomIndices: number[]}}
 */
export function resolveResidue(system, resId) {
  const { mode, sel, ff } = system;
  if (!sel || !ff) throw new Error("resolveResidue: system needs {sel, ff} (build the system first).");
  if (mode === "cg") {
    const beads = sel.beads;
    let idx = -1;
    if (typeof resId === "number") {
      idx = resId;
    } else {
      let chain = null, seq = null;
      if (typeof resId === "string") {
        const m = resId.trim().toUpperCase().replace(":", "|").split("|");
        if (m.length === 2) { chain = m[0]; seq = Number(m[1]); }
        else { seq = Number(m[0]); }
      } else if (resId && typeof resId === "object") {
        chain = resId.chain != null ? String(resId.chain).toUpperCase() : null;
        seq = Number(resId.resSeq);
      }
      if (!Number.isFinite(seq)) throw new Error(`resolveResidue: cannot parse residue "${resId}".`);
      idx = beads.findIndex((b) => b.resSeq === seq && (chain === null || (b.chain || "A") === chain));
    }
    if (!Number.isInteger(idx) || idx < 0 || idx >= ff.nProt) {
      throw new Error(`resolveResidue: residue "${resId}" not in the CG selection (nProt=${ff.nProt}).`);
    }
    const b = beads[idx];
    return {
      beadIndex: idx, chain: b.chain || "A", resSeq: b.resSeq, resName: b.resName,
      label: resLabel({ resName: b.resName, chain: b.chain || "A", resSeq: b.resSeq }),
      key: resKey(b.chain || "A", b.resSeq), atomIndices: [idx],
    };
  }
  // heavy mode
  const atoms = sel.atoms;
  let key = null;
  if (typeof resId === "number") {
    const a = atoms[resId];
    if (!a || !a.isProtein) throw new Error(`resolveResidue: atom ${resId} is not a protein atom.`);
    key = resKey(a.chain, a.resSeq);
  } else if (typeof resId === "string") {
    const m = resId.trim().toUpperCase().replace(":", "|").split("|");
    key = m.length === 2 ? resKey(m[0], Number(m[1])) : null;
    if (key === null || key.endsWith("|NaN")) {
      // fall back: treat as "SEQ" on the single/first chain
      const seq = Number(m[0]);
      const hit = atoms.find((a) => a.isProtein && a.resSeq === seq);
      if (!hit) throw new Error(`resolveResidue: residue "${resId}" not found.`);
      key = resKey(hit.chain, hit.resSeq);
    }
  } else if (resId && typeof resId === "object") {
    if (Number.isInteger(resId.beadIndex) || Number.isInteger(resId.caIndex)) {
      return resolveResidue(system, resId.beadIndex ?? resId.caIndex);
    }
    key = resKey(String(resId.chain).toUpperCase(), Number(resId.resSeq));
  } else {
    throw new Error(`resolveResidue: cannot parse residue "${resId}".`);
  }
  const idxs = [];
  for (let i = 0; i < atoms.length; i++) {
    const a = atoms[i];
    if (a.isProtein && resKey(a.chain, a.resSeq) === key) idxs.push(i);
  }
  if (!idxs.length) throw new Error(`resolveResidue: residue "${key}" not in the heavy selection.`);
  const rep = atoms[idxs[0]];
  return {
    beadIndex: idxs.find((i) => atoms[i].atomName === "CA") ?? idxs[0],
    chain: rep.chain, resSeq: rep.resSeq, resName: rep.resName,
    label: resLabel(rep), key, atomIndices: idxs,
  };
}

/**
 * Map a viewer pick (global particle index) to a residue label.
 * @param {object} sel  selection with .beads (cg) or .atoms (heavy)
 * @param {number} particleIndex
 * @param {"cg"|"heavy"} [mode="cg"]
 * @returns {string} residue label
 */
export function resIdFromBeadIndex(sel, particleIndex, mode = "cg") {
  if (mode === "cg") {
    const b = sel.beads[particleIndex];
    if (!b) throw new Error(`resIdFromBeadIndex: bead ${particleIndex} out of range.`);
    return resLabel({ resName: b.resName, chain: b.chain || "A", resSeq: b.resSeq });
  }
  const a = sel.atoms[particleIndex];
  if (!a) throw new Error(`resIdFromBeadIndex: atom ${particleIndex} out of range.`);
  return a.isProtein ? resKey(a.chain, a.resSeq) : `ligand:${a.resName}`;
}

/**
 * Rebuild parameters from a live CG force field so mutant/apo rebuilds match it.
 * @param {object} ff ForceField
 * @returns {object} par compatible with `new ForceField(sel, par, ligands)`
 */
export function parFromCgFF(ff) {
  return {
    rc: ff.rc, gamma: ff.gamma, kBond: ff.kBond, kAngle: ff.kAngle,
    epsRep: ff.epsRep, sigmaRep: ff.sigmaRep,
    mass: ff.masses && ff.masses.length ? ff.masses[0] : 110,
    binding: { on: ff.bindOn !== false, holo: ff.holoOn !== false },
    enmModel: ff.enmModel ?? "uniform",
    seqWeight: false, // spring weights are copied verbatim below; never re-derive
  };
}

/**
 * Build the alanine mutant of one residue. Particle count is unchanged in CG
 * mode; heavy mode deletes sidechain atoms beyond CB.
 *
 * @param {ScanSystem} system
 * @param {number|string|object} resId
 * @param {object} [opts]
 * @param {number} [opts.alaScale=ALA_SPRING_SCALE] CG spring rescale factor
 * @param {number} [opts.r0Shift=ALA_R0_SHIFT] CG rest-length contraction (Å)
 * @returns {{system: ScanSystem, res: object, scaledSprings: number,
 *   removedAtoms: number, note: string}} mutant system + bookkeeping
 */
export function buildAlanineMutant(system, resId, opts = {}) {
  const alaScale = opts.alaScale ?? ALA_SPRING_SCALE;
  const r0Shift = opts.r0Shift ?? ALA_R0_SHIFT;
  const res = resolveResidue(system, resId);
  const isIdentity = res.resName === "ALA" || res.resName === "GLY";

  if (system.mode === "cg") {
    const beads = system.sel.beads.map((b, i) =>
      i === res.beadIndex ? { ...b, resName: "ALA" } : b);
    const selMut = { ...system.sel, beads };
    const par = system.par ?? parFromCgFF(system.ff);
    const ffMut = new ForceField(selMut, par, system.ligands ?? []);
    // Apples-to-apples springs: copy WT weights (uniform/Tirion/seq/ML), then
    // perturb only the mutated bead's contacts — stiffness rescale PLUS rest-
    // length contraction (smaller sidechain packs closer). Spring order is
    // deterministic from identical reference coordinates, hence index-aligned.
    let scaledSprings = 0;
    if (ffMut.springK.length === system.ff.springK.length) {
      ffMut.springK.set(system.ff.springK);
      ffMut.springs.set(system.ff.springs);
      if (!isIdentity) {
        for (let q = 0; q < ffMut.springK.length; q++) {
          const i = ffMut.springs[3 * q], j = ffMut.springs[3 * q + 1];
          if (i === res.beadIndex || j === res.beadIndex) {
            ffMut.springK[q] *= alaScale;
            ffMut.springs[3 * q + 2] = Math.max(ALA_R0_MIN, ffMut.springs[3 * q + 2] - r0Shift);
            scaledSprings++;
          }
        }
      }
    }
    return {
      system: { mode: "cg", sel: selMut, ff: ffMut, par, ligands: system.ligands ?? [] },
      res, scaledSprings, removedAtoms: 0,
      note: isIdentity ? `${res.resName} identity (ΔΔG ≡ 0)` : `CG: ${scaledSprings} springs ×${alaScale}, Δr0 −${r0Shift}Å`,
    };
  }

  // heavy mode — explicit sidechain deletion
  const keepSet = new Set(res.atomIndices.filter((i) =>
    ALA_KEPT_ATOMS.has(system.sel.atoms[i].atomName)));
  const dropSet = new Set(res.atomIndices.filter((i) => !keepSet.has(i)));
  const oldToNew = new Int32Array(system.sel.atoms.length).fill(-1);
  const mutAtoms = [];
  for (let i = 0; i < system.sel.atoms.length; i++) {
    if (dropSet.has(i)) continue;
    const a = system.sel.atoms[i];
    oldToNew[i] = mutAtoms.length;
    mutAtoms.push(i === res.beadIndex || keepSet.has(i) && resKey(a.chain, a.resSeq) === res.key
      ? (a.isProtein && resKey(a.chain, a.resSeq) === res.key ? { ...a, resName: "ALA" } : a)
      : a);
  }
  const par = system.par ?? { gamma: system.ff.gamma ?? 1.0 };
  const selMut = { ...system.sel, atoms: mutAtoms };
  const ffMut = new HeavyForceField({ atoms: mutAtoms }, par, []);
  return {
    system: { mode: "heavy", sel: selMut, ff: ffMut, par, ligands: [] },
    res, scaledSprings: 0, removedAtoms: dropSet.size,
    note: isIdentity
      ? `${res.resName} identity (ΔΔG ≡ 0)`
      : `heavy: removed ${dropSet.size} sidechain atoms → CB`,
    _oldToNew: oldToNew,
  };
}

/**
 * Build the apo (ligand-free) variant of a system. Coordinates are unchanged;
 * only the ligand particles / molecules are dropped.
 * @param {ScanSystem} system
 * @returns {ScanSystem}
 */
export function apoSystemOf(system) {
  if (system.mode === "cg") {
    const par = system.par ?? parFromCgFF(system.ff);
    const ffApo = new ForceField(system.sel, par, []);
    if (ffApo.springK.length === system.ff.springK.length) ffApo.springK.set(system.ff.springK);
    return { mode: "cg", sel: system.sel, ff: ffApo, par, ligands: [] };
  }
  const atoms = system.sel.atoms.filter((a) => !a.isLigand);
  const par = system.par ?? { gamma: system.ff.gamma ?? 1.0 };
  const ffApo = new HeavyForceField({ atoms }, par, []);
  return { mode: "heavy", sel: { ...system.sel, atoms }, ff: ffApo, par, ligands: [] };
}

/**
 * Short deterministic relaxation (steepest descent with per-step clamp):
 *   x ← x + clamp(step · F, ±maxDisp),   F = −∇U in kcal/mol/Å.
 * Tracks the best energy seen so the returned value is monotonic.
 *
 * @param {object} ff  force field with .compute(pos) / .forces
 * @param {Float64Array} pos0  starting positions, length 3n
 * @param {object} [opts]
 * @param {number} [opts.steps=80]
 * @param {number} [opts.step=0.02]   Å·mol·Å/kcal gain
 * @param {number} [opts.maxDisp=0.05] Å per-coordinate clamp
 * @param {number} [opts.tol=1e-3]    kcal/mol/Å max-force convergence
 * @returns {{energy: number, pos: Float64Array, steps: number,
 *   maxForce: number, converged: boolean}}
 */
export function relaxEnergy(ff, pos0, opts = {}) {
  const steps = opts.steps ?? 80;
  const step = opts.step ?? 0.02;
  const maxDisp = opts.maxDisp ?? 0.05;
  const tol = opts.tol ?? 1e-3;
  const pos = Float64Array.from(pos0);
  let E = ff.compute(pos);
  let best = E, bestPos = pos.slice(), maxForce = Infinity, converged = false;
  if (!Number.isFinite(E)) return { energy: E, pos, steps: 0, maxForce: Infinity, converged: false };
  for (let s = 0; s < steps; s++) {
    E = ff.compute(pos);
    const F = ff.forces;
    maxForce = 0;
    for (let i = 0; i < F.length; i++) {
      const a = Math.abs(F[i]);
      if (a > maxForce) maxForce = a;
    }
    if (!Number.isFinite(E) || !Number.isFinite(maxForce)) break;
    if (E < best) { best = E; bestPos = pos.slice(); }
    if (maxForce < tol) { converged = true; break; }
    for (let i = 0; i < pos.length; i++) {
      let d = step * F[i];
      if (d > maxDisp) d = maxDisp;
      else if (d < -maxDisp) d = -maxDisp;
      pos[i] += d;
    }
  }
  return { energy: best, pos: bestPos, steps, maxForce, converged };
}

/**
 * Restrict a holo-relaxed position vector to an apo/mutant particle subset.
 * @param {Float64Array} posFull  length 3nFull
 * @param {number[]} keep  sorted global indices to keep
 * @returns {Float64Array}
 */
export function subsetPositions(posFull, keep) {
  const out = new Float64Array(keep.length * 3);
  for (let k = 0; k < keep.length; k++) {
    out[3 * k] = posFull[3 * keep[k]];
    out[3 * k + 1] = posFull[3 * keep[k] + 1];
    out[3 * k + 2] = posFull[3 * keep[k] + 2];
  }
  return out;
}

/**
 * Alanine-scan a single residue through the thermodynamic cycle.
 *
 * @param {ScanSystem} system  WT holo system
 * @param {number|string|object} resId
 * @param {object} [opts]
 * @param {number} [opts.relaxSteps=80]   steepest-descent steps per leg
 * @param {number} [opts.alaScale]        CG spring rescale (default 0.7)
 * @param {Float64Array} [opts.startPos]  WT holo start (default ff.ref)
 * @param {object} [opts.cache]           shared {wtHolo, wtApo} energies/positions
 * @returns {{resId: string, label: string, chain: string, resSeq: number,
 *   wtRes: string, dGmutHolo: number, dGmutApo: number, ddG: number,
 *   scaledSprings: number, removedAtoms: number, note: string}}
 */
export function scanResidue(system, resId, opts = {}) {
  const relaxSteps = opts.relaxSteps ?? 80;
  const cache = opts.cache ?? {};
  const res = resolveResidue(system, resId);
  const row = {
    resId: res.key, label: `${res.resName}${res.resSeq}:${res.chain}`,
    chain: res.chain, resSeq: res.resSeq, wtRes: res.resName,
    dGmutHolo: 0, dGmutApo: 0, ddG: 0,
    scaledSprings: 0, removedAtoms: 0, note: "",
  };
  if (res.resName === "ALA" || res.resName === "GLY") {
    row.note = `${res.resName} identity (ΔΔG ≡ 0)`;
    return row;
  }
  // WT legs (cached across a pocket scan)
  if (!cache.wtHolo) {
    const start = opts.startPos ?? system.ff.ref;
    const r = relaxEnergy(system.ff, start, { steps: relaxSteps });
    cache.wtHolo = r;
  }
  if (!cache.wtApo) {
    const apo = apoSystemOf(system);
    cache._apoSys = apo;
    let startApo;
    if (system.mode === "cg") {
      startApo = cache.wtHolo.pos; // same particle count
    } else {
      const keep = [];
      for (let i = 0; i < system.sel.atoms.length; i++) {
        if (!system.sel.atoms[i].isLigand) keep.push(i);
      }
      cache._apoKeep = keep;
      startApo = subsetPositions(cache.wtHolo.pos, keep);
    }
    cache.wtApo = relaxEnergy(apo.ff, startApo, { steps: relaxSteps });
  }
  const apoSys = cache._apoSys ?? apoSystemOf(system);

  // Mutant legs
  const mut = buildAlanineMutant(system, resId, opts);
  row.scaledSprings = mut.scaledSprings;
  row.removedAtoms = mut.removedAtoms;
  row.note = mut.note;
  let mutHoloStart = cache.wtHolo.pos;
  if (system.mode === "heavy" && mut._oldToNew) {
    const keep = [];
    for (let i = 0; i < mut._oldToNew.length; i++) if (mut._oldToNew[i] >= 0) keep.push(i);
    mutHoloStart = subsetPositions(cache.wtHolo.pos, keep);
    // mutant apo: drop ligand atoms from the mutant holo subset
    const mutApoKeep = [];
    for (let k = 0; k < mut.system.sel.atoms.length; k++) {
      if (!mut.system.sel.atoms[k].isLigand) mutApoKeep.push(k);
    }
    const mutApoSys = apoSystemOf(mut.system);
    const rMutHolo = relaxEnergy(mut.system.ff, mutHoloStart, { steps: relaxSteps });
    const rMutApo = relaxEnergy(mutApoSys.ff, subsetPositions(rMutHolo.pos, mutApoKeep), { steps: relaxSteps });
    row.dGmutHolo = rMutHolo.energy - cache.wtHolo.energy;
    row.dGmutApo = rMutApo.energy - cache.wtApo.energy;
    row.ddG = row.dGmutHolo - row.dGmutApo;
    return row;
  }
  // CG: identical particle counts everywhere. The mutant apo rebuild
  // re-derives springs from coordinates, so re-apply the mutant's perturbed
  // spring set (weights + shifted rest lengths) for a consistent cycle.
  const mutApoSys = apoSystemOf(mut.system);
  if (mutApoSys.ff.springK.length === mut.system.ff.springK.length) {
    mutApoSys.ff.springK.set(mut.system.ff.springK);
    mutApoSys.ff.springs.set(mut.system.ff.springs);
  }
  const rMutHolo = relaxEnergy(mut.system.ff, mutHoloStart, { steps: relaxSteps });
  const rMutApo = relaxEnergy(mutApoSys.ff, rMutHolo.pos.subarray(0, mutApoSys.ff.n * 3), { steps: relaxSteps });
  row.dGmutHolo = rMutHolo.energy - cache.wtHolo.energy;
  row.dGmutApo = rMutApo.energy - cache.wtApo.energy;
  row.ddG = row.dGmutHolo - row.dGmutApo;
  void apoSys;
  return row;
}

/**
 * Scan a pocket (or any residue list) with shared WT legs.
 * @param {ScanSystem} system  WT holo system
 * @param {Array<number|string|object>} resIds
 * @param {object} [opts]  forwarded to scanResidue
 * @returns {{rows: object[], wtHolo: number, wtApo: number, n: number}}
 */
export function scanPocket(system, resIds, opts = {}) {
  const cache = {};
  const rows = resIds.map((id) => scanResidue(system, id, { ...opts, cache }));
  rows.sort((a, b) => b.ddG - a.ddG);
  annotateScanNoise(rows, opts.noiseFloor ?? ALA_DDG_NOISE_FLOOR);
  return {
    rows,
    wtHolo: cache.wtHolo ? cache.wtHolo.energy : NaN,
    wtApo: cache.wtApo ? cache.wtApo.energy : NaN,
    n: rows.length,
  };
}

/**
 * Pocket-residue discovery: protein beads (CG) or CA atoms (heavy) within
 * rCut of the native ligand COM, sorted by distance.
 * @param {ScanSystem} system
 * @param {object} [opts]
 * @param {number} [opts.rCut=6.0]  Å
 * @param {number} [opts.maxN=12]
 * @returns {Array<{resId: (number|string), label: string, dist: number}>}
 */
export function pocketResidues(system, opts = {}) {
  const rCut = opts.rCut ?? 6.0;
  const maxN = opts.maxN ?? 12;
  const ff = system.ff;
  const ligStart = ff.ligandStart ?? ff.nProt;
  const nLig = ff.n - ligStart;
  if (nLig <= 0) return [];
  let lx = 0, ly = 0, lz = 0;
  for (let a = 0; a < nLig; a++) {
    lx += ff.ref[3 * (ligStart + a)];
    ly += ff.ref[3 * (ligStart + a) + 1];
    lz += ff.ref[3 * (ligStart + a) + 2];
  }
  lx /= nLig; ly /= nLig; lz /= nLig;
  const out = [];
  if (system.mode === "cg") {
    for (let i = 0; i < ff.nProt; i++) {
      const dx = ff.ref[3 * i] - lx, dy = ff.ref[3 * i + 1] - ly, dz = ff.ref[3 * i + 2] - lz;
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (d <= rCut) {
        const b = system.sel.beads[i];
        out.push({ resId: i, label: resLabel({ resName: b.resName, chain: b.chain || "A", resSeq: b.resSeq }), dist: d });
      }
    }
  } else {
    const seen = new Set();
    for (let i = 0; i < ff.nProt; i++) {
      const a = system.sel.atoms[i];
      if (!a || !a.isProtein || a.atomName !== "CA") continue;
      const dx = ff.ref[3 * i] - lx, dy = ff.ref[3 * i + 1] - ly, dz = ff.ref[3 * i + 2] - lz;
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const k = resKey(a.chain, a.resSeq);
      if (d <= rCut && !seen.has(k)) {
        seen.add(k);
        out.push({ resId: k, label: resLabel(a), dist: d });
      }
    }
  }
  out.sort((a, b) => a.dist - b.dist);
  return out.slice(0, maxN);
}

/**
 * Fixed-width ΔΔG table for HUD / console / download.
 * Honest display (Stage-4, Option B): rows with |ΔΔG| below the noise floor
 * carry a "~noise" flag and the table ends with a ranking-only disclaimer
 * (hotspot recall needs heavy-mode sidechains; CG magnitudes do not resolve
 * hotspots). Pure display — input energies untouched, deterministic.
 * @param {Array<object>} rows  scanPocket().rows
 * @param {object} [opts]
 * @param {number} [opts.floor=ALA_DDG_NOISE_FLOOR]  noise floor (kcal/mol)
 * @returns {string}
 */
export function formatMutationTable(rows, opts = {}) {
  const floor = opts.floor ?? ALA_DDG_NOISE_FLOOR;
  const head = "mutation      ΔGmut_holo  ΔGmut_apo   ΔΔGbind  note";
  const lines = [head, "-".repeat(head.length)];
  for (const r of rows) {
    const f = (v) => (Number.isFinite(v) ? v.toFixed(2).padStart(10) : "       NaN");
    const noise = (r.noise ?? isNoiseDdG(r.ddG, floor)) ? " ~noise" : "";
    lines.push(`${r.label.padEnd(12)}${f(r.dGmutHolo)}${f(r.dGmutApo)}${f(r.ddG)}  ${r.note}${noise}`);
  }
  lines.push(`[ala-scan: |ΔΔG| < ${floor.toFixed(2)} kcal/mol ≈ noise — ranking only, no CG hotspot resolution]`);
  return lines.join("\n");
}
