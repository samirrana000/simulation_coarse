/**
 * protonation.js — Heuristic client-side pKa / tautomer assigner (Phase 2: Chemistry).
 *
 * Assigns protonation states to titratable protein residues at a given pH using
 * a PROPKA-inspired heuristic (Olsson et al., JCTC 2011, 7, 525 — model pKa
 * values; heuristic desolvation / H-bond / Coulomb shifts computed from heavy-
 * atom geometry only, no Poisson–Boltzmann):
 *
 *   pKa_eff = pKa_model + dDesolv + dHBond + dCoulomb
 *   protonated (neutral acid / charged base)  iff  pH < pKa_eff   (acids)
 *   deprotonated (neutral base)                iff  pKa_eff < pH  (bases)
 *
 * State alphabet (AMBER convention):
 *   HIS -> HIE (H on NE2/epsilon) | HID (H on ND1/delta) | HIP (both, +1)
 *   ASP -> ASP (-1) | ASH (neutral)      GLU -> GLU (-1) | GLH (neutral)
 *   LYS -> LYS (+1) | LYN (neutral)      ARG -> ARG (+1) | ARN (neutral)
 *   CYS -> CYS (free thiol, 0) | CYX (disulfide 0 / metal-bound thiolate -1)
 *
 * Inputs are duck-typed so the assigner runs on any existing parsed structure:
 *   - heavy.js parseHeavy()/selectHeavy() output: { atoms: [...] }
 *   - flat heavy-atom array: [{ atomName|name, resName, chain, resSeq, x,y,z|pos, element|elem }]
 *   - residue array: [{ resName, resSeq, chain, atoms: [{ name, elem, pos|[x,y,z]|x,y,z }] }]
 *
 * No imports, no npm deps, vanilla ES module. Pure computation except for an
 * optional in-place rename helper (applyProtonationStates). Never throws on
 * missing atoms — residues lacking side-chain coordinates fall back to
 * backbone/exposure defaults with an explanatory reason.
 */

/** Model (intrinsic) pKa values in water (Nozaki & Tanford; Antosiewicz et al.). */
export const MODEL_PKA = Object.freeze({
  ASP: 3.9, GLU: 4.3, HIS: 6.5, CYS: 8.3, LYS: 10.5, ARG: 12.5,
});

/** Residue names this assigner can re-type. */
export const TITRATABLE = new Set(["HIS", "ASP", "GLU", "LYS", "ARG", "CYS"]);

/** Formal residue charge after assignment (CYX depends on context; see reason). */
export const STATE_CHARGE = Object.freeze({
  HID: 0, HIE: 0, HIP: 1,
  ASP: -1, ASH: 0, GLU: -1, GLH: 0,
  LYS: 1, LYN: 0, ARG: 1, ARN: 0,
  CYS: 0, CYX: 0, CYM: -1,
});

/** Side-chain heavy atoms (PDB names) used for centroids / environment probes. */
const SIDECHAIN_ATOMS = {
  HIS: ["CG", "ND1", "CD2", "CE1", "NE2"],
  ASP: ["CB", "CG", "OD1", "OD2"],
  GLU: ["CB", "CG", "CD", "OE1", "OE2"],
  LYS: ["CB", "CG", "CD", "CE", "NZ"],
  ARG: ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
  CYS: ["CB", "SG"],
};

const HBOND_DIST = 3.5;   // Å — donor/acceptor contact cutoff
const SALTBRIDGE_DIST = 4.5; // Å — N–O salt-bridge cutoff
const BURIAL_RADIUS = 8.0;   // Å — solvent-exposure proxy sphere
const DISULFIDE_DIST = 2.3;  // Å — SG–SG disulfide cutoff (CSD 2.03–2.05 + margin)
const METAL_S_DIST = 2.8;    // Å — SG–metal coordination cutoff

function clamp01(x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
function dist3(a, b) {
  const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Normalize one atom record to { name, elem, x, y, z, resName, chain, resSeq }.
 * Accepts heavy.js atoms ({atomName, element, x,y,z}), MOL2/HETATM style
 * ({atomName|name, element|elem, x,y,z}), and {name, elem, pos:[x,y,z]}.
 * @param {object} a
 * @param {object} [fallbackRes]  { resName, chain, resSeq } when atom lacks them
 * @returns {object|null}  null when coordinates are non-finite
 */
export function normalizeAtom(a, fallbackRes = {}) {
  if (!a || typeof a !== "object") return null;
  const x = a.x ?? a.pos?.[0];
  const y = a.y ?? a.pos?.[1];
  const z = a.z ?? a.pos?.[2];
  if (!Number.isFinite(x + y + z)) return null;
  const name = String(a.atomName ?? a.name ?? "").trim().toUpperCase();
  const elem = String(a.element ?? a.elem ?? guessElement(name)).trim().toUpperCase();
  return {
    name,
    elem,
    x, y, z,
    resName: String(a.resName ?? fallbackRes.resName ?? "UNK").trim().toUpperCase(),
    chain: String(a.chain ?? fallbackRes.chain ?? "_"),
    resSeq: Number.isFinite(+a.resSeq) ? +a.resSeq : (fallbackRes.resSeq ?? 0),
    ref: a,
  };
}

function guessElement(name) {
  const m = String(name || "").match(/^([A-Za-z]{1,2})/);
  if (!m) return "C";
  const el = m[1].toUpperCase();
  if (el === "CA" && name.trim() === "CA") return "C"; // alpha carbon, not calcium
  if (["CL", "BR", "ZN", "FE", "MG", "CA", "CU", "MN", "NI", "CO", "NA", "SE", "SI"].includes(el)) return el;
  return el[0];
}

/**
 * Group any supported input into residue records.
 * @param {Array|object} input  residues | flat atoms | { atoms } | { beads }
 * @returns {Array<{resName:string,chain:string,resSeq:number,atoms:Array}>}
 */
export function groupResidues(input) {
  let residues = null;
  if (Array.isArray(input) && input.length && Array.isArray(input[0]?.atoms)) {
    residues = input.map((r) => ({
      resName: String(r.resName ?? "UNK").trim().toUpperCase(),
      chain: String(r.chain ?? "_"),
      resSeq: Number.isFinite(+r.resSeq) ? +r.resSeq : 0,
      atoms: r.atoms.map((a) => normalizeAtom(a, r)).filter(Boolean),
    }));
  } else {
    const flat = Array.isArray(input) ? input
      : Array.isArray(input?.atoms) ? input.atoms
      : Array.isArray(input?.beads) ? input.beads.map((b) => ({
          atomName: "CA", element: "C", x: b.x, y: b.y, z: b.z,
          resName: b.resName, chain: b.chain, resSeq: b.resSeq,
        }))
      : [];
    const map = new Map();
    for (const a of flat) {
      const n = normalizeAtom(a);
      if (!n) continue;
      const key = `${n.chain}|${n.resSeq}`;
      if (!map.has(key)) {
        map.set(key, { resName: n.resName, chain: n.chain, resSeq: n.resSeq, atoms: [] });
      }
      map.get(key).atoms.push(n);
    }
    residues = [...map.values()];
  }
  return residues.filter((r) => r.atoms.length > 0);
}

function atomByName(res, names) {
  const set = new Set(names);
  return res.atoms.find((a) => set.has(a.name)) ?? null;
}

function centroid(atoms) {
  const c = { x: 0, y: 0, z: 0 };
  for (const a of atoms) { c.x += a.x; c.y += a.y; c.z += a.z; }
  c.x /= atoms.length; c.y /= atoms.length; c.z /= atoms.length;
  return c;
}

/**
 * Environmenthoood around a point: neighbor counts + specific contacts.
 * @param {{x:number,y:number,z:number}} p
 * @param {Array} allAtoms  normalized atoms (with resKey)
 * @param {string} ownKey   residue key to exclude
 */
function probeEnvironment(p, allAtoms, ownKey) {
  let nTotal = 0, nHydrophobic = 0, nPolar = 0;
  const contacts = []; // { atom, d }
  for (const a of allAtoms) {
    if (a.resKey === ownKey) continue;
    const dx = a.x - p.x, dy = a.y - p.y, dz = a.z - p.z;
    const d2 = dx * dx + dy * dy + dz * dz;
    if (d2 > BURIAL_RADIUS * BURIAL_RADIUS) continue;
    nTotal++;
    if (a.elem === "C" || a.elem === "S") nHydrophobic++; else nPolar++;
    const d = Math.sqrt(d2);
    if (d <= SALTBRIDGE_DIST) contacts.push({ atom: a, d });
  }
  return { nTotal, nHydrophobic, nPolar, contacts };
}

function nearestNamed(probe, allAtoms, ownKey, names, maxD) {
  let best = null;
  for (const a of allAtoms) {
    if (a.resKey === ownKey) continue;
    if (!names.has(a.name)) continue;
    const d = dist3(probe, a);
    if (d <= maxD && (!best || d < best.d)) best = { atom: a, d };
  }
  return best;
}

function nearestElement(probe, allAtoms, ownKey, elems, maxD, excludeResNames = null) {
  let best = null;
  for (const a of allAtoms) {
    if (a.resKey === ownKey) continue;
    if (!elems.has(a.elem)) continue;
    if (excludeResNames && excludeResNames.has(a.resName)) continue;
    const d = dist3(probe, a);
    if (d <= maxD && (!best || d < best.d)) best = { atom: a, d };
  }
  return best;
}

const CARBOXYL_O = new Set(["OD1", "OD2", "OE1", "OE2"]);
const CATION_N = new Set(["NZ", "NH1", "NH2", "NE"]);
const METALS = new Set(["ZN", "FE", "MG", "CA", "MN", "CU", "NI", "CO", "NA", "K"]);

function fmtRes(chain, resSeq, resName) { return `${resName} ${chain}:${resSeq}`; }
function fmtContact(c) {
  return `${c.atom.name}(${c.atom.resName} ${c.atom.chain}:${c.atom.resSeq}) ${c.d.toFixed(2)}Å`;
}

/**
 * Assign protonation states to all titratable residues in a structure.
 *
 * @param {Array|object} residuesOrAtoms  any format accepted by groupResidues()
 * @param {object} [opts]
 * @param {number} [opts.pH=7.0]  target pH for Henderson–Hasselbalch decisions
 * @returns {Array} per-residue assignments
 *   [{ chain, resSeq, resName, state, protonated, charge, pKaEff, reasons: string[] }]
 *   with a `.summary` property { pH, nTitr, byState } (array-with-metadata,
 *   same convention as parseLigands().warnings).
 */
export function assignProtonationStates(residuesOrAtoms, opts = {}) {
  const pH = Number.isFinite(opts.pH) ? opts.pH : 7.0;
  const residues = groupResidues(residuesOrAtoms);
  const allAtoms = [];
  for (const r of residues) {
    const key = `${r.chain}|${r.resSeq}`;
    for (const a of r.atoms) { a.resKey = key; allAtoms.push(a); }
  }

  /** @type {Array} */
  const out = [];
  for (const res of residues) {
    const key = `${res.chain}|${res.resSeq}`;
    const base = {
      chain: res.chain, resSeq: res.resSeq, resName: res.resName,
      state: res.resName, protonated: false, charge: 0, pKaEff: null, reasons: [],
    };
    if (!TITRATABLE.has(res.resName)) { out.push(base); continue; }
    try {
      out.push(assignOne(res, key, allAtoms, pH, base));
    } catch (e) {
      base.reasons.push(`assigner fallback (${e?.message ?? e}); kept native state`);
      out.push(base);
    }
  }
  const byState = {};
  for (const a of out) byState[a.state] = (byState[a.state] ?? 0) + 1;
  out.summary = { pH, nResidues: out.length, nTitr: out.filter((a) => TITRATABLE.has(a.resName)).length, byState };
  return out;
}

function sidechainProbe(res) {
  const names = SIDECHAIN_ATOMS[res.resName] ?? [];
  const sc = res.atoms.filter((a) => names.includes(a.name));
  if (sc.length) return { probe: centroid(sc), nSC: sc.length };
  const ca = atomByName(res, ["CA"]);
  if (ca) return { probe: { x: ca.x, y: ca.y, z: ca.z }, nSC: 0 };
  return { probe: centroid(res.atoms), nSC: 0 };
}

function assignOne(res, key, allAtoms, pH, base) {
  const { probe, nSC } = sidechainProbe(res);
  const env = probeEnvironment(probe, allAtoms, key);
  const burial = clamp01((env.nTotal - 15) / 30); // 0 exposed … 1 deeply buried
  const hydBurial = env.nTotal > 0
    ? burial * (env.nHydrophobic / env.nTotal)
    : 0;
  const tag = fmtRes(res.chain, res.resSeq, res.resName);
  const R = base.reasons;

  if (nSC === 0) {
    R.push(`${tag}: side-chain atoms missing — decision from backbone/exposure proxy only (n8Å=${env.nTotal})`);
  }

  switch (res.resName) {
    case "HIS": {
      const nd1 = atomByName(res, ["ND1"]);
      const ne2 = atomByName(res, ["NE2"]);
      const accND1 = nd1 ? nearestElement(nd1, allAtoms, key, new Set(["O", "N"]), HBOND_DIST) : null;
      const accNE2 = ne2 ? nearestElement(ne2, allAtoms, key, new Set(["O", "N"]), HBOND_DIST) : null;
      const acidBridge = nearestNamed(probe, allAtoms, key, CARBOXYL_O, 4.0);
      const cationNear = nearestNamed(probe, allAtoms, key, CATION_N, 5.0);
      let pKa = MODEL_PKA.HIS + 0.8 * burial - 1.0 * (cationNear ? 1 : 0);
      if (acidBridge) pKa += 1.2;
      base.pKaEff = +pKa.toFixed(2);
      if (pH < pKa - 0.5) {
        base.state = "HIP"; base.protonated = true; base.charge = 1;
        R.push(`${tag}→HIP (+1): pKa_eff ${base.pKaEff} > pH ${pH}` +
          (acidBridge ? `; salt-bridge/H-bond to ${fmtContact(acidBridge)} raises pKa` : "; buried, desolvated") +
          (cationNear ? `; note nearby cation ${fmtContact(cationNear)} (pKa lowered 1.0)` : ""));
      } else if (accND1 && !accNE2) {
        base.state = "HID"; base.charge = 0;
        R.push(`${tag}→HID (H on ND1): H-bond partner ${fmtContact(accND1)} near ND1 can accept; NE2 unpartnered (pKa_eff ${base.pKaEff} ≥ pH ${pH}, neutral)`);
      } else if (accNE2 && !accND1) {
        base.state = "HIE"; base.charge = 0;
        R.push(`${tag}→HIE (H on NE2): H-bond partner ${fmtContact(accNE2)} near NE2 can accept; ND1 unpartnered (pKa_eff ${base.pKaEff} ≥ pH ${pH}, neutral)`);
      } else if (accND1 && accNE2) {
        if (accND1.d <= accNE2.d) {
          base.state = "HID"; base.charge = 0;
          R.push(`${tag}→HID: both nitrogens H-bonded (${fmtContact(accND1)} vs ${fmtContact(accNE2)}); shorter ND1 contact donates`);
        } else {
          base.state = "HIE"; base.charge = 0;
          R.push(`${tag}→HIE: both nitrogens H-bonded (${fmtContact(accND1)} vs ${fmtContact(accNE2)}); shorter NE2 contact donates`);
        }
      } else {
        base.state = "HIE"; base.charge = 0;
        R.push(`${tag}→HIE: no H-bond partner within ${HBOND_DIST}Å of either nitrogen (burial ${(burial).toFixed(2)}, n8Å=${env.nTotal}); epsilon-tautomer default (most populated at pH 7)`);
      }
      return base;
    }
    case "ASP":
    case "GLU": {
      const isAsp = res.resName === "ASP";
      const oNames = isAsp ? ["OD1", "OD2"] : ["OE1", "OE2"];
      const o1 = atomByName(res, [oNames[0]]);
      const o2 = atomByName(res, [oNames[1]]);
      const bridge1 = o1 ? nearestNamed(o1, allAtoms, key, CATION_N, SALTBRIDGE_DIST) : null;
      const bridge2 = o2 ? nearestNamed(o2, allAtoms, key, CATION_N, SALTBRIDGE_DIST) : null;
      const salt = bridge1 ?? bridge2;
      const acidNear = nearestNamed(probe, allAtoms, key, CARBOXYL_O, 5.0);
      let pKa = MODEL_PKA[res.resName] + 4.0 * burial + (acidNear ? 1.0 : 0);
      if (salt) pKa -= 2.5;
      base.pKaEff = +pKa.toFixed(2);
      const protName = isAsp ? "ASH" : "GLH";
      if (pH < pKa) {
        base.state = protName; base.protonated = true; base.charge = 0;
        R.push(`${tag}→${protName} (neutral): buried (n8Å=${env.nTotal}, burial ${burial.toFixed(2)} shifts pKa +${(4.0 * burial).toFixed(1)})` +
          `, no salt bridge within ${SALTBRIDGE_DIST}Å` +
          (acidNear ? `, nearby carboxyl ${fmtContact(acidNear)} (+1.0 charge-repulsion)` : "") +
          ` → pKa_eff ${base.pKaEff} > pH ${pH}`);
      } else {
        base.state = res.resName; base.charge = -1;
        R.push(salt
          ? `${tag}→${res.resName} (−1): salt-bridged to ${fmtContact(salt)} (pKa −2.5) → pKa_eff ${base.pKaEff} ≤ pH ${pH}, deprotonated`
          : `${tag}→${res.resName} (−1): solvent-accessible (burial ${burial.toFixed(2)}, n8Å=${env.nTotal}) → pKa_eff ${base.pKaEff} ≤ pH ${pH}, deprotonated`);
      }
      return base;
    }
    case "LYS":
    case "ARG": {
      const isLys = res.resName === "LYS";
      const tipNames = isLys ? ["NZ"] : ["NH1", "NH2", "NE"];
      const tipP = atomByName(res, tipNames) ?? probe;
      const salt = nearestNamed(tipP, allAtoms, key, CARBOXYL_O, SALTBRIDGE_DIST);
      let pKa = MODEL_PKA[res.resName] - 4.5 * hydBurial;
      if (salt) pKa += 2.0;
      if (env.nPolar === 0 && env.nTotal > 10) pKa -= 1.5; // no polar solvation at all
      base.pKaEff = +pKa.toFixed(2);
      const neutName = isLys ? "LYN" : "ARN";
      if (pKa < pH) {
        base.state = neutName; base.protonated = false; base.charge = 0;
        R.push(`${tag}→${neutName} (neutral): buried hydrophobic (n8Å=${env.nTotal}, hydrophobic fraction ${(env.nTotal ? env.nHydrophobic / env.nTotal : 0).toFixed(2)})` +
          `, no counterion within ${SALTBRIDGE_DIST}Å → pKa_eff ${base.pKaEff} < pH ${pH}, deprotonated`);
      } else {
        base.state = res.resName; base.charge = 1;
        R.push(salt
          ? `${tag}→${res.resName} (+1): salt-bridged to ${fmtContact(salt)} (pKa +2.0) → pKa_eff ${base.pKaEff} ≥ pH ${pH}, charged`
          : `${tag}→${res.resName} (+1): exposed or solvated (burial ${burial.toFixed(2)}, polar neighbors ${env.nPolar}) → pKa_eff ${base.pKaEff} ≥ pH ${pH}, charged`);
      }
      return base;
    }
    case "CYS": {
      const sg = atomByName(res, ["SG"]);
      base.pKaEff = MODEL_PKA.CYS;
      if (!sg) {
        base.state = "CYS"; base.charge = 0;
        R.push(`${tag}→CYS: SG missing, assumed free thiol (no disulfide/metal evidence)`);
        return base;
      }
      let ssPartner = null, metalPartner = null;
      for (const a of allAtoms) {
        if (a.resKey === key) continue;
        if (a.name === "SG" && (a.resName === "CYS" || a.resName === "CYX")) {
          const d = dist3(sg, a);
          if (d <= DISULFIDE_DIST && (!ssPartner || d < ssPartner.d)) {
            ssPartner = { atom: a, d };
          }
        }
        if (METALS.has(a.elem)) {
          const d = dist3(sg, a);
          if (d <= METAL_S_DIST && (!metalPartner || d < metalPartner.d)) {
            metalPartner = { atom: a, d };
          }
        }
      }
      if (ssPartner) {
        base.state = "CYX"; base.charge = 0;
        R.push(`${tag}→CYX (disulfide, 0): SG–SG ${ssPartner.d.toFixed(2)}Å ≤ ${DISULFIDE_DIST}Å to ${fmtContact(ssPartner)}`);
      } else if (metalPartner) {
        base.state = "CYM"; base.charge = -1;
        R.push(`${tag}→CYM (thiolate, −1): SG–${metalPartner.atom.elem} ${metalPartner.d.toFixed(2)}Å ≤ ${METAL_S_DIST}Å to ${fmtContact(metalPartner)}`);
      } else {
        base.state = "CYS"; base.charge = 0;
        R.push(`${tag}→CYS (free thiol, 0): no SG–SG ≤ ${DISULFIDE_DIST}Å, no metal ≤ ${METAL_S_DIST}Å (pKa ${MODEL_PKA.CYS} > pH ${pH}, protonated SH)`);
        base.protonated = true;
      }
      return base;
    }
    default:
      return base;
  }
}

/**
 * Apply assignments onto a heavy.js-style flat atom array (in-place resName
 * rename). Opt-in helper for the structure-load path: with rename=false it
 * only annotates atom.protState / atom.protReasons and changes nothing.
 *
 * NOTE (Phase 3): src/physics/charges.js has no HIE/HID/HIP/ASH/GLH/LYN/ARN
 * entries yet, so renamed atoms fall back to element-based charges there.
 * Keep the returned `chargePatch` (resKey -> charge delta vs native table)
 * for the future charge-table integration.
 *
 * @param {Array} atoms  heavy-style atoms (mutated when rename=true)
 * @param {Array} assignments  output of assignProtonationStates()
 * @param {object} [opts]
 * @param {boolean} [opts.rename=true]
 * @returns {{renamed:number, annotated:number}}
 */
export function applyProtonationStates(atoms, assignments, opts = {}) {
  const rename = opts.rename ?? true;
  const byKey = new Map(assignments.map((a) => [`${a.chain}|${a.resSeq}`, a]));
  let renamed = 0, annotated = 0;
  for (const a of atoms) {
    const key = `${a.chain}|${a.resSeq}`;
    const hit = byKey.get(key);
    if (!hit || !TITRATABLE.has(a.resName)) continue;
    a.protState = hit.state;
    a.protReasons = hit.reasons;
    annotated++;
    if (rename && hit.state !== a.resName) {
      a.resName = hit.state;
      renamed++;
    }
  }
  return { renamed, annotated };
}
