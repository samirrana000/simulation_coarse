/**
 * heavy/parse.js — PDB text to a flat heavy-atom system with classified hetero
 * groups. Heavy-mode input only; nothing here is shared with CG mode.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10).
 *
 *   SOLVENT          water residue names, dropped before anything else
 *   countWarnings    A06 warning-count helper for the parse result
 *   parseHeavy       the ATOM/HETATM loop: altloc de-duplication, metal flag,
 *                    hetero-group grouping, malformed-line warnings
 *   elementFromName  PDB columns 77-78 are often blank; this recovers the
 *                    element from the atom name, including the ATOM CA (alpha
 *                    carbon) vs HETATM CA (calcium) distinction
 */
import { METAL_ELEMENT } from "../ff-params.js";
import { parseAltLoc, parseOccupancy, shouldReplaceAltloc } from "../pdb_altloc.js";

const SOLVENT = new Set([
  "HOH", "WAT", "H2O", "DOD", "HHO", "TIP", "TIP3", "TIP3P", "SPC", "SPCE", "SOL",
]);

/**
 * Count validation warnings helper (A06).
 * @param {string[]|object} warnings
 * @returns {number}
 */
export function countWarnings(warnings) {
  if (Array.isArray(warnings)) return warnings.length;
  if (warnings && Array.isArray(warnings.warnings)) return warnings.warnings.length;
  return 0;
}

/**
 * Parse a PDB text into a flat heavy-atom system with classified hetero groups.
 *
 * Input validation / warnings (A06): malformed lines are skipped with
 * `warnings` + `console.warn`; result includes `.warnings` (string[]).
 * Use `countWarnings(result.warnings)` to count.
 *
 * @param {string} pdbText
 * @returns {{atoms: Array, n: number, bySerial: Map, heteroGroups: Array, warnings: string[]}}
 */
export function parseHeavy(pdbText) {
  const atoms = [];
  // Stage-4 altloc cleaner (default-ON, bit-identical on clean files): per
  // (chain, resSeq, iCode, atom name) key keep the highest-occupancy altloc
  // ('A' on tie, else first); zero-occupancy duplicates dropped. Shared rule
  // with parseCa/parseLigands — see src/pdb_altloc.js. Fixes the S7 NaN root
  // cause (coincident duplicate-altloc atoms) at the input, in place.
  const seen = new Map();
  const bySerial = new Map();
  const heteroGroupMap = new Map();
  const warnings = [];
  let n = 0;

  for (const line of pdbText.split(/\r?\n/)) {
    const rec = line.slice(0, 6);
    if (rec !== "ATOM  " && rec !== "HETATM") continue;
    const atomName = line.slice(12, 16).trim();
    const resName = line.slice(17, 20).trim();
    if (SOLVENT.has(resName)) continue;

    const chain = (line.charAt(21) || " ").trim() || "_";
    const resSeq = parseInt(line.slice(22, 26), 10);
    const iCode = line.charAt(26).trim();
    const x = parseFloat(line.slice(30, 38));
    const y = parseFloat(line.slice(38, 46));
    const z = parseFloat(line.slice(46, 54));
    if (Number.isNaN(resSeq) || Number.isNaN(x + y + z)) {
      const msg = `parseHeavy: malformed line skipped (resSeq=${resSeq} x=${x}) line="${line.slice(0, 66).trim()}"`;
      warnings.push(msg); console.warn(msg);
      continue;
    }

    let element = line.slice(76, 78).trim().toUpperCase();
    if (!element) element = elementFromName(atomName, rec);
    if (!element || element === "H") continue;

    const serial = parseInt(line.slice(6, 11), 10);
    const key = `${chain}|${resSeq}|${iCode}|${atomName}`;
    const altLoc = parseAltLoc(line);
    const occ = parseOccupancy(line);
    if (seen.has(key)) {
      const prev = seen.get(key);
      if (shouldReplaceAltloc(prev, { occupancy: occ, altLoc })) {
        const msg = `parseHeavy: duplicate atom ${key} altloc '${prev.altLoc || " "}'→'${altLoc || " "}' replaced (occ ${prev.occupancy}→${occ})`;
        warnings.push(msg); console.warn(msg);
        const idx = prev.idx;
        const oldAtom = atoms[idx];
        bySerial.delete(oldAtom.serial);
        const isProtein = rec === "ATOM  ";
        const isHetero = !isProtein;
        const isMetal = !!METAL_ELEMENT[element];
        const heteroKey = isHetero ? `${chain}|${resSeq}|${resName}` : null;
        const a = {
          x, y, z, element, atomName, resName, chain, resSeq, serial,
          isProtein, isMetal, isHetero, isWater: false, heteroKey,
        };
        atoms[idx] = a;
        bySerial.set(serial, a);
        seen.set(key, { idx, occupancy: occ, altLoc });
        if (isHetero) {
          const g = heteroGroupMap.get(heteroKey);
          if (g) {
            const pos = g.atomIndices.indexOf(idx);
            if (pos >= 0) g.elements[pos] = element;
            if (isMetal) g.isMetal = true;
          }
        }
      } else {
        const msg = `parseHeavy: duplicate atom ${key} skipped (altLoc '${altLoc || " "}' occ ${occ})`;
        warnings.push(msg); console.warn(msg);
      }
      continue;
    }
    seen.set(key, { idx: atoms.length, occupancy: occ, altLoc });

    const isProtein = rec === "ATOM  ";
    const isHetero = !isProtein;
    const isMetal = !!METAL_ELEMENT[element];
    const heteroKey = isHetero ? `${chain}|${resSeq}|${resName}` : null;
    const a = {
      x, y, z, element, atomName, resName, chain, resSeq, serial,
      isProtein, isMetal, isHetero, isWater: false, heteroKey,
    };
    atoms.push(a);
    bySerial.set(serial, a);
    n++;

    if (isHetero) {
      let g = heteroGroupMap.get(heteroKey);
      if (!g) {
        g = { key: heteroKey, resName, chain, resSeq, atomIndices: [], isMetal: false, elements: [] };
        heteroGroupMap.set(heteroKey, g);
      }
      g.atomIndices.push(atoms.length - 1);
      g.elements.push(element);
      if (isMetal) g.isMetal = true;
    }
  }

  if (atoms.length === 0) throw new Error("No heavy atoms found in PDB file.");

  const heteroGroups = [...heteroGroupMap.values()].map((g) => ({
    key: g.key,
    resName: g.resName,
    chain: g.chain,
    resSeq: g.resSeq,
    atomIndices: g.atomIndices,
    isMetal: g.isMetal,
    element: g.elements[0] ?? null,
  }));

  if (warnings.length) console.warn(`[parseHeavy] ${warnings.length} warning(s) total`);
  return { atoms, n, bySerial, heteroGroups, warnings };
}

function elementFromName(name, rec) {
  if (!name) return null;
  const m = name.match(/^[0-9]*([A-Za-z]{1,2})/);
  if (!m) return null;
  let el = m[1].toUpperCase();
  // Guard: ATOM CA is alpha carbon, not calcium (HETATM CA is calcium)
  if (rec === "ATOM  " && name.trim() === "CA") return "C";
  if (el.length === 2) {
    const two = new Set(["CL", "BR", "ZN", "FE", "MG", "CA", "CU", "MN", "NI", "CO", "NA", "K", "SE", "SI", "AL"]);
    if (two.has(el)) return el;
    return el[0];
  }
  return el;
}
