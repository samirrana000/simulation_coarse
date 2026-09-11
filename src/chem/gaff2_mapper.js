/**
 * gaff2_mapper.js — GAFF2-lite atom typer + parameter mapper for small molecules (Phase 2).
 *
 * Covers arbitrary ligands arriving as MOL2/HETATM heavy-atom graphs
 * ({ x,y,z, element } + bond index pairs). Pipeline:
 *
 *   perceiveBondOrders()  — distance-based single/double/triple/aromatic assignment
 *   typeMolecule()         — GAFF2-lite atom types (ca/c3/c2/na/nb/n3/nh/o/oh/os/...)
 *   assignCharges()        — Gasteiger–Marsili PEOE iteration + AM1-BCC-lite
 *                            bond corrections, renormalized to integer net charge
 *   getBondParams() / getAngleParams() / getTorsionParams()
 *                          — equilibrium tables with order-class fallbacks
 *   buildMoleculeParams()  — full topology parameter assignment (convenience)
 *
 * References:
 *   GAFF: Wang et al., J. Comput. Chem. 25, 1157 (2004); GAFF2: amber.dat
 *     (Case et al., AMBER 2020 Manual) — types/values below are rounded
 *     consensus subsets, NOT a byte-exact GAFF2 port (full GAFF2 has ~80 atom
 *     types and ~6000 bond/angle/torsion entries incl. H).
 *   Gasteiger–Marsili PEOE: Gasteiger & Marsili, Tetrahedron Lett. 1979;
 *     (a,b,c) orbital parameters per element (H 7.17/6.24/−0.56 ... I 9.90/7.96/0.96).
 *   AM1-BCC: Jakalian et al., J. Comput. Chem. 21, 132 (2000); 23, 1623 (2002) —
 *     bond-charge corrections (BCC) applied here as a documented lite subset.
 *
 * Notes:
 *   - Heavy atoms only (this engine drops H at parse); implicit-H valence is
 *     accounted for in typing (e.g. degree-2 aromatic N classes).
 *   - No united-atom collapse: every input heavy atom keeps its own type.
 *   - Metals / exotic elements fall back gracefully (type = lowercase element,
 *     generic PEOE parameters, order-class bond/angle defaults) — never throws.
 */

const DEG2RAD = Math.PI / 180;

/* ------------------------------------------------------------------ */
/* Bond-order perception                                              */
/* ------------------------------------------------------------------ */

/** Reference lengths [single, double, triple] in Å per sorted element pair. */
const REF_LENGTHS = {
  "C-C": [1.54, 1.34, 1.20],
  "C-N": [1.47, 1.29, 1.16],
  "C-O": [1.43, 1.21, 1.13],
  "C-S": [1.82, 1.60, 1.53],
  "C-P": [1.84, 1.66, 1.54],
  "C-F": [1.35, 1.30, 1.28],
  "C-CL": [1.77, 1.66, 1.63],
  "C-BR": [1.94, 1.86, 1.83],
  "C-I": [2.14, 2.05, 2.02],
  "N-N": [1.45, 1.25, 1.10],
  "N-O": [1.40, 1.21, 1.06],
  "O-O": [1.48, 1.21, 1.10],
  "S-S": [2.05, 1.90, 1.85],
  "O-S": [1.57, 1.44, 1.40],
  "O-P": [1.63, 1.48, 1.43],
};

const AROM_WINDOW = [1.28, 1.52]; // Å — resonance window for aromatic candidates

function elPair(a, b) {
  const A = String(a).toUpperCase(), B = String(b).toUpperCase();
  return A < B ? `${A}-${B}` : `${B}-${A}`;
}

function bondLength(atoms, i, j) {
  const A = atoms[i], B = atoms[j];
  const ax = A.x ?? A.pos?.[0], ay = A.y ?? A.pos?.[1], az = A.z ?? A.pos?.[2];
  const bx = B.x ?? B.pos?.[0], by = B.y ?? B.pos?.[1], bz = B.z ?? B.pos?.[2];
  return Math.hypot(bx - ax, by - ay, bz - az);
}

function elemOf(a) { return String(a.element ?? a.elem ?? "C").trim().toUpperCase(); }

function adjacency(atoms, bonds) {
  const adj = Array.from({ length: atoms.length }, () => []);
  bonds.forEach(([i, j], b) => {
    if (i == null || j == null || i === j) return;
    if (i < 0 || j < 0 || i >= atoms.length || j >= atoms.length) return;
    adj[i].push({ nbr: j, bond: b });
    adj[j].push({ nbr: i, bond: b });
  });
  return adj;
}

/**
 * Perceive bond orders from heavy-atom distances.
 * Carboxylate resonance (−COO⁻: C with two O at ~1.25Å) is detected
 * explicitly and assigned order 1.5 on both C–O bonds.
 *
 * @param {Array} atoms  [{ x,y,z, element }]
 * @param {Array} bonds  [[i,j], ...] (order ignored on input)
 * @returns {{ orders:number[], aromatic:boolean[], warnings:string[] }}
 *   orders: 1 | 2 | 3 per bond; aromatic: ring-aromatic flags (refined by
 *   findAromaticRings when called via typeMolecule).
 */
export function perceiveBondOrders(atoms, bonds) {
  const orders = new Array(bonds.length).fill(1);
  const aromatic = new Array(bonds.length).fill(false);
  const warnings = [];
  const adj = adjacency(atoms, bonds);

  bonds.forEach(([i, j], b) => {
    if (i == null || j == null || i < 0 || j < 0 || i >= atoms.length || j >= atoms.length || i === j) {
      warnings.push(`bond ${b}: invalid indices [${i},${j}] — treated as single`);
      return;
    }
    const r = bondLength(atoms, i, j);
    const ref = REF_LENGTHS[elPair(elemOf(atoms[i]), elemOf(atoms[j]))];
    if (!ref) { orders[b] = 1; return; } // metal/exotic: generic single
    const dev = ref.map((x) => Math.abs(r - x));
    let best = 0;
    for (let k = 1; k < 3; k++) if (dev[k] < dev[best] - 1e-9) best = k;
    orders[b] = best + 1;
    if (Math.min(...dev) > 0.18) {
      warnings.push(`bond ${b} (${elemOf(atoms[i])}-${elemOf(atoms[j])} ${r.toFixed(2)}Å): ${(best + 1)}x by closest match, deviation ${Math.min(...dev).toFixed(2)}Å — check geometry`);
    }
  });

  // Carboxylate resonance: C with exactly two O neighbors at similar ~1.25Å.
  for (let i = 0; i < atoms.length; i++) {
    if (elemOf(atoms[i]) !== "C") continue;
    const oBonds = adj[i].filter(({ nbr }) => elemOf(atoms[nbr]) === "O");
    if (oBonds.length !== 2) continue;
    const r0 = bondLength(atoms, i, oBonds[0].nbr);
    const r1 = bondLength(atoms, i, oBonds[1].nbr);
    if (Math.abs(r0 - r1) < 0.08 && r0 > 1.18 && r0 < 1.34 && r1 > 1.18 && r1 < 1.34) {
      orders[oBonds[0].bond] = 1.5;
      orders[oBonds[1].bond] = 1.5;
    }
  }
  // Sulfone/sulfoxide S(VI): S–C/S–N bonds are single even when short
  // (hypervalent contraction; e.g. DMSO S–C 1.80Å must not read as multiple).
  for (let i = 0; i < atoms.length; i++) {
    if (elemOf(atoms[i]) !== "S") continue;
    const oDouble = adj[i].some(({ nbr, bond }) =>
      elemOf(atoms[nbr]) === "O" && orders[bond] >= 2);
    if (!oDouble) continue;
    for (const { nbr, bond } of adj[i]) {
      const e = elemOf(atoms[nbr]);
      if ((e === "C" || e === "N") && orders[bond] > 1) {
        orders[bond] = 1;
        warnings.push(`bond ${bond} (S-${e}): reset to single — hypervalent S(VI) carries single S–C/S–N bonds`);
      }
    }
  }
  return { orders, aromatic, warnings };
}

/**
 * Find planar 5-/6-membered aromatic rings (C/N/O/S) and mark ring bonds.
 * Criteria: simple ring, all lengths in AROM_WINDOW, max out-of-plane
 * deviation < 0.12 Å (Newell plane). Mutates the `aromatic` bool array and
 * returns ring atom sets.
 *
 * @param {Array} atoms
 * @param {Array} bonds
 * @param {Array<number>} orders  (unused for detection besides adjacency; kept for API clarity)
 * @param {Array<boolean>} aromatic  per-bond flags to fill
 * @returns {{ rings:number[][], aromaticAtoms:boolean[] }}
 */
export function findAromaticRings(atoms, bonds, orders, aromatic) {
  const adj = adjacency(atoms, bonds);
  const rings = [];
  const seen = new Set();
  const RING_ELEMS = new Set(["C", "N", "O", "S"]);

  for (let start = 0; start < atoms.length; start++) {
    if (!RING_ELEMS.has(elemOf(atoms[start]))) continue;
    const path = [start];
    const inPath = new Set([start]);
    const dfs = (node, parent) => {
      for (const { nbr } of adj[node]) {
        if (nbr === parent) continue;
        if (!RING_ELEMS.has(elemOf(atoms[nbr]))) continue;
        if (nbr === start) {
          if (path.length === 5 || path.length === 6) {
            const key = [...path].sort((a, b) => a - b).join("-");
            if (!seen.has(key)) {
              seen.add(key);
              if (isAromaticRing(atoms, path)) {
                rings.push([...path]);
                for (let t = 0; t < path.length; t++) {
                  const a = path[t], c = path[(t + 1) % path.length];
                  const bi = bonds.findIndex(([p, q]) => (p === a && q === c) || (p === c && q === a));
                  if (bi >= 0) aromatic[bi] = true;
                }
              }
            }
          }
        } else if (!inPath.has(nbr) && path.length < 6) {
          inPath.add(nbr); path.push(nbr); dfs(nbr, node); path.pop(); inPath.delete(nbr);
        }
      }
    };
    dfs(start, -1);
  }
  void orders;
  const aromaticAtoms = new Array(atoms.length).fill(false);
  for (const r of rings) for (const i of r) aromaticAtoms[i] = true;
  return { rings, aromaticAtoms };
}

function ringPlaneDeviation(atoms, ring) {
  // Newell plane normal + max absolute deviation.
  let nx = 0, ny = 0, nz = 0;
  const P = (i) => {
    const a = atoms[i];
    return [a.x ?? a.pos?.[0], a.y ?? a.pos?.[1], a.z ?? a.pos?.[2]];
  };
  const n = ring.length;
  for (let i = 0; i < n; i++) {
    const [x0, y0, z0] = P(ring[i]);
    const [x1, y1, z1] = P(ring[(i + 1) % n]);
    nx += (y0 - y1) * (z0 + z1);
    ny += (z0 - z1) * (x0 + x1);
    nz += (x0 - x1) * (y0 + y1);
  }
  const L = Math.hypot(nx, ny, nz) || 1;
  nx /= L; ny /= L; nz /= L;
  const c = P(ring[0]);
  let mx = 0;
  for (const i of ring) {
    const [x, y, z] = P(i);
    mx = Math.max(mx, Math.abs((x - c[0]) * nx + (y - c[1]) * ny + (z - c[2]) * nz));
  }
  return mx;
}

function isAromaticRing(atoms, ring) {
  for (let t = 0; t < ring.length; t++) {
    const r = bondLength(atoms, ring[t], ring[(t + 1) % ring.length]);
    if (r < AROM_WINDOW[0] || r > AROM_WINDOW[1]) return false;
  }
  return ringPlaneDeviation(atoms, ring) < 0.12;
}

/* ------------------------------------------------------------------ */
/* GAFF2-lite typing                                                  */
/* ------------------------------------------------------------------ */

/**
 * GAFF2-lite type rules (documented subset):
 *   C: ca (aromatic) | c1 (triple) | c2 (double-bonded / carbonyl / carboxyl)
 *      | c3 (sp3, incl. methyl/ammonium-adjacent)
 *   N: na (5-ring aromatic, pyrrole-like, incl. substituted) |
 *      nb (6-ring aromatic, pyridine-like) | nh (amide: single bond to C=O) |
 *      n4 (quaternary / formal +) | n2 (imine, double-bonded) | n3 (sp3 amine)
 *   O: o (double-bonded incl. carboxyl/sulfone) | oh (1 single bond: hydroxyl/
 *      phenol) | os (2 single bonds: ether/ester)
 *   S: s6 (hypervalent: S=O present or degree ≥ 4) | sh (1 single: thiol) |
 *      ss (2 singles: sulfide)    P: p5    Halogens: f/cl/br/i
 *   other: lowercase element (graceful fallback, e.g. "zn").
 *
 * Sybyl hints are honored when present (atom.sybylType | atom.mol2Type |
 * atom.type like "C.ar", "N.am", "O.3"): aromatic flags and amide flags from
 * the hint steer the geometry decision instead of replacing it.
 *
 * @param {Array} atoms  [{ x,y,z, element, ... }] (annotated in place: .gaffType)
 * @param {Array} bonds  [[i,j], ...]
 * @returns {{ types:string[], orders:number[], aromaticAtoms:boolean[], aromaticBonds:boolean[], rings:number[][], warnings:string[] }}
 */
export function typeMolecule(atoms, bonds) {
  const warnings = [];
  const { orders, aromatic } = perceiveBondOrders(atoms, bonds);
  const { rings, aromaticAtoms } = findAromaticRings(atoms, bonds, orders, aromatic);
  const adj = adjacency(atoms, bonds);
  const types = new Array(atoms.length);

  const sybylOf = (a) => String(a.sybylType ?? a.mol2Type ?? (String(a.type ?? "").includes(".") ? a.type : "") ?? "").toLowerCase();
  const hintArom = (a) => /\.ar|arom/.test(sybylOf(a));
  const hintAmide = (a) => /\.am|amidel|n\.am/.test(sybylOf(a));

  for (let i = 0; i < atoms.length; i++) {
    const el = elemOf(atoms[i]);
    const nbrs = adj[i].map(({ nbr }) => nbr);
    const isArom = aromaticAtoms[i] || hintArom(atoms[i]);
    const nHeavy = nbrs.length;
    const maxOrder = Math.max(0, ...adj[i].map(({ bond }) => orders[bond]));
    let t;
    if (el === "C") {
      // Exocyclic double bond (C=O, C=N out of the ring, e.g. xanthine/quinone
      // carbonyl C): sp2 c2 even when the atom sits in an aromatic ring.
      const exoDouble = adj[i].some(({ nbr, bond }) =>
        orders[bond] >= 2 && !(aromatic[bond] && aromaticAtoms[nbr]));
      if (maxOrder >= 3) t = "c1";
      else if (exoDouble) t = "c2";
      else if (isArom) t = "ca";
      else if (maxOrder >= 1.5) t = "c2";
      else if (nHeavy >= 4) t = "c3";
      else if (nHeavy === 3) t = "c2"; // carbocation-like / incomplete valence → sp2
      else t = "c3"; // methyl-like with implicit H
    } else if (el === "N") {
      const formal = Number(atoms[i].formalCharge ?? atoms[i].charge ?? 0);
      if (formal > 0.5 || nHeavy >= 4) t = "n4";
      else if (isArom) {
        const ringSize = rings.find((r) => r.includes(i))?.length ?? 6;
        t = ringSize === 5 ? "na" : "nb";
        if (nHeavy >= 3 && ringSize === 5) t = "na"; // substituted pyrrole-like (e.g. caffeine N–CH3)
      } else if (hintAmide(atoms[i]) || isAmideN(i, atoms, adj, orders)) t = "nh";
      else if (maxOrder >= 3) t = "n2";
      else if (maxOrder >= 2) t = "n2";
      else t = "n3";
    } else if (el === "O") {
      if (maxOrder >= 1.5) t = "o";
      else if (nHeavy >= 2) t = "os";
      else t = "oh";
    } else if (el === "S") {
      const hasDoubleO = adj[i].some(({ nbr, bond }) => elemOf(atoms[nbr]) === "O" && orders[bond] >= 2);
      if (hasDoubleO || nHeavy >= 4) t = "s6";
      else if (nHeavy <= 1) t = "sh";
      else t = "ss";
    } else if (el === "P") {
      t = "p5";
    } else if (["F", "CL", "BR", "I"].includes(el)) {
      t = el.toLowerCase();
    } else {
      t = el.toLowerCase(); // graceful fallback for metals / exotic elements
      warnings.push(`atom ${i} (${el}): no GAFF2-lite type — fallback "${t}" with generic parameters`);
    }
    types[i] = t;
    atoms[i].gaffType = t;
  }
  // Carboxylate O: both 1.5-order O are type "o" (perceive step already set
  // order 1.5 → maxOrder ≥ 1.5 → "o"). Nothing further needed.
  return { types, orders, aromaticAtoms, aromaticBonds: aromatic, rings, warnings };
}

function isAmideN(i, atoms, adj, orders) {
  // N single-bonded to a C that carries a double-bonded O (C=O).
  for (const { nbr, bond } of adj[i]) {
    if (elemOf(atoms[nbr]) !== "C" || orders[bond] !== 1) continue;
    const cAdj = adj[nbr];
    if (cAdj.some(({ nbr: o, bond: b }) => elemOf(atoms[o]) === "O" && orders[b] >= 2)) return true;
  }
  return false;
}

/* ------------------------------------------------------------------ */
/* Gasteiger–Marsili PEOE + AM1-BCC-lite charges                       */
/* ------------------------------------------------------------------ */

/** PEOE (a, b, c) orbital parameters per element (Gasteiger & Marsili 1979). */
export const PEOE_PARAMS = Object.freeze({
  H: [7.17, 6.24, -0.56], C: [7.98, 9.18, 1.88], N: [11.54, 10.82, 1.36],
  O: [14.18, 12.92, 1.39], F: [14.66, 13.85, 2.31], P: [8.90, 8.01, 1.06],
  S: [10.14, 9.13, 1.38], CL: [11.00, 9.69, 1.35], BR: [10.08, 8.47, 1.16],
  I: [9.90, 7.96, 0.96],
});
const PEOE_GENERIC = [8.0, 8.0, 1.0]; // fallback for metals / exotic elements

/**
 * AM1-BCC-lite bond corrections: key "typeA|typeB|orderClass" -> Δq moved
 * from the less- to the more-electronegative partner. Full AM1-BCC has ~100+
 * BCCs (Jakalian et al. 2000/2002); this lite subset covers the library
 * chemistry (carbonyl, hydroxyl, ether, amine, amide, sulfone, halide).
 */
export const BCC_LITE = Object.freeze({
  "c2|o|2": 0.12,    // carbonyl C=O
  "c2|o|1.5": 0.10,  // carboxylate resonance
  "ca|oh|1": 0.08,   // phenol
  "c3|oh|1": 0.06,   // alcohol
  "c2|os|1": 0.08,   // ester C–O
  "c3|os|1": 0.04,   // ether
  "ca|os|1": 0.05,
  "c3|n3|1": 0.03,   // amine
  "ca|na|ar": 0.04, "ca|nb|ar": 0.04, "ca|n3|1": 0.04,
  "c2|nh|1": 0.10, "c2|n2|2": 0.08,   // amide / imine
  "s6|o|2": 0.15,    // sulfone/sulfoxide S=O
  "c3|f|1": 0.06, "c3|cl|1": 0.03, "ca|cl|1": 0.03, "c3|br|1": 0.02,
});

function orderClass(order, aromatic) {
  if (aromatic) return "ar";
  if (order === 1.5) return "1.5";
  return String(Math.round(order));
}

/**
 * Assign partial charges: Gasteiger–Marsili electronegativity equalization
 * iterated to convergence, then AM1-BCC-lite bond corrections, then exact
 * renormalization to the integer net charge.
 *
 * @param {Array} atoms  typed atoms (atom.gaffType set; call typeMolecule first —
 *   untyped atoms are typed on the fly from element with a warning-free fallback)
 * @param {Array} [bonds=[]]  [[i,j], ...]
 * @param {object} [opts]
 * @param {number} [opts.totalCharge]  default: Σ atom.formalCharge ?? 0
 * @param {number} [opts.maxIter=100]  PEOE iterations (damping 2^-n, converges ≪ this)
 * @param {number} [opts.tol=1e-4]     max |Δq| stopping criterion
 * @param {boolean} [opts.writeBack=true]  store result on atom.charge
 * @returns {{ charges:Float64Array, iterations:number, maxDelta:number, netCharge:number }}
 */
export function assignCharges(atoms, bonds = [], opts = {}) {
  const n = atoms.length;
  const charges = new Float64Array(n);
  if (n === 0) return { charges, iterations: 0, maxDelta: 0, netCharge: 0 };
  const types = atoms.map((a, i) => {
    if (!a.gaffType) {
      try {
        const r = typeMolecule(atoms, bonds);
        return r.types[i];
      } catch { return elemOf(a).toLowerCase(); }
    }
    return a.gaffType;
  });
  void types;
  const params = atoms.map((a) => PEOE_PARAMS[elemOf(a)] ?? PEOE_GENERIC);
  const maxIter = opts.maxIter ?? 100;
  const tol = opts.tol ?? 1e-4;
  let iterations = 0, maxDelta = Infinity;

  const chi = (k, q) => params[k][0] + params[k][1] * q + params[k][2] * q * q;

  for (let it = 1; it <= maxIter; it++) {
    iterations = it;
    const damp = Math.pow(0.5, it); // geometric damping: later shells perturb less
    maxDelta = 0;
    for (const [rawI, rawJ] of bonds) {
      const i = rawI, j = rawJ;
      if (i == null || j == null || i < 0 || j < 0 || i >= n || j >= n || i === j) continue;
      const dChi = chi(j, charges[j]) - chi(i, charges[i]);
      if (dChi === 0) continue;
      // Charge flows toward higher electronegativity: dq > 0 moves +charge
      // from j to i (i becomes more positive, j more negative).
      const dq = (dChi / (params[i][1] + params[j][1])) * damp;
      charges[i] += dq;
      charges[j] -= dq;
      const ad = Math.abs(dq);
      if (ad > maxDelta) maxDelta = ad;
    }
    if (maxDelta < tol) break;
  }

  // AM1-BCC-lite bond corrections (antisymmetric shifts, net-conserving).
  const aromFlags = (() => {
    try {
      const adj = adjacency(atoms, bonds);
      void adj;
      const { orders, aromatic } = perceiveBondOrders(atoms, bonds);
      findAromaticRings(atoms, bonds, orders, aromatic);
      const perAtomOrder = atoms.map(() => ({ order: 1, arom: false }));
      bonds.forEach(([i, j], b) => {
        if (i < 0 || j < 0 || i >= n || j >= n) return;
        perAtomOrder[b] = { order: orders[b], arom: aromatic[b] };
      });
      return perAtomOrder;
    } catch { return bonds.map(() => ({ order: 1, arom: false })); }
  })();
  const chiOf = (k) => {
    const p = params[k];
    return p[0] + p[1] * charges[k] + p[2] * charges[k] * charges[k];
  };
  bonds.forEach(([rawI, rawJ], b) => {
    const i = rawI, j = rawJ;
    if (i == null || j == null || i < 0 || j < 0 || i >= n || j >= n || i === j) return;
    const ti = String(atoms[i].gaffType ?? elemOf(atoms[i]).toLowerCase());
    const tj = String(atoms[j].gaffType ?? elemOf(atoms[j]).toLowerCase());
    const oc = orderClass(aromFlags[b]?.order ?? 1, aromFlags[b]?.arom ?? false);
    const key = `${ti}|${tj}|${oc}`;
    const rev = `${tj}|${ti}|${oc}`;
    const corr = BCC_LITE[key] ?? BCC_LITE[rev];
    if (!corr) return;
    // Move electrons toward the more electronegative partner (higher chi).
    const s = chiOf(j) >= chiOf(i) ? 1 : -1;
    charges[i] -= s * corr * 0.5;
    charges[j] += s * corr * 0.5;
  });

  // Exact renormalization to integer net charge.
  let totalCharge = opts.totalCharge;
  if (!Number.isFinite(totalCharge)) {
    totalCharge = 0;
    for (const a of atoms) totalCharge += Number(a.formalCharge ?? 0);
    totalCharge = Math.round(totalCharge);
  }
  let sum = 0;
  for (let i = 0; i < n; i++) sum += charges[i];
  const shift = (totalCharge - sum) / n;
  for (let i = 0; i < n; i++) charges[i] += shift;

  if (opts.writeBack ?? true) {
    for (let i = 0; i < n; i++) atoms[i].charge = charges[i];
  }
  return { charges, iterations, maxDelta, netCharge: totalCharge };
}

/* ------------------------------------------------------------------ */
/* Equilibrium parameter tables                                       */
/* ------------------------------------------------------------------ */

/** Bond table: "t1|t2|orderClass" -> { r0 Å, k kcal/mol/Å² }. Symmetric. */
export const GAFF_BONDS = Object.freeze({
  "ca|ca|ar": { r0: 1.400, k: 469 }, "ca|ca|1": { r0: 1.400, k: 469 },
  "ca|c3|1": { r0: 1.510, k: 317 }, "ca|ca|2": { r0: 1.400, k: 469 },
  "c3|c3|1": { r0: 1.526, k: 310 }, "c3|c2|1": { r0: 1.510, k: 310 },
  "c2|c2|2": { r0: 1.340, k: 620 }, "c2|c2|1": { r0: 1.460, k: 410 },
  "c2|o|2": { r0: 1.229, k: 570 }, "c2|o|1.5": { r0: 1.250, k: 570 },
  "c2|oh|1": { r0: 1.360, k: 450 }, "c2|os|1": { r0: 1.360, k: 450 },
  "c3|oh|1": { r0: 1.430, k: 320 }, "c3|os|1": { r0: 1.430, k: 320 },
  "ca|oh|1": { r0: 1.360, k: 320 }, "ca|os|1": { r0: 1.360, k: 320 },
  "c3|n3|1": { r0: 1.470, k: 337 }, "c3|nh|1": { r0: 1.450, k: 337 },
  "ca|na|ar": { r0: 1.380, k: 400 }, "ca|nb|ar": { r0: 1.340, k: 400 },
  "c2|na|ar": { r0: 1.380, k: 400 }, "c2|nb|ar": { r0: 1.380, k: 400 },
  "c2|na|1": { r0: 1.400, k: 400 }, "c2|nb|1": { r0: 1.400, k: 400 },
  "c2|ca|ar": { r0: 1.400, k: 410 }, "c2|ca|1": { r0: 1.460, k: 410 },
  "na|c3|1": { r0: 1.470, k: 337 }, "nb|c3|1": { r0: 1.470, k: 337 },
  "ca|n3|1": { r0: 1.450, k: 337 }, "ca|nh|1": { r0: 1.400, k: 400 },
  "c2|nh|1": { r0: 1.335, k: 490 }, "c2|n2|2": { r0: 1.290, k: 520 },
  "c2|n2|1": { r0: 1.400, k: 400 }, "n2|n2|2": { r0: 1.250, k: 500 },
  "c3|s6|1": { r0: 1.810, k: 237 }, "c3|ss|1": { r0: 1.810, k: 237 },
  "s6|o|2": { r0: 1.440, k: 550 }, "s6|o|1": { r0: 1.570, k: 350 },
  "c3|f|1": { r0: 1.350, k: 400 }, "c3|cl|1": { r0: 1.760, k: 250 },
  "ca|cl|1": { r0: 1.740, k: 250 }, "c3|br|1": { r0: 1.940, k: 220 },
  "c3|i|1": { r0: 2.100, k: 200 }, "c2|p5|1": { r0: 1.800, k: 230 },
  "o|p5|1": { r0: 1.600, k: 300 }, "o|p5|2": { r0: 1.480, k: 450 },
});

/** Order-class bond defaults (GAFF consensus rounded). */
export const BOND_ORDER_DEFAULT = Object.freeze({
  "1": { r0: 1.500, k: 300 }, "ar": { r0: 1.400, k: 469 },
  "1.5": { r0: 1.400, k: 469 }, "2": { r0: 1.330, k: 570 }, "3": { r0: 1.200, k: 620 },
});

/** Angle table: "t1|t2|t3" (t2 apex, reversal-symmetric) -> { th0 deg, k }. */
export const GAFF_ANGLES = Object.freeze({
  "ca|ca|ca": { th0: 120.0, k: 60 }, "ca|ca|c3": { th0: 120.0, k: 60 },
  "ca|ca|oh": { th0: 120.0, k: 60 },
  "ca|ca|na": { th0: 120.0, k: 60 }, "ca|ca|nb": { th0: 120.0, k: 60 },
  "nb|ca|nb": { th0: 120.0, k: 60 }, "na|ca|ca": { th0: 120.0, k: 60 },
  "ca|na|ca": { th0: 122.0, k: 60 }, "ca|nb|ca": { th0: 118.0, k: 60 },
  "ca|na|c3": { th0: 122.0, k: 60 }, "ca|nb|c3": { th0: 120.0, k: 60 },
  "nb|c2|nb": { th0: 120.0, k: 70 }, "na|c2|na": { th0: 120.0, k: 70 },
  "o|c2|nb": { th0: 122.0, k: 70 }, "o|c2|na": { th0: 122.0, k: 70 },
  "o|c2|ca": { th0: 120.0, k: 70 }, "nb|c2|ca": { th0: 115.0, k: 60 },
  "na|c2|ca": { th0: 115.0, k: 60 },
  "ca|nb|c2": { th0: 120.0, k: 60 }, "c2|nb|c3": { th0: 120.0, k: 60 },
  "c2|nb|c2": { th0: 122.0, k: 70 }, "c2|ca|ca": { th0: 120.0, k: 60 },
  "c2|ca|na": { th0: 120.0, k: 60 }, "c2|ca|nb": { th0: 120.0, k: 60 },
  "na|ca|na": { th0: 120.0, k: 60 },
  "ca|ca|cl": { th0: 120.0, k: 50 }, "ca|ca|br": { th0: 120.0, k: 50 },
  "ca|ca|f": { th0: 120.0, k: 50 },
  "ca|ca|na": { th0: 120.0, k: 60 }, "na|ca|nb": { th0: 120.0, k: 60 },
  "ca|na|c3": { th0: 122.0, k: 60 }, "c3|c3|c3": { th0: 112.0, k: 58 },
  "c3|c3|oh": { th0: 109.5, k: 50 }, "c3|c3|os": { th0: 109.5, k: 50 },
  "c3|c3|n3": { th0: 111.0, k: 60 }, "c3|os|c3": { th0: 112.0, k: 50 },
  "c2|c2|c2": { th0: 120.0, k: 60 }, "o|c2|oh": { th0: 122.0, k: 70 },
  "o|c2|os": { th0: 122.0, k: 70 }, "o|c2|nh": { th0: 122.0, k: 70 },
  "o|c2|c3": { th0: 120.0, k: 70 }, "nh|c2|c3": { th0: 115.0, k: 60 },
  "o|c2|o": { th0: 126.0, k: 80 }, "c3|c2|nh": { th0: 115.0, k: 60 },
  "c3|s6|c3": { th0: 100.0, k: 50 }, "c3|s6|o": { th0: 106.0, k: 50 },
  "o|s6|o": { th0: 119.0, k: 50 }, "c3|c3|f": { th0: 109.5, k: 50 },
  "c3|c3|cl": { th0: 109.5, k: 50 }, "c3|ss|c3": { th0: 100.0, k: 50 },
});

/** Central-type angle defaults: sp2/aromatic 120°, sp3 109.5°. */
function angleDefault(centerType) {
  const t = String(centerType);
  if (t === "ca" || t === "c2" || t === "o" || t === "n2" || t === "nb" || t === "na" || t === "nh" || t === "o") {
    return { th0: 120.0 * DEG2RAD, k: 60.0, source: "fallback-sp2" };
  }
  return { th0: 109.5 * DEG2RAD, k: 50.0, source: "fallback-sp3" };
}

/**
 * Torsion table: "t1|t2|t3|t4" with X wildcards -> Fourier terms.
 * Energy: U = Σ (Vn/2)[1 + cos(n φ − γ)]. γ stored in radians on output.
 */
export const GAFF_TORSIONS = Object.freeze({
  "X|ca|ca|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|ca|na|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|na|ca|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|ca|nb|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|nb|ca|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|c2|ca|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|ca|c2|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|c2|na|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|na|c2|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|c2|nb|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|nb|c2|X": [{ Vn: 3.0, n: 2, gamma: 180 }],
  "X|c3|na|X": [{ Vn: 1.4, n: 3, gamma: 0 }],
  "X|na|c3|X": [{ Vn: 1.4, n: 3, gamma: 0 }],
  "X|c3|nb|X": [{ Vn: 1.4, n: 3, gamma: 0 }],
  "X|nb|c3|X": [{ Vn: 1.4, n: 3, gamma: 0 }],
  "X|c3|ss|X": [{ Vn: 1.0, n: 3, gamma: 0 }],
  "X|ss|c3|X": [{ Vn: 1.0, n: 3, gamma: 0 }],
  "X|c3|s6|X": [{ Vn: 1.0, n: 3, gamma: 0 }],
  "X|s6|c3|X": [{ Vn: 1.0, n: 3, gamma: 0 }],
  "X|nh|c2|X": [{ Vn: 8.0, n: 2, gamma: 180 }],
  "X|c3|c3|X": [{ Vn: 1.4, n: 3, gamma: 0 }],
  "X|c3|ca|X": [{ Vn: 1.4, n: 3, gamma: 0 }],
  "X|c2|c2|X": [{ Vn: 6.0, n: 2, gamma: 180 }],
  "X|c2|nh|X": [{ Vn: 8.0, n: 2, gamma: 180 }], // amide torsion barrier
  "X|c3|n3|X": [{ Vn: 1.4, n: 3, gamma: 0 }],
  "X|c3|oh|X": [{ Vn: 1.0, n: 3, gamma: 0 }],
  "X|c3|os|X": [{ Vn: 1.2, n: 3, gamma: 0 }],
  "X|ca|oh|X": [{ Vn: 2.0, n: 2, gamma: 180 }],
});
const TORSION_DEFAULT = [{ Vn: 1.4, n: 3, gamma: 0 }];

function bondTableKey(t1, t2, oc) {
  const a = `${t1}|${t2}|${oc}`, b = `${t2}|${t1}|${oc}`;
  return GAFF_BONDS[a] ? a : (GAFF_BONDS[b] ? b : null);
}

/**
 * Bond parameters for a GAFF2-lite typed pair.
 * @param {string} t1  GAFF2-lite type (or element symbol — normalized inside)
 * @param {string} t2
 * @param {number} [order=1]  1 | 1.5 | 2 | 3 (aromatic bonds: pass order=1.5 or aromatic=true)
 * @param {boolean} [aromatic=false]
 * @returns {{ r0:number, k:number, source:string }}  r0 Å, k kcal/mol/Å²
 */
export function getBondParams(t1, t2, order = 1, aromatic = false) {
  const a = String(t1 ?? "c3").toLowerCase(), b = String(t2 ?? "c3").toLowerCase();
  const oc = orderClass(order, aromatic);
  const key = bondTableKey(a, b, oc)
    ?? bondTableKey(a, b, "1") ?? bondTableKey(a, b, "ar");
  if (key) {
    const hit = GAFF_BONDS[key];
    return { r0: hit.r0, k: hit.k, source: `table:${key}` };
  }
  const d = BOND_ORDER_DEFAULT[oc] ?? BOND_ORDER_DEFAULT["1"];
  return { r0: d.r0, k: d.k, source: `fallback-order:${oc}` };
}

/**
 * Angle parameters for a GAFF2-lite typed triple (t2 = apex).
 * @returns {{ th0:number, k:number, source:string }}  th0 radians
 */
export function getAngleParams(t1, t2, t3) {
  const a = String(t1 ?? "c3").toLowerCase();
  const b = String(t2 ?? "c3").toLowerCase();
  const c = String(t3 ?? "c3").toLowerCase();
  const fwd = `${a}|${b}|${c}`, rev = `${c}|${b}|${a}`;
  const hit = GAFF_ANGLES[fwd] ?? GAFF_ANGLES[rev];
  if (hit) return { th0: hit.th0 * DEG2RAD, k: hit.k, source: `table:${GAFF_ANGLES[fwd] ? fwd : rev}` };
  return angleDefault(b);
}

/**
 * Proper-dihedral Fourier terms for a GAFF2-lite typed quadruple.
 * Always returns an array (never empty); wildcard X matching, then default.
 * @returns {{ terms:Array<{Vn:number,n:number,gamma:number}>, source:string }}  gamma radians
 */
export function getTorsionParams(t1, t2, t3, t4) {
  const q = [t1, t2, t3, t4].map((t) => String(t ?? "X").toLowerCase());
  const cands = [
    `${q[0]}|${q[1]}|${q[2]}|${q[3]}`,
    `X|${q[1]}|${q[2]}|X`,
    `X|${q[1]}|${q[2]}|${q[3]}`,
    `${q[0]}|${q[1]}|${q[2]}|X`,
  ];
  for (const key of cands) {
    if (GAFF_TORSIONS[key]) {
      return {
        terms: GAFF_TORSIONS[key].map((t) => ({ Vn: t.Vn, n: t.n, gamma: t.gamma * DEG2RAD })),
        source: `table:${key}`,
      };
    }
  }
  return {
    terms: TORSION_DEFAULT.map((t) => ({ Vn: t.Vn, n: t.n, gamma: t.gamma * DEG2RAD })),
    source: "fallback-default",
  };
}

/**
 * Assign parameters to a full molecular topology (convenience wrapper).
 * Angles are enumerated from adjacency; proper torsions from bonded paths
 * i–j–k–l (each unique central bond once, both directions deduped).
 *
 * @param {Array} atoms  typed atoms (types via typeMolecule when absent)
 * @param {Array} bonds  [[i,j], ...]
 * @returns {{ bondParams:Array, angleParams:Array, torsionParams:Array,
 *   counts:{bonds:number,angles:number,torsions:number},
 *   fallbacks:{bonds:number,angles:number,torsions:number} }}
 */
export function buildMoleculeParams(atoms, bonds) {
  let types = atoms.map((a) => a.gaffType);
  let orders = null, aromatic = null;
  if (types.some((t) => !t)) {
    const r = typeMolecule(atoms, bonds);
    types = r.types; orders = r.orders; aromatic = r.aromaticBonds;
  } else {
    const p = perceiveBondOrders(atoms, bonds);
    orders = p.orders; aromatic = p.aromatic;
    try { findAromaticRings(atoms, bonds, orders, aromatic); } catch { /* keep flags */ }
  }
  const adj = adjacency(atoms, bonds);
  const bondParams = [];
  let fbB = 0;
  bonds.forEach(([i, j], b) => {
    if (i < 0 || j < 0 || i >= atoms.length || j >= atoms.length || i === j) return;
    const p = getBondParams(types[i], types[j], orders[b], aromatic[b]);
    if (!p.source.startsWith("table:")) fbB++;
    bondParams.push({ i, j, order: orders[b], aromatic: aromatic[b], ...p });
  });

  // NOTE: angle entries use kAngle for the force constant (not k) because
  // k is already the third atom index — same flat [i,j,k,th0] spirit as
  // ligand.js, without the key collision.
  const angleParams = [];
  let fbA = 0;
  for (let j = 0; j < atoms.length; j++) {
    const nb = adj[j].map(({ nbr }) => nbr);
    for (let a = 0; a < nb.length; a++) {
      for (let b = a + 1; b < nb.length; b++) {
        const p = getAngleParams(types[nb[a]], types[j], types[nb[b]]);
        if (!p.source.startsWith("table:")) fbA++;
        angleParams.push({ i: nb[a], j, k: nb[b], th0: p.th0, kAngle: p.k, source: p.source });
      }
    }
  }

  const torsionParams = [];
  let fbT = 0;
  const seenT = new Set();
  for (let j = 0; j < atoms.length; j++) {
    for (const { nbr: i } of adj[j]) {
      for (const { nbr: k } of adj[j]) {
        if (k === i) continue;
        for (const { nbr: l } of adj[k]) {
          if (l === i || l === j) continue;
          const key = [i, j, k, l].join("-");
          const rev = [l, k, j, i].join("-");
          if (seenT.has(key) || seenT.has(rev)) continue;
          seenT.add(key);
          const p = getTorsionParams(types[i], types[j], types[k], types[l]);
          if (!p.source.startsWith("table:")) fbT++;
          torsionParams.push({ i, j, k, l, ...p });
        }
      }
    }
  }
  return {
    bondParams, angleParams, torsionParams,
    counts: { bonds: bondParams.length, angles: angleParams.length, torsions: torsionParams.length },
    fallbacks: { bonds: fbB, angles: fbA, torsions: fbT },
  };
}
