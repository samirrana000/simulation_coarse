/**
 * rotbonds.js — Rotatable-bond detection for ligands (Stage-4).
 *
 * Zero-dependency, pure functions. Input is the ligand heavy-atom graph used
 * everywhere in this repo (`{ x, y, z, element }` atoms + `[[i, j]]` bond
 * index pairs, as produced by `parseLigands`, MOL2 parsing, or the heavy
 * covalent topology). Bond entries may carry an optional third element
 * `[i, j, order]`; parallel `opts.orders` / `opts.aromatic` arrays (as from
 * `gaff2_mapper.perceiveBondOrders`) are honored when supplied.
 *
 * Rotatable-bond rule (documented):
 *  1. Single bond only: order ≈ 1 (within 0.15) and not aromatic. When no
 *     order info is available the bond is assumed single — safe here because
 *     ring bonds (incl. aromatic) are excluded by rule 2 and terminal
 *     carbonyls are excluded by rule 4 (degree-1 oxygen side).
 *  2. Non-ring: the bond must not lie on any cycle. Tested per bond by BFS:
 *     removing the bond itself, the two endpoints must become disconnected.
 *     This covers 5/6-membered aromatic rings and any other ring size.
 *  3. Non-amide: C–N bonds where the carbon carries a carbonyl oxygen are
 *     excluded (amide torsion barrier ~8 kcal/mol, not freely rotatable).
 *     Detection order: (a) gaffType 'nh' on the nitrogen when typed;
 *     (b) explicit orders — C has an order ≥ 1.5 bond to O; (c) geometry
 *     fallback — C has an oxygen neighbor within 1.35 Å (carbonyl length).
 *  4. Both sides substituted: degree(j) > 1 and degree(k) > 1 (heavy atoms).
 *     This excludes terminal methyl/methoxy/hydroxyl/halogen rotors whose
 *     rotation is degenerate or trivial.
 *  5. Symmetric-top exclusion: if all non-partner neighbors of either side
 *     are leaves (degree 1) AND there are ≥ 3 of them (tert-butyl-like
 *     C(CH₃)₃, S(=O)₃ sulfonate head, …), the bond is excluded — rotation
 *     permutes indistinguishable methyl/oxygen leaves. Methyl itself is
 *     already excluded by rule 4; this rule only fires for ≥ 3-leaf fans.
 *
 * Torsion quadruplets: one dihedral per rotatable bond (j–k), with i chosen
 * among j's neighbors (≠ k) and l among k's neighbors (≠ j) by highest
 * neighbor degree (most substituted arm), smallest index on ties —
 * deterministic across runs.
 *
 * Typical results: benzene (6 ring bonds) → 0; HEPES/EPE (4W52, 15 bonds) →
 * 5 (piperazine tails + ethanesulfonate chain; ring + S(=O)₃ head excluded).
 */

/**
 * Normalize a bond list to [i, j] pairs (drops order element if present).
 * @param {Array} bonds  [[i,j]] or [[i,j,order], ...]
 * @returns {Array<Array<number>>} [[i,j], ...]
 */
function normBonds(bonds) {
  return (bonds ?? []).map((b) => [b[0], b[1]]);
}

/**
 * Build adjacency lists.
 * @param {number} n  atom count
 * @param {Array} bonds  [[i,j], ...]
 * @returns {Array<Array<number>>} adjacency
 */
function buildAdj(n, bonds) {
  const adj = Array.from({ length: n }, () => []);
  for (const [rawI, rawJ] of bonds) {
    const i = rawI, j = rawJ;
    if (i == null || j == null || i === j) continue;
    if (i < 0 || j < 0 || i >= n || j >= n) continue;
    if (!adj[i].includes(j)) adj[i].push(j);
    if (!adj[j].includes(i)) adj[j].push(i);
  }
  return adj;
}

/**
 * Test whether bond (a, b) lies on a cycle: BFS from a to b while ignoring
 * the direct edge a–b. Reachable ⇒ alternative path exists ⇒ in a ring.
 * @param {number} a
 * @param {number} b
 * @param {Array<Array<number>>} adj
 * @returns {boolean}
 */
function bondInRing(a, b, adj) {
  const seen = new Set([a]);
  const queue = [a];
  while (queue.length) {
    const v = queue.shift();
    for (const w of adj[v]) {
      if ((v === a && w === b) || (v === b && w === a)) continue;
      if (w === b) return true;
      if (!seen.has(w)) { seen.add(w); queue.push(w); }
    }
  }
  return false;
}

/**
 * Element symbol (uppercase) for an atom record.
 * @param {object} a  { element|elem }
 * @returns {string}
 */
function elemOf(a) {
  return String(a?.element ?? a?.elem ?? "C").trim().toUpperCase();
}

/**
 * Position [x,y,z] for an atom record.
 * @param {object} a
 * @returns {Array<number>}
 */
function posOf(a) {
  return [a?.x ?? a?.pos?.[0] ?? NaN, a?.y ?? a?.pos?.[1] ?? NaN, a?.z ?? a?.pos?.[2] ?? NaN];
}

/**
 * Find rotatable bonds in a ligand graph.
 * @param {Array} atoms  [{ x,y,z, element, gaffType? }]
 * @param {Array} bonds  [[i,j]] or [[i,j,order], ...]
 * @param {object} [opts]
 * @param {Array<number>} [opts.orders]  per-bond orders parallel to bonds
 * @param {Array<boolean>} [opts.aromatic]  per-bond aromatic flags
 * @returns {{ rotBonds:Array<Array<number>>, torsions:Array<Array<number>>, count:number }}
 *   rotBonds: [[j,k], ...] central bonds; torsions: [[i,j,k,l], ...] one per
 *   rotatable bond; count = rotBonds.length.
 */
export function findRotatableBonds(atoms, bonds, opts = {}) {
  const n = atoms?.length ?? 0;
  if (!n || !bonds?.length) return { rotBonds: [], torsions: [], count: 0 };
  const pairs = normBonds(bonds);
  const adj = buildAdj(n, pairs);
  const deg = adj.map((l) => l.length);
  const orders = opts.orders ?? null;
  const arom = opts.aromatic ?? null;

  // Per-bond order lookup: inline [i,j,order] wins, then opts.orders.
  const orderOf = (bIdx) => {
    const raw = bonds[bIdx];
    if (raw && raw.length > 2 && Number.isFinite(raw[2])) return raw[2];
    if (orders && Number.isFinite(orders[bIdx])) return orders[bIdx];
    return 1; // assumed single (documented; ring/terminal filters cover the risk)
  };
  const aromOf = (bIdx) => {
    const raw = bonds[bIdx];
    if (raw && raw.length > 3 && typeof raw[3] === "boolean") return raw[3];
    if (arom && typeof arom[bIdx] === "boolean") return arom[bIdx];
    if (Number.isFinite(orderOf(bIdx)) && Math.abs(orderOf(bIdx) - 1.5) < 0.26) return true;
    return false;
  };
  // Neighbor-bond index lookup for amide order check.
  const bondIdxBetween = (a, b) => {
    for (let bi = 0; bi < pairs.length; bi++) {
      if ((pairs[bi][0] === a && pairs[bi][1] === b) || (pairs[bi][0] === b && pairs[bi][1] === a)) return bi;
    }
    return -1;
  };

  /**
   * Amide C–N test for central bond (j, k).
   * @param {number} j
   * @param {number} k
   * @returns {boolean} true → exclude as amide
   */
  function isAmide(j, k) {
    const ej = elemOf(atoms[j]), ek = elemOf(atoms[k]);
    const isCN = (ej === "C" && ek === "N") || (ej === "N" && ek === "C");
    if (!isCN) return false;
    const cIdx = ej === "C" ? j : k;
    const nIdx = ej === "N" ? j : k;
    // (a) GAFF2-lite typing wins when present.
    const nt = String(atoms[nIdx]?.gaffType ?? "").toLowerCase();
    if (nt === "nh") return true;
    // (b) Explicit orders: C double-bonded to O.
    let hasOrderInfo = false;
    for (const w of adj[cIdx]) {
      if (w === nIdx) continue;
      if (elemOf(atoms[w]) !== "O") continue;
      const bi = bondIdxBetween(cIdx, w);
      if (bi >= 0 && (orders || (bonds[bi] && bonds[bi].length > 2))) {
        hasOrderInfo = true;
        if (orderOf(bi) >= 1.5) return true;
      }
    }
    if (hasOrderInfo) return false; // orders known, no carbonyl found → not amide
    // (c) Geometry fallback: C···O neighbor within carbonyl distance.
    const pc = posOf(atoms[cIdx]);
    for (const w of adj[cIdx]) {
      if (w === nIdx) continue;
      if (elemOf(atoms[w]) !== "O") continue;
      const po = posOf(atoms[w]);
      const r = Math.hypot(po[0] - pc[0], po[1] - pc[1], po[2] - pc[2]);
      if (Number.isFinite(r) && r < 1.35) return true;
    }
    return false;
  }

  const rotBonds = [];
  for (let b = 0; b < pairs.length; b++) {
    const [j, k] = pairs[b];
    if (j == null || k == null || j < 0 || k < 0 || j >= n || k >= n || j === k) continue;
    // Rule 1: single, non-aromatic.
    const o = orderOf(b);
    if (!(Math.abs(o - 1) < 0.15 + 1e-9)) continue;
    if (aromOf(b)) continue;
    // Rule 2: non-ring.
    if (bondInRing(j, k, adj)) continue;
    // Rule 3: non-amide.
    if (isAmide(j, k)) continue;
    // Rule 4: both sides substituted.
    if (deg[j] <= 1 || deg[k] <= 1) continue;
    // Rule 5: symmetric-top fan (≥3 leaf neighbors on either side).
    let exclude = false;
    for (const [side, other] of [[j, k], [k, j]]) {
      const rest = adj[side].filter((w) => w !== other);
      if (rest.length >= 3 && rest.every((w) => deg[w] <= 1)) { exclude = true; break; }
    }
    if (exclude) continue;
    rotBonds.push([j, k]);
  }

  // One torsion quadruplet per rotatable bond (deterministic arm choice).
  const torsions = rotBonds.map(([j, k]) => {
    let bestI = -1, bestDeg = -1;
    for (const w of adj[j]) {
      if (w === k) continue;
      const d = deg[w];
      if (d > bestDeg || (d === bestDeg && (bestI < 0 || w < bestI))) { bestDeg = d; bestI = w; }
    }
    let bestL = -1; bestDeg = -1;
    for (const w of adj[k]) {
      if (w === j) continue;
      const d = deg[w];
      if (d > bestDeg || (d === bestDeg && (bestL < 0 || w < bestL))) { bestDeg = d; bestL = w; }
    }
    return [bestI, j, k, bestL];
  }).filter((q) => q[0] >= 0 && q[3] >= 0);

  return { rotBonds, torsions, count: rotBonds.length };
}

/**
 * Emit absolute-index torsion quadruplets for a ligand placed at an offset
 * inside full-system frames (e.g. ligandStart = nProt).
 * @param {Array} atoms  ligand-local atoms
 * @param {Array} bonds  ligand-local bonds
 * @param {number} [offset=0]  atom-index offset added to each quadruplet entry
 * @param {object} [opts]  passed through to findRotatableBonds
 * @returns {Array<Array<number>>} absolute-index torsions (may be empty)
 */
export function autoTorsions(atoms, bonds, offset = 0, opts = {}) {
  const { torsions } = findRotatableBonds(atoms, bonds, opts);
  if (!offset) return torsions;
  return torsions.map((q) => q.map((i) => i + offset));
}
