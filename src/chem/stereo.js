/**
 * stereo.js — Chirality + planarity checker (Phase 2: Chemistry).
 *
 *   checkStereochemistry(mol, opts) — tetrahedral stereocenter detection,
 *     chiral-volume computation, CIP-lite R/S tags, inversion flags vs. a
 *     reference volume set.
 *   validatePlanarity(mol, opts)    — aromatic / sp2 / caller-supplied group
 *     plane-deviation check (Newell plane fit).
 *
 * Conventions:
 *   - Heavy atoms only (engine drops H at parse): a tetrahedral candidate is
 *     an sp3 C/N (Si/P tolerated) with 3–4 heavy neighbors joined by single
 *     (order ≈ 1, non-aromatic) bonds. Degree-3 centers carry one implicit H,
 *     ranked lowest in the CIP-lite ordering.
 *   - Chiral volume V = (b−a)·((c−a)×(d−a)) over the four neighbor positions
 *     in deterministic (sorted-index) order; |V| < 0.05 Å³ reads as planar.
 *   - R/S tags use a CIP-lite ranking (atomic number, then substituent size,
 *     then distance — heavy atoms only). Tags are self-consistent parity
 *     labels for inversion detection, NOT certified CIP assignments (full CIP
 *     needs explicit H + sequence-rule recursion; Documented limitation).
 *   - Pyramidal degree-3 N is reported with stable:false (rapid inversion in
 *     solution) but still volume-tracked.
 *
 * Pure functions, no imports, never throws on malformed input (reports via
 * `warnings`). Input molecule shape:
 *   { atoms: [{ x,y,z | pos, element|elem }], bonds: [[i,j], ...] }
 * Bond orders accepted as [[i,j]] (assumed single unless aromatic flags given)
 * or [[i,j,order]] / { orders, aromaticBonds } via opts.
 */

const VOLUME_PLANAR = 0.05; // Å³ — |V| below this is not a stereocenter

function posOf(a) {
  return [a.x ?? a.pos?.[0] ?? NaN, a.y ?? a.pos?.[1] ?? NaN, a.z ?? a.pos?.[2] ?? NaN];
}

function elemOf(a) { return String(a.element ?? a.elem ?? "C").trim().toUpperCase(); }

function triple(ax, ay, az, bx, by, bz, cx, cy, cz) {
  return ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx);
}

/**
 * Signed chiral volume of four points (6× tetrahedron volume).
 * @param {Array<number>} a  [x,y,z] apex (usually the stereocenter)
 * @param {Array<number>} b
 * @param {Array<number>} c
 * @param {Array<number>} d
 * @returns {number}  Å³, sign = handedness for the given neighbor order
 */
export function chiralVolume(a, b, c, d) {
  const bax = b[0] - a[0], bay = b[1] - a[1], baz = b[2] - a[2];
  const cax = c[0] - a[0], cay = c[1] - a[1], caz = c[2] - a[2];
  const dax = d[0] - a[0], day = d[1] - a[1], daz = d[2] - a[2];
  // (b−a)·((c−a)×(d−a))
  const cx = cay * daz - caz * day;
  const cy = caz * dax - cax * daz;
  const cz = cax * day - cay * dax;
  return bax * cx + bay * cy + baz * cz;
}

function buildAdj(atoms, bonds, orders, aromaticBonds) {
  const adj = Array.from({ length: atoms.length }, () => []);
  bonds.forEach((bd, b) => {
    const i = bd[0], j = bd[1];
    if (i == null || j == null || i === j || i < 0 || j < 0 || i >= atoms.length || j >= atoms.length) return;
    const order = bd.length > 2 ? bd[2] : (orders?.[b] ?? 1);
    const arom = aromaticBonds?.[b] ?? (order === 1.5);
    adj[i].push({ nbr: j, order, arom });
    adj[j].push({ nbr: i, order, arom });
  });
  return adj;
}

const ZNUM = {
  H: 1, B: 5, C: 6, N: 7, O: 8, F: 9, NA: 11, MG: 12, AL: 13, SI: 14, P: 15,
  S: 16, CL: 17, K: 19, CA: 20, MN: 25, FE: 26, CO: 27, NI: 28, CU: 29, ZN: 30,
  SE: 34, BR: 35, I: 53,
};

/**
 * CIP-lite substituent rank key: [atomic number, heavy-degree, H-corrected mass proxy].
 * Implicit H sorts lowest by construction (Z=1).
 */
function rankKey(atoms, adj, idx) {
  if (idx === -1) return [1, 0, 1]; // implicit H
  const el = elemOf(atoms[idx]);
  return [ZNUM[el] ?? 6, adj[idx].length, (ZNUM[el] ?? 6) * 4 + adj[idx].length];
}

function cmpRank(a, b) {
  for (let k = 0; k < 3; k++) {
    if (a[k] !== b[k]) return b[k] - a[k]; // descending: rank 1 = highest
  }
  return 0;
}

/**
 * Detect tetrahedral stereocenters and compute their chiral volumes.
 *
 * @param {object} mol  { atoms, bonds }
 * @param {object} [opts]
 * @param {Array<number>} [opts.orders]  per-bond orders parallel to bonds
 * @param {Array<boolean>} [opts.aromaticBonds]
 * @param {object} [opts.ref]  { index -> volume } reference volumes for inversion flags
 * @param {number} [opts.planarTol=0.05]  |V| below this → planar/achiral
 * @returns {{ centers:Array, nCenters:number, nInverted:number, warnings:string[] }}
 *   centers: [{ index, element, neighbors:[i,j,k,l(-1=implicit H)], volume,
 *     tag:'R'|'S'|'planar'|'degenerate', stable:boolean, inverted:boolean, note }]
 */
export function checkStereochemistry(mol, opts = {}) {
  const warnings = [];
  const atoms = mol?.atoms ?? [];
  const bonds = mol?.bonds ?? [];
  const planarTol = opts.planarTol ?? VOLUME_PLANAR;
  const adj = buildAdj(atoms, bonds, opts.orders, opts.aromaticBonds);
  const centers = [];
  const TETRA = new Set(["C", "N", "SI", "P"]);

  for (let i = 0; i < atoms.length; i++) {
    const el = elemOf(atoms[i]);
    if (!TETRA.has(el)) continue;
    const links = adj[i].filter((e) => !e.arom && e.order <= 1.15);
    if (links.length !== adj[i].length) continue; // pi-bonded — not tetrahedral
    if (links.length < 3 || links.length > 4) continue;
    const nbrs = links.map((e) => e.nbr).sort((a, b) => a - b);
    const withH = links.length === 3 ? [...nbrs, -1] : [...nbrs];
    const P = atoms.map(posOf);
    if (P.some((p) => !Number.isFinite(p[0] + p[1] + p[2]))) {
      warnings.push(`center ${i}: non-finite coordinates — skipped`);
      continue;
    }
    // Signed tetrahedral volume det[p1−p0, p2−p0, p3−p0] over the four
    // substituent positions (implicit H placed geometrically when needed).
    const full = det4(P, i, withH);
    let entry;
    if (Math.abs(full) < planarTol) {
      entry = {
        index: i, element: el, neighbors: withH, volume: full,
        tag: "planar", stable: false, inverted: false,
        note: `|V|=${Math.abs(full).toFixed(3)}Å³ < ${planarTol}Å³ — planar/achiral arrangement`,
      };
    } else {
      const tag = cipLiteTag(P, atoms, adj, i, withH);
      const allEqual = withH.every((w) => cmpRank(rankKey(atoms, adj, w), rankKey(atoms, adj, withH[0])) === 0);
      entry = {
        index: i, element: el, neighbors: withH, volume: full,
        tag: allEqual ? "degenerate" : tag,
        stable: !(el === "N" && links.length === 3),
        inverted: false,
        note: allEqual ? "four indistinguishable substituents — achiral"
          : (el === "N" && links.length === 3) ? `pyramidal N (V=${full.toFixed(2)}Å³), fluxional — volume-tracked only`
          : `tetrahedral stereocenter (V=${full.toFixed(2)}Å³, CIP-lite ${allEqual ? "degenerate" : tag})`,
      };
    }
    if (opts.ref && Number.isFinite(opts.ref[i])) {
      const rv = opts.ref[i];
      if (Math.abs(rv) >= planarTol && Math.abs(entry.volume) >= planarTol) {
        entry.inverted = Math.sign(rv) !== Math.sign(entry.volume);
        if (entry.inverted) entry.note += ` — INVERTED vs reference (${rv.toFixed(2)}→${entry.volume.toFixed(2)}Å³)`;
      }
    }
    centers.push(entry);
  }
  return { centers, nCenters: centers.filter((c) => c.tag === "R" || c.tag === "S").length, nInverted: centers.filter((c) => c.inverted).length, warnings };
}

/** Signed volume of tetrahedron (p0,p1,p2,p3): det[p1−p0, p2−p0, p3−p0]. */
function det4(P, i, withH) {
  const pts = withH.map((w) => (w === -1 ? implicitHPos(P, i, withH) : P[w]));
  const [p0, p1, p2, p3] = pts;
  return triple(
    p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2],
    p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2],
    p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2],
  );
}

/** Place implicit H opposite the centroid of the three heavy neighbors. */
function implicitHPos(P, i, withH) {
  const c = P[i];
  const heavy = withH.filter((w) => w !== -1).map((w) => P[w]);
  const cx = heavy.reduce((s, p) => s + p[0], 0) / heavy.length;
  const cy = heavy.reduce((s, p) => s + p[1], 0) / heavy.length;
  const cz = heavy.reduce((s, p) => s + p[2], 0) / heavy.length;
  const dx = c[0] - cx, dy = c[1] - cy, dz = c[2] - cz;
  const L = Math.hypot(dx, dy, dz) || 1;
  return [c[0] + (dx / L) * 1.09, c[1] + (dy / L) * 1.09, c[2] + (dz / L) * 1.09];
}

/**
 * CIP-lite R/S tag: rank substituents, view from lowest priority, read 1→2→3.
 * det > 0 → 'R', < 0 → 'S' (self-consistent parity convention).
 */
function cipLiteTag(P, atoms, adj, i, withH) {
  const ranked = withH.map((w) => ({ w, r: rankKey(atoms, adj, w) }))
    .sort((a, b) => cmpRank(a.r, b.r) || (a.w - b.w));
  const [p1, p2, p3, p4] = ranked.map(({ w }) => (w === -1 ? implicitHPos(P, i, withH) : P[w]));
  const det = triple(
    p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2],
    p2[0] - p4[0], p2[1] - p4[1], p2[2] - p4[2],
    p3[0] - p4[0], p3[1] - p4[1], p3[2] - p4[2],
  );
  return det >= 0 ? "R" : "S";
}

/**
 * Validate planarity of aromatic / sp2 / caller-supplied atom groups.
 * Groups default to: connected `ca`-typed sets (when atoms carry gaffType)
 * plus auto-detected 5/6-membered C/N/O/S rings with short bonds.
 *
 * @param {object} mol  { atoms, bonds }
 * @param {object} [opts]
 * @param {Array<Array<number>>} [opts.groups]  explicit index groups
 * @param {number} [opts.rmsTol=0.10]   RMS tolerance (Å)
 * @param {number} [opts.maxTol=0.25]   max-deviation tolerance (Å)
 * @returns {{ groups:Array<{indices:number[],rms:number,maxDev:number,planar:boolean}>, ok:boolean, warnings:string[] }}
 */
export function validatePlanarity(mol, opts = {}) {
  const warnings = [];
  const atoms = mol?.atoms ?? [];
  const bonds = mol?.bonds ?? [];
  const rmsTol = opts.rmsTol ?? 0.10;
  const maxTol = opts.maxTol ?? 0.25;
  let groups = opts.groups ?? null;
  if (!groups) {
    groups = autoPlanarGroups(atoms, bonds);
    if (!groups.length) warnings.push("no planar groups detected (no ca types / aromatic-length rings)");
  }
  const out = groups.map((indices) => {
    const { rms, maxDev } = planeDeviation(atoms, indices);
    return { indices: [...indices], rms, maxDev, planar: rms <= rmsTol && maxDev <= maxTol };
  });
  return { groups: out, ok: out.every((g) => g.planar), warnings };
}

function autoPlanarGroups(atoms, bonds) {
  const groups = [];
  // (a) connected ca-typed components
  const caSet = new Set(atoms.map((a, i) => (String(a.gaffType ?? "").toLowerCase() === "ca" ? i : -1)).filter((i) => i >= 0));
  if (caSet.size >= 4) {
    const adj = new Map();
    for (const [i, j] of bonds) {
      if (!caSet.has(i) || !caSet.has(j)) continue;
      if (!adj.has(i)) adj.set(i, []);
      if (!adj.has(j)) adj.set(j, []);
      adj.get(i).push(j); adj.get(j).push(i);
    }
    const seen = new Set();
    for (const s of caSet) {
      if (seen.has(s)) continue;
      const comp = [];
      const stack = [s];
      seen.add(s);
      while (stack.length) {
        const v = stack.pop();
        comp.push(v);
        for (const w of (adj.get(v) ?? [])) {
          if (!seen.has(w)) { seen.add(w); stack.push(w); }
        }
      }
      if (comp.length >= 4) groups.push(comp.sort((x, y) => x - y));
    }
  }
  // (b) short-bond 5/6-rings (fallback when untyped)
  if (!groups.length) {
    const adj = Array.from({ length: atoms.length }, () => []);
    for (const [i, j] of bonds) {
      if (i < 0 || j < 0 || i >= atoms.length || j >= atoms.length || i === j) continue;
      adj[i].push(j); adj[j].push(i);
    }
    const P = atoms.map(posOf);
    const blen = (a, b) => Math.hypot(P[a][0] - P[b][0], P[a][1] - P[b][1], P[a][2] - P[b][2]);
    const seen = new Set();
    for (let start = 0; start < atoms.length; start++) {
      const path = [start];
      const inPath = new Set([start]);
      const dfs = (node, parent) => {
        for (const nbr of adj[node]) {
          if (nbr === parent) continue;
          if (nbr === start) {
            if (path.length === 5 || path.length === 6) {
              const key = [...path].sort((x, y) => x - y).join("-");
              if (seen.has(key)) continue;
              seen.add(key);
              let ok = true;
              for (let t = 0; t < path.length && ok; t++) {
                const r = blen(path[t], path[(t + 1) % path.length]);
                if (r < 1.28 || r > 1.52) ok = false;
              }
              if (ok) groups.push([...path]);
            }
          } else if (!inPath.has(nbr) && path.length < 6) {
            inPath.add(nbr); path.push(nbr); dfs(nbr, node); path.pop(); inPath.delete(nbr);
          }
        }
      };
      dfs(start, -1);
    }
  }
  return groups;
}

/**
 * RMS + max out-of-plane deviation of a point set (Newell plane).
 * @param {Array} atoms
 * @param {Array<number>} indices
 * @returns {{ rms:number, maxDev:number }}  Å
 */
export function planeDeviation(atoms, indices) {
  const P = indices.map((i) => posOf(atoms[i]));
  const n = P.length;
  if (n < 3) return { rms: 0, maxDev: 0 };
  let nx = 0, ny = 0, nz = 0;
  for (let i = 0; i < n; i++) {
    const [x0, y0, z0] = P[i];
    const [x1, y1, z1] = P[(i + 1) % n];
    if (!Number.isFinite(x0 + y0 + z0 + x1 + y1 + z1)) return { rms: NaN, maxDev: NaN };
    nx += (y0 - y1) * (z0 + z1);
    ny += (z0 - z1) * (x0 + x1);
    nz += (x0 - x1) * (y0 + y1);
  }
  const L = Math.hypot(nx, ny, nz);
  if (L < 1e-12) return { rms: 0, maxDev: 0 };
  nx /= L; ny /= L; nz /= L;
  const cx = P.reduce((s, p) => s + p[0], 0) / n;
  const cy = P.reduce((s, p) => s + p[1], 0) / n;
  const cz = P.reduce((s, p) => s + p[2], 0) / n;
  let s2 = 0, mx = 0;
  for (const [x, y, z] of P) {
    const d = Math.abs((x - cx) * nx + (y - cy) * ny + (z - cz) * nz);
    s2 += d * d;
    if (d > mx) mx = d;
  }
  return { rms: Math.sqrt(s2 / n), maxDev: mx };
}
