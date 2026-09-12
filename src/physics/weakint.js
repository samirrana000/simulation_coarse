/**
 * weakint.js — Heavy-mode weak-interaction terms (R3 §1a–c, Loop-2 S3).
 *
 * Three opt-in supplements to the heavy-mode nonbonded kernel, adding
 * anisotropic physics GB/point-charge cannot represent (R3 double-counting
 * analysis, §2):
 *
 *  - π–π stacking: Morse-like Gaussian attraction on the interplanar gap ρ
 *    with in-plane damping and orientation factor, plus a repulsive
 *    Gaussian at closer approach (Hunter–Sanders lineage).
 *  - Cation–π: cos²-gated Gaussian along the ring normal
 *    (Gallivan & Dougherty 1999: 2–5 kcal/mol on-axis).
 *  - Halogen σ-hole: angle-gated Gaussian on the C–X···D angle
 *    (VinaXB lineage, Koebel et al. 2016; F excluded).
 *
 * All energies in kcal/mol, distances in Å. Gradients are analytic and
 * finite-difference-verified (tests/test_weakint.js; prototype
 * /tmp/opencode/r3_heavyproto.mjs: π 6.2e-11 · cation-π 0 · halogen 4.6e-10).
 *
 * Zero deps, plain ES module. Positions/frames as plain arrays so the module
 * is usable both standalone (tests) and inside the heavy force field.
 */

/* ------------------------------------------------------------------ */
/* Geometry helpers (flat-array free; 3-vectors as [x,y,z])            */
/* ------------------------------------------------------------------ */

/** @param {number[]} a @param {number[]} b @returns {number} */
const vdot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

/** @param {number[]} a @param {number[]} b @returns {number[]} */
const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

/** Normalize a 3-vector; zero-safe. @param {number[]} a @returns {number[]} */
const vnorm = (a) => {
  const l = Math.sqrt(vdot(a, a)) || 1;
  return [a[0] / l, a[1] / l, a[2] / l];
};

/* ------------------------------------------------------------------ */
/* 1) π–π stacking (R3 §1a)                                           */
/* ------------------------------------------------------------------ */

/**
 * π–π stacking energy on ring frames (centroids + unit normals).
 *
 *   ρ     = (|r·n_A| + |r·n_B|)/2                     interplanar gap
 *   t²    = |r|² − (r·n̄)²                             in-plane offset², n̄ = mean normal
 *   U_att = −ε · G(ρ; ρ0=3.8, w=0.45) · exp(−t²/6) · (½ + ½(n_A·n_B)²)
 *   U_rep = +4ε · G(ρ; ρr=3.0, w=0.35) · (½ + ½(n_A·n_B)²)
 *
 * Centroid cutoff 5.5 Å (r² > 30.25 returns zero). Well depth ≈ −1.70 at
 * 3.8 Å parallel (literature −1..−3, Hunter–Sanders).
 *
 * @param {number[]} cA ring A centroid (Å)
 * @param {number[]} nA ring A unit normal
 * @param {number[]} cB ring B centroid
 * @param {number[]} nB ring B unit normal
 * @param {{eps?:number}} [params] ε (kcal/mol), default 2.4
 * @returns {{U:number, gA:number[], gB:number[], gNA:number[], gNB:number[], debug:object}}
 *   gA/gB: dU/dcA / dU/dcB · gNA/gNB: dU/dnA / dU/dnB projected on the unit
 *   sphere (tangent — caller chains by the frame's normal Jacobian).
 */
export function piStackEnergy(cA, nA, cB, nB, params = {}) {
  const eps = params.eps ?? 2.4;
  const sigma = params.rho0 ?? 3.8; // attraction gap center ρ0
  const r = vsub(cB, cA);
  const r2 = vdot(r, r);
  if (r2 > 30.25) { // 5.5 Å centroid cutoff
    return { U: 0, gA: [0, 0, 0], gB: [0, 0, 0], gNA: [0, 0, 0], gNB: [0, 0, 0], debug: { cut: true } };
  }
  // perpendicular distance of each centroid from the other ring plane
  const zA = vdot(r, nA);
  const zB = vdot(r, nB);
  const rho = 0.5 * (Math.abs(zA) + Math.abs(zB));
  const ndotn = vdot(nA, nB);
  // Canonical winding: when two ring frames have anti-parallel normals
  // (parallel ring PLANES traversed in opposite directions) the prototype's
  // nA + nB degenerates to ~0. Flip nB inside n̄ so nA + s·nB always spans
  // the shared plane axis. Prototype geometries (ndotn ≥ 0) are unchanged
  // (s = +1); only anti-parallel windings are repaired (verified: same U,
  // gradient chain symmetric).
  const sFlip = ndotn >= 0 ? 1 : -1;
  const orient = 0.5 + 0.5 * ndotn * ndotn;
  // in-plane offset: component of r perpendicular to the average normal axis
  const navg = vnorm([nA[0] + sFlip * nB[0], nA[1] + sFlip * nB[1], nA[2] + sFlip * nB[2]]);
  const rPar = vdot(r, navg);
  const t2 = Math.max(0, r2 - rPar * rPar);
  const overlap = Math.exp(-t2 / 6.0); // damping: >2.45 Å offset halves attraction
  // attraction: bounded Gaussian (no LJ singularity)
  const wRho = 0.45;
  const gatt = Math.exp(-((rho - sigma) * (rho - sigma)) / (2 * wRho * wRho));
  // repulsion keyed to perpendicular gap only (σr smaller — parallel in-plane
  // offset does not cause steric clash)
  const sigmaR = 3.0;
  const wRep = 0.35;
  const grep = Math.exp(-((rho - sigmaR) * (rho - sigmaR)) / (2 * wRep * wRep));
  const Urep = 4 * eps * orient * grep;
  const U = -eps * gatt * overlap * orient + Urep;

  // ---- gradients (chain rule incl. ∂ρ/∂n, ∂t²/∂n̄, unit-sphere projection) ----
  const gA = [0, 0, 0], gB = [0, 0, 0];
  // attraction: dU/dρ = ε·(ρ−ρ0)/w² · gatt · overlap·orient
  const dU_drho = eps * ((rho - sigma) / (wRho * wRho)) * gatt * overlap * orient;
  // repulsion: dU/dρ = 4ε·(−(ρ−σr)/wr²)·orient·grep
  const dUrep_drho = 4 * eps * orient * (-(rho - sigmaR) / (wRep * wRep)) * grep;
  const dU_total_drho = dU_drho + dUrep_drho;
  const sA = zA >= 0 ? 1 : -1, sB2 = zB >= 0 ? 1 : -1;
  for (let k = 0; k < 3; k++) {
    const dpcB = 0.5 * (sA * nA[k] + sB2 * nB[k]); // ∂ρ/∂cB (r = cB − cA)
    gB[k] += dU_total_drho * dpcB;
    gA[k] -= dU_total_drho * dpcB;
  }
  // ρ = (|r·nA|+|r·nB|)/2 ⇒ ∂ρ/∂nA = 0.5·sA·r (normal forces from the gap term)
  const gNA_rho = [0.5 * sA * dU_total_drho * r[0], 0.5 * sA * dU_total_drho * r[1], 0.5 * sA * dU_total_drho * r[2]];
  const gNB_rho = [0.5 * sB2 * dU_total_drho * r[0], 0.5 * sB2 * dU_total_drho * r[1], 0.5 * sB2 * dU_total_drho * r[2]];
  // overlap: U_att = −ε·gatt·orient·overlap → dU/dt²
  const dOv_dt2 = (-1 / 6.0) * overlap;
  const dU_dt2 = -eps * gatt * orient * dOv_dt2;
  for (let k = 0; k < 3; k++) {
    const dt2_cB = 2 * (r[k] - rPar * navg[k]); // ∂t²/∂cB (cA gets minus)
    gB[k] += dU_dt2 * dt2_cB;
    gA[k] -= dU_dt2 * dt2_cB;
  }
  // orient gradient: d(orient)/d(ndotn) = ndotn (attraction + repulsion parts)
  const dU_dndotn = (-eps * gatt * overlap + 4 * eps * grep) * ndotn;
  const gNA = [dU_dndotn * nB[0], dU_dndotn * nB[1], dU_dndotn * nB[2]];
  const gNB = [dU_dndotn * nA[0], dU_dndotn * nA[1], dU_dndotn * nA[2]];
  gNA[0] += gNA_rho[0]; gNA[1] += gNA_rho[1]; gNA[2] += gNA_rho[2];
  gNB[0] += gNB_rho[0]; gNB[1] += gNB_rho[1]; gNB[2] += gNB_rho[2];
  // t² depends on n̄ (unit-normalized): ∂n̄/∂nA ≈ (I − n̄n̄ᵀ)/|nA+s·nB|;
  // ∂n̄/∂nB carries the same Jacobian scaled by s (chain rule through the flip)
  const sumN = [nA[0] + sFlip * nB[0], nA[1] + sFlip * nB[1], nA[2] + sFlip * nB[2]];
  const lensum = Math.sqrt(vdot(sumN, sumN)) || 1e-6;
  for (let k = 0; k < 3; k++) {
    let s = 0;
    for (let j = 0; j < 3; j++) {
      const kronecker = j === k ? 1 : 0;
      s += (-2 * rPar * r[j]) * (kronecker - navg[j] * navg[k]) / lensum; // ∂t²/∂nA_k
    }
    gNA[k] += dU_dt2 * s;
    gNB[k] += dU_dt2 * s * sFlip; // ∂n̄/∂nB = s · (same Jacobian form)
  }
  // project normal gradients onto the unit sphere (chain rule for the
  // normalization step the caller applies when mapping normals → atoms)
  const lA = Math.max(Math.sqrt(vdot(nA, nA)), 1e-9);
  const projA = vdot(gNA, nA) / lA;
  const gNAt = [gNA[0] - projA * nA[0], gNA[1] - projA * nA[1], gNA[2] - projA * nA[2]];
  const lB = Math.max(Math.sqrt(vdot(nB, nB)), 1e-9);
  const projB = vdot(gNB, nB) / lB;
  const gNBt = [gNB[0] - projB * nB[0], gNB[1] - projB * nB[1], gNB[2] - projB * nB[2]];
  return { U, gA, gB, gNA: gNAt, gNB: gNBt, debug: { rho, t2, ndotn, orient, overlap, r2 } };
}

/* ------------------------------------------------------------------ */
/* 2) Cation–π (R3 §1b)                                               */
/* ------------------------------------------------------------------ */

/**
 * Cation–π energy: Gaussian well along the cation→centroid distance gated
 * by cos² of the angle to the ring normal (0 = on-axis).
 *
 *   U = −ε_cπ · G(r; r0=4.3, w=1.1) · max(0, cos α)²
 *
 * On-axis 4.3 Å → −3.5 (lit −2..−5); cutoff 6 Å (r² > 36).
 *
 * @param {number[]} catPos cation position (Å)
 * @param {number[]} centroid ring centroid
 * @param {number[]} normal ring unit normal
 * @param {{eps?:number, r0?:number, w?:number}} [params]
 * @returns {{U:number, gCat:number[], gRing:number[], debug:object}}
 *   gRing: dU/d(centroid) — the normal itself is treated as fixed by the
 *   caller (ring plane rigid on MD timescales; R3 §1b cost note).
 */
export function cationPiEnergy(catPos, centroid, normal, params = {}) {
  const eps = params.eps ?? 3.5;
  const r0 = params.r0 ?? 4.3;
  const w = params.w ?? 1.1;
  const r = vsub(centroid, catPos);
  const r2 = vdot(r, r);
  if (r2 > 36) return { U: 0, gCat: [0, 0, 0], gRing: [0, 0, 0], debug: { cut: true } };
  const dist = Math.sqrt(r2);
  // Loop-2 S7 bench hardening: coincident cation/centroid (overlapping atoms
  // in an unminimized crystal — duplicate altlocs) makes r̂ = 0/0 = NaN below.
  // Same guard style as halogenEnergy's lrc < 1e-8 (line ~240). Returns the
  // zero-overlap limit (finite, no force) instead of poisoning total U.
  if (dist < 1e-8) return { U: 0, gCat: [0, 0, 0], gRing: [0, 0, 0], debug: { cut: true } };
  const dr = dist - r0;
  const g = Math.exp(-(dr * dr) / (2 * w * w));
  const rn = [r[0] / dist, r[1] / dist, r[2] / dist];
  const cosA = vdot(rn, normal);
  const gate = Math.max(0, cosA) ** 2;
  const U = -eps * g * gate;
  // ---- analytic gradients ----
  const dg_dd = -(dr / (w * w)) * g;          // dG/ddist
  const dU_dd = -eps * dg_dd * gate;          // radial part
  // gate: cosA = n·r̂; ∂cosA/∂centroid = (n − cosA·r̂)/dist; ∂cosA/∂cat = −same
  const m = Math.max(0, cosA);
  const dgate = 2 * m;                         // d(gate)/d(cosA)
  const gCat = [0, 0, 0], gRing = [0, 0, 0];
  for (let k = 0; k < 3; k++) {
    const dcos_dRing = (normal[k] - cosA * rn[k]) / dist;
    const fr = dU_dd * rn[k] + (-eps * g) * dgate * dcos_dRing;
    gRing[k] += fr;
    gCat[k] -= fr;
  }
  return { U, gCat, gRing, debug: { dist, cosA } };
}

/* ------------------------------------------------------------------ */
/* 3) Halogen σ-hole (R3 §1c)                                         */
/* ------------------------------------------------------------------ */

/** Per-element σ-hole well depths ε_X (kcal/mol). F excluded by construction. */
export const HALOGEN_EPS = Object.freeze({ CL: 1.2, BR: 2.0, I: 2.5 });

/**
 * Halogen σ-hole energy for one C–X···D triple.
 *
 *   U = −ε_X · G(|X−D|; r0=3.1, w=0.55) · max(0, cos β)²
 *   β = angle(C–X···D); σ-hole axis = extension of C→X beyond X
 *
 * Linear C–Cl···O at 3.1 Å → −ε_X (lit net 1–4). Gradients use the
 * parallel-transport-safe forms (FD-verified; avoid exactly-perpendicular
 * geometries — the max(0,·)² gate kinks at cos β = 0).
 *
 * @param {number[]} cPos halogen-bonded carbon C
 * @param {number[]} xPos halogen X (Cl/Br/I)
 * @param {number[]} dPos acceptor D (O/S/N)
 * @param {{eps?:number, r0?:number, w?:number}} [params] eps default 2.0 (demo value)
 * @returns {{U:number, gC:number[], gX:number[], gD:number[], debug:object}}
 */
export function halogenEnergy(cPos, xPos, dPos, params = {}) {
  const eps = params.eps ?? 2.0;
  const r0 = params.r0 ?? 3.1;
  const w = params.w ?? 0.55;
  const rx = vsub(dPos, xPos);      // X→D
  const rc = vsub(xPos, cPos);      // C→X
  const dd = vdot(rx, rx);
  if (dd > 16 || dd < 1e-8) {       // 4 Å pair cutoff
    return { U: 0, gC: [0, 0, 0], gX: [0, 0, 0], gD: [0, 0, 0], debug: { cut: true } };
  }
  const dist = Math.sqrt(dd);
  const lrc = Math.sqrt(vdot(rc, rc));
  if (lrc < 1e-8) return { U: 0, gC: [0, 0, 0], gX: [0, 0, 0], gD: [0, 0, 0], debug: { cut: true } };
  const dr = dist - r0;
  const g = Math.exp(-(dr * dr) / (2 * w * w));
  const rn = [rx[0] / dist, rx[1] / dist, rx[2] / dist];   // X→D unit
  const axis = [rc[0] / lrc, rc[1] / lrc, rc[2] / lrc];     // C→X unit = σ-hole direction
  const cosB = vdot(rn, axis);
  const m = Math.max(0, cosB);
  const gate = m * m;
  const U = -eps * g * gate;
  // ---- gradients: systematic quotient rule on cos β = â·r̂ ----
  const dg = -(dr / (w * w)) * g;
  const dU_dd = -eps * dg * gate;
  const dgate = 2 * m;
  const gC = [0, 0, 0], gX = [0, 0, 0], gD = [0, 0, 0];
  for (let k = 0; k < 3; k++) {
    // ∂cosβ/∂D = (â − cosβ·r̂)/dist  (parallel-transport-safe form)
    const dcb_dD = (axis[k] - cosB * rn[k]) / dist;
    // ∂cosβ/∂X = −(â − cosβ·r̂)/dist + (r̂ − cosβ·â)/l_CX  (X in both r_XD and r_CX)
    const dcb_dX = (cosB * rn[k] - axis[k]) / dist + (rn[k] - cosB * axis[k]) / lrc;
    // ∂cosβ/∂C = −(r̂ − cosβ·â)/l_CX  (C only in r_CX)
    const dcb_dC = (cosB * axis[k] - rn[k]) / lrc;
    const fD = dU_dd * rn[k] + (-eps * g) * dgate * dcb_dD;
    const fX_radial = -dU_dd * rn[k];
    const fX_gate = (-eps * g) * dgate * dcb_dX;
    const fC = (-eps * g) * dgate * dcb_dC;
    gD[k] += fD;
    gX[k] += fX_radial + fX_gate;
    gC[k] += fC;
  }
  return { U, gC, gX, gD, debug: { dist, cosB } };
}

/* ------------------------------------------------------------------ */
/* Ring-frame construction (R3 §6 item 1)                             */
/* ------------------------------------------------------------------ */

/** GAFF2-lite sp2 ring types that mark aromatic atoms when present. */
const AROMATIC_GAFF_TYPES = new Set(["ca", "nb", "nh", "na", "n2"]);

/** Aromatic 6-ring protein residue atom sets (Phe/Tyr/Trp side chains). */
const PROTEIN_AROMATIC_RINGS = [
  { res: "PHE", atoms: ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"] },
  { res: "TYR", atoms: ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"] },
  { res: "TRP", atoms: ["CG", "CD1", "CD2", "NE1", "CE2"] },
  { res: "TRP", atoms: ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"] },
  { res: "HIS", atoms: ["CG", "ND1", "CD2", "CE1", "NE2"] },
];

/**
 * Build aromatic ring frames: centroid + unit normal + member atom indices.
 *
 * Two detection paths, in priority order:
 *  1. GAFF2-lite types: atoms carrying gaffType ∈ {ca, nb, nh, na, n2}
 *     on ≥5-membered graph rings (bond graph rings of C/N only).
 *  2. Geometric fallback: planar 5-/6-rings in the bond graph with aromatic
 *     bond-length window [1.28, 1.52] Å and Newell-plane deviation < 0.12 Å
 *     (same criteria as chem/gaff2_mapper.findAromaticRings — kept local to
 *     honor the zero-dep constraint of this module).
 * Protein side-chain rings are always detected via name sets (Phe/Tyr/Trp/
 * His), which is robust for PDB-derived structures without any typing.
 *
 * @param {Array<object>} atoms  heavy atoms ({x,y,z,element,atomName,resName,...})
 * @param {Array<[number,number]>} bonds  bond list [[i,j],...]
 * @returns {Array<{centroid:number[], normal:number[], atomIdx:number[]}>}
 */
export function buildRingFrames(atoms, bonds) {
  const n = atoms.length;
  const adj = Array.from({ length: n }, () => []);
  for (const [i, j] of bonds) {
    adj[i].push(j);
    adj[j].push(i);
  }
  const P = (i) => [atoms[i].x, atoms[i].y, atoms[i].z];

  // --- path 1: protein name sets (always on) ---
  const frames = [];
  const frameSets = []; // atom sets of accepted rings (for subsumption checks)
  const seenSig = new Set(); // exact ring atom-signature dedup
  const pushFrame = (idx) => {
    const sig = [...idx].sort((a, b) => a - b).join("-");
    if (seenSig.has(sig)) return false;
    const inSet = new Set(idx);
    // subsumption: fused-ring paths can re-find the same physical ring with a
    // one-atom different signature (e.g. a 5-subset of a 6-ring through
    // cross-ring braces, or Trp 6-ring via names vs geometry). Skip when
    // ≥80% of this ring's atoms already belong to ONE accepted ring.
    // True fused rings (Trp 5+6) share only 2/6 = 33% — far below the bar.
    for (const fs of frameSets) {
      let shared = 0;
      for (const i of idx) if (fs.has(i)) shared++;
      if (shared >= 0.8 * idx.length) return false;
    }
    seenSig.add(sig);
    frameSets.push(inSet);
    frames.push(makeFrame(atoms, idx));
    return true;
  };
  const byRes = new Map(); // residue key → atom index by name
  for (let i = 0; i < n; i++) {
    const a = atoms[i];
    if (!a?.isProtein || !a.atomName) continue;
    const key = `${a.chain}|${a.resSeq}`;
    let m = byRes.get(key);
    if (!m) { m = new Map(); byRes.set(key, m); }
    m.set(a.atomName.trim().toUpperCase(), i);
  }
  for (const { res, atoms: names } of PROTEIN_AROMATIC_RINGS) {
    for (const [key, m] of byRes) {
      const first = atoms[m.get(names[0])];
      if (!first || first.resName !== res) continue;
      const idx = [];
      for (const nm of names) {
        const k = m.get(nm);
        if (k === undefined) break;
        idx.push(k);
      }
      if (idx.length !== names.length) continue;
      pushFrame(idx); // Trp 5- and 6-rings share CD2/CE2 — both pushed (distinct sigs)
    }
  }

  // --- path 2: GAFF2 sp2 ring types on graph rings (≥5 sp2 ring atoms) ---
  const typedAromatic = new Set();
  for (let i = 0; i < n; i++) {
    const t = atoms[i]?.gaffType;
    if (typeof t === "string" && AROMATIC_GAFF_TYPES.has(t)) typedAromatic.add(i);
  }
  if (typedAromatic.size >= 5) {
    for (const ring of graphRings(adj, 5, 6, 300)) {
      if (ring.length < 5) continue;
      const sp2Count = ring.filter((i) => typedAromatic.has(i)).length;
      if (sp2Count < 5) continue;
      pushFrame(ring);
    }
  }

  // --- path 3: geometric aromatic rings (planarity + bond-length window) ---
  const RING_ELEMS = new Set(["C", "N", "O", "S"]);
  for (const ring of graphRings(adj, 5, 6, 300)) {
    if (!ring.every((i) => RING_ELEMS.has(String(atoms[i]?.element ?? "").toUpperCase()))) continue;
    if (!isGeometricAromatic(atoms, ring)) continue;
    // reject partial rings: every atom of a true 5-ring has exactly 2
    // neighbors inside the ring set; a 5-subset of a 6-ring contains atoms
    // whose third aromatic neighbor lies outside the candidate set.
    if (ring.length === 5) {
      const inSet = new Set(ring);
      let closed = true;
      for (const i of ring) {
        let deg = 0;
        for (const nbr of adj[i]) if (inSet.has(nbr)) deg++;
        if (deg !== 2) { closed = false; break; }
      }
      if (!closed) continue;
    }
    pushFrame(ring);
  }

  return frames;

  /**
   * Newell frame for an atom index set, ordered in cycle order around the
   * ring (required for a meaningful normal and the normal→atom force chain).
   * Protein name-set indices arrive unordered, so order them by walking the
   * ring adjacency; graph-ring indices are already in path order.
   */
  function makeFrame(atoms, idx) {
    const ordered = orderRingCycle(idx, adj);
    const pts = ordered.map(P);
    const centroid = [0, 0, 0];
    for (const p of pts) { centroid[0] += p[0] / pts.length; centroid[1] += p[1] / pts.length; centroid[2] += p[2] / pts.length; }
    let nx = 0, ny = 0, nz = 0;
    for (let t = 0; t < pts.length; t++) {
      const [x0, y0, z0] = pts[t];
      const [x1, y1, z1] = pts[(t + 1) % pts.length];
      nx += (y0 - y1) * (z0 + z1);
      ny += (z0 - z1) * (x0 + x1);
      nz += (x0 - x1) * (y0 + y1);
    }
    const L = Math.hypot(nx, ny, nz) || 1;
    return { centroid, normal: [nx / L, ny / L, nz / L], atomIdx: ordered };
  }
}

/**
 * Order a set of ring atom indices into cycle order via the ring adjacency.
 * Falls back to the input order when the set does not form a closed cycle in
 * the graph (still a valid Newell frame for near-planar sets).
 * @param {number[]} idx
 * @param {number[][]} adj
 * @returns {number[]}
 */
function orderRingCycle(idx, adj) {
  const inRing = new Set(idx);
  const nextIn = (i, prev) => {
    for (const nbr of adj[i]) if (inRing.has(nbr) && nbr !== prev) return nbr;
    return -1;
  };
  const ordered = [idx[0]];
  const seen = new Set([idx[0]]);
  let prev = -1, cur = idx[0];
  while (ordered.length < idx.length) {
    const nxt = nextIn(cur, prev);
    if (nxt === -1 || seen.has(nxt)) return idx; // degenerate — keep input order
    ordered.push(nxt); seen.add(nxt);
    prev = cur; cur = nxt;
  }
  return ordered;
}

/** Aromatic bond-length window (Å) and plane deviation (Å) — gaff2_mapper parity. */
const AROM_WINDOW = [1.28, 1.52];
const AROM_PLANE_DEV = 0.12;

/** Geometric aromaticity: all ring bonds in window + Newell-plane deviation. */
function isGeometricAromatic(atoms, ring) {
  const d = (i, j) => Math.hypot(atoms[i].x - atoms[j].x, atoms[i].y - atoms[j].y, atoms[i].z - atoms[j].z);
  for (let t = 0; t < ring.length; t++) {
    const r = d(ring[t], ring[(t + 1) % ring.length]);
    if (r < AROM_WINDOW[0] || r > AROM_WINDOW[1]) return false;
  }
  // Newell plane max deviation
  let nx = 0, ny = 0, nz = 0;
  for (let t = 0; t < ring.length; t++) {
    const a = atoms[ring[t]], b = atoms[ring[(t + 1) % ring.length]];
    nx += (a.y - b.y) * (a.z + b.z);
    ny += (a.z - b.z) * (a.x + b.x);
    nz += (a.x - b.x) * (a.y + b.y);
  }
  const L = Math.hypot(nx, ny, nz) || 1;
  nx /= L; ny /= L; nz /= L;
  const c0 = atoms[ring[0]];
  let mx = 0;
  for (const i of ring) {
    const a = atoms[i];
    mx = Math.max(mx, Math.abs((a.x - c0.x) * nx + (a.y - c0.y) * ny + (a.z - c0.z) * nz));
  }
  return mx < AROM_PLANE_DEV;
}

/**
 * Find simple 5- and 6-membered graph rings (DFS, deduplicated, capped).
 * @param {number[][]} adj
 * @param {number} minLen
 * @param {number} maxLen
 * @param {number} cap  max rings returned
 * @returns {number[][]}
 */
function graphRings(adj, minLen, maxLen, cap) {
  const rings = [];
  const seen = new Set();
  const n = adj.length;
  for (let start = 0; start < n && rings.length < cap; start++) {
    const path = [start];
    const inPath = new Set([start]);
    const dfs = (node, parent) => {
      if (rings.length >= cap) return;
      for (const nbr of adj[node]) {
        if (nbr === parent) continue;
        if (nbr === start) {
          if (path.length >= minLen && path.length <= maxLen) {
            const key = [...path].sort((a, b) => a - b).join("-");
            if (!seen.has(key)) {
              seen.add(key);
              rings.push([...path]);
            }
          }
        } else if (!inPath.has(nbr) && path.length < maxLen) {
          inPath.add(nbr);
          path.push(nbr);
          dfs(nbr, node);
          path.pop();
          inPath.delete(nbr);
        }
      }
    };
    dfs(start, -1);
  }
  return rings;
}

/* ------------------------------------------------------------------ */
/* Cation / halogen site lists                                        */
/* ------------------------------------------------------------------ */

/**
 * Build the cation site list for cation–π terms.
 * Protein: Lys NZ, Arg guanidinium center (CZ centroid proxy: CZ atom),
 * protonated His (both ND1+NE2 present → HIP-like; via chem/protonation
 * state names when the caller applied them). Ligands: N atoms with
 * positive parsed charge or GAFF n4/n3 type.
 * @param {Array<object>} atoms
 * @returns {Array<number>} cation atom indices
 */
export function buildCationList(atoms) {
  const out = [];
  for (let i = 0; i < atoms.length; i++) {
    const a = atoms[i];
    if (!a) continue;
    const el = String(a.element ?? "").toUpperCase();
    const name = String(a.atomName ?? "").trim().toUpperCase();
    const res = String(a.resName ?? "").trim().toUpperCase();
    if (a.isProtein) {
      if (res === "LYS" && name === "NZ") { out.push(i); continue; }
      if (res === "LYN" && name === "NZ") continue; // neutral Lys (protonation.js state)
      if (res === "ARG" && name === "CZ") { out.push(i); continue; }
      if (res === "ARN" && name === "CZ") continue; // neutral Arg
      if (res === "HIS" && name === "CZ") continue; // HIS has no CZ
      if ((res === "HIS" || res === "HIP") && (name === "ND1" || name === "NE2")) {
        // protonated only: HIP by name, or both nitrogens typed/placed
        if (res === "HIP") { out.push(i); continue; }
      }
    } else if (!a.isMetal && !a.isProtein) {
      // ligand / cofactor N cations: positive charge or GAFF2 n4/n3
      const q = Number(a.charge ?? 0);
      const t = String(a.gaffType ?? "");
      if (el === "N" && (q > 0.5 || t === "n4" || t === "n3")) out.push(i);
    }
  }
  return out;
}

/**
 * Build the halogen-bond donor list: Cl/Br/I atoms with exactly one bonded C
 * (the C partner recorded for the σ-hole axis). F is excluded by construction.
 * @param {Array<object>} atoms
 * @param {Array<[number,number]>} bonds
 * @returns {Array<{x:number, c:number}>} {x: halogen index, c: carbon index}
 */
export function buildHalogenList(atoms, bonds) {
  const adj = Array.from({ length: atoms.length }, () => []);
  for (const [i, j] of bonds) { adj[i].push(j); adj[j].push(i); }
  const out = [];
  for (let i = 0; i < atoms.length; i++) {
    const el = String(atoms[i]?.element ?? "").toUpperCase();
    if (el !== "CL" && el !== "BR" && el !== "I") continue;
    const carbons = adj[i].filter((j) => String(atoms[j]?.element ?? "").toUpperCase() === "C");
    if (carbons.length === 0) continue;
    out.push({ x: i, c: carbons[0] });
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* FD self-test helper (exported for tests)                           */
/* ------------------------------------------------------------------ */

/**
 * Central finite-difference gradient check for any weakint term.
 * Normals are renormalized after perturbation (unit-sphere chain rule,
 * matching the caller's use), so the analytic normal gradients must be
 * the projected (tangent) ones piStackEnergy returns.
 *
 * For normalized args the FD is taken on the UNIT SPHERE: the perturbed
 * point is vnorm(n ± h·e_k) and the step length is the actual chord
 * |vnorm(n+h·e_k) − vnorm(n−h·e_k)|/2 — dividing by the raw 2h would
 * silently mis-scale by 1/|n| (the prototype's own fdCheck only avoided
 * this because its demo normals were exactly unit). The analytic
 * comparison contracts the tangent gradient onto the same unit direction.
 *
 * @param {Function} fn  (...args) => {U, ...grads} — vector args perturbed
 * @param {Array} args  positional args; array args are perturbed componentwise
 * @param {Array<string|null>} gradKeys  gradient key per array arg (null = skip)
 * @param {object} [opts] {h=1e-5, tolAbs=1e-9, normalizeIdx=[] (renormalize after perturb)}
 * @returns {{maxRel:number, worst:string}}
 */
export function fdCheckWeak(fn, args, gradKeys, opts = {}) {
  const h = opts.h ?? 1e-5;
  const tolAbs = opts.tolAbs ?? 1e-9;
  const normalizeIdx = new Set(opts.normalizeIdx ?? []);
  const clone = (a) => (Array.isArray(a) ? [...a] : a);
  const base = fn(...args.map(clone));
  let maxRel = 0, worst = "";
  for (let ai = 0; ai < args.length; ai++) {
    if (!Array.isArray(args[ai])) continue;
    const gArr = base[gradKeys[ai]];
    if (!gArr) continue;
    for (let k = 0; k < 3; k++) {
      const a1 = args.map(clone); a1[ai][k] += h;
      for (const ni of normalizeIdx) a1[ni] = vnorm(a1[ni]);
      const a2 = args.map(clone); a2[ai][k] -= h;
      for (const ni of normalizeIdx) a2[ni] = vnorm(a2[ni]);
      const dU = fn(...a1).U - fn(...a2).U;
      if (normalizeIdx.has(ai)) {
        // unit-sphere FD: compare along the actual displacement direction
        const n1 = a1[ai], n2 = a2[ai];
        const dx = n1[0] - n2[0], dy = n1[1] - n2[1], dz = n1[2] - n2[2];
        const chord = Math.hypot(dx, dy, dz);
        if (chord < 1e-12) continue;
        const fdStep = dU / chord;                 // directional derivative · unit vector
        const an = (gArr[0] * dx + gArr[1] * dy + gArr[2] * dz) / chord;
        const denom = Math.abs(fdStep) + Math.abs(an);
        const rel = denom > tolAbs ? Math.abs(fdStep - an) / denom : 0;
        if (rel > maxRel) { maxRel = rel; worst = `${gradKeys[ai]}[${k}] fd=${fdStep.toPrecision(6)} an=${an.toPrecision(6)}`; }
      } else {
        const fd = dU / (2 * h);
        const denom = Math.abs(fd) + Math.abs(gArr[k]);
        const rel = denom > tolAbs ? Math.abs(fd - gArr[k]) / denom : 0;
        if (rel > maxRel) { maxRel = rel; worst = `${gradKeys[ai]}[${k}] fd=${fd.toPrecision(6)} an=${gArr[k].toPrecision(6)}`; }
      }
    }
  }
  return { maxRel, worst };
}

/* ------------------------------------------------------------------ */
/* Atom-level force accumulation (heavy.js hot-path pass)             */
/* ------------------------------------------------------------------ */

/**
 * Accumulate π-stack energy + forces of one ring pair into flat buffers.
 *
 * The term is parametrized on (centroid, normal) per ring; the chain rule
 * maps them back to atoms:
 *   ∂centroid/∂r_k = 1/m · 1
 *   normal n = Newell(ordered ring atoms) / |Newell| — its atom Jacobian is
 *   approximated by the in-plane "tangent lever" form
 *     f_k ⊖= d_k × gN / |Newell|,  d_k = r_{k+1} − r_{k−1}  (cycle order)
 *   which is exact for pure rotations of a rigid ring and conserves
 *   linear momentum (Σ d_k = 0) and the centroid (Σ r_k × d_k = 0 up to
 *   collinearity). Verified against atom-level central differences in
 *   tests/test_weakint.js (ring-pair FD < 1e-6).
 *
 * @param {ArrayLike<number>} pos flat 3n
 * @param {Float64Array} f 3n force accumulator (ENERGY gradients SUBTRACTED:
 *   f is force = −∇U)
 * @param {{atomIdx:number[]}} ringA
 * @param {{atomIdx:number[]}} ringB
 * @param {{eps?:number}} [params]
 * @returns {number} energy contribution (kcal/mol)
 */
export function piStackForces(pos, f, ringA, ringB, params = {}) {
  const idxA = ringA.atomIdx, idxB = ringB.atomIdx;
  const cA = ringCentroid(pos, idxA), cB = ringCentroid(pos, idxB);
  const nA = ringNormal(pos, idxA), nB = ringNormal(pos, idxB);
  const res = piStackEnergy(cA, nA, cB, nB, params);
  if (res.U === 0 && res.debug?.cut) return 0;
  // centroid gradient → atoms (equal split), FORCE = −gradient
  const ma = 1 / idxA.length, mb = 1 / idxB.length;
  for (let t = 0; t < idxA.length; t++) {
    const i3 = 3 * idxA[t];
    f[i3] -= res.gA[0] * ma; f[i3 + 1] -= res.gA[1] * ma; f[i3 + 2] -= res.gA[2] * ma;
  }
  for (let t = 0; t < idxB.length; t++) {
    const i3 = 3 * idxB[t];
    f[i3] -= res.gB[0] * mb; f[i3 + 1] -= res.gB[1] * mb; f[i3 + 2] -= res.gB[2] * mb;
  }
  // normal gradients → atom torques (Newell lever form), FORCE = −gradient
  applyNormalGrad(pos, f, idxA, res.gNA);
  applyNormalGrad(pos, f, idxB, res.gNB);
  return res.U;
}

/**
 * Accumulate cation–π energy + forces for one cation–ring pair.
 * The ring normal is held fixed on MD timescales (R3 §1b); centroid gradient
 * splits over ring atoms, cation gets the full gradient.
 * @param {ArrayLike<number>} pos
 * @param {Float64Array} f
 * @param {number} catIdx cation atom index
 * @param {{atomIdx:number[]}} ring
 * @param {{eps?:number, r0?:number, w?:number}} [params]
 * @returns {number} energy
 */
export function cationPiForces(pos, f, catIdx, ring, params = {}) {
  const idx = ring.atomIdx;
  const centroid = ringCentroid(pos, idx);
  const normal = ringNormal(pos, idx);
  const catPos = [pos[3 * catIdx], pos[3 * catIdx + 1], pos[3 * catIdx + 2]];
  const res = cationPiEnergy(catPos, centroid, normal, params);
  if (res.U === 0 && res.debug?.cut) return 0;
  const c3 = 3 * catIdx;
  f[c3] -= res.gCat[0]; f[c3 + 1] -= res.gCat[1]; f[c3 + 2] -= res.gCat[2];
  const m = 1 / idx.length;
  for (const i of idx) {
    const i3 = 3 * i;
    f[i3] -= res.gRing[0] * m; f[i3 + 1] -= res.gRing[1] * m; f[i3 + 2] -= res.gRing[2] * m;
  }
  return res.U;
}

/**
 * Accumulate halogen σ-hole energy + forces for one C–X···D triple.
 * @param {ArrayLike<number>} pos
 * @param {Float64Array} f
 * @param {number} cIdx carbon index
 * @param {number} xIdx halogen index
 * @param {number} dIdx acceptor index
 * @param {{eps?:number, r0?:number, w?:number}} [params]
 * @returns {number} energy
 */
export function halogenForces(pos, f, cIdx, xIdx, dIdx, params = {}) {
  const p = (i) => [pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]];
  const res = halogenEnergy(p(cIdx), p(xIdx), p(dIdx), params);
  if (res.U === 0 && res.debug?.cut) return 0;
  const c3 = 3 * cIdx, x3 = 3 * xIdx, d3 = 3 * dIdx;
  f[c3] -= res.gC[0]; f[c3 + 1] -= res.gC[1]; f[c3 + 2] -= res.gC[2];
  f[x3] -= res.gX[0]; f[x3 + 1] -= res.gX[1]; f[x3 + 2] -= res.gX[2];
  f[d3] -= res.gD[0]; f[d3 + 1] -= res.gD[1]; f[d3 + 2] -= res.gD[2];
  return res.U;
}

/** Centroid of an atom index list from flat positions. */
function ringCentroid(pos, idx) {
  const c = [0, 0, 0];
  const m = 1 / idx.length;
  for (const i of idx) {
    const i3 = 3 * i;
    c[0] += pos[i3] * m; c[1] += pos[i3 + 1] * m; c[2] += pos[i3 + 2] * m;
  }
  return c;
}

/** Newell unit normal of an atom index list (cycle order) from flat positions. */
function ringNormal(pos, idx) {
  let nx = 0, ny = 0, nz = 0;
  const m = idx.length;
  for (let t = 0; t < m; t++) {
    const i3 = 3 * idx[t], j3 = 3 * idx[(t + 1) % m];
    const x0 = pos[i3], y0 = pos[i3 + 1], z0 = pos[i3 + 2];
    const x1 = pos[j3], y1 = pos[j3 + 1], z1 = pos[j3 + 2];
    nx += (y0 - y1) * (z0 + z1);
    ny += (z0 - z1) * (x0 + x1);
    nz += (x0 - x1) * (y0 + y1);
  }
  const L = Math.hypot(nx, ny, nz) || 1;
  return [nx / L, ny / L, nz / L];
}

/**
 * Map a unit-normal gradient gN to per-atom forces via the Newell lever:
 *   F_k = − gN × d_k / |Newell|,  d_k = r_{k+1} − r_{k−1}
 * (energy-gradient subtraction folded in: forces = −∇U).
 * @param {ArrayLike<number>} pos
 * @param {Float64Array} f
 * @param {number[]} idx ring atom indices (cycle order)
 * @param {number[]} gN normal gradient
 */
function applyNormalGrad(pos, f, idx, gN) {
  const m = idx.length;
  if (m < 3) return;
  // |Newell| for this ring (recomputed from positions)
  let nx = 0, ny = 0, nz = 0;
  for (let t = 0; t < m; t++) {
    const i3 = 3 * idx[t], j3 = 3 * idx[(t + 1) % m];
    const x0 = pos[i3], y0 = pos[i3 + 1], z0 = pos[i3 + 2];
    const x1 = pos[j3], y1 = pos[j3 + 1], z1 = pos[j3 + 2];
    nx += (y0 - y1) * (z0 + z1);
    ny += (z0 - z1) * (x0 + x1);
    nz += (x0 - x1) * (y0 + y1);
  }
  const L = Math.hypot(nx, ny, nz);
  if (L < 1e-9) return;
  const s = 1 / L;
  for (let k = 0; k < m; k++) {
    const i3 = 3 * idx[k];
    const p3 = 3 * idx[(k - 1 + m) % m]; // r_{k−1}
    const n3 = 3 * idx[(k + 1) % m];    // r_{k+1}
    const dx = pos[n3] - pos[p3], dy = pos[n3 + 1] - pos[p3 + 1], dz = pos[n3 + 2] - pos[p3 + 2];
    // ∂U/∂r_k = (d_k × gN)/L  ⇒  F_k = −∂U/∂r_k = (gN × d_k)/L
    f[i3] += (gN[1] * dz - gN[2] * dy) * s;
    f[i3 + 1] += (gN[2] * dx - gN[0] * dz) * s;
    f[i3 + 2] += (gN[0] * dy - gN[1] * dx) * s;
  }
}
