/**
 * heavy/topology.js — covalent topology perception from geometry, metal
 * coordination, and the flat typed-array lists the kernels iterate.
 *
 * PART OF THE src/heavy.js SPLIT (2026-10). This is the algorithmic core of
 * the old monolith, and none of it is shared with CG mode.
 *
 *   topologyBondRows        the bond double-loop over a ROW RANGE. Range-shaped
 *                           on purpose: buildTopology (whole range, sync) and
 *                           buildTopologyChunked (row slices with event-loop
 *                           yields) share it so both are bit-identical.
 *   findTopologyRings       5- and 6-membered ring DFS
 *   finishTopology          angles + propers + ring distance braces (shared tail)
 *   buildTopology           the synchronous builder
 *   buildTopologyChunked    the cancellable, progress-reporting builder
 *   buildMetalCoordination  nearest N/O/S donors per metal, coordN of them
 *   buildLists              flat [i,j,r0] / [i,j,k,theta0] / [i,j,k,l,phi0] lists
 *   angleAt, improperAngleFlat  reference geometry at the current pose
 */
import { METAL_ELEMENT, METAL_ELEMENT_DEFAULT, COVALENT_RADIUS, BOND_SLACK } from "../ff-params.js";
import { ignore } from "../errors.js";

/**
 * Build covalent topology from geometry.
 *
 * Covalent radii source: CSD surveys (Allen et al., Cordero et al. Dalton
 * 2008) / Bondi 1964 — see src/ff-params.js:COVALENT_RADIUS table. Two heavy
 * atoms i,j are bonded iff  r_ij < BOND_SLACK*(r_cov(i)+r_cov(j))  AND
 * r_ij < 2.2 Å hard cap. BOND_SLACK=1.15 (15% slack) absorbs thermal spread;
 * cap 2.2 Å rejects spurious long contacts (Ca–N 2.9 Å) while correctly
 * capturing S–S disulfide 2.04 Å (S 1.02+1.02=2.04, 2.04*1.15=2.35→capped 2.20,
 * still passes). Validated on crambin 1CRN (3 SSBOND): tests/test_topology.js
 * asserts ≥3 S–S bonds recovered and no Ca–N 2.9 Å false bond.
 */
/**
 * Bond double-loop over the row range [i0, i1): two heavy atoms i,j are
 * bonded iff r_ij < BOND_SLACK*(r_cov(i)+r_cov(j)) AND r_ij < 2.2 Å.
 * Shared by buildTopology (full range, sync) and buildTopologyChunked
 * (row slices with event-loop yields) so both produce bit-identical output.
 * @param {Array} atoms heavy-atom records
 * @param {Array} bonds pair list to append [i, j] into
 * @param {Array<Array<number>>} nbond adjacency lists to append into
 * @param {number} i0 first row (inclusive)
 * @param {number} i1 one-past-last row
 */
function topologyBondRows(atoms, bonds, nbond, i0, i1) {
  const n = atoms.length;
  for (let i = i0; i < i1; i++) {
    const A = atoms[i];
    if (A.isMetal) continue;
    const rA = COVALENT_RADIUS[A.element] ?? 0.77;
    for (let j = i + 1; j < n; j++) {
      const B = atoms[j];
      if (B.isMetal) continue;
      const rB = COVALENT_RADIUS[B.element] ?? 0.77;
      const dx = A.x - B.x, dy = A.y - B.y, dz = A.z - B.z;
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (r < BOND_SLACK * (rA + rB) && r < 2.2) {
        bonds.push([i, j]);
        nbond[i].push(j);
        nbond[j].push(i);
      }
    }
  }
}

/**
 * Ring detection: 5- and 6-membered rings whose cross-ring pairs become
 * distance braces (shared by both topology builders).
 * @param {number} n atom count
 * @param {Array<Array<number>>} nbond adjacency lists
 * @returns {Array<Array<number>>} detected rings (atom-index paths)
 */
function findTopologyRings(n, nbond) {
  const rings = [];
  const visited = new Set();
  for (let start = 0; start < n; start++) {
    const path = [start];
    const inPath = new Set([start]);
    const dfs = (node, parent) => {
      for (const nbr of nbond[node]) {
        if (nbr === parent) continue;
        if (nbr === start) {
          if (path.length === 5 || path.length === 6) {
            const sortedKey = [...path].sort((a,b)=>a-b).join("-");
            if (!visited.has(sortedKey)) {
              visited.add(sortedKey);
              rings.push([...path]);
            }
          }
        } else if (!inPath.has(nbr) && path.length < 6) {
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

/**
 * Angles + propers + ring braces from a finished bond graph (shared tail).
 * @param {Array} atoms heavy-atom records
 * @param {Array} bonds pair list (ring braces appended here)
 * @param {Array<Array<number>>} nbond adjacency lists
 * @returns {{bonds:Array, angles:Array, propers:Array, impropers:Array}}
 */
function finishTopology(atoms, bonds, nbond) {
  const n = atoms.length;
  const angles = [];
  const angleSet = new Set();
  for (let j = 0; j < n; j++) {
    const nb = nbond[j];
    for (let a = 0; a < nb.length; a++) {
      for (let b = a + 1; b < nb.length; b++) {
        const i = nb[a], k = nb[b];
        const key = i < k ? `${i}-${j}-${k}` : `${k}-${j}-${i}`;
        if (angleSet.has(key)) continue;
        angleSet.add(key);
        angles.push([i, j, k]);
      }
    }
  }

  const propers = [];
  const properSet = new Set();
  for (let j = 0; j < n; j++) {
    for (const i of nbond[j]) {
      for (const k of nbond[j]) {
        if (k === i) continue;
        for (const l of nbond[k]) {
          if (l === i || l === j) continue;
          const key = [i, j, k, l].join("-");
          if (properSet.has(key)) continue;
          properSet.add(key);
          propers.push([i, j, k, l]);
        }
      }
    }
  }

  const impropers = [];

  const detectedRings = findTopologyRings(n, nbond);
  for (const ring of detectedRings) {
    const len = ring.length;
    // Cross-ring distance restraints to rigidly maintain planar geometry
    for (let i = 0; i < len; i++) {
      for (let j = i + 2; j < len; j++) {
        if (i === 0 && j === len - 1) continue;
        const idxA = ring[i], idxB = ring[j];
        bonds.push([idxA, idxB]);
      }
    }
  }

  return { bonds, angles, propers, impropers };
}

export function buildTopology(atoms) {
  const n = atoms.length;
  const bonds = [];
  const nbond = Array.from({ length: n }, () => []);

  topologyBondRows(atoms, bonds, nbond, 0, n);

  return finishTopology(atoms, bonds, nbond);
}

/**
 * Chunked heavy-topology builder (FP5: heavy-build progress UX).
 * Same bond loop as buildTopology, run in row slices of `chunkRows` with a
 * `setTimeout(0)` yield between slices so progress captions paint and the
 * event loop stays free (heartbeat timers fire mid-build). Angles/propers/
 * ring braces run once at the end via the shared finishTopology tail, so the
 * result is bit-identical to buildTopology. Cancellation is cooperative:
 * `isCancelled()` is polled at every slice boundary (and once before the
 * tail); on cancel it throws `Error("heavy build cancelled")`.
 * Zero deps; sync tail slices stay < ~500 ms at bundled-system sizes.
 * @param {Array} atoms heavy-atom records
 * @param {object} [opts]
 * @param {number} [opts.chunkRows=128] bond-loop rows per slice
 * @param {(done:number, total:number) => void} [opts.onProgress] per-slice callback (never throws the build)
 * @param {() => boolean} [opts.isCancelled] cooperative-cancel poll
 * @returns {Promise<{bonds:Array, angles:Array, propers:Array, impropers:Array}>}
 */
export async function buildTopologyChunked(atoms, opts = {}) {
  const chunkRows = Math.max(1, Math.floor(opts.chunkRows ?? 128));
  const onProgress = typeof opts.onProgress === "function" ? opts.onProgress : null;
  const isCancelled = typeof opts.isCancelled === "function" ? opts.isCancelled : null;
  const n = atoms.length;
  const bonds = [];
  const nbond = Array.from({ length: n }, () => []);

  for (let i0 = 0; i0 < n; i0 += chunkRows) {
    if (isCancelled && isCancelled()) throw new Error("heavy build cancelled");
    const i1 = Math.min(n, i0 + chunkRows);
    topologyBondRows(atoms, bonds, nbond, i0, i1);
    if (onProgress) {
      try { onProgress(i1, n); } catch (e) { ignore(e, "onProgress@buildTopologyChunked", "progress must never fail the build; the chunk loop owns the result"); }
    }
    if (i1 < n) await new Promise((r) => setTimeout(r, 0));
  }
  if (isCancelled && isCancelled()) throw new Error("heavy build cancelled");

  return finishTopology(atoms, bonds, nbond);
}

export function buildMetalCoordination(atoms) {
  const pairs = [];
  for (let i = 0; i < atoms.length; i++) {
    const M = atoms[i];
    if (!M.isMetal) continue;
    const p = METAL_ELEMENT[M.element] ?? METAL_ELEMENT_DEFAULT;
    const donors = [];
    for (let j = 0; j < atoms.length; j++) {
      if (j === i) continue;
      const D = atoms[j];
      if (D.isMetal) continue;
      if (!(D.element === "N" || D.element === "O" || D.element === "S")) continue;
      const dx = M.x - D.x, dy = M.y - D.y, dz = M.z - D.z;
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (r < p.coordR) donors.push({ j, r });
    }
    donors.sort((a, b) => a.r - b.r);
    const take = Math.min(p.coordN, donors.length);
    for (let k = 0; k < take; k++) pairs.push([i, donors[k].j]);
  }
  return pairs;
}

export function buildLists(atoms, topo, coordPairs) {
  const bonds = [];
  for (const [i, j] of topo.bonds) {
    const A = atoms[i], B = atoms[j];
    const r0 = Math.hypot(A.x - B.x, A.y - B.y, A.z - B.z);
    bonds.push(i, j, r0);
  }

  const angles = [];
  for (const [i, j, k] of topo.angles) {
    const th0 = angleAt(atoms, i, j, k);
    angles.push(i, j, k, th0);
  }

  const impropers = [];
  for (const [i, j, k, l] of topo.impropers) {
    const phi0 = improperAngleFlat(atoms, i, j, k, l);
    impropers.push(i, j, k, l, phi0);
  }

  const propers = [];
  for (const [i, j, k, l] of topo.propers) {
    const phi0 = improperAngleFlat(atoms, i, j, k, l);
    propers.push(i, j, k, l, phi0);
  }

  const coord = [];
  for (const [i, j] of coordPairs) {
    const A = atoms[i], B = atoms[j];
    const r0 = Math.hypot(A.x - B.x, A.y - B.y, A.z - B.z);
    coord.push(i, j, r0);
  }

  return { bonds, angles, impropers, propers, coord, holoSprings: new Float64Array(0), nHolo: 0 };
}

export function angleAt(atoms, i, j, k) {
  const A = atoms[i], B = atoms[j], C = atoms[k];
  const ax = A.x - B.x, ay = A.y - B.y, az = A.z - B.z;
  const bx = C.x - B.x, by = C.y - B.y, bz = C.z - B.z;
  const la = Math.hypot(ax, ay, az) || 1e-12;
  const lb = Math.hypot(bx, by, bz) || 1e-12;
  let c = (ax * bx + ay * by + az * bz) / (la * lb);
  c = Math.min(1, Math.max(-1, c));
  return Math.acos(c);
}

function improperAngleFlat(atoms, i, j, k, l) {
  const A = atoms[i], B = atoms[j], C = atoms[k], D = atoms[l];
  // vectors ij, jk, kl
  const ij_x = B.x - A.x, ij_y = B.y - A.y, ij_z = B.z - A.z;
  const jk_x = C.x - B.x, jk_y = C.y - B.y, jk_z = C.z - B.z;
  const kl_x = D.x - C.x, kl_y = D.y - C.y, kl_z = D.z - C.z;

  const mx = ij_y * jk_z - ij_z * jk_y;
  const my = ij_z * jk_x - ij_x * jk_z;
  const mz = ij_x * jk_y - ij_y * jk_x;

  const nx = jk_y * kl_z - jk_z * kl_y;
  const ny = jk_z * kl_x - jk_x * kl_z;
  const nz = jk_x * kl_y - jk_y * kl_x;

  const m_len = Math.hypot(mx, my, mz) || 1e-12;
  const n_len = Math.hypot(nx, ny, nz) || 1e-12;
  const jk_len = Math.hypot(jk_x, jk_y, jk_z) || 1e-12;

  const cos_phi = (mx * nx + my * ny + mz * nz) / (m_len * n_len);
  const mxn_x = my * nz - mz * ny;
  const mxn_y = mz * nx - mx * nz;
  const mxn_z = mx * ny - my * nx;
  const sin_phi = (mxn_x * jk_x + mxn_y * jk_y + mxn_z * jk_z) / (m_len * n_len * jk_len);

  return Math.atan2(sin_phi, Math.max(-1, Math.min(1, cos_phi)));
}
