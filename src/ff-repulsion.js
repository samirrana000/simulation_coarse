/**
 * ff-repulsion.js — excluded-volume (repulsive-only 12-6 LJ) kernel for the
 * ForceField class in forcefield.js (item 5 modularization).
 *
 * Verbatim body of the original ForceField._repulsion with `this` replaced
 * by an explicit `ff` parameter; the spatial-hash helpers (_cellKey,
 * _encodeCell, _decodeX/Y/Z, _pairKey) stay on the ForceField instance and
 * are invoked through ff. forcefield.js wraps this with the same
 * _repulsion(pos, f) signature, so no call site changes.
 */

/**
 * Repulsive-only 12-6 LJ for all non-bonded/non-contact pairs (implicit
 * solvent: beads cannot overlap). Spatial hash grid with numeric keys.
 *   U_rep(r) = ε[(r_e/r)¹² − 2(r_e/r)⁶ + 1],   r < r_e = 2^{1/6}σ
 *   F(r)     = 24ε[2(σ/r)¹² − (σ/r)⁶]/r² — see algebra in code body.
 */
export function repulsion(ff, pos, f) {
  const n = ff.n, cell = ff._cell, rc = ff.rcRep, eps = ff.epsRep;
  const grid = ff._grid, sig = ff._repSigma;
  // clear lists in place (avoid GC churn)
  for (const arr of grid.values()) arr.length = 0;

  // insert beads into grid cells
  for (let i = 0; i < n; i++) {
    const key = ff._cellKey(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2], cell);
    let arr = grid.get(key);
    if (!arr) grid.set(key, (arr = []));
    arr.push(i);
  }

  let U = 0;
  // iterate over cells; for each, self-pairs handled separately
  for (const [key, arr] of grid) {
    const cx = ff._decodeX(key), cy = ff._decodeY(key), cz = ff._decodeZ(key);
    // visit 14 of 27 neighbours (half-shell) to count each pair once
    for (let ox = -1; ox <= 1; ox++)
      for (let oy = -1; oy <= 1; oy++)
        for (let oz = -1; oz <= 1; oz++) {
          if (ox < 0 || (ox === 0 && oy < 0) || (ox === 0 && oy === 0 && oz < 0)) continue;
          const same = ox === 0 && oy === 0 && oz === 0;
          const nbr = same ? arr : grid.get(ff._encodeCell(cx + ox, cy + oy, cz + oz));
          if (!nbr) continue;
          for (let ai = 0; ai < arr.length; ai++) {
            const i = arr[ai], ix = 3 * i;
            for (let bi = same ? ai + 1 : 0; bi < nbr.length; bi++) {
              const j = nbr[bi];
              const pk = ff._pairKey(i, j);
              if ((i < ff.nProt) !== (j < ff.nProt)) continue;
              if (ff._excluded.has(pk)) continue;
              const jx = 3 * j;
              const dx = pos[jx] - pos[ix], dy = pos[jx + 1] - pos[ix + 1], dz = pos[jx + 2] - pos[ix + 2];
              const r2 = dx * dx + dy * dy + dz * dz;
              // per-pair σ = arithmetic mean of the two bead sizes; its own
              // r_e = 2^(1/6)σ per pair (keeps U(contact edge) = 0 exactly)
              const s = 0.5 * (sig[i] + sig[j]);
              const re = rc * (s / ff.sigmaRep);
              const re2 = re * re;
              if (r2 >= re2 || r2 < 1e-10) continue;
              const r = Math.sqrt(r2);
              const sr = s / r;
              const sr2 = sr * sr;
              const sr6 = sr2 * sr2 * sr2;              // (σ/r)^6
              const uLJ = sr6 * sr6 - sr6;              // [(σ/r)¹² − (σ/r)⁶]
              // energy: ε·(4·uLJ + 1)   (shifted so U(rcRep) = 0 exactly)
              U += eps * (4 * uLJ + 1);
              // −dU/dr = 24ε[2(σ/r)¹² − (σ/r)⁶]/r ; vector form divides by r²
              const fm = (24 * eps * (2 * sr6 * sr6 - sr6)) / r2;
              // i feels −∇_i U → away from j (repulsive out of overlap)
              f[ix] -= fm * dx; f[ix + 1] -= fm * dy; f[ix + 2] -= fm * dz;
              f[jx] += fm * dx; f[jx + 1] += fm * dy; f[jx + 2] += fm * dz;
            }
          }
        }
  }
  return U;
}
