/**
 * cg/grid.js — the CG uniform spatial grid: its numeric cell-key codec and the
 * preallocated buffers both grid kernels share. Split out of
 * src/forcefield.js; moved verbatim.
 *
 * WHAT THIS FILE OWNS
 * -------------------
 * 1. allocBuffers() — every preallocated Float64Array / Map the compute path
 *    reuses, including the EEF1-lite desolvation pair records. This is the
 *    ZERO-ALLOC contract in one place (docs/CG_HEAVY.md G65).
 * 2. The cell-key codec: _encodeCell / _cellKey / _decode{X,Y,Z}. Numeric
 *    keys, no string allocs, and the non-finite sentinel that keeps the Map
 *    from growing without bound (the 1HVR hang).
 *
 * WHY IT IS ITS OWN MODULE
 * ------------------------
 * The codec and the kernel that consumes it (ff-repulsion.js, ff-binding.js)
 * must agree on the key encoding exactly. While they lived apart inside one
 * class, a change to the encode/decode constants was a search. Here the
 * codec and the buffers it indexes are the same file, and the kernels read
 * the key type from here by import — so a mismatch is a type error at the
 * import site rather than a silently-empty neighbour bucket.
 *
 * The buffers live here rather than in each kernel because ff-repulsion.js
 * and ff-binding.js share _grid / _dens / _bpA/J/R: the binding pass writes
 * the density records the desolvation pass reads, and neither owns them.
 */

/**
 * Allocate every preallocated buffer the CG compute path uses. Runs LAST in
 * the constructor, after n / nLigAtoms / _repSigma are final.
 * @param {object} ff   the ForceField under construction
 * @param {object} par  constructor params (unused today; kept for symmetry
 *                      with the other constructors so a future par-driven
 *                      scratch option lands in the obvious place)
 */
export function allocBuffers(ff, par) {
  const nProt = ff.nProt;
  // Scratch / grid buffers (spatial hash with numeric keys — GC-free)
  ff.forces = new Float64Array(ff.n * 3);
  ff.rcRep = (2 ** (1 / 6)) * ff.sigmaRep; // r_e = 2^(1/6) σ ≈ 5.61 Å
  ff._cell = ff.rcRep;                    // cell size ≥ repulsive range
  ff._grid = new Map();                     // hashCellKey -> bead index list
  ff._gridB = new Map();                    // protein-only grid for binding pass
  ff.energy = 0; // last computed potential energy (kcal/mol)
  ff.bindingU = 0;   // protein–ligand nonbonded energy (last compute)
  ff.desolvU = 0;    // EEF1-lite burial/desolvation energy (last compute)
  // Loop-2 S4 (R4 §5 item 1, R6 §5): per-term binding-accumulator opt-in.
  // trackTerms = true fills bindLJU/bindCoulU/bindHBU/desolvU + the
  // bindU vector {lj, coul, hb, desolv} each compute() — consumed by the
  // BindLog energy channel (main.js tick). DEFAULT false: kernel skips
  // all tracker branches → bit-identical to pre-S4, zero overhead.
  ff.trackTerms = false;
  ff.bindLJU = 0; ff.bindCoulU = 0; ff.bindHBU = 0;
  ff.bindU = { lj: 0, coul: 0, hb: 0, desolv: 0 };

  // EEF1-lite desolvation scratch (reused every compute call — GC-free).
  // dens[a]  = soft protein-occupancy count n_a (Pass 1 accumulation)
  // dBdn[a]  = dB/dn = exp(−n_a/3)/3 at the current density (Pass 2 lookup)
  // bpA/bpJ/bpR = flat pair records (ligand table index, protein bead index,
  //   distance r) written in Pass 1 so Pass 2 adds forces without a re-scan.
  //   Parallel Float64Arrays are chosen over object arrays to stay allocation-
  //   free. The worst case is every ligand atom paired with every protein bead,
  //   so nProt·nLigAtoms is a hard size cap and no growth is ever needed.
  ff._dens = new Float64Array(Math.max(1, ff.nLigAtoms));
  ff._dBdn = new Float64Array(Math.max(1, ff.nLigAtoms));
  const maxBP = ff.nLigAtoms ? nProt * ff.nLigAtoms : 1;
  ff._bpA = new Float64Array(maxBP);
  ff._bpJ = new Float64Array(maxBP);
  ff._bpR = new Float64Array(maxBP);
}

/* ------------------------------------------------------------------ *
 *  Numeric spatial-hash cell keys (no string allocs)                *
 * ------------------------------------------------------------------ */

/** Pack integer cell coordinates into one numeric key. */
export function encodeCell(cx, cy, cz) { return ((cx + 2048) * 4096 + (cy + 2048)) * 4096 + (cz + 2048); }

/**
 * Cell key for a world coordinate.
 *
 * Non-finite coordinates must never reach the grid: Math.floor(NaN) is
 * NaN and the Map would grow unboundedly with a fresh garbage key on
 * every compute call (the 1HVR hang). Route everything non-finite to a
 * single sentinel cell at the grid corner — the pair distance tests will
 * then reject those particles (r2 is NaN, comparisons false).
 */
export function cellKey(x, y, z, cell) {
  if (!Number.isFinite(x + y + z)) return 0;
  return encodeCell(Math.floor(x / cell), Math.floor(y / cell), Math.floor(z / cell));
}

export function decodeX(k) { return Math.floor(k / (4096 * 4096)) - 2048; }
export function decodeY(k) { return Math.floor(k / 4096) % 4096 - 2048; }
export function decodeZ(k) { return k % 4096 - 2048; }