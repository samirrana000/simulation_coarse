/**
 * spatial-grid.js — High-performance 3D uniform spatial grid & cell list.
 *
 * Replaces O(N^2) pairwise non-bonded scans with O(N) neighbor searches.
 * Uses typed arrays with head/next linked lists to avoid any memory allocation during simulation loops.
 * G66 — Verlet skin 2Å, rebuild every 10 steps, 20% cut — aspirational target (currently rebuilds every compute() call, skin not yet implemented)
 */

export class SpatialGrid {
  /**
   * @param {number} [cellSize=8.5] Grid cell size in Å (must be >= cutoff)
   * @param {number} [maxAtoms=16384] Maximum capacity for pre-allocated arrays
   */
  constructor(cellSize = 8.5, maxAtoms = 16384) {
    this.cellSize = cellSize;
    this.invCell = 1.0 / cellSize;
    this.maxAtoms = maxAtoms;

    // Hash table size (prime for good distribution)
    this.tableSize = 4093;
    this.head = new Int32Array(this.tableSize);
    this.next = new Int32Array(maxAtoms);
    this.cellCoords = new Int32Array(maxAtoms * 3);
  }

  /**
   * Hash 3D integer cell coordinates into table index.
   */
  hashCell(cx, cy, cz) {
    let h = (cx * 73856093) ^ (cy * 19349663) ^ (cz * 83492791);
    h = h % this.tableSize;
    return h < 0 ? h + this.tableSize : h;
  }

  /**
   * Build the cell list from particle positions. Zero allocations.
   * G66 — Verlet skin 2Å, rebuild every 10 steps, 20% cut — aspirational; current implementation rebuilds every call (no skin)
   *
   * @param {Float64Array} pos Flat Cartesian coordinates [x0,y0,z0, x1,y1,z1, ...]
   * @param {number} n         Number of particles
   */
  build(pos, n) {
    if (n > this.maxAtoms) {
      this.maxAtoms = Math.max(n * 2, this.maxAtoms * 2);
      this.next = new Int32Array(this.maxAtoms);
      this.cellCoords = new Int32Array(this.maxAtoms * 3);
    }

    this.head.fill(-1);
    const inv = this.invCell;

    for (let i = 0; i < n; i++) {
      const i3 = 3 * i;
      const cx = Math.floor(pos[i3] * inv);
      const cy = Math.floor(pos[i3 + 1] * inv);
      const cz = Math.floor(pos[i3 + 2] * inv);

      this.cellCoords[i3] = cx;
      this.cellCoords[i3 + 1] = cy;
      this.cellCoords[i3 + 2] = cz;

      const cellIdx = this.hashCell(cx, cy, cz);
      this.next[i] = this.head[cellIdx];
      this.head[cellIdx] = i;
    }
  }

  /**
   * Iterate over all unique pairs (i < j) within cutoff distance.
   *
   * @param {Float64Array} pos
   * @param {number} n
   * @param {number} cutoff
   * @param {(i: number, j: number, dx: number, dy: number, dz: number, r2: number, r: number) => void} callback
   */
  forEachPair(pos, n, cutoff, callback) {
    const cut2 = cutoff * cutoff;

    // 14 unique neighbor cell offsets for half-neighborhood symmetric pair checks
    const offsets = [
      [0, 0, 0],
      [1, 0, 0], [-1, 1, 0], [0, 1, 0], [1, 1, 0],
      [-1, -1, 1], [0, -1, 1], [1, -1, 1],
      [-1, 0, 1], [0, 0, 1], [1, 0, 1],
      [-1, 1, 1], [0, 1, 1], [1, 1, 1],
    ];

    for (let i = 0; i < n; i++) {
      const i3 = 3 * i;
      const xi = pos[i3], yi = pos[i3 + 1], zi = pos[i3 + 2];
      const cx = this.cellCoords[i3];
      const cy = this.cellCoords[i3 + 1];
      const cz = this.cellCoords[i3 + 2];

      for (let o = 0; o < 14; o++) {
        const off = offsets[o];
        const ncx = cx + off[0];
        const ncy = cy + off[1];
        const ncz = cz + off[2];
        const cellIdx = this.hashCell(ncx, ncy, ncz);

        let j = this.head[cellIdx];
        while (j !== -1) {
          if (o === 0) {
            // Same cell: only take j > i
            if (j > i) {
              const j3 = 3 * j;
              const dx = pos[j3] - xi;
              const dy = pos[j3 + 1] - yi;
              const dz = pos[j3 + 2] - zi;
              const r2 = dx * dx + dy * dy + dz * dz;
              if (r2 < cut2) {
                callback(i, j, dx, dy, dz, r2, Math.sqrt(r2) || 1e-12);
              }
            }
          } else {
            // Neighbor cell: check cell coords to handle hash collisions
            const j3 = 3 * j;
            if (this.cellCoords[j3] === ncx &&
                this.cellCoords[j3 + 1] === ncy &&
                this.cellCoords[j3 + 2] === ncz) {
              const dx = pos[j3] - xi;
              const dy = pos[j3 + 1] - yi;
              const dz = pos[j3 + 2] - zi;
              const r2 = dx * dx + dy * dy + dz * dz;
              if (r2 < cut2) {
                callback(i, j, dx, dy, dz, r2, Math.sqrt(r2) || 1e-12);
              }
            }
          }
          j = this.next[j];
        }
      }
    }
  }
}
