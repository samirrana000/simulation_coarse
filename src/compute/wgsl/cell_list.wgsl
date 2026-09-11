// cell_list.wgsl — GPU spatial-hash cell list for MD non-bonded neighbor search.
//
// Purpose: bin N atoms into a uniform 3D grid with cell edge == non-bonded
// cutoff (RCUT = 12.0 A, see NONBONDED_RCUT below) so the force kernel only
// walks the 27 neighboring cells instead of all N^2 pairs. Two-pass counting
// sort (count -> CPU prefix-sum -> fill); both passes are one thread per atom.
//
// Units: positions in Angstrom. All arithmetic f32.
//
// Pass 1 (entry point `bin_count`):
//   per atom i: c = flatten(floor((pos - gridMin) / cellSize)); atomCell[i] = c;
//               atomicAdd(&cellCounts[c], 1u)
//   host then prefix-sums cellCounts -> cellStart (length numCells), zeroes
//   cellCursor, and dispatches pass 2 with the same bindings plus the
//   fill-specific buffers.
//
// Pass 2 (entry point `bin_fill`):
//   per atom i: c = atomCell[i]; s = atomicAdd(&cellCursor[c], 1u);
//               sortedIdx[cellStart[c] + s] = i
//   after dispatch, atoms in cell c occupy sortedIdx[cellStart[c] .. +count).
//
// Binding layout (group 0, shared by both entry points):
//   @binding(0) pos         : array<vec4<f32>>, read-only.
//                             xyz = position (A); w unused (padding for 16 B alignment).
//   @binding(1) params      : CellParams uniform (48 B, see struct below).
//   @binding(2) atomCell    : array<u32>, read-write. Per-atom flattened cell id.
//   @binding(3) cellCounts  : array<atomic<u32>>, read-write, length numCells.
//                             Pass 1 histogram. Must be zeroed before dispatch.
//   @binding(4) cellStart   : array<u32>, read-only, length numCells + 1
//                             (pass 2 only). Prefix-sum of counts with sentinel
//                             cellStart[numCells] = numAtoms, so the force
//                             kernel derives cell end as cellStart[c+1].
//   @binding(5) cellCursor  : array<atomic<u32>>, read-write, length numCells
//                             (pass 2 only). Must be zeroed before dispatch.
//   @binding(6) sortedIdx   : array<u32>, read-write, length numAtoms (pass 2 only).
//
// Flattened index: id = ix + nx * (iy + ny * iz), ix in [0, nx). Out-of-box
// atoms are clamped into range so no thread ever writes out of bounds.
//
// Workgroup: 64 threads (matches nonbonded_forces.wgsl and gpu.js).

// Non-bonded cutoff == cell edge (A). Matches the R_MAX/HCT cutoff used by
// src/physics/solvation/gb_obc2.js (cutoff 12 A) so Born-radii descreening
// and pair loops share one neighbor definition.
const NONBONDED_RCUT : f32 = 12.0;

struct CellParams {
  gridMin : vec3<f32>, // world-space origin of cell (0,0,0) (A)
  cellSize : f32,      // cell edge, must be >= NONBONDED_RCUT (A)
  gridDim : vec3<u32>, // cells per axis (nx, ny, nz)
  numAtoms : u32,      // N
  numCells : u32,      // nx*ny*nz
  pad : vec3<u32>,     // std140/uniform padding (keep struct 48 B)
};

@group(0) @binding(0) var<storage, read> pos : array<vec4<f32>>;
@group(0) @binding(1) var<uniform> params : CellParams;
@group(0) @binding(2) var<storage, read_write> atomCell : array<u32>;
@group(0) @binding(3) var<storage, read_write> cellCounts : array<atomic<u32>>;
@group(0) @binding(4) var<storage, read> cellStart : array<u32>;
@group(0) @binding(5) var<storage, read_write> cellCursor : array<atomic<u32>>;
@group(0) @binding(6) var<storage, read_write> sortedIdx : array<u32>;

// Clamped flattened cell id for a position.
fn cell_id(p : vec3<f32>) -> u32 {
  let rel = (p - params.gridMin) / params.cellSize;
  let nx = params.gridDim.x;
  let ny = params.gridDim.y;
  let nz = params.gridDim.z;
  let ix = clamp(u32(max(rel.x, 0.0)), 0u, nx - 1u);
  let iy = clamp(u32(max(rel.y, 0.0)), 0u, ny - 1u);
  let iz = clamp(u32(max(rel.z, 0.0)), 0u, nz - 1u);
  return ix + nx * (iy + ny * iz);
}

@compute @workgroup_size(64)
fn bin_count(@builtin(global_invocation_id) gid : vec3<u32>) {
  let i = gid.x;
  if (i >= params.numAtoms) { return; }
  let c = cell_id(pos[i].xyz);
  atomCell[i] = c;
  // c < numCells by construction (clamped), so the atomic is in bounds.
  atomicAdd(&cellCounts[c], 1u);
}

@compute @workgroup_size(64)
fn bin_fill(@builtin(global_invocation_id) gid : vec3<u32>) {
  let i = gid.x;
  if (i >= params.numAtoms) { return; }
  let c = atomCell[i];
  // c was clamped in pass 1; re-clamp defensively (host may reuse buffers).
  let cc = min(c, params.numCells - 1u);
  let s = atomicAdd(&cellCursor[cc], 1u);
  sortedIdx[cellStart[cc] + s] = i;
}
