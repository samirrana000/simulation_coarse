// nonbonded_forces.wgsl — tiled parallel non-bonded force evaluation.
//
// Physics (per pair i<j within RCUT, free boundary, kcal/mol, A, e):
//   LJ 12-6 with Lorentz-Berthelot combining rules, sigma/r clamped to 5
//     (G64 parity clamp with src/gpu.js), smooth-switched S(r) on [switchOn,
//     cutOff] with the exact heavy.js switchFunc/switchDeriv polynomial.
//   Screened Coulomb U_C = C*q_i*q_j/(epsIn*r), C = 332.06371.
//   GB-OBC reaction field (Still f_GB + Debye-Hueckel kappa screening) with
//     OBC effective Born radii streamed in params[i].w (computed on the CPU
//     via src/physics/solvation/gb_obc2.js computeBornRadii — HCT descreening
//     is gather-heavy and stays host-side). Analytic dU/dr mirrors
//     gbEnergyForces (dStillFdr/dgdf chain), fixed-radii explicit part.
//   Convention parity with the CPU kernels: LJ is switched, Coulomb+GB are
//     evaluated unswitched inside RCUT (same convention as
//     HeavyForceField._nonBondedGrid / src/physics/gb.js pairInteraction,
//     which hard-cuts GB at 12 A). Directional H-bonds, SASA, membrane and
//     bonded terms stay CPU-side by design (documented, not omitted).
//   1-2/1-3 exclusions are skipped via the sorted `excludedPairs` key list;
//     1-4 pairs in `scaledPairs` get s14 = 0.5 (heavy.js _scale14 parity).
//     Key packing is identical to the CPU: key = min*1000000u + max.
//
// Neighbor iteration: one thread per atom i; the 27 neighbor cells of i's
// own cell are walked via cellStart/sortedIdx from cell_list.wgsl. Atoms are
// stored cell-contiguously, so threads of a warp that share a cell stream
// the same sortedIdx run (coalesced); pos/params reads for i are consecutive
// across consecutive threads (coalesced); neighbor j reads are sequential
// within each cell run. Per-atom energy (LJ + Coulomb + GB pair + GB self)
// is returned in forces[i].w so the host reduces energies without a second
// dispatch.
//
// Binding layout (group 0):
//   @binding(0) pos           : array<vec4<f32>>, read. xyz (A), w unused.
//   @binding(1) params        : array<vec4<f32>>, read. (sigma, eps, q, bornR).
//   @binding(2) atomCell      : array<u32>, read. Flattened cell id per atom.
//   @binding(3) cellStart     : array<u32>, read, length numCells+1.
//                               Cell c holds sortedIdx[cellStart[c] .. cellStart[c+1]).
//   @binding(4) sortedIdx     : array<u32>, read, length numAtoms.
//   @binding(5) forces        : array<vec4<f32>>, read-write. xyz = force
//                               (kcal/mol/A), w = per-atom energy (kcal/mol).
//   @binding(6) uniforms      : NBParams (80 B).
//   @binding(7) excludedPairs : array<u32>, read, sorted skip keys, len exclCount.
//   @binding(8) scaledPairs   : array<u32>, read, sorted 1-4 keys, len scaleCount.
//
// Workgroup: 64 threads (one atom per thread).

const COULOMB_FLOOR_R2 : f32 = 1e-6; // r < 1e-3 A skipped (CPU: 1e-4 A^2 floor family)
const SR_CLAMP : f32 = 5.0;          // G64 sigma/r clamp (gpu.js parity)

struct NBParams {
  gridDim : vec3<u32>, // cells per axis (nx, ny, nz)
  numAtoms : u32,      // N
  gridMin : vec3<f32>, // cell-list origin (A)
  cellSize : f32,      // == RCUT (A)
  cutOff : f32,        // pair cutoff RCUT = 12.0 (A)
  switchOn : f32,      // LJ switch start (A)
  coulConst : f32,     // 332.06371 kcal*A/mol/e^2
  epsIn : f32,         // solute dielectric
  epsOut : f32,        // solvent dielectric
  kappa : f32,         // inverse Debye length (A^-1)
  useGB : u32,         // 1 = add Still-GB reaction field, 0 = Coulomb only
  exclCount : u32,     // len(excludedPairs)
  scaleCount : u32,    // len(scaledPairs)
  pad : vec3<u32>,     // 80 B total
};

@group(0) @binding(0) var<storage, read> pos : array<vec4<f32>>;
@group(0) @binding(1) var<storage, read> params : array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> atomCell : array<u32>;
@group(0) @binding(3) var<storage, read> cellStart : array<u32>;
@group(0) @binding(4) var<storage, read> sortedIdx : array<u32>;
@group(0) @binding(5) var<storage, read_write> forces : array<vec4<f32>>;
@group(0) @binding(6) var<uniform> uni : NBParams;
@group(0) @binding(7) var<storage, read> excludedPairs : array<u32>;
@group(0) @binding(8) var<storage, read> scaledPairs : array<u32>;

// Sorted-list membership (binary search). Lists are tiny in practice but
// binary search keeps the worst case logarithmic.
fn is_excluded(key : u32) -> bool {
  var lo = 0u;
  var hi = uni.exclCount;
  while (lo < hi) {
    let mid = (lo + hi) >> 1u;
    let v = excludedPairs[mid];
    if (v == key) { return true; }
    if (v < key) { lo = mid + 1u; } else { hi = mid; }
  }
  return false;
}

fn scale14(key : u32) -> f32 {
  var lo = 0u;
  var hi = uni.scaleCount;
  while (lo < hi) {
    let mid = (lo + hi) >> 1u;
    let v = scaledPairs[mid];
    if (v == key) { return 0.5; }
    if (v < key) { lo = mid + 1u; } else { hi = mid; }
  }
  return 1.0;
}

// Smooth switch S(r): 1 below on, 0 above off, heavy.js polynomial between.
fn switch_s(r : f32) -> f32 {
  let on = uni.switchOn;
  let off = uni.cutOff;
  if (r <= on) { return 1.0; }
  if (r >= off) { return 0.0; }
  let rsq = r * r;
  let on2 = on * on;
  let cut2 = off * off;
  let denom = (cut2 - on2) * (cut2 - on2) * (cut2 - on2);
  let num = (cut2 - rsq) * (cut2 - rsq) * (cut2 + 2.0 * rsq - 3.0 * on2);
  return num / denom;
}

fn switch_ds(r : f32) -> f32 {
  let on = uni.switchOn;
  let off = uni.cutOff;
  if (r <= on || r >= off) { return 0.0; }
  let rsq = r * r;
  let on2 = on * on;
  let cut2 = off * off;
  let denom = (cut2 - on2) * (cut2 - on2) * (cut2 - on2);
  let a = (cut2 - rsq) * (cut2 - rsq);
  let b = cut2 + 2.0 * rsq - 3.0 * on2;
  let da = -4.0 * r * (cut2 - rsq);
  let db = 4.0 * r;
  return (da * b + a * db) / denom;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid : vec3<u32>) {
  let i = gid.x;
  if (i >= uni.numAtoms) { return; }

  let pi = pos[i].xyz;
  let prm_i = params[i];
  let sig_i = prm_i.x;
  let eps_i = prm_i.y;
  let qi = prm_i.z;
  let Ri = max(prm_i.w, 0.5);

  let nx = uni.gridDim.x;
  let ny = uni.gridDim.y;
  let c = atomCell[i];
  let ix = c % nx;
  let iy = (c / nx) % ny;
  let iz = c / (nx * ny);

  let cut2 = uni.cutOff * uni.cutOff;
  var F = vec3<f32>(0.0, 0.0, 0.0);
  var eAcc : f32 = 0.0;

  for (var dz = -1; dz <= 1; dz++) {
    let nz = i32(iz) + dz;
    if (nz < 0 || nz >= i32(uni.gridDim.z)) { continue; }
    for (var dy = -1; dy <= 1; dy++) {
      let nyy = i32(iy) + dy;
      if (nyy < 0 || nyy >= i32(ny)) { continue; }
      for (var dx = -1; dx <= 1; dx++) {
        let nxx = i32(ix) + dx;
        if (nxx < 0 || nxx >= i32(nx)) { continue; }
        let nid = u32(nxx) + nx * (u32(nyy) + ny * u32(nz));
        let start = cellStart[nid];
        let endp = cellStart[nid + 1u];
        for (var k = start; k < endp; k++) {
          let j = sortedIdx[k];
          if (j == i) { continue; }
          let d = pos[j].xyz - pi;
          let r2 = dot(d, d);
          if (r2 >= cut2 || r2 < COULOMB_FLOOR_R2) { continue; }
          let r = sqrt(r2);

          // Exclusions / 1-4 scaling (CPU key parity: min*1e6+max).
          let key = select(j * 1000000u + i, i * 1000000u + j, i < j);
          if (is_excluded(key)) { continue; }
          let s14 = scale14(key);

          let prm_j = params[j];
          let qj = prm_j.z;

          // ---- Lennard-Jones (Lorentz-Berthelot, switched) ----
          let s = 0.5 * (sig_i + prm_j.x);
          let eps = sqrt(max(eps_i, 0.0) * max(prm_j.y, 0.0));
          let sr = min(s / r, SR_CLAMP);
          let sr2 = sr * sr;
          let sr6 = sr2 * sr2 * sr2;
          let ljE = 4.0 * eps * (sr6 * sr6 - sr6);
          let ljF = 4.0 * eps * (12.0 * sr6 * sr6 - 6.0 * sr6) / r; // -dU/dr
          let S = switch_s(r);
          let dS = switch_ds(r);
          let dudr = s14 * (S * (-ljF) + dS * ljE);
          eAcc += s14 * S * ljE;

          // ---- Screened Coulomb + Still GB reaction field ----
          if (qi != 0.0 && qj != 0.0) {
            let qq = qi * qj * s14;
            let uC = (uni.coulConst / uni.epsIn) * qq / r;
            var dudrE = -uC / r;
            eAcc += uC;
            if (uni.useGB == 1u) {
              let Rj = max(prm_j.w, 0.5);
              let a = max(Ri * Rj, 1e-6);
              let e = exp(-r2 / (4.0 * a));
              let f2 = r2 + a * e;
              let f = sqrt(max(f2, 1e-12));
              let eK = exp(-uni.kappa * f);
              let P = 1.0 / uni.epsIn - eK / uni.epsOut;
              let g = P / f;
              eAcc += -uni.coulConst * qq * g;
              let dfdr = (r * (1.0 - 0.25 * e)) / f;
              let dP = uni.kappa * eK / uni.epsOut;
              let dgdf = (dP * f - P) / f2;
              dudrE += -uni.coulConst * qq * dgdf * dfdr;
            }
            // F_i = dU/dr * dvec / r (dvec points i -> j).
            F += (dudrE / r) * d;
          }
          F += (dudr / r) * d;
        }
      }
    }
  }

  // ---- GB self (Born) term: per-atom, no pairs ----
  if (uni.useGB == 1u && qi != 0.0) {
    let R = Ri;
    let eK = exp(-uni.kappa * R);
    let P = 1.0 / uni.epsIn - eK / uni.epsOut;
    eAcc += -0.5 * uni.coulConst * qi * qi * P / R;
  }

  forces[i] = vec4<f32>(F, eAcc);
}
