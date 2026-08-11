/**
 * ff-binding.js — protein–ligand binding kernel (cross LJ + screened
 * electrostatics + H-bond + EEF1-lite desolvation) for the ForceField class
 * in forcefield.js (item 5 modularization).
 *
 * Verbatim body of the original ForceField._binding with `this` replaced by
 * an explicit `ff` parameter; the spatial-hash helpers (_cellKey,
 * _encodeCell, _decodeX/Y/Z, _pairKey) stay on the ForceField instance and
 * are invoked through ff. forcefield.js wraps this with the same
 * _binding(pos, f) signature, so no call site changes. The full physics
 * documentation lives on the class-method wrapper in forcefield.js.
 *
 * Two-pass structure:
 *   Pass 1 — pair terms (LJ / electrostatics / H-bond), plus EEF1-lite
 *            occupancy-density accumulation and (a, j, r) recording.
 *   Pass 2 — burial fraction/energy per ligand atom, then desolvation
 *            forces over the recorded pairs.
 */
export function binding(ff, pos, f) {
  ff.bindingU = 0;
  ff.desolvU = 0;
  if (!ff.bindOn || ff.nLigAtoms === 0) return 0;

  const nProt = ff.nProt;
  const rc = ff.bindRcut, rsw = 0.85 * rc;
  const rc2 = rc * rc, invDelta = 1 / (rc - rsw);
  const cell = rc;                  // cell = cutoff ⇒ 27-cell scan is complete
  const grid = ff._gridB;

  // EEF1-lite constants: contact distance r0, Gaussian width σ, and the
  // density→burial scale nScale (n_a ≈ nScale ≈ 3 ⇒ ~63% buried).
  const R0 = 4.5, SIG = 1.8, NS = 3.0;
  const inv2sig2 = 1 / (2 * SIG * SIG); // 1/(2σ²) for g(r)
  const invSig2 = 1 / (SIG * SIG);      // 1/σ²   for dg/dr

  // Per-atom occupancy density + burial-derivative scratch (GC-free)
  const dens = ff._dens, dBdn = ff._dBdn;
  dens.fill(0);

  // Reused flat pair records for Pass 2 (a, j, r); capped at nProt·nLigAtoms.
  let bpN = 0;
  const bpA = ff._bpA, bpJ = ff._bpJ, bpR = ff._bpR;

  // rebuild the protein-only grid in place (GC-free, like _repulsion)
  for (const arr of grid.values()) arr.length = 0;
  for (let i = 0; i < nProt; i++) {
    const key = ff._cellKey(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2], cell);
    let arr = grid.get(key);
    if (!arr) grid.set(key, (arr = []));
    arr.push(i);
  }

  // =============================================================
  // PASS 1 — pair terms (LJ / electrostatics / H-bond) + density
  // =============================================================
  let U = 0;
  const EPSHB = 0.8, HB_R0 = 3.2, HB_W = 0.6, ELC = 332.0637;
  const w2 = HB_W * HB_W;

  for (let a = 0; a < ff.nLigAtoms; a++) {
    const la = nProt + a, lx = 3 * la;
    const laX = pos[lx], laY = pos[lx + 1], laZ = pos[lx + 2];
    const ligSig = ff._ligSigma[a], ligEps = ff._ligEps[a];
    const ligQ = ff._ligQ[a], ligHB = ff._ligHB[a];
    const ckey = ff._cellKey(laX, laY, laZ, cell);
    const cx = ff._decodeX(ckey), cy = ff._decodeY(ckey), cz = ff._decodeZ(ckey);

    // 27 neighbouring cells (directed protein→ligand, no double counting)
    for (let ox = -1; ox <= 1; ox++)
      for (let oy = -1; oy <= 1; oy++)
        for (let oz = -1; oz <= 1; oz++) {
          const arr = grid.get(ff._encodeCell(cx + ox, cy + oy, cz + oz));
          if (!arr) continue;
          for (let bi = 0; bi < arr.length; bi++) {
            const i = arr[bi], ix = 3 * i;
            const dx = pos[ix] - laX, dy = pos[ix + 1] - laY, dz = pos[ix + 2] - laZ;
            const r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= rc2 || r2 < 1e-10) continue;
            const r = Math.sqrt(r2);

            // Holo contact springs already govern this native pair — skip it
            // so the binding pair pass does not double-count the contact.
            const pk = ff._pairKey(i, la);
            if (ff._excluded.has(pk)) continue;

            // --- EEF1-lite occupancy density + pair record for Pass 2 ------
            // g(r) feeds the soft count dens[a]; the (a, j, r) triple lets
            // Pass 2 add the desolvation force without re-scanning the grid.
            const gr = Math.exp(-((r - R0) * (r - R0)) * inv2sig2);
            dens[a] += gr;
            bpA[bpN] = a; bpJ[bpN] = i; bpR[bpN] = r; bpN++;

            // smooth switch + its r-derivative
            let sw, dsw;
            if (r <= rsw) { sw = 1; dsw = 0; }
            else {
              const t = (r - rsw) * invDelta, t2 = t * t;
              sw = 1 - t2 * t * (10 - 15 * t + 6 * t2);
              dsw = -(30 * invDelta) * t2 * (1 - t) * (1 - t);
            }

            let dUdr = 0;   // total dU/dr of this pair (all cross terms)

            // --- cross LJ 12-6 (attractive well) ---------------------------
            const s = 0.5 * (ff._protSigma[i] + ligSig);
            const e = Math.sqrt(ff._protEps[i] * ligEps);
            const s6 = Math.pow(s / r, 6), s12 = s6 * s6;
            const phi = s12 - s6;                       // φ = (σ/r)¹² − (σ/r)⁶
            const dphi = (6 * s6 - 12 * s12) / r;       // dφ/dr
            U += 4 * e * phi * sw;
            dUdr += 4 * e * (dphi * sw + phi * dsw);

            // --- screened electrostatics ----------------------------------
            const q1 = ff._protQ[i];
            if (q1 !== 0 && ligQ !== 0) {
              const ch = Math.cosh(r / 8);
              const epsr = 4 + 76 * Math.tanh(r / 8);
              const epsrP = 9.5 / (ch * ch);            // dεr/dr = 9.5·sech²(r/8)
              const g = 1 / (epsr * r);                 // g = 1/(εr·r)
              const gp = -(epsr + r * epsrP) / (epsr * epsr * r * r); // dg/dr
              const A = ELC * q1 * ligQ;
              U += A * g * sw;
              dUdr += A * (gp * sw + g * dsw);
            }

            // --- H-bond (donor/acceptor flags both set) -------------------
            if (ff._protHB[i] && ligHB) {
              const g = Math.exp(-((r - HB_R0) * (r - HB_R0)) / (2 * w2));
              const gp = -g * (r - HB_R0) / w2;         // dg/dr
              const B = -EPSHB;                          // U = −epsHB·g·sw
              U += B * g * sw;
              dUdr += B * (gp * sw + g * dsw);
            }

            // force on protein i = −(dU/dr)·(dx/r); on the ligand, opposite
            const F = -dUdr / r;
            const fX = F * dx, fY = F * dy, fZ = F * dz;
            f[ix] += fX; f[ix + 1] += fY; f[ix + 2] += fZ;
            f[lx] -= fX; f[lx + 1] -= fY; f[lx + 2] -= fZ;
          }
        }
  }

  // =============================================================
  // PASS 2 — EEF1-lite burial: burial fraction, energy, forces
  // =============================================================
  // B_a = 1 − exp(−n_a/3); U_desolv = Σ_a ΔG_a·B_a. dB/dn is cached per atom
  // so every pair of atom a reuses the same factor (B depends on total n_a).
  let Udesolv = 0;
  for (let a = 0; a < ff.nLigAtoms; a++) {
    const n = dens[a];
    const e = Math.exp(-n / NS);           // exp(−n_a/3)
    dBdn[a] = e / NS;                      // dB/dn = exp(−n_a/3)/3
    Udesolv += ff._ligdG[a] * (1 - e);   // ΔG_a·B_a
  }

  // Chain rule dU_desolv/dr = ΔG_a·(dB/dn)·(dg/dr); force on protein j is
  // −(dU/dr)·(dx/r) with dx = pos_j − pos_la — same convention as Pass 1.
  for (let k = 0; k < bpN; k++) {
    const a = bpA[k], j = bpJ[k], r = bpR[k];
    const g = Math.exp(-((r - R0) * (r - R0)) * inv2sig2);
    const dgdr = -g * (r - R0) * invSig2;              // dg/dr
    const dUdr = ff._ligdG[a] * dBdn[a] * dgdr;      // dU_desolv/dr
    const Fs = -dUdr / r;                              // −(dU/dr)/r, vector factor
    const la = nProt + a, lx = 3 * la, jx = 3 * j;
    const dx = pos[jx] - pos[lx], dy = pos[jx + 1] - pos[lx + 1], dz = pos[jx + 2] - pos[lx + 2];
    f[jx] += Fs * dx; f[jx + 1] += Fs * dy; f[jx + 2] += Fs * dz;
    f[lx] -= Fs * dx; f[lx + 1] -= Fs * dy; f[lx + 2] -= Fs * dz;
  }

  ff.bindingU = U + Udesolv;
  ff.desolvU = Udesolv;
  return U + Udesolv;
}
