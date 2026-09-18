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
 * Revolution 1 / Issue 3: L0 vs L1 wiring (kernels UNCHANGED — this note only
 * makes the existing opt-in explicit):
 *   L0 default (charges OFF, hbMode "off"): screened-Coulomb branch stays
 *     dead (q1==0), H-bond is the legacy isotropic bead-flag term (P/Cp/Cn
 *     beads only) — salt-bridge electrostatics + backbone directional H-bonds
 *     are missing by design (fast baseline, disclosed via
 *     ForceField.describePhysics() / bindingTermsActive()).
 *   L1 opt-in (charges ON, hbMode "directional"): CG_FORMAL_CHARGES revive the
 *     Coulomb branch; valid Cα-triplet virtual O-sites take the directional
 *     HB path (PASS 1b) with the isotropic term kept only as fallback for
 *     invalid sites. See forcefield.js:CG_PHYSICS_LEVELS/resolvePhysicsLevel.
 *
 * Two-pass structure:
 *   Pass 1 — pair terms (LJ / electrostatics / H-bond), plus EEF1-lite
 *            occupancy-density accumulation and (a, j, r) recording.
 *   Pass 2 — burial fraction/energy per ligand atom, then desolvation
 *            forces over the recorded pairs.
 */

/**
 * Revolution 1 / Issue 3: read-only wiring descriptor for the binding kernel.
 * Reports which cross terms CAN contribute given the ForceField flags — pure
 * introspection, zero energy effect (kernels below untouched).
 *
 * Revolution 2 / Issue 2: fallback-aware — in directional mode the kernel
 * keeps the legacy isotropic bead-flag term as fallback for invalid virtual
 * sites (see PASS 1 branch: `valid` site → directional, else isotropic when
 * `_protHB[i] && ligHB`). Reporting `isotropicHB:false` alone therefore lied
 * about the live fallback branch. `isotropicHBFallback` / `hasInvalidFallback`
 * is true when directional mode has at least one invalid/missing site, so the
 * descriptor is truthful. `isotropicHB` keeps its primary-only meaning for
 * backward compatibility; actual fallback energy is additionally gated by the
 * `_protHB` / `_ligHB` bead flags.
 * @param {object} ff ForceField instance (reads bindOn/nLigAtoms/chargesOn/hbMode/_vSites/nProt)
 * @returns {{lj:boolean, desolv:boolean, coulomb:boolean,
 *   directionalHB:boolean, isotropicHB:boolean,
 *   isotropicHBFallback:boolean, hasInvalidFallback:boolean}}
 */
export function bindingTermsActive(ff) {
  const hasPairs = !!(ff && ff.bindOn && ff.nLigAtoms > 0);
  const isDirectionalMode = hasPairs && ff.hbMode === "directional";
  // Revolution 3 / Issue 2: count valid vSites — directional is live only
  // with >=1 valid site (zero-valid ⇒ no dir); fallback additionally gated
  // on HB-capable beads on both sides (any(_protHB) && any(_ligHB)).
  let nValid = 0;
  let hasInvalid = false;
  if (isDirectionalMode) {
    const vs = ff._vSites;
    if (!vs) {
      hasInvalid = true;
    } else {
      const n = (ff && Number.isFinite(ff.nProt)) ? ff.nProt : vs.length;
      for (let i = 0; i < n; i++) {
        const s = vs[i];
        if (s && s.valid) nValid++;
        else hasInvalid = true;
      }
    }
  }
  const directional = isDirectionalMode && nValid >= 1;
  let hasInvalidFallback = false;
  if (hasInvalid) {
    let anyProtHB = false;
    const phb = ff ? ff._protHB : null;
    if (phb) { for (let i = 0; i < phb.length; i++) { if (phb[i]) { anyProtHB = true; break; } } }
    let anyLigHB = false;
    const lhb = ff ? ff._ligHB : null;
    if (lhb) { for (let i = 0; i < lhb.length; i++) { if (lhb[i]) { anyLigHB = true; break; } } }
    hasInvalidFallback = anyProtHB && anyLigHB;
  }
  return {
    lj: hasPairs,
    desolv: hasPairs,
    coulomb: hasPairs && ff.chargesOn === true,
    directionalHB: directional,
    // Legacy isotropic bead-flag term is the primary HB path unless directional
    // takes over; in directional mode it survives only as fallback for beads
    // whose virtual site is invalid (see PASS 1 branch below).
    isotropicHB: hasPairs && ff.hbMode !== "directional",
    isotropicHBFallback: hasInvalidFallback,
    hasInvalidFallback,
  };
}

/**
 * One-line human-readable wiring summary (for logs/regression output).
 * Revolution 2 / Issue 2: fallback-aware — reports `isoHB(fallback)` when the
 * primary isotropic term is off but the directional-mode invalid-site fallback
 * branch is live (see bindingTermsActive).
 * @param {object} ff ForceField instance
 * @returns {string} e.g. "L0: LJ+desolv+isoHB (coulomb OFF, directional OFF)"
 */
export function describeBindingTerms(ff) {
  const t = bindingTermsActive(ff);
  const lvl = (ff && ff.physicsLevel) || "L0";
  const on = [];
  if (t.lj) on.push("LJ");
  if (t.desolv) on.push("desolv");
  if (t.isotropicHB) on.push("isoHB");
  else if (t.isotropicHBFallback) on.push("isoHB(fallback)");
  if (t.coulomb) on.push("coulomb");
  if (t.directionalHB) on.push("dirHB");
  return `${lvl}: ${on.length ? on.join("+") : "no-binding"} `
    + `(coulomb ${t.coulomb ? "ON" : "OFF"}, directional ${t.directionalHB ? "ON" : "OFF"})`;
}
export function binding(ff, pos, f) {
  ff.bindingU = 0;
  ff.desolvU = 0;
  // Loop-2 S4 (R4 §5 item 1, R6 §5): per-term accumulators for the BindLog
  // energy channel. DEFAULT OFF — when trackTerms is false the branches
  // below are skipped entirely and the pass is bit-identical to pre-S4.
  const trk = ff.trackTerms === true;
  let tLJ = 0, tCoul = 0, tHB = 0;
  if (!ff.bindOn || ff.nLigAtoms === 0) {
    if (trk) {
      ff.bindLJU = 0; ff.bindCoulU = 0; ff.bindHBU = 0;
      ff.bindU = { lj: 0, coul: 0, hb: 0, desolv: 0 };
    }
    return 0;
  }

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
            const uLJ = 4 * e * phi * sw;
            U += uLJ;
            dUdr += 4 * e * (dphi * sw + phi * dsw);
            if (trk) tLJ += uLJ;

            // --- screened electrostatics ----------------------------------
            const q1 = ff._protQ[i];
            if (q1 !== 0 && ligQ !== 0) {
              const ch = Math.cosh(r / 8);
              const epsr = 4 + 76 * Math.tanh(r / 8);
              const epsrP = 9.5 / (ch * ch);            // dεr/dr = 9.5·sech²(r/8)
              const g = 1 / (epsr * r);                 // g = 1/(εr·r)
              const gp = -(epsr + r * epsrP) / (epsr * epsr * r * r); // dg/dr
              const A = ELC * q1 * ligQ;
              const uC = A * g * sw;
              U += uC;
              dUdr += A * (gp * sw + g * dsw);
              if (trk) tCoul += uC;
            }

            // --- H-bond -------------------------------------------------------
            if (ff.hbMode === "directional" && ff._vSites && ff._vSites[i] && ff._vSites[i].valid) {
              // Loop-2 S2: directional virtual-site term (R2 §2a) — all residues
              // accept via backbone O-sites (fixes 52% class-A gap), cos²/cos¹
              // gates by SS, valence-capped. Handled AFTER the pair loop over
              // pre-collected site pairs (below) — skip the isotropic bead term.
            } else if (ff._protHB[i] && ligHB) {
              const g = Math.exp(-((r - HB_R0) * (r - HB_R0)) / (2 * w2));
              const gp = -g * (r - HB_R0) / w2;         // dg/dr
              const B = -EPSHB;                          // U = −epsHB·g·sw
              const uH = B * g * sw;
              U += uH;
              dUdr += B * (gp * sw + g * dsw);
              if (trk) tHB += uH;
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
  // PASS 1b — directional H-bonds over virtual O-sites (Loop-2 S2)
  // =============================================================
  if (ff.hbMode === "directional" && ff._vSites) {
    ff._vSiteValence.fill(0);
    for (let a = 0; a < ff.nLigAtoms; a++) {
      if (!ff._ligHB[a]) continue;                    // ligand must be H-bond-capable
      const la = nProt + a, lx = 3 * la;
      for (let i = 0; i < nProt; i++) {
        const site = ff._vSites[i];
        if (!site || !site.valid) continue;
        if (ff._vSiteValence[i] >= 2) continue;      // acceptor valence cap
        // cheap distance gate on Cα first (sites sit ~1 Å off Cα)
        const cx = 3 * i;
        const dxc = pos[cx] - pos[lx], dyc = pos[cx + 1] - pos[lx + 1], dzc = pos[cx + 2] - pos[lx + 2];
        if (dxc * dxc + dyc * dyc + dzc * dzc > 36) continue; // 6 Å Cα gate
        // rebuild site position from the live triplet (sites must track motion)
        const im = Math.max(0, i - 1), ip = Math.min(nProt - 1, i + 1);
        if (im === i || ip === i) continue;
        const li = [pos[3 * im], pos[3 * im + 1], pos[3 * im + 2]];
        const ci = [pos[cx], pos[cx + 1], pos[cx + 2]];
        const ni = [pos[3 * ip], pos[3 * ip + 1], pos[3 * ip + 2]];

        let u = [ni[0] - li[0], ni[1] - li[1], ni[2] - li[2]];
        let ul = Math.hypot(...u) || 1; u = [u[0] / ul, u[1] / ul, u[2] / ul];
        let v = [2 * ci[0] - ni[0] - li[0], 2 * ci[1] - ni[1] - li[1], 2 * ci[2] - ni[2] - li[2]];
        let vl = Math.hypot(...v) || 1; v = [v[0] / vl, v[1] / vl, v[2] / vl];
        let w2v = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        let wl = Math.hypot(...w2v) || 1; w2v = [w2v[0] / wl, w2v[1] / wl, w2v[2] / wl];
        const ss = site.ss; // class from reference (cheap; motion doesn't change SS fast)
        const coef = ss === 0 ? [1.98, -0.39] : [1.75, -0.64]; // fitted (S2)
        const oX = ci[0] + coef[0] * u[0] + coef[1] * v[0];
        const oY = ci[1] + coef[0] * u[1] + coef[1] * v[1];
        const oZ = ci[2] + coef[0] * u[2] + coef[1] * v[2];
        // cone axis = normalize(oPos − Cα)
        const ca = [oX - ci[0], oY - ci[1], oZ - ci[2]];
        const cal = Math.hypot(...ca) || 1;
        const caN = [ca[0] / cal, ca[1] / cal, ca[2] / cal];
        // energy
        const dx = pos[lx] - oX, dy = pos[lx + 1] - oY, dz = pos[lx + 2] - oZ;
        const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > 16 || r2 < 1e-10) continue;
        const r = Math.sqrt(r2);
        const dr = r - 3.0;
        const g = Math.exp(-(dr * dr) / (2 * 0.5 * 0.5));
        if (g < 1e-4) continue;
        let sw = 1, dsw = 0;
        if (r > 3.5) { const t = (r - 3.5) / 0.5, t2 = t * t; sw = 1 - t2 * t * (10 - 15 * t + 6 * t2); dsw = -(30 / 0.5) * t2 * (1 - t) * (1 - t); }
        const inv = 1 / r;
        const cosT = ((dx * inv) * caN[0] + (dy * inv) * caN[1] + (dz * inv) * caN[2]); // donor beyond O along cone
        const gatePow = ss === 2 ? 1 : 2;
        const m0 = Math.max(0, cosT);
        if (m0 <= 0) continue;
        const gate = Math.pow(m0, gatePow);
        if (ff._vSiteValence[i] >= 2) continue;
        ff._vSiteValence[i]++;
          const epsHB = 2.0;
          const Uhb = -epsHB * g * sw * gate;
          U += Uhb;
          if (trk) tHB += Uhb;
        // forces (radial + angular on ligand; thirds to triplet)
        const dg = -(dr / 0.25) * g;
        const dU_dr = -epsHB * (dg * sw + g * dsw) * gate;
        const dgate = gatePow * Math.pow(m0, gatePow - 1);
        const pre = -epsHB * g * sw * dgate;
        const cxs = dx * inv, cys = dy * inv, czs = dz * inv;
        const cT_dx = -caN[0] / r - cxs * cosT / r;
        const cT_dy = -caN[1] / r - cys * cosT / r;
        const cT_dz = -caN[2] / r - czs * cosT / r;
        const Fx = -(dU_dr * cxs + pre * cT_dx);
        const Fy = -(dU_dr * cys + pre * cT_dy);
        const Fz = -(dU_dr * czs + pre * cT_dz);
        f[lx] -= Fx; f[lx + 1] -= Fy; f[lx + 2] -= Fz;
        const w3 = 1 / 3;
        f[3 * im] += Fx * w3; f[3 * im + 1] += Fy * w3; f[3 * im + 2] += Fz * w3;
        f[cx] += Fx * w3; f[cx + 1] += Fy * w3; f[cx + 2] += Fz * w3;
        f[3 * ip] += Fx * w3; f[3 * ip + 1] += Fy * w3; f[3 * ip + 2] += Fz * w3;
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
  if (trk) {
    ff.bindLJU = tLJ; ff.bindCoulU = tCoul; ff.bindHBU = tHB;
    ff.bindU = { lj: tLJ, coul: tCoul, hb: tHB, desolv: Udesolv };
  }
  return U + Udesolv;
}
