/**
 * unbinding_smd.js — steered MD (constant-velocity / constant-force pulling)
 * along a ligand exit channel, with Jarzynski free-energy reconstruction and
 * a rupture-based koff ranking surrogate.
 *
 * Pulling coordinate: the ligand-COM projection onto the exit unit vector,
 *   s(t) = (c(t) − c₀) · d̂,   c = ligand COM, c₀ = bound reference COM.
 * Exit vector default: d̂ = normalize(c₀ − pocketCOM) pointing to solvent.
 *
 * Constant-velocity (harmonic guide, stiff-spring limit):
 *   U(λ) = ½k (s − λ)²,   λ(t) = v·t,
 *   per-atom force  f_a = −k(s−λ)·d̂/nLig   (∂s/∂x_a = d̂/nLig: COM-only pull,
 *                                            no torque — cf. funnel.js),
 *   work increment dW = ∂U/∂λ·dλ = −k(s−λ)·v·dt = k(λ−s)·v·dt.
 * Constant-force:
 *   per-atom force  f_a = f0·d̂/nLig,   dW = f0·ds.
 *
 * Jarzynski equality (Jarzynski, PRL 78, 2690 (1997)):
 *   exp(−βΔF) = ⟨exp(−βW)⟩,   β = 1/(kB·T),
 * evaluated with a max-shift for floating-point safety:
 *   ΔF = Wmin − (1/β)·ln⟨exp(−β(W−Wmin))⟩.
 *
 * koff surrogate: rupture-decomposition ranking score
 *   score_i = exp(+β·(Wrupt,i − ΔF)) ≥ ~1,
 * where Wrupt,i is the nonequilibrium work accumulated up to first passage of
 * s ≥ sRupture. Near 1 ⇒ near-reversible pull; larger ⇒ more dissipative
 * escape. This is an ORDERING proxy for koff across ligands/protocols — not
 * an absolute rate; use docs/NETWORK.md (MSM) for publishable kinetics.
 *
 * Dynamics: a self-contained BAOAB Langevin loop (same SDE as
 * src/integrator.js: m·dv = (−∇U − mζv)dt + √(2mζkBT)·dW) with a seeded LCG so
 * ensembles are reproducible. Units Å / ps / kcal/mol / Da.
 *
 * DOM-free; runs under plain Node.
 */

import { KB_KCAL, KCONV } from "../units.js?v=10";

/* ------------------------------------------------------------------ */
/*  Small utilities                                                    */
/* ------------------------------------------------------------------ */

/**
 * Ligand centre of mass.
 * @param {Float64Array} pos  length 3n
 * @param {number} ligStart   first ligand global index
 * @param {number} nLig       ligand atom count
 * @returns {[number, number, number]}
 */
export function ligandCOM(pos, ligStart, nLig) {
  let x = 0, y = 0, z = 0;
  for (let a = 0; a < nLig; a++) {
    x += pos[3 * (ligStart + a)];
    y += pos[3 * (ligStart + a) + 1];
    z += pos[3 * (ligStart + a) + 2];
  }
  return [x / nLig, y / nLig, z / nLig];
}

/**
 * Exit channel: pocket COM (protein beads within rPocket of the native ligand
 * COM) and the solvent-pointing unit vector.
 * @param {object} ff  force field with .ref, .nProt, .n, .ligandStart?
 * @param {object} [opts]
 * @param {number} [opts.rPocket=8.0]
 * @param {number[]} [opts.dir]  explicit exit vector (normalized internally)
 * @returns {{dir: number[], pocketCOM: number[], ligCOM0: number[], nLig: number, ligStart: number}}
 */
export function exitVector(ff, opts = {}) {
  const nProt = ff.nProt;
  const ligStart = ff.ligandStart ?? nProt;
  const nLig = ff.n - ligStart;
  if (nLig <= 0) throw new Error("exitVector: system has no ligand atoms.");
  const ligCOM0 = ligandCOM(ff.ref, ligStart, nLig);
  let dir;
  if (opts.dir && opts.dir.length === 3) {
    const l = Math.hypot(opts.dir[0], opts.dir[1], opts.dir[2]) || 1;
    dir = [opts.dir[0] / l, opts.dir[1] / l, opts.dir[2] / l];
  } else {
    const rPocket = opts.rPocket ?? 8.0;
    let px = 0, py = 0, pz = 0, c = 0;
    for (let i = 0; i < nProt; i++) {
      const dx = ff.ref[3 * i] - ligCOM0[0];
      const dy = ff.ref[3 * i + 1] - ligCOM0[1];
      const dz = ff.ref[3 * i + 2] - ligCOM0[2];
      if (Math.sqrt(dx * dx + dy * dy + dz * dz) <= rPocket) {
        px += ff.ref[3 * i]; py += ff.ref[3 * i + 1]; pz += ff.ref[3 * i + 2]; c++;
      }
    }
    if (c === 0) { // degenerate: nearest 6 beads
      const order = [];
      for (let i = 0; i < nProt; i++) {
        const dx = ff.ref[3 * i] - ligCOM0[0];
        const dy = ff.ref[3 * i + 1] - ligCOM0[1];
        const dz = ff.ref[3 * i + 2] - ligCOM0[2];
        order.push([dx * dx + dy * dy + dz * dz, i]);
      }
      order.sort((a, b) => a[0] - b[0]);
      for (let k = 0; k < Math.min(6, nProt); k++) {
        const i = order[k][1];
        px += ff.ref[3 * i]; py += ff.ref[3 * i + 1]; pz += ff.ref[3 * i + 2]; c++;
      }
    }
    const pocketCOM = [px / c, py / c, pz / c];
    const ex = ligCOM0[0] - pocketCOM[0];
    const ey = ligCOM0[1] - pocketCOM[1];
    const ez = ligCOM0[2] - pocketCOM[2];
    const l = Math.hypot(ex, ey, ez) || 1;
    dir = [ex / l, ey / l, ez / l];
    return { dir, pocketCOM, ligCOM0, nLig, ligStart };
  }
  // Explicit dir: pocket COM still reported (native ligand COM proxy).
  return { dir, pocketCOM: [...ligCOM0], ligCOM0, nLig, ligStart };
}

/** Seeded LCG (Numerical Recipes) → uniform [0,1). */
function makeRng(seed) {
  let s = (seed >>> 0) || 1;
  return () => {
    s = (1664525 * s + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** Standard normal via Box–Muller on a uniform rng. */
function makeGaussian(rng) {
  let spare = null;
  return () => {
    if (spare !== null) { const v = spare; spare = null; return v; }
    let u = 0, v = 0;
    while (u === 0) u = rng();
    v = rng();
    const r = Math.sqrt(-2 * Math.log(u));
    spare = r * Math.sin(2 * Math.PI * v);
    return r * Math.cos(2 * Math.PI * v);
  };
}

/* ------------------------------------------------------------------ */
/*  Single pulls                                                       */
/* ------------------------------------------------------------------ */

/**
 * @typedef {object} PullResult
 * @property {number} work           total nonequilibrium work (kcal/mol)
 * @property {number} workToRupture  work at first passage s ≥ sRupture (kcal/mol)
 * @property {number} ruptureForce   max |guide force| before rupture (kcal/mol/Å)
 * @property {number} ruptureStep    step index of first passage (−1 if never)
 * @property {number} sFinal         final projection s (Å)
 * @property {object} trace          sampled {t, s, f, w} arrays
 */

/**
 * Constant-velocity SMD pull with a self-contained BAOAB Langevin loop.
 *
 * @param {object} o
 * @param {object} o.ff    force field (.compute(pos), .forces, .masses, .n)
 * @param {Float64Array} o.pos0  bound starting positions
 * @param {number[]} o.dir       exit unit vector
 * @param {number[]} o.ligCOM0   bound ligand COM
 * @param {number} [o.k=2.0]     guide spring (kcal/mol/Å²)
 * @param {number} [o.v=1.0]     pulling speed (Å/ps)
 * @param {number} [o.nSteps=2000]
 * @param {number} [o.dt=0.001]  ps
 * @param {number} [o.T=300]     K
 * @param {number} [o.zeta=8.0]  ps⁻¹
 * @param {number} [o.seed=1]
 * @param {number} [o.sRupture=8.0]  Å first-passage threshold
 * @param {number} [o.recordEvery=10]
 * @returns {PullResult}
 */
export function runConstantVelocityPull(o) {
  const { ff, pos0, dir, ligCOM0 } = o;
  const k = o.k ?? 2.0, v = o.v ?? 1.0;
  const nSteps = o.nSteps ?? 2000, dt = o.dt ?? 0.001;
  const T = o.T ?? 300, zeta = o.zeta ?? 8.0;
  const sRupture = o.sRupture ?? 8.0, recordEvery = o.recordEvery ?? 10;
  const n = ff.n;
  const ligStart = ff.ligandStart ?? ff.nProt;
  const nLig = n - ligStart;
  if (nLig <= 0) throw new Error("runConstantVelocityPull: no ligand atoms.");
  const rng = makeRng(o.seed ?? 1);
  const gauss = makeGaussian(rng);
  const mass = new Float64Array(3 * n);
  for (let i = 0; i < 3 * n; i++) mass[i] = (ff.masses && ff.masses[(i / 3) | 0]) || 110;
  const thermal = new Float64Array(3 * n);
  for (let i = 0; i < 3 * n; i++) thermal[i] = Math.sqrt(KB_KCAL * T * KCONV / mass[i]);
  const invM = new Float64Array(3 * n);
  for (let i = 0; i < 3 * n; i++) invM[i] = 1 / mass[i];

  const pos = Float64Array.from(pos0);
  const vel = new Float64Array(3 * n);
  for (let i = 0; i < 3 * n; i++) vel[i] = thermal[i] * gauss();
  const proj = (p) => {
    const c = ligandCOM(p, ligStart, nLig);
    return (c[0] - ligCOM0[0]) * dir[0] + (c[1] - ligCOM0[1]) * dir[1] + (c[2] - ligCOM0[2]) * dir[2];
  };
  const c1 = Math.exp(-zeta * dt), c2 = Math.sqrt(Math.max(0, 1 - c1 * c1));
  const h = dt / 2;
  const kick = (F, hh) => {
    for (let i = 0; i < 3 * n; i++) vel[i] += hh * KCONV * F[i] * invM[i];
  };
  // Initial force + SMD guide.
  let U = ff.compute(pos);
  if (!Number.isFinite(U)) throw new Error("runConstantVelocityPull: non-finite initial energy.");
  let s = proj(pos);
  let lam = 0;
  const applyGuide = (F, ss, ll) => {
    const g = -k * (ss - ll); // total COM force along d̂ (kcal/mol/Å)
    const per = g / nLig;
    for (let a = 0; a < nLig; a++) {
      const c = 3 * (ligStart + a);
      F[c] += per * dir[0]; F[c + 1] += per * dir[1]; F[c + 2] += per * dir[2];
    }
    return g;
  };
  let F = ff.forces;
  let g = applyGuide(F, s, lam);
  let W = 0;
  let maxF = Math.abs(g);
  let ruptured = s >= sRupture;
  let ruptureStep = ruptured ? 0 : -1;
  let workToRupture = ruptured ? 0 : NaN;
  const tr = { t: [], s: [], f: [], w: [] };
  const push = (t) => tr.t.push(t) || tr.s.push(s) || tr.f.push(g) || tr.w.push(W);
  push(0);
  for (let step = 1; step <= nSteps; step++) {
    kick(F, h);
    for (let i = 0; i < 3 * n; i++) pos[i] += h * vel[i];
    for (let i = 0; i < 3 * n; i++) vel[i] = c1 * vel[i] + c2 * thermal[i] * gauss();
    for (let i = 0; i < 3 * n; i++) pos[i] += h * vel[i];
    ff.compute(pos);
    F = ff.forces;
    s = proj(pos);
    lam = v * step * dt;
    g = applyGuide(F, s, lam);
    W += k * (lam - s) * v * dt; // dW = ∂U/∂λ·dλ
    if (!ruptured) {
      if (Math.abs(g) > maxF) maxF = Math.abs(g);
      if (s >= sRupture) { ruptured = true; ruptureStep = step; workToRupture = W; }
    }
    kick(F, h);
    if (step % recordEvery === 0 || step === nSteps) push(step * dt);
  }
  if (!ruptured) workToRupture = W;
  return {
    work: W, workToRupture, ruptureForce: maxF, ruptureStep, sFinal: s,
    trace: { t: Float64Array.from(tr.t), s: Float64Array.from(tr.s), f: Float64Array.from(tr.f), w: Float64Array.from(tr.w) },
  };
}

/**
 * Constant-force SMD pull (same loop, guide U = −f0·s).
 * @param {object} o  as above with [o.f0=5.0] kcal/mol/Å instead of k/v
 * @returns {PullResult} (ruptureForce ≡ f0)
 */
export function runConstantForcePull(o) {
  const { ff, pos0, dir, ligCOM0 } = o;
  const f0 = o.f0 ?? 5.0;
  const nSteps = o.nSteps ?? 2000, dt = o.dt ?? 0.001;
  const T = o.T ?? 300, zeta = o.zeta ?? 8.0;
  const sRupture = o.sRupture ?? 8.0, recordEvery = o.recordEvery ?? 10;
  const n = ff.n;
  const ligStart = ff.ligandStart ?? ff.nProt;
  const nLig = n - ligStart;
  if (nLig <= 0) throw new Error("runConstantForcePull: no ligand atoms.");
  const rng = makeRng(o.seed ?? 1);
  const gauss = makeGaussian(rng);
  const mass = new Float64Array(3 * n);
  for (let i = 0; i < 3 * n; i++) mass[i] = (ff.masses && ff.masses[(i / 3) | 0]) || 110;
  const thermal = new Float64Array(3 * n);
  for (let i = 0; i < 3 * n; i++) thermal[i] = Math.sqrt(KB_KCAL * T * KCONV / mass[i]);
  const invM = new Float64Array(3 * n);
  for (let i = 0; i < 3 * n; i++) invM[i] = 1 / mass[i];
  const pos = Float64Array.from(pos0);
  const vel = new Float64Array(3 * n);
  for (let i = 0; i < 3 * n; i++) vel[i] = thermal[i] * gauss();
  const proj = (p) => {
    const c = ligandCOM(p, ligStart, nLig);
    return (c[0] - ligCOM0[0]) * dir[0] + (c[1] - ligCOM0[1]) * dir[1] + (c[2] - ligCOM0[2]) * dir[2];
  };
  const c1 = Math.exp(-zeta * dt), c2 = Math.sqrt(Math.max(0, 1 - c1 * c1));
  const h = dt / 2;
  const kick = (F, hh) => {
    for (let i = 0; i < 3 * n; i++) vel[i] += hh * KCONV * F[i] * invM[i];
  };
  ff.compute(pos);
  let s = proj(pos);
  let sPrev = s;
  let W = 0;
  let ruptureStep = s >= sRupture ? 0 : -1;
  let workToRupture = ruptureStep === 0 ? 0 : NaN;
  const tr = { t: [], s: [], f: [], w: [] };
  const applyForce = (F) => {
    const per = f0 / nLig;
    for (let a = 0; a < nLig; a++) {
      const c = 3 * (ligStart + a);
      F[c] += per * dir[0]; F[c + 1] += per * dir[1]; F[c + 2] += per * dir[2];
    }
  };
  let F = ff.forces;
  applyForce(F);
  tr.t.push(0); tr.s.push(s); tr.f.push(f0); tr.w.push(0);
  for (let step = 1; step <= nSteps; step++) {
    kick(F, h);
    for (let i = 0; i < 3 * n; i++) pos[i] += h * vel[i];
    for (let i = 0; i < 3 * n; i++) vel[i] = c1 * vel[i] + c2 * thermal[i] * gauss();
    for (let i = 0; i < 3 * n; i++) pos[i] += h * vel[i];
    ff.compute(pos);
    F = ff.forces;
    sPrev = s;
    s = proj(pos);
    W += f0 * (s - sPrev); // dW = f0·ds
    applyForce(F);
    if (ruptureStep < 0 && s >= sRupture) { ruptureStep = step; workToRupture = W; }
    kick(F, h);
    if (step % recordEvery === 0 || step === nSteps) {
      tr.t.push(step * dt); tr.s.push(s); tr.f.push(f0); tr.w.push(W);
    }
  }
  if (ruptureStep < 0) workToRupture = W;
  return {
    work: W, workToRupture, ruptureForce: f0, ruptureStep, sFinal: s,
    trace: { t: Float64Array.from(tr.t), s: Float64Array.from(tr.s), f: Float64Array.from(tr.f), w: Float64Array.from(tr.w) },
  };
}

/* ------------------------------------------------------------------ */
/*  Ensemble + Jarzynski + koff surrogate                              */
/* ------------------------------------------------------------------ */

/** Switch off bound-state restraints; returns saved state for restore.
 * NOTE: CG holo springs act in compute() unconditionally once built
 * (holoOn only gates rebuilds), so the spring ARRAY is swapped out — toggling
 * the flag alone would not remove the force. */
function dropRestraints(ff) {
  const saved = {};
  if (ff && ff.holoSprings && ff.holoSprings.length) {
    saved.holoSprings = ff.holoSprings;
    ff.holoSprings = new Float64Array(0);
  }
  if (ff && typeof ff.funnelOn === "boolean") { saved.funnelOn = ff.funnelOn; ff.funnelOn = false; }
  return saved;
}

/** Restore state saved by dropRestraints(). */
function restoreRestraints(ff, saved) {
  if (!ff || !saved) return;
  if ("holoSprings" in saved) ff.holoSprings = saved.holoSprings;
  if ("funnelOn" in saved) ff.funnelOn = saved.funnelOn;
}

/**
 * Run an ensemble of independent pulls (seeds seed..seed+nPulls−1).
 *
 * Bias hygiene: holo native-pose springs (CG ForceField.holoOn) and the funnel
 * metadynamics bias (ff.funnelOn) encode/restrain the BOUND state — pulling
 * against them measures restraint-breaking, not unbinding. Unless
 * opts.keepRestraints is true they are switched OFF for the duration of each
 * pull and restored afterwards (try/finally), so the work reflects the
 * unbiased unbinding path.
 *
 * @param {object} o
 * @param {object} [o.ff]         force field (with .ref bound reference)
 * @param {Float64Array} [o.pos0] bound start (default ff.ref)
 * @param {Function} [o.makeSystem] alternative factory (i) → {ff, pos0}
 * @param {number} [o.nPulls=4]
 * @param {"velocity"|"force"} [o.mode="velocity"]
 * @param {number[]} [o.dir]      exit vector (computed once if omitted)
 * @param {number} [o.seed=42]
 * @param {boolean} [o.keepRestraints=false]  leave holo/funnel bias on
 * @param {object} [o.pull]       forwarded to the single-pull runner
 * @returns {{mode: string, pulls: PullResult[], works: Float64Array,
 *   meanWork: number, dir: number[], ligCOM0: number[]}}
 */
export function runPullingEnsemble(o = {}) {
  const nPulls = o.nPulls ?? 4;
  if (!Number.isInteger(nPulls) || nPulls < 1) {
    throw new Error("runPullingEnsemble: nPulls ≥ 1 is required.");
  }
  const mode = o.mode ?? "velocity";
  if (mode !== "velocity" && mode !== "force") {
    throw new Error(`runPullingEnsemble: unknown mode "${mode}".`);
  }
  const first = o.makeSystem ? o.makeSystem(0) : { ff: o.ff, pos0: o.pos0 ?? o.ff?.ref };
  if (!first || !first.ff) throw new Error("runPullingEnsemble: supply {ff, pos0} or makeSystem().");
  const ex = o.dir ? { dir: o.dir, ligCOM0: ligandCOM(first.pos0 ?? first.ff.ref, first.ff.ligandStart ?? first.ff.nProt, first.ff.n - (first.ff.ligandStart ?? first.ff.nProt)) }
    : exitVector(first.ff, { rPocket: o.rPocket });
  const run = mode === "velocity" ? runConstantVelocityPull : runConstantForcePull;
  const keepRestraints = o.keepRestraints === true;
  const pulls = [];
  for (let p = 0; p < nPulls; p++) {
    const sys = o.makeSystem ? o.makeSystem(p) : first;
    const sysEx = o.dir ? ex : (p === 0 ? ex : exitVector(sys.ff, { rPocket: o.rPocket }));
    // Drop bound-state restraints for the pull; restore afterwards.
    const saved = keepRestraints ? null : dropRestraints(sys.ff);
    try {
      pulls.push(run({
        ...(o.pull ?? {}), ff: sys.ff, pos0: sys.pos0 ?? sys.ff.ref,
        dir: sysEx.dir, ligCOM0: sysEx.ligCOM0, seed: (o.seed ?? 42) + p,
      }));
    } finally {
      if (saved) restoreRestraints(sys.ff, saved);
    }
  }
  const works = Float64Array.from(pulls.map((r) => r.work));
  const meanWork = works.reduce((a, w) => a + w, 0) / works.length;
  return { mode, pulls, works, meanWork, dir: ex.dir, ligCOM0: ex.ligCOM0 };
}

/**
 * Jarzynski free energy from an ensemble of nonequilibrium works:
 *   ΔF = Wmin − (1/β)·ln⟨exp(−β(W − Wmin))⟩,  β = 1/(kB·T),
 * with a seeded nonparametric bootstrap standard error.
 *
 * @param {Array<number>|Float64Array} works  kcal/mol
 * @param {object} [opts]
 * @param {number} [opts.T=300]  K
 * @param {number} [opts.bootstrap=200]  resamples (0 disables)
 * @param {number} [opts.seed=1]
 * @returns {{dF: number, meanWork: number, minWork: number, maxWork: number,
 *   dissipated: number, n: number, se: number|null, beta: number}}
 */
export function jarzynskiFreeEnergy(works, opts = {}) {
  const T = opts.T ?? 300;
  const W = Array.from(works).filter(Number.isFinite);
  if (!W.length) throw new Error("jarzynskiFreeEnergy: no finite works.");
  const beta = 1 / (KB_KCAL * T);
  const wMin = Math.min(...W);
  const avg = W.reduce((a, w) => a + Math.exp(-beta * (w - wMin)), 0) / W.length;
  const dF = wMin - Math.log(Math.max(1e-300, avg)) / beta;
  const meanWork = W.reduce((a, w) => a + w, 0) / W.length;
  let se = null;
  const B = opts.bootstrap ?? 200;
  if (B > 0 && W.length > 1) {
    const rng = makeRng(opts.seed ?? 1);
    const est = [];
    for (let b = 0; b < B; b++) {
      let s = 0;
      for (let i = 0; i < W.length; i++) s += Math.exp(-beta * (W[Math.floor(rng() * W.length)] - wMin));
      est.push(wMin - Math.log(Math.max(1e-300, s / W.length)) / beta);
    }
    const m = est.reduce((a, v) => a + v, 0) / est.length;
    se = Math.sqrt(est.reduce((a, v) => a + (v - m) ** 2, 0) / est.length);
  }
  return {
    dF, meanWork, minWork: wMin, maxWork: Math.max(...W),
    dissipated: meanWork - dF, n: W.length, se, beta,
  };
}

/**
 * Rupture-based koff ranking surrogate (unitless ordering proxy, NOT a rate).
 * Scores are reported in kT (dissipated work beyond reversible, same quantity
 * compared to itself):
 *   score_i = β·(Wrupt,i − ΔFrupt),   ΔFrupt = Jarzynski({Wrupt}),  β = 1/kT.
 * Near 0 kT ⇒ near-reversible escape; larger ⇒ more dissipative escape path.
 * kT units keep fast-protocol ensembles comparable (no exp() blow-up); the
 * ranking is identical to the exp(β·ΔW) form by monotonicity. The `dF`
 * argument (e.g. total-work Jarzynski ΔF) is reported for context as inputDF;
 * the score reference refDF is computed from rupture works.
 * @param {PullResult[]} pulls
 * @param {number} dF  Jarzynski free energy of the total-work ensemble (kcal/mol)
 * @param {object} [opts]
 * @param {number} [opts.T=300]
 * @returns {{scores: number[], meanScore: number, refDF: number, inputDF: number,
 *   meanRuptureForce: number, meanWorkToRupture: number, note: string}}
 */
export function koffSurrogate(pulls, dF, opts = {}) {
  const T = opts.T ?? 300;
  const beta = 1 / (KB_KCAL * T);
  if (!pulls || !pulls.length) throw new Error("koffSurrogate: empty pulls.");
  if (!Number.isFinite(dF)) throw new Error("koffSurrogate: dF must be finite.");
  const ruptWorks = pulls.map((p) => (Number.isFinite(p.workToRupture) ? p.workToRupture : p.work));
  if (!ruptWorks.every(Number.isFinite)) throw new Error("koffSurrogate: non-finite rupture works.");
  const refDF = jarzynskiFreeEnergy(ruptWorks, { T, bootstrap: 0 }).dF;
  const scores = ruptWorks.map((w) => beta * (w - refDF)); // kT units
  return {
    scores,
    meanScore: scores.reduce((a, s) => a + s, 0) / scores.length,
    refDF, inputDF: dF,
    meanRuptureForce: pulls.reduce((a, p) => a + p.ruptureForce, 0) / pulls.length,
    meanWorkToRupture: ruptWorks.reduce((a, w) => a + w, 0) / ruptWorks.length,
    note: "ranking surrogate only (kT dissipated past reversible) — larger score ⇒ more dissipative escape; absolute koff needs an MSM (see docs/NETWORK.md)",
  };
}
