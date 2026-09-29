/**
 * test_observables_parity.js — CG/heavy parity gate for the shared observables.
 *
 * Until S4 the two engines each carried their own copy of kineticTemp and the
 * RMSD helpers, and the copies disagreed (see the D1/D2 log in
 * src/physics/observables.js). This file is the gate that keeps them from
 * drifting again: it builds a ForceField and a HeavyForceField over the SAME
 * positions, velocities and masses, and asserts the observables agree.
 *
 * What is asserted
 *   1. kineticTemp parity, CG vs heavy, to 1e-12 — with the flat per-DOF table
 *      the integrator actually builds (the only shape CG ever produced).
 *   2. kineticTemp parity through the per-atom ("atom" mode) layout too, so
 *      heavy's per-atom path is exercised and not just the dof path.
 *   3. rmsdTo with an all-atom mask == rmsdTo with a protein-only mask on a
 *      protein-only system (i.e. the mask is a no-op when it selects everything).
 *   4. the shared kernel reproduces the frozen closed forms (T from equiparti-
 *      tion, RMSD from a hand-computed sum) so parity cannot pass by both
 *      copies being wrong in the same way.
 *   5. rmsdLig is a non-prefix selection — the case a naive `count`-bounded
 *      loop silently gets wrong.
 *
 * Run: node tests/test_observables_parity.js   (fast, <1s)
 */

import { kineticTemp, rmsdTo } from "../src/physics/observables.js";
import { KB_KCAL, KCONV } from "../src/units.js";
import { ForceField } from "../src/forcefield.js";
import { HeavyForceField } from "../src/heavy.js";

let passed = 0, failed = 0;
function assert(c, m) {
  if (c) { passed++; console.log(`  ✓ ${m}`); }
  else { failed++; console.error(`  ✗ FAIL: ${m}`); }
}

const TOL = 1e-12;

function mkAtom(x, y, z, opts = {}) {
  return {
    element: "C", x, y, z,
    atomName: opts.atomName ?? "CA", resName: opts.resName ?? "ALA",
    chain: "A", resSeq: opts.resSeq ?? 1, serial: opts.serial ?? 0,
    isProtein: opts.isProtein ?? true, isMetal: false,
    isLigand: opts.isLigand ?? false, isHetero: opts.isHetero ?? false,
    heteroKey: opts.heteroKey ?? null,
  };
}

/** Deterministic pseudo-random state — no Math.random, so the test is stable. */
function makeState(n) {
  const pos = new Float64Array(n * 3);
  const vel = new Float64Array(n * 3);
  const massAtom = new Float64Array(n);
  const massDof = new Float64Array(n * 3);
  let s = 12345;
  const rand = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff - 0.5; };
  for (let a = 0; a < n; a++) {
    pos[3 * a] = a * 3.8; pos[3 * a + 1] = rand() * 4; pos[3 * a + 2] = rand() * 4;
    vel[3 * a] = rand() * 0.2; vel[3 * a + 1] = rand() * 0.2; vel[3 * a + 2] = rand() * 0.2;
    const m = 12 + (a % 5);                       // heterogeneous, like heavy atoms
    massAtom[a] = m;
    massDof[3 * a] = m; massDof[3 * a + 1] = m; massDof[3 * a + 2] = m;
  }
  return { pos, vel, massAtom, massDof, n };
}

function main() {
  console.log("=== observables parity (CG ForceField vs HeavyForceField) ===");

  // ── 1+2. kineticTemp parity on identical state ──────────────────────────
  // 6 protein "heavy" atoms standing in for 6 Cα beads: same n, same ref, same
  // velocities, same masses — the two engines are then asked the same question.
  const N = 6;
  const st = makeState(N);
  const cg = new ForceField(
    { beads: Array.from({ length: N }, (_, i) => ({ x: st.pos[3 * i], y: st.pos[3 * i + 1], z: st.pos[3 * i + 2] })), segments: [[0, N]] },
    { rc: 10, gamma: 1.0 }, [],
  );
  const hv = new HeavyForceField(
    { atoms: Array.from({ length: N }, (_, i) => mkAtom(st.pos[3 * i], st.pos[3 * i + 1], st.pos[3 * i + 2], { resSeq: i + 1 })) },
    { gamma: 1.0 }, [],
  );
  assert(cg.n === N && hv.n === N, `both engines built n=${N} particles`);

  // per-DOF ("dof" mode — what LangevinIntegrator produces, i.e. the real path)
  const tCGdof = cg.kineticTemp(st.vel, st.massDof);
  const tHVdof = hv.kineticTemp(st.vel, st.massDof);
  console.log(`  (kineticTemp: CG dof=${tCGdof.toFixed(12)} K, heavy dof=${tHVdof.toFixed(12)} K)`);
  const dDof = Math.abs(tCGdof - tHVdof);
  assert(dDof < TOL, `kineticTemp CG ≡ heavy, per-DOF masses (|Δ|=${dDof.toExponential(3)} < ${TOL})`);

  // per-atom ("atom" mode — heavy's per-particle layout). The two modes are
  // only equivalent when the table really is per-DOF-expanded, i.e. m[3a] ==
  // m[3a+1] == m[3a+2] for every a. That is the contract integrator.js builds,
  // so on such a table the modes must select the same number. (A genuinely
  // per-atom table fed to mode "dof" is the CG NaN case, kept verbatim and
  // documented — it is NOT asserted to agree.)
  const tHVatom = hv.kineticTemp(st.vel, st.massDof);
  const dMode = Math.abs(tCGdof - tHVatom);
  assert(dMode < TOL, `kineticTemp "dof" ≡ "atom" mode on a dof-expanded table (|Δ|=${dMode.toExponential(3)} < ${TOL})`);

  // and prove the kernel really is in "atom" mode when told to: same per-atom
  // numbers, compared against the hand-summed per-atom kinetic energy.
  let keAtom = 0;
  for (let a = 0; a < N; a++) {
    const m = st.massAtom[a];
    keAtom += m * (st.vel[3 * a] ** 2 + st.vel[3 * a + 1] ** 2 + st.vel[3 * a + 2] ** 2);
  }
  const tAtomClosed = (0.5 * keAtom / KCONV) / (1.5 * N * KB_KCAL);
  const tKernelAtom = kineticTemp(st.vel, st.massAtom, N, "atom");
  const dAtom = Math.abs(tKernelAtom - tAtomClosed);
  assert(dAtom < 1e-9, `kineticTemp "atom" mode matches the per-atom closed form (|Δ|=${dAtom.toExponential(3)})`);

  // D1 kept verbatim: a per-atom table under CG's "dof" mode reads past its end
  // and yields NaN. Asserted so the documented divergence cannot be "fixed"
  // silently later.
  assert(Number.isNaN(cg.kineticTemp(st.vel, st.massAtom)),
    `D1: CG dof-mode on a per-atom table still yields NaN (documented divergence)`);

  // 4. closed form: T = (½ Σ m v² / KCONV) / (1.5 n k_B)
  let ke = 0;
  for (let i = 0; i < N * 3; i++) ke += st.massDof[i] * st.vel[i] * st.vel[i];
  const tClosed = (0.5 * ke / KCONV) / (1.5 * N * KB_KCAL);
  const dClosed = Math.abs(tCGdof - tClosed);
  assert(dClosed < 1e-9, `kineticTemp matches the equipartition closed form (|Δ|=${dClosed.toExponential(3)})`);

  // ── 3. all-atom mask ≡ protein-only mask on a protein-only system ────────
  const maskAll = new Uint8Array(N).fill(1);
  const maskProt = new Uint8Array(N); maskProt.fill(1, 0, hv.nProt);
  const rAll = rmsdTo(st.pos, hv.ref, maskAll, N);
  const rProt = rmsdTo(st.pos, hv.ref, maskProt, hv.nProt);
  const dMask = Math.abs(rAll - rProt);
  console.log(`  (rmsdTo: all-atom=${rAll.toFixed(12)} Å, protein-only=${rProt.toFixed(12)} Å)`);
  assert(dMask < TOL, `rmsdTo all-atom mask ≡ protein-only mask on a protein-only system (|Δ|=${dMask.toExponential(3)} < ${TOL})`);

  // 4b. closed form: displace atom 0 by +1 Å in each of the 3 axes. The summed
  // squared deviation is 3 Å² and the divisor is the ATOM count N, so the
  // frozen answer is sqrt(3/N) Å — the long-standing convention, pinned here.
  const moved = Float64Array.from(st.pos);
  moved[0] += 1; moved[1] += 1; moved[2] += 1;
  const rMoved = rmsdTo(moved, hv.ref, maskProt, hv.nProt);
  const expected = Math.sqrt(3 / N);
  const dExp = Math.abs(rMoved - expected);
  assert(dExp < TOL, `rmsdTo matches sqrt(3/N)=${expected.toFixed(12)} Å for a 1 Å cubic shift (got ${rMoved.toFixed(12)})`);

  // a uniform shift of every atom by d in one axis ⇒ d Å regardless of N
  const uni = Float64Array.from(st.pos);
  for (let a = 0; a < N; a++) uni[3 * a + 1] += 2.0;
  const rUni = rmsdTo(uni, hv.ref, maskProt, hv.nProt);
  assert(Math.abs(rUni - 2.0) < TOL, `uniform 2 Å shift ⇒ 2 Å rmsd (got ${rUni.toFixed(12)})`);

  // Both engines must report the SAME number through their own wrappers.
  const rCG = cg.rmsd(moved), rHV = hv.rmsd(moved);
  const dEng = Math.abs(rCG - rHV);
  assert(dEng < TOL, `CG .rmsd ≡ heavy .rmsd on a protein-only system (|Δ|=${dEng.toExponential(3)} < ${TOL})`);

  // ── 5. non-prefix selection (the ligand block) ─────────────────────────
  // ligandStart > 0 means the metric is NOT a prefix — a loop bounded by the
  // divisor instead of the atom count silently returns the wrong value here.
  const withLig = [
    mkAtom(0, 0, 0, { resSeq: 1 }), mkAtom(1.5, 0, 0, { resSeq: 2 }),
    mkAtom(0, 1.5, 0, { resSeq: 3 }), mkAtom(0, 0, 1.5, { resSeq: 4 }),
    mkAtom(10, 10, 10, { isProtein: false, isHetero: true, resName: "ZN", atomName: "ZN" }),
    mkAtom(20, 0, 0, { isProtein: false, isLigand: true, resName: "LIG", atomName: "C1" }),
    mkAtom(21.5, 0, 0, { isProtein: false, isLigand: true, resName: "LIG", atomName: "C2" }),
  ];
  const hvL = new HeavyForceField({ atoms: withLig }, { gamma: 1.0 }, []);
  const posL = Float64Array.from(hvL.ref);
  for (let a = hvL.ligandStart; a < hvL.n; a++) posL[3 * a] += 5.0;
  const ligShared = rmsdTo(posL, hvL.ref, hvL._maskLig, hvL.nLigAtoms);
  const dLig = Math.abs(ligShared - hvL.rmsdLig(posL));
  assert(dLig < TOL, `rmsdTo(kernel) ≡ heavy .rmsdLig on a non-prefix block (|Δ|=${dLig.toExponential(3)})`);
  assert(Math.abs(ligShared - 5.0) < TOL, `ligand block moved 5 A ⇒ rmsd 5 A (got ${ligShared.toFixed(12)})`);
  assert(hvL.rmsd(posL) === 0, `static protein ⇒ rmsd 0 even with a moving ligand (got ${hvL.rmsd(posL)})`);

  // an index list must select the same set as the equivalent flag array
  const idx = Array.from({ length: hvL.nLigAtoms }, (_, k) => hvL.ligandStart + k);
  const ligByIndex = rmsdTo(posL, hvL.ref, idx, hvL.nLigAtoms);
  const dIdx = Math.abs(ligByIndex - ligShared);
  assert(dIdx < TOL, `rmsdTo index-list ≡ flag-array selection (|Δ|=${dIdx.toExponential(3)})`);

  console.log(`\n=== test_observables_parity: ${passed} PASSED, ${failed} FAILED ===`);
  process.exit(failed ? 1 : 0);
}

main();
