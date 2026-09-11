# Heavy-Mode, Physics, Performance & UI — Engineering Plan

Status: **PLAN (pre-implementation)**

This plan is grounded in the current codebase. Before writing it I read the
full source tree and ran two probes under Node:

- Heavy mode parses 4W52 into **1308 heavy atoms**.
- `HeavyForceField.compute()` takes **~51 ms/call** on this machine (Node),
  with `dt = 4 fs` ⇒ about **0.08 ps of simulation per second** of wall-clock.
  Heavy mode is currently ~100× too slow to be usable, independent of GPU.
- The MOL2 path *constructs* correctly in isolation (benzene.mol2 appends 6
  atoms), but in the real app it is broken by **double-counting HETATM
  ligands** and by the heavy mode **ignoring the "Include ligands" checkbox**.

---

## 0. Root-cause diagnosis

### 0.1 Heavy mode ignores the ligand/HETATM checkbox

`buildSystem()` in `src/main.js` calls `selectHeavy(state.parsedHeavy, ...)`
and `parseHeavy()` (`src/heavy.js`) keeps **every** HETATM except water. The
`ui.includeLig` checkbox is only consulted in the CG branch. Consequence:

- Loading 4W52 (which contains crystallographic `BNZ` + `HEPES`) in heavy mode
  **always** includes those two HETATM molecules, even when the user unchecks
  "Include ligands".
- Then loading `benzene.mol2` **appends a second benzene** on top of the
  crystallographic benzene ⇒ native overlap ⇒ the 4.5 kcal/mol → blow-up /
  non-finite path the user is seeing ("MOL2 fails to load").

### 0.2 No hetero-atom / metal selection

`parseHeavy()` detects metals via `METAL_ELEMENT` but there is no UI to select
which metals / cofactors / hetero atoms are included. Everything non-water is
silently included. Requirement unmet: **detect + let the user select**.

### 0.3 `nProt` / `nLigAtoms` semantics are wrong in heavy mode

`HeavyForceField` computes:

```js
this.nProt     = atoms.filter(a => a.isProtein).length;
this.nLigAtoms = this.n - this.nProt;   // = ALL non-protein atoms
```

So metals, cofactors and the appended MOL2 ligand are all lumped into
"ligand". The funnel (`src/funnel.js`) and HUD assume ligand atoms are the
*last* `nLigAtoms` particles — but in heavy mode the metal/cofactor atoms come
*before* the appended MOL2 ligand in `state.sel.atoms`. This mis-targets the
funnel CV, `ligRMSD`, `nContacts`, and `U_bind`.

### 0.4 Heavy mode is O(n²) + finite-difference torsions

`HeavyForceField._nonBonded()` is a pure O(n²) double loop over ~1308 atoms
(~855k pairs / call). Improper and proper dihedrals use **central finite
differences** (`improperFlat`/`properFlat`), i.e. 24 `improperAngle` calls per
dihedral; 4W52 has 4160 propers + 442 impropers ⇒ ~110k dihedral evaluations
per force call. Combined with 250 steps/frame, heavy mode is unusable.

### 0.5 No multi-core / GPU path

The entire app is single-threaded and CPU-bound. There is no worker pool and
no WebGPU/WebGL compute path.

### 0.6 Physics is visualization-grade, not binding-grade

Heavy mode uses element LJ + `ε(r) = r` screened Coulomb with **mean-charge
neutralization**, no implicit-solvent solvation (GB/SASA), a radial-only
Gaussian H-bond, and no real metal chemistry beyond harmonic coordination. CG
mode has a slightly better binding term (EEF1-lite + screened Coulomb +
radial H-bond) but still not a defensible binding-physics model.

### 0.7 UI has no settings page

All controls (physics + display + future compute backend) are mixed into the
sidebar. There is no place for "basic settings" (backend, threads, GPU,
quality) vs "tuning" sliders.

---

## 1. Design direction (the "what we build")

Keep the app **browser-native, no build step, no backend**, but split it into
two complementary models and a pluggable compute layer:

1. **Cα ENM** (existing) — fast global fold dynamics. Stays the default.
2. **Heavy mode** — becomes a *binding zoom* model:
   - protein **residue-rigid-body / torsion-space** dynamics rather than
     full all-atom stiff-bond MD (fewer DOF, larger timestep, faster, more
     stable), OR keep explicit heavy atoms but with analytic forces + cell
     lists + workers/GPU so it actually runs;
   - **implicit solvent** (pairwise Generalized Born or distance-dependent
     dielectric + Debye salt screening) + **SASA** non-polar term;
   - **directional H-bonds** and explicit **metal selection**;
   - a **kinetic / chemical network mode** for binding events (Markov state
     model / milestoning-lite) that visualizes binding/unbinding as
     rate-governed transitions between metastable poses — far cheaper than
     integrating every atom, and physically interpretable.

3. **Compute backends** (pluggable): CPU-single, CPU-multi (Web Workers),
   GPU-WebGPU, GPU-WebGL2. A Settings page selects the backend, thread count,
   device, and quality preset.

**Answer to "all-atom but faster than MD / chemical network":**
Full all-atom MD at interactive speed in the browser is not realistic for
1300+ atoms. The pragmatic route is a **multi-scale hybrid**:

- ENM for global protein motion,
- heavy atoms only where binding matters (pocket side chains + ligand) with
  torsion-space / rigid-body dynamics,
- implicit solvent + electrostatics for physically grounded energies,
- and a **kinetic network (MSM)** layered on top for binding/unbinding rates
  and visualization.

This is exactly how modern visualization + biophysics tools balance speed and
realism.

---

## 2. Phased implementation plan

### Phase 0 — Correctness: heavy MOL2, hetero/metal selection, placement

**Files:** `src/heavy.js`, `src/pdb.js`, `src/main.js`, `src/ui.js`,
`src/ligand-panel.js`, `index.html`, `css/style.css`.

- **0.1** `parseHeavy()` returns an explicit hetero classification:
  `{ protein: [...], hetero: [...] }` with per-atom
  `{isProtein, isMetal, isWater, resName, chain, resSeq, element}` and a
  summary of detected hetero groups (`resName → count`). Keep metals, keep
  cofactors, but tag them (not "ligand").
- **0.2** Add a **Hetero & Metals** selection UI (panel 2, heavy mode only):
  a list of detected hetero residues (e.g. `ZN`, `HEM`, `BNZ`, `HEPES`, `SO4`,
  `GOL` …) with checkboxes. Defaults: **metals ON**, **small organic cofactors
  OFF when a MOL2/library ligand is supplied**, water/ions OFF. Heavy mode
  builds only checked hetero atoms.
- **0.3** Fix **MOL2/library override semantics** in heavy mode to mirror CG:
  when an external ligand (library or MOL2) is active, **exclude** the PDB's
  organic HETATM *ligand* groups (but keep metals + explicitly selected
  cofactors). No more double benzene.
- **0.4** Split `HeavyForceField` particle bookkeeping:
  `nProt` = protein atoms, `nHetero` = selected cofactor/metal atoms,
  `nLigAtoms` = appended external ligand atoms. Reorder particles as
  **protein → hetero → ligand** so the funnel/HUD/`ligRMSD`/`nContacts`
  correctly target the external ligand only. Keep `ligandAtoms`/`ligandBonds`
  populated for the viewer and integrator.
- **0.5** Generalize placement to heavy mode **and** to any MOL2 molecule
  (not just the built-in library): a "Place MOL2 ligand" action that runs the
  same rigid clash-relaxation (`src/placement.js`) against the heavy-atom
  protein using per-atom `_elem[i].sigma`. Verify clash-free placement is
  reported correctly.
- **0.6** End-to-end acceptance tests (Node + puppeteer):
  - protein-only PDB + `benzene.mol2` in heavy mode ⇒ 1 ligand, no duplicate;
  - 4W52 + `benzene.mol2` ⇒ 1 benzene, no double count;
  - hetero/metal selection persists across rebuilds;
  - placed ligand is clash-free and `U_bind` is finite.

**Exit gate:** heavy mode loads MOL2/library ligands correctly, metals/hetero
are user-selectable, placement works, no NaN/blow-up on the demo systems.

---

### Phase 1 — CPU performance (make heavy mode actually run)

**Files:** `src/heavy.js`, new `src/neighbors.js`, `src/ff-harmonic.js`.

- **1.1** Replace O(n²) `_nonBonded()` with a **cell-list / spatial hash**
  (reuse the pattern already proven in `forcefield.js` `_repulsion` /
  `_binding`). Cutoff 8 Å ⇒ O(N) average. This alone is ~10–30×.
- **1.2** Add **analytic gradients** for improper and proper torsions
  (shared kernel in `ff-harmonic.js` or `heavy.js`), replacing the
  finite-difference kernels. ~20–50× on the bonded terms.
- **1.3** Add a **Verlet neighbor list** with a skin (e.g. 2 Å) rebuilt every
  N steps, so the pair list is not re-derived every call; benchmark grid
  rebuild vs. Verlet and keep the faster.
- **1.4** Replace the `Set`/`Map` exclusion lookups (`pairKey` numbers) with
  typed-array/sorted-pair structures for 1-2/1-3/1-4 exclusions, or a
  per-particle excluded-neighbor bitmask for the common cases.
- **1.5** Target: `compute()` ≤ **2 ms** for ~2k heavy atoms on a desktop
  CPU (currently 51 ms for 1308), so ~1 ps/frame at 60 fps is reachable.

**Exit gate:** benchmark script shows ≥ 20× compute speedup and the 4W52
heavy-mode demo runs at interactive fps on CPU.

---

### Phase 2 — Multi-core CPU (Web Workers)

**Files:** new `src/compute-pool.js`, `src/force-worker.js`, `src/main.js`,
`src/settings-panel.js`.

- **2.1** Extract force evaluation into a **worker protocol**:
  `{type:"init", ff}` → `{type:"compute", pos}` → `{forces, energy, bindingU,
  desolvU, ...}`. Integrator stays on the main thread (it is O(N) and cheap).
- **2.2** Implement a **worker pool** of size
  `navigator.hardwareConcurrency - 1` (clamped). Two parallelization schemes,
  chosen by benchmark:
  - **term-parallel**: bonded terms on one worker, non-bonded on others, then
    reduce forces; or
  - **domain-parallel**: split the non-bonded cell list by cell ranges.
- **2.3** Use transferable `Float64Array` copies. Optionally support
  **SharedArrayBuffer** when COOP/COEP headers are present (provide the
  server-header instructions and a settings toggle).
- **2.4** Settings page exposes **CPU threads: Auto / 1 / 2 / … / N**.

**Exit gate:** multi-core path is measurably faster than single-thread and
produces forces within float tolerance of the reference.

---

### Phase 3 — GPU backend (WebGPU compute, WebGL2 fallback)

**Files:** new `src/gpu.js`, `src/gpu-wgsl.js`, `src/gpu-glsl.js`,
`src/settings-panel.js`.

- **3.1** **WebGPU** compute shader for the non-bonded term
  (LJ + screened Coulomb + switch), plus a GPU-friendly neighbor/pair list
  built on CPU and uploaded. Positions uploaded each step, forces/energy read
  back. Use WGSL strings (no build step).
- **3.2** **WebGL2** fallback (float textures / transform feedback) for
  browsers without WebGPU. This is the riskiest/largest chunk; gate it behind
  "if time", since WebGPU covers modern Chrome/Edge.
- **3.3** Settings page: **Compute backend: Auto / CPU / CPU (multi-core) /
  GPU (WebGPU) / GPU (WebGL2)**; auto-detect `navigator.gpu`.
- **3.4** Force-parity test: GPU vs CPU must agree to `< 1e-4 kcal/mol/Å`
  per component before the backend is selectable.

**Exit gate:** WebGPU backend runs the 4W52 heavy demo at 60 fps; auto
fallback to CPU/workers when unavailable.

---

### Phase 4 — Physics realism (biophysics layer)

**Files:** new `src/physics/gb.js`, `src/physics/sasa.js`,
`src/physics/hbond.js`, `src/physics/charges.js`, `src/physics/network.js`;
integrate into `forcefield.js` and `heavy.js`.

- **4.1 Charges.** Replace the "mean-neutralize every atom" hack with
  **residue/atom partial charges** (AMBER-family table for the 20 standard
  amino acids + common cofactors; ligand charges from MOL2 when present, else
  element defaults). Enforce neutrality with a background neutralizing term /
  explicit counterions where needed.
- **4.2 Electrostatics + implicit solvent.** Implement a **pairwise
  Generalized Born (GB)** term (Still/HCT/GB-OBC) for solvation + a
  **Debye–Hückel salt screening** for ionic strength, or (fallback) a
  physically motivated distance-dependent dielectric. GB is O(n²) but, with
  the Phase 1 neighbor list + cutoff, feasible for the target sizes.
- **4.3 Non-polar solvation.** Replace EEF1-lite with a proper **SASA** term
  (LCPO or Shrake–Rupley) and `U_np = γ·SASA + b`. Reuse the cell list.
- **4.4 H-bonds.** Add an **angle-dependent** donor-H-acceptor term (radial
  + donor/acceptor angular windows) instead of the radial Gaussian.
- **4.5 Metals/hetero.** Explicit, selectable metal coordination with an
  optional directional term; keep the harmonic coordination as the default
  cheap option.
- **4.6 Validation.** Add tests that recover, within an order of magnitude,
  known binding energies (e.g. 4W52 benzene ΔG ≈ −5 kcal/mol), B-factor
  correlation, and RMSIP. Document every approximation in the energy model.

**Exit gate:** heavy-mode energies are physically defensible (documented
terms, no mean-charge hack, solvation + electrostatics present) and validation
tests pass.

---

### Phase 5 — Kinetic / chemical network mode (fast binding dynamics)

**Files:** new `src/physics/network.js`, `src/network-panel.js`, integrate
with `analysis.js`.

- **5.1** Define binding CVs: ligand↔pocket COM distance, native-contact
  count, ligand RMSD, burial/SASA. (The funnel already gives the 1-D PMF.)
- **5.2** Build a **Markov state model (MSM)** over clustered ligand poses:
  - source microstates from many short CG/heavy trajectories (offline or
    in-browser), **or**
  - approximate the transfer operator from the reconstructed PMF + a
    diffusive (Smoluchowski / Kramers) model.
- **5.3** Estimate **binding/unbinding rates** (mean first-passage times) from
  the transition matrix.
- **5.4** Add a **kinetic Monte Carlo (Gillespie) visualization**: the ligand
  hops between metastable states with the MSM rates. This is the "chemical
  network" layer — it animates binding/unbinding as rate-governed events and
  is orders of magnitude faster than integrating every atom, while remaining
  physically interpretable.
- **5.5** Expose the network in the Analysis panel: state populations, rates,
  transition matrix, and a network graph (optional).

**Exit gate:** the network mode visualizes a binding/unbinding cycle with
physically sensible rates and a documented connection to the underlying
energy model.

---

### Phase 6 — UI reorganization & Settings page

**Files:** `index.html`, `css/style.css`, `src/ui.js`, new
`src/settings-panel.js`.

- **6.1** Add a **⚙ Settings** button (header or toolbar) opening a
  modal/settings view containing **basic settings**:
  - Compute backend (Auto/CPU/CPU-multi/GPU-WebGPU/GPU-WebGL2)
  - CPU thread count
  - GPU device (if multiple)
  - Quality preset (Low / Medium / High) → timestep, cutoff, skin, neighbor
    rebuild frequency, sphere quality
  - Physics toggles (implicit solvent GB, SASA, H-bonds, electrostatics,
    kinetic network)
  - Display defaults
  Persist to `localStorage`.
- **6.2** Move the remaining **tuning** controls out of the "basic" set and
  keep them in the sidebar as a compact **Tuning** group: `Rc`, `γ`, `T`, `ζ`,
  mass, motion gain, binding/holo/funnel toggles.
- **6.3** Reorganize sidebar panels:
  1 Structure, 2 Selection & Hetero, 3 Force Field & Integrator (tuning),
  Ligand Library & Placement, Recording, PMF, Analysis. Make heavy vs CG
  controls context-sensitive.
- **6.4** Responsive/accessibility pass; keep the header/HUD clean.

**Exit gate:** settings live in one place, tuning lives in the sidebar,
state persists, and switching backends/threads is a one-click operation.

---

## 3. Cross-cutting work

- **Testing** — Node unit tests for every new pure module (`gb`, `sasa`,
  `hbond`, `charges`, `network`, `neighbors`, `gpu` reference parity);
  puppeteer end-to-end for MOL2/heavy/selection/placement/backend.
- **Benchmarks** — a `bench/` script reporting `compute ms`, `steps/s`, and
  `fps` for 4W52 in both modes, before/after each phase.
- **Docs** — update `README.md` energy model, heavy-mode section, and the
  settings page after each phase.

---

## 4. Recommended execution order & dependencies

1. **Phase 0** (correctness) — no dependencies; do first. The MOL2/hetero bugs
   block everything else.
2. **Phase 1** (CPU perf) — depends on 0; makes heavy mode testable.
3. **Phase 6 (partial)** — add the Settings page shell + backend selector
   early, so Phase 2/3 have a place to land.
4. **Phase 2** (workers) and **Phase 3** (GPU) — can proceed in parallel after
   Phase 1; GPU is the higher-impact, higher-risk item.
5. **Phase 4** (physics) — can start after Phase 1 (independent of 2/3), but
   its terms should be written so they also run on the GPU later.
6. **Phase 5** (network) — after Phase 4 (needs the energy model + CVs).

**First milestone (recommended):** Phase 0 + Phase 1 + the Phase 6 settings
shell. This fixes the reported bugs and makes heavy mode fast enough to use,
before investing in GPU and advanced physics.

---

## 5. Open questions to resolve during implementation

- **GB vs distance-dependent dielectric** — GB is more physical but O(n²);
  we should benchmark GB with the Phase 1 neighbor list before committing.
- **SharedArrayBuffer** — requires COOP/COEP headers; decide whether to ship a
  small server config or stay with transferable copies.
- **WebGL2 fallback scope** — WebGPU is available in current Chrome/Edge;
  WebGL2 may be deferred unless wide support is required.
- **Network source** — MSM from short trajectories vs. PMF-derived rates; the
  latter is cheaper and sufficient for a first version.
- **Hetero default policy** — metals ON, cofactors OFF-when-ligand-supplied,
  water/ions OFF; confirm this matches user expectations.
