# Adding a physics term without touching the integrator

A **term** is one contribution to the potential energy: it reads the positions,
adds its own forces to the shared force buffer, and returns its energy. In this
repo every term is a small declarative descriptor in a registry, so a new one is
a data object — not an edit inside `compute()`.

- **The contract:** [`src/physics/terms/registry.js`](../src/physics/terms/registry.js)
  (leaf module, zero imports).
- **The built-in CG terms:** [`src/cg/terms.js`](../src/cg/terms.js).
- **The built-in heavy terms:** [`src/heavy/terms.js`](../src/heavy/terms.js) plus
  `heavy.weak` in [`src/heavy/weak.js`](../src/heavy/weak.js).
- **The walk that runs them:** [`src/cg/compute.js`](../src/cg/compute.js) and
  [`src/heavy/energy.js`](../src/heavy/energy.js).
- **The proof:** [`tests/test_physics_terms.js`](../tests/test_physics_terms.js)
  registers a term from the test file and asserts analytic numbers.
- **Deeper context:** [PHYSICS_RIGOR.md](PHYSICS_RIGOR.md) (what each term is),
  [CG_HEAVY.md](CG_HEAVY.md) (the two resolutions).

## The interface, in full

```js
{
  id:      "cg.myWell",            // "engine.name", unique, stable forever
  engine:  "cg",                   // "cg" | "heavy" — which compute() owns it
  order:   110,                    // integer; ascending = position in the U sum
  label:   "My restraint",         // one line for the HUD and the docs
  unit:    "kcal/mol",             // optional, defaults to kcal/mol

  // null (the default) ⇒ always on. Otherwise `function () { … this … }`.
  // A false predicate SKIPS the term: energy() is never called and U is untouched.
  enabled: function () { return this.myFlag === true; },

  // REQUIRED. `this` IS the force field. MUST accumulate −∇U into `f`, and MUST
  // return its own ΔU as ONE float64 (see "Grouping" below). No allocation.
  energy: function (pos, f) {
    let U = 0;
    for (let a = 0; a < this.pairs.length; a += 3) {
      const i = 3 * this.pairs[a], j = 3 * this.pairs[a + 1];
      const dx = pos[j] - pos[i], dy = pos[j + 1] - pos[i + 1], dz = pos[j + 2] - pos[i + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
      const dr = r - this.pairs[a + 2];
      U += 0.5 * this.kWell * dr * dr;
      const s = (this.kWell * dr) / r;
      f[i] += s * dx; f[i + 1] += s * dy; f[i + 2] += s * dz;
      f[j] -= s * dx; f[j + 1] -= s * dy; f[j + 2] -= s * dz;
    }
    this.myWellU = U;              // publish detail on the FIELD, not in the return
    return U;
  },

  // Optional, documentation only — never read in the hot path.
  reports: ["myWellU"],            // names of the ff accumulators this term publishes
  report:  function () { return this.myWellU; },   // how describeTerms() reads it back
}
```

Registry API: `registerTerm(t)`, `unregisterTerm(id)`, `resetTerms(engine)`,
`listTerms(engine)`, `getTerm(id)`, `describeTerms(engine, ff)`,
`validateTerm(t)`, `TERM_GENERATION`, `TERM_ENGINES`.

`energy` must be a `function (pos, f) {}`, not an arrow: the force field arrives
as `this`, and `registerTerm` **rejects** an arrow (`energy.prototype ===
undefined`) rather than let it fail silently at step 1 of a trajectory.

## Worked example — a complete, working term

Copy this; it is the test fixture, verbatim and runnable.

```js
import { registerTerm } from "./src/physics/terms/registry.js";

const K = 2.0;              // kcal/mol/Å²
const ANCHOR = [0, 0, 0];   // Å, a fixed external point

registerTerm({
  id: "cg.tetherWell", engine: "cg",
  order: 5,                              // runs BEFORE every built-in term
  label: "Tether well — harmonic to a fixed point",
  unit: "kcal/mol",
  enabled: function () { return this.n > 0; },
  energy: function (pos, f) {
    let U = 0;
    const dx = ANCHOR[0] - pos[0], dy = ANCHOR[1] - pos[1], dz = ANCHOR[2] - pos[2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-12;
    const dr = r - 4.25;
    U = 0.5 * K * dr * dr;               // U = 0 exactly at r = 4.25
    const s = (K * dr) / r;
    f[0] += s * dx; f[1] += s * dy; f[2] += s * dz;
    return U;
  },
  report: function () { return 0; },     // no per-term accumulator; see below
});

const ff = new ForceField(sel, { rc: 10, gamma: 1.0 });
ff.compute(ff.ref);                      // the term is live on the very next call
```

`tests/test_physics_terms.js` asserts the analytic values for exactly this
shape: at `r = r0` the energy is `0` and the force buffer is bit-identical to a
field with no term registered; at `r = r0 + 1.75` the term contributes
`U = 3.0625` and `f₀ₓ = −3.5`, to 1e-12. It also asserts that `unregisterTerm`
restores energy **and every force component bit-exactly**, which is what an
analysis module that perturbs the model relies on.

## Order is the physics

`U` is a float64 sum, so the order terms are added in is part of the result.
`order` ascending is the addition order. Two consequences:

1. **Give the CG term an integer order on one side of 50.** The CG assembly has
   one hard-coded term — the inlined holo compression floor — that cannot be a
   descriptor, because it accumulates `U` pair-by-pair into the running total.
   `compute()` therefore walks the plan in two segments around it: `order < 50`
   runs before the floor, `order >= 50` after. There is no order-50 slot.
2. **Group what was grouped.** If the original code summed N terms in one `+=`,
   keep them in ONE descriptor that returns that same expression. Splitting
   `U += a + b` into two descriptors inserts an extra rounding. That is why
   `heavy.covalent` (bonds + coordination), `heavy.dihedrals` (impropers +
   propers) and `heavy.nonbonded` (lj + elec + gb + hbond) are one term each.

Moving a term to a different `order`, or splitting/merging one, is a **physics
change**: `tests/test_heavy_golden.js` (24 accumulators × 6 configs × 12 poses,
bit-exact) and `tests/test_cg_ligand_golden.js` will fail. That is the point.

## Performance

The registry is DATA, not a per-step closure pipeline, because the compute paths
are inlined typed-array loops at ~14 ms (heavy) and ~0.09 ms (CG) per call. The
plan is a frozen array built once per registry generation; per `compute()` the
cost is one integer compare against a module-local, plus one `.call` per term.
Terms allocate nothing.

Measured with interleaved A/B blocks (`bench/perf.js`, and a block-timed
9 × 4000-call CG benchmark because `bench/perf.js` times 30 individual calls and
the `performance.now()` pair dominates at 0.09 ms):

| | before | after | delta |
|---|---|---|---|
| Heavy, 1308 atoms (`bench/perf.js`, 8 pairs) | 13.76 – 14.16 ms | 13.66 – 14.16 ms | **−0.3 % (noise)** |
| CG, 164 beads (block-timed, 15 pairs) | 0.09065 ms mean | 0.09118 ms mean | **+0.58 %** (+0.53 µs) |

The heavy path is unchanged. The CG path is reproducibly ~0.6 % slower
(12 of 15 pairs positive). That cost is **not** the dispatch: three alternative
implementations were built and measured, and all three landed on the same
+0.55 µs —

1. a `switch` on an integer `kind` with nine *direct* monomorphic call sites
   (the usual "descriptors as data the assembly switches on" answer),
2. an explicit `energy(ff, pos, f)` signature instead of `.call`,
3. unfrozen descriptors and plan array.

So it is the cost of the loop that replaced the straight-line statement
sequence, not of the indirection, and the simplest design is the one shipped.
If you ever find dispatch *has* become the bottleneck, the lever is (1); it was
measured here and buys nothing.

## What this does NOT support

Stated plainly, because a "plugin system" that only supports a fixed set of
shapes is worth less than one that admits it:

- **One shape only.** A term is `(pos, f) → ΔU` plus an optional predicate. There
  is no hook for a term that needs its own scratch, its own pair list build, or
  its own per-topology setup phase.
- **One float64 per term.** A term that must publish several energies publishes
  them on the force field (as every built-in does: `ff.repU`, `ff.elecU`, …) and
  returns their pre-summed total.
- **No per-force-field registry.** The registry is process-global: register once,
  and every `ForceField` / `HeavyForceField` in the process picks it up. An
  analysis module must `unregisterTerm` in a `finally`.
- **No lifecycle hooks.** There is no `setup(ff)` / `teardown(ff)`. Term state
  must be derived from the force field.
- **No sandbox and no dependency injection.** A term closes over whatever it
  likes, including module-level state, exactly like any other module in `src/`.
- **Not every built-in term is a descriptor.** The CG holo compression floor is
  inlined in `compute.js` on purpose (see *Order is the physics*), and the heavy
  NaN/∞ guard and the stage-5a binding accounting are not terms either — they
  contribute no force and (for the guard) no energy.
- **No versioning or provenance.** Replacing a built-in by re-registering the
  same id is supported and silent; there is no record of who changed what.
- **`respa.js` does not use the registry.** The r-RESPA short-range path
  (`src/physics/integrators/respa.js`) still spells out its own kernel calls.
  A term registered here is not seen by the RESPA inner loop.
- **CG reports no per-term energies.** The CG assembly never published them and
  this change did not start: adding them would add nine own instance fields to
  `ForceField`, which `tests/test_cg_class_shape.js` pins to exactly 70. So
  `describeTerms("cg", ff)` returns `energy: null` for every CG term.
