# CG ↔ Heavy — Coarse-Grained ENM & All-Atom Handshake (C24/C26)

*Source anchors: `src/forcefield.js:1`, `src/ff-params.js:1`, `src/heavy.js:1`, `src/viewer.js:163`, `src/pdb.js:74`.*

This document explains the two force-field resolutions, their single-responsibility
handshake, and the double-coverage guards (notably `nativeContacts`) that keep
`compute()` allocation-free and formally correct.

---

## 1. CG ENM — Cα Elastic Network (Tirion 1996)

**Beads:** one particle per residue at Cα (`src/pdb.js:83` `parseCa` → `beads[]`,
chaingroups `selectSystem` at `src/pdb.js:149`, segments split at >4.8 Å).

**Potential** (`src/forcefield.js:1` header):
```
U = Σ ½ k_b (r−r0)²          k_b=100 kcal/mol/Å², r0≈3.81 Å (AMBER CA)
  + Σ ½ k_θ (θ−θ0)²          k_θ=20  kcal/mol/rad² (Tirion 1996; see src/ff-params.js:10)
  + Σ ½ γ (r−r0)²·H(Rc−r0)   γ=1.0, Rc=10 Å (default), |i−j|>2, protein-only
  + WCA repulsion            ε=0.3, σ=4.0 Å, r_e=2^(1/6)σ≈5.61 Å, smooth at re
  + cross protein–ligand LJ + screened Coulomb (Lorentz–Berthelot) + H-bond + EEF1-lite burial
  + holo springs ½γ_lig(r−r0)² γ_lig=0.5, r0≤6 Å (holo-pose pin)
  + nativeContacts Gö springs ½·1.0·(r−r0)² for 1-5+ pairs with r0<re (see §4)
```

*Backbone constants provenance:* Boltzmann inversion `k = k_B T/σ²` on AMBER
ff14SB explicit-solvent MD (`src/ff-params.js:10`):
σ_b≈0.045 Å ⇒ k_b≈120→100 (stability, `src/integrator.js:145` dt=4 fs);
σ_θ≈0.14 rad ⇒ k_θ≈30→20. Validated `tests/test_bond_dist.js`: 200-step
Langevin on 1crn keeps ⟨b⟩=3.81±0.05 Å.

**ENM cutoff scan:** `notebooks/cutoff_scan.md` — default `Rc=10 Å` is near
optimal for B-factor Pearson on 4W52 (T4 lysozyme L99A). Springs count is
monotonic in Rc (verified `tests/test_cutoff.js` at 8/10/12 Å).

**Sequence-dependent ENM (Bahar-style, opt-in):**
```js
w_i = SEQ_WEIGHT[RES_CLASS_OF[bead.resName]]   // src/ff-params.js:SEQ_WEIGHT
// H:1.0  A:0.9  P:1.1  Cp:1.05  Cn:1.05  (Miyazawa–Jernigan ±10%)
K_seq(i,j) = gamma * (1 + 0.2*(w_i + w_j)/2)  // src/forcefield.js:≈235, applySeqWeights()
```
Uniform `springK.fill(gamma)` is default; call `ff.applySeqWeights(beads)` or
`new ForceField(sel,{seqWeight:true})`. Tested `tests/test_enm_seq.js` on
1ubq (ubiquitin, 76 aa): uniform vs seq Pearson both finite and seq not
worse (`ΔR ≥ −0.05`), with at least one `K_seq ≠ gamma`.

**Residue classes:** `src/ff-params.js:60` `RES_CLASS` (H/A/P/Cp/Cn) + sigma/eps/q;
ligand elements `LIG_ELEMENT` (C,N,O,S,F,CL,BR,I,P) with ΔG burial.

---

## 2. Heavy (All-Atom) — `src/heavy.js`

**Beads:** all heavy atoms (protein `ATOM` + `HETATM` hetero/metal/ligand,
solvent `HOH` dropped) via `parseHeavy` at `src/heavy.js:82`; `selectHeavy`
at `src/heavy.js:165` partitions `nProt / nHetero / nLigAtoms`.

**Potential:** covalent bonds (`buildTopology` CSD radii `src/ff-params.js:114`
`COVALENT_RADIUS` + `BOND_SLACK=1.15` hard cap 2.2 Å — disulfide S–S 2.04 Å
passes, Ca–N 2.9 Å rejected, validated `tests/test_topology.js` on 1crn),
angles, proper/improper dihedrals (analytic), metal coordination springs
(`METAL_K=40`, `METAL_ELEMENT` coordR/N), SASA burial (`SasaModel`), GB
(`GeneralizedBorn`) + screened Coulomb + directional H-bond (`DirectionalHBond`),
spatial-grid non-bonded (`SpatialGrid` RCUT=8.5 Å).

**Disulfides & PTMs (C28):** `src/heavy.js:248` covalent scan captures
`CYS SG–SG 2.03–2.05 Å` (CSD S 1.02×2=2.04 Å, 2.04×1.15=2.35→capped 2.20 still
passes with 0.16 Å margin). PTMs appear as explicit heavy atoms (e.g.
SEP phosphoserine). CG mode intentionally coalesces them to Cα; see
`src/pdb.js:74` comment. Test coverage: `tests/test_topology.js`.

---

## 3. Handshake — Switching Resolutions

| Aspect | CG (`ForceField`) | Heavy (`HeavyForceField`) |
|--------|-------------------|---------------------------|
| Parser | `parseCa` + `selectSystem` | `parseHeavy` + `selectHeavy` (+`appendHeavyLigands`) |
| FF ctor | `new ForceField(sel, par, ligands)` | `new HeavyForceField(system, par, ligands)` |
| Mass | 110 Da/Cα, ligand united-atom (12–127 Da) | `heavyMass` IUPAC/AMBER per element |
| Reference | `ff.ref` Float64Array(3n) native Cα | `ff.ref` Float64Array(3n) heavy coords |
| Units | `src/units.js:22` KCONV=418.4, KB_KCAL=0.001987 | identical (shared leaf) |
| Viewer | `viewer.setSystem(sel, ff)` → `viewer.center` protein-only | same API, same protein-only center |

Conversion is *one-shot* at `main.js:Build System` — no live inter-conversion.
Ligand MOL2/HETATM placement is concatenated after protein beads in both
modes, preserving `nProt` as the cut index. `ff.nProt` / `ff.n` / `ff.isLigand`
are the handshake invariants inspected by `src/analysis.js` and `src/viewer.js`.

---

## 4. `nativeContacts` Double-Coverage Guard (C24)

*Location:* `src/forcefield.js:340` (`this.nativeContacts`).

**Problem:** folded ligands/sugars place non-bonded 1-5+ pairs at 2.9–4.0 Å,
inside their own LJ `r_e = 2^(1/6)σ` (≈3.8–5.0 Å). The grid's WCA repulsion
would see `U≈745 kcal/mol` strain at the native geometry and launch united
atoms as projectiles.

**Fix (Gö-like):** at construction, every *non-excluded* intra-molecule pair
`(i,j)` with native distance `r0 < r_e` is *excluded from the grid* and given
a soft harmonic spring `½·1.0·(r−r0)²` (`nativeContacts`, same `[i,j,r0]` layout
as `holoSprings`). Pair is then counted exactly once.

```
native r0 ──► [already in _excluded? (bonds/angles/springs)] ──► keep excluded
     │                  │no
     │          r0 < re =2^(1/6)·(σi+σj)/2 ? ──► yes → nativeContacts spring + _excluded.add
     │                  │no → leave to WCA grid (repulsion zero at re, smooth)
     └─ inter-molecule P–L ? ──► skip (P–L already holo-pinned or binding-evaluated)
```

*Why O(N²) only at build:* the scan is construction-time; `compute()` then
uses the grid + springs without re-scanning. Forces are still `O(N)` average
(cell lists). Diagram is also inline at `src/forcefield.js:340`.

*Verification:* `grep -n "nativeContacts" src/forcefield.js` shows construction
and energy application (`src/forcefield.js:425` `_harmonicPairs` with k=1.0).

---

## 5. Handshake Center / Radius Fix (C26)

**Bug it fixed:** previously `viewer.center` / `viewer.radius` were computed over
*all* particles `n = nProt + nLigAtoms`. A single benzene at 20 Å from the
protein COM would shift `center` by ~0.3 Å per ligand and inflate `radius`
(`maxR²`) so the protein appeared small and off-screen.

**Fix:** `src/viewer.js:163`–`182`:
```js
// center: Σ_{i < nProt} r_i / nProt   (protein-only)
// radius: max_{i < nProt} |r_i - center| *1.2, floor 12 Å
```
Heavy mode inherits the same path (`src/viewer.js:99` `setSystem(sel,ff)` uses
`ff.nProt` / `ff.ref` protein prefix regardless of `sel.heavy`). Documented
`src/heavy.js:1` header handshake note.

*Center/radius are viewer-only* (no physics); fixing them avoids a subtle
multi-scale coupling where CG and heavy would otherwise show different zoom for
the same PDB+ligand.

*Related:* `src/heavy.js` itself does **not** compute a viewer radius — it
only provides `ref`; the protein-only rule lives in the single place
`src/viewer.js`, so both modes stay consistent (see this doc §3 table).

---

## 6. References & Measurable Checks

```bash
grep -n "SEQ_WEIGHT\|w_i" src/forcefield.js src/ff-params.js  # C21
node tests/test_enm_seq.js                                    # C21 PASS
ls notebooks/cutoff_scan.md                                     # C23
grep -n "nativeContacts\|CG_HEAVY" src/forcefield.js docs/CG_HEAVY.md  # C24
node tests/test_bond_dist.js                                    # C22 3.81±0.05
node tests/test_cutoff.js                                       # C23 monotonic
node tests/test_topology.js                                     # C28 disulfides
```

*See `notebooks/cutoff_scan.md`, `tests/test_enm_seq.js:1`, `tests/test_bond_dist.js:1`,
`src/forcefield.js:340` for line-cited provenance.*
