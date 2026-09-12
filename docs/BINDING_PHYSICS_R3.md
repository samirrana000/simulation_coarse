# Phase R3 — Heavy-Atom Weak-Interaction Physics (design + prototype)

Loop-1, step 3. Follows `BINDING_PHYSICS_R1.md` (ranks 6–8) and `BINDING_PHYSICS_R2.md`
(CG non-goals: no ring planes at Cα). All terms below are HEAVY-MODE ONLY.
Prototype: `/tmp/opencode/r3_heavyproto.mjs` (zero-dep, gradients FD-verified).

## 0. Current heavy-mode gaps (verified)

- Nonbonded = LJ + GB/Coulomb + H-bond only (`src/heavy.js:947-1015`).
- H-bond angle hard-coded 180°, no angular gradient (`src/physics/hbond.js:91,115-120`).
- Metals = distance springs k=40 (`src/heavy.js:383-404`); `src/chem/metals.js`
  `enforceCoordination` (cross-angle harmonics) exists but is NOT wired into heavy.
- Salt bridges: **already emergent** — AMBER united-atom charges live (Asp/Glu −1,
  Arg/Lys +1, `src/physics/charges.js:24-26`); no mean-neutralization zeroing found in
  heavy path (GB neutralizes net charge via local dipole form). No fix needed.
- GAFF2-lite types (ca/nh/nb/cl/br/i) + protonation states exist from Phase 2 —
  the ring/halogen detection hooks below build directly on them.

## 1. Term specifications

### (a) π–π stacking — anisotropic Gaussian ring-frame potential
Frames: ring centroid + unit normal from aromatic rings (typed `ca/nh/nb/c2`; ≥5 ring
atoms). One form captures parallel + displaced-parallel + (partially) T-shaped:

```
ρ  = (|r·n_A| + |r·n_B|)/2                      interplanar gap
t² = |r|² − (r·n̄)²                             in-plane offset² (n̄ = mean normal)
U_att = −ε · G(ρ; ρ0=3.8, w=0.45) · exp(−t²/6) · (½ + ½(n_A·n_B)²)
U_rep = +4ε · G(ρ; ρr=3.0, w=0.35) · (½ + ½(n_A·n_B)²)
U = U_att + U_rep        (ε=2.4 kcal/mol; centroid cutoff 5.5 Å)
```

- Well depth ≈ −1.7 at 3.8 Å parallel; orientation factor 1.0 parallel / 0.5 orthogonal;
  in-plane damping halves attraction at 2.45 Å offset (Hunter–Sanders geometry).
- Gradients: full chain rule incl. ∂ρ/∂n (gap depends on the normals), ∂t²/∂n̄, unit-sphere
  projection for normal gradients. **FD-verified max rel err 6.2e-11**.
- Cost: ~120 FLOPs/ring-pair; ring list built once per topology (aromatic rings ≪ N atoms).

### (b) Cation–π — cos²-gated Gaussian along ring normal
```
U = −ε_cπ · G(r; r0=4.3, w=1.1) · max(0, cos α)²,  ε_cπ=3.5
α = angle(cation→centroid, ring normal)
```
- Cations: Lys NZ, Arg (guanidinium C ζ-center), protonated-His (both N), any +1-typed
  ligand N (GAFF2 `n4`/`n3` + positive charge).
- On-axis 4.3 Å → −3.5 (lit 2–5, Gallivan & Dougherty 1999); off-axis 2.8 Å lateral → −2.0.
  NOTE: cos² gate is soft; if sharper discrimination wanted, switch to cos⁴ (doc'd choice).
- Gradients FD-verified exact (0.0 rel err). Cost ~40 FLOPs/pair.

### (c) Halogen σ-hole — angle-gated Gaussian (VinaXB lineage, Koebel et al. 2016)
```
U = −ε_X · G(r_XD; 3.1, w=0.55) · max(0, cos β)²
β = angle(C–X···D), σ-hole axis = extension of C→X beyond X; X ∈ {Cl, Br, I}, F EXCLUDED
ε_X: Cl 1.2 / Br 2.0 / I 2.5 (scaled from prototype ε=2.0 demo)
```
- Linear C–Cl···O → −2.0; 160° → −1.77; 120° → −0.50 (lit: 1–4 net in water, 3–5 buried).
- Acceptors D: backbone/Asn/Gln carbonyl O, Ser/Thr/Asp/Glu O, sulfur S, aromatic π.
- Gradients use the parallel-transport-safe forms:
  ∂cosβ/∂D = (â − cosβ·r̂)/d, ∂cosβ/∂X = −(â−cosβ·r̂)/d + (r̂−cosβ·â)/l_CX,
  ∂cosβ/∂C = −(r̂−cosβ·â)/l_CX — **FD-verified 4.6e-10** (careful: the max(0,·)² kink at
  cosβ=0 is non-smooth; FD must avoid exactly-perpendicular geometries).
- Cost ~60 FLOPs/pair; pairs pre-filtered to halogen atoms only (tiny set).

### (d) Salt bridge — VERDICT: no new term needed
Live ±1 charges + GB already produce it; R3 validation: Lys35 NZ–HEPES SO₃ pair
scores −0.76 kcal/mol via existing Coulomb path in the CG prototype (R2) and the heavy
GB path carries the same charges. Keep; nothing to implement.

### (e) Metal coordination upgrade — wiring swap
Replace heavy.js metal distance springs (`:383-404`) with `src/chem/metals.js`
`enforceCoordination()` (radial k=40 + cross-angle k=20, ideal angles per geometry class).
The function already exists, shares `METAL_ELEMENT` radii, accumulates into caller forces.
Swap cost: ~30 lines in heavy.js force assembly. Prevents coordination-sphere
dissociation that pure distance springs allow (octahedral → square-planar drift).

## 2. Double-counting analysis (each term vs GB/LJ)

| Term | Verdict | Rationale |
|---|---|---|
| π-stack | **supplement** (dispersion/CT GB can't represent) | aromatic LJ ε (GAFF2 ca ~0.086) badly underbinds stacks; anisotropic term adds the missing −1..−3. Reduce ring-atom pair LJ ε by ~30% when both in aromatic rings to avoid double attraction. |
| cation-π | **supplement, charge-gated OFF when GB already binds** | GB point-charge captures only ~30% (quadrupole missing). Keep full ε_cπ but subtract linear-charge component: ε_eff = ε_cπ − 0.5·q_GB(r) contribution. Simplest robust choice: ε_cπ·(1−|q_cat|/2). |
| halogen | **supplement** | σ-hole anisotropy is absent from isotropic LJ entirely; LJ Cl···O well (~0.2) ≪ σ-hole (1–4). No subtraction needed. |
| salt bridge | none | fully emergent from GB + charges. |
| metal cross-angles | **replaces** nothing | new angular constraints; radial part already exists (drop k=40 springs when enforceCoordination on to avoid double radial). |

## 3. Cost table (per pair, FLOPs, vs existing heavy nonbonded ~200/pair)

| Term | FLOPs | Pairs | Notes |
|---|---|---|---|
| π-stack | ~120 | rings×rings + ring×aromatic-lig (~10-50) | ring cache per topology |
| cation-π | ~40 | cations×rings (~5-20) | cation list from charges |
| halogen | ~60 | X×acceptors (~5-30) | halogens rare |
| metal upgrade | ~0 net | ≤ metals×coordN | swap, not addition |
| **Total added** | **< 5% of nonbonded time** | | heaviest system 4HHB ~15 rings |

## 4. Exclusion bookkeeping
All new terms run INSIDE the existing nonbonded pair loop cutoff (9 Å) with the same
`_excluded` set (heavy.js:568-583): ring internal pairs excluded automatically
(1-2/1-3); π-stack cross-terms need no exclusion (rings from different molecules or
non-bonded rings only — same-molecule adjacent rings (Phe-Phe) ARE in the nonbonded
set and SHOULD stack (arithmetically fine). 1-4 scale 0.5 applies via existing map.

## 5. Prototype numbers (see /tmp/opencode/r3_heavyproto.mjs)

| Case | Energy (kcal/mol) | Literature |
|---|---|---|
| π parallel-displaced (3.6 Å, 1.5 Å off) | +0.71 (near crossover; well at 3.8) | −1..−3 |
| π parallel scan | min −1.70 @ 3.8 Å | −1..−3 ✓ |
| π T-shape | +2.88 (partial clash at 4.0 Å) | −1..−2 (real T at 5.0 Å scores −0.5) |
| cation-π on-axis 4.3 Å | −3.50 | −2..−5 ✓ |
| cation-π off-axis 2.8 Å | −2.03 | gate softness noted |
| halogen linear | −2.00 | −1..−4 ✓ |
| halogen 160° / 120° | −1.77 / −0.50 | angular decay ✓ |
| Gradient FD max rel err | π 6.2e-11 · cation-π 0 · halogen 4.6e-10 | all < 1e-6 ✓ |

## 6. Loop-2 heavy implementation order — STATUS (Loop-2 S3, implemented)

1. **π-stack — DONE** (`src/physics/weakint.js` `piStackEnergy`/`piStackForces`;
   heavy.js opt-in `par.weak="on"`). Ring cache via `buildRingFrames` (protein
   name sets + GAFF2 sp2 types + geometric fallback, fused-ring subsumption);
   atom-level gradient chain adds the Newell-normal lever
   `F_k = (g_N × d_k)/|N|` (FD 1.2e-9) — an extension beyond the frame-level
   prototype gradients. LJ ε-downscale (§2) deliberately NOT applied yet:
   per-term accumulators (S4) land before any double-count guard is tuned.
2. **cation-π — DONE** (`cationPiEnergy`/`cationPiForces`). Cation list: Lys NZ,
   Arg CZ, HIP N (via protonation-state names), charged ligand N (GAFF n3/n4).
   Charge-gating ε·(1−|q|/2) (§2) NOT yet applied (same S4 sequencing).
3. **halogen σ-hole — DONE** (`halogenEnergy`/`halogenForces`, ε per element
   Cl 1.2 / Br 2.0 / I 2.5, F excluded; parallel-transport-safe gradients).
   Halogen list: Cl/Br/I with bonded C; acceptors from the H-bond acceptor set.
4. **metal `enforceCoordination` swap — DONE** (opt-in `par.metalAngles=true`,
   default false). Construction: detectCoordination + classifyGeometry; those
   metals' k=40 springs removed (no double radial); compute calls
   `enforceCoordination(..., {ideal:true})`. Metals without geometry keep springs.
5. H-bond angle fix — NOT in S3 scope (R3 §6 item 5, later step).

Measured (Loop-2 S3, tests/test_weakint.js 49/49; test_all 32/32):
FD π 1.8e-10 frame / 1.2e-9 atom-level · cation-π 4.7e-10 · halogen 4.6e-10;
energies π −1.70 @3.8 (min −2.03 @3.95) · cπ −3.50 on-axis · XB −ε_X linear;
4W52 rings 5 Phe + 6 Tyr + 6 Trp(2×3 fused) + 1 His + 1 benzene = 19;
default-off bit-identity exact (U = −2077.8811924685633 pre/post);
weak-on ΔU = −2.38 (π +3.68 T-clash contacts — prototype behavior, §5 caveat;
cπ −6.06; XB 0 — no Cl/Br/I in 4W52); perf 14.87 ms/step off vs 14.87 on
(ratio 1.000, budget 1.1×); NVT 100 steps RMSD 0.35 Å (4W52+benzene) /
0.94 Å (4HHB 4 Fe hemes, metalAngles exercised).

Non-goals (unchanged from R1/R2): explicit waters, PME, polarization/charge-transfer,
covalent inhibitors.

## 7. References

- Hunter, Sanders 1990 (π-system geometry); Grimme 2011 (stacking dispersion).
- Gallivan, Dougherty 1999 PNAS (cation-π energies/geometries).
- Koebel et al. 2016 J. Cheminf. (VinaXB σ-hole functional form).
- Genheden & Ryde 2015 (MM/GBSA error discipline for validating supplements).
- Wang, Hou, Xu 2007 (halogenated-ligand affinity gaps from missing σ-hole physics).
