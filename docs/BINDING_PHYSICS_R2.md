# R2 — CG Binding-Physics Upgrade: Virtual Sites + Term Set (design + prototype)

Step 2/7 of the binding-physics program, Loop 1. **Design document with
offline-validated prototype numbers.** *Loop-2 status: step 1 (term b charge
table) has landed in `src/` as a default-off opt-in — see the status note in
§2(b) and §4 item 1; everything else here is still design-only.* Builds on
`docs/BINDING_PHYSICS_R1.md` (interaction priority table, energy scales, code
citations; all `src:` line references below verified against current code).
Prototype: `/tmp/opencode/r2_cgproto.mjs` (zero-dep Node; may import repo
modules read-only), raw output `/tmp/opencode/r2_proto_output.txt`.

---

## 0. Grounding recap (from R1, verified)

| Fact | Code | Implication for R2 |
|---|---|---|
| CG H-bond isotropic Gaussian −0.8 kcal/mol @3.2 Å, flags only, no angle | `src/ff-binding.js:57,127-133` | Replace with directional term (a) |
| Protein bead charges all q=0 | `src/ff-params.js:61-67` | Coulomb path (`ff-binding.js:114-124`) is coded but dead — revive with charge table (term b) **[DONE Loop-2 S1: CG_FORMAL_CHARGES, opt-in `binding.charges`]** |
| H-bond flag = residue class P/Cp/Cn only | `src/forcefield.js:185` | **52% of backbone donors/acceptors score ZERO today (prototype §C)** — every residue gets backbone virtual sites |
| Holo springs k=0.5 hold native pose | `src/forcefield.js:116,276-303` | Keep; binding terms must not double-count (exclusion set already handles) |
| EEF1-lite ΔG_a −0.25…−0.55, B_a=1−exp(−n/3) | `src/ff-binding.js:89-155` | Term (c) multiplies this — no new loop |
| Cross cutoff 9.0 Å, switch at 7.65 Å | `src/forcefield.js:115`, `ff-binding.js:25` | All new terms live inside existing pair loop |

The three native 4w52 HEPES contacts used throughout as the test case (X-ray,
all < 3.0 Å, textbook geometries):

- **N4⁺···O=Phe104 2.73 Å** — backbone carbonyl acceptor; Phe104 is class A ⇒
  scores **0.00** today.
- **O8···H–N–Leu32 2.96 Å** — backbone NH donor; Leu32 class H ⇒ **0.00** today.
- **O1S⁻···H₃N⁺–Lys35 2.66 Å** — sidechain salt-bridge; Cp flag gives −0.53
  isotropic, Coulomb 0.

---

## 1. Virtual-site scheme: how much direction information survives Cα-only?

### 1a. Method and lineage

Per-residue local orthonormal frame from the Cα triplet (Levitt 1976 lineage;
same frame used by BBQ Gront et al. 2007 and PD2 Maupetit et al. 2006):

```
u = normalize(Cα(i-1) − Cα(i));  w = normalize(Cα(i+1) − Cα(i))
b = normalize(u + w)              bisector, in-chain
n = normalize(u × w)              ≈ peptide-plane normal
t = b × n                         third axis
```

Virtual sites = fixed offset β (in b,n,t basis) per **secondary-structure
class** per site type. SS from Cα only (DSSP-lite, ~5 FLOPs): helix if
d(Cα_i, Cα_{i+3}) < 5.8 Å; strand if bend angle > 125°; else coil.
Rebuild per force call from current bead positions: ~30 FLOPs/residue —
O(N), amortized across all ligand pairs.

Literature anchors for reconstruction fidelity: **BBQ (Gront, Kulp et al.,
Nucleic Acids Res. 2007)** 0.6–0.8 Å all-backbone RMSD; **PD2 (Maupetit,
Tuffery et al., Nucleic Acids Res. 2006)** 0.5–1.0 Å; **PULCHRA (Rohl et al.,
J. Chem. Theory Comput. 2004)** similar; original **Levitt, J. Mol. Biol.
1976** first Cα→backbone reconstruction. Our 3-parameter per-SS fits land at
the good end of that range in regular SS (below), with zero iteration cost —
no fragment threading, just 3 basis vectors.

### 1b. Measured fidelity (prototype §A; 4w52 + 1crn + 1ubq, 50/50 CV by parity within SS class)

**4w52 (main testbed):**

| SS | n | O-site RMSD | N-site RMSD | C=O-axis error (med / p90) | NH-axis (proxy) error |
|---|---|---|---|---|---|
| Helix | 102 | 0.93 Å | 0.16 Å | **6.1°** / 54° | ≤ C=O (NH is tighter) |
| Strand | 10 | 0.47 Å | 0.08 Å | **14.8°** / 21° | " |
| Coil | 50 | 1.39 Å | 0.34 Å | **40.8°** / 100° | " |

Representative regression coefficients (b,n,t basis): helix β_O=[0.56,−1.10,−1.84],
β_CO=[−0.24,−0.77,−0.59], β_N=[0.87,0.38,1.10]. 1crn/1ubq reproduce the
pattern (H 1.06/1.37 Å site, 13.7/23.8° axis — small-n classes noisier; coil
~35–60° everywhere). **Angular error vs true H-bond geometry: 6–15° in
regular SS (helix/strand = ~70% of residues), 35–60° in coils.** The N-site is
2–4× tighter than O-site (amide N is more constrained by the Cα frame than
the carbonyl O; both match BBQ/PD2 benchmarks).

Two validations of axis proxies (the part that matters for gating):
- **Frame normal tracks the true peptide plane to 19° median** (N,Cα,C plane)
  — sufficient for an out-of-plane soft penalty, not a hard gate.
- **Cα-only NH direction proxy** (anti-bisector at virtual N between the two
  flanking Cα directions) vs the actual native O8···H–N(Leu32) approach:
  **cos 0.95** (prototype §B [2]) — backbone donor direction is *knowable*
  from Cα alone at the required precision.

### 1c. Sidechain sites: statistical placement (per residue type)

Rotamer statistics (4w52, per type; Dunbrack & Karplus 1993 library
multimodality is the cause):

| Residue | n | tip dist (mean±rms) | site RMSD vs type-mean | direction spread (med) |
|---|---|---|---|---|
| THR | 11 | 2.4±0.0 Å | 2.07 Å | 79° |
| ASP | 10 | 2.9±0.3 Å | 3.19 Å | 88° |
| LYS | 9 | 5.9±0.4 Å | 5.16 Å | 80° |
| GLU | 8 | 4.0±0.6 Å | 3.55 Å | 52° |
| ASN | 11 | 3.0±0.2 Å | 3.00 Å | 65° |
| (SER/ARG/GLN/TYR/TRP similar: 2.1–6.7 Å site noise, 48–85° spread) | | | | |

**Verdict: sidechain polar tips CANNOT be sharp-gated.** Site placement noise
0.5–5 Å + 50–90° direction spread ≫ H-bond geometry tolerance (±0.3 Å,
±30°). The virtual-NZ experiment makes it concrete: type-mean Lys NZ lands
2.42 Å from crystal NZ, and single-axis cos² would kill the *native* contact
(cos 0.15 → cos² 0.02). Design rule from the data:

- **Backbone CO/NH sites** (per-SS regressed): real virtual sites with
  directional gates. These carry the specificity.
- **Sidechain charged/polar residues**: ±1 **formal charge on the bead**
  (Coulomb term b) + **soft cone gates (≥70° half-width) only** if a virtual
  donor/acceptor axis is wanted; never sharp cos² on rotamer-averaged axes.
  Matches R1 §4 rank-4 ("flip q table + GB-desolvation does the rest").
- Short sidechains (Ser/Thr/Asp/Asn, tip ≤3 Å) are borderline (2 Å site noise)
  — treat as bead-charge + soft-cone; long (Lys/Arg/Glu/Gln) charge-only.
- Cβ vs −bisector: 49° median error ⇒ naive bisector is not a sidechain
  axis; use per-type regressed vectors where used at all.

### 1d. Concrete scheme (what Loop 2 implements)

Per protein residue i, rebuilt each force call from the Cα trace:

```
frame(b,n,t) at Cα_i  (skip termini & chain breaks; skip if |u+w|<0.05)
sites:
  O_acc(i):  pos = Cα_i + β_O[SS]        axis = β_CO[SS] (C→O direction)
  N_don(i):  pos = Cα_i + β_N[SS]        axis = anti-bisector proxy (§1b)
  q(i):      +1 Lys/Arg (His per protonation heuristic, src/chem/protonation.js
             already exists for heavy) / −1 Asp/Glu — on the BEAD
  (optional, soft) sidechain vector for Ser/Thr/Asp/Asn tips
```

Arrays: 3×3 floats + 2 charges per residue = ~15 Float64/residue, allocated
once in ForceField, recomputed in the binding pass preamble (O(nProt)). Both
virtual sites participate in the existing neighbor-grid pass: H-bond pair
loop uses site positions and bead→site mapping is fixed (no new grid cells
needed — sites sit within 3 Å of their bead, well inside the 9 Å cutoff).

---

## 2. CG+ binding term set — equations, parameters, costs, expected gains

Standing notation: r = distance between interaction sites; all terms are
evaluated **inside the existing cross pair loop** (`ff-binding.js` PASS 1)
for pairs passing the 9 Å cell scan and smooth switch `sw(r)`, excluded from
holo-spring pairs (no double counting, existing `_excluded` set).

### (a) Directional H-bond (backbone virtual sites) — RANK 1

```
E_hb = −ε_hb · G(r) · A_D · A_A
G(r)  = exp(−(r − r₀)²/(2σ²))                r₀ = 3.0 Å, σ = 0.5 Å
A_D   = max(0, cosθ_D)²                      sp2 donors (backbone NH)
A_D   = cone(θ_D; 110°)                      sp3 donors (Lys/Arg/amines, OH)
A_A   = max(0, cosα)                          Baker–Hubbard hemisphere:
        cosα = (donor − acc) · C=O-axis       at-O 90–180° (McDonald–Thornton)
cone(x;θc) = clip((cosθ − cosθc)/(1 − cosθc), 0, 1)
```

Parameters (from R1 §1a/§2 #1 + calibration anchors):
- ε_hb = **2.0 kcal/mol** (HYDE/Chemgauss net-HB window 1–3; Böhm LUDI fit
  −1.1 was pre-desolvation-balanced; today's 0.8 underweights buried H-bonds
  ~2–3× per R1 §5 rank-1).
- r₀ = 3.0, σ = 0.5 (D···A 2.7–3.2 survey range; Vina uses r₀≈3 with two
  Gaussians, single Gaussian suffices per Vinardo ablation).
- **Valence caps**: each backbone O accepts ≤2 (bifurcated allowed), each NH
  donates 1, sidechain donors capped at 2 for amines — running per-site
  counters in the pair pass, greedy by |E| (Chemgauss4/HYDE rule). Cost: one
  increment + compare per accepted pair.
- Coil residues: widen A_A to the max(0,cosα) plus 0.5·(1−oop²) in-plane
  fallback (axis error 40°+ there); helix/strand keep the sharp form.

**Cost**: +2 dot3 + 2 max + 3 mul + 1 exp (vs isotropic's 1 exp) per flagged
pair ⇒ measured **21.1 vs 15.7 ns/pair = 1.3×** on the tiny H-bond-flagged
subset (prototype §F); frames ~30 FLOPs/residue amortized per call. Net frame
impact estimated <5% (H-bond pairs are ≪1% of all cross pairs).

**Expected ranking gain** (R1 §3b evidence): H-bond carries **26–28% of
Vina/Vinardo score mass** (Quiroga & Villarreal 2016 audit); isotropic
versions are the known false-positive source — Glide-XP's H-bond motifs +
HYDE's desolvation balance are what separate winners. The prototype shows the
mechanism: isotropic scores all orientation-decoys identically; directional
kills carbon-side approaches to 0.000 kcal/mol at identical distance.

### (b) Salt bridges via assigned formal charges — RANK 2 (cheapest Δ per line)

```
E_coul = 332.0637 · q_i · q_a / (ε(r) · r) · sw(r)     [ALREADY CODED]
ε(r) = 4 + 76·tanh(r/8)
q: Lys/Arg +1, Asp/Glu −1, His ±per protonation heuristic on bead
```

**Zero new math** — the term is live in `ff-binding.js:114-124`, dead only
because `RES_CLASS` has q:0 for all five classes (`ff-params.js:61-67`).
Revival = set q in that table + flow through `forcefield.js:184`. Prototype §D:
−2.2 kcal/mol at 2.66 Å (ε=28), −1.06 at 6 Å, −0.56 at the 9 Å edge — a sane
screened profile. Desolvation counterweight already exists (ΔG_a burial),
which is the Hendsch & Tidor lesson (bare Coulomb over-praises salt bridges;
burial penalty rescues it). Ligand side already carries q (N −0.30, O −0.50,
P +0.40 — `ff-params.js:76-87`), so ionic ligands light up immediately.

> **Loop-2 S1 status: DONE (2026-09-12, default-off opt-in).** Implemented
> as `CG_FORMAL_CHARGES` in `src/ff-params.js` (ASP/GLU −1, LYS/ARG +1;
> HIS deliberately absent ⇒ neutral 0 — HIP +1 hookup to
> `chem/protonation.js` deferred) wired into `forcefield.js` `_protQ` behind
> `par.binding.charges === true` (DEFAULT OFF — `RES_CLASS.q` stays 0, so
> pre-S1 trajectories are bit-identical; rollback = flag flip). Per the R2
> spec, formal charges are kept as-is — **no mean-neutralization** (ε(r) +
> the existing EEF1-lite burial/desolvation counterweight handle net-charge
> effects, the Hendsch–Tidor balance). Measured through the real CG pipeline
> (4W52 + ligandLib acetate with O⁻ at the native O1S⁻···Lys35 geometry,
> tests/test_charges.js): **ΔU_coul = −2.30 kcal/mol at 2.66 Å** (ε=28.4;
> prototype predicted −2.20), screened profile −2.30 → −1.42 @4 Å →
> −0.84 @6 Å → −0.03 @8 Å (sw cutoff 9 Å), analytic closed-form match
> <1e-9 at every r. Default path bit-identical to pre-S1 over 10 steps
> (energy/posHash/velHash exact); `test_all` 32/32, golden/negative/bindlog
> (18)/bindviz (21)/R1-validation/wikiskill-gate all PASS; CG default
> 0.24 ms/step unchanged (charges-on 0.23 — the `q1!==0` branch is ~free).
> Charged census on 4W52: 18 Asp/Glu, 26 Lys/Arg beads.

**Cost**: ~0 (condition `q1!==0` already in the loop; table edit only).
**Gain**: unblocks an entire dead interaction channel; charged-ligand
rankings (phosphate/sulfonate/carboxylate/ammonium drugs — ~40% of druglike
molecules) currently scored by LJ alone.

### (c) Hydrophobic-enclosure multiplier on existing burial — RANK 3

```
E_burial = Σ_a ΔG_a · B_a · (1 + λ_enc · enc_a)
enc_a   = clip((n_H(a) − 4)/8, 0, 1)      n_H = H/A-class protein beads within 7 Å
B_a     = 1 − exp(−n_a/3)                  [existing, ff-binding.js:150-155]
λ_enc   = 1.0
```

Glide-XP's signature super-linearity (Friesner et al., J. Med. Chem. 2006:
lipophilic atoms enclosed on multiple sides rewarded beyond linear SASA).
Prototype §E: benzene in the L99A cavity sits amid **7/7 H-class Cα** (enc =
0.38 ⇒ ×1.38 multiplier) vs a solvent-exposed control placement with enc = 0
(×1.00). The density accumulation `dens[a]` is already computed in PASS 1 —
n_H(a) needs only a second class-masked accumulation into the same loop
(+1 compare+add per pair), then the multiplier applies in PASS 2 (zero new
loop).

**Cost**: ~1 add+cmp per pair in PASS 1 + 1 mul per ligand atom in PASS 2.
**Gain**: Glide-XP's motif with near-zero code; separates buried-cavity
poses from surface-snorkeling decoys that burial alone scores equally.

### (d) Halogen-lite: C–X···O angle gate — RANK 4

```
E_xb = −ε_X · G_x(r) · max(0, cosθ_XB)²        X ∈ {Cl, Br, I}, F excluded
G_x  = exp(−(r − r_X)²/(2·0.3²)),  r_X = 3.1 Å (X···O)
θ_XB = angle(C–X vector, X···acceptor vector); ideal 180°, gate >120°
ε_X  = 0.6 (Cl) / 1.2 (Br) / 2.0 (I)  [net; R1 §2 #5: 1–4 kcal/mol gas]
acceptors: backbone O virtual sites (term a machinery) + ligand O/N
```

Ligand side already knows C–X bonds (united-atom internal FF, `ligand.js`);
the C–X direction is a rigid-body property cached per placement. Protein side
uses term (a)'s virtual O. VinaXB precedent (Koebel et al., J. Cheminf. 2016:
adding one XB term to Vina fixes halogenated-ligand docking at negligible
cost). PDB-REDO 2025 medians θ₁: Cl 143°, Br 154°, I 169° — the max(0,cos)²
gate with the 120° threshold captures the survey spread.

**Cost**: 1 extra dot + reuse of G(r) per halogen pair (halogens ≪ ligand
atoms; ~25% of lead molecules carry one). **Gain**: corrects currently
*LJ-blind* halogen interactions (σ/ε treat Cl/Br/I as fat carbons; Vina
misses XB entirely — same failure mode we'd inherit).

### (e) Unmet-buried-polar penalty (HYDE-lite) — RANK 5

```
E_pen = Σ_a w_pen · B_a · flag_pol(a) · max(0, 1 − s_a/2)
s_a   = # satisfied H-bonds of atom a (from term (a) counters, this frame)
w_pen = 1.5 kcal/mol  [R1 §1b: +1–3 per unmet buried polar; HYDE/LeScore]
```

HYDE (Reulecke et al. 2008; Schneider et al. 2013) and Glide-XP agree this
is the consensus false-positive killer (LeScore 2025: AUC 0.71 *with* the
penalty). Needs term (a) live (satisfaction counters), then it's a per-ligand-
atom post-pass over existing arrays — no new pair loop.

**Cost**: O(nLigAtoms) post-pass. **Gain**: kills the classic decoy class —
buried ligand N/O with no partner — that all radial terms over-reward.

### Term table summary (Pareto: gain/cost)

| # | Term | Equation (core) | Key params | Per-pair extra cost | Ranking-power evidence | Priority |
|---|---|---|---|---|---|---|
| a | Directional H-bond | −ε·G(r)·max(0,cosθ_D)²·max(0,cosα) | ε=2.0, r₀=3.0, σ=0.5 | 2 dot3 + 2 max (+1.3× on hbond subset, meas.) | 26–28% of Vina/Vinardo score mass; HYDE minimal-set | **1** |
| b | Salt bridge (live Coulomb) | 332·q_i·q_a/(ε(r)·r) | q=±1 beads; ε(r) existing | ≈0 (table edit) | dead channel; R1 §5 rank 4; screened −2.2@2.66 Å | **2** |
| c | Hydrophobic enclosure | ΔG_a·B_a·(1+enc) | enc=clip((n_H−4)/8), λ=1 | 1 add+cmp (PASS 1) | Glide-XP super-linear motif, prototype pocket ×1.38 | **3** |
| d | Halogen-lite | −ε_X·G_x·max(0,cosθ)² | ε Cl/Br/I 0.6/1.2/2.0, r_X=3.1 | 1 dot per halogen pair | VinaXB (Koebel 2016); 25% lead prevalence | **4** |
| e | Unmet-polar penalty | w_pen·B_a·(1−s_a/2) | w_pen=1.5 | O(nLig) post-pass | HYDE/Glide-XP consensus FP-killer; needs (a) | **5** |

Calibration discipline (R1 §3b Vinardo lesson): fewer, calibrated terms beat
more terms — weights get re-fit on PDBbind after each addition (R4/R5 loop),
not hand-tuned per-term.

---

## 3. Prototype validation numbers (offline, /tmp/opencode/)

Full output: `/tmp/opencode/r2_proto_output.txt`. Headlines (all through the
CG pipeline = Cα triplets → per-SS regression → virtual sites → gates):

**Angular discrimination (the feasibility question):**

| Contact | True r | Virtual r | Isotropic today | Directional CG | Verdict |
|---|---|---|---|---|---|
| HEPES N4⁺···O=Phe104 | 2.73 Å | 2.97 Å | **0.00** (class A, flag off; −0.59 if on) | **−1.87** (A_A=0.94) | native kept |
| HEPES O8···H–N–Leu32 | 2.96 Å | 3.14 Å | **0.00** (class H; −0.74 if on) | **−1.68** (A_D=0.87) | native kept |
| HEPES O1S⁻···NZ–Lys35 | 2.66 Å | 4.85 Å (virtual NZ) | −0.53 (flag on, no Coulomb) | −0.76 Coulomb + 0 H-bond (cone) | charged pair via term (b) |
| — decoy: C-side approach (same r) | 2.97 | — | −0.74 (blind) | **0.000** | killed |
| — decoy: θ_D=90° (same r) | 3.14 | — | −0.74 (blind) | **0.000** | killed |
| — decoy: H anti-oriented (same r) | 3.14 | — | −0.74 (blind) | **0.000** | killed |

Site placement honesty: virtual-O(Phe104) sits 0.30 Å from crystal O (helix
fit); the Lys35 type-mean NZ lands 2.42 Å off (rotamer noise — hence
charge-first design for sidechains). The NH-proxy validation: cos 0.95 vs the
true approach.

**Cost microbench (§F)**: directional 21.1 ns vs isotropic 15.7 ns per pair
= 1.3× on the flagged subset only; ≤5% projected frame impact.

**Enclosure (§E)**: benzene pocket enc=0.38 (7/7 H-class Cα < 7 Å) vs
solvent-exposed control enc=0.

---

## 4. Loop-2 implementation order (Pareto-ranked, per term table)

1. **Charge table** (term b): set q=±1 in `RES_CLASS` (or a parallel
   `_protQ` fill in `forcefield.js:184`) + His protonation hookup
   (reuse `chem/protonation.js` heuristic). ~10 lines. Biggest Δ-per-line;
   unblocks docking-screen validation of charged ligands immediately.
   **✅ DONE (Loop-2 S1, 2026-09-12)** — as the parallel `_protQ` fill:
   `CG_FORMAL_CHARGES` map (ff-params.js) + opt-in `par.binding.charges`
   (default OFF, preserves legacy bit-exactly; review-doc S1 spec).
   His hookup deferred (neutral 0 default, documented). No
   mean-neutralization (spec: ε(r) handles net charge).
   Measured salt-bridge ΔU: **−2.30 kcal/mol @ 2.66 Å** (Lys35(+1)···
   acetate O⁻), −1.42 @4, −0.84 @6, −0.03 @8 — matches the §2(b)/§D
   screened profile. Tests: tests/test_charges.js 20/20; test_all 32/32;
   default-off bit-identical (10-step golden hashes exact).
2. **Virtual-site arrays + per-SS regression constants** (§1d): precompute
   the 3×3 β tables offline (prototype values above as seeds; re-fit on a
   PDB survey in R4); rebuild frames in binding-pass preamble. ~60 lines.
3. **Directional H-bond gates** (term a) replacing the isotropic block
   `ff-binding.js:126-133`, incl. valence counters + coil widening. ~40 lines.
   *Simultaneously fixes the 52%-invisible-backbone bug (every residue gets
   sites — the flag becomes site-existence, not class).*
4. **Flag-coverage fix** (if 3 delayed): make `_protHB` per-residue always-1
   (backbone) as a stopgap — 1 line, immediately recovers Phe104-O-class
   contacts isotropically.
5. **Enclosure multiplier** (term c): class-masked density accumulation +
   PASS-2 multiplier. ~15 lines.
6. **Halogen gate** (term d): ligand C–X bond vectors (already in ligand
   graph) × protein O_acc virtual sites. ~25 lines.
7. **Unmet-polar penalty** (term e) once (a)'s counters exist. ~15 lines.

Sequencing rationale: 1+2 are prerequisites for everything (charges
independent; sites feed a/d/e); 3 is the flagship; 4 is a 1-line stopgap if
integration slips; 5–7 are independent increments validated by the R4/R5
rescoring loop. All inside `ff-binding.js` + `forcefield.js` parameter
plumbing; no new modules; all gated behind existing `bindOn`.

## 5. Explicit non-goals (stay heavy-only, per R1 CG-feasibility 0)

- **Cation-π and π-stacking** — need ring planes/normals; at Cα level no
  ring exists (R1 §3d/§4: feasibility 0). Heavy mode gets them (R3) from
  precomputed ring normals (ligand rings + Phe/Tyr/Trp/His reconstructible
  from heavy atoms). CG keeps aromatic-class LJ only.
- **Metal CG bead re-typing** — heavy already correct (`heavy.js:383-404`);
  CG metal sites are rare in the target set and springs would just re-create
  holo-restraints. Defer unless a metal-dependent test system enters.
- **Chalcogen/S–π standalone** — share the (d) σ-hole path at downscaled ε
  if ever needed; never standalone (R1 §4).
- **Explicit water / ordered-water displacement** — no water in the engine;
  enclosure (c) is the sanctioned proxy.
- **Absolute ΔG** — ranking only; entropy/GB error bars per Genheden & Ryde
  2015 make absolute claims from these terms meaningless (R1 §1c/§6).

## 6. References for this document's numbers

Levitt, J. Mol. Biol. 104:59 (1976) — first Cα→backbone reconstruction.
Maupetit et al. (PD2), Nucleic Acids Res. 34:e129 (2006) — 0.5–1.0 Å.
Gront/Kulp (BBQ), Nucleic Acids Res. 35:W327 (2007) — 0.6–0.8 Å.
Rohl et al. (PULCHRA), JCTC 5:1833 (2004). Dunbrack & Karplus, J. Mol. Biol.
231:543 (1993) — rotamer library multimodality. Baker & Hubbard, Protein Sci.
2:1810 (1993); McDonald & Thornton, J. Mol. Biol. 238:777 (1994) —
at-acceptor 90–180° hemisphere geometry. Trott & Olson (Vina) J. Comput.
Chem. 31:455 (2010); Quiroga & Villarreal (Vinardo) PLoS ONE 11:e0155180
(2016) — term-mass + minimal-set evidence. Friesner et al. (Glide XP), J.
Med. Chem. 49:6177 (2006) — enclosure. Reulecke et al. (HYDE) ChemMedChem
3:885 (2008); Schneider et al. JCIM 53:2063 (2013). Koebel et al. (VinaXB),
J. Cheminform. 8:27 (2016). Böhm, J. Comput.-Aided Mol. Des. 8:243 (1994) —
LUDI H-bond scale. Hendsch & Tidor, Protein Sci. 3:211 (1994) — salt-bridge
desolvation balance.

## 7. Validation

- Prototype: `node /tmp/opencode/r2_cgproto.mjs` (exit 0; output archived at
  `/tmp/opencode/r2_proto_output.txt`). No npm deps; reads repo PDBs only.
- **Loop-2 S1 landed**: `node tests/test_charges.js` (20/20 — charge census,
  salt-bridge ΔU vs closed form, screening profile, default-off bit-identity
  vs pre-S1 reference); `node tests/test_all.js` 32/32; golden / negative /
  bindlog 18/18 / bindviz 21/21 / `validate_binding_physics_r1.mjs` /
  `wikiskill_gate.js` OPEN all re-run PASS after the change.
- `src/` changes from S1: `ff-params.js` (CG_FORMAL_CHARGES + docs),
  `forcefield.js` (`chargesOn` flag + `_protQ` fill — default-off opt-in);
  `ff-binding.js` Coulomb kernel untouched (it was already correct, just
  starved of charges).
- `src/` untouched by this step (the 3 modified files in `git status` —
  main.js/ml-tier.js/network-panel.js — pre-date this work, from the prior
  commit's session; verified `git diff` shows no physics changes from R2).
- This file: `docs/BINDING_PHYSICS_R2.md`.
