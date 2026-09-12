# R1 — Binding Biophysics Literature Survey (grounded to this repo)

Step 1/7 of the binding-physics program. Scope: what to model, at what energy
scale, in which resolution. No `src/` physics was modified for this report.

## 0. Grounding: what the code actually does today (verified by reading)

CG binding kernel (`src/ff-binding.js`, wrapped by `ForceField._binding`
documented in `src/forcefield.js:723-772`):

- Cross 12-6 LJ, Lorentz–Berthelot combining, smoothstep-switched at
  `bindRcut = 9.0 Å` (`src/forcefield.js:115`, `src/ff-binding.js:104-111`).
- Screened Coulomb `U = 332.0637·q1·q2/(ε(r)·r)·sw` with
  `ε(r) = 4 + 76·tanh(r/8)` (`src/ff-binding.js:115-124`; doc
  `src/forcefield.js:728-730`). Protein bead charges are all `q: 0`
  (`src/ff-params.js:61-67`), so this term is **identically zero** in CG
  today (charges flow through `src/forcefield.js:184` unchanged).
- ISOTROPIC H-bond Gaussian, `EPSHB = 0.8 kcal/mol`, `HB_R0 = 3.2 Å`,
  `HB_W = 0.6 Å`, gated only on donor/acceptor **flags**
  (`src/ff-binding.js:57-58,127-133`). No angular term. Protein flag is
  residue-class based: P/Cp/Cn = 1 (`src/forcefield.js:185`); ligand flag is
  element based: N/O/F = true (`src/ff-params.js:76-87`).
- EEF1-lite burial desolvation: `g(r) = exp(−(r−4.5)²/2·1.8²)`,
  `B_a = 1 − exp(−n_a/3)`, `U = Σ ΔG_a·B_a`
  (`src/ff-binding.js:32,89-91,150-155`; doc `src/forcefield.js:748-764`).
  Per-atom transfer energies `ΔG_a ∈ [−0.55, −0.25]` kcal/mol
  (`src/ff-params.js:76-87`).
- Native holo springs `k = holoGamma = 0.5 kcal/mol/Å²`
  (`src/forcefield.js:116`), pairs with `r0 ≤ 6.0 Å`
  (`src/forcefield.js:276-303`), **excluded** from the binding pair pass
  (`src/ff-binding.js:83-84`) so no double counting.
- CG protein beads carry only 5 residue classes H/A/P/Cp/Cn with
  `σ ∈ [3.6, 4.1] Å, ε ∈ [0.10, 0.18]` kcal/mol
  (`src/ff-params.js:61-74`).

Heavy mode (`src/heavy.js`):

- Non-bonded grid kernel is LJ + GB/screened-Coulomb + directional-H-bond
  **only** (`src/heavy.js:947-1015`). No PME, no explicit water, no
  halogen / cation-π / π-stack / chalcogen terms anywhere (absence verified
  by reading the full kernel).
- GB-HCT default (`epsIn 4.0, epsOut 78.5, salt 0.15 M`,
  `src/heavy.js:501`; Still/Hawkins formalism `src/physics/gb.js:9-15,134-190`);
  OBC2/LCPO/membrane-slab are opt-in (`src/heavy.js:632-655,779-836`).
- H-bond module has a `cos²θ` angular gate with `r_eq = 2.9 Å, σ_r = 0.5 Å,
  epsHB = 2.5–3.0` (`src/physics/hbond.js:12-22,76-99`), but the angular
  gradient is **not propagated** ("stub", `src/physics/hbond.js:91`) and the
  production path `evaluatePair()` hard-codes ideal 180°
  (`src/physics/hbond.js:115-120`). Effective heavy H-bond today ≈ radial.
- Covalent topology: `r < 1.15·(r_cov(i)+r_cov(j))` and `r < 2.2 Å` cap
  (`src/heavy.js:278-299`; radii `src/ff-params.js:150-155`); metal
  coordination = distance springs to N/O/S donors within `coordR`
  (`src/heavy.js:383-404`; `coordR/coordN` table `src/ff-params.js:102-113`).
- Charges: united-atom AMBER-ish protein table + element fallbacks for
  ligand (`src/physics/charges.js:96-132`); GAFF2-lite typer opt-in
  (`src/heavy.js:215-249`); heuristic PROPKA-style protonation
  (`src/chem/protonation.js:1-33`, pKa shifts `HBOND_DIST 3.5 Å`,
  `SALTBRIDGE_DIST 4.5 Å` at lines 56-57).
- Non-polar: SASA pair-overlap model, `γ = 0.0072 kcal/mol/Å²`
  (`src/heavy.js:502`, `src/physics/sasa.js:19-28`).

Implication for R2/R3: the single biggest CG gaps vs literature are (1) no
H-bond directionality, (2) dead Coulomb channel (q=0), (3) no aromatic /
cation-π / halogen / metal-ligand distinction beyond σ/ε classes, (4) no
unmet-polar / ordered-water penalty. All four are cheap to add — see §4–5.

---

## 1. Binding free-energy decomposition and typical scales

```
ΔG_bind = ΔH − TΔS
        = (ΔE_vdW + ΔE_elec + ΔE_hb + ΔG_polar-solv + ΔG_nonpolar-solv)
          − T(ΔS_trans+rot + ΔS_conf,lig + ΔS_conf,prot + ΔS_solv)
```

Reference frame: mean drug-like ΔG° ≈ −36.5 kJ/mol ≈ **−8.7 kcal/mol**
(≈70% of 3025 affinities fall −11 to −6 kcal/mol; Forouzesh et al./Olsson ITC
compilation, reviewed in MDPI Biophysica 2024). Individual weak contacts are
each 1–5 kcal/mol, so net affinity is a **small difference of large
opposing terms** — the central fact behind enthalpy–entropy compensation
(Chodera & Mobley, Annu. Rev. Biophys. 2013, 42:121).

### 1a. Enthalpic drivers (gas-phase → net-in-water)

| Driver | Gas-phase / intrinsic | Net in water (what a scorer feels) | Source |
|---|---|---|---|
| vdW packing (per heavy-atom contact) | LJ well ε ≈ 0.1–0.2 kcal/mol/pair | Sums over buried surface; effective **−0.02 to −0.03 kcal/mol/Å²** buried non-polar surface | Chothia 1974; Eisenberg & McLachlan 1986; Sharp–Honig; our code γ=0.0072 is the low (MM-GBSA) end — see §1c |
| Neutral H-bond | 3–7 kcal/mol | **1–3 kcal/mol net** (desolvation cancels ~half+) | Böhm LUDI fitted −4.7 kJ/mol ≈ −1.1 kcal/mol per H-bond; Fersht; Schneider/HYDE |
| Charge-assisted H-bond / salt bridge | 80–120 kcal/mol bare ions | **0–2 kcal/mol net**, occasionally 3–4 buried | Sheinerman & Honig; Hendsch & Tidor (desolvation ≈ cancels Coulomb); frightening MM/PBSA −80 kcal/mol artefacts without entropy (Woo & Roux, PNAS 2005) |
| Cation-π | Li⁺–benzene −38, K⁺–benzene −19, NH₄⁺–benzene −19 gas | **2–5 kcal/mol in proteins** | Ma & Dougherty, Chem. Rev. 1997; Gallivan & Dougherty, PNAS 1999 (1 per 77 residues; E_es ≤ −2 kcal/mol cutoff); Dougherty, Acc. Chem. Res. 2013 |
| π–π stacking | benzene dimer −2 to −3 gas | **1–3 kcal/mol**, geometry-sensitive | Hunter & Sanders 1990 (electrostatic model); McGaughey et al. 1998 (vdW dominates stacked, electrostatics repulsive); Burley & Petsko 1985 |
| Halogen bond C–X···O | F ≈ 0, Cl 1–2, Br 2–3.5, I 3–5 (tunable to ~7 with strong EWG) | **1–4 kcal/mol net, low desolvation penalty** (key advantage over H-bonds) | Riley et al. 2017; Wilcken et al.; PDB-REDO survey 2025 (θ₁ medians F 132°/Cl 143°/Br 154°/I 169°); Verteramo et al., iScience 2024; VinaXB (Koebel et al., J. Cheminf. 2016) |
| Chalcogen (S/Se···O/N) | S ≈ 1–3, Se stronger | **1–3 kcal/mol**, σ-hole like halogen | Adhav et al. 2023 (S chalcogen vs H-bond in proteins); Politzer σ-hole reviews |
| Sulfur–π / S–aromatic | — | **~1–2 kcal/mol** (Met–aromatic enrichment) | Reid et al.; Vernon et al. 2018 (π-contacts underestimated by force fields) |
| Metal coordination (Zn²⁺/Mg²⁺/Fe) | 50–200 kcal/mol bare | **Net large but conditionally screened**; harmonic-spring + formal-charge treatment standard; chelate effect mostly entropic (+ΔS from released waters) | Irving–Williams; Ataie et al. 2008 (Zn geometry–affinity); chelate textbooks (ΔG° ≈ −54 kJ/mol for [Co(en)₃]³⁺ vs hexammine formation) |
| Hydrophobic transfer | — | **≈ 20–30 cal/mol/Å²** (0.8–1.2 kT per methyl ≈ 0.5–0.7 kcal/mol) | Tanford/Chothia ~25 cal/mol/Å²; Kyte 2003 (counts H–C bonds, not SASA); Sharp et al. |

### 1b. Entropic penalties (positive ΔG contributions at 300 K)

| Term | Scale | Source |
|---|---|---|
| Rigid-body (trans+rot) loss, 1 M std state | **+5 to +8 kcal/mol** (cratic ≈ +2.4; full RRHO 7–12 before solvent release offsets) | Janin; Murray & Verdonk; Gohlke & Case (3RT enthalpy correction, 1 M std); Woo & Roux decomposition (G_c+G_o ≈ +5–9) |
| Ligand configurational (per frozen rotor) | **+0.3 to +0.6 kcal/mol/rotor** (Böhm +1.4 kJ/mol ≈ +0.33; AutoDock +0.3; Vina Nrot nonlinear) | Böhm 1994 LUDI; Morris et al. AutoDock4; Trott & Olson Vina 2010; Veber (rotors+TPSA ↔ bioavailability) |
| Protein conformational (side-chain freezing, loop ordering) | **0 to +5 kcal/mol**, system-dependent; loop-state errors alone cause 1.5–2.5 kcal/mol enthalpy misses | Çınaroğlu & Biggin, Chem. Sci. 2023 (ZA-loop ⇒ RMSE 2.49→0.90); Deng et al. |
| Desolvation, polar (unmet buried H-bond) | **+1 to +3 kcal/mol per unmet donor/acceptor** — the Glide-XP/HYDE/LeScore penalty term | Friesner et al. Glide XP 2006; Reulecke et al. HYDE 2008; Schneider et al. 2013; LeScore 2025 (AUC 0.71 with H-bond penalty) |
| Desolvation, hydrophobic (solvent entropy gain) | **Favorable, −TΔS ≈ −20 to −50 cal/mol/Å²**; displacing "unhappy" waters adds −1 to −3 kcal/mol each | Dunitz; Young et al. WaterMap; Glide XP hydrophobic-enclosure term |
| Net experimental entropy −TΔS | Typically **−2 to +6 kcal/mol**; ITC shows ΔH ∈ [−15,+5], −TΔS ∈ [−5,+10] with near-cancellation | Chodera & Mobley 2013; Olsson ITC set (171 complexes); Du et al. 2016 review |

### 1c. Solvation-model numbers the code should respect

- Non-polar γ: literature spans **0.005 (MM-PBSA default) → 0.0072 (our
  `SasaModel`) → 0.02–0.03 (microscopic hydrophobic γ)** kcal/mol/Å².
  Genheden & Ryde review (Expert Opin. Drug Discov. 2015): γ choice barely
  moves congeneric rankings but absolute ΔG needs the larger value; fine for
  ranking, wrong for absolute ΔG. Our EEF1-lite per-atom ΔG_a (−0.25…−0.55)
  implicitly encodes ~0.01–0.02/Å² — consistent, keep.
- PB/GB without explicit-water entropy can err by **up to ~75 kJ/mol in
  relative ΔG** across buried-water cases (Genheden & Ryde §entropy). Never
  promise absolute ΔG from GB alone — R5 must state this.
- Normal-mode entropy dominates MM-PBSA cost and routinely **omits the true
  binding entropy**; SASA/rotor-count proxies perform "as good as normal
  modes" (Genheden & Ryde citing Kongsted/Ylilauri). License for our cheap
  rotor+burial entropy proxy in R4.

---

## 2. Noncovalent interaction catalog (numbers-first)

Geometry conventions: `r` = heavy–heavy distance unless noted; angles in
degrees. "Cheapest faithful form" = minimal term that reproduces ranking
signal at 60 fps in JS (radial part reuses the existing pair loop; angular
part ≤ 1 cos evaluation).

| # | Interaction | Strength (net in protein) | Key geometry (structural surveys) | Cheapest faithful functional form | Notes |
|---|---|---|---|---|---|
| 1 | H-bond (neutral N–H···O, O–H···O) | 1–3 kcal/mol each; Böhm fit 1.1; dense-buried up to ~4 | H···A 1.7–2.2, D···A **2.7–3.2**; D–H···A **>140°** (ideal 180); H···A–AA >120°; local-mode force const 0.2–0.4 mDyn/Å | Radial Gaussian at 2.9–3.2 Å × **cos²(D–H···A) × cos²(H···A–AA)** gate; cap donors/acceptors at valence (Chemgauss4/HYDE rule) | Single most cost-effective term. Current code has radial half only (§0) |
| 2 | Salt bridge (Asp/Glu···Arg/Lys/lig NH⁺) | 0–2 net (buried 3–4); strongest NCIs but most desolvated | N···O **2.8–3.5** (strong <2.3 H···O per JCIM 2018 survey); N–O–C angle ~90–120°; bridging waters common | Screened Coulomb (already present) + buried-polar penalty when unsatisfied; **no extra pair term needed** — fix is charges + GB, not new physics | Code gap is q=0 + missing unmet-polar penalty, not missing term |
| 3 | Cation-π (Arg/Lys/lig NH⁺···Phe/Tyr/Trp) | 2–5; E_es ≤ −2 cutoff (Gallivan–Dougherty); survives water better than salt bridge (only cation desolvated) | Cation-to-ring-centroid **3.5–5.5** (Na⁺–benzene optimum ~2.4 above plane, protein 4–6); Arg parallel-stacked (θ≈0°) preferred in PDB, T-shaped strongest in gas | Point-charge–quadrupole approx: `E = −k·q·Θ·(3cos²θ−1)/2r³` damped, or Gaussian-well × `cos²θ_plane` toward ring normal; needs ring centroid+normal (precomputable per aromatic residue/ligand ring) | Highest-value aromatic term; Vina/Vinardo miss it entirely |
| 4 | π–π stacking (Phe/Tyr/Trp···lig aryl) | 1–3; parallel-displaced ≈ T-shaped | Centroid–centroid **3.5–5.5**; inter-plane **3.3–3.8**; offset 1–2 Å (displaced) or 70–90° dihedral (T) | Anisotropic LJ: distance Gaussian × orientation factor from precomputed normals (`f = cos²φ_stack + cos²φ_T` blend, Hunter–Sanders-lite); falls back to isotropic LJ when normals unavailable | McGaughey: vdW dominates stacked ⇒ CG LJ already half-captures; angle term adds specificity |
| 5 | Halogen bond C–X···O/N/π (X=Cl/Br/I; F ≈ none) | Cl 1–2, Br 2–3.5, I 3–5 | X···O **3.0–3.5** (< ΣvdW); C–X···O **θ₁ 150–180°** (I 169° median, Br 154°, Cl 143°, F 132° — PDB-REDO 2025); X···O=C θ₂ 90–120°; X···π 3.3–3.8 to centroid | Extra-point (σ-hole) charge: +0.1–0.3 e point at 1.2–1.5 Å beyond X along C–X + standard Coulomb; or angle-gated Gaussian `E_XB·exp(−(r−r₀)²/2σ²)·cos²(θ₁−180°)` | VinaXB (Koebel 2016) proves empirical-XB-fixes-docking; ~25% of leads halogenated ⇒ high prevalence |
| 6 | Hydrophobic burial (C/S···C/S contacts) | 0.02–0.03/Å²; methyl ~0.5–0.7 | C···C 3.5–4.5; burial = loss of SASA | Keep EEF1-lite/SASA (already present); add unhappy-water bonus −1–2 kcal/mol per displaced buried water proxy (grid occupancy) if R4 budget allows | Already our best term; Glide-XP enclosure is the upgrade path |
| 7 | Metal coordination (Zn²⁺/Mg²⁺/Ca²⁺/Fe···N/O/S + ligand chelators) | Net 2–10+ per chelate; chelate effect entropic | M···N/O **1.9–2.5** (Zn–N ~2.0, Zn–O ~1.9–2.1, Mg–O ~2.1); coordination numbers 4 (tetrahedral Zn) / 6 (octahedral Mg/Ca/Fe) | Harmonic spring at ideal M–donor distance + formal-charge Coulomb (heavy already does this, `heavy.js:383-404`); CG: typed metal-bead LJ + charge | Heavy done; CG needs only a metal bead type, not new math |
| 8 | Chalcogen S/Se···O/N (Met/Cys–lig) | 1–3 (Se > S) | S···O **3.0–3.6**; C–S···O ~160–180° (σ-hole like halogen, weaker) | Same σ-hole machinery as #5 with smaller E and σ; share one code path | Cheap iff #5 built; never build standalone |
| 9 | Sulfur–π / S–aromatic (Met···aryl) | 1–2 | S···centroid **4–6**; S above ring plane | Isotropic LJ + small extra well toward aromatic centroid; subsumed by #4 machinery | Lowest standalone priority; free rider on #4 |
| 10 | (Reference) vdW sterics / clash | Repulsion >> attraction below ΣvdW | Clash < 0.8·ΣvdW | Keep 12-6 repulsion / Vina-style quadratic repulsion (already present both resolutions) | Baseline; Vinardo shows repulsion weight dominates (46% of score mass) |

Strength ladder (net, per-contact, kcal/mol, for R2 weighting):

```
salt-bridge-buried(3-4) > cation-pi(2-5) ≈ H-bond-buried(2-4) ≈ XB-iodo(3-5)
  > H-bond-typical(1-3) ≈ pi-stack(1-3) ≈ XB-bromo(2-3.5) ≈ chalcogen(1-3)
  > hydrophobic-methyl(0.5-0.7) ≈ S-pi(1-2) >> per-atom-vdW(0.1-0.2)
```

---

## 3. What CG / scoring-function literature says we can vs cannot capture

### 3a. Knowledge-based pair potentials (the CG precedent)

- **Miyazawa–Jernigan (1985, 1996)**: 20×20 quasi-chemical contact energies
  from PDB frequencies. Proof that **residue-level pair counts rank
  folds/binds** with zero angles. Direct license for our 5-class CG table —
  and diagnosis: our H/A/P/Cp/Cn ε spread (0.10–0.18) is flatter than MJ
  (several kT spread); R2 should re-derive MJ-like cross terms for
  ligand-element × residue-class.
- **Sippl PMF / Samudrala–Moult / DFIRE (Zhou & Zhou 2002)**: distance
  distributions beat single-cutoff contacts; DFIRE's finite-ideal-gas
  reference fixes the reference-state pathology. Lesson: if R2 adds one CG
  upgrade, make it a **distance-dependent** (multi-bin or Gaussian-mixture)
  pair term, not a deeper single well. DLIGAND2 (DFIRE for ligands) is still
  a top VS rescorer — distance-binned statistics work for ligands too.
- **Martini protein–ligand (Souza et al., Nat. Commun. 2020)**: CG with ~4
  heavy atoms/bead reproduces binding poses/pathways with calibrated LJ +
  charges. Our CG is coarser (1 bead/residue) so pose RMSD will be worse,
  but **ranking by burial + contacts transfers**.

### 3b. Vina / Vinardo term anatomy (the minimal-term-set experiment)

- **Vina (Trott & Olson 2010)**: 5 steric/contact terms + Nrot:
  `gauss1 (−0.0356) + gauss2 (−0.00516) + repulsion (+0.840) + hydrophobic
  (−0.0351) + H-bond (−0.587)`, Nrot nonlinear. Score-mass audit (Quiroga &
  Villarreal 2016): Gauss2-long-range 58%, H-bond 26%, hydrophobic 11%,
  Gauss1/repulsion ~5%. I.e. **three physics (sterics + hydrophobic + H-bond)
  carry essentially all ranking power**.
- **Vinardo (Quiroga & Villarreal, PLoS ONE 2016)**: drops Gauss2 entirely,
  re-optimizes radii, and **improves** scoring+docking+VS. Score mass shifts
  to repulsion 46% / H-bond 28% / hydrophobic 26%. Lesson for R2: fewer,
  better-calibrated terms beat more terms; do not add exotic terms before
  the big three are calibrated. Our CG already has exactly the Vinardo
  triple (LJ≈steric, EEF1-lite≈hydrophobic, Gaussian≈H-bond) — architecture
  is right, calibration + directionality are missing.
- **VinaXB (Koebel et al., J. Cheminf. 2016)**: adds one XB term to Vina,
  fixes halogenated-ligand docking with negligible cost. Template for how
  R3 should add #5.

### 3c. Glide XP / Chemgauss / HYDE motifs (what separates winners)

- **Glide XP (Friesner et al., J. Med. Chem. 2006)**: base GlideScore +
  (a) **hydrophobic-enclosure** rewards (lipophilic ligand atoms enclosed on
  multiple sides score super-linearly), (b) **special H-bond motifs**
  (neutral–neutral, charged, water-bridged each with own geometry/weight),
  (c) **desolvation penalties** for buried unsatisfied polars, (d) large
  hydrophobic-contact base. Takeaway: the delta over Vina is not new
  elements — it is **context (enclosure), H-bond typing, and penalties**.
  All three are implementable on our grids without new pair loops.
- **Chemgauss3/4 (OpenEye)**: Gaussian-smoothed shape + directional H-bond
  via lone-pair/polar-H positions + metal-chelator positions on precomputed
  grids. Validates the "precompute direction vectors, evaluate cheap
  Gaussians" engineering pattern R2/R3 should copy.
- **HYDE (Reulecke 2008; Schneider 2013)**: only H-bonds + hydrophobic +
  desolvation, atom-decomposed, no training set — yet competitive. Strongest
  evidence that the **minimal term set (H-bond + burial ± desolvation
  penalty) suffices** when balanced.

### 3d. CAN vs CANNOT (CG Cα+ligand vs heavy-atom)

| Capability | CG Cα + ligand heavy atoms | Heavy-atom (all heavy) |
|---|---|---|
| Shape complementarity / clash ranking | CAN (backbone of Vina/Vinardo/MJ signal) | CAN, better (true ΣvdW) |
| Hydrophobic burial ranking | CAN (EEF1-lite/SASA proxy; MJ precedent) | CAN (SASA/LCPO already) |
| Salt-bridge / Coulomb ranking | PARTIAL (needs residue charges on beads; geometry coarse) | CAN (GB + charges already) |
| H-bond counting | CAN (flags exist) | CAN |
| H-bond directionality / networks | CANNOT faithfully (no donor-H/acceptor-lone-pair vectors; backbone direction only infers ~50%) | CAN (needs H reconstruction or vector stubs — code stub exists, unwired) |
| Cation-π / π-stack specificity | CANNOT (no ring planes at Cα level; aromatic-class LJ only) | CAN with precomputed normals (ligand rings known; protein Phe/Tyr/Trp/His planes reconstructible from heavy atoms) |
| Halogen / chalcogen σ-hole | EFFECTIVELY CANNOT in CG (needs C–X bond vector — ligand side HAS it, protein side rarely matters; ligand-XB-to-backbone-carbonyl is the one CG-feasible sub-case) | CAN via extra-point or angle gate (ligand C–X vectors known) |
| Unmet-polar / ordered-water penalties | PARTIAL (burial × polarity proxy works; true water sites need heavy geometry) | CAN (burial + H-bond-satisfaction counting à la HYDE/XP) |
| Polarization, charge transfer, covalent inhibition | CANNOT (both resolutions; needs QM/MM) | CANNOT — declare out of scope |
| Absolute ΔG | CANNOT (no rigorous entropy; GB error bars ≫ signal) | CANNOT without FEP/TI — rank, don't promise |

**Minimal term set that gives ranking power** (consensus Vina+Vinardo+HYDE+
Glide-base): `sterics + hydrophobic-burial + directional-H-bond +
desolvation-penalty + rotor-entropy-penalty`. Everything else is a
single-digit-percent upgrade conditioned on the big five being right.

---

## 4. Resolution mapping with browser-feasibility scores (0–2)

Scale: 0 = infeasible/wrong-resolution · 1 = feasible with scaffolding
(precompute vectors, extra arrays, <10% frame cost) · 2 = trivially feasible
(fits existing pair loop, ~0 marginal cost at 60 fps).

| Interaction (§2 #) | CG-feasible? | Heavy-feasible? | Cheapest form in our engine | Priority driver |
|---|---|---|---|---|
| H-bond directionality (#1) | 1 (backbone N–H/O=C vectors inferable from Cα trace ± H-stub; side chains not) | 2 (vectors from heavy topology; stub `hbond.js:76-99` already written, needs wiring through `_nonBondedGrid`) | cos² gates on existing Gaussian | Highest relevance × feasibility |
| Unmet-polar / desolvation penalty (#2-adjacent) | 1 (burial × flag proxy) | 2 (HYDE-lite counting on grid) | Per-atom satisfaction counter, no new loop | Largest false-positive killer per Glide/HYDE |
| Hydrophobic enclosure upgrade (#6) | 1 (occupancy anisotropy from existing dens pairs) | 2 (SASA-neighbor count already looped) | Multi-side enclosure multiplier on ΔG_a·B_a | Super-linear ranking win, ~free |
| Salt bridge via real charges (#2) | 1 (assign ±1 to Cp/Cn beads; ε(r) channel already coded but dormant) | 2 (already live through GB) | Flip q table + GB-desolvation does the rest | Unblocks a dead code path |
| Rotor entropy penalty (config) | 2 (ligand graph already parsed; count rotors once per molecule) | 2 (same) | +0.3–0.5 kcal/mol × Nrot added to score (Vina-style, not dynamics) | Cheapest ΔG correction in the survey |
| Cation-π (#3) | 0 (no ring planes at Cα; aromatic-class LJ is the ceiling) | 1 (precompute 4 residue ring normals + ligand ring normals per frame or per placement) | Quadrupole-lite well × cos² to normal | Top aromatic term; heavy-only |
| π-stack (#4) | 0 (same reason) | 1 (same normals as cation-π, shared code path) | Anisotropic LJ blend | Free rider on cation-π normals |
| Halogen σ-hole (#5) | 1 (ligand C–X vector known; acceptors = backbone carbonyl proxy; narrow sub-case only) | 1 (extra-point charge or cos² gate; piggybacks pair loop) | Angle-gated Gaussian, Cl/Br/I typed (F excluded) | 25% of leads halogenated; VinaXB precedent |
| Metal–ligand (#7) | 1 (typed metal bead + charge; springs exist in heavy only) | 2 (already shipped: springs + formal charges) | CG: new bead type, no new math | Heavy done; CG small |
| Chalcogen (#8) / S–π (#9) | 0 standalone | 1 as free rider on #5/#4 paths | Shared σ-hole / stack code with downscaled E | Never standalone |

Cost honesty: angular terms need per-frame normals (aromatic planes,
backbone vectors). Protein ring normals change slowly ⇒ cache per
placement/minimization step, not per MD substep; ligand normals update with
rigid-body transform only. Estimated <5% frame budget. Extra-point charges
add zero new loops (fold into Coulomb pass).

---

## 5. Deliverable: interaction priority table (top 8 by relevance × feasibility)

Ranked by (fraction of PDBbind complexes where it decides ranking) ×
(browser cost)⁻¹, heaviest weight to terms Vina/Vinardo/Glide-XP/HYDE agree on.

| Rank | Interaction | Strength (net) | Geometry gate | Cheapest functional form | CG? | Heavy? | Why this rank |
|---|---|---|---|---|---|---|---|
| 1 | Directional H-bond | 1–3 kcal/mol | D···A 2.7–3.2 Å, D–H···A >140° | Gaussian(r; 2.9–3.2) × cos²θ_D × cos²θ_A, valence-capped | 1 | 2 | Carries ~26–28% of Vina/Vinardo score mass; current isotropic −0.8 underweights buried H-bonds ~2–3× |
| 2 | Unmet-polar / desolvation penalty | +1–3 per buried unsatisfied polar | burial > 0.5 + no H-bond partner | HYDE-lite per-atom penalty, reuses burial+H-bond loops | 1 | 2 | Glide-XP/HYDE/LeScore consensus false-positive killer; zero new loops |
| 3 | Hydrophobic enclosure (upgrade of existing burial) | 0.02–0.03/Å², super-linear when enclosed | lipophilic atom contacted on ≥2–3 sides | enclosure multiplier × (ΔG_a·B_a) from existing pair anisotropy | 1 | 2 | Glide-XP's signature win; free on current grids |
| 4 | Salt bridge via live charges | 0–2 net (buried 3–4) | N···O 2.8–3.5 Å | activate Cp/Cn ±1 + existing ε(r)/GB (no new term) | 1 | 2 | Unblocks dead Coulomb channel (§0); biggest Δ-per-line-changed |
| 5 | Rotor configurational-entropy penalty | +0.3–0.5 / rotor | Nrot from ligand graph (static per molecule) | score += w·Nrot (Vina-style, outside dynamics) | 2 | 2 | Cheapest absolute-ΔG correction; Böhm/AutoDock/Vina agree on scale |
| 6 | Cation-π | 2–5 | cation···centroid 3.5–5.5 Å, approach along ring normal | quadrupole-lite well × cos²θ_normal, cached normals | 0 | 1 | Strongest aromatic term; invisible to Vina — differentiation opportunity |
| 7 | Halogen σ-hole (Cl/Br/I) | 1–4 (I>Br>Cl; F=0) | X···O 3.0–3.5 Å, C–X···O 150–180° | angle-gated Gaussian or +0.1–0.3 e extra point | 1-sub | 1 | VinaXB precedent; 25% lead prevalence; tiny code |
| 8 | π-stack (+ S–π / chalcogen as free riders) | 1–3 (S–π 1–2, chalcogen 1–3) | centroid 3.5–5.5 Å, plane 3.3–3.8 Å | anisotropic-LJ blend on shared normals from #6 | 0 | 1 | Shares #6 normals ⇒ marginal cost ~0 once #6 exists |

Explicitly deferred: metal CG re-typing (heavy already correct; CG bead type
is R2 filler), polarization/charge-transfer/covalency (QM territory),
explicit/ordered waters beyond a displacement proxy (engine has no water),
PME (no periodicity; screened Coulomb + GB is the documented substitute).

---

## 6. Key energy scales (pocket card for R2/R3)

- Typical drug ΔG° ≈ **−6 to −11 kcal/mol**; design battles are won/lost in
  **±1–2 kcal/mol** ⇒ every term below matters at that resolution.
- Per-contact net ladder: buried salt bridge 3–4 ≈ cation-π 2–5 ≈ buried
  H-bond 2–4 ≈ I-XB 3–5 > typical H-bond 1–3 ≈ stack 1–3 ≈ Br-XB 2–3.5 ≈
  chalcogen 1–3 > methyl burial 0.5–0.7 >> single vdW pair 0.1–0.2.
- Penalties: rigid-body +5–8; rotor +0.3–0.5 each; unmet buried polar +1–3;
  unhappy-water displacement −1–3 (favorable).
- Compensation warning (Chodera–Mobley): optimizing ΔH alone is
  systematically cancelled by −TΔS ⇒ R4/R5 must score **ΔG proxies**, never
  raw interaction sums.

## 7. Three most important papers for R2/R3

1. **Trott & Olson, J. Comput. Chem. 2010, 31:455 — "AutoDock Vina:
   improving the speed and accuracy of docking."** The minimal-term-set
   blueprint (gauss/hydrophobic/H-bond/Nrot) and the weight table R2 should
   start from. Companion: **Quiroga & Villarreal, PLoS ONE 2016 —
   "Vinardo"** (drop Gauss2, re-fit radii; proves fewer-calibrated-beats-more).
2. **Friesner et al., J. Med. Chem. 2006, 49:6177 — "Extra Precision Glide:
   docking and scoring incorporating a model of hydrophobic enclosure."**
   The enclosure + special-H-bond + desolvation-penalty motifs that separate
   winners from Vina-class scorers; directly maps to priority ranks 2–3.
3. **Genheden & Ryde, Expert Opin. Drug Discov. 2015, 10:449 — "The
   MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities."**
   The error-bar bible: γ-range, ~75 kJ/mol buried-water risk, normal-mode
   entropy ≈ SASA/rotor proxies, absolute-vs-relative ΔG discipline R5 must
   enforce. (Theory anchor behind it: **Chodera & Mobley, Annu. Rev.
   Biophys. 2013** on compensation; **Gallivan & Dougherty, PNAS 1999** for
   cation-π energetics/geometry; **Koebel et al., J. Cheminf. 2016 (VinaXB)**
   for the halogen implementation template.)

Supporting: Böhm, J. Comput.-Aided Mol. Des. 1994/1998 (LUDI scales:
H-bond −1.1, ionic −2.0, lipo −0.04/Å², rotor +0.33 kcal/mol, ΔG₀ +1.3);
Zhou & Zhou 2002 (DFIRE reference state); Souza et al., Nat. Commun. 2020
(Martini protein–ligand CG precedent); McGaughey et al. 1998 + Hunter &
Sanders 1990 (stack physics); Riley et al./PDB-REDO 2025 (XB geometry);
Schneider et al. 2013 (HYDE desolvation balance); Woo & Roux, PNAS 2005
(absolute-ΔG decomposition discipline).

---

## 8. Validation

- This file: `docs/BINDING_PHYSICS_R1.md` (written, this step).
- Code-citation check script: `scripts/validate_binding_physics_r1.mjs`
  (zero-dep node): asserts every `src:` line-number claim in §0
  (EPSHB/HB_R0/HB_W, ε(r), R0/SIG/NS, ΔG_a range, q=0 protein charges,
  bindRcut/holoGamma, holo-exclusion, heavy kernel term inventory, GB/SASA
  constants, BOND_SLACK/cap, salt-bridge cutoffs). Run:
  `node scripts/validate_binding_physics_r1.mjs`.
- Physics untouched: `git status --porcelain src/` is clean (docs + script
  only), no npm deps added.
