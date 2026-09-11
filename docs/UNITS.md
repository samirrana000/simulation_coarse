# Units — simulation_coarse

> **Canonical base units** for all physics, integrator, and I/O. No hidden
> conversions — every module imports from `src/units.js`.

## Base Units

| Quantity | Symbol | Unit | Notes |
|----------|--------|------|-------|
| Length | L | **Ångström (Å)** | PDB coordinates are Å; all positions, cutoffs, RMSD in Å |
| Time | T | **picosecond (ps)** | Integrator `dt` in ps; trajectory times in ps; 1 ns = 1000 ps |
| Energy | E | **kcal/mol** | Force field energies, kT, bond springs, PMF, ΔG |
| Mass | M | **Dalton (Da)** = g/mol | Per-bead masses: Cα ≈110 Da, heavy atoms 12–65 Da |
| Temperature | Θ | **Kelvin (K)** | Baxter T in K; `KB * T` in kcal/mol |

## Derived & Conversion Constants

All defined in `src/units.js:1` (single source of truth).

| Constant | Symbol | Value | Unit | Role |
|----------|--------|-------|------|------|
| `KCONV` | KCONV | **418.4** | Da·Å²/ps² per kcal/mol | Force→acceleration: `a[Å/ps²] = KCONV * F[kcal/mol/Å] / m[Da]` |
| `KCAL_TO_DA_A2_PS2` | — | **418.4** | same as KCONV | Alias for backward compat |
| `KB` | `KB_KCAL` | **0.001987** (0.001987204) | kcal/mol/K | Boltzmann: `k_B T` in kcal/mol |
| `STANDARD_VOLUME` | V° | **1660.54** | Å³ | Standard-state volume per molecule at 1 M (for K_D) |
| `Coulomb const` | `K_ELEC` | 332.0 | kcal·Å/mol/e² | `E = K * q1 q2 / (ε r)` (heavy.js) |
| Screening | `SCREEN_LEN` | 8.0 | Å | Debye–Hückel screening (heavy.js) |

> **Why 418.4?** 1 kcal/mol = 4184 J/mol = 4184/(N_A)·J per particle = 418.4 Da·Å²/ps²
> (via `1 Da = 1.6605e-27 kg`, `1 Å = 1e-10 m`, `1 ps = 1e-12 s`).

## Time Scales

| Parameter | Typical | Unit | Description |
|-----------|---------|------|-------------|
| `dt` (Cα, no ligand) | 0.004 | ps | =4 fs, BAOAB stability limit for 110 Da + k≈160 kcal/mol/Å² |
| `dt` (Cα + ligand) | 0.0017 | ps | 1.7 fs, light 12 Da benzene requires smaller |
| `dt` (heavy) | 0.001 | ps | 1 fs, all-atom bonds 200–300 kcal/mol/Å² |
| `stridePs` | 2.0 | ps | Recording stride in ps of *simulation time* (not wall time) |
| `friction ζ` | 8.0 | ps⁻¹ | Langevin friction; c1=exp(-ζ dt) |
| `simSpeed` | 1.0 | ps/frame | `advance(stepsPerFrame)` where steps ≈ simSpeed/dt |

## Masses

| Bead type | Mass | Unit | Source |
|-----------|------|------|--------|
| Cα bead | 110 (slider 50–200) | Da | Mean residue mass (slider `mass`) |
| C (heavy) | 12.011 | Da | `heavy.js:heavyMass()` |
| N | 14.007 | Da | — |
| O | 15.999 | Da | — |
| S | 32.06 | Da | — |
| ZN | 65.38 | Da | — |
| Other | 14.0 fallback | Da | — |

Unified-atom ligand masses match heavy masses (e.g. benzene C=12.01).

## Energy Terms (all in kcal/mol)

- Bonds: ½ k_b (r−r0)², k_b≈160–300 kcal/mol/Å²
- Angles: ½ k_θ (θ−θ0)², k_θ=40 kcal/mol/rad² (heavy)
- Dihedrals: ½ k_φ (φ−φ0)², k∈2–20 kcal/mol/rad²
- LJ: 4ε[(σ/r)^12−(σ/r)^6], ε≈0.05–0.2 kcal/mol, σ≈3.2–3.8 Å
- Coulomb/GB: screened + GB ~ kcal/mol
- SASA: γ·ΔA, γ=0.0072 kcal/mol/Å²

## Temperature & Kinetic Energy

`T_inst = ( Σ ½ m v² * 0.5 / KCONV ) / (1.5 N k_B)` via equipartition.
Thermal velocity scale: `σ_v = sqrt(k_B T * KCONV / m)` Å/ps, validated as
`sqrt(KB_KCAL*T*KCONV/m)` (see `src/integrator.js:_rebuildThermal`).

## File I/O Units

- PDB: Å (ATOM x/y/z columns 31–54)
- MOL2: Å (ATOM x/y/z)
- XYZ: Å, time in comment line `time = X ps`
- PDB export (recorder): Å, `REMARK time = X ps`, `REMARK simulation_coarse vX ...`
- Analysis: all Å, ps, kcal/mol

## References

- `src/units.js:22` — KCONV=418.4
- `src/integrator.js:22` — derivation KCONV, thermal scale
- `src/heavy.js:28` — K_ELEC, screening lengths
- `src/physics/charges.js` — partial charges in e
- `src/version.js` — VERSION, BUILD_DATE for provenance

*Last updated: 2026-08-29. Keep this table in sync with `src/units.js` — any change to constants must update both.*
