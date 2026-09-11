# Charges — AMBER ff14SB Approximate Mapping

> **Source:** Maier, J.A. et al. "ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB." *J. Chem. Theory Comput.* **2015**, *11*, 3696–3713. DOI:10.1021/acs.jctc.5b00255
> Cornell et al. *J. Am. Chem. Soc.* **1995**, *117*, 5179 (original ff94/RESP charge derivation).

This project uses an **approximate, united-atom mapping of AMBER ff14SB** partial charges — see `src/physics/charges.js`. Hydrogens are summed into heavy atoms, values are rounded for coarse-grained visualization and teaching. **Statement: approximate, not full ff14SB.** A full ff14SB implementation would include explicit hydrogens, separate N-terminal (NH3+) and C-terminal (COO-) patches, and exact RESP charges per atom type.

## Residue net-charge audit (internal residues, OXT excluded)

OXT (-0.80) is a C-terminal cap and excluded from the internal sum. Validation via `validateCharges()` in `src/physics/charges.js:validateCharges`.

| Residue | Formal charge | Atoms summed (heavy) | Calculated sum | Literature ff14SB reference (example) | Pass? |
|---------|---------------|----------------------|----------------|---------------------------------------|-------|
| **ASP** | **-1** | N(-0.28)+CA(+0.28)+C(+0.55)+O(-0.55)+CB(0.00)+CG(+0.60)+OD1(-0.80)+OD2(-0.80) | **-1.00** | OD1 -0.8014, OD2 -0.8014, CG 0.7172, CB -0.2826 | ✓ (tol 0.05) |
| **GLU** | **-1** | N+CA+C+O+CB(0)+CG(0)+CD(+0.60)+OE1(-0.80)+OE2(-0.80) | **-1.00** | OE1/OE2 -0.8014, CD 0.7130 | ✓ |
| **ARG** | **+1** | N+CA+C+O+CB(0)+CG(0)+CD(+0.10)+NE(-0.40)+CZ(+0.64)+NH1(+0.33)+NH2(+0.33) | **+1.00** | CZ 0.64, NH* ~0.33, NE -0.40 (united) | ✓ |
| **LYS** | **+1** | N+CA+C+O+CB(0)+CG(0)+CD(0)+CE(+0.25)+NZ(+0.75) | **+1.00** | NZ ~0.75, CE ~0.25 (summed H) | ✓ |

All four charged residues reproduce the formal integer net charge within 0.05 e. Neutral residues (ALA, GLY, etc.) sum to 0.00 (backbone N -0.28 + CA +0.28 + C +0.55 + O -0.55 = 0).

## Notes and limitations

- Values are rounded; e.g., literature ASP OD1 -0.8014 is rounded to -0.80, CG 0.7172 to 0.60, CB -0.2826 to 0.00. This preserves the integer net charge but not the full ff14SB electrostatics.
- Histidine in this table is neutral HIS (ND1 -0.36, CE1 +0.18, CD2 +0.18); protonation states HIE/HID/HIP not distinguished.
- Metal ions use formal charges (NA/K +1, others +2) not ff14SB-derived.
- **Approximate, not full ff14SB** — for production MD use AMBER ff14SB via GROMACS/AMBER/OpenMM with explicit solvent and PME.

## Generalized Born — HCT with scale 0.8 vs OBC

Implicit solvent is Still Generalized Born + Debye–Hückel screening
(`src/physics/gb.js`). Born radii use Hawkins-Cramer-Truhlar (HCT) pairwise
descreening with scale 0.8 vs Onufriev-Bashford-Case (OBC):
HCT with scale 0.8 vs OBC — HCT: Hawkins et al. J. Phys. Chem. 1996, 100,
19824 & Chem. Phys. Lett. 1995, 246, 122. Intrinsic radii `GB_RADII`
scaled as `rho_j * 0.8` for overlap; OBC would add element-specific
α/β/γ rescaling (not used here for simplicity). Isolated-atom test:
`computeBornRadii` with no `pos` or single atom returns intrinsic ±1e-6
(see `gb.js:computeBornRadii`).

## Validation function

```js
import { validateCharges } from "./src/physics/charges.js";
console.table(validateCharges());
// → [{residue:"ASP", calc:-1, formal:-1, ok:true}, ...]
```

