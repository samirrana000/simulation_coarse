# Phase 2 — Chemistry, Topology & Protonation (opt-in, default OFF)

- `src/chem/protonation.js` — `assignProtonationStates(residues|atoms,{pH})`, PROPKA-style
  `pKa_eff=model+desolvation+H-bond+Coulomb`; HIS→HIE/HID/HIP, ASP/GLU→ASH/GLH when buried w/o salt bridge,
  LYS/ARG→LYN/ARN when buried-hydrophobic, CYS→CYX (SG–SG≤2.3 Å) / CYM (metal≤2.8 Å); every call carries `reasons[]`.
  UI: `protAssign` checkbox in Structure → Advanced chemistry.
- `src/chem/gaff2_mapper.js` — `typeMolecule/assignCharges (Gasteiger+AM1-BCC-lite, integer-renormalized)/getBondParams/getAngleParams/getTorsionParams/buildMoleculeParams`.
  Types ca/c3/c2/na/nb/n/o/oh/os/s/f/cl/br/i; metals fall back to element; never throws. UI: `gaffLig` checkbox.
- `src/chem/stereo.js` — `checkStereochemistry` (tetrahedral volumes, CIP-lite R/S, inversion flags), `validatePlanarity` (Newell plane).
- `src/chem/metals.js` — `detectCoordination/classifyGeometry (linear→square-antiprismatic)/assignIdealAngles/enforceCoordination`
  (radial k=40 + cross-angle k=20).

Validation (4W52+benzene/indole/caffeine): HIS A31→HIP (salt bridge ASP70 3.86 Å); 1CRN 6 CYS→CYX correct pairs;
benzene 6×ca net 0, indole 8×ca+na, caffeine full table coverage; 2000-step relax max bond dev 0.0001 Å, planarity rms 0.009.
`test_all` 32/32, `test_mol2_fidelity` 15/15 PASS.
