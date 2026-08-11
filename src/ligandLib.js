/**
 * ligandLib.js — Built-in small-molecule library for pocket placement (item 2).
 *
 * Each entry is a complete Tripos MOL2 file (inline text) plus a display
 * name and SMILES for identification. The same strings feed parseMol2()
 * (mol2.js), so everything that applies to user-supplied MOL2 input applies
 * here too:
 *   - hydrogens are dropped by the parser (united-atom model),
 *   - elements must be covered by LIG_ELEMENT (ff-params.js) — this library
 *     deliberately stays inside C / N / O / S / F / Cl,
 *   - ring aromaticity is re-derived from the geometry at build time
 *     (ligand.js), so ring bond lengths below sit in the 1.30–1.48 Å
 *     resonance window where applicable.
 *
 * Atom ids are contiguous 1..n within each file; bond entries use `ar`
 * for aromatic bonds, but only connectivity is used downstream.
 *
 * Shape: { id, name, smiles, mol2 }  — `mol2` parsed lazily by the caller.
 */

const BENZENE = `@<TRIPOS>MOLECULE
BNZ
 6 6 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 C1 1.3960 0.0000 0.0000 C.ar 1 BNZ 0.0000
 2 C2 0.6980 1.2090 0.0000 C.ar 1 BNZ 0.0000
 3 C3 -0.6980 1.2090 0.0000 C.ar 1 BNZ 0.0000
 4 C4 -1.3960 0.0000 0.0000 C.ar 1 BNZ 0.0000
 5 C5 -0.6980 -1.2090 0.0000 C.ar 1 BNZ 0.0000
 6 C6 0.6980 -1.2090 0.0000 C.ar 1 BNZ 0.0000
@<TRIPOS>BOND
 1 1 2 ar
 2 2 3 ar
 3 3 4 ar
 4 4 5 ar
 5 5 6 ar
 6 6 1 ar
`;

const PHENOL = `@<TRIPOS>MOLECULE
PHN
 7 7 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 C1 1.3950 0.0000 0.0000 C.ar 1 PHN 0.0000
 2 C2 0.6980 1.2080 0.0000 C.ar 1 PHN 0.0000
 3 C3 -0.6980 1.2080 0.0000 C.ar 1 PHN -0.1000
 4 C4 -1.3950 0.0000 0.0000 C.ar 1 PHN 0.0000
 5 C5 -0.6980 -1.2080 0.0000 C.ar 1 PHN 0.0000
 6 C6 0.6980 -1.2080 0.0000 C.ar 1 PHN 0.0000
 7 O1 -2.0900 1.7080 0.0000 O.3 1 PHN -0.4000
@<TRIPOS>BOND
 1 1 2 ar
 2 2 3 ar
 3 3 4 ar
 4 4 5 ar
 5 5 6 ar
 6 6 1 ar
 7 3 7 1
`;

const TOLUENE = `@<TRIPOS>MOLECULE
TOL
 7 7 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 C1 1.3950 0.0000 0.0000 C.ar 1 TOL 0.0000
 2 C2 0.6980 1.2080 0.0000 C.ar 1 TOL 0.0000
 3 C3 -0.6980 1.2080 0.0000 C.ar 1 TOL 0.0000
 4 C4 -1.3950 0.0000 0.0000 C.ar 1 TOL 0.0000
 5 C5 -0.6980 -1.2080 0.0000 C.ar 1 TOL 0.0000
 6 C6 0.6980 -1.2080 0.0000 C.ar 1 TOL 0.0000
 7 C7 2.9000 0.0000 0.0000 C.3 1 TOL 0.0000
@<TRIPOS>BOND
 1 1 2 ar
 2 2 3 ar
 3 3 4 ar
 4 4 5 ar
 5 5 6 ar
 6 6 1 ar
 7 1 7 1
`;

const CHLOROBENZENE = `@<TRIPOS>MOLECULE
CLB
 7 7 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 C1 1.3950 0.0000 0.0000 C.ar 1 CLB 0.0000
 2 C2 0.6980 1.2080 0.0000 C.ar 1 CLB 0.0000
 3 C3 -0.6980 1.2080 0.0000 C.ar 1 CLB 0.0000
 4 C4 -1.3950 0.0000 0.0000 C.ar 1 CLB 0.0000
 5 C5 -0.6980 -1.2080 0.0000 C.ar 1 CLB 0.0000
 6 C6 0.6980 -1.2080 0.0000 C.ar 1 CLB 0.0000
 7 CL1 -3.1350 0.0000 0.0000 Cl 1 CLB 0.0000
@<TRIPOS>BOND
 1 1 2 ar
 2 2 3 ar
 3 3 4 ar
 4 4 5 ar
 5 5 6 ar
 6 6 1 ar
 7 4 7 1
`;

const INDOLE = `@<TRIPOS>MOLECULE
IND
 9 10 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 C1 1.3960 0.0000 0.0000 C.ar 1 IND 0.0000
 2 C2 0.6980 1.2090 0.0000 C.ar 1 IND 0.0000
 3 C3 -0.6980 1.2090 0.0000 C.ar 1 IND 0.0000
 4 C4 -1.3960 0.0000 0.0000 C.ar 1 IND 0.0000
 5 C5 -0.6980 -1.2090 0.0000 C.ar 1 IND 0.0000
 6 C6 0.6980 -1.2090 0.0000 C.ar 1 IND 0.0000
 7 C7 -1.1286 -2.5348 0.0000 C.ar 1 IND 0.0000
 8 N1 0.0000 -3.3548 0.0000 N.ar 1 IND -0.3000
 9 C8 1.1286 -2.5348 0.0000 C.ar 1 IND 0.0000
@<TRIPOS>BOND
 1 1 2 ar
 2 2 3 ar
 3 3 4 ar
 4 4 5 ar
 5 5 6 ar
 6 6 1 ar
 7 5 7 ar
 8 7 8 ar
 9 8 9 ar
 10 9 6 ar
`;

const IMIDAZOLE = `@<TRIPOS>MOLECULE
IMI
 5 5 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 N1 1.1867 0.0000 0.0000 N.ar 1 IMI -0.3500
 2 C1 0.3667 1.1286 0.0000 C.ar 1 IMI 0.1500
 3 N2 -0.9592 0.6975 0.0000 N.ar 1 IMI -0.3500
 4 C2 -0.9592 -0.6975 0.0000 C.ar 1 IMI 0.0000
 5 C3 0.3667 -1.1286 0.0000 C.ar 1 IMI 0.1500
@<TRIPOS>BOND
 1 1 2 ar
 2 2 3 ar
 3 3 4 ar
 4 4 5 ar
 5 5 1 ar
`;

const ACETATE = `@<TRIPOS>MOLECULE
ACT
 4 3 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 C1 0.0000 0.0000 0.0000 C.2 1 ACT 0.4000
 2 O1 1.2300 0.0000 0.0000 O.co2 1 ACT -0.7000
 3 O2 -0.6200 1.0700 0.0000 O.co2 1 ACT -0.7000
 4 C2 -0.8500 -1.2400 0.0000 C.3 1 ACT 0.0000
@<TRIPOS>BOND
 1 1 2 2
 2 1 3 1
 3 1 4 1
`;

const ETHANOLAMINE = `@<TRIPOS>MOLECULE
ETA
 4 3 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 C1 0.0000 0.0000 0.0000 C.3 1 ETA 0.1000
 2 C2 1.5200 0.1000 0.0000 C.3 1 ETA 0.1000
 3 O1 2.2400 -1.0200 0.4000 O.3 1 ETA -0.4000
 4 N1 -0.5200 -0.8000 1.0500 N.3 1 ETA -0.3000
@<TRIPOS>BOND
 1 1 2 1
 2 2 3 1
 3 1 4 1
`;

const DMSO = `@<TRIPOS>MOLECULE
DMS
 4 3 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 S1 0.0000 0.0000 0.0000 S.o 1 DMS 0.3000
 2 O1 0.0000 1.4800 0.0000 O.2 1 DMS -0.5000
 3 C1 1.0500 -0.6200 0.9000 C.3 1 DMS 0.1000
 4 C2 -1.0500 -0.6200 -0.9000 C.3 1 DMS 0.1000
@<TRIPOS>BOND
 1 1 2 2
 2 1 3 1
 3 1 4 1
`;

const CAFFEINE = `@<TRIPOS>MOLECULE
CAF
 14 15 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
 1 N1 1.3950 0.0000 0.0000 N.ar 1 CAF -0.2000
 2 C1 0.6975 1.2081 0.0000 C.ar 1 CAF 0.3000
 3 N2 -0.6975 1.2081 0.0000 N.ar 1 CAF -0.2000
 4 C2 -1.3950 0.0000 0.0000 C.ar 1 CAF 0.3000
 5 C3 -0.6975 -1.2081 0.0000 C.ar 1 CAF 0.3000
 6 C4 0.6975 -1.2081 0.0000 C.ar 1 CAF 0.3000
 7 N3 1.1286 -2.5348 0.0000 N.ar 1 CAF -0.2000
 8 C5 0.0000 -3.3548 0.0000 C.ar 1 CAF 0.3000
 9 N4 -1.1286 -2.5348 0.0000 N.ar 1 CAF -0.2000
 10 O1 0.6975 2.4381 0.0000 O.2 1 CAF -0.4000
 11 O2 -2.6250 0.0000 0.0000 O.2 1 CAF -0.4000
 12 C6 2.8950 0.0000 0.0000 C.3 1 CAF 0.1000
 13 C7 -0.6975 2.7081 0.0000 C.3 1 CAF 0.1000
 14 C8 -2.4086 -3.3048 0.0000 C.3 1 CAF 0.1000
@<TRIPOS>BOND
 1 1 2 ar
 2 2 3 ar
 3 3 4 ar
 4 4 5 ar
 5 5 6 ar
 6 6 1 ar
 7 6 7 ar
 8 7 8 ar
 9 8 9 ar
 10 9 5 ar
 11 2 10 2
 12 4 11 2
 13 1 12 1
 14 3 13 1
 15 9 14 1
`;

/** Built-in ligand library (item 2). `mol2` strings feed parseMol2() directly. */
export const LIGAND_LIBRARY = [
  { id: "benzene",        name: "Benzene",        smiles: "c1ccccc1",                mol2: BENZENE },
  { id: "phenol",         name: "Phenol",         smiles: "Oc1ccccc1",               mol2: PHENOL },
  { id: "toluene",        name: "Toluene",        smiles: "Cc1ccccc1",               mol2: TOLUENE },
  { id: "chlorobenzene",  name: "Chlorobenzene",  smiles: "Clc1ccccc1",              mol2: CHLOROBENZENE },
  { id: "indole",         name: "Indole",         smiles: "c1ccc2[nH]ccc2c1",        mol2: INDOLE },
  { id: "imidazole",      name: "Imidazole",      smiles: "c1c[nH]cn1",              mol2: IMIDAZOLE },
  { id: "acetate",        name: "Acetate",        smiles: "CC(=O)[O-]",              mol2: ACETATE },
  { id: "ethanolamine",   name: "Ethanolamine",   smiles: "NCCO",                    mol2: ETHANOLAMINE },
  { id: "dmso",           name: "DMSO",           smiles: "CS(C)=O",                 mol2: DMSO },
  { id: "caffeine",       name: "Caffeine",       smiles: "Cn1c(=O)c2[nH]cnc2n(C)c1=O", mol2: CAFFEINE },
];
