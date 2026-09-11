# Exclusions — 1-2 / 1-3 / 1-4 non-bonded policy

> Audit: `src/heavy.js:524-542` (`pairKey i*1e6+j`) and `src/forcefield.js:250-309`
> ⇒ documented here for B13.

## Summary

| Class | Atoms | Coulomb / LJ | GB | Reason |
|-------|-------|-------------|----|--------|
| **1-2** (bonded) | i–j bonded per `buildTopology` or `bonds` segment list | **excluded** (`_excluded`) | excluded | Harmonic bond already holds distance; non-bonded would double-count and fight bond |
| **1-3** (angle) | i–k two bonds apart (`angles` list) | **excluded** | excluded | 1-3 distance restrained by angle term |
| **1-4** (dihedral / proper) | i–l three bonds apart (`propers` list) | **scaled 0.5** (`_scale14` = 0.5) | scaled 0.5 | AMBER/CHARMM convention — 1-4 non-bonded partially retained (electrostatics & LJ scaled) |
| **intra-ligand** | all pairs i<j with i,j ≥ `ligandStart` | **fully excluded** (`_excluded`) | excluded | Small-molecule geometry is governed by its own harmonic/angle/improper FF; grid LJ would otherwise crush rings / folded conformers (para carbons at 2.8 Å ≪ r_e 5 Å) |
| **coarse ENM / holo / native contacts** | ENM springs, holo springs, `nativeContacts` (`forcefield.js`) | excluded | excluded | Spring already encodes native contact |

## Heavy mode (`src/heavy.js:524-542`)

```js
this._excluded = new Set();
this._scale14 = new Map();
for (b in bonds)   _excluded.add(pairKey(bonds[2a], bonds[2a+1]));      // 1-2
for (a in angles)  _excluded.add(pairKey(angles[4a], angles[4a+2]));    // 1-3
for (p in propers) _scale14.set(pairKey(propers[5a], propers[5a+3]), 0.5); // 1-4 scaled 0.5
for (i=ligandStart; i<n; i++) for (j=i+1; j<n; j++) _excluded.add(pairKey(i,j)); // intra-ligand
```

`_nonBondedGrid` lookup:

```js
if (excluded.has(pairKey(i,j))) return;
const s14 = _scale14.get(pairKey(i,j)) ?? 1.0; // 0.5 for 1-4, else 1
```

## Coarse-grained (`src/forcefield.js:250-309`)

Same `_pairKey(i,j)=i<j ? i*1e6+j : j*1e6+i`, but also excludes ENM springs, holo springs, nativeContacts, and ligand internal FF pairs (bonds/angles/impropers) plus a transitive 1-4 walk for ligands (any `x–a–b–c` chain adds `pairKey(x,c)`). All 1-4 proper exclusions via that walk are excluded; protein 1-4 via ENM/bonds are excluded, not scaled — coarse model uses repulsive-only LJ (no 1-4 partial electrostatics).

## pairKey capacity and collision audit — `i*1e6 + j`

Choice `pairKey(i,j)= min*1e6 + max` encodes ordered pair in one integer key for `Set`/`Map`. With factor 1e6, supports `n < 1_000_000` atoms without collision (2 distinct pairs never map to same integer because `j < 1e6`). Real workloads: crambin 46 aa ≈ 327 heavy atoms, T4 lysozyme 164 aa ≈ 1308 heavy, 4HHB ~ 4300 heavy, largest PDBBind coreset < 20k heavy — three orders of magnitude below limit. If `n ≥ 1e6` collisions would occur (e.g. 0*1e6+1_000_001 collides with 1*1e6+1); such systems are out-of-scope (would need BigInt or string key or 1e9 factor). `forcefield.js:_pairKey` and `heavy.js:pairKey` both use 1e6; see `_pairKey` definition at `heavy.js:848` and `forcefield.js:383`.

Validation: `tests/test_exclusions.js` asserts 1-4 scale 0.5 on a 4-atom chain and that `intra-ligand` excluded count = `nLig*(nLig-1)/2`.
