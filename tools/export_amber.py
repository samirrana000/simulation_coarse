#!/usr/bin/env python3
"""
export_amber.py — OpenMM XML interchange stub for 4W52 heavy (I85)

Generates an OpenMM XML skeleton for the 4W52 all-atom heavy system
(Cα→heavy per src/heavy.js, 1308 protein heavy atoms + ligand) as a
placeholder for AMBER/GAFF → OpenMM FF interchange.

Usage:
    python tools/export_amber.py --pdb 4w52.pdb --out system.xml
    python tools/export_amber.py --fetch 4W52 --out system.xml  # not yet: needs internet

Status: STUB / skeleton — placeholders only, not a full serialized System.
Real OpenMM inter-op would:
  1) Assign AMBER ff14SB + GAFF2 atom types/charges (antechamber/GAFF or openmmforcefields),
  2) Build topology with Modeller, add H, run 200 steps minimization,
  3) Serialize System to XML via openmm.XmlSerializer and compare
     browser HeavyForceField vs OpenMM GBSA energies on the same frames.

Energy fidelity target (aspirational): energy diff <2% aspirational between
browser Heavy GB/SA (src/heavy.js HeavyForceField + src/physics/gb.js GeneralizedBorn)
and OpenMM GBSA (GBn2/OBC2) on 20 4W52 NVT frames — see docs/OPENMM_REF.md R>0.75.
Absolute energies will differ by a systematic offset (radii/cutoffs), so the
<2% is rank-/scale-relative after mean-subtraction or per-component (LJ/electrostatics).
This gate is aspirational until real force-field assignment lands.

Output is a minimal XML skeleton with correct element masses/charges/positions
for 4W52 heavy atoms so that `openmm` / `openmmforcefields` can be wired later
without changing the CLI contract.
"""

import argparse
import sys
import xml.etree.ElementTree as ET
from xml.dom import minidom
import datetime

def _pretty_xml(elem):
    rough = ET.tostring(elem, encoding="unicode")
    return minidom.parseString(rough).toprettyxml(indent="  ")

AA3_TO_ELEMENT_GUESS = {
    "ALA":"C","ARG":"N","ASN":"N","ASP":"O","CYS":"S","GLN":"N","GLU":"O",
    "GLY":"C","HIS":"N","ILE":"C","LEU":"C","LYS":"N","MET":"S","PHE":"C",
    "PRO":"C","SER":"O","THR":"O","TRP":"N","TYR":"O","VAL":"C",
}

def parse_heavy_atoms(pdb_text):
    """Count heavy atoms (non-H, non-water) — close to src/heavy.js parseHeavy for 4W52."""
    atoms = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM","HETATM")):
            continue
        # PDB columns 77-78 element, else atom name
        elem = line[76:78].strip().upper() if len(line)>=78 else ""
        name = line[12:16].strip()
        if not elem:
            first = name.lstrip("0123456789")[0:1].upper() if name else "C"
            elem = "C" if first=="C" else ("N" if first=="N" else ("O" if first=="O" else ("S" if first=="S" else "C")))
        if elem == "H":
            continue
        res = line[17:20].strip().upper()
        if res in ("HOH","WAT","SOL"):
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except: continue
        atoms.append((elem, x, y, z))
    return atoms

def build_skeleton(atoms, pdb_id="4W52"):
    """
    Build OpenMM XML skeleton:
      <ForceField><System id=4W52 heavy><Particles>...
    Placeholder bonds/angles omit true FF params — marked TODO(amber).
    """
    ff = ET.Element("ForceField")
    # Provenance
    info = ET.SubElement(ff, "Info")
    gen = ET.SubElement(info, "Source")
    gen.text = f"simulation_coarse tools/export_amber.py stub — {pdb_id} heavy — {datetime.datetime.now(datetime.timezone.utc).isoformat()}"
    note = ET.SubElement(info, "Note")
    note.text = "STUB: placeholder masses/charges only; wire AMBER ff14SB/GAFF2 + GBSA for real interchange. Energy diff <2% aspirational (see docs/OPENMM_REF.md)."
    # System scaffold
    system = ET.SubElement(ff, "System", id=pdb_id, heavy=str(len(atoms)))
    particles = ET.SubElement(system, "Particles")
    for i, (elem, x, y, z) in enumerate(atoms, start=1):
        # placeholder LJ/charge — use mean-neutralized / ff-params defaults
        mass = {"C":12.01,"N":14.01,"O":16.00,"S":32.06,"P":30.97,"F":19.0,"CL":35.45,"BR":79.9}.get(elem, 12.01)
        ET.SubElement(particles, "Particle", id=str(i), element=elem, mass=f"{mass:.4f}", charge="0.0000",
                      x=f"{x:.3f}", y=f"{y:.3f}", z=f"{z:.3f}")
    # Placeholder force blocks — not valid FF, just skeletons so downstream parsers find sections
    forces = ET.SubElement(system, "Forces")
    ET.SubElement(forces, "HarmonicBondForce", placeholder="true", note="TODO(amber): replace with AMBER covalent bonds; src/heavy.js BOND_SLACK 1.15 cap 2.2A")
    ET.SubElement(forces, "HarmonicAngleForce", placeholder="true")
    ET.SubElement(forces, "PeriodicTorsionForce", placeholder="true")
    ET.SubElement(forces, "NonbondedForce", placeholder="true", coulombConst="332.06371", note="GB/SASA placeholder — use OpenMM GBSA OBC2; energy diff <2% aspirational")
    ET.SubElement(forces, "GBSAOBCForce", placeholder="true", solventDielectric="78.5", soluteDielectric="4.0")
    return ff

def main(argv=None):
    ap = argparse.ArgumentParser(description="Export 4W52 heavy OpenMM XML skeleton (stub). Energy diff <2% aspirational target.")
    src = ap.add_mutually_exclusive_group(required=False)
    src.add_argument("--pdb", help="input PDB file (default: tries 4w52.pdb)")
    src.add_argument("--fetch", help="PDB ID to fetch from RCSB (requires network; stub will warn)")
    ap.add_argument("--out", default="system.xml", help="output OpenMM XML path (default: system.xml)")
    ap.add_argument("--id", default="4W52", help="PDB ID for XML metadata")
    args = ap.parse_args(argv)

    pdb_text = None
    if args.pdb:
        try:
            with open(args.pdb, "r", encoding="utf-8") as fh:
                pdb_text = fh.read()
        except OSError as e:
            print(f"ERROR: cannot read --pdb {args.pdb}: {e}", file=sys.stderr)
            return 1
    elif args.fetch:
        print(f"WARNING: --fetch {args.fetch} stub — fetching requires network; using local 4w52.pdb if present", file=sys.stderr)
        try:
            with open("4w52.pdb", "r", encoding="utf-8") as fh:
                pdb_text = fh.read()
        except OSError:
            # minimal synthetic 4W52 heavy placeholder so XML still writes
            pdb_text = "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00           C\nATOM      2  CA  GLY A   2      13.800  10.000  10.000  1.00  0.00           C\n"
    else:
        try:
            with open("4w52.pdb", "r", encoding="utf-8") as fh:
                pdb_text = fh.read()
        except OSError:
            try:
                with open("data/4w52.pdb", "r", encoding="utf-8") as fh:
                    pdb_text = fh.read()
            except OSError:
                print("No PDB found — writing 2-atom synthetic placeholder skeleton (stub).", file=sys.stderr)
                pdb_text = "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00           C\nATOM      2  CA  GLY A   2      13.800  10.000  10.000  1.00  0.00           C\n"

    atoms = parse_heavy_atoms(pdb_text)
    if not atoms:
        print("WARNING: no heavy atoms parsed — skeleton will be empty (stub)", file=sys.stderr)
    ff = build_skeleton(atoms, pdb_id=args.id)
    xml_str = _pretty_xml(ff)
    # Pretty adds extra header; ensure output
    with open(args.out, "w", encoding="utf-8") as out:
        out.write(xml_str)
    print(f"Wrote OpenMM XML skeleton for {args.id} heavy ({len(atoms)} atoms) → {args.out}")
    print("Energy diff <2% aspirational — real AMBER/GBSA assignment needed to hit target (see docs/OPENMM_REF.md).")
    return 0

if __name__ == "__main__":
    sys.exit(main())
