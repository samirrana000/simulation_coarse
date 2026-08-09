#!/usr/bin/env python3
"""
export_esm_contacts.py
======================

Generate a residue-pair contact prior (`contacts.json`) for a coarse-grained
protein–ligand simulator. The output drives an Elastic Network Model where
per-pair spring constants are scaled by an ML contact prior.

USAGE
-----
    python3 ml/export_esm_contacts.py \\
        --pdb  data/4W52.pdb   |   --fetch 4W52 \\
        --out  data/4W52_contacts.json \\
        [--model esm2_t33_650M_UR50D] [--min-p 0.05] [--max-pairs 20000]

    --pdb FILE   read structure from a PDB file (exactly one of --pdb/--fetch)
    --fetch ID   download https://files.rcsb.org/download/{ID}.pdb (needs net)
    --out FILE   output JSON path (default: contacts.json)
    --model      ESM-2 model name, only used in ESM mode
    --min-p      keep only pairs with p_ij > min-p (default 0.05)
    --max-pairs  cap the number of kept pairs (default 20000)

TWO MODES
---------
1. ESM mode (preferred): requires `esm` (fair-esm) + PyTorch.
   Loads an ESM-2 model, tokenizes the single sequence, runs a forward pass
   (no grad), and derives residue-pair contact probability from the mean
   attention over the last 4 layers averaged over all heads (standard trick).
   The BOS/EOS special tokens are stripped, the matrix is symmetrized, and the
   diagonal is masked.

2. Heuristic fallback (default on machines without `esm`): a deterministic,
   physically-motivated pseudo prior so the pipeline runs end-to-end without
   heavy dependencies. This is a documented PLACEHOLDER that ESM replaces:
       p_ij = 0.6 * exp(-(i-j)^2 / (2 * 30^2))
                * (1 + 0.4 * hydroScore(res_i, res_j))
   where hydroScore is the Kyte–Doolittle scale product for hydrophobic pairs
   plus a bonus if both residues are aromatic. Sequence separation does NOT
   dominate; the 30-residue Gaussian only weakly favours local contacts.

NOTES
-----
* Sequence extraction is a coarse approximation for a Cα model: only the FIRST
  chain is used, one residue per (chain, resSeq) ATOM record, standard
  3-letter->1-letter mapping (unknown residues -> 'X'), ATOM records only
  (HETATM skipped), and parsing stops at the first TER record.
* To enable full ML mode:
      pip install fair-esm biopython
* If `esm` is missing the script prints
      [esm unavailable - using heuristic prior]
  and continues deterministically.
"""

import argparse
import datetime
import json
import math
import sys
import urllib.request
from collections import OrderedDict

# ---------------------------------------------------------------------------
# 3-letter -> 1-letter amino acid map
# ---------------------------------------------------------------------------

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # common variants
    "SEC": "C", "PYL": "K", "ASX": "B", "GLX": "Z", "XAA": "X",
    "UNK": "X",
}


def aa3_to_1(three_letter):
    """Map a 3-letter residue code to 1-letter; unknown -> 'X'."""
    return AA3_TO_1.get(three_letter.strip().upper(), "X")


# ---------------------------------------------------------------------------
# PDB parsing (coarse Cα model: first chain, one residue per (chain, resSeq))
# ---------------------------------------------------------------------------

def parse_pdb_text(text):
    """Parse a PDB blob. Returns (chain, sequence, n_res).

    Coarse approximation: only the FIRST chain is considered, one residue per
    (chain, resSeq) ATOM record, parsing stops at the first TER record,
    HETATM records are skipped.
    """
    residues = OrderedDict()  # (chain, resSeq, iCode) -> 1-letter code
    chain = None
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith("TER"):
            break
        if not line.startswith("ATOM"):
            continue
        # Fixed-width PDB columns.
        resname = line[17:20]
        chain_id = line[21]
        res_seq = line[22:26].strip()
        i_code = line[26] if len(line) > 26 else " "
        if chain is None:
            chain = chain_id
        if chain_id != chain:
            continue  # first chain only
        try:
            res_seq_int = int(res_seq)
        except ValueError:
            continue
        key = (chain_id, res_seq_int, i_code)
        if key not in residues:
            residues[key] = aa3_to_1(resname)
    if chain is None:
        raise ValueError("no ATOM records found in PDB text")
    sequence = "".join(residues.values())
    return chain, sequence, len(sequence)


def read_pdb_file(path):
    with open(path, "r", encoding="utf-8") as fh:
        return parse_pdb_text(fh.read())


def fetch_pdb(pdb_id):
    url = "https://files.rcsb.org/download/{}.pdb".format(pdb_id)
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8")


# ---------------------------------------------------------------------------
# Heuristic pseudo contact prior (deterministic, documented placeholder)
# ---------------------------------------------------------------------------

# Kyte-Doolittle hydropathy scale.
KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
    "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5,
    "K": -3.9, "R": -4.5,
}
AROMATIC = {"F", "W", "Y"}

SEQ_SEP_SIGMA = 30.0   # Gaussian width (residues) for local-contact smoothing
HYDRO_WEIGHT = 0.4     # weight of the hydrophobic/aromatic term


def hydro_score(a, b):
    """Kyte-Doolittle product for hydrophobic pairs, plus aromatic bonus."""
    kda = max(KD.get(a, 0.0), 0.0)
    kdb = max(KD.get(b, 0.0), 0.0)
    score = 0.0
    if kda > 0 and kdb > 0:
        score = min(1.0, (kda * kdb) / 12.0)
    if a in AROMATIC and b in AROMATIC:
        score = min(1.0, score + 0.3)
    return score


def heuristic_contacts(sequence):
    """Deterministic pseudo contact prior.

    p_ij = 0.6 * exp(-(i-j)^2 / (2 * 30^2)) * (1 + 0.4 * hydroScore(i, j))
    """
    n = len(sequence)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            d = j - i
            base = 0.6 * math.exp(-(d * d) / (2.0 * SEQ_SEP_SIGMA * SEQ_SEP_SIGMA))
            p = base * (1.0 + HYDRO_WEIGHT * hydro_score(sequence[i], sequence[j]))
            p = min(1.0, max(0.0, p))
            pairs.append((i, j, p))
    return pairs


# ---------------------------------------------------------------------------
# ESM mode (fair-esm, optional)
# ---------------------------------------------------------------------------

def esm_contacts(sequence, model_name):
    """Residue-pair contact probability from ESM-2 mean attention.

    Mean attention over the last 4 layers, averaged over all heads. Special
    BOS/EOS tokens are stripped, then the matrix is symmetrized and the
    diagonal masked.
    """
    import torch  # noqa: F401  (only importable when fair-esm present)
    import esm as esm_module

    model, alphabet = esm_module.pretrained.load_model_and_alphabet(model_name)
    model.eval()
    if hasattr(model, "transformer"):
        model.transformer.eval()
    if hasattr(model, "contact_head"):
        model.contact_head = None

    batch_converter = alphabet.get_batch_converter()
    batch_labels, batch_strs, batch_tokens = batch_converter([("protein", sequence)])
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    tokens = batch_tokens.to(device)

    with torch.no_grad():
        out = model(tokens, repr_layers=[], return_contacts=False)
    attentions = out["attentions"]  # list over layers

    # Collect per-layer attention tensors, dropping BOS/EOS token positions.
    layer_tensors = []
    for layer_attn in attentions[-4:]:
        attn_weights = layer_attn[0] if isinstance(layer_attn, tuple) else layer_attn
        # shape (batch, heads, seq_len, seq_len)
        layer_tensors.append(attn_weights[0].cpu().numpy())  # first batch item
    # (n_layers, heads, L, L) -> mean over layers and heads
    n_layers = len(layer_tensors)
    attn = sum(layer_tensors) / float(n_layers)
    attn = attn.mean(axis=0)  # (L, L)

    # Strip BOS/EOS: residues live at token positions 1..L (ESM-2 alphabets).
    L = len(sequence)
    attn = attn[1 : L + 1, 1 : L + 1]

    # Symmetrize + mask diagonal.
    attn = (attn + attn.T) / 2.0
    for i in range(L):
        attn[i, i] = 0.0

    pairs = []
    for i in range(L):
        for j in range(i + 1, L):
            p = min(1.0, max(0.0, float(attn[i, j])))
            pairs.append((i, j, p))
    return pairs


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def dump_json(out_path, source, chain, sequence, pairs, min_p, max_pairs):
    pairs = sorted(pairs, key=lambda t: t[2], reverse=True)
    if max_pairs > 0:
        pairs = pairs[:max_pairs]
    payload = {
        "source": source,
        "chain": chain,
        "sequence": sequence,
        "nRes": len(sequence),
        "meta": {
            "created": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "min_p": min_p,
        },
        "contacts": [[int(i), int(j), float(p)] for (i, j, p) in pairs],
    }
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")
    return payload


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="Export an ESM-derived (or heuristic) contact prior as JSON."
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb", metavar="FILE", help="read structure from a PDB file")
    src.add_argument("--fetch", metavar="PDBID", help="download PDB from RCSB")
    p.add_argument("--out", default="contacts.json", help="output JSON path")
    p.add_argument(
        "--model",
        default="esm2_t33_650M_UR50D",
        help="ESM-2 model name (ESM mode only)",
    )
    p.add_argument("--min-p", type=float, default=0.05, help="p_ij threshold")
    p.add_argument(
        "--max-pairs", type=int, default=20000, help="cap on kept pairs (0 = no cap)"
    )
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.pdb:
        try:
            chain, sequence, n_res = read_pdb_file(args.pdb)
        except (OSError, ValueError) as exc:
            print("ERROR: could not read PDB '{}': {}".format(args.pdb, exc))
            return 1
    else:
        try:
            text = fetch_pdb(args.fetch)
            chain, sequence, n_res = parse_pdb_text(text)
        except Exception as exc:  # network or parse failure
            print(
                "ERROR: could not fetch PDB '{}': {}".format(args.fetch, exc)
            )
            return 1

    try:
        import esm  # noqa: F401
        esm_available = True
    except ImportError:
        esm_available = False

    if esm_available:
        source = args.model
        pairs = esm_contacts(sequence, args.model)
    else:
        print("[esm unavailable - using heuristic prior]")
        source = "heuristic"
        pairs = heuristic_contacts(sequence)

    if not pairs:
        pairs = [(0, 1, 0.0)] if len(sequence) > 1 else []

    dump_json(args.out, source, chain, sequence, pairs, args.min_p, args.max_pairs)
    print(
        "wrote {} contacts (source={}, chain={}, nRes={}, {} pairs)".format(
            args.out, source, chain, len(sequence), len(pairs)
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
