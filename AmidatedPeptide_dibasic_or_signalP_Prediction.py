#!/usr/bin/env python3
"""
Prediction of candidate C-terminally amidated peptides from a secretome.

Accompanies Wiggenhorn et al., "Global profiling of peptide amidation to
discover bioactive fragments of the secretome."

Method
------
Each classically secreted protein is scanned from its N-terminus to its
C-terminus for a glycine-dibasic tripeptide (GKR, GKK, GRK or GRR). This motif
indicates proprotein convertase cleavage at the dibasic site followed by
PAM-mediated amidation of the resulting C-terminal glycine.

On finding such a motif, an N-terminal cleavage site is sought upstream:
  1. the most 3' dibasic pair (KR, KK, RK, RR) lying downstream of the signal
     peptide, otherwise
  2. the annotated signal peptide cleavage site.

The predicted amidated peptide is the sequence between that N-terminal site and
the residue preceding the glycine.

Usage
-----
    python predict_amidated_peptides.py \\
        --secretome data/mouse_secretome.tsv \\
        --output    results/mouse_amidated_peptides.csv

    # restrict to the subset examined experimentally (< 15 residues)
    python predict_amidated_peptides.py \\
        --secretome data/mouse_secretome.tsv \\
        --output    results/mouse_predicted_under15.csv \\
        --max-length 15

Input
-----
A UniProt TSV export with, in column order:
    0 Entry (accession)   2 Gene names        4 Modified residue
    7 Signal peptide      8 Peptide           9 Sequence
    10 Propeptide         11 Chain
Rows lacking a SIGNAL annotation are skipped.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterator, NamedTuple

# ── Motif definitions ────────────────────────────────────────────────────────
# Glycine-dibasic tripeptides: cleavage at the dibasic site leaves a C-terminal
# glycine, which PAM converts to the amide.
AMIDATION_MOTIFS = ("GRR", "GKR", "GKK", "GRK")

# Dibasic pairs marking proprotein convertase cleavage at the N-terminus.
DIBASIC_MOTIFS = ("RR", "KR", "KK", "RK")

# A motif at the very start of a sequence cannot represent a processed site.
MIN_MOTIF_POSITION = 2

# Residues of flanking sequence reported either side of the peptide.
FLANK = 5

OUTPUT_COLUMNS = [
    "peptide name", "peptide sequence", "uniprot accession", "GN",
    "peptide length", "start position", "stop position",
    "around N-term", "around C-term", "signal start", "signal peptide end",
    "peptide info", "propeptide info", "chain info", "mod info",
]


class SecretedProtein(NamedTuple):
    accession: str
    name: str
    sequence: str
    signal_end: str          # kept as text: UniProt may give '?'
    peptide_info: str
    propeptide_info: str
    chain_info: str
    mod_info: str


# ── Input ────────────────────────────────────────────────────────────────────
def read_secretome(path: Path) -> Iterator[SecretedProtein]:
    """Yield proteins carrying a SIGNAL annotation from a UniProt TSV export."""
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= 1 or fields[0] == "Entry":
                continue

            ptm = fields[7].replace('"', "")
            if not ptm.startswith("SIGNAL"):
                continue

            # 'SIGNAL 1..22; /evidence=...' -> '22'
            signal_end = fields[7].split(";")[0].split("..")[1]

            yield SecretedProtein(
                accession       = fields[0],
                name            = fields[2] or fields[0],
                sequence        = fields[9],
                signal_end      = str(signal_end),
                peptide_info    = fields[8],
                propeptide_info = fields[10] if len(fields) > 10 else "",
                chain_info      = fields[11] if len(fields) > 11 else "",
                mod_info        = fields[4],
            )


# ── Prediction ───────────────────────────────────────────────────────────────
def find_n_terminal_site(fragment: str, signal_end: int,
                         max_upstream: int | None) -> int:
    """
    Return the 0-based start index of the predicted peptide.

    Preference is for the most 3' dibasic pair downstream of the signal peptide;
    failing that, the signal peptide cleavage site; failing that, position 0.

    If ``max_upstream`` is set, a dibasic pair is only accepted when it lies
    within that many residues of the amidation motif.
    """
    search_from = 0
    if max_upstream is not None:
        search_from = max(0, len(fragment) - max_upstream - len(DIBASIC_MOTIFS[0]))

    positions = [fragment.rfind(motif, search_from) for motif in DIBASIC_MOTIFS]
    dibasic_pos = max(positions)

    if dibasic_pos > signal_end:
        return dibasic_pos + 2          # cleave immediately after the pair
    if signal_end != -1:
        return signal_end
    return 0


def predict_peptides(protein: SecretedProtein, min_length: int,
                     max_length: int | None, max_upstream: int | None) -> list[dict]:
    """Predict all amidated peptides for one protein, N-terminus to C-terminus."""
    seq = protein.sequence
    signal_end = int(protein.signal_end) if protein.signal_end.isdigit() else -1

    rows: list[dict] = []
    site_count = 0

    for position in range(len(seq) - 1):
        # NOTE: the parentheses around the motif comparison are required.
        # Without them Python's precedence (and binds tighter than or) would
        # apply the position test to the final motif only, admitting motifs at
        # positions 0-1 for the other three. Reported by a reviewer; corrected
        # here. Re-running produced no change to the predicted peptide set, as
        # all previously identified candidates lay beyond position 1.
        is_motif = seq[position:position + 3] in AMIDATION_MOTIFS
        if not (is_motif and position > MIN_MOTIF_POSITION - 1):
            continue

        fragment = seq[:position]
        start = find_n_terminal_site(fragment, signal_end, max_upstream)
        peptide = fragment[start:]

        if len(peptide) < min_length:
            continue
        if max_length is not None and len(peptide) >= max_length:
            continue

        # Flanking context. max(0, ...) guards the sequence start: a bare
        # negative index would silently wrap to the C-terminus.
        upstream   = seq[max(0, start - FLANK - 1):start + FLANK]
        downstream = seq[max(0, position - 4):position + FLANK + 2]

        site_count += 1
        gene = protein.name.split(" ")[0]

        # Peptides are numbered in the order encountered, i.e. by position in
        # the precursor: the second site in a protein becomes PEP-GENE-2.
        # NOTE: this counter is per UniProt entry. Where one gene appears under
        # multiple accessions, peptides from each are numbered independently
        # and must be reconciled downstream.
        name = f"PEP-{gene}" if site_count == 1 else f"PEP-{gene}-{site_count}"

        rows.append({
            "peptide name":        name,
            "peptide sequence":    peptide,
            "uniprot accession":   protein.accession,
            "GN":                  protein.name,
            "peptide length":      len(peptide),
            "start position":      start + 1,          # report 1-based
            "stop position":       start + len(peptide),
            "around N-term":       upstream,
            "around C-term":       downstream,
            "signal start":        "yes" if start == signal_end else "no",
            "signal peptide end":  protein.signal_end,
            "peptide info":        protein.peptide_info,
            "propeptide info":     protein.propeptide_info,
            "chain info":          protein.chain_info,
            "mod info":            protein.mod_info,
        })

    return rows


# ── Entry point ──────────────────────────────────────────────────────────────
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Predict C-terminally amidated peptides from a secretome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--secretome", type=Path, required=True,
                        help="UniProt TSV export of secreted proteins")
    parser.add_argument("--output", type=Path, required=True,
                        help="destination CSV")
    parser.add_argument("--min-length", type=int, default=3,
                        help="minimum peptide length, inclusive")
    parser.add_argument("--max-length", type=int, default=None,
                        help="exclusive upper bound on peptide length. Omit for "
                             "no cutoff; use 15 for the subset examined "
                             "experimentally")
    parser.add_argument("--max-upstream", type=int, default=None,
                        help="only accept a dibasic pair within this many "
                             "residues of the amidation motif. Omit to search "
                             "the whole upstream sequence")
    args = parser.parse_args(argv)

    if not args.secretome.exists():
        parser.error(f"secretome file not found: {args.secretome}")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    seen: set[str] = set()
    n_proteins = n_written = 0

    with args.output.open("w", newline="") as handle:
        # csv.writer quotes fields containing commas. UniProt annotation
        # columns frequently do, so manual string concatenation would shift
        # columns for those rows.
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()

        for protein in read_secretome(args.secretome):
            n_proteins += 1
            for row in predict_peptides(protein, args.min_length,
                                        args.max_length, args.max_upstream):
                # Identical sequences predicted from more than one protein are
                # reported once, on first encounter.
                if row["peptide sequence"] in seen:
                    continue
                seen.add(row["peptide sequence"])
                writer.writerow(row)
                n_written += 1

    print(f"proteins scanned : {n_proteins}", file=sys.stderr)
    print(f"peptides written : {n_written}", file=sys.stderr)
    print(f"output           : {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
