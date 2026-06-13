"""
Extract peptides from *_deepimmuno.fasta, deduplicate, and write {sample}_tcr_input.csv
in the same directory.

Input FASTA naming: {sample}_deepimmuno.fasta
Output CSV naming: {sample}_tcr_input.csv (same directory as FASTA)
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import OrderedDict

from Bio import SeqIO

FASTA_SUFFIX = "_deepimmuno.fasta"
CSV_SUFFIX = "_tcr_input.csv"


def extract_sample_name(fasta_filename: str) -> str:
    """Extract sample name from {sample}_deepimmuno.fasta."""
    stem = os.path.splitext(fasta_filename)[0]
    lower = stem.lower()
    suffix = FASTA_SUFFIX[:-len(".fasta")]  # _deepimmuno
    if lower.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def find_deepimmuno_fastas(root_dir: str) -> list[str]:
    """Recursively find all *_deepimmuno.fasta under a directory."""
    matches: list[str] = []
    for dirpath, _, filenames in os.walk(root_dir):
        for name in filenames:
            if name.lower().endswith(FASTA_SUFFIX):
                matches.append(os.path.join(dirpath, name))
    return sorted(matches)


def dedupe_peptides_from_fasta(fasta_path: str) -> list[str]:
    """Read FASTA and deduplicate by peptide sequence (preserve first-seen order)."""
    seen: OrderedDict[str, None] = OrderedDict()
    for record in SeqIO.parse(fasta_path, "fasta"):
        peptide = str(record.seq).strip().upper()
        if not peptide:
            continue
        if peptide not in seen:
            seen[peptide] = None
    return list(seen.keys())


def write_tcr_input_csv(peptides: list[str], csv_path: str) -> None:
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["peptide"])
        writer.writeheader()
        for peptide in peptides:
            writer.writerow({"peptide": peptide})


def process_fasta(fasta_path: str) -> tuple[str, int]:
    sample_name = extract_sample_name(os.path.basename(fasta_path))
    out_dir = os.path.dirname(os.path.abspath(fasta_path))
    csv_path = os.path.join(out_dir, f"{sample_name}{CSV_SUFFIX}")

    peptides = dedupe_peptides_from_fasta(fasta_path)
    write_tcr_input_csv(peptides, csv_path)
    return csv_path, len(peptides)


def parse_arguments():
    p = argparse.ArgumentParser(
        description="Deduplicate peptides from *_deepimmuno.fasta and write *_tcr_input.csv in the same directory"
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="Root folder containing {sample}_deepimmuno.fasta files (searched recursively)",
    )
    return p.parse_args()


def main():
    args = parse_arguments()
    root_dir = os.path.abspath(args.input)
    if not os.path.isdir(root_dir):
        raise SystemExit(f"Directory not found: {root_dir}")

    fasta_files = find_deepimmuno_fastas(root_dir)
    if not fasta_files:
        raise SystemExit(
            f"No *{FASTA_SUFFIX} files found; check directory: {root_dir}"
        )

    print(f"Found {len(fasta_files)} FASTA file(s) under {root_dir}")
    total_peptides = 0

    for fasta_path in fasta_files:
        csv_path, n = process_fasta(fasta_path)
        total_peptides += n
        print(
            f"✅ {os.path.basename(fasta_path)} → "
            f"{n} deduplicated peptides → {csv_path}"
        )

    print(f"\nAll done; wrote {total_peptides} peptides total (across samples).")


if __name__ == "__main__":
    main()
