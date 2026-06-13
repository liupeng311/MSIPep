"""
Extract Class I peptides from MSFragger TSV output and export FASTA.

Default thresholds (aligned with common peptide search practice and Comet post-processing):
- MHC-I length 8–11 aa (same as comet-pep);
- expectation / expectscore (common MSFragger column names) and E-value-like columns ≤ threshold (default 0.01) when present;
- hyperscore ≥ 15 when column exists (common conservative lower bound for closed search);
- Remove typical decoy protein matches; include target accessions only in headers.
"""

import argparse
import glob
import os
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

_MS_DIR = Path(__file__).resolve().parent
if str(_MS_DIR) not in sys.path:
    sys.path.insert(0, str(_MS_DIR))

from peptide_filters import (
    in_class1_length_range,
    row_has_target_hit,
    split_protein_identifiers,
    target_proteins_only,
)


def _get_float(row, *keys, default=None):
    """
    Try multiple column names in order; return the first value convertible to float.
    Compatible with pandas Series (NaN) and plain dict (None / empty string).
    Same logic as the same-named function in comet-pep.py for consistent E-value lookup.
    """
    for key in keys:
        v = row.get(key)
        if v is None:
            continue
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        import math
        if math.isnan(fv):
            continue
        return fv
    return default


def _norm_col(c: str) -> str:
    return str(c).strip().lower().replace(" ", "_")


def _rename_df_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [_norm_col(c) for c in df.columns]
    return df


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract peptides from MSFragger TSV, apply filters, and write FASTA (decoy removal, length and score thresholds)."
    )
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing .tsv files")
    parser.add_argument("-o", "--output_fasta", required=True, help="Output FASTA")
    parser.add_argument(
        "--min_length",
        type=int,
        default=8,
        help="Class I minimum peptide length (default 8, same as Comet post-processing)",
    )
    parser.add_argument("--max_length", type=int, default=11, help="Class I maximum peptide length (default 11)")
    parser.add_argument(
        "--max_expectation",
        type=float,
        default=0.01,
        help="Maximum expectation (E-value scale, default 0.01). ≤0 disables.",
    )
    parser.add_argument(
        "--min_hyperscore",
        type=float,
        default=15.0,
        help="Minimum hyperscore (default 15). ≤0 disables.",
    )
    parser.add_argument(
        "--min_nextscore_ratio",
        type=float,
        default=0.0,
        help="If >0 and nextscore exists: require hyperscore/nextscore ≥ this ratio (reduce ambiguous hits).",
    )
    parser.add_argument(
        "--keep_decoy",
        action="store_true",
        help="Keep decoy matches (default: remove DECOY/REV etc.)",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    exclude_decoy = not args.keep_decoy
    input_dir = args.input_dir
    output_fasta = args.output_fasta

    use_exp = args.max_expectation > 0
    use_hyp = args.min_hyperscore > 0
    use_ratio = args.min_nextscore_ratio > 0

    tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
    if not tsv_files:
        print(f"No .tsv files found in folder {input_dir}; exiting.")
        return

    print(f"Found {len(tsv_files)} TSV file(s) in folder {input_dir}.")
    print(f"Results will be written to {output_fasta}")
    if exclude_decoy:
        print(
            f"Filters: decoy removal; peptide length {args.min_length}-{args.max_length} aa; "
            f"expectation≤{args.max_expectation if use_exp else '(disabled)'}; "
            f"hyperscore≥{args.min_hyperscore if use_hyp else '(disabled)'}"
        )
    else:
        print("Warning: keeping decoy matches.")

    peptide_protein_map = defaultdict(set)
    warn_exp = False
    warn_hyp = False
    warn_next = False

    for tsv_file in tsv_files:
        print(f"Processing file: {tsv_file}")
        try:
            df = pd.read_csv(tsv_file, sep="\t", dtype=str, low_memory=False)
        except Exception as e:
            print(f"Error reading file {tsv_file}: {e}")
            continue

        df = _rename_df_columns(df)

        if "protein" not in df.columns or "peptide" not in df.columns:
            print(f"File {tsv_file} missing protein / peptide columns; skipping.")
            continue

        has_exp = any(
            c in df.columns
            for c in ("expectation", "expectscore", "e-value", "expect", "e_value")
        )
        has_hyp = "hyperscore" in df.columns
        has_next = "nextscore" in df.columns

        if use_exp and not has_exp:
            if not warn_exp:
                print(
                    "Note: expectation / expectscore / e-value column not found in TSV; skipping E-value-like filtering."
                )
                warn_exp = True
        if use_hyp and not has_hyp:
            if not warn_hyp:
                print("Note: hyperscore column not found in TSV; skipping hyperscore filtering.")
                warn_hyp = True
        if use_ratio and not has_next:
            if not warn_next:
                print("Note: nextscore column not found; skipping hyperscore/nextscore ratio filtering.")
                warn_next = True

        for _, row in df.iterrows():
            protein = row.get("protein")
            peptide = row.get("peptide")
            if pd.isna(protein) or pd.isna(peptide):
                continue
            peptide = str(peptide).strip()
            if not peptide:
                continue

            if exclude_decoy and not row_has_target_hit(protein):
                continue

            if not in_class1_length_range(
                peptide, args.min_length, args.max_length
            ):
                continue

            if use_exp and has_exp:
                evf = _get_float(
                    row,
                    "expectation", "expectscore", "e-value", "expect", "e_value",
                    default=None,
                )
                if evf is None:
                    continue
                if evf > args.max_expectation:
                    continue

            if use_hyp and has_hyp:
                hy = _get_float(row, "hyperscore", default=None)
                if hy is None or hy < args.min_hyperscore:
                    continue

            if use_ratio and has_hyp and has_next:
                hy = _get_float(row, "hyperscore", default=None)
                nx = _get_float(row, "nextscore", default=None)
                if hy is None or nx is None or nx <= 0:
                    continue
                if (hy / nx) < args.min_nextscore_ratio:
                    continue

            if exclude_decoy:
                for acc in target_proteins_only(str(protein)):
                    peptide_protein_map[peptide].add(acc)
            else:
                for acc in split_protein_identifiers(str(protein)):
                    peptide_protein_map[peptide].add(acc)

    with open(output_fasta, "w", encoding="utf-8") as fasta_out:
        for peptide, protein_set in sorted(peptide_protein_map.items()):
            if not protein_set:
                continue
            merged = "|".join(sorted(protein_set))
            fasta_out.write(f">{merged}\n{peptide}\n")

    print(f"Done; FASTA file written: {output_fasta}")


if __name__ == "__main__":
    main()
