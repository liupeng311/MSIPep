#!/usr/bin/env python3
"""
7.protcr_filter.py
------------------
Step 7: Filter peptides with plabels ≥ threshold from ProTCR output CSVs.

Usage:
  python 7.protcr_filter.py /path/to/protcr_results --threshold 0.9 -o filtered/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def filter_protcr_folder(
    folder: Path,
    output_dir: Path,
    threshold: float = 0.9,
    recursive: bool = True,
    output_format: str = "csv",
) -> tuple[int, int]:
    pattern = "**/*.csv" if recursive else "*.csv"
    csv_files = list(folder.glob(pattern))
    if not csv_files:
        print(f"❌ No CSV files found in {folder}")
        return 0, 0

    output_dir.mkdir(parents=True, exist_ok=True)
    total_count = 0
    processed = 0

    for csv_file in csv_files:
        if csv_file.name.endswith("_tcr_input.csv"):
            continue
        if "_filtered" in csv_file.stem:
            continue

        try:
            df = pd.read_csv(csv_file, sep=None, engine="python", encoding="utf-8")
            df.columns = [col.strip() for col in df.columns]
            lower_map = {c.lower(): c for c in df.columns}
            pep_col = lower_map.get("peptide")
            plab_col = lower_map.get("plabels")
            if not pep_col or not plab_col:
                print(f"⚠️  Skipping (missing peptide/plabels): {csv_file.name}")
                continue

            df[plab_col] = pd.to_numeric(df[plab_col], errors="coerce")
            filtered = df[df[plab_col] >= threshold].copy()
            if filtered.empty:
                print(f"ℹ️  No peptides with plabels≥{threshold}: {csv_file.name}")
                continue

            count = len(filtered)
            total_count += count
            processed += 1
            out_base = f"{csv_file.stem}_filtered"

            if output_format == "xlsx":
                out_path = output_dir / f"{out_base}.xlsx"
                filtered[[pep_col, plab_col]].to_excel(out_path, index=False)
            else:
                out_path = output_dir / f"{out_base}.csv"
                filtered[[pep_col, plab_col]].to_csv(out_path, index=False)

            print(f"✅ {csv_file.name} → {count} rows → {out_path.name}")

        except Exception as e:
            print(f"❌ Failed to process {csv_file.name}: {e}")

    return processed, total_count


def main():
    parser = argparse.ArgumentParser(
        description="Filter ProTCR results: extract peptides with plabels ≥ threshold"
    )
    parser.add_argument("folder_path", help="Directory containing ProTCR result CSVs")
    parser.add_argument("--output_dir", default="filtered_results", help="Output directory for filtered results")
    parser.add_argument("--threshold", type=float, default=0.9, help="Minimum plabels value")
    parser.add_argument("--recursive", action="store_true", help="Search subdirectories recursively")
    parser.add_argument("--format", choices=("csv", "xlsx"), default="csv")
    args = parser.parse_args()

    folder = Path(args.folder_path)
    if not folder.is_dir():
        print(f"❌ Directory not found: {folder}")
        return

    print(f"Filtering directory: {folder} | plabels ≥ {args.threshold}")
    print("=" * 60)

    n_files, n_pep = filter_protcr_folder(
        folder,
        Path(args.output_dir),
        threshold=args.threshold,
        recursive=args.recursive,
        output_format=args.format,
    )

    print("=" * 60)
    print(f"🎉 Done: processed {n_files} file(s), {n_pep} peptide(s) total")
    print(f"   Output directory: {Path(args.output_dir).resolve()}")


if __name__ == "__main__":
    main()
