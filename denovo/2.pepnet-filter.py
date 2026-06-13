"""
PepNet output TSV → FASTA (length + confidence filtering).

Two-stage filtering (Methods):
- Length filter: 8–11 aa
- Confidence filter: Score ≥ 0.85 and |PPM Difference| ≤ 10 ppm

Usage examples:
  python 2.pepnet-filter.py -i sample.tsv -o sample-i.fasta
  python 2.pepnet-filter.py -i /path/to/denovo/*.tsv -o merged.fasta
  python 2.pepnet-filter.py -i denovo_dir/ --per-file --output-dir out/
  python 2.pepnet-filter.py -i sample.tsv -o out.fasta --length-only
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import pandas as pd


def peptide_aa_length(peptide: str) -> int:
    if not peptide:
        return 0
    return sum(1 for c in peptide if c.isalpha())


def in_class1_length_range(peptide: str, min_len: int, max_len: int) -> bool:
    ln = peptide_aa_length(peptide.strip())
    return min_len <= ln <= max_len


def detect_encoding(file_path: str) -> str:
    try:
        import chardet

        with open(file_path, "rb") as f:
            result = chardet.detect(f.read())
        enc = result.get("encoding")
        if enc:
            return enc
    except ImportError:
        pass
    return "utf-8"


def read_tsv(file_path: str) -> pd.DataFrame:
    encoding = detect_encoding(file_path)
    try:
        return pd.read_csv(file_path, sep="\t", encoding=encoding, dtype=str)
    except UnicodeDecodeError:
        return pd.read_csv(file_path, sep="\t", encoding="utf-8", errors="replace", dtype=str)


def _norm_col(c: str) -> str:
    return str(c).strip().lower().replace(" ", "_")


def _rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [_norm_col(c) for c in df.columns]
    return df


def _to_float(val, default=None):
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return default
    s = str(val).strip()
    if not s or s.lower() in ("nan", "none"):
        return default
    try:
        return float(s)
    except (TypeError, ValueError):
        return default


def parse_positional_scores(raw) -> list[float]:
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return []
    if isinstance(raw, (list, tuple)):
        out = []
        for x in raw:
            v = _to_float(x)
            if v is not None:
                out.append(v)
        return out

    s = str(raw).strip()
    if not s or s.lower() in ("nan", "none", "[]"):
        return []
    if s.startswith("[") and s.endswith("]"):
        s = s[1:-1]
    out = []
    for part in s.split(","):
        v = _to_float(part.strip())
        if v is not None:
            out.append(v)
    return out


def _avg_positional_score(raw) -> float:
    scores = parse_positional_scores(raw)
    if not scores:
        return 0.0
    return sum(scores) / len(scores)


def _make_header(row: pd.Series, file_stem: str, row_idx: int) -> str:
    title = row.get("title")
    if title is not None and str(title).strip() and str(title).strip().lower() != "nan":
        return f"{str(title).strip()}_{row_idx + 1}"
    return f"{file_stem}_{row_idx + 1}"


def filter_pepnet_rows(
    df: pd.DataFrame,
    file_stem: str,
    *,
    min_length: int = 8,
    max_length: int = 11,
    min_score: float = 0.85,
    max_ppm_difference: float = 10.0,
    min_avg_positional_score: float = 0.0,
    length_only: bool = False,
) -> list[tuple[str, str]]:
    df = _rename_columns(df)
    if "denovo" not in df.columns:
        raise ValueError("TSV missing DENOVO column")

    required_score_cols = ("score", "ppm_difference")
    use_positional = (not length_only) and (min_avg_positional_score > 0)
    if not length_only:
        missing = [c for c in required_score_cols if c not in df.columns]
        if missing:
            raise ValueError(
                "Confidence filtering requires columns: Score, PPM Difference; "
                f"missing: {', '.join(missing)}. "
                "Use --length-only for length-only filtering"
            )
        if use_positional and "positional_score" not in df.columns:
            raise ValueError(
                "Positional Score filtering enabled but column missing; "
                "set --min_avg_positional_score 0 to disable"
            )

    records: list[tuple[str, str]] = []
    for i, row in df.iterrows():
        peptide = row.get("denovo")
        if peptide is None or str(peptide).strip() == "" or str(peptide).lower() == "nan":
            continue
        peptide = str(peptide).strip()

        if not in_class1_length_range(peptide, min_length, max_length):
            continue

        if not length_only:
            score = _to_float(row.get("score"), default=-1.0)
            ppm = _to_float(row.get("ppm_difference"), default=1e9)

            if score < min_score:
                continue
            if abs(ppm) > max_ppm_difference:
                continue
            if use_positional:
                avg_pos = _avg_positional_score(row.get("positional_score"))
                if avg_pos < min_avg_positional_score:
                    continue

        records.append((_make_header(row, file_stem, int(i)), peptide))

    return records


def write_fasta(records: list[tuple[str, str]], output_file: str) -> int:
    out_path = os.path.abspath(output_file)
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")
    return len(records)


def collect_tsv_paths(inputs: list[str]) -> list[str]:
    paths: list[str] = []
    seen: set[str] = set()

    for item in inputs:
        if any(ch in item for ch in "*?[]"):
            matched = sorted(glob.glob(item))
            if not matched:
                raise FileNotFoundError(f"Glob matched no files: {item}")
            for p in matched:
                ap = os.path.abspath(p)
                if os.path.isfile(ap) and ap not in seen:
                    seen.add(ap)
                    paths.append(ap)
            continue

        p = os.path.abspath(item)
        if os.path.isfile(p):
            if p not in seen:
                seen.add(p)
                paths.append(p)
        elif os.path.isdir(p):
            for name in sorted(os.listdir(p)):
                if name.lower().endswith(".tsv"):
                    fp = os.path.join(p, name)
                    if fp not in seen:
                        seen.add(fp)
                        paths.append(fp)
        else:
            raise FileNotFoundError(f"Path not found: {item}")

    if not paths:
        raise FileNotFoundError("No .tsv input files found")
    return paths


def parse_arguments():
    p = argparse.ArgumentParser(
        description="Filter PepNet TSV (length + confidence) and export FASTA"
    )
    p.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="PepNet output .tsv, directory, or glob (multiple allowed)",
    )
    p.add_argument(
        "-o",
        "--output",
        help="Merged output FASTA (mutually exclusive with --per-file)",
    )
    p.add_argument(
        "--per-file",
        action="store_true",
        help="Write one {filename}-i.fasta per input TSV",
    )
    p.add_argument(
        "--output-dir",
        help="With --per-file: specify output directory (default: same as input TSV)",
    )
    p.add_argument("--min_length", type=int, default=8, help="Minimum peptide length (default 8)")
    p.add_argument("--max_length", type=int, default=11, help="Maximum peptide length (default 11; legacy 2.1 script used 12)")
    p.add_argument(
        "--length-only",
        action="store_true",
        help="Length filter only (skip Score/PPM/Positional filtering)",
    )
    p.add_argument(
        "--min_score",
        type=float,
        default=0.85,
        help="Minimum Score (Methods: 0.85)",
    )
    p.add_argument(
        "--max_ppm_difference",
        type=float,
        default=10.0,
        help="Maximum absolute PPM Difference in ppm (Methods: 10)",
    )
    p.add_argument(
        "--min_avg_positional_score",
        type=float,
        default=0.0,
        help="Minimum average Positional Score (0 = disabled; not in Methods)",
    )
    return p.parse_args()


def main():
    args = parse_arguments()
    if args.per_file and args.output:
        print("⚠️ Both --output and --per-file specified: --output will be ignored; use --output-dir or input directory.")
    if not args.per_file and not args.output:
        raise SystemExit("Specify -o/--output, or use --per-file")

    tsv_files = collect_tsv_paths(args.input)
    filter_kwargs = dict(
        min_length=args.min_length,
        max_length=args.max_length,
        min_score=args.min_score,
        max_ppm_difference=args.max_ppm_difference,
        min_avg_positional_score=args.min_avg_positional_score,
        length_only=args.length_only,
    )

    mode = "length only" if args.length_only else "length + confidence"
    print(f"{len(tsv_files)} PepNet TSV file(s) | Mode: {mode}")
    if not args.length_only:
        pos_msg = (
            f"Positional avg≥{args.min_avg_positional_score}; "
            if args.min_avg_positional_score > 0
            else ""
        )
        print(
            f"  Score≥{args.min_score}; |PPM|≤{args.max_ppm_difference}; "
            f"{pos_msg}"
            f"peptide length {args.min_length}-{args.max_length} aa"
        )
    else:
        print(f"  peptide length {args.min_length}-{args.max_length} aa")

    total_written = 0

    if args.per_file:
        for tsv_file in tsv_files:
            stem = os.path.splitext(os.path.basename(tsv_file))[0]
            if args.output_dir:
                out_path = os.path.join(os.path.abspath(args.output_dir), f"{stem}-i.fasta")
            else:
                out_path = os.path.join(os.path.dirname(tsv_file), f"{stem}-i.fasta")

            df = read_tsv(tsv_file)
            records = filter_pepnet_rows(df, stem, **filter_kwargs)
            n = write_fasta(records, out_path)
            total_written += n
            print(f"✅ {os.path.basename(tsv_file)} → {n} entries → {out_path}")
    else:
        all_records: list[tuple[str, str]] = []
        for tsv_file in tsv_files:
            stem = os.path.splitext(os.path.basename(tsv_file))[0]
            df = read_tsv(tsv_file)
            records = filter_pepnet_rows(df, stem, **filter_kwargs)
            all_records.extend(records)
            print(f"  {os.path.basename(tsv_file)}: {len(records)} passed")

        n = write_fasta(all_records, args.output)
        total_written = n
        print(f"✅ Wrote {n} entries total → {os.path.abspath(args.output)}")

    if total_written == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
