"""
Step 3: IEDB immunogenicity filter (Methods §2.5).

Mask positions 2, 3, and C-terminal; retain peptides with score > 0.
Runs the bundled IEDB tool from immunogenicity.zip (see start.sh).

Usage:
  python 3.iedb_filter.py -i sample_SWB.fasta -o sample_IEDB.fasta \\
      --iedb-script ./immunogenicity/new-predict_immunogenicity.py
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _mask_for_length(length: int) -> str:
    """Positions 2, 3, and C-terminal (Methods)."""
    if length < 3:
        return ""
    return f"2,3,{length}"


def run_iedb_predict(
    peptides: list[str],
    length: int,
    script: str,
    work_dir: str,
) -> dict[str, float]:
    """Return peptide -> immunogenicity score."""
    if not peptides:
        return {}

    script_path = os.path.abspath(script)
    if not os.path.isfile(script_path):
        raise FileNotFoundError(f"IEDB script not found: {script_path}")

    inp = os.path.join(work_dir, f"peptides_{length}mer.txt")
    out = os.path.join(work_dir, f"scores_{length}mer.tsv")

    with open(inp, "w", encoding="utf-8", newline="\n") as f:
        for pep in peptides:
            f.write(f"{pep}\n")

    cmd = [sys.executable, script_path, "--output", out, inp]
    mask = _mask_for_length(length)
    if mask:
        cmd.insert(-1, f"--custom_mask={mask}")

    print(f"[IEDB] length={length} n={len(peptides)} mask={mask or 'default'}")
    result = subprocess.run(cmd, cwd=os.path.dirname(script_path), check=False)
    if result.returncode != 0:
        raise RuntimeError(f"IEDB prediction failed (exit {result.returncode}) for {length}mers")

    scores: dict[str, float] = {}
    if not os.path.isfile(out):
        raise FileNotFoundError(f"IEDB output not found: {out}")

    with open(out, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.lower().startswith("peptide"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            pep = parts[0].strip().upper()
            try:
                score = float(parts[-1])
            except ValueError:
                continue
            scores[pep] = score
    return scores


def filter_fasta(
    input_fasta: str,
    output_fasta: str,
    *,
    iedb_script: str,
    min_score: float = 0.0,
    scores_out: str | None = None,
) -> tuple[int, int]:
    """
    Filter FASTA: keep peptides with IEDB score > min_score (Methods: > 0).
    """
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        SeqIO.write([], output_fasta, "fasta")
        return 0, 0

    by_len: dict[int, list[SeqRecord]] = defaultdict(list)
    for rec in records:
        seq = str(rec.seq).strip().upper()
        if seq:
            by_len[len(seq)].append(
                SeqRecord(rec.seq, id=rec.id, description=rec.description)
            )

    all_scores: dict[str, float] = {}
    with tempfile.TemporaryDirectory(prefix="iedb_") as tmp:
        for length, recs in sorted(by_len.items()):
            unique_peps = sorted({str(r.seq).upper() for r in recs})
            all_scores.update(run_iedb_predict(unique_peps, length, iedb_script, tmp))

    kept: list[SeqRecord] = []
    for rec in records:
        seq = str(rec.seq).strip().upper()
        score = all_scores.get(seq)
        if score is not None and score > min_score:
            kept.append(rec)

    out_dir = os.path.dirname(os.path.abspath(output_fasta))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    SeqIO.write(kept, output_fasta, "fasta")

    if scores_out:
        with open(scores_out, "w", encoding="utf-8", newline="\n") as f:
            f.write("peptide\tlength\tscore\tpass\n")
            for rec in records:
                seq = str(rec.seq).strip().upper()
                sc = all_scores.get(seq)
                passed = sc is not None and sc > min_score
                f.write(f"{seq}\t{len(seq)}\t{sc if sc is not None else ''}\t{passed}\n")
        print(f"[Scores] {scores_out}")

    print(
        f"[Done] Input {len(records)} | kept {len(kept)} (score > {min_score}) | "
        f"removed {len(records) - len(kept)}"
    )
    print(f"[Output] {output_fasta}")
    return len(records), len(kept)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="IEDB immunogenicity filter (mask 2,3,C-term; score > threshold)"
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA (e.g. NetMHCpan SB/WB)")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA")
    parser.add_argument(
        "--iedb-script",
        default="./immunogenicity/new-predict_immunogenicity.py",
        help="Path to IEDB predict_immunogenicity script",
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=0.0,
        help="Keep peptides with score strictly greater than this (Methods: > 0)",
    )
    parser.add_argument("--scores", default=None, help="Optional TSV of all scores")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        parser.error(f"Input not found: {args.input}")

    filter_fasta(
        args.input,
        args.output,
        iedb_script=args.iedb_script,
        min_score=args.min_score,
        scores_out=args.scores,
    )


if __name__ == "__main__":
    main()
