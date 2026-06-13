"""
Screen 9/10mer peptides in a FASTA with DeepImmuno immunogenicity scores (Step 4):
remove 9/10mers below threshold; keep all other lengths unchanged.

Requires: deepimmuno-cnn.py (and data/, models/) from
https://github.com/frankligy/DeepImmuno

DeepImmuno needs peptide + HLA pairs. By default HLA is read from the first
field of the FASTA header, e.g.:
  >HLA-B*08:01|peptide=DSSNVVHLL|...
Use --hla to assign a uniform HLA to sequences missing header HLA.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# HLA-A02:01 / HLA-A0201 / HLA-A*02:01 etc.
_HLA_RE = re.compile(
    r"(HLA-[ABC]\*?\d{2}:?\d{2,4}|HLA-[ABC]\d{4,6})",
    re.IGNORECASE,
)


@dataclass
class PepRecord:
    record_id: str
    description: str
    sequence: str
    hla: str | None = None

    @property
    def length(self) -> int:
        return len(self.sequence)


def normalize_hla(hla: str) -> str:
    h = hla.strip().replace(":", "").upper()
    if not h.startswith("HLA-"):
        h = "HLA-" + h
    return h


def hla_from_header(header: str) -> str | None:
    """
    Parse HLA from FASTA header.

    Prefer the field before the first '|', e.g.:
      HLA-B*08:01|peptide=DSSNVVHLL|...  ->  HLA-B*0801
    Otherwise search the full line for HLA- patterns.
    """
    first_field = header.split("|", 1)[0].strip()
    if first_field.upper().startswith("HLA-"):
        return normalize_hla(first_field)

    m = _HLA_RE.search(header)
    if not m:
        return None
    return normalize_hla(m.group(1))


def read_fasta(path: str) -> list[PepRecord]:
    records: list[PepRecord] = []
    for rec in SeqIO.parse(path, "fasta"):
        seq = str(rec.seq).strip().upper()
        desc = rec.description or rec.id
        records.append(
            PepRecord(
                record_id=rec.id,
                description=desc,
                sequence=seq,
                hla=hla_from_header(desc),
            )
        )
    return records


def write_fasta(path: str, items: list[PepRecord]) -> None:
    seq_records = [
        SeqRecord(
            id=item.record_id,
            description=item.description,
            seq=item.sequence,
        )
        for item in items
    ]
    SeqIO.write(seq_records, path, "fasta")


def build_immuno_csv(
    items: list[PepRecord],
    default_hla: str | None,
    csv_path: str,
) -> list[PepRecord]:
    """Write DeepImmuno input CSV for 9/10mers; return entries actually used for prediction."""
    rows: list[list[str]] = []
    used: list[PepRecord] = []

    for item in items:
        if item.length not in (9, 10):
            continue
        hla = item.hla or (normalize_hla(default_hla) if default_hla else None)
        if not hla:
            print(
                f"[Warning] Skipping (no HLA): {item.record_id}  "
                f"Use --hla or annotate HLA in the header line",
                file=sys.stderr,
            )
            continue
        item.hla = hla
        rows.append([item.sequence, hla])
        used.append(item)

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(rows)

    return used


def run_deepimmuno(
    csv_path: str,
    work_dir: str,
    deepimmuno_script: str,
) -> str:
    script = os.path.abspath(deepimmuno_script)
    if not os.path.isfile(script):
        raise FileNotFoundError(f"deepimmuno-cnn.py not found: {script}")

    deepimmuno_root = os.path.dirname(script)
    cmd = [
        sys.executable,
        script,
        "--mode",
        "multiple",
        "--intdir",
        os.path.abspath(csv_path),
        "--outdir",
        os.path.abspath(work_dir),
    ]
    print(f"[DeepImmuno] {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=deepimmuno_root, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"deepimmuno-cnn.py exit code {result.returncode}")

    result_path = os.path.join(work_dir, "deepimmuno-cnn-result.txt")
    if not os.path.isfile(result_path):
        raise FileNotFoundError(
            f"Result file not generated: {result_path} (verify DeepImmuno data/ and models/ are complete)"
        )
    return result_path


def load_scores(result_path: str) -> dict[str, float]:
    """peptide -> immunogenicity score"""
    scores: dict[str, float] = {}
    with open(result_path, encoding="utf-8") as f:
        header = f.readline()
        if not header.strip():
            return scores
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                pep = cols[0].strip().upper()
                scores[pep] = float(cols[2])
    return scores


def build_output_records(
    all_records: list[PepRecord],
    scores: dict[str, float],
    predicted_sequences: set[str],
    threshold: float,
) -> list[PepRecord]:
    """Keep all non-9/10mers; keep 9/10mers only if above threshold (unpredicted ones removed)."""
    output: list[PepRecord] = []
    for rec in all_records:
        if rec.length not in (9, 10):
            output.append(rec)
        elif rec.sequence in predicted_sequences:
            if scores.get(rec.sequence, -1.0) >= threshold:
                output.append(rec)
        # 9/10mer without HLA / not predicted: omit from output
    return output


def filter_fasta(
    input_fasta: str,
    output_fasta: str,
    *,
    default_hla: str | None,
    threshold: float,
    deepimmuno_script: str,
    scores_out: str | None = None,
) -> tuple[int, int, int]:
    """
    Read FASTA → run DeepImmuno on 9/10mers → write FASTA.

    Output = all non-9/10mers + 9/10mers passing threshold (preserves input order).

    Returns:
        (input count, predicted count, output count)
    """
    all_records = read_fasta(input_fasta)
    candidates = [r for r in all_records if r.length in (9, 10)]
    other_length = [r for r in all_records if r.length not in (9, 10)]

    scores: dict[str, float] = {}
    predicted: list[PepRecord] = []

    if candidates:
        with tempfile.TemporaryDirectory(prefix="deepimmuno_") as tmp:
            csv_path = os.path.join(tmp, "input.csv")
            predicted = build_immuno_csv(candidates, default_hla, csv_path)

            if predicted:
                result_path = run_deepimmuno(csv_path, tmp, deepimmuno_script)
                scores = load_scores(result_path)

            if scores_out:
                predicted_seqs = {r.sequence for r in predicted}
                with open(scores_out, "w", newline="", encoding="utf-8") as f:
                    w = csv.writer(f, delimiter="\t")
                    w.writerow(["peptide", "HLA", "immunogenicity", "pass", "note"])
                    for item in candidates:
                        if item.sequence not in predicted_seqs:
                            w.writerow([item.sequence, "", "", "no", "missing_hla"])
                            continue
                        sc = scores.get(item.sequence, float("nan"))
                        w.writerow(
                            [
                                item.sequence,
                                item.hla,
                                sc,
                                "yes" if sc >= threshold else "no",
                                "",
                            ]
                        )
                print(f"[Scores] {scores_out}")

    predicted_sequences = {r.sequence for r in predicted}
    output = build_output_records(
        all_records, scores, predicted_sequences, threshold
    )
    write_fasta(output_fasta, output)

    passed_9_10 = sum(
        1
        for r in candidates
        if r.sequence in predicted_sequences
        and scores.get(r.sequence, -1.0) >= threshold
    )
    removed_9_10 = len(candidates) - passed_9_10

    print(
        f"[Done] Input {len(all_records)} | "
        f"Other lengths kept {len(other_length)} | "
        f"9/10mer {len(candidates)} (predicted {len(predicted)}) | "
        f"9/10mer passed {passed_9_10}, removed {removed_9_10} | "
        f"Output total {len(output)}"
    )
    print(f"[Output] {output_fasta}")
    return len(all_records), len(predicted), len(output)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run DeepImmuno on 9/10mers and remove low-scoring peptides; keep other lengths unchanged"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input FASTA path",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output FASTA path (all non-9/10mers + passing 9/10mers)",
    )
    parser.add_argument(
        "--hla",
        default=None,
        help="Uniform HLA (e.g. HLA-A0201); default parse from first header field, e.g. >HLA-B*08:01|peptide=...",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Minimum immunogenicity score, default 0.5",
    )
    parser.add_argument(
        "--deepimmuno-script",
        default="deepimmuno-cnn.py",
        help="Full path to deepimmuno-cnn.py",
    )
    parser.add_argument(
        "--scores",
        default=None,
        help="Optional: write score table for all 9/10mers (TSV)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        parser.error(f"Input file not found: {args.input}")

    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    filter_fasta(
        args.input,
        args.output,
        default_hla=args.hla,
        threshold=args.threshold,
        deepimmuno_script=args.deepimmuno_script,
        scores_out=args.scores,
    )


if __name__ == "__main__":
    main()
