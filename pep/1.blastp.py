"""
pep module — immunopeptide filtering pipeline (Methods §2.5)

  Step 1  1.blastp.py              Remove 100% human proteome matches
  Step 2  2.netmhcpan_filter.py    NetMHCpan SB/WB (%Rank ≤0.5 / ≤2)
  Step 3  3.iedb_filter.py         IEDB immunogenicity (mask 2,3,C-term; score > 0)
  Step 4  4.deepimmuno_filter.py   DeepImmuno 9/10mer (score ≥ 0.5)
  Step 5  5.protcr_input.py        ProTCR input CSV
  Step 6  6.protcr_run.py          ProTCR prediction
  Step 7  7.protcr_filter.py       plabels ≥ 0.9

Run all:  python run_pep_pipeline.py --config pep_config.yaml
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO

DEFAULT_INPUT_FASTA = ""
DEFAULT_BLAST_DB = "./software/ncbi-blast/test/human_uniprot"
BATCH_INPUT_SUFFIX = "_merged_classI_dedup.fasta"


def _resolve_blastp_executable(db_path: str | None, explicit: str | None) -> str:
    """Resolve blastp: --blastp → PATH → bin/ next to --db."""
    if explicit:
        exp = os.path.abspath(os.path.expanduser(explicit))
        if os.path.isfile(exp):
            return exp
        if os.name == "nt" and not exp.lower().endswith(".exe"):
            exp_exe = exp + ".exe"
            if os.path.isfile(exp_exe):
                return exp_exe
        raise FileNotFoundError(f"blastp executable not found: {explicit} (also tried .exe)")

    for name in ("blastp", "blastp.exe"):
        w = shutil.which(name)
        if w:
            return w

    if db_path:
        cur = os.path.abspath(os.path.dirname(db_path))
        for _ in range(10):
            for leaf in ("blastp.exe", "blastp"):
                cand = os.path.join(cur, "bin", leaf)
                if os.path.isfile(cand):
                    return os.path.abspath(cand)
            parent = os.path.dirname(cur)
            if parent == cur:
                break
            cur = parent

    raise FileNotFoundError(
        "blastp not found. Add NCBI BLAST bin to PATH or pass --blastp, e.g. "
        r'--blastp "D:\Blast\ncbi-blast-2.15.0+\bin\blastp.exe"'
    )


def _run_blastp(blastp_exe: str, query_fasta: str, db_path: str, blast_out: str):
    with open(blast_out, "w", encoding="utf-8", newline="\n") as out_f:
        return subprocess.run(
            [
                blastp_exe,
                "-query",
                query_fasta,
                "-db",
                db_path,
                "-outfmt",
                "6 qacc sacc qlen length pident qseq sseq",
                "-evalue",
                "1e9",
                "-gapopen",
                "11",
                "-gapextend",
                "1",
            ],
            stdout=out_f,
            stderr=subprocess.PIPE,
            text=True,
        )


def _load_query_records(query_fasta: str):
    """Load FASTA; key = record.id (matches blast qacc)."""
    records = {}
    qlen_map = {}
    for rec in SeqIO.parse(query_fasta, "fasta"):
        seq_str = str(rec.seq).replace("-", "").replace("*", "").upper()
        qlen_map[rec.id] = len(seq_str)
        records[rec.id] = rec
    return records, qlen_map


def _qacc_has_full_alignment(blast_tsv: str, qlen_map: dict) -> set:
    """Full match: no gaps, pident==100, alignment length equals query length."""
    fully = set()
    if not os.path.exists(blast_tsv) or os.path.getsize(blast_tsv) == 0:
        return fully

    with open(blast_tsv, encoding="utf-8") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            qacc, qlen_s, length_s, pident_s, qseq, sseq = (
                cols[0],
                cols[2],
                cols[3],
                cols[4],
                cols[5],
                cols[6],
            )
            if "-" in qseq or "-" in sseq:
                continue
            try:
                blast_qlen = int(qlen_s)
                aln_len = int(length_s)
                pident_f = float(pident_s)
            except ValueError:
                continue
            ref_len = qlen_map.get(qacc)
            if ref_len is None:
                continue
            if blast_qlen != ref_len or aln_len != ref_len:
                continue
            if round(pident_f, 2) == 100.0:
                fully.add(qacc)
    return fully


def filter_fasta_remove_full_matches(
    input_fasta: str,
    output_fasta: str,
    db_path: str,
    blastp_exe: str,
    keep_blast_tsv=None,
):
    records, qlen_map = _load_query_records(input_fasta)
    if not records:
        print("Input FASTA contains no sequences.")
        return 0

    tmpdir = None
    if keep_blast_tsv:
        blast_txt = os.path.abspath(keep_blast_tsv)
    else:
        tmpdir = tempfile.mkdtemp(prefix="blastp_pep_")
        blast_txt = os.path.join(tmpdir, "blastp_result.tsv")

    try:
        result = _run_blastp(blastp_exe, input_fasta, db_path, blast_txt)
        if result.returncode != 0:
            err = (result.stderr or "").strip()
            raise RuntimeError(
                f"blastp failed (exit {result.returncode})." + (f" stderr: {err}" if err else "")
            )

        remove_ids = _qacc_has_full_alignment(blast_txt, qlen_map)
        kept = [records[qid] for qid in records if qid not in remove_ids]

        out_dir = os.path.dirname(os.path.abspath(output_fasta))
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        SeqIO.write(kept, output_fasta, "fasta")

        n_in = len(records)
        n_out = len(kept)
        print(
            f"Input {n_in}; removed {len(remove_ids)} full matches; "
            f"wrote {n_out} to {output_fasta}"
        )
        return n_out
    finally:
        if tmpdir is not None:
            shutil.rmtree(tmpdir, ignore_errors=True)


def run_batch(batch_dir: str, db_path: str, blastp_exe: str | None, keep_blast: bool):
    """Process all *_merged_classI_dedup.fasta files in a directory."""
    pattern = os.path.join(batch_dir, f"*{BATCH_INPUT_SUFFIX}")
    fasta_files = sorted(glob.glob(pattern))
    if not fasta_files:
        print(f"No files matching *{BATCH_INPUT_SUFFIX} in {batch_dir}")
        return

    bp = _resolve_blastp_executable(os.path.abspath(db_path), blastp_exe)
    print(f"Found {len(fasta_files)} files; starting batch blastp...\n")

    for file_path in fasta_files:
        sample = os.path.basename(file_path).replace(BATCH_INPUT_SUFFIX, "")
        out_file = os.path.join(batch_dir, f"{sample}_classI_mut.fasta")
        tsv_file = os.path.join(batch_dir, f"{sample}_classI_mut.tsv") if keep_blast else None
        print(f"Sample: {sample}")
        filter_fasta_remove_full_matches(
            file_path,
            out_file,
            os.path.abspath(db_path),
            bp,
            tsv_file,
        )
        print("=" * 60)

    print("Batch processing complete.")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "blastp against a reference proteome; remove query sequences with "
            "gap-free 100% full-length matches to the database."
        )
    )
    parser.add_argument("--input", "-i", default=None, help="Query FASTA (single file)")
    parser.add_argument(
        "--batch-dir",
        "-d",
        default=None,
        help=f"Batch mode: process all *{BATCH_INPUT_SUFFIX} in this directory",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output FASTA (default: <input>_nofullhit.fasta)",
    )
    parser.add_argument(
        "--db",
        dest="db_path",
        default=DEFAULT_BLAST_DB,
        help="blastp -db prefix (no extension)",
    )
    parser.add_argument(
        "--blastp",
        default=None,
        metavar="EXE",
        help="Path to blastp executable (required on Windows if not in PATH)",
    )
    parser.add_argument(
        "--keep-blast",
        metavar="PATH",
        default=None,
        help="Save blast tabular output to this path",
    )
    parser.add_argument(
        "--keep-blast-tsv",
        action="store_true",
        help="In batch mode, save {sample}_classI_mut.tsv per sample",
    )

    args = parser.parse_args()
    db = args.db_path or DEFAULT_BLAST_DB

    if args.batch_dir:
        run_batch(args.batch_dir, db, args.blastp, args.keep_blast or args.keep_blast_tsv)
        return

    inp = args.input or DEFAULT_INPUT_FASTA
    if not inp:
        parser.error("Provide --input/-i or --batch-dir/-d")
    inp = os.path.abspath(inp)
    if not os.path.isfile(inp):
        raise SystemExit(f"Input file not found: {inp}")

    out_fa = args.output
    if not out_fa:
        stem, ext = os.path.splitext(inp)
        out_fa = stem + "_nofullhit" + (ext if ext else ".fasta")

    bp = _resolve_blastp_executable(os.path.abspath(db), args.blastp)
    filter_fasta_remove_full_matches(
        inp,
        os.path.abspath(out_fa),
        os.path.abspath(db),
        bp,
        args.keep_blast,
    )
    print("Done.")


if __name__ == "__main__":
    main()
