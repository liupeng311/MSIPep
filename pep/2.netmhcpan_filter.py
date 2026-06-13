#!/usr/bin/env python3
"""
2.netmhcpan_filter.py
---------------------
Step 2: NetMHCpan peptide-MHC binding affinity prediction and Strong/Weak Binder (SB/WB) filtering.

Input: Step 1 blastp-filtered FASTA (peptides identical to human proteome removed)
Output:
  - {prefix}_netmhcpan.txt     merged NetMHCpan raw results (with header)
  - {prefix}_SWB.fasta          SB/WB peptide FASTA (headers include HLA, for Step 3 DeepImmuno)
  - {prefix}_SWB.tsv            SB/WB peptide table (peptide, HLA, %Rank_EL, BindLevel, etc.)

HLA typing source (choose one):
  --hla           direct specification, e.g. HLA-A*02:01,HLA-B*07:02,...
  --optitype-dir  OptiType output directory (auto-find *_result.tsv)

FASTA header format (compatible with 3.deepimmuno_filter.py):
  >HLA-A*02:01|peptide=SIINFEKL|bind=SB|rank_el=0.32|rank_ba=1.20

Usage:
  python 2.netmhcpan_filter.py -i sample_classI_mut.fasta -o sample_SWB.fasta \\
      --optitype-dir /path/to/hla_folder
  python 2.netmhcpan_filter.py -i sample_classI_mut.fasta -o sample_SWB.fasta \\
      --hla HLA-A*02:01,HLA-A*11:01,HLA-B*07:02,HLA-B*44:02,HLA-C*07:02,HLA-C*12:03
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

_BIND_RE = re.compile(r"<=\s*(SB|WB)\b")
_HLA_RE = re.compile(r"^HLA-[ABC]", re.I)


def load_hla_from_optitype(optitype_dir: str) -> str:
    """Read 6 HLA alleles from OptiType *_result.tsv."""
    folder = Path(optitype_dir)
    if not folder.is_dir():
        raise FileNotFoundError(f"OptiType directory not found: {optitype_dir}")

    candidates = sorted(folder.glob("*_result.tsv"))
    if not candidates:
        raise FileNotFoundError(f"No *_result.tsv found in {optitype_dir}")

    path = candidates[0]
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if len(lines) < 2:
        raise ValueError(f"OptiType file format error (insufficient rows): {path}")

    # Typical format: row 1 header, row 2+ typing results; take first data row A1–C2 (cols 1–6)
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) >= 7:
            alleles = [p.strip() for p in parts[1:7] if p.strip()]
            if len(alleles) >= 4 and all(_HLA_RE.match(a) for a in alleles):
                print(f"[HLA] Using OptiType result: {path.name} → {', '.join(alleles)}")
                return ",".join(alleles)

    raise ValueError(f"Unable to parse HLA from OptiType file: {path}")


def group_records_by_length(fasta_path: str) -> dict[int, list]:
    groups: dict[int, list] = defaultdict(list)
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).strip().upper()
        if not seq:
            continue
        groups[len(seq)].append(rec)
    return dict(groups)


def run_netmhcpan_per_length(
    records: list,
    length: int,
    hla_string: str,
    work_dir: Path,
    prefix: str,
    netmhcpan_bin: str,
) -> Path:
    """Run NetMHCpan for peptides of a given length; return result txt path."""
    work_dir.mkdir(parents=True, exist_ok=True)
    in_fa = work_dir / f"{prefix}_netmhcpan_in_{length}.fasta"
    out_txt = work_dir / f"{prefix}-pepBA_{length}.txt"

    SeqIO.write(records, str(in_fa), "fasta")

    cmd = [
        netmhcpan_bin,
        "-f", str(in_fa),
        "-a", hla_string,
        "-l", str(length),
        "-BA",
    ]
    print(f"[NetMHCpan] l={length} | {len(records)} peptides")
    print(f"  CMD: {' '.join(cmd)} > {out_txt}")

    with open(out_txt, "w", encoding="utf-8") as fout:
        result = subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        err = (result.stderr or "").strip()
        raise RuntimeError(f"NetMHCpan failed (l={length}, exit code {result.returncode}): {err}")

    return out_txt


def _parse_netmhcpan_line(line: str) -> dict | None:
    """Parse NetMHCpan data line; keep SB/WB only."""
    if not _BIND_RE.search(line):
        return None
    cols = line.split()
    if len(cols) < 12:
        return None

    peptide = cols[2].strip().upper()
    mhc = cols[1].strip()
    bind = "SB" if "<= SB" in line or cols[-1] == "SB" else "WB"

    def _flt(idx: int, default: float = 999.0) -> float:
        try:
            return float(cols[idx])
        except (ValueError, IndexError):
            return default

    return {
        "peptide": peptide,
        "mhc": mhc,
        "rank_el": _flt(-4),
        "rank_ba": _flt(-3),
        "aff_nm": _flt(-2),
        "bind": bind,
        "raw_line": line.rstrip("\n"),
    }


def parse_netmhcpan_outputs(txt_files: list[Path]) -> list[dict]:
    """Merge NetMHCpan outputs across lengths and parse SB/WB lines."""
    hits: list[dict] = []
    for path in txt_files:
        if not path.is_file():
            continue
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                if line.lower().startswith("pos") or "peptide" in line.lower() and "mhc" in line.lower():
                    continue
                parsed = _parse_netmhcpan_line(line)
                if parsed:
                    hits.append(parsed)
    return hits


def dedupe_best_per_peptide_mhc(hits: list[dict]) -> list[dict]:
    """For each (peptide, HLA) pair, keep the row with best (lowest) %Rank_EL."""
    best: dict[tuple[str, str], dict] = {}
    for h in hits:
        key = (h["peptide"], h["mhc"])
        if key not in best or h["rank_el"] < best[key]["rank_el"]:
            best[key] = h
    return sorted(best.values(), key=lambda x: (x["peptide"], x["rank_el"]))


def write_outputs(
    hits: list[dict],
    output_fasta: str,
    output_tsv: str | None,
    merged_txt: str | None,
    txt_sources: list[Path],
):
    out_dir = os.path.dirname(os.path.abspath(output_fasta))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Merge raw NetMHCpan text
    if merged_txt:
        header_written = False
        with open(merged_txt, "w", encoding="utf-8") as mf:
            for src in txt_sources:
                if not src.is_file():
                    continue
                with open(src, encoding="utf-8", errors="replace") as sf:
                    for line in sf:
                        mf.write(line)

    # TSV
    if output_tsv:
        with open(output_tsv, "w", encoding="utf-8", newline="\n") as tf:
            tf.write("peptide\tmhc\tbind\trank_el\trank_ba\taff_nm\n")
            for h in hits:
                tf.write(
                    f"{h['peptide']}\t{h['mhc']}\t{h['bind']}\t"
                    f"{h['rank_el']}\t{h['rank_ba']}\t{h['aff_nm']}\n"
                )

    # FASTA (headers include HLA for DeepImmuno parsing)
    with open(output_fasta, "w", encoding="utf-8", newline="\n") as ff:
        for i, h in enumerate(hits, 1):
            header = (
                f">{h['mhc']}|peptide={h['peptide']}|bind={h['bind']}"
                f"|rank_el={h['rank_el']}|rank_ba={h['rank_ba']}"
            )
            ff.write(f"{header}\n{h['peptide']}\n")

    print(f"[Output] {len(hits)} SB/WB peptides")
    print(f"  FASTA → {output_fasta}")
    if output_tsv:
        print(f"  TSV   → {output_tsv}")
    if merged_txt:
        print(f"  Raw   → {merged_txt}")


def filter_netmhcpan(
    input_fasta: str,
    output_fasta: str,
    *,
    hla: str | None = None,
    optitype_dir: str | None = None,
    netmhcpan_bin: str = "netMHCpan",
    work_dir: str | None = None,
    class1_lengths: tuple[int, ...] = (8, 9, 10, 11),
    output_tsv: str | None = None,
    merged_txt: str | None = None,
) -> int:
    if not os.path.isfile(input_fasta):
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    if hla:
        hla_string = hla.strip()
    elif optitype_dir:
        hla_string = load_hla_from_optitype(optitype_dir)
    else:
        raise ValueError("Must specify --hla or --optitype-dir")

    prefix = Path(output_fasta).stem.replace("_SWB", "").replace(".fasta", "")
    if work_dir is None:
        work_dir = str(Path(output_fasta).resolve().parent / "netmhcpan_tmp")
    tmp = Path(work_dir)

    length_groups = group_records_by_length(input_fasta)
    if not length_groups:
        print("[Warning] Input FASTA has no valid sequences")
        open(output_fasta, "w").close()
        return 0

    txt_files: list[Path] = []
    for length, records in sorted(length_groups.items()):
        if length not in class1_lengths:
            print(f"[Skip] Peptide length {length} not in Class I range {class1_lengths}")
            continue
        txt_files.append(
            run_netmhcpan_per_length(
                records, length, hla_string, tmp, prefix, netmhcpan_bin
            )
        )

    if not txt_files:
        raise RuntimeError(f"No peptide lengths in {class1_lengths}; NetMHCpan not run")

    hits = parse_netmhcpan_outputs(txt_files)
    if not hits:
        print("[Warning] No SB/WB peptides found in NetMHCpan results")
        open(output_fasta, "w").close()
        return 0

    hits = dedupe_best_per_peptide_mhc(hits)
    write_outputs(hits, output_fasta, output_tsv, merged_txt, txt_files)
    return len(hits)


def parse_args():
    p = argparse.ArgumentParser(
        description="NetMHCpan binding affinity prediction + SB/WB filter (pep pipeline Step 2)"
    )
    p.add_argument("-i", "--input", required=True, help="Input FASTA (Step 1 output)")
    p.add_argument("-o", "--output", required=True, help="Output SB/WB FASTA path")
    p.add_argument("--hla", default=None, help="HLA alleles, comma-separated")
    p.add_argument("--optitype-dir", default=None, help="OptiType result directory")
    p.add_argument("--netmhcpan", default="netMHCpan", help="netMHCpan executable")
    p.add_argument("--work-dir", default=None, help="Intermediate files directory")
    p.add_argument(
        "--lengths",
        default="8,9,10,11",
        help="Peptide lengths to predict (comma-separated, default 8,9,10,11)",
    )
    p.add_argument("--tsv", default=None, help="Optional: output SB/WB summary TSV")
    p.add_argument("--merged-txt", default=None, help="Optional: merge raw NetMHCpan output")
    return p.parse_args()


def main():
    args = parse_args()
    lengths = tuple(int(x) for x in args.lengths.split(",") if x.strip())

    out_fa = os.path.abspath(args.output)
    tsv = args.tsv or str(Path(out_fa).with_suffix(".tsv"))
    merged = args.merged_txt or str(
        Path(out_fa).parent / f"{Path(out_fa).stem.replace('_SWB', '')}_netmhcpan.txt"
    )

    n = filter_netmhcpan(
        os.path.abspath(args.input),
        out_fa,
        hla=args.hla,
        optitype_dir=args.optitype_dir,
        netmhcpan_bin=args.netmhcpan,
        work_dir=args.work_dir,
        class1_lengths=lengths,
        output_tsv=tsv,
        merged_txt=merged,
    )
    if n == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
