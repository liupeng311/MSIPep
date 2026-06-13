#!/usr/bin/env python3
"""
5.coding_mut_pep.py
-------------------
Coding-region VCF → mutant peptide FASTA (reimplements pVACtools pVACseq FastaGenerator)

Workflow:
  1. Parse VEP-annotated VCF (requires WildtypeProtein, FrameshiftSequence fields)
  2. Filter missense / inframe / frameshift variants per pVACtools VcfConverter rules
  3. Generate WT / MT protein subsequences (flanking + downstream frameshift logic matches pVACtools)
  4. Output FASTA with headers compatible with downstream peptide annotation scripts:
     >MT.{index}.{gene}.{ENST}.{MS|FS|IF}.{amino_acid_change}
     >WT.{index}.{gene}.{ENST}.{MS|FS|IF}.{amino_acid_change}

Dependencies: Python standard library only (gzip optional)

Usage:
  python 5.coding_mut_pep.py coding.vcf output_prefix \\
      [--flanking-length 24] [--downstream-length 0] [--epitope-length 0]
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass
from typing import Iterator


SUPPORTED_AAS = set("ARNDCEQGHILKMFPSTWYV")
INVALID_AAS = {"*", "X", "?"}

VARIANT_CLASS_MAP = {
    "missense": "MS",
    "FS": "FS",
    "inframe_ins": "IF",
    "inframe_del": "IF",
}


# ---------------------------------------------------------------------------
# VCF / CSQ parsing
# ---------------------------------------------------------------------------

def open_text(path: str, mode: str = "r"):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def extract_csq_fields(header_lines: list[str]) -> list[str]:
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            marker = "Format: "
            if marker in line:
                tail = line.split(marker, 1)[1].rstrip('">')
                return [x.strip() for x in tail.split("|")]
    return []


def parse_info(info_str: str) -> dict:
    info: dict = {}
    if info_str == ".":
        return info
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def is_insertion(ref: str, alt: str) -> bool:
    return len(alt) > len(ref)


def is_deletion(ref: str, alt: str) -> bool:
    return len(alt) < len(ref)


def resolve_alleles(ref: str, alt_field: str) -> dict[str, str]:
    alleles: dict[str, str] = {}
    for alt in alt_field.split(","):
        alt = alt.strip()
        if not alt or alt == ".":
            continue
        if is_insertion(ref, alt) or is_deletion(ref, alt):
            if alt[0:1] != ref[0:1]:
                alleles[alt] = alt
            elif alt[1:] == "":
                alleles[alt] = "-"
            else:
                alleles[alt] = alt[1:]
        else:
            alleles[alt] = alt
    return alleles


def parse_csq_for_allele(csq_value: str, csq_fields: list[str], allele: str) -> list[dict]:
    transcripts = []
    for entry in csq_value.split(","):
        parts = entry.split("|")
        record = {csq_fields[i]: (parts[i] if i < len(parts) else "") for i in range(len(csq_fields))}
        if record.get("Allele") == allele:
            transcripts.append(record)
    return transcripts


def resolve_consequence(consequence_string: str, ref: str, alt: str) -> str | None:
    if "&" in consequence_string:
        consequences = {c.lower() for c in consequence_string.split("&")}
    elif "." in consequence_string:
        consequences = {c.lower() for c in consequence_string.split(".")}
    else:
        consequences = {consequence_string.lower()}

    if "start_lost" in consequences or "stop_retained_variant" in consequences:
        return None
    if "frameshift_variant" in consequences:
        return "FS"
    if "missense_variant" in consequences:
        return "missense"
    if "inframe_insertion" in consequences:
        return "inframe_ins"
    if "inframe_deletion" in consequences:
        return "inframe_del"
    if "protein_altering_variant" in consequences:
        if len(ref) > len(alt) and (len(ref) - len(alt)) % 3 == 0:
            return "inframe_del"
        if len(alt) > len(ref) and (len(alt) - len(ref)) % 3 == 0:
            return "inframe_ins"
    return None


def decode_hex_url(s: str) -> str:
    return re.sub(
        r"%[0-9A-Fa-f]{2}",
        lambda m: bytes.fromhex(m.group(0)[1:]).decode("utf-8"),
        s,
    )


def parse_protein_position(raw: str) -> str | None:
    if "/" in raw:
        pos = raw.split("/")[0]
        if pos == "-":
            pos = raw.split("/")[1]
    else:
        pos = raw
    if pos in ("", "-", "."):
        return None
    return pos


def genotype_has_alt(fmt: str, sample: str, alt: str) -> bool:
    if not sample or sample == ".":
        return True
    fmt_keys = fmt.split(":")
    sample_vals = sample.split(":")
    fmt_map = dict(zip(fmt_keys, sample_vals))
    gt = fmt_map.get("GT", "")
    if gt in ("", ".", "./."):
        return True
    alleles = re.split(r"[|/]", gt)
    alt_idx = None
    # In GT, 1 means first ALT
    for i, a in enumerate(alleles):
        if a not in ("0", "0/0", "."):
            alt_idx = a
            break
    return alt_idx is not None


# ---------------------------------------------------------------------------
# pVACtools FastaGenerator core algorithm
# ---------------------------------------------------------------------------

@dataclass
class VariantRecord:
    index: str
    variant_type: str
    gene: str
    transcript: str
    protein_position: str
    amino_acid_change: str
    wildtype_sequence: str
    mutant_sequence: str


def distance_from_start(position: int) -> int:
    return position


def distance_from_end(position: int, sequence: str) -> int:
    return len(sequence) - 1 - position


def get_wildtype_subsequence(
    position: int,
    full_wt: str,
    wt_aa_len: int,
    flanking: int,
) -> tuple[int, str]:
    peptide_len = min(2 * flanking + wt_aa_len, len(full_wt))
    if distance_from_start(position) < flanking:
        subseq = full_wt[:peptide_len]
        mut_pos = position
    elif distance_from_end(position, full_wt) < flanking:
        start = len(full_wt) - peptide_len
        subseq = full_wt[start:]
        mut_pos = peptide_len - distance_from_end(position, full_wt) - 1
    else:
        start = position - flanking
        subseq = full_wt[start : start + peptide_len]
        mut_pos = flanking
    return mut_pos, subseq


def get_frameshift_subsequences(
    position: int,
    full_wt: str,
    full_mt: str,
    flanking: int,
    downstream_len: int | None,
) -> tuple[str, str]:
    start = 0 if position < flanking else position - flanking
    wt_stop = position + flanking
    wt_sub = full_wt[start:wt_stop]
    if downstream_len:
        mt_sub = full_mt[start : position + downstream_len]
    else:
        mt_sub = full_mt[start:]
    return wt_sub, mt_sub


def trim_unsupported_ends(sequence: str) -> str:
    seq = sequence
    if seq and seq[0] not in SUPPORTED_AAS:
        seq = seq[1:]
    if seq and seq[-1] not in SUPPORTED_AAS:
        seq = seq[:-1]
    return seq


def contains_invalid(sequence: str) -> bool:
    return any(c in INVALID_AAS for c in sequence) or not all(c in SUPPORTED_AAS for c in sequence)


def determine_neoepitopes(sequence: str, length: int) -> dict[int, str]:
    return {i + 1: sequence[i : i + length] for i in range(max(0, len(sequence) - length + 1))}


def has_novel_epitope(wt_sub: str, mt_sub: str, epitope_len: int) -> bool:
    if epitope_len <= 0:
        return mt_sub not in wt_sub
    neo = determine_neoepitopes(mt_sub, epitope_len)
    return not all(pep in wt_sub for pep in neo.values())


def build_mutant_sequence(
    variant_type: str,
    wt_full: str,
    mt_full_fs: str,
    protein_position: str,
    amino_acids: str,
    flanking: int,
    downstream_len: int | None,
    epitope_len: int,
) -> tuple[str, str] | None:
    if variant_type == "FS":
        position = int(protein_position.split("-", 1)[0]) - 1
        if mt_full_fs.startswith(wt_full):
            return None
        wt_sub, mt_sub = get_frameshift_subsequences(position, wt_full, mt_full_fs, flanking, downstream_len)
    else:
        if "/" not in amino_acids:
            return None
        wt_aa, mt_aa = amino_acids.split("/", 1)
        stop_added = False
        for tag in ("*", "X"):
            if tag in wt_aa:
                wt_aa = wt_aa.split(tag)[0]
            if tag in mt_aa:
                mt_aa = mt_aa.split(tag)[0]
                stop_added = True
        if mt_aa == "-":
            mt_aa = ""

        if wt_aa == "-":
            position = int(protein_position.split("-", 1)[0])
            wt_aa_len = 0
        elif "-" in protein_position:
            position = int(protein_position.split("-", 1)[0]) - 1
            wt_aa_len = len(wt_aa)
        else:
            position = int(protein_position) - 1
            wt_aa_len = len(wt_aa)

        if position > len(wt_full) - 1:
            return None

        if variant_type == "inframe_ins":
            mut_start, wt_sub = get_wildtype_subsequence(position, wt_full, len(mt_aa), flanking)
        else:
            mut_start, wt_sub = get_wildtype_subsequence(position, wt_full, wt_aa_len, flanking)

        mut_end = mut_start + wt_aa_len
        if wt_aa != "-" and wt_sub[mut_start:mut_end] != wt_aa:
            return None

        if stop_added:
            mt_sub = wt_sub[:mut_start] + mt_aa
        else:
            mt_sub = wt_sub[:mut_start] + mt_aa + wt_sub[mut_end:]

    wt_sub = trim_unsupported_ends(wt_sub)
    mt_sub = trim_unsupported_ends(mt_sub)
    if not wt_sub or not mt_sub or contains_invalid(wt_sub) or contains_invalid(mt_sub):
        return None

    min_len = epitope_len if epitope_len > 0 else 1
    if len(wt_sub) < min_len or len(mt_sub) < min_len:
        return None
    if not has_novel_epitope(wt_sub, mt_sub, epitope_len):
        return None
    return wt_sub, mt_sub


def construct_index(count: int, gene: str, transcript: str, vclass: str, change: str) -> str:
    return f"{count}.{gene}.{transcript}.{vclass}.{change}"


def transcript_priority(tx: dict) -> tuple:
    canonical = 0 if tx.get("CANONICAL") == "YES" else 1
    tsl_raw = tx.get("TSL", "99") or "99"
    try:
        tsl = int(tsl_raw)
    except ValueError:
        tsl = 99
    return (canonical, tsl)


def iter_coding_variants(
    vcf_path: str,
    pass_only: bool = True,
    canonical_only: bool = True,
    biotypes: tuple[str, ...] = ("protein_coding",),
) -> Iterator[dict]:
    headers: list[str] = []
    csq_fields: list[str] = []
    count = 1
    seen_indexes: set[str] = set()

    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                headers.append(line.rstrip("\n"))
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            if not csq_fields:
                csq_fields = extract_csq_fields(headers)
                required = {"WildtypeProtein", "FrameshiftSequence", "Consequence", "Feature"}
                missing = required - set(csq_fields)
                if missing:
                    sys.exit(
                        f"❌ VCF CSQ missing required pVACtools fields: {missing}\n"
                        "Re-annotate with VEP + Wildtype/Frameshift plugins."
                    )

            chrom, pos, _id, ref, alt_field, flt, fmt, info_str = (
                cols[0], cols[1], cols[2], cols[3], cols[4], cols[6], cols[8], cols[7]
            )
            sample = cols[9] if len(cols) > 9 else "."

            if pass_only and flt != "PASS":
                continue

            info = parse_info(info_str)
            if "CSQ" not in info:
                continue

            allele_map = resolve_alleles(ref, alt_field)
            for alt, csq_allele in allele_map.items():
                if not genotype_has_alt(fmt, sample, alt):
                    continue

                transcripts = parse_csq_for_allele(info["CSQ"], csq_fields, csq_allele)
                if not transcripts and is_deletion(ref, alt):
                    transcripts = parse_csq_for_allele(info["CSQ"], csq_fields, "deletion")
                if not transcripts:
                    transcripts = parse_csq_for_allele(info["CSQ"], csq_fields, alt)

                if canonical_only:
                    canonical = [t for t in transcripts if t.get("CANONICAL") == "YES"]
                    pool = canonical if canonical else transcripts
                else:
                    pool = transcripts
                pool = sorted(pool, key=transcript_priority)

                for tx in pool:
                    flags = (tx.get("FLAGS") or "").lower()
                    if "cds_start_nf" in flags or "cds_end_nf" in flags:
                        continue

                    biotype = tx.get("BIOTYPE") or ""
                    if biotype and biotype not in biotypes:
                        continue

                    vtype = resolve_consequence(tx.get("Consequence", ""), ref, alt)
                    if vtype is None:
                        continue

                    protein_pos = parse_protein_position(tx.get("Protein_position", ""))
                    if protein_pos is None:
                        continue

                    wt_protein = tx.get("WildtypeProtein", "")
                    if not wt_protein or "*" in wt_protein:
                        continue

                    fs_protein = tx.get("FrameshiftSequence", "")
                    if vtype == "FS":
                        if not fs_protein:
                            continue
                        aa_change = f"{protein_pos}{ref}/{alt}"
                    else:
                        aa_field = tx.get("Amino_acids", "")
                        if not aa_field or "/" not in aa_field:
                            continue
                        aa_change = protein_pos + aa_field

                    gene = tx.get("SYMBOL") or tx.get("Gene") or "UNKNOWN"
                    transcript = tx.get("Feature") or "UNKNOWN"
                    vclass = VARIANT_CLASS_MAP.get(vtype, vtype)
                    index = construct_index(count, gene, transcript, vclass, aa_change)
                    if index in seen_indexes:
                        continue
                    seen_indexes.add(index)

                    yield {
                        "index": index,
                        "variant_type": vtype,
                        "gene": gene,
                        "transcript": transcript,
                        "protein_position": protein_pos,
                        "amino_acid_change": aa_change,
                        "wildtype_protein": wt_protein,
                        "frameshift_protein": fs_protein,
                        "amino_acids": tx.get("Amino_acids", ""),
                    }
                    count += 1
                    break  # one highest-priority transcript per ALT


def generate_fasta_records(
    variants: Iterator[dict],
    flanking: int,
    downstream_len: int | None,
    epitope_len: int,
) -> OrderedDict[tuple[str, str], str]:
    """Return {(kind, index): sequence}; kind is WT or MT."""
    records: OrderedDict[tuple[str, str], str] = OrderedDict()
    skipped = 0
    kept = 0

    for var in variants:
        result = build_mutant_sequence(
            var["variant_type"],
            var["wildtype_protein"],
            var["frameshift_protein"],
            var["protein_position"],
            var["amino_acids"],
            flanking,
            downstream_len,
            epitope_len,
        )
        if result is None:
            skipped += 1
            continue

        wt_sub, mt_sub = result
        idx = var["index"]
        records[("WT", idx)] = wt_sub
        records[("MT", idx)] = mt_sub
        kept += 1

    print(f"[INFO] Coding variants → peptides: kept {kept}, skipped {skipped}")
    return records


def write_fasta(records: OrderedDict[tuple[str, str], str], out_mt: str, out_wt: str, out_combined: str):
    with open(out_combined, "w", encoding="utf-8") as fc, \
         open(out_mt, "w", encoding="utf-8") as fm, \
         open(out_wt, "w", encoding="utf-8") as fw:
        for (kind, index), seq in records.items():
            header = f">{kind}.{index}"
            fc.write(f"{header}\n{seq}\n")
            if kind == "MT":
                fm.write(f"{header}\n{seq}\n")
            else:
                fw.write(f"{header}\n{seq}\n")


def parse_args():
    p = argparse.ArgumentParser(
        description="Coding-region VCF → pVACtools-style MT/WT peptide FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input_vcf", help="VEP-annotated coding-region VCF")
    p.add_argument("output_prefix", help="Output prefix; generates .fasta / _MT.fasta / _WT.fasta")
    p.add_argument("--flanking-length", type=int, default=24, help="Amino acids on each side of mutation (Methods: 24)")
    p.add_argument("--downstream-length", type=int, default=0, help="Frameshift MT downstream aa (0 = to protein end)")
    p.add_argument("--epitope-length", type=int, default=0, help="If >0, keep only peptides containing neo-epitope")
    p.add_argument("--no-canonical-only", action="store_true", help="Do not limit to canonical transcripts")
    p.add_argument("--include-filtered", action="store_true", help="Include variants with FILTER != PASS")
    return p.parse_args()


def main():
    args = parse_args()
    variants = iter_coding_variants(
        args.input_vcf,
        pass_only=not args.include_filtered,
        canonical_only=not args.no_canonical_only,
    )
    records = generate_fasta_records(
        variants,
        flanking=args.flanking_length,
        downstream_len=args.downstream_length if args.downstream_length > 0 else None,
        epitope_len=args.epitope_length,
    )

    prefix = args.output_prefix
    out_combined = f"{prefix}.fasta"
    out_mt = f"{prefix}_MT.fasta"
    out_wt = f"{prefix}_WT.fasta"
    write_fasta(records, out_mt, out_wt, out_combined)

    n_mt = sum(1 for k in records if k[0] == "MT")
    print(f"✅ Coding peptide FASTA written")
    print(f"   Combined: {out_combined}  ({len(records)} sequences)")
    print(f"   MT library: {out_mt}  ({n_mt} entries)")
    print(f"   WT library: {out_wt}  ({n_mt} entries)")


if __name__ == "__main__":
    main()
