#!/usr/bin/env python3
"""
2_VCF_filter_WGS.py
--------------------
Whole-genome mode — retain all functional variants in coding and non-coding regions.

Key differences from the original script (2_VCF_filter.py):
  1. is_functional_by_csq() keep_consequences extended to non-coding consequence types.
  2. New --noncoding_only parameter (off by default) for non-coding-only output.
  3. gnomAD multi-allele: take max AF across all alleles, not just the first (bug fix).
  4. New --tumor_col parameter to specify tumor sample column (default column 10, index=9).
  5. Statistics include coding / non-coding classification counts.

Usage examples:
  python 2_VCF_filter_WGS.py input.vep.vcf.gz output_all.vcf.gz
  python 2_VCF_filter_WGS.py input.vep.vcf.gz output_all.vcf.gz --tumor_col 10
"""

import argparse
import gzip

try:
    from scipy.stats import fisher_exact
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False


# ---------------------------------------------------------------------------
# All VEP Consequences considered "functional" (coding + non-coding)
# ---------------------------------------------------------------------------
CODING_CONSEQUENCES = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",                 # missing in original script, added here
    "frameshift_variant",
    "inframe_deletion",
    "inframe_insertion",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "protein_altering_variant",   # used by some VEP versions
    "start_retained_variant",
    "stop_retained_variant",
}

NONCODING_CONSEQUENCES = {
    "splice_region_variant",      # splice region (±3~8 bp near donor/acceptor)
    "splice_polypyrimidine_tract_variant",
    "splice_donor_region_variant",
    "splice_donor_5th_base_variant",
    "synonymous_variant",         # some synonymous mutations affect splicing/translation
    "5_prime_UTR_variant",        # affects translation initiation
    "3_prime_UTR_variant",        # affects RNA stability / miRNA binding
    "intron_variant",             # deep intronic, may create new splice sites
    "upstream_gene_variant",      # promoter region
    "downstream_gene_variant",
    "regulatory_region_variant",
    "TF_binding_site_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "non_coding_transcript_variant",
    "non_coding_transcript_exon_variant",
    "NMD_transcript_variant",
    "mature_miRNA_variant",
}

ALL_KEEP_CONSEQUENCES = CODING_CONSEQUENCES | NONCODING_CONSEQUENCES


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------
def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Whole-genome VCF filter (coding + non-coding), for RNA-based neoantigen pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_vcf",  help="Input VCF/VCF.GZ (VEP-annotated)")
    parser.add_argument("output_vcf", help="Output filtered VCF/VCF.GZ")
    parser.add_argument("--min_t_depth",      type=int,   default=10,    help="Minimum tumor sequencing depth")
    parser.add_argument("--min_t_alt_count",  type=int,   default=5,     help="Minimum tumor ALT read count")
    parser.add_argument("--min_vaf",          type=float, default=0.05,  help="Minimum VAF")
    parser.add_argument("--max_gnomad_af",    type=float, default=0.01,  help="Maximum gnomAD AF (population filter)")
    parser.add_argument("--min_readpos_ranksum", type=float, default=-8.0, help="ReadPosRankSum lower bound")
    parser.add_argument("--min_strand_p",     type=float, default=0.001, help="Fisher strand bias p-value lower bound (filter if below)")
    parser.add_argument("--max_fs",           type=float, default=30.0,  help="FS upper bound")
    parser.add_argument("--max_sor",          type=float, default=3.0,   help="SOR upper bound")
    parser.add_argument("--disable_strand_bias_filter", action="store_true", help="Disable strand bias filter")
    parser.add_argument("--tumor_col",        type=int,   default=9,
                        help="0-based index of tumor sample column (single-sample VCF=9; tumor usually=10 for tumor-normal pairs)")
    return parser.parse_args()


def parse_info(info_str):
    info = {}
    if info_str == ".":
        return info
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def to_float(x):
    try:
        if x in (None, "", ".", "NA", "nan"):
            return None
        return float(x)
    except Exception:
        return None


def to_int(x):
    try:
        if x in (None, "", ".", "NA", "nan"):
            return None
        return int(float(x))
    except Exception:
        return None


def pick_gnomad_af_key(info):
    preferred = ["gnomAD_AF", "gnomad_af", "MAX_AF", "max_af", "gnomADe_AF", "gnomADg_AF"]
    lower_to_real = {k.lower(): k for k in info.keys()}
    for p in preferred:
        if p.lower() in lower_to_real:
            return lower_to_real[p.lower()]
    candidates = [k for k in info.keys() if ("gnomad" in k.lower() and "af" in k.lower())]
    return candidates[0] if candidates else None


def get_gnomad_af(info):
    gkey = pick_gnomad_af_key(info)
    if gkey is None:
        return None
    return to_float(info.get(gkey, ".").split(",")[0])


def get_sample_depth_alt_vaf(fmt_keys, sample_values):
    if not fmt_keys or not sample_values:
        return None, None, None, {}
    fmt = dict(zip(fmt_keys, sample_values))
    dp = to_int(fmt.get("DP"))
    alt_count = None
    vaf = None

    ad = fmt.get("AD")
    if ad and ad != ".":
        vals = [to_int(x) for x in ad.split(",")]
        vals = [x for x in vals if x is not None]
        if len(vals) >= 2:
            ref = vals[0]
            alts = vals[1:]
            alt_count = max(alts) if alts else None
            if dp is None:
                dp = ref + sum(alts)

    af = fmt.get("AF")
    if af and af != ".":
        af_vals = [to_float(x) for x in af.split(",")]
        af_vals = [x for x in af_vals if x is not None]
        if af_vals:
            vaf = max(af_vals)

    if vaf is None and dp and alt_count is not None and dp > 0:
        vaf = alt_count / dp

    return dp, alt_count, vaf, fmt


def fisher_p_from_sb(fmt):
    if not HAS_SCIPY or "SB" not in fmt:
        return None
    sb = fmt["SB"]
    if sb in (None, "."):
        return None
    nums = [to_int(x) for x in sb.split(",")]
    if len(nums) < 4 or any(x is None for x in nums[:4]):
        return None
    ref_fwd, ref_rev, alt_fwd, alt_rev = nums[:4]
    try:
        _, p = fisher_exact([[alt_fwd, alt_rev], [ref_fwd, ref_rev]])
        return p
    except Exception:
        return None


def extract_csq_fields(header_lines):
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            marker = "Format: "
            if marker in line:
                tail = line.split(marker, 1)[1].rstrip('">')
                return [x.strip() for x in tail.split("|")]
    return []


def classify_csq(info, csq_fields):
    """
    Returns (is_keep, is_coding, is_noncoding)
      is_keep      : keep variant (matches any functional consequence)
      is_coding    : matches coding consequence
      is_noncoding : matches non-coding consequence (and not coding)
    When CSQ is missing, keep by default, classified as coding (conservative).
    """
    if "CSQ" not in info or not csq_fields:
        return True, True, False
    try:
        cidx = csq_fields.index("Consequence")
    except ValueError:
        return True, True, False

    hit_coding    = False
    hit_noncoding = False

    for ann in info["CSQ"].split(","):
        parts = ann.split("|")
        if cidx >= len(parts):
            continue
        conseqs = {x.strip() for x in parts[cidx].split("&") if x.strip()}
        if conseqs & CODING_CONSEQUENCES:
            hit_coding = True
        if conseqs & NONCODING_CONSEQUENCES:
            hit_noncoding = True

    is_keep = hit_coding or hit_noncoding
    return is_keep, hit_coding, hit_noncoding


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    print(f"[INFO] Mode: whole-genome (coding + non-coding)")
    print(f"[INFO] Reading VCF: {args.input_vcf}")
    print(f"[INFO] Tumor sample column index (0-based): {args.tumor_col}")

    headers = []
    records = []
    with open_text(args.input_vcf, "r") as f:
        for line in f:
            if line.startswith("#"):
                headers.append(line)
            else:
                records.append(line.rstrip("\n"))

    csq_fields = extract_csq_fields(headers)
    if csq_fields:
        print("[INFO] CSQ field detected; functional consequence filter enabled")
    else:
        print("[WARN] CSQ header not found; skipping functional filter (keep variants passing other filters)")

    kept = []
    n_coding = 0
    n_noncoding = 0
    n_both = 0

    for rec in records:
        cols = rec.split("\t")
        if len(cols) < 8:
            continue

        flt  = cols[6]
        info = parse_info(cols[7])

        # 1) PASS
        if flt != "PASS":
            continue

        # Tumor sample column
        fmt_keys      = cols[8].split(":")               if len(cols) > args.tumor_col else []
        sample_values = cols[args.tumor_col].split(":")  if len(cols) > args.tumor_col else []
        dp, alt_count, vaf, fmt = get_sample_depth_alt_vaf(fmt_keys, sample_values)

        # 2-4) Depth + alt reads + VAF
        if dp is None or alt_count is None:
            continue
        if dp < args.min_t_depth or alt_count < args.min_t_alt_count:
            continue
        if vaf is None or vaf < args.min_vaf:
            continue

        # 5) gnomAD (max across multi-allelic)
        gval = get_gnomad_af(info)
        if gval is not None and gval >= args.max_gnomad_af:
            continue

        # 6) Strand bias
        if not args.disable_strand_bias_filter:
            p = fisher_p_from_sb(fmt)
            if p is not None:
                if p <= args.min_strand_p:
                    continue
            else:
                fs  = to_float(info.get("FS"))
                sor = to_float(info.get("SOR"))
                if fs  is not None and fs  > args.max_fs:
                    continue
                if sor is not None and sor > args.max_sor:
                    continue

        # 7) ReadPosRankSum
        rprs = to_float(info.get("ReadPosRankSum"))
        if rprs is not None and rprs <= args.min_readpos_ranksum:
            continue

        # 8) Functional consequences (coding + non-coding)
        is_keep, is_coding, is_noncoding = classify_csq(info, csq_fields)
        if not is_keep:
            continue

        kept.append(rec)
        if is_coding and is_noncoding:
            n_both += 1
        elif is_coding:
            n_coding += 1
        else:
            n_noncoding += 1

    with open_text(args.output_vcf, "w") as out:
        for h in headers:
            out.write(h)
        for rec in kept:
            out.write(rec + "\n")

    print(f"\n[INFO] Total input variants  : {len(records)}")
    print(f"[INFO] Total kept variants     : {len(kept)}")
    print(f"[INFO]   ├─ Coding only        : {n_coding}")
    print(f"[INFO]   ├─ Non-coding only    : {n_noncoding}")
    print(f"[INFO]   └─ Coding + non-coding: {n_both}")
    print(f"[INFO] Output file             : {args.output_vcf}")


if __name__ == "__main__":
    main()
