#!/usr/bin/env python3
"""
generate_noncoding_peptides.py
------------------------------
Non-coding VCF → mutant neoantigen candidate peptide FASTA

Output FASTA header format:
  >NC|chr1|123456|A>T|GENE_NAME|ENSG00000xxxxxx|splice_region_variant|ENST00000xxxxxx|strand:+|mode:gtf|pep:9aa|idx:1

Field description:
  NC            fixed prefix, marks non-coding source
  chr1          chromosome
  123456        genomic position (1-based, same as VCF)
  A>T           ref>alt
  GENE_NAME     gene symbol (from GTF gene_name field)
  ENSGxxxxxxx   Ensembl gene ID
  consequence   VEP mutation consequence (from CSQ field)
  ENSTxxxxxxx   transcript ID
  strand        transcript strand
  mode          gtf (exon splicing) or fallback (genomic window)
  pep           peptide length
  idx           peptide index at same locus

Dependencies:
  pip install pyfaidx gtfparse pandas

Usage:
  python generate_noncoding_peptides.py \\
      noncoding.vcf \\
      hg38.fa \\
      gencode.v47.basic.annotation.gtf \\
      output.fasta \\
      [--pep_len 8,9,10,11] \\
      [--window 100]
"""

import sys
import argparse
from pathlib import Path
from collections import defaultdict

try:
    from pyfaidx import Fasta
except ImportError:
    sys.exit("❌ Install dependencies first: pip install pyfaidx")

try:
    import gtfparse
    import pandas as pd
    HAS_GTF = True
except ImportError:
    HAS_GTF = False
    print("[WARN] gtfparse not installed; all variants will use genomic window six-frame fallback")
    print("       Recommended: pip install gtfparse pandas")


# ---------------------------------------------------------------------------
# Genetic code table
# ---------------------------------------------------------------------------
CODON_TABLE = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGA':'*',
}


def reverse_complement(seq: str) -> str:
    comp = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    return ''.join(comp.get(b.upper(), 'N') for b in reversed(seq))


def translate_dna(dna: str) -> str:
    """Translate until stop codon or end of sequence; return protein without '*'."""
    protein = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def extract_peptides(protein: str, lengths=(8, 9, 10, 11)) -> set:
    """Sliding-window extract peptides of given lengths; filter low-quality fragments with X."""
    peps = set()
    for L in lengths:
        for i in range(len(protein) - L + 1):
            pep = protein[i:i+L]
            if 'X' not in pep:
                peps.add(pep)
    return peps


# ---------------------------------------------------------------------------
# VCF CSQ field parsing (extract gene, gene_id, consequence, transcript)
# ---------------------------------------------------------------------------

def parse_csq_header(header_lines: list) -> list:
    """Parse CSQ column order from VCF header."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            marker = "Format: "
            if marker in line:
                tail = line.split(marker, 1)[1].rstrip('">')
                return [x.strip() for x in tail.split("|")]
    return []


def extract_csq_annotations(info_str: str, csq_fields: list) -> list:
    """
    Parse CSQ annotations from INFO field; return list of all transcript annotations.
    Each element is a dict with SYMBOL, Gene, Feature, Consequence, etc.
    """
    if not csq_fields:
        return []
    annotations = []
    for item in info_str.split(";"):
        if item.startswith("CSQ="):
            csq_val = item[4:]
            for ann in csq_val.split(","):
                parts = ann.split("|")
                d = {}
                for i, field in enumerate(csq_fields):
                    d[field] = parts[i] if i < len(parts) else ""
                annotations.append(d)
            break
    return annotations


def best_csq_annotation(annotations: list) -> dict:
    """
    Pick the most representative transcript from multiple annotations:
    prefer CANONICAL=YES, then most severe consequence,
    finally fall back to first or empty dict.
    """
    if not annotations:
        return {}
    # Prefer canonical transcripts
    canonical = [a for a in annotations if a.get("CANONICAL", "") == "YES"]
    pool = canonical if canonical else annotations
    # Sort by consequence severity (string length proxy: splice > UTR > intron)
    severity_order = [
        "splice_acceptor", "splice_donor", "splice_region",
        "5_prime_utr", "3_prime_utr", "intron", "upstream", "downstream",
        "non_coding", "regulatory", "tf_binding",
    ]
    def severity_score(ann):
        csq = ann.get("Consequence", "").lower()
        for i, kw in enumerate(severity_order):
            if kw in csq:
                return i
        return len(severity_order)
    pool.sort(key=severity_score)
    return pool[0]


# ---------------------------------------------------------------------------
# GTF parsing and transcript structure cache
# ---------------------------------------------------------------------------

def load_transcripts(gtf_path: str):
    """
    Return transcripts[chrom] = [ {transcript_id, gene_id, gene_name,
                                   strand, exons, tx_start, tx_end}, ... ]
    exons: [(start, end), ...] 0-based half-open intervals, sorted by genomic coordinate.
    """
    print(f"[INFO] Reading GTF: {gtf_path}")
    df = gtfparse.read_gtf(gtf_path)
    df = df[df["feature"] == "exon"].copy()
    df["start"] = df["start"].astype(int) - 1   # GTF 1-based → 0-based
    df["end"]   = df["end"].astype(int)

    # Extract gene_name (GENCODE GTF has this field)
    has_gene_name = "gene_name" in df.columns

    transcripts = defaultdict(list)
    grouped = df.groupby(["seqname", "transcript_id"])
    for (chrom, tx_id), grp in grouped:
        strand    = grp["strand"].iloc[0]
        gene_id   = grp["gene_id"].iloc[0]   if "gene_id"   in grp.columns else ""
        gene_name = grp["gene_name"].iloc[0] if has_gene_name               else gene_id
        exons     = sorted(zip(grp["start"].tolist(), grp["end"].tolist()))
        tx_start  = exons[0][0]
        tx_end    = exons[-1][1]
        transcripts[chrom].append({
            "transcript_id": tx_id,
            "gene_id":       gene_id,
            "gene_name":     gene_name,
            "strand":        strand,
            "exons":         exons,
            "tx_start":      tx_start,
            "tx_end":        tx_end,
        })
    print(f"[INFO] Loaded {sum(len(v) for v in transcripts.values())} transcripts")
    return transcripts


def get_overlapping_transcripts(transcripts, chrom, pos):
    """Return transcripts overlapping pos (±500 bp upstream/downstream for splice region)."""
    MARGIN = 500
    return [
        tx for tx in transcripts.get(chrom, [])
        if tx["tx_start"] - MARGIN <= pos <= tx["tx_end"] + MARGIN
    ]


def build_mrna(tx: dict, ref: Fasta, chrom: str) -> tuple:
    """
    Splice exons into mRNA and return coordinate mapping table.
    coords: [(genome_start, genome_end, mrna_offset), ...]
    """
    chrom_len  = len(ref[chrom])
    mrna_parts = []
    coords     = []
    offset     = 0
    exons = tx["exons"] if tx["strand"] == "+" else list(reversed(tx["exons"]))

    for (ex_s, ex_e) in exons:
        ex_s = max(0, ex_s)
        ex_e = min(chrom_len, ex_e)
        seq  = str(ref[chrom][ex_s:ex_e])
        if tx["strand"] == "-":
            seq = reverse_complement(seq)
        coords.append((ex_s, ex_e, offset))
        mrna_parts.append(seq)
        offset += ex_e - ex_s

    return ''.join(mrna_parts), coords


def genome_to_mrna_offset(pos: int, coords: list, strand: str):
    """
    Genomic coordinate (0-based) → mRNA offset.
    Returns (offset, is_intronic).
    Intronic sites return nearest exon boundary offset with is_intronic=True.
    """
    if strand == "+":
        for (ex_s, ex_e, mrna_off) in coords:
            if ex_s <= pos < ex_e:
                return mrna_off + (pos - ex_s), False
        for (ex_s, ex_e, mrna_off) in sorted(coords, key=lambda x: x[0]):
            if pos < ex_s:
                return mrna_off, True
        last = sorted(coords, key=lambda x: x[0])[-1]
        return last[2] + (last[1] - last[0]) - 1, True
    else:
        for (ex_s, ex_e, mrna_off) in coords:
            if ex_s <= pos < ex_e:
                return mrna_off + (ex_e - 1 - pos), False
        for (ex_s, ex_e, mrna_off) in sorted(coords, key=lambda x: -x[1]):
            if pos >= ex_e:
                return mrna_off, True
        last = sorted(coords, key=lambda x: -x[1])[-1]
        return last[2] + (last[1] - last[0]) - 1, True


# ---------------------------------------------------------------------------
# Mutation introduction + differential peptide extraction
# ---------------------------------------------------------------------------

def apply_snv_to_mrna(mrna: str, mrna_offset: int,
                      ref_allele: str, alt_allele: str) -> str:
    return mrna[:mrna_offset] + alt_allele + mrna[mrna_offset + len(ref_allele):]


def get_diff_peptides(ref_mrna: str, mut_mrna: str,
                      mrna_offset: int, pep_lengths: tuple,
                      context_nt: int = 100) -> set:
    """
    Local window ±context_nt nt around mutation; return mutant-only peptides (set difference).
    Methods: 100 nt upstream and downstream of each mutation site.
    """
    win_s      = max(0, mrna_offset - context_nt)
    win_e_ref  = min(len(ref_mrna), mrna_offset + context_nt + 1)
    delta      = len(mut_mrna) - len(ref_mrna)
    win_e_mut  = min(len(mut_mrna), win_e_ref + max(0, delta))

    ref_window = ref_mrna[win_s:win_e_ref]
    mut_window = mut_mrna[win_s:win_e_mut]

    ref_peps = set()
    mut_peps = set()
    for frame in range(3):
        ref_peps |= extract_peptides(translate_dna(ref_window[frame:]), pep_lengths)
        mut_peps |= extract_peptides(translate_dna(mut_window[frame:]), pep_lengths)

    return mut_peps - ref_peps


# ---------------------------------------------------------------------------
# Fallback: genomic window six-frame translation
# ---------------------------------------------------------------------------

def fallback_six_frame(ref: Fasta, chrom: str, pos: int,
                       ref_allele: str, alt_allele: str,
                       pep_lengths: tuple, window: int) -> set:
    chrom_len = len(ref[chrom])
    start     = max(0, pos - window)
    end       = min(chrom_len, pos + window + len(ref_allele))
    ref_seq   = str(ref[chrom][start:end])
    off       = pos - start
    mut_seq   = ref_seq[:off] + alt_allele + ref_seq[off + len(ref_allele):]

    ref_peps = set()
    mut_peps = set()
    for frame in range(3):
        ref_peps |= extract_peptides(translate_dna(ref_seq[frame:]), pep_lengths)
        mut_peps |= extract_peptides(translate_dna(mut_seq[frame:]), pep_lengths)
    rev_ref = reverse_complement(ref_seq)
    rev_mut = reverse_complement(mut_seq)
    for frame in range(3):
        ref_peps |= extract_peptides(translate_dna(rev_ref[frame:]), pep_lengths)
        mut_peps |= extract_peptides(translate_dna(rev_mut[frame:]), pep_lengths)

    return mut_peps - ref_peps


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Non-coding VCF → mutant candidate peptide FASTA (transcript structure + full header)"
    )
    parser.add_argument("vcf",    help="Non-coding VCF file (VEP-annotated)")
    parser.add_argument("fasta",  help="Reference genome FASTA (requires .fai index)")
    parser.add_argument("gtf",    help="Genome annotation GTF (e.g. gencode.v47.basic.annotation.gtf)")
    parser.add_argument("output", help="Output FASTA file")
    parser.add_argument("--pep_len", default="8,9,10,11",
                        help="Peptide lengths, comma-separated (default 8,9,10,11)")
    parser.add_argument("--window", type=int, default=100,
                        help="Genomic window (bp) each side of mutation when no transcript overlap (Methods: 100)")
    args = parser.parse_args()

    pep_lengths = tuple(int(x) for x in args.pep_len.split(","))
    ref         = Fasta(args.fasta)

    transcripts = {}
    if HAS_GTF and Path(args.gtf).exists():
        transcripts = load_transcripts(args.gtf)
    else:
        print("[WARN] GTF unavailable; using genomic window fallback for all variants")

    total_variants = 0
    total_peptides = 0
    fallback_count = 0

    # Read VCF header (for CSQ column names)
    header_lines = []
    with open(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line.rstrip())
            else:
                break
    csq_fields = parse_csq_header(header_lines)
    if csq_fields:
        print(f"[INFO] CSQ field detected with {len(csq_fields)} annotation columns")
    else:
        print("[WARN] CSQ header not detected; gene/consequence will be labeled unknown")

    with open(args.vcf) as f, open(args.output, 'w') as out:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            chrom      = fields[0]
            pos_1based = int(fields[1])
            ref_allele = fields[3]
            alt_field  = fields[4]
            info_str   = fields[7] if len(fields) > 7 else "."
            pos        = pos_1based - 1   # 0-based

            if chrom not in ref:
                print(f"[WARN] Chromosome {chrom} not in reference genome; skipping")
                continue

            # Parse CSQ: extract gene / consequence / transcript
            annotations  = extract_csq_annotations(info_str, csq_fields)
            best_ann     = best_csq_annotation(annotations)
            gene_name    = best_ann.get("SYMBOL",      "unknown") or "unknown"
            gene_id      = best_ann.get("Gene",        "unknown") or "unknown"
            consequence  = best_ann.get("Consequence", "unknown") or "unknown"
            transcript   = best_ann.get("Feature",     "unknown") or "unknown"
            # Exon/intron number, e.g. "3/12" (3 of 12); empty if absent
            exon_no      = best_ann.get("EXON",   "") or ""
            intron_no    = best_ann.get("INTRON", "") or ""
            # Location tag: prefer exon number, then intron, else "."
            if exon_no:
                location_tag = f"exon{exon_no}"
            elif intron_no:
                location_tag = f"intron{intron_no}"
            else:
                location_tag = "."

            # Process each alt allele
            alt_alleles = [a.strip() for a in alt_field.split(',')
                           if a.strip() not in ('.', '*')]

            for alt_allele in alt_alleles:
                total_variants += 1
                diff_peps  = set()
                mode       = "gtf"
                tx_used    = "NA"
                strand_used = "."

                # ── GTF mode: exon splicing ─────────────────────────────────
                if transcripts:
                    overlapping = get_overlapping_transcripts(transcripts, chrom, pos)
                    for tx in overlapping:
                        try:
                            ref_mrna, coords = build_mrna(tx, ref, chrom)
                            mrna_offset, _   = genome_to_mrna_offset(
                                pos, coords, tx["strand"]
                            )
                            mut_mrna = apply_snv_to_mrna(
                                ref_mrna, mrna_offset, ref_allele, alt_allele
                            )
                            peps = get_diff_peptides(
                                ref_mrna, mut_mrna, mrna_offset, pep_lengths
                            )
                            if peps:
                                diff_peps   |= peps
                                tx_used      = tx["transcript_id"]
                                strand_used  = tx["strand"]
                                # Fill gene info from GTF if CSQ missing
                                if gene_name == "unknown":
                                    gene_name = tx.get("gene_name", "unknown")
                                if gene_id == "unknown":
                                    gene_id   = tx.get("gene_id",   "unknown")
                        except Exception as e:
                            print(f"[WARN] Failed to process transcript {tx['transcript_id']}: {e}")
                            continue

                # ── Fallback: genomic window six-frame translation ──────────
                if not diff_peps:
                    diff_peps = fallback_six_frame(
                        ref, chrom, pos, ref_allele, alt_allele,
                        pep_lengths, args.window
                    )
                    mode        = "fallback"
                    tx_used     = "NA"
                    strand_used = "."
                    fallback_count += 1

                # ── Write FASTA ─────────────────────────────────────────────
                # header format:
                # >NC|chrom|pos|ref>alt|gene_name|gene_id|consequence|location|transcript|strand|mode|pepLen|idx
                for idx, pep in enumerate(sorted(diff_peps), 1):
                    header = (
                        f">NC"
                        f"|{chrom}"
                        f"|{pos_1based}"
                        f"|{ref_allele}>{alt_allele}"
                        f"|{gene_name}"
                        f"|{gene_id}"
                        f"|{consequence}"
                        f"|{location_tag}"
                        f"|{tx_used}"
                        f"|strand:{strand_used}"
                        f"|mode:{mode}"
                        f"|pep:{len(pep)}aa"
                        f"|idx:{idx}"
                    )
                    out.write(f"{header}\n{pep}\n")
                    total_peptides += 1

    print(f"\n✅ Processing complete")
    print(f"   Variant sites (allele count): {total_variants}")
    print(f"   Total candidate peptides    : {total_peptides}")
    print(f"   Fallback to genomic window  : {fallback_count} sites")
    print(f"   Output file                 : {args.output}")


if __name__ == "__main__":
    main()
