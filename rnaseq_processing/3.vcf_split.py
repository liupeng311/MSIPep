#!/usr/bin/env python3
"""
Split a VEP-annotated VCF into coding and non-coding region files
"""

import sys
import gzip
from pathlib import Path


def is_coding_consequence(consequence_str: str) -> bool:
    """Determine whether the consequence is coding-region related"""
    coding_terms = {
        'missense_variant', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion',
        'stop_gained', 'stop_lost', 'start_lost', 'synonymous_variant',
        'splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant',
        'protein_altering_variant', 'coding_sequence_variant', 'incomplete_terminal_codon_variant'
    }

    consequences = consequence_str.upper().split('&')  # handle multiple consequences separated by &
    for cons in consequences:
        if cons in coding_terms or any(term in cons for term in ['MISSENSE', 'FRAMESHIFT', 'SPLICE', 'STOP']):
            return True
    return False


def split_vep_vcf(input_vcf: str, output_coding: str, output_noncoding: str):
    open_func = gzip.open if str(input_vcf).endswith('.gz') else open
    mode = 'rt'

    coding_count = noncoding_count = 0

    with open_func(input_vcf, mode) as fin, \
            open(output_coding, 'w') as f_coding, \
            open(output_noncoding, 'w') as f_noncoding:

        for line in fin:
            if line.startswith('#'):  # keep header
                f_coding.write(line)
                f_noncoding.write(line)
                continue

            # Find CSQ field
            if 'CSQ=' not in line:
                # unannotated variants go to non-coding by default
                f_noncoding.write(line)
                noncoding_count += 1
                continue

            # Extract CSQ info
            csq_part = None
            for field in line.strip().split('\t')[7].split(';'):  # INFO column
                if field.startswith('CSQ='):
                    csq_part = field[4:]  # strip CSQ=
                    break

            if not csq_part:
                f_noncoding.write(line)
                noncoding_count += 1
                continue

            # Use Consequence from first transcript (usually most important)
            first_csq = csq_part.split(',')[0]
            consequence = first_csq.split('|')[1]  # Consequence is column 2 (index 1)

            if is_coding_consequence(consequence):
                f_coding.write(line)
                coding_count += 1
            else:
                f_noncoding.write(line)
                noncoding_count += 1

    print(f"✅ Split complete!")
    print(f"   Coding variants    : {coding_count} → {output_coding}")
    print(f"   Non-coding variants: {noncoding_count} → {output_noncoding}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python split_vep_vcf.py <input.vcf> <output_coding.vcf> <output_noncoding.vcf>")
        sys.exit(1)

    split_vep_vcf(sys.argv[1], sys.argv[2], sys.argv[3])
