#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import logging
import sys

# Set logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define codon translation table
codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

# Define manual translation function: translate every 3 nucleotides into an amino acid (ignoring stop codons)
def translate_sequence(seq, codon_table):
    protein = []
    if len(seq) < 3:
        return ""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        amino_acid = codon_table.get(codon, "X")
        if amino_acid != '*':
            protein.append(amino_acid)
    return ''.join(protein)

# Reverse complement sequence
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

# Six-frame translation: 3 forward frames + 3 reverse-complement frames
def six_frame_translation(seq, codon_table):
    # Retain only A/T/C/G
    seq = ''.join([base for base in seq if base in 'ATCG'])

    # Forward 3 frames
    frame1 = translate_sequence(seq, codon_table)
    frame2 = translate_sequence(seq[1:], codon_table)
    frame3 = translate_sequence(seq[2:], codon_table)

    # Reverse complement frames
    rev_seq = reverse_complement(seq)
    rev_frame1 = translate_sequence(rev_seq, codon_table)
    rev_frame2 = translate_sequence(rev_seq[1:], codon_table)
    rev_frame3 = translate_sequence(rev_seq[2:], codon_table)

    return [frame1, frame2, frame3, rev_frame1, rev_frame2, rev_frame3]

# Extract sequence from reference genome around pos±flank
def extract_sequence_from_genome(genome, chrom, start, flank=100):
    if chrom not in genome:
        raise ValueError(f"Chromosome {chrom} not found in the reference genome.")
    seq_len = len(genome[chrom])
    # Note: Depending on your coordinate system (0-based), you may need to use pos-1
    start_idx = max(0, (start - 1) - flank)
    end_idx = min(seq_len, (start - 1) + flank)
    ref_seq = genome[chrom].seq[start_idx:end_idx]
    return str(ref_seq)

# Generate FASTA for a single mutation: only translate ALT sequence (mutated)
def process_mutation(row, genome, codon_table, flank):
    """
    Only output the 6-frame translation of the ALT (mutated) sequence
    """
    chrom = row['Chr']
    pos   = int(row['Start'])  # Ensure this is 1-based or 0-based as per input
    ref   = row['Ref']
    alt   = row['Alt']
    gene  = row.get('Gene.refGene', 'NA')  # Default to 'NA' if column is missing
    func = row.get('Func.refGene', 'NA')
    # 1) Get reference sequence (±100bp)
    try:
        ref_seq = extract_sequence_from_genome(genome, chrom, pos, flank)
    except ValueError as e:
        logging.error(e)
        return None

    # 2) Replace ref with alt in the sequence
    #    ref_seq[:flank] + alt + ref_seq[flank + len(ref):]
    mutated_seq = ref_seq[:flank] + alt + ref_seq[flank + len(ref):]
    mutated_seq = mutated_seq.upper()

    # 3) Translate only mutated sequence in six frames
    alt_translation = six_frame_translation(mutated_seq, codon_table)

    # 4) Write 6-frame translations in FASTA format: Frame1_ALT ~ Frame6_ALT
    fasta_entries = []
    for i, frame_seq in enumerate(alt_translation):
        frame_type = 'pl' if i < 3 else 'nl'  # pl = plus, nl = minus
        header = f">{chrom}|{pos}|{ref}|{alt}|{frame_type}|{gene}|{func}|Frame{i + 1}_ALT"
        fasta_entries.append(f"{header}\n{frame_seq}\n")

    # Return 6 translated sequences
    return ''.join(fasta_entries)

# Multithreaded FASTA generation: only ALT sequence translations
def generate_fasta(df, output_fasta, genome_file, codon_table, flank=100, threads=32):
    try:
        genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
        logging.info("Reference genome loaded successfully.")
    except Exception as e:
        logging.error(f"Failed to load reference genome: {e}")
        sys.exit(1)

    fasta_entries = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for index, row in df.iterrows():
            futures.append(
                executor.submit(process_mutation, row, genome, codon_table, flank)
            )

        for future in as_completed(futures):
            result = future.result()
            if result:
                fasta_entries.append(result)

    # Write all sequences to FASTA
    with open(output_fasta, 'w') as fasta:
        fasta.write(''.join(fasta_entries))
    logging.info(f"FASTA file generated: {output_fasta}")

# Read txt input file
def read_txt(file):
    try:
        df = pd.read_csv(file, sep="\t")
        logging.info(f"Input file {file} loaded successfully.")
        return df
    except Exception as e:
        logging.error(f"Failed to read input file: {e}")
        sys.exit(1)

# Main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a FASTA file from mutations (only ALT) with six-frame translation.')
    parser.add_argument('-i', '--input', required=True, help='Input txt file with mutation data')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file path')
#    parser.add_argument('-g', '--genome', required=True, help='Reference genome file path (FASTA format)')
    parser.add_argument('-t', '--threads', type=int, default=16, help='Number of threads (default: 16)')
    parser.add_argument('-f', '--flank', type=int, default=100, help='Flank size around mutation site (default: 100)')

    args = parser.parse_args()
    reference_genome_path = "./reference/hg38.fa"

    mutation_df = read_txt(args.input)
    generate_fasta(
        mutation_df,
        args.output,
        reference_genome_path,
        codon_table,
        flank=args.flank,
        threads=args.threads
    )

# Example usage:
# python script_alt_only.py -i input_file.txt -o output_file.fasta
