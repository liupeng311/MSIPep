import os
import glob
import pandas as pd
import argparse
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract (protein, peptide) pairs from multiple TSV files in a folder, merge and output into a single FASTA file, keeping only peptides of length 9–11, removing duplicates and merging protein info for identical peptides."
    )
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to input directory containing .tsv files")
    parser.add_argument("-o", "--output_fasta", required=True,
                        help="Path to output FASTA file")
    return parser.parse_args()

def main():
    args = parse_arguments()
    input_dir = args.input_dir
    output_fasta = args.output_fasta
    # Search for all .tsv files in the folder
    tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
    if not tsv_files:
        print(f"No .tsv files found in directory {input_dir}, exiting.")
        return
    print(f"Found {len(tsv_files)} TSV files in directory {input_dir}.")
    print(f"Results will be written to {output_fasta}")
    # Dictionary to store peptide as key and corresponding set of proteins as value
    peptide_protein_map = defaultdict(set)
    # Process each file
    for tsv_file in tsv_files:
        print(f"Processing file: {tsv_file}")
        try:
            # Read as string type to avoid confusion
            df = pd.read_csv(tsv_file, sep="\t", dtype=str)
        except Exception as e:
            print(f"Error reading file {tsv_file}: {e}")
            continue
        # Check if required columns exist
        if "protein" not in df.columns or "peptide" not in df.columns:
            print(f"File {tsv_file} does not contain 'protein' or 'peptide' columns, skipping.")
            continue
        # Iterate through each record
        for _, row in df.iterrows():
            protein = row["protein"]
            peptide = row["peptide"]
            # Skip empty or missing values
            if pd.isna(protein) or pd.isna(peptide):
                continue
            # Keep only peptides of length 9 to 11
            if 9 <= len(peptide) <= 11:
                # Add the protein to the set corresponding to this peptide
                peptide_protein_map[peptide].add(protein)

    # Write to FASTA file
    with open(output_fasta, "w") as fasta_out:
        for peptide, protein_set in peptide_protein_map.items():
            # protein_set contains all proteins matching this peptide
            # Merge into one header line, using '|' as delimiter
            merged_proteins = "|".join(sorted(protein_set))
            fasta_header = f">{merged_proteins}"
            fasta_out.write(f"{fasta_header}\n{peptide}\n")

    print(f"Processing complete. FASTA file written to: {output_fasta}")

if __name__ == "__main__":
    main()

# Example usage:
# python extract_classI_peptides.py -i /path/to/tsv_folder -o /path/to/output.fasta
