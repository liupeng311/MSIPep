import csv
import argparse
from concurrent.futures import ThreadPoolExecutor

def comet_to_fasta(input_file, output_fasta, output_combined_fasta, xcorr_threshold=None, evalue_threshold=None):
    """
    Process Comet result file to extract peptides into a main FASTA file and an optional
    filtered FASTA file for Class I peptides (8-11 length).

    Parameters:
    - input_file: Path to the Comet result file (tab-separated format).
    - output_fasta: Path to save all extracted peptides in FASTA format.
    - output_combined_fasta: Path to save filtered peptides (8-11 length) in FASTA format.
    - xcorr_threshold: Optional threshold for xcorr filtering.
    - evalue_threshold: Optional threshold for e-value filtering.
    """
    try:
        with open(input_file, 'r') as infile, open(output_fasta, 'w') as main_fasta, \
                open(output_combined_fasta, 'w') as combined_fasta:

            # Skip the first line
            infile.readline()

            # Use DictReader to read data rows after the headers
            reader = csv.DictReader(infile, delimiter='\t', skipinitialspace=True)
            sequence_count = 1

            for row in reader:
                try:
                    peptide_sequence = row['plain_peptide']
                    protein_name = row['protein']
                    peptide_length = len(peptide_sequence)

                    # Apply filtering if thresholds are specified
                    if xcorr_threshold is not None and evalue_threshold is not None:
                        xcorr = float(row['xcorr'])
                        evalue = float(row['e-value'])
                        if xcorr < xcorr_threshold or evalue > evalue_threshold:
                            continue

                    # Write to the main FASTA file
                    main_fasta.write(f">{protein_name}\n{peptide_sequence}\n")

                    # Write peptides with length 8-11 to the combined FASTA
                    if 8 <= peptide_length <= 11:
                        combined_fasta.write(f">Sequence_{sequence_count}\n{peptide_sequence}\n")
                        sequence_count += 1

                except KeyError as e:
                    print(f"Missing expected column: {e} in row {row}")
                except ValueError as e:
                    print(f"Error converting xcorr or e-value: {e} in row {row}")

        print(f"FASTA files for {input_file} created successfully.")
    except Exception as e:
        print(f"Error processing {input_file}: {e}")

def process_files_in_parallel(file_paths, output_fasta_paths, output_combined_fasta_paths, xcorr_threshold,
                              evalue_threshold, threads):
    """
    Process multiple Comet result files in parallel.
    """
    if len(file_paths) != len(output_fasta_paths) or len(file_paths) != len(output_combined_fasta_paths):
        print("Error: The number of input files must match the number of output FASTA paths.")
        return

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for input_file, output_fasta, output_combined_fasta in zip(file_paths, output_fasta_paths,
                                                                   output_combined_fasta_paths):
            futures.append(executor.submit(
                comet_to_fasta, input_file, output_fasta, output_combined_fasta,
                xcorr_threshold, evalue_threshold
            ))

        for future in futures:
            future.result()  # Wait for all threads to complete

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Comet search result files to extract peptides.")
    parser.add_argument('--input_files', type=str, nargs='+', required=True,
                        help="Paths to the Comet result input files.")
    parser.add_argument('--output_fasta', type=str, nargs='+', required=True,
                        help="Paths to the output FASTA files for all peptides.")
    parser.add_argument('--output_combined_fasta', type=str, nargs='+', required=True,
                        help="Paths to the output FASTA files for Class I peptides (8-11 length).")
    parser.add_argument('--xcorr_threshold', type=float, default=None, help="Optional xcorr threshold for filtering.")
    parser.add_argument('--evalue_threshold', type=float, default=None,
                        help="Optional e-value threshold for filtering.")
    parser.add_argument('--threads', type=int, default=1, help="Number of threads to use for parallel processing.")

    args = parser.parse_args()

    process_files_in_parallel(
        file_paths=args.input_files,
        output_fasta_paths=args.output_fasta,
        output_combined_fasta_paths=args.output_combined_fasta,
        xcorr_threshold=args.xcorr_threshold,
        evalue_threshold=args.evalue_threshold,
        threads=args.threads
    )
