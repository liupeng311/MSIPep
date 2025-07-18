import argparse
import os
import subprocess
import logging
import csv

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    logging.info(f"Running command: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        logging.info("Command executed successfully")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command execution failed: {e}")
        raise

def run_rna_mutation_pipeline(rna_seq_1P, rna_seq_2P, threads):
    run_command(f"python rna_mutation_pipeline.py -rna_seq_1P {rna_seq_1P} -rna_seq_2P {rna_seq_2P} -threads {threads}")

def filter_annovar(input_file, col1_index, col2_index, threshold):
    base_dir = os.path.dirname(input_file)
    output_coding = os.path.join(base_dir, "coding.txt")
    output_non_coding = os.path.join(base_dir, "noncoding-pep.txt")
    coding_framstop = os.path.join(base_dir, "coding_framstop.txt")
    coding_nonfrasnv = os.path.join(base_dir, "coding_snvnonfra.txt")

    with open(input_file, 'r') as infile, \
         open(output_coding, 'w', newline='') as coding_file, \
         open(output_non_coding, 'w', newline='') as non_coding_file:

        reader = csv.reader(infile, delimiter='\t')
        coding_writer = csv.writer(coding_file, delimiter='\t')
        non_coding_writer = csv.writer(non_coding_file, delimiter='\t')

        header = next(reader)
        coding_writer.writerow(header)
        non_coding_writer.writerow(header)

        otherinfo10_idx = header.index('Otherinfo10')

        for row in reader:
            try:
                value1 = float(row[col1_index]) if row[col1_index] != "." else 0
                value2 = float(row[col2_index]) if row[col2_index] != "." else 0
            except ValueError:
                continue

            if row[otherinfo10_idx] != 'PASS':
                continue

            if value1 <= threshold and value2 <= threshold:
                if row[5] == 'exonic':
                    if row[8] != 'synonymous SNV':
                        coding_writer.writerow(row)
                else:
                    non_coding_writer.writerow(row)

    with open(output_coding, 'r') as coding_file, \
         open(coding_framstop, 'w', newline='') as framstop_file, \
         open(coding_nonfrasnv, 'w', newline='') as nonfrasnv_file:

        reader = csv.reader(coding_file, delimiter='\t')
        framstop_writer = csv.writer(framstop_file, delimiter='\t')
        nonfrasnv_writer = csv.writer(nonfrasnv_file, delimiter='\t')

        header = next(reader)
        framstop_writer.writerow(header)
        nonfrasnv_writer.writerow(header)

        for row in reader:
            if row[8] in ['stoploss', 'frameshift deletion', 'frameshift insertion']:
                framstop_writer.writerow(row)
            else:
                nonfrasnv_writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Integrated RNA mutation analysis pipeline")
    parser.add_argument("-1", "--rna_seq_1P", required=True, help="Path to first RNA-Seq read file")
    parser.add_argument("-2", "--rna_seq_2P", required=True, help="Path to second RNA-Seq read file")
    parser.add_argument("-t", "--threads", type=int, default=16, help="Number of threads to use")
    parser.add_argument("--filter_col1", type=int, required=True, help="Column 1 for filtering (1-based)")
    parser.add_argument("--filter_col2", type=int, required=True, help="Column 2 for filtering (1-based)")
    parser.add_argument("--threshold", type=float, default=0.05, help="Filtering threshold")
    args = parser.parse_args()

    # Step 1: Run rna_mutation_pipeline.py
    run_rna_mutation_pipeline(args.rna_seq_1P, args.rna_seq_2P, args.threads)

    input_dir = os.path.dirname(args.rna_seq_1P)
    annovar_txt = os.path.join(input_dir, "mutadddb.hg38_multianno.txt")

    # Step 2: Filter ANNOVAR output
    col1_index = args.filter_col1 + 10 - 1
    col2_index = args.filter_col2 + 10 - 1
    filter_annovar(annovar_txt, col1_index, col2_index, args.threshold)

    # Step 3: Run downstream peptide generation scripts
    run_command(f"python coding-pep.py -i {os.path.join(input_dir, 'noncoding-pep.txt')} -o {os.path.join(input_dir, 'nocoding-pep.fasta')}")
    run_command(f"python coding_snvnonfra.py {os.path.join(input_dir, 'coding_snvnonfra.txt')} {os.path.join(input_dir, 'coding_snvnonfra.fasta')}")
    run_command(f"python coding_framstop.py {os.path.join(input_dir, 'coding_framstop.txt')} {os.path.join(input_dir, 'coding_framstop.fasta')}")

    logging.info("Full pipeline execution completed successfully!")

if __name__ == "__main__":
    main()

# Example usage:
# python run_full_pipeline.py -1 /path/sample_1.fastq -2 /path/sample_2.fastq -t 16 --filter_col1 12 --filter_col2 13 --threshold 0.05
