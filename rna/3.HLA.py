import subprocess
import argparse
import gzip
import shutil
import os
import pathlib

# path
OPTITYPE_SCRIPT = "./software/OptiType/OptiTypePipeline.py"
REFERENCE_HLA = "./OptiType/data/hla_reference_rna.fasta"

def run_command(command):
    """run shell """
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)

def main(input_fastq, output_dir, output_prefix):
    os.makedirs(output_dir, exist_ok=True)

    # index
    uncompressed_fastq = os.path.join(output_dir, f"{output_prefix}.fastq")
    fished_bam = os.path.join(output_dir, f"{output_prefix}_RNA.bam")
    sorted_bam = os.path.join(output_dir, f"{output_prefix}_RNA_sorted.bam")
    fastq_output = os.path.join(output_dir, f"{output_prefix}_fished.fastq")
    bam_index = sorted_bam + ".bai"

    # Step 0: ungzip .fastq.gz
    with gzip.open(input_fastq, 'rb') as f_in, open(uncompressed_fastq, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    # Step 1: razers3 
    run_command(f"razers3 -i 95 -m 1 -dr 0 -o {fished_bam} {REFERENCE_HLA} {uncompressed_fastq}")

    # Step 1.1: BAM sort
    run_command(f"samtools sort {fished_bam} -o {sorted_bam}")
    run_command(f"samtools index {sorted_bam}")

    # Step 2: BAM to FASTQ
    run_command(f"samtools bam2fq {sorted_bam} > {fastq_output}")

    # Step 3: run OptiType
    run_command(f"python {OPTITYPE_SCRIPT} -i {fastq_output} -r -o {output_dir} -p {output_prefix} -v")

    # Step 4: remove 
    intermediate_files = [
        uncompressed_fastq,
        fished_bam,
        sorted_bam,
        bam_index,
        fastq_output,
    ]
    for file in intermediate_files:
        if os.path.exists(file):
            print(f"Deleting intermediate file: {file}")
            os.remove(file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run HLA typing pipeline for single-end FASTQ.gz data")
    parser.add_argument("-i", "--input_fastq", required=True, help="Path to input single-end FASTQ.gz file")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("-p", "--output_prefix", required=True, help="Prefix for output files")

    args = parser.parse_args()
    main(args.input_fastq, args.output_dir, args.output_prefix)

#python run_optitype_pipeline.py -i /path/to/sample.fastq.gz -o /path/to/output_dir  -p sample1