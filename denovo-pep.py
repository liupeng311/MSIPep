#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import shutil
import subprocess
import pandas as pd
import chardet


def detect_encoding(file_path):
    with open(file_path, 'rb') as f:
        result = chardet.detect(f.read())
        return result['encoding']


def read_tsv(file_path):
    encoding = detect_encoding(file_path)
    return pd.read_csv(file_path, sep='\t', encoding=encoding)


def write_fasta(sequences, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as file:
        for identifier, seq in sequences:
            file.write(f">{identifier}\n{seq}\n")


def filter_peptides(tsv_data, min_score=50, max_ppm_difference=5, min_avg_positional_score=0.7):
    filtered_sequences = []
    for i, row in tsv_data.iterrows():
        try:
            positional_scores = row['Positional Score']
            if isinstance(positional_scores, str):
                positional_scores = [float(x) for x in positional_scores.strip('[]').split(',')]
            else:
                positional_scores = []
            avg_positional_score = sum(positional_scores) / len(positional_scores) if positional_scores else 0
            if (
                row['Score'] >= min_score and
                abs(row['PPM Difference']) <= max_ppm_difference and
                avg_positional_score >= min_avg_positional_score and
                8 <= len(row['DENOVO']) <= 12
            ):
                filtered_sequences.append((f"{row['TITLE']}_{i + 1}", row['DENOVO']))
        except Exception as e:
            print(f"Warning: skipping row due to error: {e}")
    return filtered_sequences


def run_command(cmd):
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"✅ Finished: {cmd}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed: {cmd}\nError: {e}")


def process_file(input_file, output_dir):
    input_dir = os.path.dirname(input_file)
    base_name = os.path.splitext(os.path.basename(input_file))[0]

    if input_file.lower().endswith(".raw"):
        # Step 1: Convert raw -> mgf
        run_command(
            f"singularity exec --cleanenv --bind {input_dir}:/data ./software/msconvert.sif "
            f"wine msconvert /data/{base_name}.raw --mgf "
            f"--filter \"peakPicking true 1-\" --filter \"zeroSamples removeExtra\" --filter \"chargeStatePredictor\" "
            f"--outdir /data 2>/dev/null"
        )
        mgf_file = os.path.join(input_dir, f"{base_name}.mgf")
    elif input_file.lower().endswith(".mgf"):
        mgf_file = input_file
    else:
        print(f"Skipping unsupported file: {input_file}")
        return

    # Step 2: Add charge prediction again
    charged_mgf = os.path.join(input_dir, f"{base_name}_charge.mgf")
    run_command(
        f"singularity exec --cleanenv --bind {input_dir}:/data ./software/msconvert.sif "
        f"wine msconvert /data/{os.path.basename(mgf_file)} --mgf "
        f"--filter \"chargeStatePredictor\" "
        f"--outfile {base_name}_charge.mgf --outdir /data 2>/dev/null"
    )

    # Step 3: Run PepNet
    tsv_path = os.path.join(output_dir, f"{base_name}.tsv")
    run_command(
        f"python ./software/PepNet-master/denovo.py "
        f"--input {charged_mgf} --model ./software/PepNet-master/model.h5 "
        f"--output {tsv_path}"
    )


def main():
    parser = argparse.ArgumentParser(description="Batch pipeline for PepNet from raw/mgf to merged FASTA.")
    parser.add_argument("--input_dir", type=str, required=True, help="Input directory containing .raw or .mgf files.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory for TSV and final FASTA.")
    parser.add_argument("--output_fasta", type=str, required=True, help="Path to output merged FASTA file.")
    args = parser.parse_args()

    input_files = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir)
                   if f.lower().endswith((".raw", ".mgf"))]

    os.makedirs(args.output_dir, exist_ok=True)

    # Step 1-3: Process each raw/mgf file
    for file in input_files:
        process_file(file, args.output_dir)

    # Step 4: Filter all TSVs and merge to one FASTA
    all_sequences = []
    for f in os.listdir(args.output_dir):
        if f.endswith(".tsv"):
            tsv_path = os.path.join(args.output_dir, f)
            try:
                tsv_data = read_tsv(tsv_path)
                if "DENOVO" in tsv_data.columns:
                    filtered = filter_peptides(tsv_data)
                    all_sequences.extend(filtered)
            except Exception as e:
                print(f"Failed to parse {tsv_path}: {e}")

    # Write to merged FASTA
    write_fasta(all_sequences, args.output_fasta)
    print(f"✅ Final merged FASTA saved to {args.output_fasta}")


if __name__ == "__main__":
    main()

#python run_pepnet_batch.py  --input_dir /path/to/raw_or_mgf_files  --output_dir /path/to/output_folder  --output_fasta /path/to/output_folder/final_filtered_peptides.fasta
