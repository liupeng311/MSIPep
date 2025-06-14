import os
import sys
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
    return filtered_sequences

def run_command(cmd):
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"✅ Finished: {cmd}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed: {cmd}\nError: {e}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python run_pipeline_from_raw.py /path/to/data/ms/filename.raw")
        sys.exit(1)

    raw_file = sys.argv[1]
    if not os.path.isfile(raw_file):
        print(f"Error: file not found: {raw_file}")
        sys.exit(1)

    raw_dir = os.path.dirname(raw_file)
    base_name = os.path.splitext(os.path.basename(raw_file))[0]

    mgf1_file = os.path.join(raw_dir, f"{base_name}.mgf")
    mgf2_file = os.path.join(raw_dir, f"{base_name}_charge.mgf")
    tsv_file = os.path.join(raw_dir, f"{base_name}.tsv")
    fasta_file = os.path.join(raw_dir, f"{base_name}.fasta")

    # Step 1: Convert raw to mgf
    run_command(
        f"singularity exec --cleanenv --bind {raw_dir}:/data ../../software/msconvert.sif "
        f"wine msconvert /data/{base_name}.raw --mgf "
        f"--filter \"peakPicking true 1-\" --filter \"zeroSamples removeExtra\" --filter \"chargeStatePredictor\" "
        f"--outdir /data 2>/dev/null"
    )

    # Step 2: Re-process mgf for charge state prediction again
    run_command(
        f"singularity exec --cleanenv --bind {raw_dir}:/data ../../software/msconvert.sif "
        f"wine msconvert /data/{base_name}.mgf --mgf "
        f"--filter \"chargeStatePredictor\" "
        f"--outfile {base_name}_charge.mgf --outdir /data 2>/dev/null"
    )

    # Step 3: Run PepNet to generate TSV
    run_command(
        f"python ../../software/PepNet-master/denovo.py "
        f"--input {mgf2_file} --model ../../software/PepNet-master/model.h5 "
        f"--output {tsv_file}"
    )

    # Step 4: Filter TSV and write FASTA
    try:
        tsv_data = read_tsv(tsv_file)
        if 'DENOVO' in tsv_data.columns:
            filtered_sequences = filter_peptides(tsv_data)
            write_fasta(filtered_sequences, fasta_file)
            print(f"✅ FASTA written: {fasta_file}")
        else:
            print(f" DENOVO column not found in {tsv_file}")
    except Exception as e:
        print(f" Error processing {tsv_file}: {e}")

if __name__ == "__main__":
    main()

#python run_pipeline_from_raw.py /path/to/data/ms/yourfile.raw
