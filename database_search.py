import os
import glob
import shutil
import argparse
import subprocess
from pathlib import Path

# Hardcoded paths to external scripts
COMET_PEP_SCRIPT = "./ms/comet-pep-I.py"
MSFRAGGER_PEP_SCRIPT = "./ms/msfragger-pep-I.py"

def merge_fasta_files(fasta_list, output_path):
    with open(output_path, 'w') as outfile:
        for fasta in fasta_list:
            if not os.path.exists(fasta):
                print(f"FASTA file does not exist: {fasta}")
                continue
            with open(fasta, 'r') as infile:
                shutil.copyfileobj(infile, outfile)
    print(f"FASTA merge completed: {output_path}")

def update_param_file(param_path, db_path, new_param_path):
    if not os.path.exists(param_path):
        raise FileNotFoundError(f"Parameter file not found: {param_path}")
    with open(param_path, 'r') as f:
        lines = f.readlines()

    with open(new_param_path, 'w') as f:
        for line in lines:
            if line.strip().startswith('database_name'):
                f.write(f"database_name = {db_path}\n")
            else:
                f.write(line)
    print(f"Parameter file updated: {new_param_path}")

def run_msfragger(jar_path, param_file, mgf_file):
    if not all(os.path.exists(f) for f in [jar_path, param_file, mgf_file]):
        raise FileNotFoundError("MSFragger executable, parameter file, or MGF file does not exist")
    cmd = f"java -jar {jar_path} {param_file} {mgf_file}"
    print(f"Running MSFragger: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def run_comet(exe_path, param_file, mgf_file, output_dir):
    if not all(os.path.exists(f) for f in [exe_path, param_file, mgf_file]):
        raise FileNotFoundError("Comet executable, parameter file, or MGF file does not exist")
    cmd = f"{exe_path} -P{param_file} {mgf_file}"
    print(f"Running Comet: {cmd}")
    subprocess.run(cmd, shell=True, check=True, cwd=output_dir)  # Ensure outputs are written to output_dir

def main():
    parser = argparse.ArgumentParser(description="Merge FASTA files, run MSFragger and Comet, and post-process results.")
    parser.add_argument("--fasta_dir1", required=True, help="Directory containing the first set of FASTA files")
    parser.add_argument("--fasta_list2", nargs='+', required=True, help="List of additional FASTA files")
    parser.add_argument("--mgf", required=True, help="Path to MGF file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    args = parser.parse_args()

    # Check if the hardcoded scripts exist
    for script in [COMET_PEP_SCRIPT, MSFRAGGER_PEP_SCRIPT]:
        if not os.path.exists(script):
            raise FileNotFoundError(f"Script file not found: {script}")

    # Create the output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Merge FASTA files
    fasta1_files = glob.glob(os.path.join(args.fasta_dir1, "*.fasta"))
    if not fasta1_files:
        raise FileNotFoundError(f"No FASTA files found in {args.fasta_dir1}")
    all_fasta_files = fasta1_files + args.fasta_list2
    merged_fasta = os.path.join(args.output_dir, "combined_protein_db.fasta")
    merge_fasta_files(all_fasta_files, merged_fasta)

    # MSFragger configuration
    msfragger_dir = "./software/MSFragger-3.8"
    msfragger_param = os.path.join(msfragger_dir, "closed_fragger.params")
    msfragger_param_new = os.path.join(msfragger_dir, "closed_fragger-n.params")
    update_param_file(msfragger_param, merged_fasta, msfragger_param_new)

    # Comet configuration
    comet_dir = "./software/comet"
    comet_param = os.path.join(comet_dir, "comet.params")
    comet_param_new = os.path.join(comet_dir, "comet-n.params")
    update_param_file(comet_param, merged_fasta, comet_param_new)

    # Run MSFragger
    run_msfragger(os.path.join(msfragger_dir, "MSFragger-3.8.jar"), msfragger_param_new, args.mgf)

    # Run Comet
    run_comet(os.path.join(comet_dir, "comet.2021010.linux.exe"), comet_param_new, args.mgf, args.output_dir)

    # Post-process MSFragger results
    fragger_output_dir = os.path.dirname(args.mgf)
    fragger_tsv_dir = os.path.join(args.output_dir, "fragger_tsv")
    os.makedirs(fragger_tsv_dir, exist_ok=True)
    for f in glob.glob(os.path.join(fragger_output_dir, "*.tsv")):
        shutil.move(f, fragger_tsv_dir)
    msfragger_cmd = ["python", MSFRAGGER_PEP_SCRIPT, "-i", fragger_tsv_dir,
                     "-o", os.path.join(args.output_dir, "fragger_classI.fasta")]
    print(f"Running MSFragger post-processing: {' '.join(msfragger_cmd)}")
    subprocess.run(msfragger_cmd, check=True)

    # Post-process Comet results
    comet_txt_files = glob.glob(os.path.join(args.output_dir, "*.txt"))
    if not comet_txt_files:
        raise FileNotFoundError(f"No Comet output files found in {args.output_dir}")
    comet_all = os.path.join(args.output_dir, "comet_all.fasta")
    comet_classI = os.path.join(args.output_dir, "comet_classI.fasta")
    comet_cmd = [
        "python", COMET_PEP_SCRIPT,
        "--input_files", *comet_txt_files,
        "--output_fasta", comet_all,
        "--output_combined_fasta", comet_classI,
        "--xcorr_threshold", "2.0",
        "--evalue_threshold", "0.01",
        "--threads", "4"
    ]
    print(f"Running Comet post-processing: {' '.join(comet_cmd)}")
    subprocess.run(comet_cmd, check=True)

if __name__ == "__main__":
    main()

# Example usage:
# python database_search.py --fasta_dir1 /path/to/fasta_group1 --fasta_list2 /path/to/extra1.fasta --mgf /path/to/input.mgf --output_dir /path/to/output_folder
