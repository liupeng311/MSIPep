import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor

def merge_mgf_files(mgf_files, merged_mgf_path):
    with open(merged_mgf_path, 'w') as outfile:
        for mgf_file in mgf_files:
            with open(mgf_file, 'r') as infile:
                outfile.write(infile.read())
    print(f"Merged {len(mgf_files)} MGF files into {merged_mgf_path}")

def run_comet(comet_path, params_file, mgf_path, output_path):
    try:
        result = subprocess.run([comet_path, f"-P{params_file}", mgf_path], check=True, capture_output=True, text=True)
        print(f"Comet ran successfully for {mgf_path}. Output written to {output_path}")
        print("Output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running Comet for {mgf_path}:")
        print(e.stderr)

def process_folder(comet_path, params_file, mgf_paths, output_dir):
    # Merge MGF files
    merged_mgf_path = os.path.join(output_dir, "merged.mgf")
    merge_mgf_files(mgf_paths, merged_mgf_path)

    # Set output file path
    output_path = os.path.join(output_dir, "comet_output.txt")

    # Run Comet
    run_comet(comet_path, params_file, merged_mgf_path, output_path)

def process_all_folders(comet_path, root_folder, params_file, num_threads):
    tasks = []
    # Iterate through all subdirectories in the root folder
    for subdir in os.listdir(root_folder):
        subdir_path = os.path.join(root_folder, subdir)
        if os.path.isdir(subdir_path):
            mgf_files = [os.path.join(subdir_path, f) for f in os.listdir(subdir_path) if f.endswith(".mgf")]
            if mgf_files:
                output_dir = subdir_path
                tasks.append((comet_path, params_file, mgf_files, output_dir))

    # Run tasks in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_folder, *task) for task in tasks]
        for future in futures:
            future.result()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process MGF files using Comet")
    parser.add_argument('--comet_path', required=True, help="Path to Comet executable")
    parser.add_argument('--root_folder', required=True, help="Root folder containing subdirectories with MGF files")
    parser.add_argument('--params_file', required=True, help="Path to Comet parameters file")
    parser.add_argument('--num_threads', type=int, default=16, help="Number of parallel threads (default: 16)")

    # Parse arguments
    args = parser.parse_args()

    # Process all folders with provided arguments
    process_all_folders(args.comet_path, args.root_folder, args.params_file, args.num_threads)

if __name__ == "__main__":
    main()

#python comet-run.py --comet_path ./software/comet/comet.2021010.linux.exe --root_folder /path/to/mgf --params_file ./software/comet/comet-TNBC.params --num_threads 16