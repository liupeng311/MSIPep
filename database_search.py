import os
import glob
import shutil
import argparse
import subprocess
from pathlib import Path

def merge_fasta_files(fasta_list, output_path):
    with open(output_path, 'w') as outfile:
        for fasta in fasta_list:
            with open(fasta, 'r') as infile:
                shutil.copyfileobj(infile, outfile)
    print(f"FASTA 合并完成: {output_path}")

def update_param_file(param_path, db_path, new_param_path):
    with open(param_path, 'r') as f:
        lines = f.readlines()

    with open(new_param_path, 'w') as f:
        for line in lines:
            if line.strip().startswith('database_name'):
                f.write(f"database_name = {db_path}\n")
            else:
                f.write(line)
    print(f"参数文件已保存: {new_param_path}")

def run_msfragger(jar_path, param_file, mgf_file):
    cmd = f"java -jar {jar_path} {param_file} {mgf_file}"
    print(f"运行 MSFragger: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def run_comet(exe_path, param_file, mgf_file):
    cmd = f"{exe_path} -P{param_file} {mgf_file}"
    print(f"运行 Comet: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="合并 FASTA 并配置 MSFragger 和 Comet 并执行结果后处理")
    parser.add_argument("--fasta_dir1", required=True, help="第一类 fasta 文件所在目录")
    parser.add_argument("--fasta_list2", nargs='+', required=True, help="第二类 fasta 文件列表")
    parser.add_argument("--mgf", required=True, help="MGF 文件路径")
    parser.add_argument("--output_dir", required=True, help="结果输出目录")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 求合并的 fasta 文件
    fasta1_files = glob.glob(os.path.join(args.fasta_dir1, "*.fasta"))
    fasta2_files = args.fasta_list2
    all_fasta_files = fasta1_files + fasta2_files
    merged_fasta = os.path.join(args.output_dir, "combined_protein_db.fasta")
    merge_fasta_files(all_fasta_files, merged_fasta)

    # MSFragger 参数配置
    msfragger_dir = os.path.join("../../software/MSFragger-3.8")
    msfragger_param = os.path.join(msfragger_dir, "closed_fragger.params")
    msfragger_param_new = os.path.join(msfragger_dir, "closed_fragger-n.params")
    update_param_file(msfragger_param, merged_fasta, msfragger_param_new)

    # Comet 参数配置
    comet_dir = os.path.join("../../software/comet")
    comet_param = os.path.join(comet_dir, "comet.params")
    comet_param_new = os.path.join(comet_dir, "comet-n.params")
    update_param_file(comet_param, merged_fasta, comet_param_new)

    # 执行 MSFragger
    run_msfragger(os.path.join(msfragger_dir, "MSFragger-3.8.jar"), msfragger_param_new, args.mgf)

    # 执行 Comet
    run_comet(os.path.join(comet_dir, "comet.2021010.linux.exe"), comet_param_new, args.mgf)

    # 处理 MSFragger 结果
    fragger_output_dir = os.path.dirname(args.mgf)  # 结果和 mgf 在同一目录
    fragger_tsv_dir = os.path.join(fragger_output_dir, "fragger_tsv")
    os.makedirs(fragger_tsv_dir, exist_ok=True)
    for f in glob.glob(os.path.join(fragger_output_dir, "*.tsv")):
        shutil.move(f, fragger_tsv_dir)
    subprocess.run(["python", "extract_classI_peptides.py", "-i", fragger_tsv_dir,
                    "-o", os.path.join(args.output_dir, "fragger_classI.fasta")])

    # 处理 Comet 结果
    comet_txt_files = glob.glob(os.path.join(fragger_output_dir, "*.txt"))
    comet_all = os.path.join(args.output_dir, "comet_all.fasta")
    comet_classI = os.path.join(args.output_dir, "comet_classI.fasta")
    subprocess.run([
        "python", "extract_comet_peptides.py",
        "--input_files", *comet_txt_files,
        "--output_fasta", comet_all,
        "--output_combined_fasta", comet_classI,
        "--xcorr_threshold", "2.0",
        "--evalue_threshold", "0.01",
        "--threads", "4"
    ])

if __name__ == "__main__":
    main()

#python database.py --fasta_dir1 /path/to/fasta_group1  --fasta_list2 /path/to/extra1.fasta /path/to/extra2.fasta  --mgf /path/to/input.mgf  --output_dir /path/to/output_folder
