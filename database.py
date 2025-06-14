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
    print(f"合并完成: {output_path}")

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
    parser = argparse.ArgumentParser(description="合并FASTA并配置MSFragger和Comet")
    parser.add_argument("--ms_dir", required=True, help="MGF和参考数据库目录")
    parser.add_argument("--fasta1", required=True, help="noncoding-pep.fasta 路径")
    parser.add_argument("--fasta2", required=True, help="coding_snvnonfra.fasta 路径")
    parser.add_argument("--fasta3", required=True, help="coding_framstop.fasta 路径")
    parser.add_argument("--mgf", required=True, help="MGF文件路径")
    args = parser.parse_args()

    # 合并FASTA
    all_fasta_files = [args.fasta1, args.fasta2, args.fasta3] + glob.glob(os.path.join(args.ms_dir, "*.fasta"))
    merged_fasta = os.path.join(args.ms_dir, "combined_protein_db.fasta")
    merge_fasta_files(all_fasta_files, merged_fasta)

    # 配置MSFragger参数
    msfragger_dir = os.path.join("../../software/MSFragger-3.8")
    msfragger_param = os.path.join(msfragger_dir, "closed_fragger.params")
    msfragger_param_new = os.path.join(msfragger_dir, "closed_fragger-n.params")
    update_param_file(msfragger_param, merged_fasta, msfragger_param_new)

    # 配置Comet参数
    comet_dir = os.path.join("../../software/comet")
    comet_param = os.path.join(comet_dir, "comet.params")
    comet_param_new = os.path.join(comet_dir, "comet-n.params")
    update_param_file(comet_param, merged_fasta, comet_param_new)

    # 运行搜索
    run_msfragger(os.path.join(msfragger_dir, "MSFragger-3.8.jar"), msfragger_param_new, args.mgf)
    run_comet(os.path.join(comet_dir, "comet.2021010.linux.exe"), comet_param_new, args.mgf)

if __name__ == "__main__":
    main()