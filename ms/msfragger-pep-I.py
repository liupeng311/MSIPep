import os
import glob
import pandas as pd
import argparse
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="从文件夹内多个 TSV 文件提取 (protein, peptide)，合并输出到单个 FASTA 文件，仅保留 9~11 长度肽段，去重并合并相同肽段的蛋白信息。"
    )
    parser.add_argument("-i", "--input_dir", required=True,
                        help="输入文件夹路径，包含若干 .tsv 文件")
    parser.add_argument("-o", "--output_fasta", required=True,
                        help="输出的 FASTA 文件路径")
    return parser.parse_args()

def main():
    args = parse_arguments()
    input_dir = args.input_dir
    output_fasta = args.output_fasta
    # 搜索文件夹下所有 .tsv 文件
    tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
    if not tsv_files:
        print(f"在文件夹 {input_dir} 中未找到任何 .tsv 文件，程序退出。")
        return
    print(f"在文件夹 {input_dir} 中找到 {len(tsv_files)} 个 TSV 文件。")
    print(f"结果将写入 {output_fasta}")
    # 使用一个字典, key=peptide, value=对应的一组 protein
    peptide_protein_map = defaultdict(set)
    # 逐个文件读取
    for tsv_file in tsv_files:
        print(f"处理文件: {tsv_file}")
        try:
            # 读取为字符串类型，避免混淆
            df = pd.read_csv(tsv_file, sep="\t", dtype=str)
        except Exception as e:
            print(f"读取文件 {tsv_file} 时出错：{e}")
            continue
        # 检查是否存在所需的列
        if "protein" not in df.columns or "peptide" not in df.columns:
            print(f"文件 {tsv_file} 不包含 'protein' 或 'peptide' 列，跳过。")
            continue
        # 遍历每条记录
        for _, row in df.iterrows():
            protein = row["protein"]
            peptide = row["peptide"]
            # 过滤空值或缺失
            if pd.isna(protein) or pd.isna(peptide):
                continue
            # 仅保留 9~11 长度的肽段
            if 9 <= len(peptide) <= 11:
                # 将该蛋白加入该肽段对应的集合
                peptide_protein_map[peptide].add(protein)

    # 现在开始写 FASTA
    with open(output_fasta, "w") as fasta_out:
        for peptide, protein_set in peptide_protein_map.items():
            # protein_set 是所有与该肽段匹配的蛋白# 合并到一个标题行，假设用 '|' 分隔蛋白# 也可根据需要改用其他分隔符
            merged_proteins = "|".join(sorted(protein_set))
            fasta_header = f">{merged_proteins}"
            fasta_out.write(f"{fasta_header}\n{peptide}\n")

    print(f"处理完成，FASTA 文件已生成：{output_fasta}")

if __name__ == "__main__":
    main()

#python extract_classI_peptides.py -i /path/to/tsv_folder -o /path/to/output.fasta
