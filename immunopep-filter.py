import os
import re
import sys
import csv
import subprocess
from Bio import SeqIO

def run_netmhcpan(input_folder, hla_folder):
    folder_name = os.path.basename(input_folder)
    print("[NetMHCpan] Running peptide-HLA binding prediction...")

    # Step 1: 聚合所有候选肽段 fasta 文件（假设为 blastp-mut.fasta）
    fasta_file = os.path.join(input_folder, f"{folder_name}-blastp-mut.fasta")
    if not os.path.exists(fasta_file):
        print(f"Fasta file not found: {fasta_file}")
        return

    # Step 2: 按长度分类保存到新 fasta 文件中
    length_to_records = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        l = len(record.seq)
        length_to_records.setdefault(l, []).append(record)

    # Step 3: 获取 optitype 结果中的 HLA 类型
    optitype_result = next((f for f in os.listdir(hla_folder) if f.endswith("_result.tsv")), None)
    if not optitype_result:
        print(f"No OptiType result file found in {hla_folder}")
        return
    optitype_path = os.path.join(hla_folder, optitype_result)
    with open(optitype_path) as f:
        hla_line = f.readlines()[1].strip().split('\t')[1:7]  # 假设前 6 列是 HLA
    hla_string = ','.join(hla_line)

    # Step 4: 调用 netMHCpan 按长度运行
    for length, records in length_to_records.items():
        out_fa = os.path.join(input_folder, f"{folder_name}_output_{length}.fasta")
        SeqIO.write(records, out_fa, "fasta")
        out_txt = os.path.join(input_folder, f"{folder_name}-pepBA_{length}.txt")

        cmd = f"netMHCpan -f {out_fa} -a {hla_string} -l {length} -BA > {out_txt}"
        print(f"Running: {cmd}")
        subprocess.run(cmd, shell=True)

    print("[NetMHCpan] Finished.")

def step1_extract_swb(input_folder):
    folder_name = os.path.basename(input_folder)
    output_txt = os.path.join(input_folder, f"{folder_name}-SWB-pMHC.txt")
    output_fasta = os.path.join(input_folder, f"{folder_name}-SWB-pMHC.fasta")
    unique_sequences = set()
    line_number = 1

    with open(output_txt, 'w') as outfile, open(output_fasta, 'w') as fasta_file:
        custom_header = " Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL Score_BA %Rank_BA  Aff(nM) BindLevel\n"
        separator = '-' * len(custom_header)
        outfile.write(separator + "\n" + custom_header + separator + "\n")

        for filename in os.listdir(input_folder):
            if re.match(rf"^{re.escape(folder_name)}-pepBA.*_.*\\.txt$", filename):
                with open(os.path.join(input_folder, filename), 'r') as infile:
                    infile.readline()
                    for line in infile:
                        if re.search(r'\s*<=\s*(SB|WB)', line):
                            columns = line.split()
                            if len(columns) >= 3 and columns[2].strip():
                                peptide = columns[2]
                                if peptide not in unique_sequences:
                                    unique_sequences.add(peptide)
                                    fasta_file.write(f">Sequence_{line_number}\n{peptide}\n")
                                    outfile.write(line)
                                    line_number += 1
    print(f"[Step1] Extracted to: {output_txt} and {output_fasta}")

def step2_blastp_classify(input_folder):
    folder_name = os.path.basename(input_folder)
    input_file = os.path.join(input_folder, f"{folder_name}-SWB-pMHC.fasta")
    length_to_seqs = {}
    for record in SeqIO.parse(input_file, "fasta"):
        l = len(record.seq)
        length_to_seqs.setdefault(l, []).append(record)

    for length, records in length_to_seqs.items():
        tmp_fa = os.path.join(input_folder, f"{folder_name}-SWB-pMHC_{length}.fasta")
        SeqIO.write(records, tmp_fa, "fasta")
        blast_out = os.path.join(input_folder, f"{folder_name}-pMHC-blastp_{length}.txt")
        os.system(
            f"blastp -query {tmp_fa} -db /data/liup/lost_found/software/ncbi-blast/test/human_uniprot "
            f"-outfmt '6 qacc qseq sacc sseq evalue length pident' -evalue 1e9 -gapopen 11 -gapextend 1 > {blast_out}"
        )
        mut_out_txt = os.path.join(input_folder, f"{folder_name}-blastp-mut_{length}.txt")
        mut_out_fa = os.path.join(input_folder, f"{folder_name}-blastp-mut_{length}.fasta")
        extract_mutations(blast_out, mut_out_txt, mut_out_fa, length)

    merge_blastp_files(input_folder, folder_name)

def extract_mutations(blast_file, out_txt, out_fa, X):
    groups = {}
    with open(blast_file, 'r') as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) >= 7 and '-' not in cols[1]:
                key = cols[0]
                score = float(cols[5])
                if key not in groups or score > float(groups[key][5]):
                    groups[key] = cols

    selected = [v for v in groups.values() if str(v[5]) == str(X) and v[6] != '100.000']

    with open(out_txt, 'w') as f_txt, open(out_fa, 'w') as f_fa:
        for i, cols in enumerate(selected, 1):
            f_txt.write(' '.join(cols) + '\n')
            f_fa.write(f">Sequence_{i}\n{cols[1]}\n")

def merge_blastp_files(folder, prefix):
    with open(os.path.join(folder, f"{prefix}-blastp-mut.txt"), 'w') as final_txt:
        for file in os.listdir(folder):
            if file.startswith(f"{prefix}-blastp-mut_") and file.endswith(".txt"):
                with open(os.path.join(folder, file)) as f:
                    final_txt.write(f.read())
    with open(os.path.join(folder, f"{prefix}-blastp-mut.fasta"), 'w') as final_fa:
        for file in os.listdir(folder):
            if file.startswith(f"{prefix}-blastp-mut_") and file.endswith(".fasta"):
                with open(os.path.join(folder, file)) as f:
                    final_fa.write(f.read())

def step3_hla_binding(folder):
    fasta_file = os.path.join(folder, f"{os.path.basename(folder)}-blastp-mut.fasta")
    swb_txt = os.path.join(folder, f"{os.path.basename(folder)}-SWB-pMHC.txt")
    output_csv = os.path.join(folder, f"{os.path.basename(folder)}-mut_9_10-immuno.csv")

    peptides = []
    with open(fasta_file, 'r') as f:
        pep = ''
        for line in f:
            if line.startswith('>'):
                if pep:
                    peptides.append(pep)
                pep = ''
            else:
                pep += line.strip()
        if pep:
            peptides.append(pep)
    peptides = [p for p in peptides if len(p) in [9, 10]]

    hla_dict = {}
    with open(swb_txt, 'r') as f:
        for line in f.readlines()[3:]:
            cols = line.strip().split()
            if len(cols) >= 3:
                hla_dict[cols[2]] = cols[1].replace(':', '')

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        for pep in peptides:
            if pep in hla_dict:
                writer.writerow([pep, hla_dict[pep]])

    print(f"[Step3] Saved HLA info to {output_csv}")
    command = f"python3 deepimmuno-cnn.py --mode 'multiple' --intdir '{output_csv}' --outdir '{folder}'"
    subprocess.run(command, shell=True)

def step4_immunogenicity(folder):
    folder_name = os.path.basename(folder)
    fasta_file = os.path.join(folder, f"{folder_name}-blastp-mut.fasta")

    def fasta_to_txt(fa, out, lens):
        with open(fa) as f, open(out, 'w') as out_f:
            seq = ''
            for line in f:
                if line.startswith('>'):
                    if seq and len(seq) in lens:
                        out_f.write(seq + '\n')
                    seq = ''
                else:
                    seq += line.strip()
            if seq and len(seq) in lens:
                out_f.write(seq + '\n')

    file_8 = os.path.join(folder, f"{folder_name}_8-mut-immunogenicity.txt")
    file_9_11 = os.path.join(folder, f"{folder_name}_9_11-mut-immunogenicity.txt")
    fasta_to_txt(fasta_file, file_8, [8])
    fasta_to_txt(fasta_file, file_9_11, [9, 10, 11])

    out_8 = os.path.join(folder, f"{folder_name}-immunogenicity-pep_8.txt")
    out_9_11 = os.path.join(folder, f"{folder_name}-immunogenicity-pep_9_11.txt")

    subprocess.run(['python', '/data/liup/software/A-IEDB-tools/immunogenicity/new-predict_immunogenicity.py', '--output', out_8, file_8])
    subprocess.run(['python', '/data/liup/software/A-IEDB-tools/immunogenicity/new-predict_immunogenicity.py', '--custom_mask=2,3,9', '--output', out_9_11, file_9_11])

    merged = os.path.join(folder, f"{folder_name}-immunogenicity-pep.txt")
    with open(merged, 'w') as out:
        for file in [out_8, out_9_11]:
            with open(file) as f:
                lines = f.readlines()[4:]
                out.writelines(lines)
    print(f"[Step4] Immunogenicity results merged to {merged}")

def step5_tcr(folder):
    def read_cnn(file):
        with open(file, 'r') as f:
            next(f)
            return [(row.split('\t')[0], float(row.split('\t')[2])) for row in f if row.strip()]

    def read_iedb(file):
        with open(file, 'r') as f:
            reader = csv.reader(f)
            for _ in range(3): next(reader)
            header = next(reader)
            if header[0].lower() != 'peptide':
                yield (header[0], int(header[1]), float(header[2]))
            for row in reader:
                yield (row[0], int(row[1]), float(row[2]))

    folder_name = os.path.basename(folder)
    cnn = read_cnn(os.path.join(folder, "deepimmuno-cnn-result.txt"))
    pep9_11 = list(read_iedb(os.path.join(folder, f"{folder_name}-immunogenicity-pep_9_11.txt")))
    pep8 = list(read_iedb(os.path.join(folder, f"{folder_name}-immunogenicity-pep_8.txt")))

    out_csv = os.path.join(folder, f"{folder_name}-TCR.csv")
    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        for p, score in cnn:
            for q, l, s in pep9_11:
                if p == q and ((l in [9, 10] and score >= 0.5 and s >= 0) or (l == 11 and s >= 0)):
                    writer.writerow([p])
        for q, l, s in pep8:
            if s >= 0:
                writer.writerow([q])
    print(f"[Step5] TCR peptides saved to {out_csv}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python integrated_pipeline.py <input_folder> <hla_result_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    hla_folder = sys.argv[2]

    run_netmhcpan(input_folder, hla_folder)
    step1_extract_swb(input_folder)
    step2_blastp_classify(input_folder)
    step3_hla_binding(input_folder)
    step4_immunogenicity(input_folder)
    step5_tcr(input_folder)


#python integrated_pipeline.py /path/to/input_folder /path/to/hla_result_folder