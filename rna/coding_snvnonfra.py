import os
import sys
import pandas as pd
import re
from Bio import SeqIO

REFERENCE_PROTEOME = "../../reference/Reference_protein_seq-RefSeq_ID.fasta"  # 替换为实际路径

def load_proteome(proteome_fasta):
    proteome_dict = {}
    for record in SeqIO.parse(proteome_fasta, "fasta"):
        header = record.description
        transcript_ids = re.findall(r'Transcript:\s*([^|)]+)', header)
        if transcript_ids:
            for transcript_id_group in transcript_ids:
                individual_ids = re.findall(r'(NM_[\d\.]+|XM_[\d\.]+)', transcript_id_group)
                for transcript_id in individual_ids:
                    base_transcript_id = transcript_id.split('.')[0]
                    proteome_dict[base_transcript_id] = str(record.seq)
    return proteome_dict

def generate_peptide(protein_seq, pos, upstream=24, downstream=24):
    mutation_index = pos - 1
    start = max(mutation_index - upstream, 0)
    end = mutation_index + downstream
    peptide = protein_seq[start:end]
    return peptide
def parse_AAChange_refGene(aa_change_field):
    mutations = []
    if pd.isna(aa_change_field):
        return mutations
    aa_changes = aa_change_field.split(',')
    for aa_change_info in aa_changes:
        parts = aa_change_info.split(':')
        if len(parts) < 5:
            continue
        gene, transcript_id, exon, c_change, p_change = parts[:5]
        mutations.append((gene, transcript_id, p_change.strip()))
    return mutations

def handle_mutation(gene, chr, start, end, transcript_id, exonic_func, ref, alt, p_change, proteome_dict):
    fasta_records = []
    if transcript_id not in proteome_dict:
        return fasta_records
    protein_seq = proteome_dict[transcript_id]
    p_change = p_change.strip()

    match_substitution = re.match(r'^p\.([A-Z\*])(\d+)([A-Z\*])$', p_change)
    match_single_del = re.match(r'^p\.([A-Z\*])(\d+)del$', p_change)
    match_multiple_del = re.match(r'^p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)del$', p_change)
    match_frameshift = re.match(r'^p\.([A-Z\*])(\d+)(?:([A-Z\*]+))?fs(?:\*(\d+))?$', p_change)
    match_delins_multiple = re.match(r'^p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)delins([A-Z\*\d]+)$', p_change)
    match_delins_single = re.match(r'^p\.([A-Z\*])(\d+)delins([A-Z\*\d]+)$', p_change)
    match_ins = re.match(r'^p\.([A-Z\*])(\d+)(?:_([A-Z\*])(\d+))?ins([A-Z\*\d]+)$', p_change)

    if match_substitution:
        ref_aa, pos, alt_aa = match_substitution.groups()
        pos = int(pos)
        if pos > len(protein_seq):
            return fasta_records
        peptide = generate_peptide(protein_seq, pos, upstream=24, downstream=24)
        peptide_start_pos = max(pos - 31, 0)  # 0基准
        mutation_relative_pos = pos - peptide_start_pos - 1  # 0基准索引
        if mutation_relative_pos < 0 or mutation_relative_pos >= len(peptide):
            return fasta_records
        peptide_list = list(peptide)
        if peptide_list[mutation_relative_pos] != ref_aa:
            return fasta_records
        peptide_list[mutation_relative_pos] = alt_aa
        peptide = ''.join(peptide_list)
        fasta_header = f">{gene}_{chr}_{start}_{end}_{transcript_id}_{exonic_func}_{ref}/{alt}"
        fasta_records.append(fasta_header)
        fasta_records.append(peptide)

    elif match_single_del:
        ref_aa, pos = match_single_del.groups()
        pos = int(pos)
        if pos > len(protein_seq):
            return fasta_records
        peptide = generate_peptide(protein_seq, pos, upstream=24, downstream=24)
        peptide_start_pos = max(pos - 31, 0)  # 0基准
        mutation_relative_pos = pos - peptide_start_pos - 1  # 0基准索引
        if mutation_relative_pos < 0 or mutation_relative_pos >= len(peptide):
            return fasta_records
        peptide_list = list(peptide)
        if peptide_list[mutation_relative_pos] != ref_aa:
            return fasta_records
        del peptide_list[mutation_relative_pos]
        peptide = ''.join(peptide_list)
        fasta_header = f">{gene}_{chr}_{start}_{end}_{transcript_id}_{exonic_func}_{ref}/{alt}"
        fasta_records.append(fasta_header)
        fasta_records.append(peptide)

    elif match_multiple_del:
        ref_aa1, pos1, ref_aa2, pos2 = match_multiple_del.groups()
        pos1 = int(pos1)
        pos2 = int(pos2)
        if pos1 > len(protein_seq) or pos2 > len(protein_seq):
            return fasta_records
        start_peptide_pos = max(pos1 - 31, 0)
        end_peptide_pos = pos2 + 30
        peptide = protein_seq[start_peptide_pos:end_peptide_pos]
        mutation_relative_start = pos1 - start_peptide_pos - 1  # 0基准索引
        mutation_relative_end = pos2 - start_peptide_pos  # 不包含在切片中
        if mutation_relative_start < 0 or mutation_relative_end > len(peptide):
            return fasta_records
        expected_segment = protein_seq[pos1 - 1:pos2]
        actual_segment = peptide[mutation_relative_start:mutation_relative_end]
        if expected_segment != actual_segment:
            return fasta_records
        peptide = peptide[:mutation_relative_start] + peptide[mutation_relative_end:]
        fasta_header = f">{gene}_{chr}_{start}_{end}_{transcript_id}_{exonic_func}_{ref}/{alt}"
        fasta_records.append(fasta_header)
        fasta_records.append(peptide)

    elif match_delins_multiple:
        ref_aa1, pos1, ref_aa2, pos2, alt_aa = match_delins_multiple.groups()
        pos1 = int(pos1)
        pos2 = int(pos2)
        if pos1 > len(protein_seq) or pos2 > len(protein_seq):
            return fasta_records
        start_peptide_pos = max(pos1 - 31, 0)
        end_peptide_pos = pos2 + 30
        peptide = protein_seq[start_peptide_pos:end_peptide_pos]
        mutation_relative_start = pos1 - start_peptide_pos - 1  # 0基准索引
        mutation_relative_end = pos2 - start_peptide_pos  # 不包含在切片中
        if mutation_relative_start < 0 or mutation_relative_end > len(peptide):
            return fasta_records
        expected_segment = protein_seq[pos1 - 1:pos2]
        actual_segment = peptide[mutation_relative_start:mutation_relative_end]
        if expected_segment != actual_segment:
            return fasta_records
        peptide = peptide[:mutation_relative_start] + alt_aa + peptide[mutation_relative_end:]
        if '*' in alt_aa:
            peptide = peptide.split('*')[0] + '*'
        fasta_header = f">{gene}_{chr}_{start}_{end}_{transcript_id}_{exonic_func}_{ref}/{alt}"
        fasta_records.append(fasta_header)
        fasta_records.append(peptide)

    elif match_delins_single:
        ref_aa, pos, alt_aa = match_delins_single.groups()
        pos = int(pos)
        if pos > len(protein_seq):
            return fasta_records
        peptide = generate_peptide(protein_seq, pos, upstream=24, downstream=24)
        peptide_start_pos = max(pos - 31, 0)
        mutation_relative_pos = pos - peptide_start_pos - 1
        if mutation_relative_pos < 0 or mutation_relative_pos >= len(peptide):
            return fasta_records
        peptide_list = list(peptide)
        if peptide_list[mutation_relative_pos] != ref_aa:
            return fasta_records
        peptide = peptide[:mutation_relative_pos + 1] + alt_aa + peptide[mutation_relative_pos + 1:]
        if '*' in alt_aa:
            peptide = peptide.split('*')[0] + '*'
        fasta_header = f">{gene}_{chr}_{start}_{end}_{transcript_id}_{exonic_func}_{ref}/{alt}"
        fasta_records.append(fasta_header)
        fasta_records.append(peptide)

    elif match_ins :
        try:
            if match_ins:
                ref_aa1, pos1, ref_aa2, pos2, inserted_aa = match_ins.groups()
                pos1 = int(pos1)
                pos2 = int(pos2)
                inserted_aa = inserted_aa.rstrip('*')
                has_stop_codon = p_change.endswith('*')
            if pos1 > len(protein_seq) or pos2 > len(protein_seq):
                return fasta_records
            start_peptide_pos = max(pos1 - 31, 0)
            end_peptide_pos = pos2 + 30
            peptide = protein_seq[start_peptide_pos:end_peptide_pos]
            insertion_relative_pos = pos2 - start_peptide_pos  # 插入在 pos2 之后

            if insertion_relative_pos < 0 or insertion_relative_pos > len(peptide):
                return fasta_records
            expected_ref_aa = protein_seq[pos2 - 1]
            if insertion_relative_pos == 0:
                actual_ref_aa = ''
            else:
                actual_ref_aa = peptide[insertion_relative_pos - 1]
            if actual_ref_aa != expected_ref_aa:
                return fasta_records
            peptide = peptide[:insertion_relative_pos] + inserted_aa + peptide[insertion_relative_pos:]
            if has_stop_codon:
                peptide = peptide.split('*')[0] + '*'

            fasta_header = f">{gene}_{chr}_{start}_{end}_{transcript_id}_{exonic_func}_{ref}/{alt}"
            fasta_records.append(fasta_header)
            fasta_records.append(peptide)

        except Exception as e:
            return fasta_records

    elif match_frameshift:
        ref_aa, pos, alt_aa, stop_pos = match_frameshift.groups()
        pos = int(pos)
        if pos > len(protein_seq):
            return fasta_records
        peptide_upstream = generate_peptide(protein_seq, pos, upstream=24, downstream=0)
        if alt_aa:
            peptide = peptide_upstream + alt_aa + '*'
        else:
            peptide = peptide_upstream + '*'
        fasta_header = f">{gene}_{chr}_{start}_{end}_{transcript_id}_{exonic_func}_{ref}/{alt}"
        fasta_records.append(fasta_header)
        fasta_records.append(peptide)

    return fasta_records

def main():
    if len(sys.argv) != 3:
        print("用法: python script2.py <annovar_annotation_file> <output_fasta>")
        sys.exit(1)
    annovar_file = sys.argv[1]
    output_fasta = sys.argv[2]
    output_dir = os.path.dirname(output_fasta)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    print("加载参考蛋白质组...")
    proteome_dict = load_proteome(REFERENCE_PROTEOME)
    print("读取ANNOVAR注释文件...")
    df = pd.read_csv(annovar_file, sep='\t', header=0)
    if 'AAChange.refGene' not in df.columns:
        print("错误: 注释文件中缺少 'AAChange.refGene' 列。请检查文件格式。")
        sys.exit(1)
    fasta_records = []
    print("处理突变信息...")
    for index, row in df.iterrows():
        exonic_func = row['ExonicFunc.refGene']
        aa_change_field = row['AAChange.refGene']
        mutations = parse_AAChange_refGene(aa_change_field)
        for gene, transcript_id, p_change in mutations:
            if p_change in ['UNKNOWN', '']:
                continue
            if p_change.startswith('p.M1?'):
                continue
            chr = row['Chr']
            start = row['Start']
            end = row['End']
            ref = row['Ref']
            alt = row['Alt']
            fasta_records.extend(
                handle_mutation(gene, chr, start, end, transcript_id, exonic_func, ref, alt, p_change, proteome_dict))
    with open(output_fasta, 'w') as f:
        for record in fasta_records:
            f.write(record + '\n')
    #print("输出文件已保存到", output_fasta)

if __name__ == "__main__":
    main()
