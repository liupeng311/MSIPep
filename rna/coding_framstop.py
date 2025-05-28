import sys
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

"""
用法:
  python script.py <annovar.txt> <output.fa>

文件路径已在脚本中硬编码：
  TRANSCRIPT_MAP = "Seq_Ens.tsv"
  CDS_FILE        = "cds.fa"
参数:
  annovar.txt    ANNOVAR 注释结果（txt或tsv），需包含 "AAChange.refGene"、Chr、Start、End、Ref、Alt 列
  output.fa      输出 FASTA 文件，包含变异肽段
"""

# 硬编码文件路径
TRANSCRIPT_MAP = "../../reference/Seq_Ens.tsv"  # RefSeq->Ensembl 转录本对应表
CDS_FILE = "../../reference/Homo_sapiens.GRCh38.cds.all.fa"  # CDS FASTA 文件
UPSTREAM_AA = 24

# 1. 读取 RefSeq->Ensembl 转录本映射
def load_transcript_map(path):
    df = pd.read_csv(path, sep='\t', dtype=str)
    return {
        row['RefSeq mRNA ID'].split('.')[0]: row['Transcript stable ID'].split('.')[0]
        for _, row in df.iterrows()
        if pd.notna(row.get('RefSeq mRNA ID')) and pd.notna(row.get('Transcript stable ID'))
    }

# 2. 读取 CDS fasta，ID 去版本号
def load_cds(path):
    return {rec.id.split('.')[0]: str(rec.seq) for rec in SeqIO.parse(path, 'fasta')}

# 3. 解析 AAChange.refGene 字段
def parse_changes(field):
    if pd.isna(field): return []
    out = []
    for seg in field.split(','):
        parts = seg.split(':')
        if len(parts) < 5: continue
        gene = parts[0]
        refseq = parts[1].split('.')[0]
        cchg = parts[3].strip()
        pchg = parts[4].strip()
        out.append((gene, refseq, cchg, pchg))
    return out

# 4. 在 CDS 上应用 c. 突变
def apply_cds(seq, cchg):
    s = seq
    # SNV
    m = re.match(r'^c\.(\d+)([ACGT])>([ACGT])$', cchg)
    if m:
        i = int(m.group(1)) - 1
        return s[:i] + m.group(3) + s[i+1:]
    # 删除
    m = re.match(r'^c\.(\d+)(?:_(\d+))?del$', cchg)
    if m:
        start = int(m.group(1)) - 1
        end = int(m.group(2)) if m.group(2) else int(m.group(1))
        return s[:start] + s[end:]
    # 插入
    m = re.match(r'^c\.(?:\d+_)?(\d+)ins([ACGT]+)$', cchg)
    if m:
        pos = int(m.group(1))
        ins = m.group(2)
        return s[:pos] + ins + s[pos:]
    # 删除-插入
    m = re.match(r'^c\.(\d+)(?:_(\d+))?delins([ACGT]+)$', cchg)
    if m:
        start = int(m.group(1)) - 1
        end = int(m.group(2)) if m.group(2) else int(m.group(1))
        ins = m.group(3)
        return s[:start] + ins + s[end:]
    # 复制
    m = re.match(r'^c\.(\d+)(?:_(\d+))?dup([ACGT]*)$', cchg)
    if m:
        start = int(m.group(1)) - 1
        end = int(m.group(2)) if m.group(2) else int(m.group(1))
        dup_seq = m.group(3) if m.group(3) else s[start:end]
        return s[:end] + dup_seq + s[end:]
    return None

# 5. 翻译 CDS 到蛋白，至首个终止

def translate_cds(mut_seq):
    length = len(mut_seq) - (len(mut_seq) % 3)
    cds = mut_seq[:length]
    if len(cds) < 3:
        return ''
    prot = Seq(cds).translate(to_stop=True)
    return str(prot) + '*'

# 6. 从蛋白序列中截取含突变位点的片段
def extract_by_p(prot, pchg):
    m = re.match(r'^p\.\D+(\d+)', pchg)
    if not m:
        return prot
    aa_pos = int(m.group(1))
    idx = aa_pos - 1
    start = max(idx - UPSTREAM_AA, 0)
    return prot[start:]

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    annovar_file, out_fa = sys.argv[1:]

    # 加载映射表和CDS序列
    tmap = load_transcript_map(TRANSCRIPT_MAP)
    cds_dict = load_cds(CDS_FILE)
    df = pd.read_csv(annovar_file, sep='\t', dtype=str)

    results = []
    for _, row in df.iterrows():
        for gene, refseq, cchg, pchg in parse_changes(row.get('AAChange.refGene')):
            enst = tmap.get(refseq)
            if not enst:
                continue
            cds_seq = cds_dict.get(enst)
            if not cds_seq:
                continue
            mutated = apply_cds(cds_seq, cchg)
            if mutated is None:
                continue
            prot = translate_cds(mutated)
            if not prot or prot == '*':
                continue
            pep = extract_by_p(prot, pchg)
            header = (
                f">{gene}_{row['Chr']}_{row['Start']}_{row['End']}_{enst}_"
                f"{cchg}_{pchg}_{row['Ref']}/{row['Alt']}_len{len(pep)}"
            )
            results.extend([header, pep])

    if results:
        with open(out_fa, 'w') as fw:
            fw.write("\n".join(results))
        print("完成输出到", out_fa)
    else:
        print("未生成任何肽段，请检查输入数据。")
