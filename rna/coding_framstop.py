import sys
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

"""
Usage:
  python script.py <annovar.txt> <output.fa>

Note: File paths are hardcoded in the script:
  TRANSCRIPT_MAP = "Seq_Ens.tsv"
  CDS_FILE       = "cds.fa"

Arguments:
  annovar.txt   ANNOVAR annotation result (txt or tsv), must contain columns: "AAChange.refGene", Chr, Start, End, Ref, Alt
  output.fa     Output FASTA file containing mutated peptide sequences
"""

# Hardcoded file paths
TRANSCRIPT_MAP = "../../reference/Seq_Ens.tsv"  # RefSeq to Ensembl transcript mapping table
CDS_FILE = "../../reference/Homo_sapiens.GRCh38.cds.all.fa"  # CDS FASTA file
UPSTREAM_AA = 24  # Number of upstream amino acids to keep when extracting peptide

# 1. Load RefSeq to Ensembl transcript mapping
def load_transcript_map(path):
    df = pd.read_csv(path, sep='\t', dtype=str)
    return {
        row['RefSeq mRNA ID'].split('.')[0]: row['Transcript stable ID'].split('.')[0]
        for _, row in df.iterrows()
        if pd.notna(row.get('RefSeq mRNA ID')) and pd.notna(row.get('Transcript stable ID'))
    }

# 2. Load CDS FASTA, removing version numbers from IDs
def load_cds(path):
    return {rec.id.split('.')[0]: str(rec.seq) for rec in SeqIO.parse(path, 'fasta')}

# 3. Parse the AAChange.refGene field
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

# 4. Apply c. mutation on CDS
def apply_cds(seq, cchg):
    s = seq
    # SNV
    m = re.match(r'^c\.(\d+)([ACGT])>([ACGT])$', cchg)
    if m:
        i = int(m.group(1)) - 1
        return s[:i] + m.group(3) + s[i+1:]
    # Deletion
    m = re.match(r'^c\.(\d+)(?:_(\d+))?del$', cchg)
    if m:
        start = int(m.group(1)) - 1
        end = int(m.group(2)) if m.group(2) else int(m.group(1))
        return s[:start] + s[end:]
    # Insertion
    m = re.match(r'^c\.(?:\d+_)?(\d+)ins([ACGT]+)$', cchg)
    if m:
        pos = int(m.group(1))
        ins = m.group(2)
        return s[:pos] + ins + s[pos:]
    # Deletion-insertion
    m = re.match(r'^c\.(\d+)(?:_(\d+))?delins([ACGT]+)$', cchg)
    if m:
        start = int(m.group(1)) - 1
        end = int(m.group(2)) if m.group(2) else int(m.group(1))
        ins = m.group(3)
        return s[:start] + ins + s[end:]
    # Duplication
    m = re.match(r'^c\.(\d+)(?:_(\d+))?dup([ACGT]*)$', cchg)
    if m:
        start = int(m.group(1)) - 1
        end = int(m.group(2)) if m.group(2) else int(m.group(1))
        dup_seq = m.group(3) if m.group(3) else s[start:end]
        return s[:end] + dup_seq + s[end:]
    return None

# 5. Translate CDS into protein, stopping at first stop codon
def translate_cds(mut_seq):
    length = len(mut_seq) - (len(mut_seq) % 3)
    cds = mut_seq[:length]
    if len(cds) < 3:
        return ''
    prot = Seq(cds).translate(to_stop=True)
    return str(prot) + '*'

# 6. Extract a peptide fragment containing the mutation site from the protein
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

    # Load transcript mapping and CDS sequences
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
        print("Output completed:", out_fa)
    else:
        print("No peptide generated. Please check your input data.")
