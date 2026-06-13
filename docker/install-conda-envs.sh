#!/bin/bash
# Create MSIPep conda environments inside the Docker image (no Zenodo download).
set -euo pipefail

eval "$(conda shell.bash hook)"

create_env() {
    local name="$1"
    shift
    if conda env list | grep -qE "^${name}[[:space:]]"; then
        echo "[SKIP] env ${name} exists"
        return
    fi
    echo "[CREATE] env ${name}"
    conda create -n "${name}" -y "$@"
}

# Module 1 / 3 / 4 (blastp, netmhcpan orchestration)
create_env msipep python=3.10
conda activate msipep
conda install -y -c bioconda -c conda-forge \
    bwa samtools gatk4 picard star trimmomatic bcftools tabix \
    perl-yaml perl-json openjdk=11
pip install --no-cache-dir pyyaml biopython pandas numpy pysam openpyxl scipy
conda deactivate

# Module 2 — PepNet
create_env pepnet python=3.8
conda activate pepnet
pip install --no-cache-dir \
    "tensorflow>=2.5.0" "pandas>=0.20" pyteomics numba chardet openpyxl
conda deactivate

# Module 4 Step 6 — ProTCR (CPU wheels; mount GPU host if needed)
create_env protcr python=3.10
conda activate protcr
pip install --no-cache-dir torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install --no-cache-dir \
    transformers==4.30.2 accelerate pandas numpy scikit-learn tqdm matplotlib openpyxl biopython
conda install -y -c conda-forge sentencepiece scikit-learn
conda deactivate

# OptiType (optional HLA typing)
create_env optitype python=3.8
conda activate optitype
pip install --no-cache-dir numpy pandas pysam matplotlib future
conda install -y -c bioconda optitype || pip install --no-cache-dir optitype
conda deactivate

conda clean -afy
echo "[OK] All conda environments ready."
