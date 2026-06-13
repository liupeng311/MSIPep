#!/bin/bash
# =============================================================================
# MSIPep — software download, reference setup, and conda environment guide
#
# Usage (from repository root):
#   chmod +x start.sh
#   bash start.sh
#
# Recommended OS: Linux (CentOS 7 / Ubuntu 20.04+), 16+ CPU cores, 64+ GB RAM
# Prerequisites: conda/miniconda, wget, unzip, Java 8+, Perl 5.26+
# =============================================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

echo "=============================================="
echo " MSIPep setup — ${REPO_ROOT}"
echo "=============================================="

# -----------------------------------------------------------------------------
# 1. Check basic tools
# -----------------------------------------------------------------------------
for cmd in wget unzip conda python3; do
    if ! command -v "${cmd}" &>/dev/null; then
        echo "[WARN] '${cmd}' not found in PATH. Install it before running pipelines."
    fi
done

# -----------------------------------------------------------------------------
# 2. Download bundled software and reference data (Zenodo)
# -----------------------------------------------------------------------------
SOFTWARE_DIR="${REPO_ROOT}/software"
mkdir -p "${SOFTWARE_DIR}"
cd "${SOFTWARE_DIR}"

ZENODO_BASE="https://zenodo.org/records/15960309/files"

download_if_missing() {
    local name="$1"
    local url="$2"
    if [[ -f "${name}" ]]; then
        echo "[SKIP] ${name} already exists"
    else
        echo "[DOWNLOAD] ${name}"
        wget -O "${name}" "${url}"
    fi
}

download_if_missing "annovar.zip"              "${ZENODO_BASE}/annovar.zip?download=1"
download_if_missing "comet.zip"                "${ZENODO_BASE}/comet.zip?download=1"
download_if_missing "DeepImmuno.zip"           "${ZENODO_BASE}/DeepImmuno.zip?download=1"
download_if_missing "DIA-NN-2.2.0.zip"       "${ZENODO_BASE}/DIA-NN-2.2.0.zip?download=1"
download_if_missing "humandb.zip"              "${ZENODO_BASE}/humandb.zip?download=1"
download_if_missing "immunogenicity.zip"       "${ZENODO_BASE}/immunogenicity.zip?download=1"
download_if_missing "msconvert.zip"           "${ZENODO_BASE}/msconvert.zip?download=1"
download_if_missing "MSFragger-3.8.zip"        "${ZENODO_BASE}/MSFragger-3.8.zip?download=1"
download_if_missing "ncbi-blast.zip"         "${ZENODO_BASE}/ncbi-blast.zip?download=1"
download_if_missing "netMHCpan-4.1.zip"      "${ZENODO_BASE}/netMHCpan-4.1.zip?download=1"
download_if_missing "OptiType.zip"             "${ZENODO_BASE}/OptiType.zip?download=1"
download_if_missing "PepNet.zip"               "${ZENODO_BASE}/PepNet.zip?download=1"
download_if_missing "reference.zip"            "${ZENODO_BASE}/reference.zip?download=1"
download_if_missing "Trimmomatic.zip"          "${ZENODO_BASE}/Trimmomatic.zip?download=1"

if [[ ! -f "TPP_7.3.0-src.tgz" ]]; then
    echo "[DOWNLOAD] TPP_7.3.0-src.tgz"
    wget -O TPP_7.3.0-src.tgz \
        "https://sourceforge.net/projects/sashimi/files/Trans-Proteomic%20Pipeline%20%28TPP%29/TPP%20v7.3%20%28Trade%20Wind%29%20rev%200/TPP_7.3.0-src.tgz/download"
fi

echo "[UNZIP] software archives"
for z in *.zip; do
    [[ -f "${z}" ]] && unzip -o -q "${z}" || true
done

# Move reference databases to repo root
if [[ -d humandb ]]; then
    rm -rf "${REPO_ROOT}/humandb"
    mv humandb "${REPO_ROOT}/"
fi
if [[ -d reference ]]; then
    rm -rf "${REPO_ROOT}/reference"
    mv reference "${REPO_ROOT}/"
fi

cd "${REPO_ROOT}/reference" 2>/dev/null && unzip -o -q '*.zip' 2>/dev/null || true
cd "${REPO_ROOT}/humandb" 2>/dev/null && unzip -o -q '*.zip' 2>/dev/null || true

# -----------------------------------------------------------------------------
# 3. Conda: base bioinformatics environment (Module 1 — RNA-seq)
# -----------------------------------------------------------------------------
# Creates env 'msipep' with aligners, GATK, and Python deps for pipeline runners.
echo ""
echo "=============================================="
echo " Conda environment: msipep (RNA-seq / common)"
echo "=============================================="

if conda env list | grep -qE '^msipep\s'; then
    echo "[SKIP] conda env 'msipep' already exists"
else
    conda create -n msipep python=3.10 -y
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate msipep

conda install -y -c bioconda \
    bwa samtools gatk4 picard star trimmomatic bcftools tabix \
    perl-yaml perl-json 2>/dev/null || \
conda install -y bwa samtools gatk picard 2>/dev/null || true

pip install --upgrade pip
pip install pyyaml biopython pandas numpy pysam

# Build BWA index if reference genome is present
REF_FA="${REPO_ROOT}/reference/hg38.fa"
if [[ -f "${REF_FA}" ]] && [[ ! -f "${REF_FA}.bwt" ]]; then
    echo "[INDEX] BWA index for ${REF_FA}"
    bwa index "${REF_FA}"
fi

conda deactivate

# -----------------------------------------------------------------------------
# 4. Conda: PepNet environment (Module 2 — de novo)
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo " Conda environment: pepnet (de novo sequencing)"
echo "=============================================="

if conda env list | grep -qE '^pepnet\s'; then
    echo "[SKIP] conda env 'pepnet' already exists"
else
    conda create -n pepnet python=3.8 -y
fi

conda activate pepnet
pip install --upgrade pip
pip install tensorflow>=2.5.0 pandas>=0.20 pyteomics numba chardet openpyxl
conda deactivate

# -----------------------------------------------------------------------------
# 5. Conda: ProTCR environment (Module 4 — pep Step 5)
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo " Conda environment: protcr (pMHC-TCR binding)"
echo "=============================================="

if conda env list | grep -qE '^protcr\s'; then
    echo "[SKIP] conda env 'protcr' already exists"
else
    conda create -n protcr python=3.10 -y
fi

conda activate protcr
pip install --upgrade pip
# CUDA 11.8 wheels; adjust index URL for your GPU / CPU-only setup
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118 || \
    pip install torch torchvision torchaudio
pip install transformers==4.30.2 accelerate pandas numpy scikit-learn tqdm matplotlib openpyxl
conda install -y sentencepiece scikit-learn
pip install biopython
conda deactivate

# -----------------------------------------------------------------------------
# 6. Conda: OptiType environment (HLA typing, optional)
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo " Conda environment: optitype (HLA typing)"
echo "=============================================="

if conda env list | grep -qE '^optitype\s'; then
    echo "[SKIP] conda env 'optitype' already exists"
else
    conda create -n optitype python=3.8 -y
fi

conda activate optitype
pip install numpy pandas pysam matplotlib future
conda install -y -c bioconda optitype 2>/dev/null || pip install optitype
conda deactivate

# -----------------------------------------------------------------------------
# 7. External tools (install separately, not in Zenodo bundle)
# -----------------------------------------------------------------------------
cat <<'EOF'

==============================================
 External tools (user-provided)
==============================================
  - Ensembl VEP + Wildtype/Frameshift plugins  (Module 1 coding peptides)
  - Singularity / Apptainer for msconvert.sif    (Module 2 RAW conversion)
  - ProTCR source repo with run_protcr.py        (Module 4 Step 5)
  - SpectraST / DIA-NN                           (DIA data only)

==============================================
 Pipeline entry points (after editing *.yaml)
==============================================
  Module 1  RNA-seq:
    conda activate msipep
    python run_rnaseq_pipeline.py --config rnaseq_processing/pipeline_config.yaml

  Module 2  De novo:
    conda activate pepnet
    python run_denovo_pipeline.py --config denovo/denovo_config.yaml

  Module 3  Database search:
    conda activate msipep
    python run_database_search_pipeline.py --config database_search/database_search_config.yaml

  Module 4  Immunopeptide filter:
    conda activate msipep          # Steps 1-4
    python run_pep_pipeline.py --config pep/pep_config.yaml
    conda activate protcr          # Step 5 (or set protcr.python in config)
    python run_pep_pipeline.py --config pep/pep_config.yaml --from-step protcr_run

==============================================
 Config templates (copy and edit paths)
==============================================
  rnaseq_processing/pipeline_config.example.yaml
  denovo/denovo_config.example.yaml
  database_search/database_search_config.example.yaml
  pep/pep_config.example.yaml

Setup complete.
EOF
