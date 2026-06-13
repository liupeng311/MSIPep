#!/bin/bash
# MSIPep Docker entrypoint — route commands to the correct conda environment.
set -euo pipefail

MSIPEP_ROOT="${MSIPEP_ROOT:-/opt/msipep}"
cd "${MSIPEP_ROOT}"

eval "$(conda shell.bash hook)"

usage() {
    cat <<'EOF'
MSIPep Docker container

Usage:
  docker run --rm -v /path/to/data:/data -v /path/to/ref:/ref IMAGE MODULE [args...]

Modules (new pipeline entry points):
  rnaseq          python run_rnaseq_pipeline.py --config CONFIG
  denovo          python run_denovo_pipeline.py --config CONFIG
  database_search python run_database_search_pipeline.py --config CONFIG
  pep             python run_pep_pipeline.py --config CONFIG
  pep-protcr      ProTCR steps only (protcr conda env)
  shell           interactive bash in msipep env
  help            show this message

Examples:
  docker run --rm -v /data:/data IMAGE rnaseq \
    --config /data/SAMPLE01/pipeline_config.yaml

  docker run --rm -v /data:/data IMAGE pep \
    --config /data/SAMPLE01/pep_config.yaml

  docker run --rm -v /data:/data IMAGE pep-protcr \
    --config /data/SAMPLE01/pep_config.yaml --from-step protcr_run

Volume layout (recommended):
  /data          FASTQ, MS raw/mgf, per-sample configs and outputs
  /ref           reference genome, STAR index, VEP cache (large; mount, do not bake in)
  /software      Zenodo bundles from start.sh (MSFragger, PepNet, netMHCpan, etc.)

Environment variables:
  MSIPEP_DATA=/data
  MSIPEP_REF=/ref
  MSIPEP_SOFTWARE=/software
  CUDA_VISIBLE_DEVICES   for ProTCR GPU on host
EOF
}

run_in_env() {
    local env_name="$1"
    shift
    conda activate "${env_name}"
    exec "$@"
}

MODULE="${1:-help}"
shift || true

case "${MODULE}" in
    help|-h|--help)
        usage
        ;;
    rnaseq)
        run_in_env msipep python "${MSIPEP_ROOT}/run_rnaseq_pipeline.py" "$@"
        ;;
    denovo)
        run_in_env pepnet python "${MSIPEP_ROOT}/run_denovo_pipeline.py" "$@"
        ;;
    database_search|db_search|search)
        run_in_env msipep python "${MSIPEP_ROOT}/run_database_search_pipeline.py" "$@"
        ;;
    pep|immunopep)
        run_in_env msipep python "${MSIPEP_ROOT}/run_pep_pipeline.py" "$@"
        ;;
    pep-protcr|protcr)
        # Steps 6–7 need protcr env (PyTorch); set protcr.python in config or run this module
        run_in_env protcr python "${MSIPEP_ROOT}/run_pep_pipeline.py" "$@"
        ;;
    shell|bash)
        run_in_env msipep bash "$@"
        ;;
    *)
        echo "Unknown module: ${MODULE}" >&2
        usage >&2
        exit 1
        ;;
esac
