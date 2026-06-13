#!/usr/bin/env python3
"""
run_rnaseq_pipeline.py
----------------------
One-click orchestration: RNA-seq somatic mutations → neoantigen peptides

Execution order:
  Step 1  Mutation calling (one of three, determined by config.mode)
          - tumor_only_pe : 1.rna_mutation_pipeline.py (tumor paired-end STAR+Mutect2)
          - tumor_only_se : 1.rna_SE_mutation_pipeline.py (tumor single-end STAR+Mutect2)
          - tumor_normal  : 1.rna_normal_tumor_mutation.py (paired BWA+Mutect2+ANNOVAR)
  Step 2  VEP annotation (Wildtype + Frameshift plugins, for coding translation)
  Step 3  VCF filtering (2_VCF_filter_WGS.py)
  Step 4  VCF split coding/non-coding (3.vcf_split.py)
  Step 5a Coding peptides (5.coding_mut_pep.py, pVACtools algorithm)
  Step 5b Non-coding peptides (4.nocoding_mut_pep.py)

Usage:
  python run_rnaseq_pipeline.py --config pipeline_config.yaml
  python run_rnaseq_pipeline.py --config pipeline_config.yaml --from-step vep
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent

STEP_ORDER = [
    "mutation_calling",
    "vep",
    "vcf_filter",
    "vcf_split",
    "coding_peptides",
    "noncoding_peptides",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(path: str) -> dict[str, Any]:
    p = Path(path)
    if not p.is_file():
        sys.exit(f"❌ Config file not found: {path}")

    text = p.read_text(encoding="utf-8")
    suffix = p.suffix.lower()

    if suffix in (".yaml", ".yml"):
        try:
            import yaml  # type: ignore
        except ImportError:
            sys.exit("❌ PyYAML required to read YAML: pip install pyyaml")
        return yaml.safe_load(text) or {}

    if suffix == ".json":
        import json
        return json.loads(text)

    sys.exit(f"❌ Unsupported config format: {suffix} (use .yaml or .json)")


def cfg_get(cfg: dict, *keys, default=None):
    cur: Any = cfg
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def ensure_dir(path: str | Path):
    Path(path).mkdir(parents=True, exist_ok=True)


def file_ok(path: str | Path) -> bool:
    p = Path(path)
    return p.is_file() and p.stat().st_size > 0


# ---------------------------------------------------------------------------
# Subprocess execution
# ---------------------------------------------------------------------------

def run_python(script: str, args: list[str], desc: str):
    cmd = [sys.executable, str(SCRIPT_DIR / script)] + args
    logger.info("▶ %s", desc)
    logger.info("  CMD: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    logger.info("✔ Done: %s", desc)


def run_shell(command: str, desc: str):
    logger.info("▶ %s", desc)
    logger.info("  CMD: %s", command)
    shell = True
    executable = None
    if os.name != "nt":
        executable = "/bin/bash"
    subprocess.run(command, shell=shell, check=True, executable=executable)
    logger.info("✔ Done: %s", desc)


# ---------------------------------------------------------------------------
# Path derivation
# ---------------------------------------------------------------------------

def output_paths(cfg: dict) -> dict[str, Path]:
    mode = cfg["mode"]
    sample = cfg["sample_name"]
    out = Path(cfg["output_dir"]).resolve()
    ensure_dir(out)

    paths: dict[str, Path] = {"out_dir": out, "sample": Path(sample)}

    if mode in ("tumor_only_pe", "tumor_only_se"):
        paths["raw_vcf"] = out / f"{sample}_no_editing.vcf"
    elif mode == "tumor_normal":
        paths["raw_vcf"] = out / "somatic_mutations_pass.vcf"
    else:
        sys.exit(f"❌ Unknown mode: {mode}")

    vep_override = cfg_get(cfg, "vep", "input_vcf", default="")
    paths["vep_vcf"] = Path(vep_override) if vep_override else out / f"{sample}_vep_annotated.vcf"
    paths["filtered_vcf"] = out / f"{sample}_vep_annotated_filter.vcf"
    paths["coding_vcf"] = out / f"{sample}_vep_filter_coding.vcf"
    paths["noncoding_vcf"] = out / f"{sample}_vep_filter_noncoding.vcf"
    paths["coding_fasta_prefix"] = out / f"{sample}_coding_mut_pep"
    paths["noncoding_fasta"] = out / f"{sample}_noncoding_mut_pep.fasta"
    return paths


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def step_mutation_calling(cfg: dict):
    mode = cfg["mode"]
    sample = cfg["sample_name"]
    out = cfg["output_dir"]
    threads = str(cfg.get("threads", 16))
    ref = cfg_get(cfg, "reference", "genome_fasta")
    inp = cfg.get("input", {})

    common = [
        "--sample_name", sample,
        "--output_dir", out,
        "--threads", threads,
        "--ref", ref,
        "--star_index", cfg_get(cfg, "reference", "star_index"),
        "--dbsnp", cfg_get(cfg, "databases", "dbsnp"),
        "--gnomad", cfg_get(cfg, "databases", "gnomad"),
        "--common_biallelic", cfg_get(cfg, "databases", "common_biallelic"),
        "--pon", cfg_get(cfg, "databases", "pon"),
        "--redi_bed", cfg_get(cfg, "databases", "redi_bed"),
    ]

    if mode == "tumor_only_pe":
        run_python(
            "1.rna_mutation_pipeline.py",
            [
                "--fq1", inp["fq1"],
                "--fq2", inp["fq2"],
                "--trimmomatic_jar", cfg_get(cfg, "tools", "trimmomatic_jar"),
                "--trimmomatic_adapters", cfg_get(cfg, "tools", "trimmomatic_adapters_pe"),
            ] + common,
            "Step 1: Tumor paired-end RNA-seq mutation calling (STAR + Mutect2)",
        )
    elif mode == "tumor_only_se":
        run_python(
            "1.rna_SE_mutation_pipeline.py",
            [
                "--fq1", inp["fq1"],
                "--trimmomatic_jar", cfg_get(cfg, "tools", "trimmomatic_jar"),
                "--trimmomatic_adapters", cfg_get(cfg, "tools", "trimmomatic_adapters_se"),
            ] + common,
            "Step 1: Tumor single-end RNA-seq mutation calling (STAR + Mutect2)",
        )
    elif mode == "tumor_normal":
        run_python(
            "1.rna_normal_tumor_mutation.py",
            [
                "-tumor_rna_1P", inp["tumor_rna_1p"],
                "-tumor_rna_2P", inp["tumor_rna_2p"],
                "-normal_rna_1P", inp["normal_rna_1p"],
                "-normal_rna_2P", inp["normal_rna_2p"],
                "-output_dir", out,
                "-threads", threads,
                "--ref", ref,
                "--trimmomatic_jar", cfg_get(cfg, "tools", "trimmomatic_jar"),
                "--trimmomatic_adapters", cfg_get(cfg, "tools", "trimmomatic_adapters_pe"),
                "--dbsnp", cfg_get(cfg, "databases", "dbsnp"),
                "--mills_indels", cfg_get(cfg, "databases", "mills_indels"),
                "--annovar_dir", cfg_get(cfg, "tools", "annovar_dir"),
                "--humandb", cfg_get(cfg, "tools", "annovar_humandb"),
            ],
            "Step 1: Tumor-normal paired RNA-seq mutation calling (BWA + Mutect2 + ANNOVAR)",
        )
    else:
        sys.exit(f"❌ Unknown mode: {mode}")


def step_vep(cfg: dict, paths: dict[str, Path]):
    if cfg_get(cfg, "vep", "skip", default=False):
        logger.info("⏭ Skipping VEP (config.vep.skip=true)")
        return

    if file_ok(paths["vep_vcf"]):
        logger.info("⏭ Skipping VEP (already exists: %s)", paths["vep_vcf"])
        return

    if not file_ok(paths["raw_vcf"]):
        sys.exit(f"❌ VEP input VCF not found: {paths['raw_vcf']}")

    vep_bin = cfg_get(cfg, "tools", "vep", default="vep")
    cache = cfg_get(cfg, "tools", "vep_cache_dir")
    plugins = cfg_get(cfg, "tools", "vep_plugins_dir")
    species = cfg_get(cfg, "vep", "species", default="homo_sapiens")
    assembly = cfg_get(cfg, "vep", "assembly", default="GRCh38")
    ref = cfg_get(cfg, "reference", "genome_fasta")

    plugin_args = ""
    if plugins:
        wt_plugin = Path(plugins) / "Wildtype.pm"
        fs_plugin = Path(plugins) / "Frameshift.pm"
        if wt_plugin.is_file():
            plugin_args += f" --plugin Wildtype"
        if fs_plugin.is_file():
            plugin_args += f" --plugin Frameshift"

    cmd = (
        f'"{vep_bin}"'
        f' --input_file "{paths["raw_vcf"]}"'
        f' --output_file "{paths["vep_vcf"]}"'
        f' --format vcf --vcf --force_overwrite'
        f' --species {species} --assembly {assembly}'
        f' --offline --cache --dir_cache "{cache}"'
        f' --fasta "{ref}"'
        f' --symbol --canonical --tsl --biotype --hgvs --pick_order canonical,tsl,biotype'
        f'{plugin_args}'
    )
    run_shell(cmd, "Step 2: VEP annotation (Wildtype/Frameshift plugins)")


def step_vcf_filter(cfg: dict, paths: dict[str, Path]):
    inp = paths["vep_vcf"]
    out = paths["filtered_vcf"]
    if file_ok(out):
        logger.info("⏭ Skipping VCF filter (already exists: %s)", out)
        return
    if not file_ok(inp):
        sys.exit(f"❌ VCF filter input not found: {inp}")

    vf = cfg.get("vcf_filter", {})
    args = [
        str(inp), str(out),
        "--min_t_depth", str(vf.get("min_t_depth", 10)),
        "--min_t_alt_count", str(vf.get("min_t_alt_count", 5)),
        "--min_vaf", str(vf.get("min_vaf", 0.05)),
        "--max_gnomad_af", str(vf.get("max_gnomad_af", 0.01)),
        "--tumor_col", str(vf.get("tumor_col", 9)),
    ]
    run_python("2_VCF_filter_WGS.py", args, "Step 3: VCF quality and functional filtering")


def step_vcf_split(cfg: dict, paths: dict[str, Path]):
    inp = paths["filtered_vcf"]
    coding = paths["coding_vcf"]
    noncoding = paths["noncoding_vcf"]
    if file_ok(coding) and file_ok(noncoding):
        logger.info("⏭ Skipping VCF split (coding/non-coding files already exist)")
        return
    if not file_ok(inp):
        sys.exit(f"❌ VCF split input not found: {inp}")

    run_python(
        "3.vcf_split.py",
        [str(inp), str(coding), str(noncoding)],
        "Step 4: Split VCF into coding / non-coding regions",
    )


def step_coding_peptides(cfg: dict, paths: dict[str, Path]):
    inp = paths["coding_vcf"]
    prefix = paths["coding_fasta_prefix"]
    out_combined = Path(f"{prefix}.fasta")
    if file_ok(out_combined):
        logger.info("⏭ Skipping coding peptide generation (already exists: %s)", out_combined)
        return
    if not file_ok(inp):
        logger.info("⚠ Coding VCF empty or missing; skipping coding peptide generation")
        return

    pep = cfg.get("peptide", {})
    args = [
        str(inp), str(prefix),
        "--flanking-length", str(pep.get("flanking_length", 24)),
        "--downstream-length", str(pep.get("downstream_length", 0)),
        "--epitope-length", str(pep.get("epitope_length", 0)),
    ]
    if not pep.get("canonical_only", True):
        args.append("--no-canonical-only")

    run_python("5.coding_mut_pep.py", args, "Step 5a: Coding mutation peptides (pVACtools style)")


def step_noncoding_peptides(cfg: dict, paths: dict[str, Path]):
    inp = paths["noncoding_vcf"]
    out = paths["noncoding_fasta"]
    if file_ok(out):
        logger.info("⏭ Skipping non-coding peptide generation (already exists: %s)", out)
        return
    if not file_ok(inp):
        logger.info("⚠ Non-coding VCF empty or missing; skipping non-coding peptide generation")
        return

    pep = cfg.get("peptide", {})
    ref = cfg_get(cfg, "reference", "genome_fasta")
    gtf = cfg_get(cfg, "reference", "gtf")
    args = [
        str(inp), ref, gtf, str(out),
        "--pep_len", str(pep.get("noncoding_lengths", "8,9,10,11")),
        "--window", str(pep.get("noncoding_window", 100)),
    ]
    run_python("4.nocoding_mut_pep.py", args, "Step 5b: Non-coding mutation peptides")


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

STEP_FUNCS = {
    "mutation_calling": lambda c, p: step_mutation_calling(c),
    "vep": step_vep,
    "vcf_filter": step_vcf_filter,
    "vcf_split": step_vcf_split,
    "coding_peptides": step_coding_peptides,
    "noncoding_peptides": step_noncoding_peptides,
}

STEP_FLAG = {
    "mutation_calling": ("steps", "run_mutation_calling"),
    "vep": ("steps", "run_vep"),
    "vcf_filter": ("steps", "run_vcf_filter"),
    "vcf_split": ("steps", "run_vcf_split"),
    "coding_peptides": ("steps", "run_coding_peptides"),
    "noncoding_peptides": ("steps", "run_noncoding_peptides"),
}


def parse_args():
    p = argparse.ArgumentParser(description="RNA-seq neoantigen peptide one-click pipeline")
    p.add_argument("--config", required=True, help="Config file path (.yaml / .json)")
    p.add_argument(
        "--from-step",
        choices=STEP_ORDER,
        default="mutation_calling",
        help="Start from this step (for resume from checkpoint)",
    )
    p.add_argument("--to-step", choices=STEP_ORDER, default="noncoding_peptides")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)

    required = ["mode", "sample_name", "output_dir"]
    for key in required:
        if key not in cfg:
            sys.exit(f"❌ Config missing required field: {key}")

    mode = cfg["mode"]
    logger.info("=" * 60)
    logger.info("RNA-seq neoantigen pipeline started")
    logger.info("  Mode       : %s", mode)
    logger.info("  Sample     : %s", cfg["sample_name"])
    logger.info("  Output dir : %s", cfg["output_dir"])
    logger.info("=" * 60)

    paths = output_paths(cfg)
    start_idx = STEP_ORDER.index(args.from_step)
    end_idx = STEP_ORDER.index(args.to_step)

    for step in STEP_ORDER[start_idx : end_idx + 1]:
        flag_keys = STEP_FLAG[step]
        if not cfg_get(cfg, *flag_keys, default=True):
            logger.info("⏭ Skipping step %s (config steps.%s=false)", step, flag_keys[1])
            continue
        STEP_FUNCS[step](cfg, paths)

    logger.info("=" * 60)
    logger.info("✅ Pipeline completed! Main outputs:")
    logger.info("  VEP VCF        : %s", paths["vep_vcf"])
    logger.info("  Filtered VCF   : %s", paths["filtered_vcf"])
    logger.info("  Coding VCF     : %s", paths["coding_vcf"])
    logger.info("  Non-coding VCF : %s", paths["noncoding_vcf"])
    logger.info("  Coding FASTA   : %s.fasta", paths["coding_fasta_prefix"])
    logger.info("  Non-coding FASTA: %s", paths["noncoding_fasta"])
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
