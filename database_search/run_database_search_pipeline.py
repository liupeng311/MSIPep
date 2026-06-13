#!/usr/bin/env python3

"""

run_database_search_pipeline.py

-------------------------------

One-click orchestration script for mass spectrometry database search



Pipeline logic:

  Step 1  Merge protein database FASTAs → combined_protein_db.fasta

  Step 2  Update parameters and run search engines (database_search.py)

          - MSFragger: process .raw → .tsv

          - Comet: process .mgf → .txt

  Step 3  Post-process and extract Class I peptides

          - msfragger-pep-I.py: TSV → fragger_classI.fasta

          - comet-pep.py: TXT → comet_classI.fasta / comet_all.fasta

  Step 4  Merge and deduplicate (built into database_search.py)

          → merged_classI_dedup.fasta



Scripts:

  database_search.py      Main pipeline (search + post-process + merge)

  msfragger-pep-I.py      MSFragger TSV post-processing

  comet-pep.py            Comet TXT post-processing

  peptide_filters.py      Shared filtering logic (decoy, peptide length, etc.)



Usage:

  python run_database_search_pipeline.py --config database_search_config.yaml

  python run_database_search_pipeline.py --config database_search_config.yaml --mode postprocess_only

"""



from __future__ import annotations



import argparse

import logging

import subprocess

import sys

from pathlib import Path

from typing import Any



SCRIPT_DIR = Path(__file__).resolve().parent



logging.basicConfig(

    level=logging.INFO,

    format="%(asctime)s [%(levelname)s] %(message)s",

    handlers=[logging.StreamHandler(sys.stdout)],

)

logger = logging.getLogger(__name__)





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

            sys.exit("❌ PyYAML is required to read YAML: pip install pyyaml")

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





def build_database_search_args(cfg: dict) -> list[str]:

    mode = cfg.get("mode", "full")

    inp = cfg.get("input", {})

    tools = cfg.get("tools", {})

    msf = cfg.get("msfragger_filter", {})

    comet = cfg.get("comet_filter", {})

    merge = cfg.get("merge", {})

    opts = cfg.get("options", {})

    steps = cfg.get("steps", {})



    args = ["--output_dir", str(Path(cfg["output_dir"]).resolve())]



    fasta_list = inp.get("fasta_list") or []

    if fasta_list:

        args.extend(["--fasta_list", *[str(x) for x in fasta_list]])



    raw_files = inp.get("raw") or []

    if raw_files:

        args.extend(["--raw", *[str(x) for x in raw_files]])



    mgf = inp.get("mgf")

    if mgf:

        args.extend(["--mgf", str(mgf)])



    args.extend([

        "--java_memory", str(tools.get("java_memory", "64G")),

        "--msfragger_dir", str(tools.get("msfragger_dir", "./software/MSFragger-3.8")),

        "--comet_dir", str(tools.get("comet_dir", "./software/comet")),

        "--comet_exe", str(tools.get("comet_exe", "")),

        "--msf_min_length", str(msf.get("min_length", 8)),

        "--msf_max_length", str(msf.get("max_length", 11)),

        "--msf_max_expectation", str(msf.get("max_expectation", 0.01)),

        "--msf_min_hyperscore", str(msf.get("min_hyperscore", 15)),

        "--comet_xcorr_threshold", str(comet.get("xcorr_threshold", 2.0)),

        "--comet_evalue_threshold", str(comet.get("evalue_threshold", 0.01)),

        "--comet_min_deltacn", str(comet.get("min_deltacn", 0.05)),

        "--comet_min_length", str(comet.get("min_length", 8)),

        "--comet_max_length", str(comet.get("max_length", 11)),

        "--class1-merge-headers", merge.get("header_mode", "join"),

    ])



    if opts.get("keep_decoy", False):

        args.append("--keep_decoy")

    if comet.get("strict_deltacn", False):

        args.append("--comet_strict_deltacn")

    if not merge.get("enabled", True) or not steps.get("run_merge_class1", True):

        args.append("--no-merge-class1")



    # Mode mapping and input validation

    if mode != "postprocess_only" and steps.get("run_search", True) and not fasta_list:

        sys.exit("❌ Search mode requires protein database FASTAs in input.fasta_list")



    if mode == "search_only" or not steps.get("run_postprocess", True):

        args.append("--search-only")

    elif mode == "postprocess_only":

        args.append("--postprocess-only")

    elif mode == "msfragger_only":

        if not raw_files:

            sys.exit("❌ msfragger_only mode requires RAW files in input.raw")

    elif mode == "comet_only":

        if not mgf:

            sys.exit("❌ comet_only mode requires MGF file in input.mgf")

    elif mode == "full":

        if not raw_files and not mgf:

            sys.exit("❌ full mode requires at least one of input.raw or input.mgf")

    else:

        sys.exit(f"❌ Unknown mode: {mode}")



    if not steps.get("run_search", True) and mode != "postprocess_only":

        if "--postprocess-only" not in args:

            args.append("--postprocess-only")



    return args





def parse_args():

    p = argparse.ArgumentParser(description="One-click mass spectrometry database search pipeline")

    p.add_argument("--config", required=True, help="Config file path (.yaml / .json)")

    p.add_argument(

        "--mode",

        choices=["full", "search_only", "postprocess_only", "msfragger_only", "comet_only"],

        default=None,

        help="Override mode in config file",

    )

    return p.parse_args()





def main():

    args = parse_args()

    cfg = load_config(args.config)



    if args.mode:

        cfg["mode"] = args.mode



    for key in ("sample_name", "output_dir"):

        if key not in cfg:

            sys.exit(f"❌ Required config field missing: {key}")



    mode = cfg.get("mode", "full")

    ds_args = build_database_search_args(cfg)



    logger.info("=" * 60)

    logger.info("Database search pipeline started")

    logger.info("  Mode       : %s", mode)

    logger.info("  Sample     : %s", cfg["sample_name"])

    logger.info("  Output dir : %s", cfg["output_dir"])

    logger.info("=" * 60)



    cmd = [sys.executable, str(SCRIPT_DIR / "database_search.py")] + ds_args

    logger.info("▶ Running main pipeline database_search.py")

    logger.info("  CMD: %s", " ".join(cmd))

    subprocess.run(cmd, check=True)



    out = Path(cfg["output_dir"]).resolve()

    logger.info("=" * 60)

    logger.info("✅ Pipeline complete! Main outputs:")

    logger.info("  Merged protein DB : %s/combined_protein_db.fasta", out)

    logger.info("  MSFragger         : %s/fragger_classI.fasta", out)

    logger.info("  Comet             : %s/comet_classI.fasta", out)

    logger.info("  Merged Class I    : %s/merged_classI_dedup.fasta", out)

    logger.info("=" * 60)





if __name__ == "__main__":

    main()

