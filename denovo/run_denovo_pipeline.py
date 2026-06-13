#!/usr/bin/env python3

"""

run_denovo_pipeline.py

----------------------

One-click orchestration script for De Novo peptide identification



Pipeline logic:

  Step 1  PepNet de novo sequencing (denovo-pep.py)

          RAW/MGF → msconvert (format conversion + charge prediction) → PepNet → TSV

  Step 2  Filter TSV and export FASTA (2.pepnet-filter.py)

          Length filtering + confidence filtering (Score, PPM, Positional Score)



Usage:

  python run_denovo_pipeline.py --config denovo_config.yaml

  python run_denovo_pipeline.py --config denovo_config.yaml --from-step filter

  python run_denovo_pipeline.py --config denovo_config.yaml --mode filter_only

"""



from __future__ import annotations



import argparse

import logging

import subprocess

import sys

from pathlib import Path

from typing import Any



SCRIPT_DIR = Path(__file__).resolve().parent



STEP_ORDER = ["pepnet", "filter"]



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





def ensure_dir(path: str | Path):

    Path(path).mkdir(parents=True, exist_ok=True)





def file_ok(path: str | Path) -> bool:

    p = Path(path)

    return p.is_file() and p.stat().st_size > 0





def has_tsv(output_dir: Path) -> bool:

    return output_dir.is_dir() and any(output_dir.glob("*.tsv"))





def run_python(script: str, args: list[str], desc: str):

    cmd = [sys.executable, str(SCRIPT_DIR / script)] + args

    logger.info("▶ %s", desc)

    logger.info("  CMD: %s", " ".join(cmd))

    subprocess.run(cmd, check=True)

    logger.info("✔ Done: %s", desc)





def resolve_merged_fasta(cfg: dict) -> Path:

    out = Path(cfg["output_dir"]).resolve()

    name_tpl = cfg_get(cfg, "filter", "merged_fasta_name", default="{sample}_denovo_filtered.fasta")

    return out / name_tpl.format(sample=cfg["sample_name"])





def step_pepnet(cfg: dict):

    mode = cfg.get("mode", "full")

    if mode == "filter_only":

        logger.info("⏭ Skipping PepNet (mode=filter_only)")

        return



    out_dir = Path(cfg["output_dir"]).resolve()

    ensure_dir(out_dir)



    if has_tsv(out_dir):

        logger.info("⏭ Skipping PepNet (output directory already contains TSV: %s)", out_dir)

        return



    raw_mgf_dir = cfg_get(cfg, "input", "raw_mgf_dir")

    if not raw_mgf_dir or not Path(raw_mgf_dir).is_dir():

        sys.exit(f"❌ PepNet input directory not found: {raw_mgf_dir}")



    tools = cfg.get("tools", {})

    args = [

        "--input_dir", str(Path(raw_mgf_dir).resolve()),

        "--output_dir", str(out_dir),

        "--pepnet-only",

        "--msconvert_sif", tools.get("msconvert_sif", "./software/msconvert.sif"),

        "--pepnet_script", tools.get("pepnet_script", "./software/PepNet-master/denovo.py"),

        "--pepnet_model", tools.get("pepnet_model", "./software/PepNet-master/model.h5"),

        "--pepnet_python", tools.get("pepnet_python", "python"),

    ]

    run_python("denovo-pep.py", args, "Step 1: msconvert + PepNet de novo sequencing")





def step_filter(cfg: dict):

    out_dir = Path(cfg["output_dir"]).resolve()

    filt = cfg.get("filter", {})

    mode = cfg.get("mode", "full")



    if mode == "filter_only":

        tsv_inputs = cfg_get(cfg, "input", "tsv", default=[])

        if not tsv_inputs:

            sys.exit("❌ filter_only mode requires TSV paths in input.tsv")

        input_args = [str(x) for x in tsv_inputs]

    else:

        if not has_tsv(out_dir):

            sys.exit(f"❌ No TSV files found for filter step: {out_dir}")

        input_args = [str(out_dir)]



    output_mode = filt.get("output_mode", "merged")

    filter_args = [

        "-i", *input_args,

        "--min_length", str(filt.get("min_length", 8)),

        "--max_length", str(filt.get("max_length", 11)),

        "--min_score", str(filt.get("min_score", 0.85)),

        "--max_ppm_difference", str(filt.get("max_ppm_difference", 10.0)),

        "--min_avg_positional_score", str(filt.get("min_avg_positional_score", 0)),

    ]

    if filt.get("length_only", False):

        filter_args.append("--length-only")



    merged_fasta = resolve_merged_fasta(cfg)



    if output_mode == "per_file":

        filter_args.extend(["--per-file", "--output-dir", str(out_dir)])

        if file_ok(merged_fasta):

            logger.info("⏭ per_file mode will overwrite/update *-i.fasta under %s", out_dir)

    else:

        if file_ok(merged_fasta):

            logger.info("⏭ Skipping filter (already exists: %s)", merged_fasta)

            return

        filter_args.extend(["-o", str(merged_fasta)])



    run_python("2.pepnet-filter.py", filter_args, "Step 2: PepNet TSV filter → FASTA")





STEP_FUNCS = {

    "pepnet": step_pepnet,

    "filter": step_filter,

}



STEP_FLAG = {

    "pepnet": ("steps", "run_pepnet"),

    "filter": ("steps", "run_filter"),

}





def parse_args():

    p = argparse.ArgumentParser(description="One-click De Novo peptide identification pipeline")

    p.add_argument("--config", required=True, help="Config file path (.yaml / .json)")

    p.add_argument(

        "--from-step",

        choices=STEP_ORDER,

        default=None,

        help="Start from the specified step (resume from checkpoint)",

    )

    p.add_argument(

        "--to-step",

        choices=STEP_ORDER,

        default=None,

        help="Run through the specified step",

    )

    p.add_argument(

        "--mode",

        choices=["full", "filter_only"],

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

    start_step = args.from_step

    end_step = args.to_step



    if mode == "filter_only":

        steps_to_run = ["filter"] if start_step is None else STEP_ORDER[STEP_ORDER.index(start_step) :]

    else:

        start_idx = STEP_ORDER.index(start_step) if start_step else 0

        end_idx = STEP_ORDER.index(end_step) if end_step else len(STEP_ORDER) - 1

        steps_to_run = STEP_ORDER[start_idx : end_idx + 1]



    logger.info("=" * 60)

    logger.info("De Novo peptide pipeline started")

    logger.info("  Mode       : %s", mode)

    logger.info("  Sample     : %s", cfg["sample_name"])

    logger.info("  Output dir : %s", cfg["output_dir"])

    logger.info("  Steps      : %s", " → ".join(steps_to_run))

    logger.info("=" * 60)



    for step in steps_to_run:

        flag_keys = STEP_FLAG[step]

        if not cfg_get(cfg, *flag_keys, default=True):

            logger.info("⏭ Skipping step %s (config steps.%s=false)", step, flag_keys[1])

            continue

        STEP_FUNCS[step](cfg)



    out_dir = Path(cfg["output_dir"]).resolve()

    merged = resolve_merged_fasta(cfg)

    logger.info("=" * 60)

    logger.info("✅ Pipeline complete! Main outputs:")

    logger.info("  PepNet TSV dir   : %s", out_dir)

    if cfg_get(cfg, "filter", "output_mode", default="merged") == "merged":

        logger.info("  Filtered FASTA   : %s", merged)

    else:

        logger.info("  Filtered FASTA   : %s/*-i.fasta", out_dir)

    logger.info("=" * 60)





if __name__ == "__main__":

    main()

