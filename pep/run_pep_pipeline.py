#!/usr/bin/env python3
"""
run_pep_pipeline.py
-------------------
One-click immunopeptide screening pipeline (pep module).

Execution order (Methods §2.5):
  Step 1  blastp (1.blastp.py) → {sample}_classI_mut.fasta
  Step 2  NetMHCpan (2.netmhcpan_filter.py) → {sample}_SWB.fasta
  Step 3  IEDB (3.iedb_filter.py) → {sample}_IEDB.fasta
  Step 4  DeepImmuno (4.deepimmuno_filter.py) → {sample}_deepimmuno.fasta
  Step 5  ProTCR input (5.protcr_input.py) → {sample}_tcr_input.csv
  Step 6  ProTCR prediction (6.protcr_run.py)
  Step 7  ProTCR filter (7.protcr_filter.py) → plabels ≥ 0.9

Usage:
  python run_pep_pipeline.py --config pep_config.yaml
  python run_pep_pipeline.py --config pep_config.yaml --from-step protcr_run
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent

STEP_ORDER = [
    "blastp",
    "netmhcpan",
    "iedb",
    "deepimmuno",
    "protcr_input",
    "protcr_run",
    "protcr_filter",
]

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
            sys.exit("❌ Reading YAML requires: pip install pyyaml")
        return yaml.safe_load(text) or {}
    if suffix == ".json":
        import json
        return json.loads(text)
    sys.exit(f"❌ Unsupported config format: {suffix}")


def cfg_get(cfg: dict, *keys, default=None):
    cur: Any = cfg
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def paths_for_sample(cfg: dict) -> dict[str, Path]:
    sample = cfg["sample_name"]
    work = Path(cfg["work_dir"]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    protcr_cfg = cfg.get("protcr", {})
    results_dir = protcr_cfg.get("results_dir") or str(work)

    return {
        "work_dir": work,
        "classI_in": Path(cfg_get(cfg, "input", "classI_fasta", default="")),
        "classI_mut": work / f"{sample}_classI_mut.fasta",
        "blast_tsv": work / f"{sample}_classI_mut.tsv",
        "swb_fasta": work / f"{sample}_SWB.fasta",
        "iedb_fasta": work / f"{sample}_IEDB.fasta",
        "iedb_scores": work / f"{sample}_iedb_scores.tsv",
        "deepimmuno_fasta": work / f"{sample}_deepimmuno.fasta",
        "tcr_input": work / f"{sample}_tcr_input.csv",
        "deepimmuno_scores": work / f"{sample}_deepimmuno_scores.tsv",
        "protcr_results_dir": Path(results_dir).resolve(),
        "protcr_filtered_dir": work / "protcr_filtered",
        "protcr_final": work / "protcr_filtered" / f"{sample}_tcr_input_filtered.csv",
    }


def run_python(script: str, args: list[str], desc: str):
    cmd = [sys.executable, str(SCRIPT_DIR / script)] + args
    logger.info("▶ %s", desc)
    logger.info("  CMD: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    logger.info("✔ Done: %s", desc)


def file_ok(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _has_protcr_results(results_dir: Path) -> bool:
    if not results_dir.is_dir():
        return False
    try:
        import pandas as pd
    except ImportError:
        return any(results_dir.glob("*.csv"))

    for csv in results_dir.rglob("*.csv"):
        if csv.name.endswith("_tcr_input.csv") or "_filtered" in csv.stem:
            continue
        try:
            df = pd.read_csv(csv, sep=None, engine="python", nrows=3)
            cols = {c.strip().lower() for c in df.columns}
            if "peptide" in cols and "plabels" in cols:
                return True
        except Exception:
            continue
    return False


def step_blastp(cfg: dict, p: dict[str, Path]):
    if file_ok(p["classI_mut"]):
        logger.info("⏭ Skipping blastp (already exists: %s)", p["classI_mut"])
        return
    if not file_ok(p["classI_in"]):
        sys.exit(f"❌ Input FASTA not found: {p['classI_in']}")

    tools = cfg.get("tools", {})
    blastp_cfg = cfg.get("blastp", {})
    args = [
        "-i", str(p["classI_in"]),
        "-o", str(p["classI_mut"]),
        "--db", str(tools.get("blastp_db", "./software/ncbi-blast/test/human_uniprot")),
    ]
    if tools.get("blastp"):
        args.extend(["--blastp", str(tools["blastp"])])
    if blastp_cfg.get("keep_blast_tsv", True):
        args.extend(["--keep-blast", str(p["blast_tsv"])])
    run_python("1.blastp.py", args, "Step 1: blastp — remove peptides identical to human proteome")


def step_netmhcpan(cfg: dict, p: dict[str, Path]):
    if file_ok(p["swb_fasta"]):
        logger.info("⏭ Skipping NetMHCpan (already exists: %s)", p["swb_fasta"])
        return
    if not file_ok(p["classI_mut"]):
        sys.exit(f"❌ NetMHCpan input not found: {p['classI_mut']}")

    inp = cfg.get("input", {})
    tools = cfg.get("tools", {})
    nm = cfg.get("netmhcpan", {})
    hla = (inp.get("hla_alleles") or "").strip()

    args = [
        "-i", str(p["classI_mut"]),
        "-o", str(p["swb_fasta"]),
        "--netmhcpan", str(tools.get("netmhcpan", "netMHCpan")),
        "--lengths", str(nm.get("lengths", "8,9,10,11")),
        "--work-dir", str(p["work_dir"] / "netmhcpan_tmp"),
    ]
    if hla:
        args.extend(["--hla", hla])
    elif inp.get("optitype_dir"):
        args.extend(["--optitype-dir", str(inp["optitype_dir"])])
    else:
        sys.exit("❌ Specify HLA typing in input.hla_alleles or input.optitype_dir")
    run_python("2.netmhcpan_filter.py", args, "Step 2: NetMHCpan SB/WB filter (%Rank ≤0.5 / ≤2)")


def step_iedb(cfg: dict, p: dict[str, Path]):
    if file_ok(p["iedb_fasta"]):
        logger.info("⏭ Skipping IEDB (already exists: %s)", p["iedb_fasta"])
        return
    if not file_ok(p["swb_fasta"]):
        sys.exit(f"❌ IEDB input not found: {p['swb_fasta']}")

    iedb = cfg.get("iedb", {})
    tools = cfg.get("tools", {})
    args = [
        "-i", str(p["swb_fasta"]),
        "-o", str(p["iedb_fasta"]),
        "--iedb-script", str(tools.get("iedb_script", "./immunogenicity/new-predict_immunogenicity.py")),
        "--min-score", str(iedb.get("min_score", 0.0)),
        "--scores", str(p["iedb_scores"]),
    ]
    run_python("3.iedb_filter.py", args, "Step 3: IEDB immunogenicity (mask 2,3,C-term; score > 0)")


def step_deepimmuno(cfg: dict, p: dict[str, Path]):
    if file_ok(p["deepimmuno_fasta"]):
        logger.info("⏭ Skipping DeepImmuno (already exists: %s)", p["deepimmuno_fasta"])
        return
    di_input = p["iedb_fasta"] if file_ok(p["iedb_fasta"]) else p["swb_fasta"]
    if not file_ok(di_input):
        sys.exit(f"❌ DeepImmuno input not found: {di_input}")

    di = cfg.get("deepimmuno", {})
    tools = cfg.get("tools", {})
    args = [
        "-i", str(di_input),
        "-o", str(p["deepimmuno_fasta"]),
        "--threshold", str(di.get("threshold", 0.5)),
        "--deepimmuno-script", str(tools.get("deepimmuno_script", "deepimmuno-cnn.py")),
        "--scores", str(p["deepimmuno_scores"]),
    ]
    if di.get("default_hla"):
        args.extend(["--hla", str(di["default_hla"])])
    run_python("4.deepimmuno_filter.py", args, "Step 4: DeepImmuno filter (9/10mer score ≥ 0.5)")


def step_protcr_input(cfg: dict, p: dict[str, Path]):
    if file_ok(p["tcr_input"]):
        logger.info("⏭ Skipping ProTCR input prep (already exists: %s)", p["tcr_input"])
        return
    if not file_ok(p["deepimmuno_fasta"]):
        sys.exit(f"❌ ProTCR input dependency not found: {p['deepimmuno_fasta']}")
    run_python("5.protcr_input.py", ["-i", str(p["work_dir"])], "Step 5: generate ProTCR input CSV")


def step_protcr_run(cfg: dict, p: dict[str, Path]):
    protcr = cfg.get("protcr", {})
    install_dir = protcr.get("install_dir")
    if not install_dir:
        sys.exit("❌ Specify ProTCR install path in protcr.install_dir")

    results_dir = p["protcr_results_dir"]
    if _has_protcr_results(results_dir):
        logger.info("⏭ Skipping ProTCR prediction (result CSV with plabels already detected)")
        return

    if not file_ok(p["tcr_input"]):
        sys.exit(f"❌ ProTCR input CSV not found: {p['tcr_input']}")

    mode = protcr.get("mode", "single")
    args = ["--protcr-dir", str(install_dir)]
    if protcr.get("python"):
        args.extend(["--python", str(protcr["python"])])

    search_dirs = [str(p["work_dir"]), str(install_dir)]
    extra = protcr.get("result_search_dirs") or []
    search_dirs.extend(str(d) for d in extra)
    args.extend(["--result-search-dirs", *search_dirs])

    if mode == "batch":
        args.extend(["--batch-dir", str(p["work_dir"])])
        desc = "Step 6: ProTCR batch prediction (batch_run_protcr.sh)"
    else:
        args.extend(["-i", str(p["tcr_input"])])
        desc = "Step 6: ProTCR single-sample prediction (run_protcr.py)"

    run_python("6.protcr_run.py", args, desc)


def step_protcr_filter(cfg: dict, p: dict[str, Path]):
    protcr = cfg.get("protcr", {})
    threshold = protcr.get("plabels_threshold", 0.9)
    out_fmt = protcr.get("output_format", "csv")
    results_dir = p["protcr_results_dir"]

    if file_ok(p["protcr_final"]):
        logger.info("⏭ Skipping ProTCR filter (already exists: %s)", p["protcr_final"])
        return

    if not _has_protcr_results(results_dir):
        sys.exit(
            f"❌ No ProTCR result CSV found in {results_dir} (must contain peptide and plabels columns). "
            "Complete Step 6 first or check protcr.results_dir."
        )

    args = [
        str(results_dir),
        "--output_dir", str(p["protcr_filtered_dir"]),
        "--threshold", str(threshold),
        "--recursive",
        "--format", out_fmt,
    ]
    run_python("7.protcr_filter.py", args, f"Step 7: ProTCR filter (plabels≥{threshold})")


STEP_FUNCS = {
    "blastp": step_blastp,
    "netmhcpan": step_netmhcpan,
    "iedb": step_iedb,
    "deepimmuno": step_deepimmuno,
    "protcr_input": step_protcr_input,
    "protcr_run": step_protcr_run,
    "protcr_filter": step_protcr_filter,
}

STEP_FLAG = {
    "blastp": ("steps", "run_blastp"),
    "netmhcpan": ("steps", "run_netmhcpan"),
    "iedb": ("steps", "run_iedb"),
    "deepimmuno": ("steps", "run_deepimmuno"),
    "protcr_input": ("steps", "run_protcr_input"),
    "protcr_run": ("steps", "run_protcr_run"),
    "protcr_filter": ("steps", "run_protcr_filter"),
}


def parse_args():
    p = argparse.ArgumentParser(description="One-click immunopeptide screening pipeline (pep)")
    p.add_argument("--config", required=True, help="Config file (.yaml / .json)")
    p.add_argument("--from-step", choices=STEP_ORDER, default="blastp")
    p.add_argument("--to-step", choices=STEP_ORDER, default="protcr_filter")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)

    for key in ("sample_name", "work_dir"):
        if key not in cfg:
            sys.exit(f"❌ Missing config key: {key}")

    p = paths_for_sample(cfg)
    start = STEP_ORDER.index(args.from_step)
    end = STEP_ORDER.index(args.to_step)
    steps_run = STEP_ORDER[start : end + 1]

    logger.info("=" * 60)
    logger.info("Immunopeptide screening pipeline started")
    logger.info("  Sample : %s", cfg["sample_name"])
    logger.info("  Dir    : %s", p["work_dir"])
    logger.info("  Steps  : %s", " → ".join(steps_run))
    logger.info("=" * 60)

    for step in steps_run:
        flag = STEP_FLAG[step]
        if not cfg_get(cfg, *flag, default=True):
            logger.info("⏭ Skipping %s (steps.%s=false)", step, flag[1])
            continue
        STEP_FUNCS[step](cfg, p)

    logger.info("=" * 60)
    logger.info("✅ Pipeline complete! Main outputs:")
    logger.info("  Step1 blastp       : %s", p["classI_mut"])
    logger.info("  Step2 NetMHCpan    : %s", p["swb_fasta"])
    logger.info("  Step3 IEDB         : %s", p["iedb_fasta"])
    logger.info("  Step4 DeepImmuno   : %s", p["deepimmuno_fasta"])
    logger.info("  Step5 ProTCR input : %s", p["tcr_input"])
    logger.info("  Step6 ProTCR results: %s", p["protcr_results_dir"])
    logger.info("  Step7 ProTCR filter: %s", p["protcr_filtered_dir"])
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
