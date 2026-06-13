#!/usr/bin/env python3
"""
6.protcr_run.py
---------------
Step 6: Invoke ProTCR for pMHC-TCR binding prediction.

Install ProTCR and dependencies in conda env protcr (see pep_config.example.yaml comments).

ProTCR internal pipeline (orchestrated by run_protcr.py):
  b1 → b2 → b3 → b4 → b5 → b6 → inference-mhc1.py → b7 → last.py

Usage (single file, same as official):
  python 5.protcr_run.py -i SAMPLE01_tcr_input.csv \\
      --protcr-dir /data/liup/software/protcr

Usage (batch, via batch_run_protcr.sh):
  python 5.protcr_run.py --batch-dir /data/output/SAMPLE01/pep \\
      --protcr-dir /data/liup/software/protcr
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path


def _resolve_python(protcr_cfg: dict) -> str:
    return protcr_cfg.get("python") or sys.executable


def run_protcr_single(
    input_csv: str,
    protcr_dir: str,
    python_exe: str,
    run_script: str = "run_protcr.py",
) -> None:
    install = Path(protcr_dir).resolve()
    script = install / run_script
    if not script.is_file():
        raise FileNotFoundError(f"ProTCR entry script not found: {script}")

    inp = Path(input_csv).resolve()
    if not inp.is_file():
        raise FileNotFoundError(f"Input CSV not found: {inp}")

    cmd = [python_exe, str(script), "-i", str(inp)]
    print(f"[ProTCR] Single-file mode")
    print(f"  CWD : {install}")
    print(f"  CMD : {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(install), check=False)
    if result.returncode != 0:
        raise RuntimeError(f"ProTCR run failed, exit code {result.returncode}")


def run_protcr_batch(
    batch_dir: str,
    protcr_dir: str,
    batch_script: str = "batch_run_protcr.sh",
) -> None:
    install = Path(protcr_dir).resolve()
    script = install / batch_script
    if not script.is_file():
        raise FileNotFoundError(f"Batch script not found: {script}")

    folder = Path(batch_dir).resolve()
    if not folder.is_dir():
        raise FileNotFoundError(f"Batch directory not found: {folder}")

    cmd = ["bash", str(script), str(folder)]
    print(f"[ProTCR] Batch mode")
    print(f"  CWD : {install}")
    print(f"  CMD : {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(install), check=False)
    if result.returncode != 0:
        raise RuntimeError(f"ProTCR batch run failed, exit code {result.returncode}")


def find_result_csvs(
    search_dirs: list[str],
    input_csv: str | None,
    patterns: list[str],
) -> list[str]:
    """Find result CSVs with plabels column in given directories (exclude input file itself)."""
    import pandas as pd

    input_abs = os.path.abspath(input_csv) if input_csv else None
    found: list[str] = []
    seen: set[str] = set()

    for d in search_dirs:
        if not d or not os.path.isdir(d):
            continue
        candidates: list[str] = []
        for pat in patterns:
            candidates.extend(glob.glob(os.path.join(d, pat)))
            candidates.extend(glob.glob(os.path.join(d, "**", pat), recursive=True))

        for path in sorted(set(candidates)):
            ap = os.path.abspath(path)
            if ap in seen or not os.path.isfile(ap):
                continue
            if input_abs and ap == input_abs:
                continue
            if ap.endswith("_tcr_input.csv"):
                continue
            try:
                df = pd.read_csv(ap, sep=None, engine="python", nrows=5)
                cols = {c.strip().lower() for c in df.columns}
                if "peptide" in cols and "plabels" in cols:
                    seen.add(ap)
                    found.append(ap)
            except Exception:
                continue
    return found


def parse_args():
    p = argparse.ArgumentParser(description="Run ProTCR prediction (pep pipeline Step 5)")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("-i", "--input", help="Single *_tcr_input.csv path")
    g.add_argument("--batch-dir", help="Directory with multiple tcr_input.csv files (calls batch_run_protcr.sh)")

    p.add_argument(
        "--protcr-dir",
        required=True,
        help="ProTCR install directory (contains run_protcr.py, batch_run_protcr.sh)",
    )
    p.add_argument("--python", default=None, help="Python from protcr env (default: current interpreter)")
    p.add_argument("--run-script", default="run_protcr.py", help="Single-file entry script name")
    p.add_argument("--batch-script", default="batch_run_protcr.sh", help="Batch shell script name")
    p.add_argument(
        "--result-search-dirs",
        nargs="*",
        default=None,
        help="Search these directories for result CSVs with plabels after run",
    )
    p.add_argument(
        "--result-patterns",
        default="*.csv",
        help="Result file glob patterns, comma-separated (default *.csv)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    python_exe = args.python or sys.executable
    patterns = [x.strip() for x in args.result_patterns.split(",") if x.strip()]

    if args.input:
        run_protcr_single(args.input, args.protcr_dir, python_exe, args.run_script)
        search_dirs = args.result_search_dirs or [
            os.path.dirname(os.path.abspath(args.input)),
            args.protcr_dir,
        ]
        input_csv = args.input
    else:
        run_protcr_batch(args.batch_dir, args.protcr_dir, args.batch_script)
        search_dirs = args.result_search_dirs or [args.batch_dir, args.protcr_dir]
        input_csv = None

    results = find_result_csvs(search_dirs, input_csv, patterns)
    if results:
        print(f"[ProTCR] Detected {len(results)} result CSV(s) with peptide/plabels:")
        for r in results:
            print(f"  → {r}")
    else:
        print(
            "[ProTCR] Run finished, but no result CSV with plabels was auto-detected.\n"
            "  Check ProTCR output path and set protcr.results_dir in config."
        )


if __name__ == "__main__":
    main()
