#!/usr/bin/env python3
"""Root entry point — delegates to pep/run_pep_pipeline.py."""

from __future__ import annotations

import runpy
from pathlib import Path

if __name__ == "__main__":
    runpy.run_path(
        str(Path(__file__).resolve().parent / "pep" / "run_pep_pipeline.py"),
        run_name="__main__",
    )
