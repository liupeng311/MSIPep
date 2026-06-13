#!/usr/bin/env python3
"""Root entry point — delegates to denovo/run_denovo_pipeline.py."""

from __future__ import annotations

import runpy
from pathlib import Path

if __name__ == "__main__":
    runpy.run_path(
        str(Path(__file__).resolve().parent / "denovo" / "run_denovo_pipeline.py"),
        run_name="__main__",
    )
