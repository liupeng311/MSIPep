#!/usr/bin/env python3
"""Root entry point — delegates to database_search/run_database_search_pipeline.py."""

from __future__ import annotations

import runpy
from pathlib import Path

if __name__ == "__main__":
    runpy.run_path(
        str(Path(__file__).resolve().parent / "database_search" / "run_database_search_pipeline.py"),
        run_name="__main__",
    )
