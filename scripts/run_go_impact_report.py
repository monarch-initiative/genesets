#!/usr/bin/env python3
"""Compatibility wrapper for genesets-workflows go-impact."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_SRC = ROOT / "python" / "genesets-workflows" / "src"
if PACKAGE_SRC.exists():
    sys.path.insert(0, str(PACKAGE_SRC))

from genesets_workflows.reports.go_impact import main


if __name__ == "__main__":
    raise SystemExit(main())
