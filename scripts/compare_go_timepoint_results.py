#!/usr/bin/env python3
"""Compare two generated GO timepoint result tables with genesets-rs."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


DEFAULT_LEFT = Path("evals/go_timepoints/generated/five_years_ago_2021_05_01/all/results.tsv")
DEFAULT_RIGHT = Path("evals/go_timepoints/generated/now_2026_03_25/all/results.tsv")
DEFAULT_OUTPUT = Path("notebooks/generated/go_2021_vs_2026.diff.parquet")
DEFAULT_METADATA = Path("notebooks/generated/go_2021_vs_2026.diff.yaml")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left", default=DEFAULT_LEFT, type=Path)
    parser.add_argument("--right", default=DEFAULT_RIGHT, type=Path)
    parser.add_argument("--output", default=DEFAULT_OUTPUT, type=Path)
    parser.add_argument("--metadata-output", default=DEFAULT_METADATA, type=Path)
    parser.add_argument("--p-adjust-cutoff", default=0.05, type=float)
    parser.add_argument("--binary", default="genesets-rs")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    missing = [path for path in [args.left, args.right] if not path.exists()]
    if missing:
        missing_text = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing timepoint result table(s): {missing_text}. "
            "Run scripts/run_go_timepoints_eval.py --skip-existing first, or regenerate the GO timepoint eval."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        args.binary,
        "compare",
        "--left",
        str(args.left),
        "--right",
        str(args.right),
        "--p-adjust-cutoff",
        str(args.p_adjust_cutoff),
        "--output-format",
        "parquet",
        "--output",
        str(args.output),
        "--metadata-output",
        str(args.metadata_output),
    ]
    subprocess.run(command, check=True)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.metadata_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
