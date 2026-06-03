#!/usr/bin/env python3
"""Run Expression20 queries against a prepared official GO snapshot."""

from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import sys
import time
from pathlib import Path

import prepare_go_eval


DEFAULT_GO_DIR = Path("evals/go_timepoints/generated/now_2026_03_25")
DEFAULT_QUERIES = Path("evals/expression20/generated/queries.gmt")
DEFAULT_OUTPUT = Path("notebooks/generated/expression20_go_now.parquet")
DEFAULT_METADATA = Path("notebooks/generated/expression20_go_now.yaml")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--go-dir", default=DEFAULT_GO_DIR, type=Path)
    parser.add_argument("--queries", default=DEFAULT_QUERIES, type=Path)
    parser.add_argument("--output", default=DEFAULT_OUTPUT, type=Path)
    parser.add_argument("--metadata-output", default=DEFAULT_METADATA, type=Path)
    parser.add_argument("--binary", default="genesets-rs")
    parser.add_argument("--min-overlap", default=2, type=int)
    parser.add_argument("--max-p-adjust", default=0.05, type=float)
    parser.add_argument("--skip-existing", action="store_true")
    return parser.parse_args()


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def require(path: Path, description: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {description}: {path}")


def main() -> int:
    args = parse_args()
    required = {
        "GO annotations": args.go_dir / "all/gene_terms.tsv",
        "GO term names": args.go_dir / "terms.tsv",
        "GO closure": args.go_dir / "closure.tsv",
        "GO background": args.go_dir / "background_all_goa_symbols.txt",
        "Expression20 queries": args.queries,
    }
    for description, path in required.items():
        require(path, description)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        args.binary,
        "matrix",
        "--annotations",
        str(args.go_dir / "all/gene_terms.tsv"),
        "--terms",
        str(args.go_dir / "terms.tsv"),
        "--closure",
        str(args.go_dir / "closure.tsv"),
        "--queries",
        str(args.queries),
        "--query-format",
        "gmt",
        "--background",
        str(args.go_dir / "background_all_goa_symbols.txt"),
        "--min-overlap",
        str(args.min_overlap),
        "--max-p-adjust",
        str(args.max_p_adjust),
        "--output-format",
        "parquet",
        "--output",
        str(args.output),
    ]

    started = time.perf_counter()
    if args.skip_existing and args.output.exists():
        status = "reused_existing"
        runtime_seconds = 0.0
    else:
        subprocess.run(command, check=True)
        status = "generated"
        runtime_seconds = round(time.perf_counter() - started, 3)

    metadata = {
        "generated_at_utc": utc_now(),
        "status": status,
        "runtime_seconds": runtime_seconds,
        "query_panel": "Expression20",
        "target": "GOA human all-evidence annotations with GO closure",
        "inputs": {
            "go_dir": str(args.go_dir),
            "annotations": str(args.go_dir / "all/gene_terms.tsv"),
            "terms": str(args.go_dir / "terms.tsv"),
            "closure": str(args.go_dir / "closure.tsv"),
            "background": str(args.go_dir / "background_all_goa_symbols.txt"),
            "queries": str(args.queries),
        },
        "result_filter": {
            "min_overlap": args.min_overlap,
            "max_p_adjust": args.max_p_adjust,
            "correction": "bonferroni",
        },
        "output": str(args.output),
        "command": command,
    }
    with args.metadata_output.open("w") as handle:
        prepare_go_eval.write_yaml_value(handle, metadata)

    print(f"{status}: {args.output}")
    print(f"metadata: {args.metadata_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
