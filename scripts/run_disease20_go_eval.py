#!/usr/bin/env python3
"""Run a MyGeneset manifest against official GOA human GO targets."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import subprocess
import sys
import time
from pathlib import Path

import prepare_go_eval


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=Path("evals/disease20/sets.tsv"), type=Path)
    parser.add_argument(
        "--out-dir", default=Path("evals/disease20_vs_go/generated"), type=Path
    )
    parser.add_argument("--binary", default=Path("target/release/genesets-rs"), type=Path)
    parser.add_argument("--variants", default="all,no_contributes_to,iba,iba_iea")
    parser.add_argument("--eval-name", default="disease20_vs_go")
    parser.add_argument(
        "--description",
        default=(
            "Twenty MyGeneset.info Disease Ontology human gene sets enriched "
            "against official GOA human GO annotations."
        ),
    )
    parser.add_argument(
        "--max-p-adjust",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff passed to genesets-rs configs.",
    )
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--skip-fetch", action="store_true")
    parser.add_argument("--skip-prepare", action="store_true")
    parser.add_argument(
        "--ontology-url",
        default=prepare_go_eval.DEFAULT_ONTOLOGY_URL,
        help="GO ontology URL passed to prepare_go_eval.py.",
    )
    parser.add_argument(
        "--gaf-url",
        default=prepare_go_eval.DEFAULT_GAF_URL,
        help="GOA GAF URL passed to prepare_go_eval.py.",
    )
    return parser.parse_args()


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def run_command(command: list[str], cwd: Path) -> dict:
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    elapsed = time.perf_counter() - started
    if completed.returncode != 0:
        print(completed.stdout, file=sys.stdout)
        print(completed.stderr, file=sys.stderr)
        raise SystemExit(f"command failed with exit code {completed.returncode}: {' '.join(command)}")
    return {
        "command": command,
        "runtime_seconds": round(elapsed, 3),
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
    }


def summarize_results(path: Path) -> dict:
    rows = 0
    significant = 0
    top: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            p_adjust = float(row["p_adjust_bonferroni"])
            if p_adjust < 0.05:
                significant += 1
            if len(top) < 10:
                top.append(
                    {
                        "query_id": row["query_id"],
                        "query_name": row["query_name"],
                        "target_id": row["target_id"],
                        "target_name": row["target_name"],
                        "overlap": int(row["overlap"]),
                        "query_size": int(row["query_size"]),
                        "target_size": int(row["target_size"]),
                        "background_size": int(row["background_size"]),
                        "p_value": row["p_value"],
                        "p_adjust_bonferroni": row["p_adjust_bonferroni"],
                    }
                )
    return {
        "path": str(path),
        "rows": rows,
        "bonferroni_significant_rows": significant,
        "top_rows": top,
    }


def main() -> int:
    args = parse_args()
    if not 0 <= args.max_p_adjust <= 1:
        raise SystemExit("--max-p-adjust must be between 0 and 1")
    repo = Path.cwd()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    variants = [name.strip() for name in args.variants.split(",") if name.strip()]

    commands: list[dict] = []
    if not args.skip_fetch:
        commands.append(
            run_command(
                [
                    sys.executable,
                    "scripts/fetch_mygeneset_eval.py",
                    "--manifest",
                    str(args.manifest),
                    "--out-dir",
                    str(args.out_dir),
                ],
                repo,
            )
        )

    if not args.skip_prepare:
        prepare_command = [
            sys.executable,
            "scripts/prepare_go_eval.py",
            "--out-dir",
            str(args.out_dir),
            "--variants",
            ",".join(variants),
            "--max-p-adjust",
            str(args.max_p_adjust),
            "--ontology-url",
            args.ontology_url,
            "--gaf-url",
            args.gaf_url,
        ]
        if args.force_download:
            prepare_command.append("--force-download")
        commands.append(run_command(prepare_command, repo))

    if not args.binary.exists():
        commands.append(run_command(["cargo", "build", "--release"], repo))

    result_summaries = {}
    for variant in variants:
        config = args.out_dir / variant / "config.yaml"
        result = args.out_dir / variant / "results.tsv"
        if result.exists():
            result.unlink()
        command_result = run_command([str(args.binary), "run", str(config)], repo)
        commands.append(command_result)
        result_summaries[variant] = summarize_results(result)
        result_summaries[variant]["timing"] = {
            "runtime_seconds": command_result["runtime_seconds"],
            "command": command_result["command"],
        }

    metadata = {
        "generated_at_utc": utc_now(),
        "eval": args.eval_name,
        "description": args.description,
        "result_filter": {
            "max_p_adjust": args.max_p_adjust,
            "correction": "bonferroni",
            "policy": "Result TSVs contain rows with Bonferroni-adjusted p-value <= max_p_adjust.",
        },
        "queries": {
            "manifest": str(args.manifest),
            "generated_gmt": str(args.out_dir / "queries.gmt"),
            "source": "MyGeneset.info",
        },
        "go_prep_metadata": str(args.out_dir / "metadata.yaml"),
        "go_inputs": {
            "ontology_url": args.ontology_url,
            "gaf_url": args.gaf_url,
        },
        "variants": result_summaries,
        "commands": commands,
    }

    with (args.out_dir / "run_metadata.yaml").open("w") as handle:
        prepare_go_eval.write_yaml_value(handle, metadata)
    with (args.out_dir / "run_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")

    for variant, summary in result_summaries.items():
        print(
            f"{variant}: {summary['rows']} rows, "
            f"{summary['bonferroni_significant_rows']} Bonferroni-significant"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
