#!/usr/bin/env python3
"""Run Disease20-vs-GO across pinned GO archive timepoints."""

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


DEFAULT_TIMEPOINTS = [
    {
        "label": "now_2026_03_25",
        "release_date": "2026-03-25",
        "description": "Current GO release used for this eval run.",
    },
    {
        "label": "five_years_ago_2021_05_01",
        "release_date": "2021-05-01",
        "description": "Closest pinned release to five years before 2026-05-11.",
    },
    {
        "label": "ten_years_ago_2016_05_01",
        "release_date": "2016-05-01",
        "description": "Closest pinned release to ten years before 2026-05-11.",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default=Path("evals/go_timepoints/generated"), type=Path)
    parser.add_argument("--manifest", default=Path("evals/disease20/sets.tsv"), type=Path)
    parser.add_argument("--variants", default="all,no_contributes_to,iba,iba_iea")
    parser.add_argument("--max-p-adjust", type=float, default=0.05)
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    return parser.parse_args()


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def run_command(command: list[str]) -> dict:
    started = time.perf_counter()
    completed = subprocess.run(command, text=True, capture_output=True)
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


def release_url(date: str, path: str) -> str:
    return f"https://release.geneontology.org/{date}/{path}"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def summarize_variant(result_path: Path) -> dict:
    rows = 0
    queries: set[str] = set()
    terms: set[str] = set()
    by_query: dict[str, int] = {}
    with result_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            queries.add(row["query_id"])
            terms.add(row["target_id"])
            by_query[row["query_id"]] = by_query.get(row["query_id"], 0) + 1
    return {
        "rows": rows,
        "query_count": len(queries),
        "unique_target_terms": len(terms),
        "rows_per_query": by_query,
    }


def build_pair_index(path: Path) -> dict[tuple[str, str], dict]:
    index: dict[tuple[str, str], dict] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for rank, row in enumerate(reader, start=1):
            p_adjust = float(row["p_adjust_bonferroni"]) if row["p_adjust_bonferroni"] != "0" else 0.0
            p_value = float(row["p_value"]) if row["p_value"] != "0" else 0.0
            index[(row["query_id"], row["target_id"])] = {
                "rank": rank,
                "query_name": row["query_name"],
                "target_name": row["target_name"],
                "overlap": int(row["overlap"]),
                "target_size": int(row["target_size"]),
                "p_value": p_value,
                "p_adjust": p_adjust,
            }
    return index


def neg_log10(value: float) -> float:
    if value <= 0:
        return 320.0
    import math

    return -math.log10(value)


def diff_pair_indexes(left: dict, right: dict, limit: int = 25) -> dict:
    shared = set(left).intersection(right)
    left_only = set(left).difference(right)
    right_only = set(right).difference(left)
    deltas = []
    for key in shared:
        l = left[key]
        r = right[key]
        deltas.append(
            {
                "query_id": key[0],
                "target_id": key[1],
                "query_name": r["query_name"] or l["query_name"],
                "target_name": r["target_name"] or l["target_name"],
                "delta_neg_log10_p_adjust": round(neg_log10(r["p_adjust"]) - neg_log10(l["p_adjust"]), 6),
                "left_p_adjust": l["p_adjust"],
                "right_p_adjust": r["p_adjust"],
                "left_overlap": l["overlap"],
                "right_overlap": r["overlap"],
                "delta_overlap": r["overlap"] - l["overlap"],
                "left_target_size": l["target_size"],
                "right_target_size": r["target_size"],
                "delta_target_size": r["target_size"] - l["target_size"],
            }
        )
    deltas.sort(key=lambda item: abs(item["delta_neg_log10_p_adjust"]), reverse=True)
    lost_examples = sorted(left_only)[:limit]
    gained_examples = sorted(right_only)[:limit]
    return {
        "shared_significant_pairs": len(shared),
        "lost_significant_pairs": len(left_only),
        "gained_significant_pairs": len(right_only),
        "lost_significant_examples": [
            {
                "query_id": query_id,
                "target_id": target_id,
                "query_name": left[(query_id, target_id)]["query_name"],
                "target_name": left[(query_id, target_id)]["target_name"],
                "p_adjust": left[(query_id, target_id)]["p_adjust"],
                "overlap": left[(query_id, target_id)]["overlap"],
            }
            for query_id, target_id in lost_examples
        ],
        "gained_significant_examples": [
            {
                "query_id": query_id,
                "target_id": target_id,
                "query_name": right[(query_id, target_id)]["query_name"],
                "target_name": right[(query_id, target_id)]["target_name"],
                "p_adjust": right[(query_id, target_id)]["p_adjust"],
                "overlap": right[(query_id, target_id)]["overlap"],
            }
            for query_id, target_id in gained_examples
        ],
        "top_abs_delta_pairs": deltas[:limit],
    }


def write_yaml(path: Path, value: object) -> None:
    with path.open("w") as handle:
        prepare_go_eval.write_yaml_value(handle, value)


def main() -> int:
    args = parse_args()
    if not 0 <= args.max_p_adjust <= 1:
        raise SystemExit("--max-p-adjust must be between 0 and 1")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    variants = [variant.strip() for variant in args.variants.split(",") if variant.strip()]

    started = time.perf_counter()
    commands = []
    timepoint_summaries = {}

    for timepoint in DEFAULT_TIMEPOINTS:
        label = timepoint["label"]
        release_date = timepoint["release_date"]
        out_dir = args.out_dir / label
        run_metadata = out_dir / "run_metadata.json"
        if args.skip_existing and run_metadata.exists():
            command_result = {
                "command": ["<skipped-existing>", str(out_dir)],
                "runtime_seconds": 0,
                "stdout": "",
                "stderr": "",
            }
        else:
            command = [
                sys.executable,
                "scripts/run_disease20_go_eval.py",
                "--manifest",
                str(args.manifest),
                "--out-dir",
                str(out_dir),
                "--variants",
                ",".join(variants),
                "--max-p-adjust",
                str(args.max_p_adjust),
                "--ontology-url",
                release_url(release_date, "ontology/go-basic.obo"),
                "--gaf-url",
                release_url(release_date, "annotations/goa_human.gaf.gz"),
            ]
            if args.force_download:
                command.append("--force-download")
            command_result = run_command(command)
        commands.append(command_result)

        prep = read_json(out_dir / "metadata.json")
        run = read_json(out_dir / "run_metadata.json")
        variant_summaries = {}
        for variant in variants:
            variant_summaries[variant] = summarize_variant(out_dir / variant / "results.tsv")
            variant_summaries[variant]["runtime_seconds"] = run["variants"][variant]["timing"][
                "runtime_seconds"
            ]
            variant_summaries[variant]["annotation_pairs"] = prep["annotations"]["stats"][
                "variant_unique_pairs"
            ][variant]
            variant_summaries[variant]["annotation_genes"] = prep["annotations"]["stats"][
                "variant_gene_counts"
            ][variant]

        timepoint_summaries[label] = {
            **timepoint,
            "out_dir": str(out_dir),
            "ontology": {
                "term_count": prep["ontology"]["tables"]["term_count"],
                "closure_rows": prep["ontology"]["tables"]["closure_rows"],
                "relations_for_closure": prep["ontology"]["relations_for_closure"],
                "sha256": prep["ontology"]["file"]["sha256"],
            },
            "annotations": {
                "background_gene_count": prep["background"]["gene_count"],
                "gaf_version": prep["annotations"]["gaf_header"].get("gaf_version"),
                "date_generated": prep["annotations"]["gaf_header"].get("date_generated"),
                "go_version": prep["annotations"]["gaf_header"].get("go_version"),
                "not_qualified_lines": prep["annotations"]["stats"]["not_qualified_lines"],
                "contributes_to_lines": prep["annotations"]["stats"]["contributes_to_lines"],
                "sha256": prep["annotations"]["file"]["sha256"],
            },
            "variants": variant_summaries,
        }

    diffs = {}
    labels = [timepoint["label"] for timepoint in DEFAULT_TIMEPOINTS]
    comparisons = [
        (labels[2], labels[1]),
        (labels[1], labels[0]),
        (labels[2], labels[0]),
    ]
    for variant in variants:
        diffs[variant] = {}
        for left_label, right_label in comparisons:
            left = build_pair_index(args.out_dir / left_label / variant / "results.tsv")
            right = build_pair_index(args.out_dir / right_label / variant / "results.tsv")
            diffs[variant][f"{left_label}__to__{right_label}"] = diff_pair_indexes(left, right)

    metadata = {
        "generated_at_utc": utc_now(),
        "runtime_seconds": round(time.perf_counter() - started, 3),
        "eval": "go_timepoints_disease20",
        "scope": "20 Disease Ontology query gene sets from evals/disease20/sets.tsv against GOA human GO targets.",
        "result_filter": {
            "max_p_adjust": args.max_p_adjust,
            "correction": "bonferroni",
        },
        "timepoints": timepoint_summaries,
        "diffs": diffs,
        "commands": commands,
    }
    with (args.out_dir / "summary.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_yaml(args.out_dir / "summary.yaml", metadata)

    print(f"Wrote {args.out_dir / 'summary.yaml'}")
    for label, summary in timepoint_summaries.items():
        all_variant = summary["variants"].get("all", {})
        print(
            f"{label}: {summary['ontology']['term_count']} terms, "
            f"{summary['annotations']['background_gene_count']} background genes, "
            f"all rows={all_variant.get('rows')}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
