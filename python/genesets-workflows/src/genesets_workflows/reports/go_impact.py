from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Any

import yaml

from genesets_workflows.runtime import display_path, duckdb_json, resolve_path, run_command, utc_now
from genesets_workflows.sources import mygeneset
from genesets_workflows.yaml_io import write_yaml_value


DEFAULT_SOURCE_DIR = Path("evals/expression_like/generated/msigdb_gse_500")
DEFAULT_OLD_GO = Path("evals/go_timepoints/generated/five_years_ago_2021_05_01")
DEFAULT_NEW_GO = Path("evals/go_timepoints/generated/now_2026_03_25")
DEFAULT_LIMIT = 500


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="genesets-workflows go-impact",
        description="Build a GO impact report over expression-derived gene sets."
    )
    parser.add_argument(
        "config",
        nargs="?",
        type=Path,
        help="YAML report configuration. Relative paths resolve from the config file.",
    )
    parser.add_argument("--source-dir", default=None, type=Path)
    parser.add_argument("--out-dir", default=None, type=Path)
    parser.add_argument("--old-go-dir", default=None, type=Path)
    parser.add_argument("--new-go-dir", default=None, type=Path)
    parser.add_argument("--limit", default=None, type=int)
    parser.add_argument("--binary", default=None)
    parser.add_argument("--min-overlap", default=None, type=int)
    parser.add_argument("--max-p-adjust", default=None, type=float)
    parser.add_argument("--fetch-query", default=None)
    parser.add_argument("--fetch-limit", default=None, type=int)
    parser.add_argument("--source-filter", default=None)
    return parser.parse_args(argv)


def default_config() -> dict[str, Any]:
    return {
        "analysis": "go_impact_5y_expression500",
        "question": (
            "Across 500 expression-derived MSigDB/GSE gene sets, which enrichment "
            "calls cross a Bonferroni-adjusted p-value threshold between the 2021 "
            "and 2026 GO/GOA snapshots?"
        ),
        "binary": "genesets-rs",
        "query_sets": {
            "source_dir": str(DEFAULT_SOURCE_DIR),
            "limit": DEFAULT_LIMIT,
            "selection": "first N emitted sets from the source query snapshot",
            "fetch": {
                "query": "GSE*",
                "source_filter": "msigdb",
                "limit": DEFAULT_LIMIT,
                "min_genes": 20,
                "max_genes": 1000,
            },
        },
        "snapshots": {
            "left": {
                "label": "2021-05-01",
                "dir": str(DEFAULT_OLD_GO),
                "annotation_variant": "all",
                "output_prefix": "go_2021_05_01",
            },
            "right": {
                "label": "2026-03-25",
                "dir": str(DEFAULT_NEW_GO),
                "annotation_variant": "all",
                "output_prefix": "go_2026_03_25",
            },
        },
        "statistics": {
            "min_overlap": 2,
            "max_p_adjust": 0.05,
            "correction": "bonferroni",
        },
        "compare": {
            "p_adjust_cutoff": 0.05,
            "crossings_only": False,
        },
        "post_processing": {
            "ranked_term_scope": {
                "subset": "goslim_generic",
                "snapshot": "right",
                "include_descendants": True,
                "max_target_size": 1000,
                "mode": "include",
            }
        },
        "output": {
            "dir": "notebooks/generated/go_impact_5y_expression500",
            "result_stem": "expression500",
        },
    }


def load_config(args: argparse.Namespace) -> tuple[dict[str, Any], Path | None]:
    if args.config is None:
        config = default_config()
        base = Path(".")
    else:
        with args.config.open() as handle:
            config = yaml.safe_load(handle) or {}
        if not isinstance(config, dict):
            raise SystemExit(f"{args.config} must contain a YAML object")
        base = args.config.parent if args.config.parent != Path("") else Path(".")

    apply_cli_overrides(config, args)
    return normalize_config(config, base), args.config


def apply_cli_overrides(config: dict[str, Any], args: argparse.Namespace) -> None:
    if args.binary is not None:
        config["binary"] = args.binary
    if args.source_dir is not None:
        config.setdefault("query_sets", {})["source_dir"] = str(args.source_dir)
    if args.limit is not None:
        config.setdefault("query_sets", {})["limit"] = args.limit
        config.setdefault("output", {})["result_stem"] = f"expression{args.limit}"
    if args.out_dir is not None:
        config.setdefault("output", {})["dir"] = str(args.out_dir)
    if args.old_go_dir is not None:
        config.setdefault("snapshots", {}).setdefault("left", {})["dir"] = str(args.old_go_dir)
    if args.new_go_dir is not None:
        config.setdefault("snapshots", {}).setdefault("right", {})["dir"] = str(args.new_go_dir)
    if args.min_overlap is not None:
        config.setdefault("statistics", {})["min_overlap"] = args.min_overlap
    if args.max_p_adjust is not None:
        config.setdefault("statistics", {})["max_p_adjust"] = args.max_p_adjust
        config.setdefault("compare", {})["p_adjust_cutoff"] = args.max_p_adjust
    fetch = config.setdefault("query_sets", {}).setdefault("fetch", {})
    if args.fetch_query is not None:
        fetch["query"] = args.fetch_query
    if args.fetch_limit is not None:
        fetch["limit"] = args.fetch_limit
    if args.source_filter is not None:
        fetch["source_filter"] = args.source_filter


def normalize_config(config: dict[str, Any], base: Path) -> dict[str, Any]:
    merged = merge_defaults(config, default_config())

    query_sets = merged["query_sets"]
    query_sets["source_dir"] = resolve_path(base, query_sets["source_dir"])
    query_sets["limit"] = int(query_sets["limit"])
    query_sets["fetch"]["limit"] = int(query_sets["fetch"].get("limit", query_sets["limit"]))
    query_sets["fetch"]["min_genes"] = int(query_sets["fetch"].get("min_genes", 20))
    query_sets["fetch"]["max_genes"] = int(query_sets["fetch"].get("max_genes", 1000))

    for side in ("left", "right"):
        snapshot = merged["snapshots"][side]
        snapshot["dir"] = resolve_path(base, snapshot["dir"])
        snapshot.setdefault("annotation_variant", "all")
        snapshot.setdefault("output_prefix", f"go_{slugify(snapshot['label'])}")

    statistics = merged["statistics"]
    statistics["min_overlap"] = int(statistics.get("min_overlap", 2))
    statistics["max_p_adjust"] = float(statistics.get("max_p_adjust", 0.05))
    if statistics.get("correction", "bonferroni") != "bonferroni":
        raise SystemExit("go-impact currently supports correction: bonferroni")

    compare = merged.setdefault("compare", {})
    compare["p_adjust_cutoff"] = float(compare.get("p_adjust_cutoff", statistics["max_p_adjust"]))
    compare["crossings_only"] = bool(compare.get("crossings_only", False))

    output = merged["output"]
    output["dir"] = resolve_path(base, output["dir"])
    output.setdefault("result_stem", f"expression{query_sets['limit']}")

    if not 0 <= statistics["max_p_adjust"] <= 1:
        raise SystemExit("statistics.max_p_adjust must be between 0 and 1")
    if not 0 <= compare["p_adjust_cutoff"] <= 1:
        raise SystemExit("compare.p_adjust_cutoff must be between 0 and 1")
    if query_sets["limit"] <= 0:
        raise SystemExit("query_sets.limit must be positive")
    return merged


def merge_defaults(config: dict[str, Any], defaults: dict[str, Any]) -> dict[str, Any]:
    merged = json.loads(json.dumps(defaults))
    deep_update(merged, config)
    return merged


def deep_update(target: dict[str, Any], updates: dict[str, Any]) -> None:
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            deep_update(target[key], value)
        else:
            target[key] = value


def slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()


def ensure_source_queries(config: dict[str, Any]) -> tuple[dict[str, Any], Path, Path]:
    query_sets = config["query_sets"]
    fetch = query_sets["fetch"]
    source_dir = query_sets["source_dir"]
    queries = source_dir / "queries.gmt"
    metadata = source_dir / "metadata.json"
    if queries.exists() and metadata.exists():
        return (
            {
                "status": "reused_existing",
                "queries": display_path(queries),
                "metadata": display_path(metadata),
                "command": None,
                "runtime_seconds": 0,
            },
            queries,
            metadata,
        )

    fetch_command = [
        "--query",
        str(fetch["query"]),
        "--limit",
        str(fetch["limit"]),
        "--out-dir",
        str(source_dir),
        "--min-genes",
        str(fetch["min_genes"]),
        "--max-genes",
        str(fetch["max_genes"]),
    ]
    if fetch.get("source_filter") is not None:
        fetch_command.extend(["--source-filter", str(fetch["source_filter"])])
    args = mygeneset.parse_args(fetch_command)
    started = time.perf_counter()
    fetched_metadata = mygeneset.fetch_query(args)
    return (
        {
            "status": "fetched",
            "queries": display_path(queries),
            "metadata": display_path(metadata),
            "command": ["genesets-workflows", "fetch-mygeneset", *fetch_command],
            "runtime_seconds": round(time.perf_counter() - started, 3),
            "emitted": fetched_metadata["emitted"],
        },
        queries,
        metadata,
    )


def select_queries(source_gmt: Path, source_metadata: Path, out_dir: Path, limit: int) -> dict[str, Any]:
    started = time.perf_counter()
    out_dir.mkdir(parents=True, exist_ok=True)
    selected_gmt = out_dir / "queries.gmt"
    selected_metadata = out_dir / "queries.metadata.json"

    rows: list[str] = []
    ids: list[str] = []
    with source_gmt.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            rows.append(line)
            ids.append(line.split("\t", 1)[0])
            if len(rows) >= limit:
                break
    if len(rows) < limit:
        raise SystemExit(f"{source_gmt} has only {len(rows)} non-empty sets; requested {limit}")
    selected_gmt.write_text("".join(rows))

    metadata = json.loads(source_metadata.read_text())
    source_sets = {item["id"]: item for item in metadata.get("sets", [])}
    selected_sets = [source_sets.get(set_id, {"id": set_id}) for set_id in ids]
    selected_metadata.write_text(
        json.dumps(
            {
                "generated_at_utc": utc_now(),
                "source_gmt": display_path(source_gmt),
                "source_metadata": display_path(source_metadata),
                "selected_count": len(ids),
                "selection": "first N emitted sets from the source query snapshot",
                "sets": selected_sets,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return {
        "runtime_seconds": round(time.perf_counter() - started, 3),
        "queries": display_path(selected_gmt),
        "metadata": display_path(selected_metadata),
        "selected_count": len(ids),
    }


def require_go_tables(go_dir: Path, annotation_variant: str) -> None:
    required = [
        go_dir / annotation_variant / "gene_terms.tsv",
        go_dir / "terms.tsv",
        go_dir / "closure.tsv",
        go_dir / "background_all_goa_symbols.txt",
    ]
    missing = [path for path in required if not path.exists()]
    if missing:
        raise SystemExit(f"Missing prepared GO table(s): {', '.join(str(path) for path in missing)}")


def matrix_command(
    config: dict[str, Any],
    snapshot: dict[str, Any],
    queries: Path,
    output: Path,
) -> list[str]:
    statistics = config["statistics"]
    go_dir = snapshot["dir"]
    annotation_variant = snapshot["annotation_variant"]
    return [
        config["binary"],
        "matrix",
        "--annotations",
        str(go_dir / annotation_variant / "gene_terms.tsv"),
        "--terms",
        str(go_dir / "terms.tsv"),
        "--closure",
        str(go_dir / "closure.tsv"),
        "--queries",
        str(queries),
        "--query-format",
        "gmt",
        "--background",
        str(go_dir / "background_all_goa_symbols.txt"),
        "--min-overlap",
        str(statistics["min_overlap"]),
        "--max-p-adjust",
        str(statistics["max_p_adjust"]),
        "--output-format",
        "parquet",
        "--output",
        str(output),
    ]


def summarize_with_duckdb(path: Path) -> tuple[int, int, int]:
    try:
        import duckdb  # type: ignore
    except ImportError:
        rows = duckdb_json(
            f"""
            SELECT
              count(*) AS rows,
              count(DISTINCT query_id) AS queries,
              count(DISTINCT target_id) AS terms
            FROM read_parquet('{path}')
            """
        )
        if not rows:
            raise SystemExit("go-impact requires Python duckdb or the duckdb CLI for summaries")
        row = rows[0]
        return int(row["rows"]), int(row["queries"]), int(row["terms"])

    con = duckdb.connect()
    row = con.execute(
        """
        SELECT
          count(*) AS rows,
          count(DISTINCT query_id) AS queries,
          count(DISTINCT target_id) AS terms
        FROM read_parquet(?)
        """,
        [str(path)],
    ).fetchone()
    return int(row[0]), int(row[1]), int(row[2])


def diff_counts_with_duckdb(path: Path) -> dict[str, int]:
    try:
        import duckdb  # type: ignore
    except ImportError:
        rows = duckdb_json(
            f"""
            SELECT class, count(*) AS rows
            FROM read_parquet('{path}')
            GROUP BY class
            """
        )
        if rows is None:
            raise SystemExit("go-impact requires Python duckdb or the duckdb CLI for summaries")
        return {str(row["class"]): int(row["rows"]) for row in rows}

    con = duckdb.connect()
    rows = con.execute(
        """
        SELECT class, count(*) AS rows
        FROM read_parquet(?)
        GROUP BY class
        """,
        [str(path)],
    ).fetchall()
    return {class_name: int(count) for class_name, count in rows}


def run_report(args: argparse.Namespace) -> dict[str, Any]:
    config, config_path = load_config(args)

    started = time.perf_counter()
    out_dir = config["output"]["dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    left_snapshot = config["snapshots"]["left"]
    right_snapshot = config["snapshots"]["right"]
    require_go_tables(left_snapshot["dir"], left_snapshot["annotation_variant"])
    require_go_tables(right_snapshot["dir"], right_snapshot["annotation_variant"])

    source, source_gmt, source_metadata = ensure_source_queries(config)
    selected = select_queries(
        source_gmt,
        source_metadata,
        out_dir,
        config["query_sets"]["limit"],
    )
    queries = out_dir / "queries.gmt"

    result_stem = config["output"]["result_stem"]
    left_prefix = left_snapshot["output_prefix"]
    right_prefix = right_snapshot["output_prefix"]
    old_results = out_dir / f"{left_prefix}_{result_stem}.parquet"
    new_results = out_dir / f"{right_prefix}_{result_stem}.parquet"
    diff_results = out_dir / f"{left_prefix}_vs_{right_prefix}_{result_stem}.diff.parquet"
    diff_metadata = out_dir / f"{left_prefix}_vs_{right_prefix}_{result_stem}.diff.yaml"

    old_run = run_command(matrix_command(config, left_snapshot, queries, old_results))
    new_run = run_command(matrix_command(config, right_snapshot, queries, new_results))
    compare_command = [
        config["binary"],
        "compare",
        "--left",
        str(old_results),
        "--right",
        str(new_results),
        "--p-adjust-cutoff",
        str(config["compare"]["p_adjust_cutoff"]),
        "--output-format",
        "parquet",
        "--output",
        str(diff_results),
        "--metadata-output",
        str(diff_metadata),
    ]
    if config["compare"]["crossings_only"]:
        compare_command.append("--crossings-only")
    compare_run = run_command(compare_command)

    old_rows, old_queries, old_terms = summarize_with_duckdb(old_results)
    new_rows, new_queries, new_terms = summarize_with_duckdb(new_results)
    diff_counts = diff_counts_with_duckdb(diff_results)

    report = {
        "generated_at_utc": utc_now(),
        "analysis": config["analysis"],
        "question": config["question"],
        "config_path": display_path(config_path) if config_path else None,
        "report_parameters": {
            "query_sets": {
                "source_dir": display_path(config["query_sets"]["source_dir"]),
                "limit": config["query_sets"]["limit"],
                "selection": config["query_sets"].get("selection"),
                "fetch": config["query_sets"]["fetch"],
            },
            "snapshots": {
                "left": {
                    **left_snapshot,
                    "dir": display_path(left_snapshot["dir"]),
                },
                "right": {
                    **right_snapshot,
                    "dir": display_path(right_snapshot["dir"]),
                },
            },
            "statistics": config["statistics"],
            "compare": config["compare"],
            "post_processing": config.get("post_processing", {}),
        },
        "query_sets": {
            **selected,
            "source": source,
        },
        "go_snapshots": {
            "old": {
                "label": left_snapshot["label"],
                "path": display_path(left_snapshot["dir"]),
                "annotation_variant": left_snapshot["annotation_variant"],
            },
            "new": {
                "label": right_snapshot["label"],
                "path": display_path(right_snapshot["dir"]),
                "annotation_variant": right_snapshot["annotation_variant"],
            },
        },
        "result_filter": {
            "min_overlap": config["statistics"]["min_overlap"],
            "max_p_adjust": config["statistics"]["max_p_adjust"],
            "correction": config["statistics"]["correction"],
        },
        "outputs": {
            "old_results": display_path(old_results),
            "new_results": display_path(new_results),
            "diff_results": display_path(diff_results),
            "diff_metadata": display_path(diff_metadata),
            "summary_json": display_path(out_dir / "summary.json"),
            "summary_yaml": display_path(out_dir / "summary.yaml"),
        },
        "result_counts": {
            "old": {
                "rows": old_rows,
                "queries_with_hits": old_queries,
                "terms": old_terms,
            },
            "new": {
                "rows": new_rows,
                "queries_with_hits": new_queries,
                "terms": new_terms,
            },
            "diff": diff_counts,
        },
        "timings": {
            "select_queries_seconds": selected["runtime_seconds"],
            "old_matrix_seconds": old_run["runtime_seconds"],
            "new_matrix_seconds": new_run["runtime_seconds"],
            "compare_seconds": compare_run["runtime_seconds"],
            "total_seconds": round(time.perf_counter() - started, 3),
        },
        "commands": {
            "old_matrix": old_run["command"],
            "new_matrix": new_run["command"],
            "compare": compare_run["command"],
        },
    }

    (out_dir / "summary.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    with (out_dir / "summary.yaml").open("w") as handle:
        write_yaml_value(handle, report)
    return report


def main(argv: list[str] | None = None) -> int:
    report = run_report(parse_args(argv))
    print(f"Wrote {report['outputs']['summary_yaml']}")
    print(
        "Diff rows: "
        + ", ".join(f"{key}={value}" for key, value in sorted(report["result_counts"]["diff"].items()))
    )
    print(
        "Timings: "
        f"old={report['timings']['old_matrix_seconds']}s, "
        f"new={report['timings']['new_matrix_seconds']}s, "
        f"compare={report['timings']['compare_seconds']}s, "
        f"total={report['timings']['total_seconds']}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
