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
from genesets_workflows.sources import mygeneset_stratified
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
        "term_coverage": {
            "enabled": True,
            "example_limit": 10,
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
    query_sets["include_regex"] = normalize_regex_list(
        query_sets.get("include_regex"),
        "query_sets.include_regex",
    )
    query_sets["exclude_regex"] = normalize_regex_list(
        query_sets.get("exclude_regex"),
        "query_sets.exclude_regex",
    )
    fetch = query_sets["fetch"]
    fetch.setdefault("mode", "query")
    fetch["limit"] = int(fetch.get("limit", query_sets["limit"]))
    fetch["min_genes"] = int(fetch.get("min_genes", 20))
    fetch["max_genes"] = int(fetch.get("max_genes", 1000))
    if fetch["mode"] == "stratified":
        if "config" not in fetch:
            raise SystemExit("query_sets.fetch.config is required when fetch.mode is stratified")
        fetch["config"] = resolve_path(base, fetch["config"])
    elif fetch["mode"] != "query":
        raise SystemExit(f"unsupported query_sets.fetch.mode: {fetch['mode']}")

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

    term_coverage = merged.setdefault("term_coverage", {})
    term_coverage["enabled"] = bool(term_coverage.get("enabled", True))
    term_coverage["example_limit"] = int(term_coverage.get("example_limit", 10))

    output = merged["output"]
    output["dir"] = resolve_path(base, output["dir"])
    output.setdefault("result_stem", f"expression{query_sets['limit']}")

    if not 0 <= statistics["max_p_adjust"] <= 1:
        raise SystemExit("statistics.max_p_adjust must be between 0 and 1")
    if not 0 <= compare["p_adjust_cutoff"] <= 1:
        raise SystemExit("compare.p_adjust_cutoff must be between 0 and 1")
    if query_sets["limit"] <= 0:
        raise SystemExit("query_sets.limit must be positive")
    if term_coverage["example_limit"] < 0:
        raise SystemExit("term_coverage.example_limit must be non-negative")
    return merged


def normalize_regex_list(value: Any, field: str) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, list) and all(isinstance(item, str) for item in value):
        return value
    raise SystemExit(f"{field} must be a string or list of strings")


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


def metadata_value(value: Any) -> Any:
    if isinstance(value, Path):
        return display_path(value)
    if isinstance(value, dict):
        return {key: metadata_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [metadata_value(item) for item in value]
    return value


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

    if fetch.get("mode") == "stratified":
        fetch_command = [
            str(fetch["config"]),
            "--out-dir",
            str(source_dir),
            "--limit",
            str(query_sets["limit"]),
        ]
        args = mygeneset_stratified.parse_args(fetch_command)
        raw_plan = mygeneset_stratified.load_plan(args.config)
        plan = mygeneset_stratified.normalize_plan(raw_plan, args.config.parent, args)
        started = time.perf_counter()
        fetched_metadata = mygeneset_stratified.fetch_stratified(plan)
        return (
            {
                "status": "fetched",
                "queries": display_path(queries),
                "metadata": display_path(metadata),
                "command": ["genesets-workflows", "fetch-mygeneset-stratified", *fetch_command],
                "runtime_seconds": round(time.perf_counter() - started, 3),
                "emitted": fetched_metadata["emitted"],
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


def compile_patterns(patterns: list[str], field: str) -> list[re.Pattern[str]]:
    compiled: list[re.Pattern[str]] = []
    for pattern in patterns:
        try:
            compiled.append(re.compile(pattern))
        except re.error as err:
            raise SystemExit(f"invalid regex in {field}: {pattern!r}: {err}") from err
    return compiled


def query_match_text(query_id: str, metadata: dict[str, Any]) -> str:
    fields = [
        query_id,
        metadata.get("name", ""),
        metadata.get("description", ""),
        metadata.get("source", ""),
        metadata.get("source_class", ""),
        metadata.get("stratum", ""),
        metadata.get("search_query", ""),
    ]
    return "\n".join(str(value) for value in fields if value is not None)


def select_queries(
    source_gmt: Path,
    source_metadata: Path,
    out_dir: Path,
    limit: int,
    include_regex: list[str] | None = None,
    exclude_regex: list[str] | None = None,
    selection: str | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    out_dir.mkdir(parents=True, exist_ok=True)
    selected_gmt = out_dir / "queries.gmt"
    selected_metadata = out_dir / "queries.metadata.json"

    include_regex = include_regex or []
    exclude_regex = exclude_regex or []
    include_patterns = compile_patterns(include_regex, "query_sets.include_regex")
    exclude_patterns = compile_patterns(exclude_regex, "query_sets.exclude_regex")
    metadata = json.loads(source_metadata.read_text())
    source_sets = {item["id"]: item for item in metadata.get("sets", [])}

    rows: list[str] = []
    ids: list[str] = []
    considered = 0
    skipped_include = 0
    skipped_exclude = 0
    with source_gmt.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            considered += 1
            set_id = line.split("\t", 1)[0]
            match_text = query_match_text(set_id, source_sets.get(set_id, {"id": set_id}))
            if include_patterns and not any(pattern.search(match_text) for pattern in include_patterns):
                skipped_include += 1
                continue
            if exclude_patterns and any(pattern.search(match_text) for pattern in exclude_patterns):
                skipped_exclude += 1
                continue
            rows.append(line)
            ids.append(set_id)
            if len(rows) >= limit:
                break
    if len(rows) < limit:
        raise SystemExit(
            f"{source_gmt} has only {len(rows)} selectable non-empty sets after filters; "
            f"requested {limit}"
        )
    selected_gmt.write_text("".join(rows))

    selected_sets = [source_sets.get(set_id, {"id": set_id}) for set_id in ids]
    selected_metadata.write_text(
        json.dumps(
            {
                "generated_at_utc": utc_now(),
                "source_gmt": display_path(source_gmt),
                "source_metadata": display_path(source_metadata),
                "selected_count": len(ids),
                "selection": selection or "first N emitted sets from the source query snapshot",
                "include_regex": include_regex,
                "exclude_regex": exclude_regex,
                "considered_count": considered,
                "skipped_include": skipped_include,
                "skipped_exclude": skipped_exclude,
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
        "selection": selection or "first N emitted sets from the source query snapshot",
        "include_regex": include_regex,
        "exclude_regex": exclude_regex,
        "considered_count": considered,
        "skipped_include": skipped_include,
        "skipped_exclude": skipped_exclude,
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


def sql_literal(value: Path | str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def write_term_coverage(
    snapshot: dict[str, Any],
    results: Path,
    output: Path,
    example_limit: int,
) -> tuple[dict[str, Any], float]:
    started = time.perf_counter()
    go_dir = snapshot["dir"]
    annotation_variant = snapshot["annotation_variant"]
    terms = go_dir / "terms.tsv"
    closure = go_dir / "closure.tsv"
    annotations = go_dir / annotation_variant / "gene_terms.tsv"
    background = go_dir / "background_all_goa_symbols.txt"
    output.parent.mkdir(parents=True, exist_ok=True)

    coverage_query = f"""
    WITH
    terms AS (
      SELECT term_id, name
      FROM read_csv_auto({sql_literal(terms)}, delim='\\t')
    ),
    background AS (
      SELECT gene_id
      FROM read_csv_auto({sql_literal(background)}, delim='\\t')
    ),
    direct_annotations AS (
      SELECT DISTINCT annotations.gene_id, annotations.term_id AS child
      FROM read_csv_auto({sql_literal(annotations)}, delim='\\t') AS annotations
      JOIN background USING (gene_id)
    ),
    propagated_targets AS (
      SELECT closure.ancestor AS term_id, count(DISTINCT direct_annotations.gene_id) AS target_size
      FROM direct_annotations
      JOIN read_csv_auto({sql_literal(closure)}, delim='\\t') AS closure
        ON direct_annotations.child = closure.child
      GROUP BY closure.ancestor
    ),
    significant_hits AS (
      SELECT
        target_id AS term_id,
        count(*) AS significant_rows,
        count(DISTINCT query_id) AS significant_queries,
        min(p_adjust_bonferroni) AS best_p_adjust,
        max(overlap) AS max_overlap
      FROM read_parquet({sql_literal(results)})
      GROUP BY target_id
    )
    SELECT
      terms.term_id,
      terms.name,
      coalesce(propagated_targets.target_size, 0)::UBIGINT AS target_size,
      coalesce(significant_hits.significant_rows, 0)::UBIGINT AS significant_rows,
      coalesce(significant_hits.significant_queries, 0)::UBIGINT AS significant_queries,
      significant_hits.best_p_adjust,
      coalesce(significant_hits.max_overlap, 0)::UBIGINT AS max_overlap,
      CASE
        WHEN coalesce(propagated_targets.target_size, 0) = 0 THEN 'unscorable'
        WHEN coalesce(significant_hits.significant_rows, 0) = 0 THEN 'scorable_never_significant'
        ELSE 'significant'
      END AS coverage_status
    FROM terms
    LEFT JOIN propagated_targets USING (term_id)
    LEFT JOIN significant_hits USING (term_id)
    """
    copy_sql = f"COPY ({coverage_query}) TO {sql_literal(output)} (FORMAT PARQUET)"

    try:
        import duckdb  # type: ignore
    except ImportError:
        run_command(["duckdb", "-c", copy_sql])
        summary_rows = duckdb_json(term_coverage_summary_sql(output))
        example_rows = duckdb_json(term_coverage_examples_sql(output, example_limit))
        if summary_rows is None or example_rows is None:
            raise SystemExit("go-impact term coverage requires Python duckdb or the duckdb CLI")
        summary = summary_rows[0]
        examples = example_rows
    else:
        con = duckdb.connect()
        con.execute(copy_sql)
        summary = con.execute(term_coverage_summary_sql(output)).fetchone()
        summary = {
            "total_terms": int(summary[0]),
            "scorable_terms": int(summary[1]),
            "significant_terms": int(summary[2]),
            "scorable_never_significant_terms": int(summary[3]),
            "unscorable_terms": int(summary[4]),
        }
        examples = [
            {
                "term_id": row[0],
                "name": row[1],
                "target_size": int(row[2]),
            }
            for row in con.execute(term_coverage_examples_sql(output, example_limit)).fetchall()
        ]

    if "total_terms" not in summary:
        summary = {
            "total_terms": int(summary["total_terms"]),
            "scorable_terms": int(summary["scorable_terms"]),
            "significant_terms": int(summary["significant_terms"]),
            "scorable_never_significant_terms": int(summary["scorable_never_significant_terms"]),
            "unscorable_terms": int(summary["unscorable_terms"]),
        }
    return (
        {
            **summary,
            "examples_largest_scorable_never_significant": examples,
            "path": display_path(output),
        },
        round(time.perf_counter() - started, 3),
    )


def term_coverage_summary_sql(path: Path) -> str:
    return f"""
    SELECT
      count(*) AS total_terms,
      count(*) FILTER (WHERE coverage_status != 'unscorable') AS scorable_terms,
      count(*) FILTER (WHERE coverage_status = 'significant') AS significant_terms,
      count(*) FILTER (WHERE coverage_status = 'scorable_never_significant') AS scorable_never_significant_terms,
      count(*) FILTER (WHERE coverage_status = 'unscorable') AS unscorable_terms
    FROM read_parquet({sql_literal(path)})
    """


def term_coverage_examples_sql(path: Path, limit: int) -> str:
    return f"""
    SELECT term_id, name, target_size
    FROM read_parquet({sql_literal(path)})
    WHERE coverage_status = 'scorable_never_significant'
    ORDER BY target_size DESC, term_id
    LIMIT {limit}
    """


def write_term_coverage_comparison(
    left_coverage: Path,
    right_coverage: Path,
    output: Path,
) -> tuple[dict[str, int], float]:
    started = time.perf_counter()
    output.parent.mkdir(parents=True, exist_ok=True)
    comparison_query = f"""
    SELECT
      coalesce(left_cov.term_id, right_cov.term_id) AS term_id,
      coalesce(left_cov.name, right_cov.name) AS name,
      left_cov.coverage_status AS left_status,
      right_cov.coverage_status AS right_status,
      left_cov.target_size AS left_target_size,
      right_cov.target_size AS right_target_size,
      left_cov.significant_rows AS left_significant_rows,
      right_cov.significant_rows AS right_significant_rows,
      CASE
        WHEN left_cov.coverage_status = 'significant' OR right_cov.coverage_status = 'significant'
          THEN 'significant_either'
        WHEN left_cov.coverage_status = 'scorable_never_significant'
          AND right_cov.coverage_status = 'scorable_never_significant'
          THEN 'scorable_never_significant_both'
        WHEN left_cov.coverage_status = 'scorable_never_significant'
          OR right_cov.coverage_status = 'scorable_never_significant'
          THEN 'scorable_never_significant_one_side'
        ELSE 'unscorable_both_or_absent'
      END AS comparison_status
    FROM read_parquet({sql_literal(left_coverage)}) AS left_cov
    FULL OUTER JOIN read_parquet({sql_literal(right_coverage)}) AS right_cov
      USING (term_id)
    """
    copy_sql = f"COPY ({comparison_query}) TO {sql_literal(output)} (FORMAT PARQUET)"

    try:
        import duckdb  # type: ignore
    except ImportError:
        run_command(["duckdb", "-c", copy_sql])
        rows = duckdb_json(term_coverage_comparison_summary_sql(output))
        if rows is None:
            raise SystemExit("go-impact term coverage comparison requires Python duckdb or the duckdb CLI")
        summary = {str(row["comparison_status"]): int(row["terms"]) for row in rows}
    else:
        con = duckdb.connect()
        con.execute(copy_sql)
        rows = con.execute(term_coverage_comparison_summary_sql(output)).fetchall()
        summary = {str(status): int(count) for status, count in rows}
    return (
        {
            **summary,
            "path": display_path(output),
        },
        round(time.perf_counter() - started, 3),
    )


def term_coverage_comparison_summary_sql(path: Path) -> str:
    return f"""
    SELECT comparison_status, count(*) AS terms
    FROM read_parquet({sql_literal(path)})
    GROUP BY comparison_status
    """


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
        config["query_sets"]["include_regex"],
        config["query_sets"]["exclude_regex"],
        config["query_sets"].get("selection"),
    )
    queries = out_dir / "queries.gmt"

    result_stem = config["output"]["result_stem"]
    left_prefix = left_snapshot["output_prefix"]
    right_prefix = right_snapshot["output_prefix"]
    old_results = out_dir / f"{left_prefix}_{result_stem}.parquet"
    new_results = out_dir / f"{right_prefix}_{result_stem}.parquet"
    diff_results = out_dir / f"{left_prefix}_vs_{right_prefix}_{result_stem}.diff.parquet"
    diff_metadata = out_dir / f"{left_prefix}_vs_{right_prefix}_{result_stem}.diff.yaml"
    left_term_coverage = out_dir / f"{left_prefix}_{result_stem}.term_coverage.parquet"
    right_term_coverage = out_dir / f"{right_prefix}_{result_stem}.term_coverage.parquet"
    term_coverage_comparison = out_dir / (
        f"{left_prefix}_vs_{right_prefix}_{result_stem}.term_coverage.parquet"
    )

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
    term_coverage_summary: dict[str, Any] | None = None
    term_coverage_timings: dict[str, float] = {}
    if config["term_coverage"]["enabled"]:
        old_term_coverage, old_term_coverage_seconds = write_term_coverage(
            left_snapshot,
            old_results,
            left_term_coverage,
            config["term_coverage"]["example_limit"],
        )
        new_term_coverage, new_term_coverage_seconds = write_term_coverage(
            right_snapshot,
            new_results,
            right_term_coverage,
            config["term_coverage"]["example_limit"],
        )
        coverage_comparison, coverage_comparison_seconds = write_term_coverage_comparison(
            left_term_coverage,
            right_term_coverage,
            term_coverage_comparison,
        )
        term_coverage_summary = {
            "old": old_term_coverage,
            "new": new_term_coverage,
            "comparison": coverage_comparison,
        }
        term_coverage_timings = {
            "old_term_coverage_seconds": old_term_coverage_seconds,
            "new_term_coverage_seconds": new_term_coverage_seconds,
            "term_coverage_comparison_seconds": coverage_comparison_seconds,
        }

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
                "include_regex": config["query_sets"]["include_regex"],
                "exclude_regex": config["query_sets"]["exclude_regex"],
                "fetch": metadata_value(config["query_sets"]["fetch"]),
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
            "term_coverage": config["term_coverage"],
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
            "old_term_coverage": display_path(left_term_coverage) if term_coverage_summary else None,
            "new_term_coverage": display_path(right_term_coverage) if term_coverage_summary else None,
            "term_coverage_comparison": (
                display_path(term_coverage_comparison) if term_coverage_summary else None
            ),
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
        "term_coverage": term_coverage_summary,
        "timings": {
            "select_queries_seconds": selected["runtime_seconds"],
            "old_matrix_seconds": old_run["runtime_seconds"],
            "new_matrix_seconds": new_run["runtime_seconds"],
            "compare_seconds": compare_run["runtime_seconds"],
            **term_coverage_timings,
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
    if report.get("term_coverage"):
        coverage = report["term_coverage"]
        comparison = coverage["comparison"]
        print(
            "Term coverage: "
            f"old_never={coverage['old']['scorable_never_significant_terms']}, "
            f"new_never={coverage['new']['scorable_never_significant_terms']}, "
            f"never_both={comparison.get('scorable_never_significant_both', 0)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
