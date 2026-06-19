from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from genesets_workflows.runtime import repo_root


DEFAULT_LIMIT = 100
MAX_LIMIT = 1000
SPECIFIC_TARGET_SIZE = 1000


def json_value(value: Any) -> Any:
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def row_dicts(cursor: Any) -> list[dict[str, Any]]:
    columns = [item[0] for item in cursor.description]
    return [
        {column: json_value(value) for column, value in zip(columns, row)}
        for row in cursor.fetchall()
    ]


def slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "-", value.strip()).strip("-").lower()
    return slug or "bundle"


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a YAML mapping")
    return data


def resolve_artifact(bundle_dir: Path, root: Path, value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if path.is_absolute():
        return path if path.exists() else path
    candidates = [
        root / path,
        bundle_dir / path,
        bundle_dir / path.name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return (root / path).resolve()


def discover_bundle_paths(paths: list[Path]) -> list[Path]:
    if paths:
        candidates = paths
    else:
        cwd = Path.cwd()
        if (cwd / "summary.yaml").exists():
            candidates = [cwd]
        else:
            candidates = sorted(Path("notebooks/generated").glob("*/summary.yaml"))

    bundle_dirs: list[Path] = []
    for candidate in candidates:
        path = candidate.resolve()
        if path.is_file() and path.name == "summary.yaml":
            bundle_dirs.append(path.parent)
        elif path.is_dir() and (path / "summary.yaml").exists():
            bundle_dirs.append(path)
        else:
            raise FileNotFoundError(f"{candidate} is not a report directory or summary.yaml")
    return bundle_dirs


def load_query_genes(path: Path | None) -> dict[str, list[str]]:
    if path is None or not path.exists():
        return {}
    genes: dict[str, list[str]] = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) >= 3:
                genes[row[0]] = [gene for gene in row[2:] if gene]
    return genes


def load_query_metadata(path: Path | None) -> list[dict[str, Any]]:
    if path is None or not path.exists():
        return []
    data = json.loads(path.read_text())
    sets = data.get("sets", [])
    if not isinstance(sets, list):
        return []
    return [item for item in sets if isinstance(item, dict) and item.get("id")]


def require_duckdb() -> Any:
    try:
        import duckdb  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "The explorer requires Python duckdb. Install with "
            "`python3 -m pip install -e python/genesets-workflows[explorer]` "
            "or `uv run --project python/genesets-workflows --extra explorer ...`."
        ) from exc
    return duckdb


@dataclass
class ReportBundle:
    path: Path
    ordinal: int
    summary: dict[str, Any] = field(init=False)
    id: str = field(init=False)
    root: Path = field(init=False)
    outputs: dict[str, Path] = field(init=False)
    query_sets: dict[str, dict[str, Any]] = field(init=False)
    query_genes: dict[str, list[str]] = field(init=False)
    rows: list[dict[str, Any]] = field(init=False)

    def __post_init__(self) -> None:
        self.path = self.path.resolve()
        self.root = repo_root(self.path)
        self.summary = read_yaml(self.path / "summary.yaml")
        analysis = str(self.summary.get("analysis") or self.path.name)
        self.id = slugify(f"{self.ordinal + 1}-{analysis}")
        self.outputs = self._resolve_outputs()

        query_meta_path = resolve_artifact(
            self.path,
            self.root,
            self.summary.get("query_sets", {}).get("metadata"),
        )
        queries_path = resolve_artifact(
            self.path,
            self.root,
            self.summary.get("query_sets", {}).get("queries"),
        )
        metadata_rows = load_query_metadata(query_meta_path)
        self.query_sets = {str(item["id"]): item for item in metadata_rows}
        self.query_genes = load_query_genes(queries_path)
        self.rows = self._build_query_rows()

    def _resolve_outputs(self) -> dict[str, Path]:
        outputs: dict[str, Path] = {}
        raw_outputs = self.summary.get("outputs", {})
        if not isinstance(raw_outputs, dict):
            return outputs
        for key, value in raw_outputs.items():
            if not isinstance(value, str):
                continue
            path = resolve_artifact(self.path, self.root, value)
            if path is not None:
                outputs[key] = path
        return outputs

    def overview(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "analysis": self.summary.get("analysis") or self.path.name,
            "question": self.summary.get("question"),
            "path": str(self.path),
            "selected_count": self.summary.get("query_sets", {}).get("selected_count"),
            "result_counts": self.summary.get("result_counts", {}),
            "timings": self.summary.get("timings", {}),
            "snapshots": self.summary.get("go_snapshots", {}),
            "parameters": self.summary.get("report_parameters", {}),
            "outputs": {key: str(path) for key, path in self.outputs.items()},
        }

    def facets(self) -> dict[str, Any]:
        source_counts: dict[str, int] = {}
        stratum_counts: dict[str, int] = {}
        for row in self.rows:
            source = str(row.get("source_class") or "unknown")
            stratum = str(row.get("stratum") or "unknown")
            source_counts[source] = source_counts.get(source, 0) + 1
            stratum_counts[stratum] = stratum_counts.get(stratum, 0) + 1
        return {
            "source_class": sorted(
                [{"value": key, "count": value} for key, value in source_counts.items()],
                key=lambda item: (-item["count"], item["value"]),
            ),
            "stratum": sorted(
                [{"value": key, "count": value} for key, value in stratum_counts.items()],
                key=lambda item: (-item["count"], item["value"]),
            ),
        }

    def _build_query_rows(self) -> list[dict[str, Any]]:
        rows: dict[str, dict[str, Any]] = {}
        for query_id, item in self.query_sets.items():
            rows[query_id] = {
                "id": query_id,
                "name": item.get("name") or query_id,
                "description": item.get("description"),
                "source": item.get("source"),
                "source_class": item.get("source_class") or item.get("source") or "unknown",
                "stratum": item.get("stratum") or "unknown",
                "search_query": item.get("search_query"),
                "gene_count": item.get("gene_count") or len(self.query_genes.get(query_id, [])),
                "old_result_count": 0,
                "new_result_count": 0,
                "old_best_p_adjust": None,
                "new_best_p_adjust": None,
                "lost_count": 0,
                "gained_count": 0,
                "shared_count": 0,
                "lost_specific_count": 0,
                "gained_specific_count": 0,
                "largest_abs_delta": 0.0,
                "top_changed_term": None,
                "top_crossing_term": None,
            }

        for query_id, genes in self.query_genes.items():
            rows.setdefault(
                query_id,
                {
                    "id": query_id,
                    "name": query_id,
                    "description": None,
                    "source": None,
                    "source_class": "unknown",
                    "stratum": "unknown",
                    "search_query": None,
                    "gene_count": len(genes),
                    "old_result_count": 0,
                    "new_result_count": 0,
                    "old_best_p_adjust": None,
                    "new_best_p_adjust": None,
                    "lost_count": 0,
                    "gained_count": 0,
                    "shared_count": 0,
                    "lost_specific_count": 0,
                    "gained_specific_count": 0,
                    "largest_abs_delta": 0.0,
                    "top_changed_term": None,
                    "top_crossing_term": None,
                },
            )

        self._merge_result_metrics(rows, "old_results", "old")
        self._merge_result_metrics(rows, "new_results", "new")
        self._merge_diff_metrics(rows)
        return sorted(rows.values(), key=lambda row: str(row["id"]))

    def _merge_result_metrics(self, rows: dict[str, dict[str, Any]], output_key: str, prefix: str) -> None:
        path = self.outputs.get(output_key)
        if path is None or not path.exists():
            return
        duckdb = require_duckdb()
        con = duckdb.connect()
        result = con.execute(
            """
            SELECT
              query_id,
              count(*) AS result_count,
              min(p_adjust_bonferroni) AS best_p_adjust,
              max(overlap) AS max_overlap
            FROM read_parquet(?)
            GROUP BY query_id
            """,
            [str(path)],
        )
        for item in row_dicts(result):
            row = rows.setdefault(str(item["query_id"]), {"id": str(item["query_id"])})
            row[f"{prefix}_result_count"] = item["result_count"]
            row[f"{prefix}_best_p_adjust"] = item["best_p_adjust"]
            row[f"{prefix}_max_overlap"] = item["max_overlap"]

    def _merge_diff_metrics(self, rows: dict[str, dict[str, Any]]) -> None:
        path = self.outputs.get("diff_results")
        if path is None or not path.exists():
            return
        duckdb = require_duckdb()
        con = duckdb.connect()
        metrics = con.execute(
            f"""
            SELECT
              query_id,
              sum(CASE WHEN class = 'lost_significant' THEN 1 ELSE 0 END) AS lost_count,
              sum(CASE WHEN class = 'gained_significant' THEN 1 ELSE 0 END) AS gained_count,
              sum(CASE WHEN class = 'shared_significant' THEN 1 ELSE 0 END) AS shared_count,
              sum(CASE WHEN class = 'lost_significant'
                         AND coalesce(left_target_size, 0) <= {SPECIFIC_TARGET_SIZE}
                       THEN 1 ELSE 0 END) AS lost_specific_count,
              sum(CASE WHEN class = 'gained_significant'
                         AND coalesce(right_target_size, 0) <= {SPECIFIC_TARGET_SIZE}
                       THEN 1 ELSE 0 END) AS gained_specific_count,
              max(abs(coalesce(delta_neg_log10_p_adjust, 0))) AS largest_abs_delta
            FROM read_parquet(?)
            GROUP BY query_id
            """,
            [str(path)],
        )
        for item in row_dicts(metrics):
            row = rows.setdefault(str(item["query_id"]), {"id": str(item["query_id"])})
            for key in (
                "lost_count",
                "gained_count",
                "shared_count",
                "lost_specific_count",
                "gained_specific_count",
                "largest_abs_delta",
            ):
                row[key] = item.get(key) or 0

        changed = con.execute(
            """
            SELECT
              query_id,
              class,
              target_id,
              target_name,
              delta_neg_log10_p_adjust,
              left_p_adjust,
              right_p_adjust,
              left_target_size,
              right_target_size
            FROM (
              SELECT
                *,
                row_number() OVER (
                  PARTITION BY query_id
                  ORDER BY abs(coalesce(delta_neg_log10_p_adjust, 0)) DESC,
                           coalesce(left_p_adjust, right_p_adjust, 1.0) ASC
                ) AS rank
              FROM read_parquet(?)
            )
            WHERE rank = 1
            """,
            [str(path)],
        )
        for item in row_dicts(changed):
            row = rows.setdefault(str(item["query_id"]), {"id": str(item["query_id"])})
            row["top_changed_term"] = {
                "class": item["class"],
                "target_id": item["target_id"],
                "target_name": item["target_name"],
                "delta": item["delta_neg_log10_p_adjust"],
                "left_p_adjust": item["left_p_adjust"],
                "right_p_adjust": item["right_p_adjust"],
                "left_target_size": item["left_target_size"],
                "right_target_size": item["right_target_size"],
            }

        crossing = con.execute(
            f"""
            SELECT
              query_id,
              class,
              target_id,
              target_name,
              delta_neg_log10_p_adjust,
              left_p_adjust,
              right_p_adjust,
              left_target_size,
              right_target_size
            FROM (
              SELECT
                *,
                row_number() OVER (
                  PARTITION BY query_id
                  ORDER BY abs(coalesce(delta_neg_log10_p_adjust, 0)) DESC,
                           coalesce(left_p_adjust, right_p_adjust, 1.0) ASC
                ) AS rank
              FROM read_parquet(?)
              WHERE class IN ('lost_significant', 'gained_significant')
                AND CASE
                      WHEN class = 'lost_significant' THEN coalesce(left_target_size, 0)
                      ELSE coalesce(right_target_size, 0)
                    END <= {SPECIFIC_TARGET_SIZE}
            )
            WHERE rank = 1
            """,
            [str(path)],
        )
        for item in row_dicts(crossing):
            row = rows.setdefault(str(item["query_id"]), {"id": str(item["query_id"])})
            row["top_crossing_term"] = {
                "class": item["class"],
                "target_id": item["target_id"],
                "target_name": item["target_name"],
                "delta": item["delta_neg_log10_p_adjust"],
                "left_p_adjust": item["left_p_adjust"],
                "right_p_adjust": item["right_p_adjust"],
                "left_target_size": item["left_target_size"],
                "right_target_size": item["right_target_size"],
            }

    def list_query_sets(
        self,
        *,
        search: str | None,
        source_class: str | None,
        stratum: str | None,
        sort: str,
        order: str,
        offset: int,
        limit: int,
    ) -> dict[str, Any]:
        rows = self.rows
        if search:
            needle = search.casefold()
            rows = [
                row
                for row in rows
                if needle in str(row.get("id") or "").casefold()
                or needle in str(row.get("name") or "").casefold()
                or needle in str(row.get("description") or "").casefold()
            ]
        if source_class:
            rows = [row for row in rows if row.get("source_class") == source_class]
        if stratum:
            rows = [row for row in rows if row.get("stratum") == stratum]

        reverse = order == "desc"

        def sort_key(row: dict[str, Any]) -> tuple[int, Any]:
            value = row.get(sort)
            if isinstance(value, dict):
                value = value.get("delta") or value.get("target_name")
            if value is None:
                return (1, "")
            if isinstance(value, str):
                return (0, value.casefold())
            return (0, value)

        rows = sorted(rows, key=sort_key, reverse=reverse)
        limit = max(1, min(limit, MAX_LIMIT))
        offset = max(0, offset)
        return {
            "total": len(rows),
            "offset": offset,
            "limit": limit,
            "rows": rows[offset : offset + limit],
        }

    def query_set(self, query_id: str) -> dict[str, Any]:
        for row in self.rows:
            if row.get("id") == query_id:
                item = dict(row)
                item["genes"] = self.query_genes.get(query_id, [])
                return item
        raise KeyError(query_id)

    def result_rows(self, query_id: str, run: str, limit: int) -> dict[str, Any]:
        output_key = {"old": "old_results", "new": "new_results"}.get(run)
        if output_key is None:
            raise ValueError("run must be old or new")
        path = self.outputs.get(output_key)
        if path is None or not path.exists():
            return {"rows": [], "total": 0}
        duckdb = require_duckdb()
        con = duckdb.connect()
        limit = max(1, min(limit, MAX_LIMIT))
        result = con.execute(
            """
            SELECT *
            FROM read_parquet(?)
            WHERE query_id = ?
            ORDER BY p_adjust_bonferroni ASC, overlap DESC, target_name ASC
            LIMIT ?
            """,
            [str(path), query_id, limit],
        )
        rows = row_dicts(result)
        count = con.execute(
            """
            SELECT count(*)
            FROM read_parquet(?)
            WHERE query_id = ?
            """,
            [str(path), query_id],
        ).fetchone()[0]
        return {"rows": rows, "total": count}

    def diff_rows(self, query_id: str, class_filter: str | None, limit: int) -> dict[str, Any]:
        path = self.outputs.get("diff_results")
        if path is None or not path.exists():
            return {"rows": [], "total": 0}
        duckdb = require_duckdb()
        con = duckdb.connect()
        limit = max(1, min(limit, MAX_LIMIT))
        params: list[Any] = [str(path), query_id]
        where = "query_id = ?"
        if class_filter:
            where += " AND class = ?"
            params.append(class_filter)
        query = f"""
            SELECT *
            FROM read_parquet(?)
            WHERE {where}
            ORDER BY
              CASE class
                WHEN 'lost_significant' THEN 0
                WHEN 'gained_significant' THEN 1
                ELSE 2
              END,
              abs(coalesce(delta_neg_log10_p_adjust, 0)) DESC,
              coalesce(left_p_adjust, right_p_adjust, 1.0) ASC,
              target_name ASC
            LIMIT ?
        """
        result = con.execute(query, [*params, limit])
        rows = row_dicts(result)
        count = con.execute(
            f"""
            SELECT count(*)
            FROM read_parquet(?)
            WHERE {where}
            """,
            params,
        ).fetchone()[0]
        return {"rows": rows, "total": count}


@dataclass
class ExplorerState:
    bundles: dict[str, ReportBundle]

    @classmethod
    def from_paths(cls, paths: list[Path]) -> "ExplorerState":
        bundles = [ReportBundle(path=path, ordinal=index) for index, path in enumerate(paths)]
        if not bundles:
            raise ValueError("No report bundles found")
        return cls({bundle.id: bundle for bundle in bundles})

    def bundle(self, bundle_id: str) -> ReportBundle:
        try:
            return self.bundles[bundle_id]
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=f"Unknown bundle: {bundle_id}") from exc


def create_app(bundle_paths: list[Path]) -> FastAPI:
    state = ExplorerState.from_paths(discover_bundle_paths(bundle_paths))
    app = FastAPI(
        title="genesets-rs Explorer",
        description="Local browser for genesets-rs workflow result bundles.",
        version="0.1.0",
    )
    static_dir = Path(__file__).resolve().parent / "static"
    app.mount("/assets", StaticFiles(directory=str(static_dir)), name="assets")

    @app.get("/", include_in_schema=False)
    def index() -> FileResponse:
        return FileResponse(str(static_dir / "index.html"))

    @app.get("/api/bundles")
    def bundles() -> dict[str, Any]:
        return {"bundles": [bundle.overview() for bundle in state.bundles.values()]}

    @app.get("/api/bundles/{bundle_id}/facets")
    def facets(bundle_id: str) -> dict[str, Any]:
        return state.bundle(bundle_id).facets()

    @app.get("/api/bundles/{bundle_id}/query-sets")
    def query_sets(
        bundle_id: str,
        search: str | None = None,
        source_class: str | None = None,
        stratum: str | None = None,
        sort: str = "lost_specific_count",
        order: str = "desc",
        offset: int = 0,
        limit: int = DEFAULT_LIMIT,
    ) -> dict[str, Any]:
        if order not in {"asc", "desc"}:
            raise HTTPException(status_code=400, detail="order must be asc or desc")
        return state.bundle(bundle_id).list_query_sets(
            search=search,
            source_class=source_class,
            stratum=stratum,
            sort=sort,
            order=order,
            offset=offset,
            limit=limit,
        )

    @app.get("/api/bundles/{bundle_id}/query-set")
    def query_set(bundle_id: str, query_id: str = Query(...)) -> dict[str, Any]:
        try:
            return state.bundle(bundle_id).query_set(query_id)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=f"Unknown query set: {query_id}") from exc

    @app.get("/api/bundles/{bundle_id}/query-set/results")
    def query_results(
        bundle_id: str,
        query_id: str = Query(...),
        run: str = "old",
        limit: int = 200,
    ) -> dict[str, Any]:
        try:
            return state.bundle(bundle_id).result_rows(query_id, run, limit)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.get("/api/bundles/{bundle_id}/query-set/diffs")
    def query_diffs(
        bundle_id: str,
        query_id: str = Query(...),
        class_filter: str | None = None,
        limit: int = 200,
    ) -> dict[str, Any]:
        return state.bundle(bundle_id).diff_rows(query_id, class_filter, limit)

    return app
