from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from pathlib import Path
from typing import Any

import yaml

from genesets_workflows.sources.mygeneset import (
    classify_msigdb_set,
    normalize_gene_values,
    request_json,
    utc_now,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="genesets-workflows fetch-mygeneset-stratified",
        description="Fetch a stratified MyGeneSet query collection into GMT for eval workloads.",
    )
    parser.add_argument("config", type=Path, help="YAML stratification plan.")
    parser.add_argument("--out-dir", type=Path, default=None, help="Override output.dir from the plan.")
    parser.add_argument("--limit", type=int, default=None, help="Override output.limit from the plan.")
    parser.add_argument("--dry-run", action="store_true", help="Validate and print the plan without fetching.")
    return parser.parse_args(argv)


def load_plan(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        plan = yaml.safe_load(handle) or {}
    if not isinstance(plan, dict):
        raise SystemExit(f"{path} must contain a YAML mapping")
    if not isinstance(plan.get("strata"), list) or not plan["strata"]:
        raise SystemExit(f"{path} must define a non-empty strata list")
    return plan


def plan_value(plan: dict[str, Any], key: str, default: Any) -> Any:
    return plan.get(key, default)


def normalize_plan(plan: dict[str, Any], base: Path, args: argparse.Namespace) -> dict[str, Any]:
    output = dict(plan.get("output") or {})
    if args.out_dir is not None:
        output["dir"] = args.out_dir
    if args.limit is not None:
        output["limit"] = args.limit
    if "dir" not in output:
        raise SystemExit("plan output.dir is required")

    output["dir"] = Path(output["dir"])
    if not output["dir"].is_absolute():
        output["dir"] = (base / output["dir"]).resolve()
    output["limit"] = int(output.get("limit", 0))

    defaults = {
        "species": plan_value(plan, "species", "human"),
        "gene_field": plan_value(plan, "gene_field", "symbol"),
        "source_filter": plan.get("source_filter"),
        "min_genes": int(plan_value(plan, "min_genes", 20)),
        "max_genes": int(plan_value(plan, "max_genes", 1000)),
        "size": int(plan_value(plan, "size", 1000)),
        "sleep": float(plan_value(plan, "sleep", 0.0)),
        "fetch_limit_per_query": int(plan_value(plan, "fetch_limit_per_query", 2000)),
    }

    normalized_strata = []
    for item in plan["strata"]:
        if not isinstance(item, dict):
            raise SystemExit("each stratum must be a YAML mapping")
        label = str(item.get("label") or "").strip()
        if not label:
            raise SystemExit("each stratum requires label")
        queries = item.get("queries", item.get("query"))
        if isinstance(queries, str):
            queries = [queries]
        if not isinstance(queries, list) or not queries:
            raise SystemExit(f"stratum {label} requires query or queries")
        quota = int(item.get("quota", 0))
        if quota <= 0:
            raise SystemExit(f"stratum {label} requires positive quota")
        normalized = dict(defaults)
        include_regex = item.get("include_regex", item.get("require_regex"))
        exclude_regex = item.get("exclude_regex")
        if isinstance(include_regex, str):
            include_regex = [include_regex]
        if isinstance(exclude_regex, str):
            exclude_regex = [exclude_regex]
        normalized.update(
            {
                "label": label,
                "quota": quota,
                "queries": [str(query) for query in queries],
                "source_filter": item.get("source_filter", defaults["source_filter"]),
                "min_genes": int(item.get("min_genes", defaults["min_genes"])),
                "max_genes": int(item.get("max_genes", defaults["max_genes"])),
                "fetch_limit_per_query": int(item.get("fetch_limit_per_query", defaults["fetch_limit_per_query"])),
                "include_regex": [str(pattern) for pattern in include_regex or []],
                "exclude_regex": [str(pattern) for pattern in exclude_regex or []],
            }
        )
        normalized_strata.append(normalized)

    return {
        "description": plan.get("description"),
        "output": output,
        "defaults": defaults,
        "strata": normalized_strata,
    }


def iter_hits(
    query: str,
    *,
    species: str,
    gene_field: str,
    size: int,
    fetch_limit: int,
    sleep: float,
) -> tuple[int, int, list[dict[str, Any]]]:
    fields = f"_id,name,source,taxid,description,genes.{gene_field}"
    fetched = 0
    total = 0
    hits_out: list[dict[str, Any]] = []
    while True:
        remaining = fetch_limit - fetched if fetch_limit else size
        page_size = min(size, remaining) if fetch_limit else size
        if page_size <= 0:
            break
        data = request_json(
            {
                "q": query,
                "species": species,
                "fields": fields,
                "size": page_size,
                "from": fetched,
            }
        )
        hits = data.get("hits", [])
        total = int(data.get("total", fetched + len(hits)))
        if not hits:
            break
        hits_out.extend(hits)
        fetched += len(hits)
        if fetched >= total or (fetch_limit and fetched >= fetch_limit):
            break
        if sleep:
            time.sleep(sleep)
    return fetched, total, hits_out


def hit_text(hit: dict[str, Any]) -> str:
    return " ".join(
        str(hit.get(key) or "")
        for key in ("_id", "name", "description", "source")
    )


def fetch_stratified(plan: dict[str, Any]) -> dict[str, Any]:
    out_dir: Path = plan["output"]["dir"]
    global_limit = int(plan["output"].get("limit", 0))
    out_dir.mkdir(parents=True, exist_ok=True)

    seen_ids: set[str] = set()
    background: set[str] = set()
    emitted_sets: list[dict[str, Any]] = []
    strata_metadata: list[dict[str, Any]] = []
    skipped_duplicate = 0
    skipped_source = 0
    skipped_include = 0
    skipped_exclude = 0
    skipped_small = 0
    skipped_large = 0
    fetched_total = 0

    with (out_dir / "queries.gmt").open("w", newline="") as gmt:
        writer = csv.writer(gmt, delimiter="\t", lineterminator="\n")
        for stratum in plan["strata"]:
            stratum_emitted = 0
            stratum_fetched = 0
            query_metadata = []
            for query in stratum["queries"]:
                if stratum_emitted >= stratum["quota"]:
                    break
                if global_limit and len(emitted_sets) >= global_limit:
                    break
                fetched, total, hits = iter_hits(
                    query,
                    species=stratum["species"],
                    gene_field=stratum["gene_field"],
                    size=stratum["size"],
                    fetch_limit=stratum["fetch_limit_per_query"],
                    sleep=stratum["sleep"],
                )
                fetched_total += fetched
                stratum_fetched += fetched
                query_emitted = 0
                for hit in hits:
                    if stratum_emitted >= stratum["quota"]:
                        break
                    if global_limit and len(emitted_sets) >= global_limit:
                        break
                    hit_id = str(hit.get("_id", "")).strip()
                    if not hit_id:
                        continue
                    if hit_id in seen_ids:
                        skipped_duplicate += 1
                        continue
                    if stratum["source_filter"] and hit.get("source") != stratum["source_filter"]:
                        skipped_source += 1
                        continue
                    if stratum["include_regex"] and not any(
                        re.search(pattern, hit_text(hit), flags=re.IGNORECASE)
                        for pattern in stratum["include_regex"]
                    ):
                        skipped_include += 1
                        continue
                    if stratum["exclude_regex"] and any(
                        re.search(pattern, hit_text(hit), flags=re.IGNORECASE)
                        for pattern in stratum["exclude_regex"]
                    ):
                        skipped_exclude += 1
                        continue
                    genes = normalize_gene_values(hit.get("genes"), stratum["gene_field"])
                    if len(genes) < stratum["min_genes"]:
                        skipped_small += 1
                        continue
                    if len(genes) > stratum["max_genes"]:
                        skipped_large += 1
                        continue

                    seen_ids.add(hit_id)
                    background.update(genes)
                    name = hit.get("name") or hit_id
                    description = hit.get("description")
                    writer.writerow([hit_id, name, *genes])
                    item = {
                        "id": hit_id,
                        "name": name,
                        "description": description,
                        "source": hit.get("source"),
                        "source_class": classify_msigdb_set(hit_id, name, description)
                        if hit.get("source") == "msigdb"
                        else hit.get("source"),
                        "taxid": hit.get("taxid"),
                        "gene_count": len(genes),
                        "stratum": stratum["label"],
                        "search_query": query,
                    }
                    emitted_sets.append(item)
                    stratum_emitted += 1
                    query_emitted += 1
                query_metadata.append(
                    {
                        "query": query,
                        "total_hits_reported": total,
                        "fetched": fetched,
                        "emitted": query_emitted,
                        "include_regex": stratum["include_regex"],
                        "exclude_regex": stratum["exclude_regex"],
                    }
                )
                if stratum["sleep"]:
                    time.sleep(stratum["sleep"])
            strata_metadata.append(
                {
                    "label": stratum["label"],
                    "quota": stratum["quota"],
                    "emitted": stratum_emitted,
                    "fetched": stratum_fetched,
                    "queries": query_metadata,
                }
            )
            if global_limit and len(emitted_sets) >= global_limit:
                break

    with (out_dir / "background.txt").open("w") as handle:
        handle.write("gene_id\n")
        for gene in sorted(background):
            handle.write(f"{gene}\n")

    metadata = {
        "generated_at_utc": utc_now(),
        "service": "MyGeneSet.info",
        "strategy": "stratified_query_quota",
        "description": plan.get("description"),
        "output_limit": global_limit,
        "emitted": len(emitted_sets),
        "fetched": fetched_total,
        "background_gene_count": len(background),
        "skipped_duplicate": skipped_duplicate,
        "skipped_source": skipped_source,
        "skipped_include": skipped_include,
        "skipped_exclude": skipped_exclude,
        "skipped_small": skipped_small,
        "skipped_large": skipped_large,
        "defaults": plan["defaults"],
        "strata": strata_metadata,
        "sets": emitted_sets,
    }
    with (out_dir / "metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return metadata


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    raw_plan = load_plan(args.config)
    plan = normalize_plan(raw_plan, args.config.parent, args)
    if args.dry_run:
        print(json.dumps(plan, indent=2, default=str, sort_keys=True))
        return 0
    metadata = fetch_stratified(plan)
    print(
        f"Fetched {metadata['fetched']} hits, emitted {metadata['emitted']} sets "
        f"to {plan['output']['dir'] / 'queries.gmt'}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
