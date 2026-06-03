from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


QUERY_URL = "https://mygeneset.info/v1/query"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="genesets-workflows fetch-mygeneset",
        description="Fetch MyGeneset.info query results into GMT for eval workloads."
    )
    parser.add_argument("--query", required=True, help="MyGeneset query string.")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--size", type=int, default=1000, help="Page size.")
    parser.add_argument("--limit", type=int, default=0, help="Maximum sets; 0 means all hits.")
    parser.add_argument("--gene-field", default="symbol")
    parser.add_argument("--species", default="human")
    parser.add_argument(
        "--source-filter",
        default=None,
        help="Optional client-side exact source filter.",
    )
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--max-genes", type=int, default=2000)
    parser.add_argument("--sleep", type=float, default=0.0)
    return parser.parse_args(argv)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def request_json(params: dict[str, Any]) -> dict[str, Any]:
    url = f"{QUERY_URL}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(url, headers={"User-Agent": "genesets-workflows/0.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        data = json.load(response)
    data["_source_url"] = url
    return data


def normalize_gene_values(raw_genes: object, gene_field: str) -> list[str]:
    if raw_genes is None:
        return []
    if isinstance(raw_genes, dict):
        raw_genes = [raw_genes]
    if not isinstance(raw_genes, list):
        return []
    values: set[str] = set()
    for gene in raw_genes:
        if not isinstance(gene, dict):
            continue
        value = gene.get(gene_field)
        if isinstance(value, list):
            values.update(str(item).strip() for item in value if str(item).strip())
        elif value is not None and str(value).strip():
            values.add(str(value).strip())
    return sorted(values)


def fetch_query(args: argparse.Namespace) -> dict[str, Any]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = f"_id,name,source,taxid,genes.{args.gene_field}"
    fetched = 0
    emitted = 0
    skipped_small = 0
    skipped_large = 0
    skipped_source = 0
    total = 0
    background: set[str] = set()
    sets: list[dict[str, Any]] = []

    with (args.out_dir / "queries.gmt").open("w", newline="") as gmt:
        writer = csv.writer(gmt, delimiter="\t", lineterminator="\n")
        while True:
            remaining = args.limit - fetched if args.limit else args.size
            page_size = min(args.size, remaining) if args.limit else args.size
            if page_size <= 0:
                break
            data = request_json(
                {
                    "q": args.query,
                    "species": args.species,
                    "fields": fields,
                    "size": page_size,
                    "from": fetched,
                }
            )
            hits = data.get("hits", [])
            if not hits:
                total = data.get("total", fetched)
                break
            total = data.get("total", fetched + len(hits))
            for hit in hits:
                fetched += 1
                if args.source_filter and hit.get("source") != args.source_filter:
                    skipped_source += 1
                    continue
                genes = normalize_gene_values(hit.get("genes"), args.gene_field)
                if len(genes) < args.min_genes:
                    skipped_small += 1
                    continue
                if len(genes) > args.max_genes:
                    skipped_large += 1
                    continue
                emitted += 1
                background.update(genes)
                geneset_id = hit["_id"]
                name = hit.get("name") or geneset_id
                writer.writerow([geneset_id, name, *genes])
                sets.append(
                    {
                        "id": geneset_id,
                        "name": name,
                        "source": hit.get("source"),
                        "taxid": hit.get("taxid"),
                        "gene_count": len(genes),
                    }
                )
                if args.limit and fetched >= args.limit:
                    break
            if fetched >= total or (args.limit and fetched >= args.limit):
                break
            if args.sleep:
                time.sleep(args.sleep)

    with (args.out_dir / "background.txt").open("w") as handle:
        handle.write("gene_id\n")
        for gene in sorted(background):
            handle.write(f"{gene}\n")

    metadata = {
        "generated_at_utc": utc_now(),
        "service": "MyGeneset.info",
        "query": args.query,
        "species": args.species,
        "gene_field": args.gene_field,
        "limit": args.limit,
        "total_hits_reported": total,
        "fetched": fetched,
        "emitted": emitted,
        "skipped_small": skipped_small,
        "skipped_large": skipped_large,
        "skipped_source": skipped_source,
        "source_filter": args.source_filter,
        "min_genes": args.min_genes,
        "max_genes": args.max_genes,
        "background_gene_count": len(background),
        "sets": sets,
    }
    with (args.out_dir / "metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return metadata


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    metadata = fetch_query(args)
    print(
        f"Fetched {metadata['fetched']} hits, emitted {metadata['emitted']} sets "
        f"to {args.out_dir / 'queries.gmt'}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
