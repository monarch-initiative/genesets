#!/usr/bin/env python3
"""Fetch MyGeneset.info gene sets listed in a TSV manifest.

The script uses only the Python standard library so eval refreshes do not
require a notebook or package environment.
"""

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


BASE_URL = "https://mygeneset.info/v1/geneset"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--gene-field",
        default="symbol",
        help="Gene field to emit from MyGeneset records, for example symbol or ncbigene.",
    )
    parser.add_argument("--sleep", type=float, default=0.0, help="Delay between API calls.")
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Reuse an existing queries.gmt, background.txt, and metadata.json snapshot.",
    )
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip() and not line.startswith("#")),
            delimiter="\t",
        )
        rows = list(reader)
    required = {"geneset_id", "name", "source", "taxid", "expected_min_symbols"}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"{path} is missing required columns: {sorted(missing)}")
    return rows


def fetch_geneset(geneset_id: str, gene_field: str) -> dict:
    params = urllib.parse.urlencode(
        {"fields": f"_id,name,source,taxid,genes.{gene_field}"}
    )
    url = f"{BASE_URL}/{urllib.parse.quote(geneset_id)}?{params}"
    with urllib.request.urlopen(url, timeout=30) as response:
        if response.status != 200:
            raise RuntimeError(f"GET {url} failed with status {response.status}")
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
        if value is None:
            continue
        if isinstance(value, list):
            values.update(str(item).strip() for item in value if str(item).strip())
        else:
            value = str(value).strip()
            if value:
                values.add(value)
    return sorted(values)


def write_outputs(
    out_dir: Path,
    gene_field: str,
    fetched: list[dict],
    manifest_rows: list[dict[str, str]],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_by_id = {row["geneset_id"]: row for row in manifest_rows}

    metadata = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "service": "MyGeneset.info",
        "base_url": BASE_URL,
        "gene_field": gene_field,
        "sets": [],
    }

    background: set[str] = set()
    with (out_dir / "queries.gmt").open("w", newline="") as gmt:
        for item in fetched:
            geneset_id = item["_id"]
            manifest = manifest_by_id[geneset_id]
            name = item.get("name") or manifest["name"]
            genes = item["_normalized_genes"]
            background.update(genes)
            gmt.write("\t".join([geneset_id, name, *genes]) + "\n")
            metadata["sets"].append(
                {
                    "id": geneset_id,
                    "name": name,
                    "source": item.get("source"),
                    "taxid": item.get("taxid"),
                    "gene_count": len(genes),
                    "selection_note": manifest.get("selection_note", ""),
                    "category": manifest.get("category", ""),
                    "contrast_type": manifest.get("contrast_type", ""),
                    "direction": manifest.get("direction", ""),
                    "source_url": item["_source_url"],
                }
            )

    with (out_dir / "background.txt").open("w", newline="") as background_file:
        background_file.write("gene_id\n")
        for gene in sorted(background):
            background_file.write(f"{gene}\n")

    with (out_dir / "metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    args = parse_args()
    expected_outputs = [
        args.out_dir / "queries.gmt",
        args.out_dir / "background.txt",
        args.out_dir / "metadata.json",
    ]
    if args.skip_existing and all(path.exists() for path in expected_outputs):
        print(f"Using existing MyGeneset snapshot at {args.out_dir}")
        return 0

    rows = read_manifest(args.manifest)
    fetched: list[dict] = []
    failures: list[str] = []

    for row in rows:
        geneset_id = row["geneset_id"]
        try:
            item = fetch_geneset(geneset_id, args.gene_field)
            genes = normalize_gene_values(item.get("genes"), args.gene_field)
            item["_normalized_genes"] = genes
            expected = int(row["expected_min_symbols"])
            if len(genes) < expected:
                failures.append(
                    f"{geneset_id} has {len(genes)} {args.gene_field} values; "
                    f"expected at least {expected}"
                )
            fetched.append(item)
        except Exception as exc:  # noqa: BLE001 - script reports all failed IDs.
            failures.append(f"{geneset_id}: {exc}")
        if args.sleep:
            time.sleep(args.sleep)

    if failures:
        print("Fetch completed with validation failures:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    write_outputs(args.out_dir, args.gene_field, fetched, rows)
    print(f"Wrote {len(fetched)} gene sets to {args.out_dir / 'queries.gmt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
