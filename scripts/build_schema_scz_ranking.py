#!/usr/bin/env python3
"""Build the SCHEMA schizophrenia exome ranking for the curated corpus.

Source: Singh et al. 2022, "Rare coding variants in ten genes confer
substantial risk for schizophrenia", Nature 604:509-516 (PMID:35396579), via
the SCHEMA browser's gene-results download. Each gene's score is
-log10(`P meta`), the meta-analysis of case-control and de novo rare-variant
burden. That p-value is one-sided (excess damaging variation in cases), so the
ranking is unsigned.

SCHEMA identifies genes by Ensembl id; the curated corpus and GOA use HGNC
symbols, so ids are mapped through HGNC's complete set (approved symbols
only). Genes with no `P meta`, or no approved HGNC symbol, are left out; if
two Ensembl ids map to one symbol, the more significant is kept.

Usage::

    uv run --project python/genesets-workflows --extra curation \\
        python scripts/build_schema_scz_ranking.py \\
        --out curation/rankings/SCHIZOPHRENIA_SCHEMA_EXOME.rnk.gz
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python" / "genesets-workflows" / "src"))

from genesets_workflows.curation.ranking import write_ranking  # noqa: E402

SCHEMA_URL = "https://atgu-exome-browser-data.s3.amazonaws.com/SCHEMA/SCHEMA_gene_results.tsv.bgz"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cache-dir", type=Path, default=Path(".cache/schema_scz"))
    return parser.parse_args()


def download(url: str, path: Path) -> Path:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(url, headers={"User-Agent": "genesets-curation/0.1"})
        with urllib.request.urlopen(request, timeout=120) as response:
            path.write_bytes(response.read())
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    print(f"{url}\n  sha256 {digest}")
    return path


def ensembl_to_symbol(path: Path) -> dict[str, str]:
    mapping = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] == "Approved" and row["ensembl_gene_id"]:
                mapping[row["ensembl_gene_id"]] = row["symbol"]
    return mapping


def main() -> None:
    args = parse_args()
    schema = download(SCHEMA_URL, args.cache_dir / "SCHEMA_gene_results.tsv.bgz")
    hgnc = download(HGNC_URL, args.cache_dir / "hgnc_complete_set.txt")
    symbols = ensembl_to_symbol(hgnc)

    scores: dict[str, float] = {}
    rows = no_p = unmapped = 0
    with gzip.open(schema, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows += 1
            raw = row["P meta"].strip()
            if raw in ("", "NA"):
                no_p += 1
                continue
            symbol = symbols.get(row["gene_id"])
            if symbol is None:
                unmapped += 1
                continue
            score = -math.log10(float(raw))
            scores[symbol] = max(score, scores.get(symbol, score))

    digest = write_ranking(scores, args.out)
    print(
        f"{rows} SCHEMA genes: {no_p} without P meta, {unmapped} without an approved "
        f"HGNC symbol, {len(scores)} ranked"
    )
    print(f"wrote {args.out}\n  sha256 {digest}")


if __name__ == "__main__":
    main()
