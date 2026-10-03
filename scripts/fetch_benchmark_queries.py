#!/usr/bin/env python3
"""Build queries.gmt for the curated benchmark: member genes of every evaluable set.

MSigDB sets in ``curation/genesets/*.yaml`` are fetched from MyGeneset.info by
their MSigDB name (``_id``); ``LIT:`` sets come from
``curation/genesets/lit_members.gmt``. The first GMT field is the YAML
``gene_set_name`` so ``score_method_vs_benchmark.py`` can join on it.

    python3 scripts/fetch_benchmark_queries.py --out /tmp/eval/queries.gmt
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
GENESETS = ROOT / "curation" / "genesets"
BASE_URL = "https://mygeneset.info/v1/geneset"


def fetch_symbols(name: str) -> list[str] | None:
    url = f"{BASE_URL}/{urllib.parse.quote(name)}?{urllib.parse.urlencode({'fields': 'genes.symbol,source'})}"
    request = urllib.request.Request(url, headers={"User-Agent": "genesets-benchmark/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            data = json.load(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise
    if data.get("source") != "msigdb":
        return None
    genes = data.get("genes") or []
    if isinstance(genes, dict):
        genes = [genes]
    return sorted({g["symbol"] for g in genes if isinstance(g, dict) and g.get("symbol")})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    lines: list[str] = []
    missing: list[str] = []
    for path in sorted(GENESETS.glob("*.yaml")):
        with path.open() as handle:
            doc = yaml.safe_load(handle) or {}
        set_id = str(doc.get("gene_set_id", ""))
        name = doc.get("gene_set_name")
        if not set_id.startswith("MSIGDB:") or not name:
            continue
        symbols = fetch_symbols(set_id.split(":", 1)[1])
        if not symbols:
            missing.append(name)
            continue
        lines.append("\t".join([name, set_id, *symbols]))

    lit = (GENESETS / "lit_members.gmt").read_text().splitlines()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines + [line for line in lit if line.strip()]) + "\n")
    print(f"{len(lines)} MSigDB + {len(lit)} LIT sets -> {args.out}; {len(missing)} MSigDB not on MyGeneset", file=sys.stderr)
    for name in missing:
        print(f"  missing: {name}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
