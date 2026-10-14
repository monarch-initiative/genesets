#!/usr/bin/env python3
"""Build ``queries.gmt`` for the curated benchmark (step 1 of the benchmark evals).

MSigDB sets (``gene_set_id: MSIGDB:<NAME>``) are fetched from mygeneset.info by
exact ``_id``; ``LIT:`` sets are taken from ``curation/genesets/lit_members.gmt``
when they have a captured membership. Sets absent from both are reported and
skipped. Standard library only.
"""

from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QUERY_URL = "https://mygeneset.info/v1/query"


def fetch_msigdb(name: str) -> list[str]:
    params = {"q": f"_id:{name}", "fields": "_id,genes.symbol", "size": 1}
    req = urllib.request.Request(
        f"{QUERY_URL}?{urllib.parse.urlencode(params)}", headers={"User-Agent": "genesets-rs-eval/0.1"}
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = json.load(resp)
    for hit in data.get("hits", []):
        if hit.get("_id") == name:
            genes = hit.get("genes") or []
            if isinstance(genes, dict):
                genes = [genes]
            return sorted({g["symbol"] for g in genes if isinstance(g, dict) and g.get("symbol")})
    return []


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genesets-dir", type=Path, default=ROOT / "curation" / "genesets")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    lit = {}
    for line in (args.genesets_dir / "lit_members.gmt").read_text().splitlines():
        parts = line.split("\t")
        if len(parts) > 2:
            lit[parts[0]] = parts[2:]

    rows, missing = [], []
    for path in sorted(args.genesets_dir.glob("*.yaml")):
        sid = name = None
        for line in path.read_text().splitlines()[:5]:
            if line.startswith("gene_set_id:"):
                sid = line.split(":", 1)[1].strip()
            elif line.startswith("gene_set_name:"):
                name = line.split(":", 1)[1].strip()
        if not sid or not name:
            continue
        genes = lit.get(name, []) if sid.startswith("LIT:") else fetch_msigdb(name)
        if genes:
            rows.append("\t".join([name, sid, *genes]))
        else:
            missing.append(name)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(rows) + "\n")
    print(f"wrote {len(rows)} sets to {args.out}; {len(missing)} without membership")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
