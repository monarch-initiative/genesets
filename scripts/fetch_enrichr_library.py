#!/usr/bin/env python3
"""Fetch an Enrichr gene-set library as a genesets-rs query GMT.

Enrichr (maayanlab.cloud) redistributes many external libraries, including
recent single-cell marker collections (Tabula_Sapiens, Azimuth_2023,
CellMarker_2024, ...). Set names are kept; the GMT description column records
the library name.

    python3 scripts/fetch_enrichr_library.py Tabula_Sapiens --out /tmp/eval/breadth/sc_tabula_sapiens/queries.gmt
    python3 scripts/fetch_enrichr_library.py CellMarker_2024 --include Human --out ...
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("library", help="Enrichr library name, e.g. Tabula_Sapiens")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--include", help="keep only sets whose name matches this regex (e.g. Human)")
    parser.add_argument("--exclude", help="drop sets whose name matches this regex (e.g. Mouse)")
    parser.add_argument("--min-genes", type=int, default=10)
    parser.add_argument("--max-genes", type=int, default=1000)
    args = parser.parse_args()

    query = urllib.parse.urlencode({"mode": "text", "libraryName": args.library})
    request = urllib.request.Request(f"{URL}?{query}", headers={"User-Agent": "genesets-eval/0.1"})
    with urllib.request.urlopen(request, timeout=300) as response:
        text = response.read().decode("utf-8")

    kept = skipped = 0
    lines = []
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name = parts[0].strip()
        if args.include and not re.search(args.include, name):
            skipped += 1
            continue
        if args.exclude and re.search(args.exclude, name):
            skipped += 1
            continue
        # Enrichr entries may carry a ",weight" suffix
        genes = sorted({g.split(",")[0].strip() for g in parts[2:] if g.strip()})
        if not args.min_genes <= len(genes) <= args.max_genes:
            skipped += 1
            continue
        lines.append("\t".join([name, args.library, *genes]))
        kept += 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    print(f"{args.library}: kept {kept} sets, skipped {skipped} -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
