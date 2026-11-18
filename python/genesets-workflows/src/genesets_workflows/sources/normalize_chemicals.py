"""Normalize a metabolite list to the canonical ChEBI ids used by prepare-chem2go.

Run this on both the query (changed metabolites) and the background (every
metabolite the platform measured) so they share one identifier space with the
annotation tables. Inputs are one id per line, or a TSV whose first column is
the id (a header row is skipped when its first cell does not parse).
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

from genesets_workflows.chemid import ChemicalNormalizer
from genesets_workflows.sources.chem2go import DEFAULT_OUT_DIR, write_tsv


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="genesets-workflows normalize-chemicals",
        description="Map metabolite ids (ChEBI, HMDB, KEGG, LIPID MAPS, CAS, ...) to canonical ChEBI.",
    )
    parser.add_argument("input", type=Path, help="List or TSV of chemical ids.")
    parser.add_argument("--prep-dir", type=Path, default=DEFAULT_OUT_DIR, help="prepare-chem2go output directory.")
    parser.add_argument("--out", type=Path, required=True, help="Canonical id list for genesets-rs.")
    parser.add_argument("--report", type=Path, help="Per-input mapping report TSV.")
    parser.add_argument(
        "--extra-xrefs",
        type=Path,
        help="Supplementary two-column TSV (xref, chebi_id) for ids ChEBI does not cross-reference.",
    )
    return parser.parse_args(argv)


def read_ids(path: Path) -> list[str]:
    ids = []
    with path.open() as handle:
        for line in handle:
            value = line.rstrip("\n").split("\t")[0].strip()
            if value and not value.startswith("#"):
                ids.append(value)
    return ids


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    normalizer = ChemicalNormalizer.load(args.prep_dir, args.extra_xrefs)
    results = [normalizer.normalize(raw) for raw in read_ids(args.input)]
    if results and results[0].status == "unparsed":
        results = results[1:]  # header row

    seen: dict[str, str] = {}
    for result in results:
        if result.canonical_id and result.canonical_id not in seen:
            seen[result.canonical_id] = result.input_id
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as handle:
        for canonical in seen:
            handle.write(f"{canonical}\n")

    if args.report:
        write_tsv(
            args.report,
            ["input_id", "curie", "canonical_id", "canonical_name", "status", "annotated", "collapsed_with", "alternatives"],
            (
                (
                    r.input_id,
                    r.curie or "",
                    r.canonical_id or "",
                    normalizer.names.get(r.canonical_id or "", ""),
                    r.status,
                    str(r.canonical_id in normalizer.annotated).lower() if r.canonical_id else "",
                    seen.get(r.canonical_id or "", "") if seen.get(r.canonical_id or "") != r.input_id else "",
                    "|".join(r.alternatives),
                )
                for r in results
            ),
        )

    statuses = Counter(r.status for r in results)
    summary = ", ".join(f"{k}={v}" for k, v in sorted(statuses.items()))
    unannotated = sum(1 for c in seen if c not in normalizer.annotated)
    print(
        f"{len(results)} inputs -> {len(seen)} canonical chemicals ({summary}); "
        f"{unannotated} mapped but carry no GO annotation",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
