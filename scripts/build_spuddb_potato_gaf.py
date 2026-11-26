#!/usr/bin/env python3
"""Build a gene-level potato (DM v6.1) GAF from SpudDB's GO annotation.

The SOLTU_* lunar-regolith sets in curation/plants/genesets use DM v6.1 gene ids
(``Soltu.DM.01G000010``), as in Handy et al. 2026 (PMID:43058136). UniProt's
potato proteomes are built on older assemblies and map to only about half of the
v6.1 genes by sequence, so the eval uses SpudDB's own v6.1 GO file instead (the
annotation the authors used). That file is GAF-shaped but transcript-level, with
the DB column set to ``MSU``. Every row is IEA, projected from the gene's
Arabidopsis best hit (``TAIR:`` in the with/from column).

This script collapses transcripts to genes, sets the DB column to ``SpudDB``
and drops duplicate rows. Feed the result to prepare_go_eval.py with::

    python scripts/prepare_go_eval.py --out-dir OUT \\
        --gaf-url file:///abs/path/spuddb_potato.gaf.gz \\
        --gaf-db SpudDB --gene-id-field db_object_id --variants all

Only the ``all`` variant is meaningful: there is no IBA, and every row is IEA.

>>> transcript_to_gene("Soltu.DM.09G009360.1")
'Soltu.DM.09G009360'
"""

from __future__ import annotations

import argparse
import gzip
import urllib.request
from pathlib import Path

SPUDDB_GO_URL = (
    "https://spuddb.uga.edu/data/dm_v61/"
    "DM_1-3_516_R44_potato.v6.1.working_models.go.txt.gz"
)


def transcript_to_gene(transcript_id: str) -> str:
    """Strip the trailing ``.N`` isoform suffix from a DM v6.1 transcript id."""
    return transcript_id.rsplit(".", 1)[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, required=True, help="Output .gaf.gz path")
    parser.add_argument("--source-url", default=SPUDDB_GO_URL)
    args = parser.parse_args()

    request = urllib.request.Request(args.source_url, headers={"User-Agent": "genesets-rs-eval/0.1"})
    with urllib.request.urlopen(request, timeout=300) as response:
        raw = gzip.decompress(response.read()).decode()

    seen: set[tuple[str, ...]] = set()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out, "wt") as out:
        out.write("!gaf-version: 2.2\n")
        out.write(f"!source: {args.source_url} (TAIR best-hit IEA), collapsed to gene\n")
        for line in raw.splitlines():
            if not line or line.startswith("!"):
                continue
            fields = line.split("\t")
            fields += [""] * (17 - len(fields))
            gene = transcript_to_gene(fields[1])
            fields[0], fields[1], fields[2] = "SpudDB", gene, gene
            key = (gene, fields[3], fields[4], fields[5], fields[6])
            if key in seen:
                continue
            seen.add(key)
            out.write("\t".join(fields[:17]) + "\n")
    genes = len({k[0] for k in seen})
    print(f"wrote {args.out}: {len(seen)} annotations on {genes} genes")


if __name__ == "__main__":
    main()
