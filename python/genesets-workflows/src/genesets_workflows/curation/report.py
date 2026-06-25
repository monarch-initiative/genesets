"""Build precision/recall/F1 report rows from curated interpretation files."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from genesets_workflows.curation import metrics, model

COLUMNS = [
    "gene_set_id",
    "gene_set_name",
    "collection",
    "n_associations",
    "recovered_adjudicated",
    "core_total",
    "precision",
    "recall",
    "f1",
]


def report_rows(paths: Iterable[Path]) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        interp = model.load_interpretation(path)
        scored = metrics.score(interp.associations)
        rows.append(
            {
                "gene_set_id": interp.gene_set_id,
                "gene_set_name": interp.gene_set_name,
                "collection": interp.collection or "",
                "n_associations": len(interp.associations),
                **scored,
            }
        )
    return rows


def write_report(paths: Iterable[Path], out: Path) -> Path:
    rows = report_rows(paths)
    with Path(out).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return Path(out)


def discover(curation_dir: Path) -> list[Path]:
    return sorted(Path(curation_dir).glob("*.yaml"))
