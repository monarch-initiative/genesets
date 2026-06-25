"""Parse genesets-rs enrichment TSV and run the Rust enrichment kernel."""

from __future__ import annotations

import csv
import io
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class EnrichmentRow:
    query_id: str
    target_id: str
    target_name: str
    overlap: int
    query_size: int
    target_size: int
    p_value: float | None
    p_adjust: float | None
    overlap_genes: list[str]


def _to_float(value: str) -> float | None:
    value = value.strip()
    if value == "":
        return None
    return float(value)


def parse_enrichment_tsv(text: str) -> list[EnrichmentRow]:
    """Parse genesets-rs enrich/matrix TSV output into EnrichmentRow records.

    The adjusted-p column may be named ``p_adjust`` or ``p_adjust_bonferroni``.
    ``overlap_genes`` is optional and ``;``-separated.
    """
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    p_adjust_key = "p_adjust_bonferroni" if "p_adjust_bonferroni" in (reader.fieldnames or []) else "p_adjust"
    rows: list[EnrichmentRow] = []
    for record in reader:
        genes_raw = record.get("overlap_genes") or ""
        rows.append(
            EnrichmentRow(
                query_id=record["query_id"],
                target_id=record["target_id"],
                target_name=record.get("target_name", ""),
                overlap=int(record["overlap"]),
                query_size=int(record["query_size"]),
                target_size=int(record["target_size"]),
                p_value=_to_float(record.get("p_value", "")),
                p_adjust=_to_float(record.get(p_adjust_key, "")),
                overlap_genes=[g for g in genes_raw.split(";") if g],
            )
        )
    return rows


def run_enrichment(
    *,
    sample: Path,
    annotations: Path,
    terms: Path,
    closure: Path,
    background: Path,
    min_overlap: int = 3,
    max_p_adjust: float = 0.05,
    binary: str = "genesets-rs",
) -> str:
    """Run `genesets-rs enrich` and return its TSV stdout."""
    cmd = [
        binary,
        "enrich",
        "--annotations", str(annotations),
        "--terms", str(terms),
        "--closure", str(closure),
        "--sample", str(sample),
        "--background", str(background),
        "--min-overlap", str(min_overlap),
        "--max-p-adjust", str(max_p_adjust),
        "--overlap-genes",
    ]
    completed = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return completed.stdout
