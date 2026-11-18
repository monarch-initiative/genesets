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
    p_adjust_method: str | None = None


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
    bonferroni = "p_adjust_bonferroni" in (reader.fieldnames or [])
    p_adjust_key = "p_adjust_bonferroni" if bonferroni else "p_adjust"
    # genesets-rs names the column p_adjust when run with --correction none.
    p_adjust_method = "bonferroni" if bonferroni else "none"
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
                p_adjust_method=p_adjust_method,
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


@dataclass
class GseaRow:
    query_id: str
    target_id: str
    target_name: str
    es: float
    nes: float
    p_value: float
    p_adjust: float
    set_size: int
    target_size: int
    leading_edge_size: int
    peak_rank: int
    leading_edge_genes: list[str]


def parse_gsea_tsv(text: str) -> list[GseaRow]:
    """Parse `genesets-rs gsea` TSV output into GseaRow records.

    ``leading_edge_genes`` is optional (``--leading-edge-genes``) and
    ``;``-separated.

    >>> text = (
    ...     "query_id\\tquery_name\\ttarget_id\\ttarget_name\\tes\\tnes\\tp_value\\tp_adjust\\t"
    ...     "set_size\\teliminated\\ttarget_size\\tranked_size\\tleading_edge_size\\tpeak_rank\\t"
    ...     "leading_edge_genes\\tleading_edge_gene_names\\n"
    ...     "r\\t\\tGO:1\\tx\\t0.5\\t1.9\\t1e-3\\t2e-2\\t20\\t0\\t25\\t100\\t2\\t7\\tA;B\\t\\n"
    ... )
    >>> row = parse_gsea_tsv(text)[0]
    >>> (row.target_id, row.nes, row.set_size, row.leading_edge_genes)
    ('GO:1', 1.9, 20, ['A', 'B'])
    """
    rows: list[GseaRow] = []
    for record in csv.DictReader(io.StringIO(text), delimiter="\t"):
        genes_raw = record.get("leading_edge_genes") or ""
        rows.append(
            GseaRow(
                query_id=record["query_id"],
                target_id=record["target_id"],
                target_name=record.get("target_name", ""),
                es=float(record["es"]),
                nes=float(record["nes"]),
                p_value=float(record["p_value"]),
                p_adjust=float(record["p_adjust"]),
                set_size=int(record["set_size"]),
                target_size=int(record["target_size"]),
                leading_edge_size=int(record["leading_edge_size"]),
                peak_rank=int(record["peak_rank"]),
                leading_edge_genes=[g for g in genes_raw.split(";") if g],
            )
        )
    return rows
