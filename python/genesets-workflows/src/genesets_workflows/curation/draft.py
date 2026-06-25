"""Assemble a draft GeneSetInterpretation from enrichment output."""

from __future__ import annotations

from typing import Callable, Iterable

from genesets_workflows.curation import model
from genesets_workflows.curation.enrichment import EnrichmentRow


def assemble_draft(
    *,
    gene_set_id: str,
    gene_set_name: str,
    collection: str,
    rows: Iterable[EnrichmentRow],
    contexts: list[model.BiologicalContext] | None = None,
    aspect_lookup: Callable[[str], str | None] | None = None,
    msigdb_release: str | None = None,
    n_genes: int | None = None,
) -> model.GeneSetInterpretation:
    """Build an un-adjudicated draft: every enrichment hit becomes a
    TermAssociation with seed_source=enrichment_recovered and no category."""
    associations: list[model.TermAssociation] = []
    for row in rows:
        aspect = aspect_lookup(row.target_id) if aspect_lookup else None
        associations.append(
            model.TermAssociation(
                term=model.Term(id=row.target_id, label=row.target_name),
                seed_source="enrichment_recovered",
                aspect=aspect,
                enrichment_stats=model.EnrichmentStats(
                    p_value=row.p_value,
                    p_adjust=row.p_adjust,
                    overlap=row.overlap,
                    query_size=row.query_size,
                    target_size=row.target_size,
                    overlap_genes=list(row.overlap_genes),
                ),
            )
        )
    return model.GeneSetInterpretation(
        gene_set_id=gene_set_id,
        gene_set_name=gene_set_name,
        collection=collection,
        msigdb_release=msigdb_release,
        n_genes=n_genes,
        contexts=contexts or [],
        associations=associations,
    )
