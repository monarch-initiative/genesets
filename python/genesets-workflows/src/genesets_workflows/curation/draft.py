"""Assemble a draft GeneSetInterpretation from enrichment output."""

from __future__ import annotations

from typing import Callable, Iterable

from genesets_workflows.curation import model
from genesets_workflows.curation.enrichment import EnrichmentRow, GseaRow


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
                    method="fisher",
                    p_adjust_method=row.p_adjust_method,
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


P_ADJUST_METHOD = {"gsea_ks": "gsea_fdr", "wilcoxon": "benjamini_hochberg"}


def assemble_ranked_draft(
    *,
    gene_set_id: str,
    gene_set_name: str,
    ranking: model.RankedSignature,
    rows: Iterable[GseaRow],
    method: str,
    algorithm: str = "classic",
    ks_weight: float = 1.0,
    contexts: list[model.BiologicalContext] | None = None,
    aspect_lookup: Callable[[str], str | None] | None = None,
) -> model.GeneSetInterpretation:
    """Build an un-adjudicated draft of a ranked signature from `gsea` rows.

    Each row becomes an association whose direction follows the sign of its
    NES. An unsigned ranking has no low pole, so its negative-NES rows are
    dropped: they say only that a term's genes are not associated.
    """
    associations: list[model.TermAssociation] = []
    for row in rows:
        if not ranking.signed and row.nes < 0:
            continue
        associations.append(
            model.TermAssociation(
                term=model.Term(id=row.target_id, label=row.target_name),
                seed_source="enrichment_recovered",
                aspect=aspect_lookup(row.target_id) if aspect_lookup else None,
                direction="UP" if row.nes >= 0 else "DN",
                enrichment_stats=model.EnrichmentStats(
                    method=method,
                    algorithm=algorithm,
                    ks_weight=ks_weight if method == "gsea_ks" else None,
                    p_adjust_method=P_ADJUST_METHOD[method],
                    p_value=row.p_value,
                    p_adjust=row.p_adjust,
                    target_size=row.target_size,
                    es=row.es,
                    nes=row.nes,
                    set_size=row.set_size,
                    leading_edge_size=row.leading_edge_size,
                    leading_edge_genes=list(row.leading_edge_genes),
                    peak_rank=row.peak_rank,
                ),
            )
        )
    return model.GeneSetInterpretation(
        gene_set_id=gene_set_id,
        gene_set_name=gene_set_name,
        ranking=ranking,
        contexts=contexts or [],
        associations=associations,
    )
