from genesets_workflows.curation import draft, model
from genesets_workflows.curation.enrichment import EnrichmentRow


def _rows():
    return [
        EnrichmentRow("dc", "GO:0019882", "antigen processing and presentation", 12, 100, 340, 1e-9, 3e-7, ["HLA-DRA"]),
        EnrichmentRow("dc", "GO:0006412", "translation", 8, 100, 1500, 2e-2, 8e-1, ["RPL3"]),
    ]


def test_assemble_draft_basic():
    interp = draft.assemble_draft(
        gene_set_id="MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL",
        gene_set_name="HAY_BONE_MARROW_DENDRITIC_CELL",
        collection="C8",
        rows=_rows(),
        provenance=model.EnrichmentProvenance(tool="genesets-rs", go_version="snap-X"),
    )
    assert interp.gene_set_id == "MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"
    assert interp.curation_status == "draft"
    assert len(interp.associations) == 2
    a0 = interp.associations[0]
    assert a0.term.id == "GO:0019882"
    assert a0.seed_source == "enrichment_recovered"
    assert a0.category is None  # un-adjudicated
    assert a0.enrichment_stats.overlap == 12
    assert a0.enrichment_stats.p_adjust == 3e-7


def test_assemble_draft_applies_aspect_lookup():
    interp = draft.assemble_draft(
        gene_set_id="MSIGDB:X",
        gene_set_name="X",
        collection="C8",
        rows=_rows()[:1],
        provenance=model.EnrichmentProvenance(),
        aspect_lookup=lambda go_id: "biological_process",
    )
    assert interp.associations[0].aspect == "biological_process"
