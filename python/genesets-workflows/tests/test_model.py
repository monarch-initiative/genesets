from pathlib import Path

from genesets_workflows.curation import model

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "curation" / "genesets" / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"


def test_load_example():
    interp = model.load_interpretation(EXAMPLE)
    assert interp.gene_set_id == "MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"
    assert len(interp.associations) == 6
    first = interp.associations[0]
    assert first.term.id == "GO:0019882"
    assert first.category == "core_process"
    assert first.seed_source == "enrichment_recovered"


def test_roundtrip_drops_none(tmp_path):
    interp = model.load_interpretation(EXAMPLE)
    out = tmp_path / "rt.yaml"
    model.dump_interpretation(interp, out)
    text = out.read_text()
    assert "null" not in text  # None-valued optional fields are omitted
    reloaded = model.load_interpretation(out)
    assert reloaded.gene_set_id == interp.gene_set_id
    assert len(reloaded.associations) == len(interp.associations)


def test_build_minimal_interpretation():
    interp = model.GeneSetInterpretation(
        gene_set_id="MSIGDB:X",
        gene_set_name="X",
        associations=[
            model.TermAssociation(
                term=model.Term(id="GO:0006955", label="immune response"),
                seed_source="enrichment_recovered",
            )
        ],
    )
    assert interp.to_dict()["associations"][0]["term"]["id"] == "GO:0006955"


PUBLISHED = ROOT / "curation" / "genesets" / "LIT_MANALO_2005_HYPOXIA_HIF1_UP.yaml"


def test_published_analyses_load():
    interp = model.load_interpretation(PUBLISHED)
    roles = [a.analysis_role for a in interp.published_analyses]
    assert roles[0] == "original"
    assert "reanalysis" in roles
    ballouz = next(a for a in interp.published_analyses if a.tool == "ErmineJ 3.0")
    assert ballouz.membership_size == 202
    by_label = {r.reported_label: r for r in ballouz.reported_terms}
    hypoxia = by_label["cellular response to hypoxia"]
    assert hypoxia.term.id == "GO:0071456"
    assert hypoxia.p_adjust_corrected < 0.05
    # A term no longer in GO keeps its printed id and label but has no mapping.
    retired = by_label[
        "regulation of transcription from RNA polymerase II promoter in response to hypoxia"
    ]
    assert retired.reported_id == "GO:0061418"
    assert retired.term is None


def test_published_analyses_roundtrip(tmp_path):
    interp = model.load_interpretation(PUBLISHED)
    out = tmp_path / "rt.yaml"
    model.dump_interpretation(interp, out)
    reloaded = model.load_interpretation(out)
    assert reloaded.published_analyses == interp.published_analyses
    assert "null" not in out.read_text()
