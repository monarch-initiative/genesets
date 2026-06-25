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
