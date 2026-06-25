from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "curation" / "genesets" / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"


def test_example_loads_and_has_required_fields():
    data = yaml.safe_load(EXAMPLE.read_text())
    assert data["gene_set_id"] == "MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"
    assert data["contexts"], "must record at least one biological context"
    assert data["associations"], "must record at least one association"
    for assoc in data["associations"]:
        assert assoc["term"]["id"].startswith("GO:")
        assert assoc["term"]["label"]
        assert assoc["seed_source"] in {"enrichment_recovered", "curator_added"}


def test_example_categories_are_in_enum():
    data = yaml.safe_load(EXAMPLE.read_text())
    allowed = {
        "core_process",
        "core_component",
        "supporting_process",
        "marker_driven_plausible",
        "nonspecific",
        "false_association",
    }
    for assoc in data["associations"]:
        if "category" in assoc:
            assert assoc["category"] in allowed
