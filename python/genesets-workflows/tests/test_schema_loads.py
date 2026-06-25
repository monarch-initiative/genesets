from pathlib import Path

from linkml_runtime import SchemaView

SCHEMA = Path(__file__).resolve().parents[3] / "curation" / "schema" / "genesets_interpretation.yaml"


def test_schema_has_expected_classes_and_enum():
    view = SchemaView(str(SCHEMA))
    classes = set(view.all_classes())
    assert {"GeneSetInterpretation", "TermAssociation", "Term", "BiologicalContext"} <= classes
    category = view.get_enum("AssociationCategoryEnum")
    values = set(category.permissible_values)
    assert values == {
        "core_process",
        "core_component",
        "supporting_process",
        "marker_driven_plausible",
        "nonspecific",
        "false_association",
    }


def test_tree_root_is_interpretation():
    view = SchemaView(str(SCHEMA))
    roots = [c.name for c in view.all_classes().values() if c.tree_root]
    assert roots == ["GeneSetInterpretation"]
