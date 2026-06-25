from genesets_workflows.curation import aspect


def test_aspect_from_ancestors_bp():
    assert aspect.aspect_from_ancestors("GO:0019882", {"GO:0008150", "GO:0019882"}) == "biological_process"


def test_aspect_from_ancestors_cc():
    assert aspect.aspect_from_ancestors("GO:0042613", {"GO:0005575", "GO:0032991"}) == "cellular_component"


def test_aspect_from_ancestors_unknown():
    assert aspect.aspect_from_ancestors("GO:9999999", set()) is None


def test_aspect_for_term_uses_lookup():
    def fake_ancestors(go_id):
        return {"GO:0003674", go_id}

    assert aspect.aspect_for_term("GO:0004672", fake_ancestors) == "molecular_function"
