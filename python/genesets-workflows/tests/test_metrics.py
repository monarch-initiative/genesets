from genesets_workflows.curation import metrics, model


def _assoc(go_id, category, seed):
    return model.TermAssociation(
        term=model.Term(id=go_id, label=go_id),
        seed_source=seed,
        category=category,
    )


def test_precision_recall_f1():
    assocs = [
        _assoc("GO:1", "core_process", "enrichment_recovered"),       # good, recovered, core
        _assoc("GO:2", "supporting_process", "enrichment_recovered"), # good, recovered
        _assoc("GO:3", "false_association", "enrichment_recovered"),  # bad, recovered
        _assoc("GO:4", "nonspecific", "enrichment_recovered"),        # bad, recovered
        _assoc("GO:5", "core_component", "curator_added"),            # core, missed
        _assoc("GO:6", None, "enrichment_recovered"),                 # un-adjudicated -> ignored
    ]
    m = metrics.score(assocs)
    assert m["recovered_adjudicated"] == 4
    assert m["precision"] == 0.5      # 2 good of 4 recovered+adjudicated
    assert m["core_total"] == 2       # GO:1 (recovered) + GO:5 (missed)
    assert m["recall"] == 0.5         # 1 of 2 core recovered
    assert round(m["f1"], 3) == 0.5


def test_score_empty_is_none():
    m = metrics.score([])
    assert m["precision"] is None
    assert m["recall"] is None
    assert m["f1"] is None
