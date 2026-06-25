import pytest

from genesets_workflows.curation import metrics, model


def _assoc(go_id, category, seed, recovery=None):
    return model.TermAssociation(
        term=model.Term(id=go_id, label=go_id),
        seed_source=seed,
        category=category,
        recovery_status=recovery,
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
    assert m["recall_supportable"] is None


def test_recovery_status_decomposition():
    # category encodes biology; recovery_status records the gap, orthogonally.
    assocs = [
        _assoc("GO:1", "core_process", "enrichment_recovered", "annotation_supported"),
        _assoc("GO:2", "core_process", "curator_added", "annotation_gap"),   # genes present, GO shallow
        _assoc("GO:3", "core_process", "curator_added", "membership_gap"),   # genes absent from set
    ]
    m = metrics.score(assocs)
    assert m["core_total"] == 3
    assert m["core_annotation_gap"] == 1
    assert m["core_membership_gap"] == 1
    # biology-complete recall: 1 of 3 core recovered
    assert m["recall"] == pytest.approx(1 / 3)
    # tool-fair recall excludes the membership_gap term the set can never support: 1 of 2
    assert m["recall_supportable"] == 0.5
