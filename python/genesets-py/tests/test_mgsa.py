import pytest

import genesets

# Genes g1-g4 sit under "child"; "parent" adds g5-g8; "other" is unrelated.
CLOSURE = [("child", "child"), ("child", "parent"), ("parent", "parent"), ("other", "other")]
ANNOTATIONS = (
    [(f"g{i}", t) for i in range(1, 5) for t in ("child", "parent")]
    + [(f"g{i}", "parent") for i in range(5, 9)]
    + [(f"g{i}", "other") for i in range(9, 21)]
)


def test_specific_term_explains_away_its_parent():
    result = genesets.mgsa(
        ANNOTATIONS, ["g1", "g2", "g3", "g4", "g9"], closure=CLOSURE, samples_per_term=5000
    )
    by_id = {term["id"]: term for term in result["terms"]}
    assert result["terms"][0]["id"] == "child"
    assert by_id["child"]["posterior"] > by_id["parent"]["posterior"]
    assert by_id["child"]["genes"] == ["g1", "g2", "g3", "g4"]
    assert result["study_size"] == 5
    assert result["background_size"] == 20
    assert result["stats"]["modeled_terms"] == 3


def test_seed_makes_runs_reproducible():
    first = genesets.mgsa(ANNOTATIONS, ["g1", "g2", "g9"], seed=7)
    second = genesets.mgsa(ANNOTATIONS, ["g1", "g2", "g9"], seed=7)
    assert first == second


def test_background_limits_the_universe():
    result = genesets.mgsa(ANNOTATIONS, ["g1", "g2"], background=[f"g{i}" for i in range(1, 9)])
    assert result["background_size"] == 8
    assert {term["id"] for term in result["terms"]} == {"child", "parent"}


def test_rejects_bad_priors():
    with pytest.raises(ValueError):
        genesets.mgsa(ANNOTATIONS, ["g1"], term_prob=2.0)
