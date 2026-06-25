"""Precision/recall/F1 of enrichment against a curated gold standard.

Precision: of the recovered terms that a curator adjudicated, how many are
biologically good (core/supporting)?
Recall: of the curator-marked core terms, how many did enrichment recover?
"""

from __future__ import annotations

from typing import Iterable

GOOD_CATEGORIES = {"core_process", "core_component", "supporting_process"}
CORE_CATEGORIES = {"core_process", "core_component"}
RECOVERED = "enrichment_recovered"


def precision(associations: Iterable) -> float | None:
    """Fraction of recovered, adjudicated terms that are good.

    >>> from types import SimpleNamespace as S
    >>> a = S(seed_source="enrichment_recovered", category="core_process")
    >>> b = S(seed_source="enrichment_recovered", category="false_association")
    >>> precision([a, b])
    0.5
    """
    recovered = [a for a in associations if a.seed_source == RECOVERED and a.category is not None]
    if not recovered:
        return None
    good = [a for a in recovered if a.category in GOOD_CATEGORIES]
    return len(good) / len(recovered)


def recall(associations: Iterable) -> float | None:
    """Fraction of core terms that enrichment recovered.

    >>> from types import SimpleNamespace as S
    >>> a = S(seed_source="enrichment_recovered", category="core_process")
    >>> b = S(seed_source="curator_added", category="core_component")
    >>> recall([a, b])
    0.5
    """
    core = [a for a in associations if a.category in CORE_CATEGORIES]
    if not core:
        return None
    recovered = [a for a in core if a.seed_source == RECOVERED]
    return len(recovered) / len(core)


def f1(p: float | None, r: float | None) -> float | None:
    if p is None or r is None or (p + r) == 0:
        return None
    return 2 * p * r / (p + r)


def score(associations: Iterable) -> dict:
    associations = list(associations)
    p = precision(associations)
    r = recall(associations)
    recovered_adjudicated = sum(1 for a in associations if a.seed_source == RECOVERED and a.category is not None)
    core_total = sum(1 for a in associations if a.category in CORE_CATEGORIES)
    return {
        "precision": p,
        "recall": r,
        "f1": f1(p, r),
        "recovered_adjudicated": recovered_adjudicated,
        "core_total": core_total,
    }
