"""Precision/recall/F1 of enrichment against a curated gold standard.

Precision: of the recovered terms that a curator adjudicated, how many are
biologically good (core/supporting)?
Recall: of the curator-marked core terms, how many did enrichment recover?

``recovery_status`` decomposes the non-recovered core terms into the two kinds
of gap so they can be scored without weakening the biological category:
``annotation_gap`` (genes present but GO annotation too shallow -- a GO curation
target) and ``membership_gap`` (genes absent -- the gene set is incomplete).
``recall_supportable`` is the tool-fair recall that excludes membership gaps the
set could never support, alongside the biology-complete ``recall``.
"""

from __future__ import annotations

from typing import Iterable

GOOD_CATEGORIES = {"core_process", "core_component", "supporting_process"}
CORE_CATEGORIES = {"core_process", "core_component"}
RECOVERED = "enrichment_recovered"
MEMBERSHIP_GAP = "membership_gap"
ANNOTATION_GAP = "annotation_gap"


def effective_recovery_status(assoc) -> str | None:
    """The association's recovery_status, defaulting an enrichment-recovered
    term to ``annotation_supported`` when the field is unset (back-compat).

    >>> from types import SimpleNamespace as S
    >>> effective_recovery_status(S(seed_source="enrichment_recovered", recovery_status=None))
    'annotation_supported'
    >>> effective_recovery_status(S(seed_source="curator_added", recovery_status="membership_gap"))
    'membership_gap'
    """
    status = getattr(assoc, "recovery_status", None)
    if status:
        return status
    if assoc.seed_source == RECOVERED:
        return "annotation_supported"
    return None


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
    core = [a for a in associations if a.category in CORE_CATEGORIES]
    core_total = len(core)
    core_annotation_gap = sum(1 for a in core if effective_recovery_status(a) == ANNOTATION_GAP)
    core_membership_gap = sum(1 for a in core if effective_recovery_status(a) == MEMBERSHIP_GAP)
    recovered_core = sum(1 for a in core if a.seed_source == RECOVERED)
    # Tool-fair recall: exclude core terms the set can never support (genes absent).
    supportable = core_total - core_membership_gap
    recall_supportable = (recovered_core / supportable) if supportable > 0 else None
    return {
        "precision": p,
        "recall": r,
        "f1": f1(p, r),
        "recovered_adjudicated": recovered_adjudicated,
        "core_total": core_total,
        "core_annotation_gap": core_annotation_gap,
        "core_membership_gap": core_membership_gap,
        "recall_supportable": recall_supportable,
    }
