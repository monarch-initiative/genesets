"""Derive the GO aspect (BP/CC/MF) of a term from its ancestors."""

from __future__ import annotations

from typing import Callable, Iterable

ASPECT_ROOTS = {
    "GO:0008150": "biological_process",
    "GO:0005575": "cellular_component",
    "GO:0003674": "molecular_function",
}


def aspect_from_ancestors(go_id: str, ancestors: Iterable[str]) -> str | None:
    """Return the GO aspect implied by an ancestor set, or None if undetermined.

    >>> aspect_from_ancestors("GO:0019882", {"GO:0008150"})
    'biological_process'
    """
    ancestor_set = set(ancestors)
    for root, name in ASPECT_ROOTS.items():
        if root in ancestor_set:
            return name
    return None


def aspect_for_term(go_id: str, ancestors_lookup: Callable[[str], Iterable[str]]) -> str | None:
    """Derive aspect using a caller-supplied ancestors lookup (e.g. an OAK adapter)."""
    return aspect_from_ancestors(go_id, ancestors_lookup(go_id))


def oak_ancestors_lookup(adapter):  # pragma: no cover - exercised in integration runs
    """Build an ancestors lookup backed by an oaklib adapter (IS_A only)."""
    from oaklib.datamodels.vocabulary import IS_A

    def lookup(go_id: str):
        return set(adapter.ancestors(go_id, predicates=[IS_A]))

    return lookup
