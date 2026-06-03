"""GO term filtering helpers for tutorial notebooks."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GoSubsetDescendants:
    subset: str
    subset_terms: frozenset[str]
    allowed_terms: frozenset[str]


def go_subset_descendants(go_dir: Path | str, subset: str = "goslim_generic") -> GoSubsetDescendants:
    """Return terms in a GO subset plus all descendants from the prepared closure table."""
    go_dir = Path(go_dir)
    subset_terms = _subset_terms(go_dir / "downloads/go-basic.obo", subset)
    allowed_terms = set(subset_terms)

    with (go_dir / "closure.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["ancestor"] in subset_terms:
                allowed_terms.add(row["child"])

    return GoSubsetDescendants(
        subset=subset,
        subset_terms=frozenset(subset_terms),
        allowed_terms=frozenset(allowed_terms),
    )


def _subset_terms(obo_path: Path, subset: str) -> set[str]:
    terms: set[str] = set()
    current_id: str | None = None
    current_subsets: set[str] = set()
    current_obsolete = False

    def commit() -> None:
        if current_id and not current_obsolete and subset in current_subsets:
            terms.add(current_id)

    with obo_path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                commit()
                current_id = None
                current_subsets = set()
                current_obsolete = False
                continue
            if line.startswith("["):
                commit()
                current_id = None
                current_subsets = set()
                current_obsolete = False
                continue
            if ": " not in line:
                continue

            key, value = line.split(": ", 1)
            if key == "id":
                current_id = value
            elif key == "subset":
                current_subsets.add(value)
            elif key == "is_obsolete" and value == "true":
                current_obsolete = True

    commit()
    return terms
