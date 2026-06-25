"""Dataclasses mirroring the genesets_interpretation LinkML schema, with
YAML round-trip helpers. Optional (None / empty) fields are omitted on dump."""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass
class Term:
    id: str
    label: str


@dataclass
class EnrichmentStats:
    p_value: float | None = None
    p_adjust: float | None = None
    overlap: int | None = None
    query_size: int | None = None
    target_size: int | None = None
    overlap_genes: list[str] = field(default_factory=list)


@dataclass
class EvidenceItem:
    reference: str | None = None
    reference_title: str | None = None
    supports: str | None = None
    snippet: str | None = None
    explanation: str | None = None


@dataclass
class TermAssociation:
    term: Term
    seed_source: str
    aspect: str | None = None
    category: str | None = None
    confidence: str | None = None
    confidence_score: float | None = None
    specificity: str | None = None
    specificity_score: float | None = None
    enrichment_stats: EnrichmentStats | None = None
    evidence: list[EvidenceItem] = field(default_factory=list)
    curator_note: str | None = None


@dataclass
class BiologicalContext:
    term: Term
    context_type: str
    role_note: str | None = None


@dataclass
class GeneSetInterpretation:
    gene_set_id: str
    gene_set_name: str
    collection: str | None = None
    msigdb_release: str | None = None
    taxon: Term | None = None
    direction: str | None = None
    n_genes: int | None = None
    description: str | None = None
    contexts: list[BiologicalContext] = field(default_factory=list)
    associations: list[TermAssociation] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return _to_dict(self)


def _to_dict(obj: Any) -> Any:
    """Recursively convert dataclasses to plain dicts, dropping None and empties."""
    if is_dataclass(obj):
        result: dict[str, Any] = {}
        for f in fields(obj):
            value = _to_dict(getattr(obj, f.name))
            if value is None:
                continue
            if isinstance(value, (list, dict)) and len(value) == 0:
                continue
            result[f.name] = value
        return result
    if isinstance(obj, list):
        return [_to_dict(item) for item in obj]
    return obj


def _term(value: dict[str, Any] | None) -> Term | None:
    return None if value is None else Term(**value)


def _association(value: dict[str, Any]) -> TermAssociation:
    value = dict(value)
    value["term"] = Term(**value["term"])
    if "enrichment_stats" in value and value["enrichment_stats"] is not None:
        value["enrichment_stats"] = EnrichmentStats(**value["enrichment_stats"])
    if value.get("evidence"):
        value["evidence"] = [EvidenceItem(**e) for e in value["evidence"]]
    return TermAssociation(**value)


def _context(value: dict[str, Any]) -> BiologicalContext:
    value = dict(value)
    value["term"] = Term(**value["term"])
    return BiologicalContext(**value)


def from_dict(data: dict[str, Any]) -> GeneSetInterpretation:
    data = dict(data)
    data["taxon"] = _term(data.get("taxon"))
    data["contexts"] = [_context(c) for c in data.get("contexts", [])]
    data["associations"] = [_association(a) for a in data.get("associations", [])]
    return GeneSetInterpretation(**data)


def load_interpretation(path: str | Path) -> GeneSetInterpretation:
    data = yaml.safe_load(Path(path).read_text())
    return from_dict(data)


def dump_interpretation(interp: GeneSetInterpretation, path: str | Path) -> None:
    text = yaml.safe_dump(interp.to_dict(), sort_keys=False, allow_unicode=True)
    Path(path).write_text(text)
