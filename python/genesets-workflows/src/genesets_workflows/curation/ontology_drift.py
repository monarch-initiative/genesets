"""Detect and migrate ontology drift in curated interpretations.

GO (and every other ontology the gold standard cites) changes under us: terms
are obsoleted, merged into other terms, or relabelled. Each change silently
degrades the benchmark (GOA never annotates to an obsolete term, so a core term
that went obsolete can never be recovered) or breaks validation outright. This
module classifies every term an interpretation asserts against the current
ontology and rewrites the changes that the ontology itself declares:

- ``label_changed``: the id is live but its primary label differs. Auto-fixed.
- ``merged``: the id is now an alternative id of a live term. Auto-fixed.
- ``obsolete_replaced``: obsolete with exactly one live ``replaced_by``
  (IAO:0100001). Auto-fixed: the ontology asserts the replacement.
- ``obsolete_consider``: obsolete with only ``consider`` suggestions. Needs a
  curator, who records the decision with an explicit remap.
- ``obsolete_no_replacement``: obsolete with nothing suggested. Needs a curator.
- ``missing``: the id is unknown to the ontology. Needs a curator.

An auto-fix is also held for review when its target is already asserted in the
same file, because applying it would duplicate an association.

Every applied change is appended to a migration log (TSV) recording the old and
new id/label, the reason, the ontology version and the tracker issue, so the
gold standard's history stays auditable. Migration never touches ``category``,
``recovery_status`` or ``insight``: it keeps the meaning a curator asserted and
only follows the ontology's own declarations (see the guardrail in
``evals/iba_vs_benchmark/README.md``).
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Protocol

import yaml

AUTO_FIXABLE = frozenset({"label_changed", "merged", "obsolete_replaced"})
MANUAL = frozenset({"obsolete_consider", "obsolete_no_replacement", "missing", "needs_review"})

REPLACED_BY = "IAO:0100001"
CONSIDER = "oio:consider"
TRACKER = "IAO:0000233"
OBSOLETE_LABEL_PREFIX = "obsolete "

LOG_FIELDS = [
    "date", "file", "old_id", "old_label", "new_id", "new_label", "kind",
    "ontology_version", "tracker", "note",
]

# NCBITaxon ids are stable and its adapter is large; it is skipped, as in the
# obsolescence gate.
SKIP_PREFIXES = frozenset({"NCBITaxon"})


class OntologyView(Protocol):
    """The slice of an ontology the classifier needs (an OAK adapter fits)."""

    def label(self, curie: str) -> str | None: ...

    def is_obsolete(self, curie: str) -> bool: ...

    def metadata(self, curie: str) -> dict[str, list]: ...

    def primary_for_alt_id(self, curie: str) -> str | None: ...

    def version(self) -> str: ...


@dataclass
class TermDrift:
    """One asserted term that no longer matches the current ontology."""

    file: str
    term_id: str
    label: str
    kind: str
    new_id: str | None = None
    new_label: str | None = None
    candidates: list[str] = field(default_factory=list)
    reason: str = ""
    tracker: str = ""
    ontology_version: str = ""

    @property
    def auto_fixable(self) -> bool:
        return self.kind in AUTO_FIXABLE

    def describe(self) -> str:
        """One human-readable line for reports and validator output.

        >>> TermDrift("a.yaml", "GO:1", "x", "obsolete_replaced", "GO:2", "y").describe()
        "GO:1 'x': obsolete_replaced -> GO:2 'y'"
        >>> TermDrift("a.yaml", "GO:1", "x", "obsolete_consider",
        ...           candidates=["GO:3"]).describe()
        "GO:1 'x': obsolete_consider; consider GO:3"
        """
        text = f"{self.term_id} {self.label!r}: {self.kind}"
        if self.new_id:
            text += f" -> {self.new_id} {self.new_label!r}"
        if self.candidates:
            text += f"; consider {', '.join(self.candidates)}"
        return text


def _live_label(view: OntologyView, curie: str) -> str | None:
    label = view.label(curie)
    if label is None or view.is_obsolete(curie):
        return None
    return label


def classify(term_id: str, label: str, view: OntologyView, file: str = "") -> TermDrift | None:
    """Classify one asserted (id, label) against the ontology; None if current.

    >>> view = FakeView(
    ...     labels={"GO:1": "obsolete old", "GO:2": "new", "GO:3": "kept",
    ...             "GO:4": "obsolete gone", "GO:5": "obsolete vague"},
    ...     obsolete={"GO:1", "GO:4", "GO:5"},
    ...     meta={"GO:1": {REPLACED_BY: ["GO:2"]}, "GO:5": {CONSIDER: ["GO:3"]}},
    ...     alt={"GO:9": "GO:3"})
    >>> classify("GO:3", "kept", view) is None
    True
    >>> classify("GO:3", "old name", view).kind
    'label_changed'
    >>> d = classify("GO:1", "old", view); (d.kind, d.new_id, d.new_label)
    ('obsolete_replaced', 'GO:2', 'new')
    >>> classify("GO:4", "gone", view).kind
    'obsolete_no_replacement'
    >>> d = classify("GO:5", "vague", view); (d.kind, d.candidates)
    ('obsolete_consider', ['GO:3'])
    >>> d = classify("GO:9", "kept", view); (d.kind, d.new_id)
    ('merged', 'GO:3')
    >>> classify("GO:404", "nothing", view).kind
    'missing'
    """
    version = view.version()
    current = view.label(term_id)
    if current is None:
        primary = view.primary_for_alt_id(term_id)
        primary_label = _live_label(view, primary) if primary else None
        if primary and primary_label:
            return TermDrift(file, term_id, label, "merged", primary, primary_label,
                             ontology_version=version)
        return TermDrift(file, term_id, label, "missing", ontology_version=version)

    if view.is_obsolete(term_id):
        meta = view.metadata(term_id)
        reason = " ".join(str(v) for v in meta.get("rdfs:comment", []))
        tracker = " ".join(str(v) for v in meta.get(TRACKER, []))
        replacements = [str(r) for r in meta.get(REPLACED_BY, [])]
        live = [(r, _live_label(view, r)) for r in replacements]
        live = [(r, lab) for r, lab in live if lab]
        candidates = [str(c) for c in meta.get(CONSIDER, [])]
        if len(replacements) == 1 and len(live) == 1:
            new_id, new_label = live[0]
            return TermDrift(file, term_id, label, "obsolete_replaced", new_id, new_label,
                             reason=reason, tracker=tracker, ontology_version=version)
        kind = "obsolete_consider" if candidates or replacements else "obsolete_no_replacement"
        return TermDrift(file, term_id, label, kind,
                         candidates=replacements + candidates, reason=reason,
                         tracker=tracker, ontology_version=version)

    if current != label:
        return TermDrift(file, term_id, label, "label_changed", term_id, current,
                         ontology_version=version)
    return None


class FakeView:
    """In-memory OntologyView for tests and doctests."""

    def __init__(self, labels, obsolete=(), meta=None, alt=None, version="test"):
        self._labels = dict(labels)
        self._obsolete = set(obsolete)
        self._meta = meta or {}
        self._alt = alt or {}
        self._version = version

    def label(self, curie):
        return self._labels.get(curie)

    def is_obsolete(self, curie):
        return curie in self._obsolete

    def metadata(self, curie):
        return self._meta.get(curie, {})

    def primary_for_alt_id(self, curie):
        return self._alt.get(curie)

    def version(self):
        return self._version


class OakView:
    """OntologyView over an OAK adapter (the same adapters the validators use)."""

    def __init__(self, adapter_string: str):
        from oaklib import get_adapter

        self.adapter_string = adapter_string
        self._adapter = get_adapter(adapter_string)
        self._obsolete = set(self._adapter.obsoletes())
        self._alt: dict[str, str] | None = None
        self._version: str | None = None

    def label(self, curie):
        return self._adapter.label(curie)

    def is_obsolete(self, curie):
        return curie in self._obsolete

    def metadata(self, curie):
        return {k: list(v) for k, v in self._adapter.entity_metadata_map(curie).items()}

    def _sql(self, query: str) -> list[tuple]:
        """Run SQL against a semsql-backed adapter; [] for other backends."""
        session = getattr(self._adapter, "session", None)
        if session is None:
            return []
        from sqlalchemy import text

        return list(session.execute(text(query)).fetchall())

    def primary_for_alt_id(self, curie):
        if self._alt is None:
            rows = self._sql(
                "SELECT subject, value FROM statements "
                "WHERE predicate = 'oio:hasAlternativeId'"
            )
            self._alt = {value: subject for subject, value in rows if value}
        return self._alt.get(curie)

    def version(self):
        if self._version is None:
            rows = self._sql(
                "SELECT value, object FROM statements "
                "WHERE predicate IN ('owl:versionInfo', 'owl:versionIRI')"
            )
            values = [v for v, _ in rows if v] + [o for _, o in rows if o]
            self._version = values[0] if values else self.adapter_string
        return self._version


def asserted_terms(path: Path) -> list[tuple[str, str]]:
    """(id, label) of every validated term asserted in an interpretation.

    Covers context and association terms, and the current-term mapping
    (``term``) of each published analysis's reported terms. A reported term's
    ``reported_id`` / ``reported_label`` are a historical record copied as
    printed and are deliberately not included (nor rewritten by ``apply``:
    ``rewrite_text`` only matches plain ``id:`` / ``label:`` lines).

    Read straight from the YAML so drift can be scanned even in files the
    pydantic model would reject.

    >>> import tempfile, textwrap
    >>> body = textwrap.dedent('''
    ...     associations:
    ...       - term: {id: GO:1, label: a}
    ...     published_analyses:
    ...       - reference: PMID:1
    ...         reported_terms:
    ...           - reported_id: GO:9
    ...             reported_label: retired name
    ...             term: {id: GO:2, label: b}
    ...           - reported_label: unmapped
    ... ''')
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = Path(d, "x.yaml"); _ = p.write_text(body)
    ...     asserted_terms(p)
    [('GO:1', 'a'), ('GO:2', 'b')]
    """
    data = yaml.safe_load(Path(path).read_text()) or {}
    items = [*(data.get("contexts") or []), *(data.get("associations") or [])]
    for analysis in data.get("published_analyses") or []:
        items.extend(analysis.get("reported_terms") or [])
    terms = []
    for item in items:
        term = item.get("term") or {}
        if term.get("id"):
            terms.append((str(term["id"]), str(term.get("label", ""))))
    return terms


def load_views(oak_config: Path, prefixes: Iterable[str]) -> dict[str, OntologyView]:
    adapters = (yaml.safe_load(Path(oak_config).read_text()) or {}).get(
        "ontology_adapters", {}
    )
    return {
        prefix: OakView(adapters[prefix])
        for prefix in sorted(set(prefixes))
        if prefix in adapters and prefix not in SKIP_PREFIXES
    }


def scan(paths: Iterable[Path], views: dict[str, OntologyView]) -> list[TermDrift]:
    """Classify every asserted term in ``paths``; return only the drifted ones.

    Auto-fixes whose target is already asserted in the same file are
    downgraded to ``needs_review``.

    >>> import tempfile, textwrap
    >>> view = FakeView({"GO:1": "obsolete a", "GO:2": "b"}, obsolete={"GO:1"},
    ...                 meta={"GO:1": {REPLACED_BY: ["GO:2"]}})
    >>> body = textwrap.dedent('''
    ...     associations:
    ...       - term: {id: GO:1, label: a}
    ...       - term: {id: GO:2, label: b}
    ... ''')
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = Path(d, "x.yaml"); _ = p.write_text(body)
    ...     [x.kind for x in scan([p], {"GO": view})]
    ['needs_review']
    """
    drifts: list[TermDrift] = []
    for path in paths:
        terms = asserted_terms(path)
        asserted_ids = {tid for tid, _ in terms}
        seen: set[tuple[str, str]] = set()
        for tid, label in terms:
            if (tid, label) in seen:
                continue
            seen.add((tid, label))
            view = views.get(tid.split(":", 1)[0])
            if view is None:
                continue
            drift = classify(tid, label, view, file=str(path))
            if drift is None:
                continue
            if drift.auto_fixable and drift.new_id != tid and drift.new_id in asserted_ids:
                drift.reason = (f"auto-fix target {drift.new_id} is already asserted in "
                                f"this file; merge the associations by hand. {drift.reason}")
                drift.candidates = [drift.new_id]
                drift.kind = "needs_review"
            drifts.append(drift)
    return drifts


def remap(drift: TermDrift, new_id: str, view: OntologyView) -> TermDrift:
    """Turn a manual drift into a fix using a curator-chosen replacement id.

    >>> view = FakeView({"GO:4": "obsolete gone", "GO:7": "chosen"}, obsolete={"GO:4"})
    >>> d = remap(classify("GO:4", "gone", view), "GO:7", view)
    >>> (d.kind, d.new_id, d.new_label)
    ('curator_remap', 'GO:7', 'chosen')
    >>> remap(classify("GO:4", "gone", view), "GO:4", view)
    Traceback (most recent call last):
    ...
    ValueError: remap target GO:4 is obsolete or unknown
    """
    label = _live_label(view, new_id)
    if label is None:
        raise ValueError(f"remap target {new_id} is obsolete or unknown")
    return TermDrift(drift.file, drift.term_id, drift.label, "curator_remap", new_id, label,
                     candidates=drift.candidates, reason=drift.reason,
                     tracker=drift.tracker, ontology_version=drift.ontology_version)


def _yaml_scalar(value: str) -> str:
    """Render a label as a YAML scalar, quoting only when plain style would misparse.

    >>> _yaml_scalar("cytoplasmic stress granule assembly")
    'cytoplasmic stress granule assembly'
    >>> _yaml_scalar("ion: binding")
    '"ion: binding"'
    >>> _yaml_scalar("yes")
    '"yes"'
    """
    if value and value == value.strip() and yaml.safe_load(value) == value and "#" not in value:
        return value
    return json.dumps(value)


_ID_LINE = re.compile(r"^(?P<indent>\s*)id:\s*(?P<id>\S+)\s*$")
_LABEL_LINE = re.compile(r"^(?P<indent>\s*)label:\s*(?P<label>.*?)\s*$")
_FLOW_TERM = re.compile(r"\{\s*id:\s*(?P<id>[^,}\s]+)\s*,\s*label:\s*(?P<label>[^}]*?)\s*\}")


def rewrite_text(text: str, old_id: str, old_label: str, new_id: str, new_label: str) -> tuple[str, int]:
    """Rewrite every ``id: old_id`` + ``label: old_label`` pair; return (text, count).

    Only the two lines of each matching term change, so comments, folding and
    key order elsewhere in a hand-curated file survive.

    >>> src = "  - term:\\n      id: GO:1\\n      label: old name\\n    aspect: biological_process\\n"
    >>> out, n = rewrite_text(src, "GO:1", "old name", "GO:2", "new name")
    >>> n
    1
    >>> print(out, end="")
      - term:
          id: GO:2
          label: new name
        aspect: biological_process
    >>> rewrite_text("- term: {id: GO:1, label: old name}\\n", "GO:1", "old name", "GO:2", "x: y")
    ('- term: {id: GO:2, label: "x: y"}\\n', 1)
    """
    lines = text.splitlines(keepends=True)
    count = 0
    for i, line in enumerate(lines):
        m = _ID_LINE.match(line.rstrip("\r\n"))
        if not m or m.group("id") != old_id:
            continue
        j = i + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        if j >= len(lines):
            continue
        lm = _LABEL_LINE.match(lines[j].rstrip("\r\n"))
        if not lm or lm.group("indent") != m.group("indent"):
            continue
        if yaml.safe_load(lm.group("label") or "''") != old_label:
            continue
        eol_i = lines[i][len(lines[i].rstrip("\r\n")):]
        eol_j = lines[j][len(lines[j].rstrip("\r\n")):]
        lines[i] = f"{m.group('indent')}id: {new_id}{eol_i}"
        lines[j] = f"{lm.group('indent')}label: {_yaml_scalar(new_label)}{eol_j}"
        count += 1
    text = "".join(lines)

    def flow(match: re.Match) -> str:
        nonlocal count
        if match.group("id") != old_id or yaml.safe_load(match.group("label") or "''") != old_label:
            return match.group(0)
        count += 1
        return f"{{id: {new_id}, label: {_yaml_scalar(new_label)}}}"

    return _FLOW_TERM.sub(flow, text), count


def apply(drifts: Iterable[TermDrift], log_path: Path | None, note: str = "") -> list[TermDrift]:
    """Rewrite each fixable drift in place, verify it, and log it; return the applied ones.

    A rewrite is verified by re-reading the file: the old (id, label) must be
    gone and the new one present, otherwise the change is reverted and an
    error raised (never a silent partial migration).
    """
    applied: list[TermDrift] = []
    rows = []
    for drift in drifts:
        if not drift.new_id or drift.kind in MANUAL:
            continue
        path = Path(drift.file)
        before = path.read_text()
        after, count = rewrite_text(before, drift.term_id, drift.label, drift.new_id, drift.new_label)
        if count == 0:
            raise RuntimeError(f"{path}: could not locate {drift.term_id} {drift.label!r} to rewrite")
        path.write_text(after)
        terms = asserted_terms(path)
        if (drift.term_id, drift.label) in terms or (drift.new_id, drift.new_label) not in terms:
            path.write_text(before)
            raise RuntimeError(f"{path}: rewrite of {drift.term_id} did not verify; reverted")
        applied.append(drift)
        rows.append({
            "date": dt.date.today().isoformat(),
            "file": path.name,
            "old_id": drift.term_id,
            "old_label": drift.label,
            "new_id": drift.new_id,
            "new_label": drift.new_label,
            "kind": drift.kind,
            "ontology_version": drift.ontology_version,
            "tracker": drift.tracker,
            "note": note or drift.reason,
        })
    if log_path is not None and rows:
        new_file = not Path(log_path).exists()
        with Path(log_path).open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=LOG_FIELDS, delimiter="\t",
                                    lineterminator="\n")
            if new_file:
                writer.writeheader()
            writer.writerows(rows)
    return applied
