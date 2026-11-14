"""Ontology-drift scan and migration, against an in-memory ontology."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

import pytest

from genesets_workflows.curation import model
from genesets_workflows.curation import ontology_drift as od

REPO = Path(__file__).resolve().parents[3]
EXAMPLE = REPO / "curation" / "genesets" / "LIT_FERROPTOSIS_SUPPRESSORS.yaml"


def _view_for(path: Path, **changes) -> od.FakeView:
    """A view where every asserted term is current, then apply `changes`."""
    labels = {tid: label for tid, label in od.asserted_terms(path)}
    labels.update(changes.get("labels", {}))
    return od.FakeView(labels, obsolete=changes.get("obsolete", ()),
                       meta=changes.get("meta", {}), alt=changes.get("alt", {}),
                       version="2099-01-01")


@pytest.fixture
def curated(tmp_path: Path) -> Path:
    path = tmp_path / EXAMPLE.name
    shutil.copy(EXAMPLE, path)
    return path


def test_current_corpus_file_has_no_drift_against_itself(curated):
    assert od.scan([curated], {"GO": _view_for(curated)}) == []


def test_replaced_by_is_migrated_verified_and_logged(curated, tmp_path):
    view = _view_for(
        curated,
        labels={"GO:0097707": "obsolete ferroptosis", "GO:9999999": "ferroptotic cell death"},
        obsolete={"GO:0097707"},
        meta={"GO:0097707": {od.REPLACED_BY: ["GO:9999999"], od.TRACKER: ["https://x/1"]}},
    )
    drifts = od.scan([curated], {"GO": view})
    assert [(d.kind, d.term_id, d.new_id) for d in drifts] == [
        ("obsolete_replaced", "GO:0097707", "GO:9999999")]

    before = curated.read_text()
    log = tmp_path / "migrations.tsv"
    od.apply(drifts, log)
    after = curated.read_text()

    # Only the id and label lines change; everything a curator wrote survives.
    changed = [(a, b) for a, b in zip(before.splitlines(), after.splitlines()) if a != b]
    assert changed == [("      id: GO:0097707", "      id: GO:9999999"),
                       ("      label: ferroptosis", "      label: ferroptotic cell death")]
    interp = model.load_interpretation(curated)
    migrated = [a for a in interp.associations if a.term.id == "GO:9999999"][0]
    assert migrated.category == "core_process" and migrated.insight == "confirmatory"

    rows = list(csv.DictReader(log.open(), delimiter="\t"))
    assert rows[0]["old_id"] == "GO:0097707" and rows[0]["new_id"] == "GO:9999999"
    assert rows[0]["ontology_version"] == "2099-01-01" and rows[0]["tracker"] == "https://x/1"


def test_label_change_and_merge_are_auto_fixable(curated):
    view = _view_for(curated, labels={"GO:0006979": "cellular response to oxidative stress"},
                     alt={"GO:0006749": "GO:0006751"})
    view._labels.pop("GO:0006749")
    view._labels["GO:0006751"] = "glutathione catabolic process"
    kinds = {d.term_id: d.kind for d in od.scan([curated], {"GO": view})}
    assert kinds == {"GO:0006979": "label_changed", "GO:0006749": "merged"}


def test_obsolete_without_replacement_needs_a_curator_remap(curated):
    view = _view_for(curated, labels={"GO:0019369": "obsolete arachidonate metabolic process",
                                      "GO:0001676": "long-chain fatty acid metabolic process"},
                     obsolete={"GO:0019369"})
    (drift,) = od.scan([curated], {"GO": view})
    assert drift.kind == "obsolete_no_replacement" and not drift.auto_fixable
    assert od.apply([drift], None) == []  # never applied without a decision
    od.apply([od.remap(drift, "GO:0001676", view)], None)
    assert ("GO:0001676", "long-chain fatty acid metabolic process") in od.asserted_terms(curated)


def test_rewrite_preserves_crlf_line_endings():
    text = "  - term:\r\n      id: GO:1\r\n      label: a\r\n"
    out, n = od.rewrite_text(text, "GO:1", "a", "GO:2", "b")
    assert n == 1 and out == "  - term:\r\n      id: GO:2\r\n      label: b\r\n"


def test_rewrite_ignores_same_id_with_a_different_label():
    text = "      id: GO:1\n      label: other\n"
    assert od.rewrite_text(text, "GO:1", "a", "GO:2", "b") == (text, 0)
