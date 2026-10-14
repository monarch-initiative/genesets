#!/usr/bin/env python3
"""Append blinded pooled judgments to the curated gene-set interpretations.

Used for the regulation_closure recuration (see
``evals/regulation_closure/README.md``): candidate terms were pooled from the
top hits of every GO variant, stripped of which variant proposed them, and
judged on the biology of each set's membership. This script only *writes* those
judgments; it never reads any variant's results, so it cannot refit the gold.

Each judgment is appended to the end of the set's ``associations`` list (the
last key in every file), with ``seed_source: enrichment_recovered`` and a
curator_note prefixed by ``NOTE_TAG`` so the recuration can be audited or
reverted. Terms already present in a file are skipped.
"""

from __future__ import annotations

import argparse
import json
import textwrap
from pathlib import Path

import yaml

NOTE_TAG = "[pooled regulation-closure recuration, 2026-10]"
NAMESPACES = {"biological_process", "molecular_function", "cellular_component"}
CATEGORIES = {
    "core_process",
    "core_component",
    "supporting_process",
    "marker_driven_plausible",
    "nonspecific",
    "false_association",
}
RECOVERY = {"annotation_supported", "annotation_gap", "membership_gap"}


def go_terms(obo: Path) -> dict[str, tuple[str, str]]:
    out: dict[str, tuple[str, str]] = {}
    cur: dict[str, str] = {}
    for line in obo.read_text(encoding="utf-8").splitlines():
        if line.startswith("["):
            if cur.get("id") and not cur.get("obsolete"):
                out[cur["id"]] = (cur.get("name", ""), cur.get("namespace", ""))
            cur = {}
        elif line.startswith("id: "):
            cur["id"] = line[4:]
        elif line.startswith("name: "):
            cur["name"] = line[6:]
        elif line.startswith("namespace: "):
            cur["namespace"] = line[11:]
        elif line == "is_obsolete: true":
            cur["obsolete"] = "1"
    if cur.get("id") and not cur.get("obsolete"):
        out[cur["id"]] = (cur.get("name", ""), cur.get("namespace", ""))
    return out


def render(j: dict, label: str, aspect: str) -> str:
    q = json.dumps  # JSON strings are valid YAML scalars
    lines = [
        "  - term:",
        f"      id: {j['term_id']}",
        f"      label: {q(label)}",
        f"    aspect: {aspect}",
        f"    category: {j['category']}",
        f"    confidence: {j['confidence']}",
        f"    specificity: {j['specificity']}",
        "    seed_source: enrichment_recovered",
    ]
    if j.get("recovery_status"):
        lines.append(f"    recovery_status: {j['recovery_status']}")
    if j.get("insight"):
        lines.append(f"    insight: {j['insight']}")
    note = f"{NOTE_TAG} {j['curator_note'].strip()}"
    lines.append("    curator_note: >-")
    lines += ["      " + w for w in textwrap.wrap(note, 72)]
    return "\n".join(lines) + "\n"


def check(j: dict) -> None:
    assert j["category"] in CATEGORIES, j
    core_like = j["category"] in {"core_process", "core_component", "supporting_process"}
    if not core_like:
        j["insight"] = None
    if j["category"] in {"nonspecific", "false_association"}:
        j["recovery_status"] = None
    if j.get("recovery_status"):
        assert j["recovery_status"] in RECOVERY, j
    if j.get("insight"):
        assert j["insight"] in {"confirmatory", "mechanistic"}, j


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("judgments", nargs="+", type=Path)
    ap.add_argument("--obo", required=True, type=Path)
    ap.add_argument("--genesets-dir", type=Path, default=Path("curation/genesets"))
    args = ap.parse_args()

    terms = go_terms(args.obo)
    by_set: dict[str, list[dict]] = {}
    for path in args.judgments:
        for j in json.loads(path.read_text()):
            check(j)
            by_set.setdefault(j["gene_set_name"], []).append(j)

    files = {}
    for p in sorted(args.genesets_dir.glob("*.yaml")):
        d = yaml.safe_load(p.read_text())
        files[d["gene_set_name"]] = (p, {a["term"]["id"] for a in d.get("associations", [])})
    added = skipped = 0
    for name, js in sorted(by_set.items()):
        path, present = files[name]
        chunks = []
        for j in sorted(js, key=lambda x: x["term_id"]):
            if j["term_id"] in present or j["term_id"] not in terms:
                skipped += 1
                continue
            label, ns = terms[j["term_id"]]
            assert ns in NAMESPACES, j
            chunks.append(render(j, label, ns))
            present.add(j["term_id"])
        if chunks:
            text = path.read_text()
            if not text.endswith("\n"):
                text += "\n"
            path.write_text(text + "".join(chunks))
            added += len(chunks)
    print(f"added {added} associations to {len(by_set)} sets; skipped {skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
