#!/usr/bin/env python3
"""How much of the judged regulation signal is transcriptional, or carried by TFs?

For the pooled candidates that are regulation terms, split by judged category
(core/supporting vs the rest), report:

- the share whose regulated process is gene expression/transcription
  (an is_a/part_of descendant of GO:0010467);
- the share of contributing set genes (members annotated to the term) that
  have DNA-binding TF activity (GO:0003700 or a descendant) in GOA human;
- the most common positive regulation terms.

Contributing genes come from the judging packets (``batch_*.json``), which
``build_pool.py`` regenerates deterministically.

Usage: tf_share.py <go.obo> <goa_human.gaf.gz> <packet-dir> [judgments.yaml]
"""

from __future__ import annotations

import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
from prepare_regulation_variants import REG_RELS, parse_obo, regulation_terms  # noqa: E402

POS = {"core_process", "core_component", "supporting_process"}


def descendants(root: str, parents: dict[str, set[str]]) -> set[str]:
    children: dict[str, set[str]] = defaultdict(set)
    for c, ps in parents.items():
        for p in ps:
            children[p].add(c)
    out, stack = {root}, [root]
    while stack:
        for c in children.get(stack.pop(), ()):
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def main() -> int:
    obo, gaf, packets = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    jpath = Path(sys.argv[4]) if len(sys.argv) > 4 else Path(__file__).with_name("judgments.yaml")
    names, _, edges = parse_obo(obo)
    reg = regulation_terms(names, edges)
    isa_po: dict[str, set[str]] = defaultdict(set)
    for rel in ("is_a", "part_of"):
        for c, ps in edges[rel].items():
            isa_po[c] |= ps
    expr = descendants("GO:0010467", isa_po)
    tf_terms = descendants("GO:0003700", {c: ps for c, ps in edges["is_a"].items()})
    tf_genes = set()
    with gzip.open(gaf, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if line[0] != "!" and f[4] in tf_terms and "NOT" not in f[3]:
                tf_genes.add(f[2])
    contrib = {}
    for b in sorted(packets.glob("batch_*.json")):
        for pk in json.loads(b.read_text()):
            for c in pk["candidates"]:
                contrib[(pk["gene_set_name"], c["term_id"])] = c["set_genes_annotated_isa_partof"]

    def targets(t: str) -> set[str]:
        return set().union(*(edges[r].get(t, set()) for r in REG_RELS))

    judged = [j for j in yaml.safe_load(jpath.read_text()) if j["term_id"] in reg]
    for label, rows in (
        ("positive", [j for j in judged if j["category"] in POS]),
        ("other", [j for j in judged if j["category"] not in POS]),
    ):
        tx = sum(bool(targets(j["term_id"]) & expr) for j in rows)
        genes = [g for j in rows for g in contrib[(j["gene_set_name"], j["term_id"])]]
        tfg = sum(g in tf_genes for g in genes)
        print(
            f"{label}: {len(rows)} regulation terms; transcriptional {tx} ({tx / len(rows):.0%}); "
            f"contributing gene instances {len(genes)}, TF {tfg} ({tfg / len(genes):.0%})"
        )
        if label == "positive":
            print("  most common:", Counter(j["term_label"] for j in rows).most_common(10))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
