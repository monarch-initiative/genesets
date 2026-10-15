#!/usr/bin/env python3
"""Precision of each GO variant on the blinded pooled judgments.

For every variant, take its significant (set, term) hits that fall in the judged
pool and tally the curators' categories. ``pos`` = core + supporting,
``neg`` = nonspecific + false_association. Only judged hits count, so the rates
describe the contested region where the variants disagree, not all output.

Usage: score_pooled_precision.py <eval-dir> [judgments.yaml]
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import yaml

VARIANTS = ["current", "drop_regulation", "regulates_closure", "collapse_to_target"]
POS = {"core_process", "core_component", "supporting_process"}
NEG = {"nonspecific", "false_association"}


def main() -> int:
    eval_dir = Path(sys.argv[1])
    jpath = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).with_name("judgments.yaml")
    judged = {(j["gene_set_name"], j["term_id"]): j["category"] for j in yaml.safe_load(jpath.read_text())}
    reg = set((eval_dir / "regulation_terms.txt").read_text().split())
    print("variant\tjudged_hits\tcore\tpos\tplausible\tneg\tpos_rate\tneg_rate\treg_term_hits\treg_pos\treg_neg")
    for v in VARIANTS:
        with (eval_dir / v / "results.tsv").open() as fh:
            hits = {(r["query_id"], r["target_id"]) for r in csv.DictReader(fh, delimiter="\t")}
        c = Counter(judged[k] for k in hits if k in judged)
        rc = Counter(judged[k] for k in hits if k in judged and k[1] in reg)
        n = sum(c.values())
        pos, neg = sum(c[x] for x in POS), sum(c[x] for x in NEG)
        core = c["core_process"] + c["core_component"]
        print(
            f"{v}\t{n}\t{core}\t{pos}\t{c['marker_driven_plausible']}\t{neg}\t"
            f"{pos / n:.3f}\t{neg / n:.3f}\t{sum(rc.values())}\t"
            f"{sum(rc[x] for x in POS)}\t{sum(rc[x] for x in NEG)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
