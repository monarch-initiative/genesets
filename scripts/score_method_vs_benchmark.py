#!/usr/bin/env python3
"""Score an enrichment method against the curated GO-interpretation benchmark.

Unlike the GO-impact analyses (which compare annotation variants to *each other*),
this scores each variant's enriched GO terms against external ground truth: the
curator-validated CORE terms in ``curation/genesets/*.yaml``, bucketed by
``recovery_status``.

Pipeline (see evals/iba_vs_benchmark/README.md):
  1. Build queries.gmt of the benchmark sets' member genes (mygeneset.info).
  2. scripts/prepare_go_eval.py --variants all,iba,iba_iea,no_contributes_to
  3. genesets-rs run <variant>/config.yaml  ->  <variant>/results.tsv
  4. this script

The PRIMARY metric is recall of the curated CORE biology (`category` in
core_process/core_component) — this is tool-independent of the gold's
`recovery_status` labels, so the eval never pressures us to refit the gold to a
method's output (which would make recall-vs-gold circular).

Metrics per variant:
  recall_core      recovered / total of all CORE terms (the biology).
  recall_confirm   recovered / total of `insight: confirmatory` terms (core +
                   supporting). Measures the annotation/method machinery on the
                   "known answer" signal.
  recall_mechan    recovered / total of `insight: mechanistic` terms. Measures
                   whether a method surfaces genuine INSIGHT, not just the
                   confirmatory plumbing. `mechanistic` is the small evaluable
                   denominator, so read it alongside that count.
  gap_recovered    count of annotation_gap CORE terms a method recovered — a
                   DISAGREEMENT between the curator's prediction and a tool run,
                   a REVIEW item adjudicated against GOA facts, NOT auto-applied
                   to the gold (the tail must not wag the dog).
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_SRC = ROOT / "python" / "genesets-workflows" / "src"
if PACKAGE_SRC.exists():
    sys.path.insert(0, str(PACKAGE_SRC))

from genesets_workflows.curation import model  # noqa: E402

CORE = {"core_process", "core_component"}


def load_gold(genesets_dir: Path) -> dict[str, dict[str, dict[str, set[str]]]]:
    """name -> {'rec': {recovery_status: CORE term ids},
               'ins': {insight: core+supporting term ids}}.

    `rec` drives the recovery_status metrics (CORE terms only). `ins` drives the
    confirmatory-vs-mechanistic split over every asserted-real term (all
    insight-tagged = core + supporting; nonspecific/false carry no insight tag).
    """
    gold: dict[str, dict[str, dict[str, set[str]]]] = {}
    for path in sorted(genesets_dir.glob("*.yaml")):
        interp = model.load_interpretation(path)
        rec: dict[str, set[str]] = defaultdict(set)
        ins: dict[str, set[str]] = defaultdict(set)
        for assoc in interp.associations:
            if assoc.category in CORE:
                rec[assoc.recovery_status or "unset"].add(assoc.term.id)
            if assoc.insight:
                ins[assoc.insight].add(assoc.term.id)
        gold[interp.gene_set_name] = {"rec": rec, "ins": ins}
    return gold


def load_hits(results_tsv: Path) -> dict[str, set[str]]:
    """query_id -> set of significant target term ids."""
    hits: dict[str, set[str]] = defaultdict(set)
    with results_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            hits[row["query_id"]].add(row["target_id"])
    return hits


def evaluable_sets(queries_gmt: Path, gold: dict) -> list[str]:
    names = [line.split("\t", 1)[0] for line in queries_gmt.read_text().splitlines() if line.strip()]
    return [n for n in names if n in gold]


def _tally(hits: dict[str, set[str]], gold: dict, names: list[str], axis: str, key: str) -> tuple[int, int]:
    """(total, recovered) over gold[name][axis][key] term ids for the given variant hits."""
    total = recov = 0
    for name in names:
        ids = gold[name][axis].get(key, set())
        total += len(ids)
        recov += len(ids & hits.get(name, set()))
    return total, recov


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--eval-dir", required=True, type=Path, help="Prepared dir: <variant>/results.tsv + queries.gmt")
    parser.add_argument("--genesets-dir", type=Path, default=ROOT / "curation" / "genesets")
    parser.add_argument("--variants", default="all,no_contributes_to,iba_iea,iba")
    parser.add_argument("--baseline", default="all", help="Variant used for unique-vs-baseline diagnostics.")
    parser.add_argument("--out", type=Path, default=None, help="Optional summary TSV path.")
    args = parser.parse_args()

    variants = [v.strip() for v in args.variants.split(",") if v.strip()]
    gold = load_gold(args.genesets_dir)
    names = evaluable_sets(args.eval_dir / "queries.gmt", gold)
    hits = {v: load_hits(args.eval_dir / v / "results.tsv") for v in variants}

    # only count sets that produced any enrichment under the baseline
    scored = [n for n in names if hits.get(args.baseline, {}).get(n)]

    def core_total(n):  # all CORE terms regardless of recovery_status
        return sum(len(b) for b in gold[n]["rec"].values())

    rows = []
    for v in variants:
        h = hits[v]
        core_t = sum(core_total(n) for n in scored)
        core_r = sum(len(set().union(*gold[n]["rec"].values()) & h.get(n, set()))
                     for n in scored if gold[n]["rec"])
        _, gap_r = _tally(h, gold, scored, "rec", "annotation_gap")
        conf_t, conf_r = _tally(h, gold, scored, "ins", "confirmatory")
        mech_t, mech_r = _tally(h, gold, scored, "ins", "mechanistic")
        rows.append({
            "variant": v,
            "sets": len(scored),
            "recall_core": round(core_r / core_t, 3) if core_t else 0.0,
            "confirmatory": conf_t,
            "recall_confirm": round(conf_r / conf_t, 3) if conf_t else 0.0,
            "mechanistic": mech_t,
            "recall_mechan": round(mech_r / mech_t, 3) if mech_t else 0.0,
            "gap_recovered": gap_r,
        })

    header = ["variant", "sets", "recall_core", "confirmatory", "recall_confirm",
              "mechanistic", "recall_mechan", "gap_recovered"]
    width = {h: max(len(h), max(len(str(r[h])) for r in rows)) for h in header}
    print("  ".join(h.rjust(width[h]) for h in header))
    for r in rows:
        print("  ".join(str(r[h]).rjust(width[h]) for h in header))

    if args.out:
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
