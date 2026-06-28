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
  recall_core       recovered / total of all CORE terms (the biology). HEADLINE.
  recall_supported  recovered / total of the annotation_supported CORE subset
                    (descriptive: of terms the curator predicted recoverable)
  gap_recovered     count of annotation_gap CORE terms a method recovered. This
                    is a DISAGREEMENT signal between the curator's prediction and
                    a tool run -> a curator REVIEW item, adjudicated against GOA
                    facts. NOT auto-applied to the gold (tail must not wag dog).
  unique_vs_all     terms a variant recovers that the `all` variant does not
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


def load_gold(genesets_dir: Path) -> dict[str, dict[str, set[str]]]:
    """name -> recovery_status bucket -> set of CORE term ids."""
    gold: dict[str, dict[str, set[str]]] = {}
    for path in sorted(genesets_dir.glob("*.yaml")):
        interp = model.load_interpretation(path)
        buckets: dict[str, set[str]] = defaultdict(set)
        for assoc in interp.associations:
            if assoc.category in CORE:
                buckets[assoc.recovery_status or "unset"].add(assoc.term.id)
        gold[interp.gene_set_name] = buckets
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


def recovered(variant_hits: dict[str, set[str]], gold: dict, names: list[str], bucket: str) -> set[tuple[str, str]]:
    out: set[tuple[str, str]] = set()
    for name in names:
        hit = variant_hits.get(name, set())
        for tid in gold[name].get(bucket, set()):
            if tid in hit:
                out.add((name, tid))
    return out


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
    base_sup = recovered(hits[args.baseline], gold, scored, "annotation_supported")

    rows = []
    for v in variants:
        sup = recovered(hits[v], gold, scored, "annotation_supported")
        gap = recovered(hits[v], gold, scored, "annotation_gap")
        sup_t = sum(len(gold[n].get("annotation_supported", set())) for n in scored)
        # all CORE terms regardless of recovery_status = the curated biology
        core_t = sum(sum(len(b) for b in gold[n].values()) for n in scored)
        core_r = sum(len(set().union(*gold[n].values()) & hits[v].get(n, set()))
                     for n in scored if gold[n])
        rows.append({
            "variant": v,
            "sets": len(scored),
            "core": core_t,
            "core_recovered": core_r,
            "recall_core": round(core_r / core_t, 3) if core_t else 0.0,
            "recall_supported": round(len(sup) / sup_t, 3) if sup_t else 0.0,
            "gap_recovered": len(gap),
            "unique_vs_baseline": len(sup - base_sup),
        })

    header = ["variant", "sets", "core", "core_recovered", "recall_core",
              "recall_supported", "gap_recovered", "unique_vs_baseline"]
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
