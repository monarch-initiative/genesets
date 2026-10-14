#!/usr/bin/env python3
"""Score the regulation-handling GO variants against the curated benchmark.

Companion to ``prepare_regulation_variants.py`` (see
``evals/regulation_closure/README.md``). Beyond the standard recall metrics of
``score_method_vs_benchmark.py`` it reports what changes when "regulation of X"
is handled differently:

- ``hits`` / ``reg_hit_frac``: significant terms, and the fraction that are
  regulation terms (output size and how much of it is regulation).
- ``recall_core``: strict recall of curated CORE terms.
- ``recall_core_reg`` / ``recall_core_part``: the same, split by whether the
  gold term is itself a regulation term or a participant term.
- ``recall_core_lenient``: a gold regulation-of-X term also counts as recovered
  when the variant hits X itself. This is the fair comparison for variants that
  cannot report regulation terms at all.
- ``neg_hits``: significant terms the curators judged ``nonspecific`` or
  ``false_association`` for that set (a precision proxy; most hits are
  uncurated, so this is a floor, not a false-positive rate).
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python" / "genesets-workflows" / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from genesets_workflows.curation import model  # noqa: E402
from prepare_regulation_variants import parse_obo, regulated_targets, regulation_terms  # noqa: E402
from score_method_vs_benchmark import evaluable_sets, load_hits  # noqa: E402

CORE = {"core_process", "core_component"}
NEGATIVE = {"nonspecific", "false_association"}
VARIANTS = "current,drop_regulation,regulates_closure,collapse_to_target"


def load_gold(genesets_dir: Path) -> dict[str, dict[str, set[str]]]:
    gold: dict[str, dict[str, set[str]]] = {}
    for path in sorted(genesets_dir.glob("*.yaml")):
        interp = model.load_interpretation(path)
        g: dict[str, set[str]] = defaultdict(set)
        for a in interp.associations:
            if a.category in CORE:
                g["core"].add(a.term.id)
            if a.category in NEGATIVE:
                g["neg"].add(a.term.id)
            if a.insight:
                g[a.insight].add(a.term.id)
        gold[interp.gene_set_name] = g
    return gold


def frac(n: int, d: int) -> float:
    return round(n / d, 3) if d else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eval-dir", required=True, type=Path)
    ap.add_argument("--obo", required=True, type=Path, help="Full go.obo used by the prepare step.")
    ap.add_argument("--genesets-dir", type=Path, default=ROOT / "curation" / "genesets")
    ap.add_argument("--variants", default=VARIANTS)
    ap.add_argument("--baseline", default="current")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--changes-out", type=Path, help="Per-set gold terms gained/lost vs baseline.")
    args = ap.parse_args()

    names_, _, edges = parse_obo(args.obo)
    reg = regulation_terms(names_, edges)
    variants = [v for v in args.variants.split(",") if v]
    gold = load_gold(args.genesets_dir)
    sets = evaluable_sets(args.eval_dir / "queries.gmt", gold)
    hits = {v: load_hits(args.eval_dir / v / "results.tsv") for v in variants}
    scored = [n for n in sets if hits[args.baseline].get(n)]

    rows = []
    for v in variants:
        h = hits[v]
        tot = defaultdict(int)
        rec = defaultdict(int)
        n_hits, n_reg_hits, per_set = 0, 0, []
        for n in scored:
            hs = h.get(n, set())
            n_hits += len(hs)
            n_reg_hits += len(hs & reg)
            per_set.append(len(hs))
            for t in gold[n]["core"]:
                kind = "reg" if t in reg else "part"
                tot[kind] += 1
                ok = t in hs
                rec[kind] += ok
                lenient = ok or (t in reg and bool(regulated_targets(t, edges, reg, names_) & hs))
                rec["lenient"] += lenient
            for key in ("confirmatory", "mechanistic"):
                tot[key] += len(gold[n][key])
                rec[key] += len(gold[n][key] & hs)
            rec["neg"] += len(gold[n]["neg"] & hs)
        core_t = tot["reg"] + tot["part"]
        rows.append(
            {
                "variant": v,
                "sets": len(scored),
                "hits": n_hits,
                "median_hits_per_set": statistics.median(per_set) if per_set else 0,
                "reg_hit_frac": frac(n_reg_hits, n_hits),
                "core_terms": core_t,
                "recall_core": frac(rec["reg"] + rec["part"], core_t),
                "core_reg_terms": tot["reg"],
                "recall_core_reg": frac(rec["reg"], tot["reg"]),
                "recall_core_part": frac(rec["part"], tot["part"]),
                "recall_core_lenient": frac(rec["lenient"], core_t),
                "recall_confirm": frac(rec["confirmatory"], tot["confirmatory"]),
                "recall_mechan": frac(rec["mechanistic"], tot["mechanistic"]),
                "neg_hits": rec["neg"],
            }
        )

    header = list(rows[0])
    width = {k: max(len(k), *(len(str(r[k])) for r in rows)) for k in header}
    print("  ".join(k.rjust(width[k]) for k in header))
    for r in rows:
        print("  ".join(str(r[k]).rjust(width[k]) for k in header))
    if args.out:
        with args.out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    if args.changes_out:
        base = hits[args.baseline]
        with args.changes_out.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["variant", "gene_set", "change", "term_id", "term_name", "is_regulation"])
            for v in variants:
                if v == args.baseline:
                    continue
                for n in scored:
                    core = gold[n]["core"]
                    for change, terms in (
                        ("gained", (hits[v].get(n, set()) - base.get(n, set())) & core),
                        ("lost", (base.get(n, set()) - hits[v].get(n, set())) & core),
                    ):
                        for t in sorted(terms):
                            w.writerow([v, n, change, t, names_.get(t, ""), "yes" if t in reg else "no"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
