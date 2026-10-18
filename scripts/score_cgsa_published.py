#!/usr/bin/env python3
"""Score cGSA's *published* outputs against its own ground truth — no model, no API.

cGSA (Wang et al., Bioinformatics 2026, btag214) emits free-text pathway names,
so its authors evaluated it two ways: an automatic MedCPT cosine-similarity match
and a blinded expert annotation. The expert annotation assigns every ground-truth
pathway and every cGSA output pathway an integer *functional category*; a shared
category is a match. That makes the expert HIT/ACC reproducible here with nothing
but the released spreadsheet — no embedding model, no OpenAI key.

This reproduces the paper's Fig. 4a expert numbers (ACC 0.828 / HIT 0.873 /
harmonic 0.850) from ``evals/cgsa/data/cgsa_{ground_truth,predictions}.tsv``
(built by ``evals/cgsa/data/build_cgsa_data.py``).

Definitions, per result set:
  ACC  = cGSA pathways whose category matches a ground-truth category / all cGSA pathways
  HIT  = ground-truth pathways whose category matches a cGSA category / all ground-truth pathways
  harmonic = 2·ACC·HIT / (ACC + HIT)

Category 0 = "unrelated to any ground-truth functional category" (the expert's
reject bucket). Two accounting variants are reported:
  matched  category 0 excluded from both sides (categories only compared among real ones)
  strict   category-0 cGSA outputs kept in the ACC denominator as misses
The strict ACC is the honest precision including off-context outputs; `matched`
mirrors the paper's headline. Aggregates are macro-averages over sets.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "evals" / "cgsa" / "data"


def _load(path: Path, cat_col: int) -> dict[str, list[int | None]]:
    by_set: dict[str, list[int | None]] = defaultdict(list)
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)  # header
        for row in reader:
            raw = row[cat_col].strip()
            cat = int(raw) if raw.lstrip("-").isdigit() else None
            by_set[row[0]].append(cat)
    return by_set


def _acc_hit(gt: list[int | None], pred: list[int | None], strict: bool) -> tuple[float, float] | None:
    gt_real = [c for c in gt if c is not None and c != 0]
    pred_kept = [c for c in pred if c is not None and (strict or c != 0)]
    if not gt_real or not pred_kept:
        return None
    gt_cats = set(gt_real)
    pred_cats = {c for c in pred_kept if c != 0}
    acc = sum(1 for c in pred_kept if c in gt_cats) / len(pred_kept)
    hit = sum(1 for c in gt_real if c in pred_cats) / len(gt_real)
    return acc, hit


def score(gt_by_set: dict, pred_by_set: dict, strict: bool) -> dict:
    ids = sorted(set(gt_by_set) | set(pred_by_set))
    accs, hits, harms = [], [], []
    for sid in ids:
        res = _acc_hit(gt_by_set.get(sid, []), pred_by_set.get(sid, []), strict)
        if res is None:
            continue
        acc, hit = res
        accs.append(acc)
        hits.append(hit)
        harms.append(0.0 if acc + hit == 0 else 2 * acc * hit / (acc + hit))
    ma, mh = statistics.mean(accs), statistics.mean(hits)
    return {
        "sets": len(accs),
        "ACC": round(ma, 3),
        "HIT": round(mh, 3),
        "harmonic_of_means": round(2 * ma * mh / (ma + mh), 3) if ma + mh else 0.0,
        "mean_per_set_harmonic": round(statistics.mean(harms), 3),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--out", type=Path, default=None, help="Optional summary TSV.")
    args = parser.parse_args()

    gt = _load(args.data_dir / "cgsa_ground_truth.tsv", 2)
    pred = _load(args.data_dir / "cgsa_predictions.tsv", 2)

    n_pred = sum(len(v) for v in pred.values())
    n_unrelated = sum(1 for v in pred.values() for c in v if c == 0)
    print(f"result sets: {len(set(gt) | set(pred))}   "
          f"ground-truth pathways: {sum(len(v) for v in gt.values())}   "
          f"cGSA pathways: {n_pred}   off-context (cat 0): {n_unrelated}\n")

    rows = []
    for variant, strict in (("matched", False), ("strict", True)):
        r = {"variant": variant, **score(gt, pred, strict)}
        rows.append(r)

    header = ["variant", "sets", "ACC", "HIT", "harmonic_of_means", "mean_per_set_harmonic"]
    width = {h: max(len(h), max(len(str(r[h])) for r in rows)) for h in header}
    print("  ".join(h.rjust(width[h]) for h in header))
    for r in rows:
        print("  ".join(str(r[h]).rjust(width[h]) for h in header))
    print("\npaper Fig. 4a (expert): ACC 0.828  HIT 0.873  harmonic 0.850")

    if args.out:
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
