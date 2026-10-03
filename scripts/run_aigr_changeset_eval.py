#!/usr/bin/env python3
"""Run and score the AIGR change-set eval: does applying AI Gene Review
decisions to GOA improve enrichment of external gene sets?

Expects a directory prepared by ``prepare_aigr_changeset_eval.py`` plus:

- ``queries.gmt``          curated benchmark sets (``fetch_benchmark_queries.py``);
                           scored against the CORE / nonspecific terms in
                           ``curation/genesets/*.yaml``.
- ``breadth/<name>/queries.gmt`` (optional) uncurated external sets, e.g. MSigDB
                           CGP chemical-exposure or C6 oncogenic signatures
                           (``genesets_workflows.sources.mygeneset``). No gold:
                           scored on size, specificity and redundancy of the result.

Writes ``<eval-dir>/report/``: per-variant summaries, per-set deltas, the null
comparison, and the core terms gained/lost with the gene edits responsible.
The gold is never modified (see evals/aigr_changeset/README.md).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "python" / "genesets-workflows" / "src"))

import prepare_go_eval  # noqa: E402
from score_method_vs_benchmark import load_gold  # noqa: E402
from genesets_workflows.curation import model  # noqa: E402

MAIN_VARIANTS = ["all", "aigr_remove", "aigr_prune", "aigr_full", "aigr_core"]
CORE = {"core_process", "core_component"}
GENERAL_TERM_SIZE = 1000  # terms annotating >1000 genes are "general"


# --------------------------------------------------------------------- running
def run_genesets(
    binary: Path, eval_dir: Path, variant: str, queries: Path, out: Path, max_p: float | None,
    max_raw_p: float | None = None,
) -> None:
    if out.exists():  # results are cached per run dir; delete runs/ to recompute
        return
    config_path = out.with_suffix(".config.yaml")
    config = {
        "mode": "matrix",
        "ontology": {
            "terms": str(eval_dir / "terms.tsv"),
            "closure": str(eval_dir / "closure.tsv"),
            "annotations": str(eval_dir / variant / "gene_terms.tsv"),
        },
        "input": {"queries": str(queries), "query_format": "gmt"},
        "background": {"file": str(eval_dir / "background_all_goa_symbols.txt")},
        "min_overlap": 2,
        "correction": "bonferroni",
        "output": str(out),
    }
    if max_p is not None:
        config["max_p_adjust"] = max_p
    if max_raw_p is not None:
        config["max_p_value"] = max_raw_p
    out.parent.mkdir(parents=True, exist_ok=True)
    with config_path.open("w") as handle:
        prepare_go_eval.write_yaml_value(handle, config)
    subprocess.run([str(binary), "run", str(config_path)], check=True, stdout=subprocess.DEVNULL)


def read_results(path: Path, keep: dict[str, set[str]] | None = None) -> dict[str, dict[str, float]]:
    """query -> {term: raw p}. With ``keep``, only rows for those (query, term)."""
    out: dict[str, dict[str, float]] = defaultdict(dict)
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            q, t = sys.intern(row["query_id"]), sys.intern(row["target_id"])
            if keep is not None and t not in keep.get(q, ()):
                continue
            out[q][t] = float(row["p_value"])
    return out


# --------------------------------------------------------------------- helpers
def read_gmt(path: Path) -> dict[str, set[str]]:
    sets = {}
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) > 2:
            sets[parts[0]] = set(parts[2:])
    return sets


def read_ids(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines()[1:] if line.strip()}


def read_closure(path: Path) -> dict[str, set[str]]:
    anc: dict[str, set[str]] = defaultdict(set)
    with path.open() as handle:
        next(handle)
        for line in handle:
            child, ancestor = line.rstrip("\n").split("\t")
            anc[child].add(ancestor)
    return anc


def read_pairs(path: Path) -> dict[str, set[str]]:
    by_gene: dict[str, set[str]] = defaultdict(set)
    with path.open() as handle:
        next(handle)
        for line in handle:
            gene, term = line.rstrip("\n").split("\t")
            by_gene[gene].add(term)
    return by_gene


def propagate(terms: set[str], anc: dict[str, set[str]]) -> set[str]:
    out = set(terms)
    for t in terms:
        out |= anc.get(t, set())
    return out


def neglog10(p: float | None) -> float:
    if p is None or p >= 1:
        return 0.0
    return -math.log10(max(p, 1e-300))


def sign_test(deltas: list[float]) -> dict:
    up = sum(1 for d in deltas if d > 0)
    down = sum(1 for d in deltas if d < 0)
    n = up + down
    if n == 0:
        return {"up": 0, "down": 0, "same": len(deltas), "p_two_sided": 1.0}
    k = min(up, down)
    p = min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2**n)
    return {"up": up, "down": down, "same": len(deltas) - n, "p_two_sided": round(p, 4)}


def bootstrap_ci(deltas: list[float], reps: int = 2000, seed: int = 0) -> tuple[float, float]:
    if not deltas:
        return (0.0, 0.0)
    rng = random.Random(seed)
    means = sorted(statistics.fmean(rng.choices(deltas, k=len(deltas))) for _ in range(reps))
    return (round(means[int(0.025 * reps)], 4), round(means[int(0.975 * reps)], 4))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# --------------------------------------------------------------- curated arm
def load_nonspecific(genesets_dir: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for path in sorted(genesets_dir.glob("*.yaml")):
        interp = model.load_interpretation(path)
        for assoc in interp.associations:
            if assoc.category == "nonspecific":
                out[interp.gene_set_name].add(assoc.term.id)
    return out


def curated_set_metrics(name, gold, nonspec, sig, pvals) -> dict:
    core = set().union(*gold[name]["rec"].values()) if gold[name]["rec"] else set()
    hits = sig.get(name, set())
    ps = pvals.get(name, {})
    return {
        "core_total": len(core),
        "core_recovered": len(core & hits),
        "confirm_total": len(gold[name]["ins"].get("confirmatory", set())),
        "confirm_recovered": len(gold[name]["ins"].get("confirmatory", set()) & hits),
        "mechan_total": len(gold[name]["ins"].get("mechanistic", set())),
        "mechan_recovered": len(gold[name]["ins"].get("mechanistic", set()) & hits),
        "nonspecific_hits": len(nonspec.get(name, set()) & hits),
        "sig_terms": len(hits),
        "core_mean_neglog10p": round(statistics.fmean(neglog10(ps.get(t)) for t in core), 4) if core else 0.0,
    }


def curated_arm(args, eval_dir: Path, variants: list[str], report: Path, reviewed: set[str]) -> dict:
    queries = eval_dir / "queries.gmt"
    members = read_gmt(queries)
    gold = load_gold(args.genesets_dir)
    nonspec = load_nonspecific(args.genesets_dir)
    names = [n for n in members if n in gold]
    keep = {n: set().union(*gold[n]["rec"].values(), nonspec.get(n, set())) if gold[n]["rec"] else set() for n in names}

    per_variant: dict[str, dict[str, dict]] = {}
    sig_by_variant: dict[str, dict[str, set[str]]] = {}
    for v in variants:
        run_dir = eval_dir / "runs" / "curated" / v
        sig_path = run_dir / "results.tsv"
        gold_p_path = run_dir / "gold_term_pvalues.tsv"
        run_genesets(args.binary, eval_dir, v, queries, sig_path, args.max_p_adjust)
        if not gold_p_path.exists():
            # the unfiltered matrix is large: keep only the gold terms' p-values
            all_path = run_dir / "results_unfiltered.tsv"
            run_genesets(args.binary, eval_dir, v, queries, all_path, None)
            pvals = read_results(all_path, keep)
            write_tsv(gold_p_path, [{"query_id": q, "target_id": t, "p_value": p}
                                    for q, ts in pvals.items() for t, p in ts.items()])
            all_path.unlink()
        sig = {q: set(ts) for q, ts in read_results(sig_path).items()}
        pvals = read_results(gold_p_path)
        sig_by_variant[v] = sig
        per_variant[v] = {n: curated_set_metrics(n, gold, nonspec, sig, pvals) for n in names}

    scored = [n for n in names if sig_by_variant["all"].get(n)]
    exposure = {n: len(members[n] & reviewed) / len(members[n]) for n in names}

    def totals(v: str, subset: list[str]) -> dict:
        m = per_variant[v]
        s = lambda k: sum(m[n][k] for n in subset)  # noqa: E731
        return {
            "variant": v,
            "sets": len(subset),
            "recall_core": round(s("core_recovered") / max(1, s("core_total")), 4),
            "core_recovered": s("core_recovered"),
            "core_total": s("core_total"),
            "recall_confirm": round(s("confirm_recovered") / max(1, s("confirm_total")), 4),
            "recall_mechan": round(s("mechan_recovered") / max(1, s("mechan_total")), 4),
            "nonspecific_hits": s("nonspecific_hits"),
            "sig_terms": s("sig_terms"),
            "core_mean_neglog10p": round(statistics.fmean(m[n]["core_mean_neglog10p"] for n in subset), 4) if subset else 0,
        }

    tiers = {
        "all_sets": scored,
        "exposure_ge_25pct": [n for n in scored if exposure[n] >= 0.25],
        "exposure_lt_10pct": [n for n in scored if exposure[n] < 0.10],
    }
    summary_rows = []
    for tier, subset in tiers.items():
        for v in variants:
            if v.startswith("null_"):
                continue
            summary_rows.append({"tier": tier, **totals(v, subset)})
    nulls = [v for v in variants if v.startswith("null_")]
    if nulls:
        for tier, subset in tiers.items():
            null_totals = [totals(v, subset) for v in nulls]
            row = {"tier": tier, "variant": f"null_prune (mean of {len(nulls)})", "sets": len(subset)}
            for k in ("recall_core", "core_recovered", "core_total", "recall_confirm", "recall_mechan",
                      "nonspecific_hits", "sig_terms", "core_mean_neglog10p"):
                row[k] = round(statistics.fmean(t[k] for t in null_totals), 4)
            summary_rows.append(row)
    write_tsv(report / "curated_summary.tsv", summary_rows)

    # paired per-set deltas vs all
    delta_rows, tests = [], {}
    for v in variants:
        if v == "all" or v.startswith("null_"):
            continue
        for metric in ("core_recovered", "core_mean_neglog10p", "nonspecific_hits", "sig_terms"):
            d = [per_variant[v][n][metric] - per_variant["all"][n][metric] for n in scored]
            tests[f"{v}:{metric}"] = {"mean_delta": round(statistics.fmean(d), 4), "ci95": bootstrap_ci(d), **sign_test(d)}
        for n in scored:
            delta_rows.append({
                "variant": v, "gene_set": n, "exposure": round(exposure[n], 3),
                **{f"d_{k}": round(per_variant[v][n][k] - per_variant["all"][n][k], 4)
                   for k in ("core_recovered", "core_mean_neglog10p", "nonspecific_hits", "sig_terms")},
            })
    write_tsv(report / "curated_per_set_deltas.tsv", delta_rows)

    # aigr_prune vs random removal of the same rows
    null_cmp = {}
    if nulls:
        for tier, subset in tiers.items():
            prune = totals("aigr_prune", subset)
            null_cmp[tier] = {}
            for k in ("core_recovered", "core_mean_neglog10p", "nonspecific_hits", "sig_terms"):
                dist = [totals(v, subset)[k] for v in nulls]
                null_cmp[tier][k] = {
                    "aigr_prune": prune[k],
                    "null_mean": round(statistics.fmean(dist), 4),
                    "null_min": min(dist), "null_max": max(dist),
                    "frac_null_ge_prune": round(sum(x >= prune[k] for x in dist) / len(dist), 3),
                    "frac_null_le_prune": round(sum(x <= prune[k] for x in dist) / len(dist), 3),
                }

    # attribution: which CORE terms were gained/lost and which gene edits did it
    anc = read_closure(eval_dir / "closure.tsv")
    base_pairs = read_pairs(eval_dir / "all" / "gene_terms.tsv")
    attr_rows = []
    for v in ("aigr_prune", "aigr_full", "aigr_core"):
        if v not in variants:
            continue
        var_pairs = read_pairs(eval_dir / v / "gene_terms.tsv")
        for n in scored:
            core = set().union(*gold[n]["rec"].values()) if gold[n]["rec"] else set()
            for t in sorted(core):
                was, now = t in sig_by_variant["all"].get(n, ()), t in sig_by_variant[v].get(n, ())
                if was == now:
                    continue
                lost_genes, gained_genes = [], []
                for g in sorted(members[n]):
                    before = t in propagate(base_pairs.get(g, set()), anc)
                    after = t in propagate(var_pairs.get(g, set()), anc)
                    if before and not after:
                        lost_genes.append(g)
                    elif after and not before:
                        gained_genes.append(g)
                attr_rows.append({
                    "variant": v, "gene_set": n, "term_id": t, "change": "gained" if now else "lost",
                    "genes_losing_term": ",".join(lost_genes), "genes_gaining_term": ",".join(gained_genes),
                })
    write_tsv(report / "curated_core_term_changes.tsv", attr_rows)

    return {"summary": summary_rows, "paired_tests": tests, "null_comparison": null_cmp,
            "core_term_changes": len(attr_rows), "scored_sets": len(scored),
            "exposure_median": round(statistics.median(exposure[n] for n in scored), 3)}


# --------------------------------------------------------------- breadth arm
def breadth_arm(args, eval_dir: Path, variants: list[str], report: Path, reviewed: set[str]) -> dict:
    breadth_root = eval_dir / "breadth"
    if not breadth_root.exists():
        return {}
    anc = read_closure(eval_dir / "closure.tsv")
    desc_of: dict[str, set[str]] = defaultdict(set)
    for child, ancestors in anc.items():
        for a in ancestors:
            if a != child:
                desc_of[a].add(child)
    # fixed (baseline) term sizes so pruning cannot inflate specificity
    size: dict[str, int] = defaultdict(int)
    for _g, terms in read_pairs(eval_dir / "all" / "gene_terms.tsv").items():
        for t in propagate(terms, anc):
            size[t] += 1
    n_bg = len(read_ids(eval_dir / "background_all_goa_symbols.txt"))

    def ic(t: str) -> float:
        return -math.log2(max(size.get(t, 1), 1) / n_bg)

    # genesets-rs Bonferroni is matrix-wide (queries x terms), so a collection's
    # stringency would depend on how many sets it holds. Breadth collections
    # instead use a per-set Bonferroni over the baseline's non-empty GO terms,
    # identical for every collection and variant.
    n_terms = sum(1 for count in size.values() if count > 0)
    raw_p_cutoff = args.max_p_adjust / n_terms
    tiers = {"exposure_lt_10pct": (0.0, 0.10), "exposure_10_25pct": (0.10, 0.25), "exposure_ge_25pct": (0.25, 1.01)}

    def lost_profile_for(sig_all, sig_v, names):
        lost = [t for n in names for t in sig_all.get(n, set()) - sig_v.get(n, set())]
        return {
            "lost_hits": len(lost),
            "lost_mean_ic": round(statistics.fmean(ic(t) for t in lost), 4) if lost else None,
            "lost_general_frac": round(sum(size.get(t, 0) > GENERAL_TERM_SIZE for t in lost) / len(lost), 4) if lost else None,
        }

    out: dict = {"per_set_raw_p_cutoff": raw_p_cutoff, "n_terms": n_terms}
    rows, tier_rows = [], []
    pooled: dict[str, dict[str, dict]] = defaultdict(dict)  # variant -> "coll/set" -> metrics
    for collection in sorted(p.name for p in breadth_root.iterdir() if (p / "queries.gmt").exists()):
        queries = breadth_root / collection / "queries.gmt"
        members = read_gmt(queries)
        exposure = {n: len(m & reviewed) / len(m) for n, m in members.items()}
        per_variant, sig_sets = {}, {}
        for v in variants:
            path = eval_dir / "runs" / f"breadth_{collection}" / v / "results_perset.tsv"
            run_genesets(args.binary, eval_dir, v, queries, path, None, raw_p_cutoff)
            sig = {q: set(ts) for q, ts in read_results(path).items()}
            m = {}
            for n in members:
                hits = sig.get(n, set())
                redundant = sum(1 for t in hits if desc_of.get(t, set()) & hits)
                m[n] = {
                    "sig_terms": len(hits),
                    "mean_ic": statistics.fmean(ic(t) for t in hits) if hits else 0.0,
                    "general_terms": sum(1 for t in hits if size.get(t, 0) > GENERAL_TERM_SIZE),
                    "redundant_frac": redundant / len(hits) if hits else 0.0,
                    "any_hit": 1 if hits else 0,
                }
            per_variant[v] = m
            sig_sets[v] = sig
        base = per_variant["all"]
        lost_profile = {v: lost_profile_for(sig_sets["all"], sig_sets[v], members) for v in variants}
        for v in variants:
            for n in members:
                pooled[v][f"{collection}/{n}"] = {
                    "sig": sig_sets[v].get(n, set()), "exposure": exposure[n],
                }
        with_hits = [n for n in members if base[n]["any_hit"]]
        for v in variants:
            if v.startswith("null_"):
                continue
            m = per_variant[v]
            row = {"collection": collection, "variant": v, "sets": len(members), "sets_with_hits": sum(m[n]["any_hit"] for n in members)}
            for k in ("sig_terms", "general_terms"):
                row[k] = sum(m[n][k] for n in members)
            for k in ("mean_ic", "redundant_frac"):
                row[k] = round(statistics.fmean(m[n][k] for n in with_hits), 4) if with_hits else 0
            if v != "all":
                d = [m[n]["sig_terms"] - base[n]["sig_terms"] for n in with_hits]
                d_ic = [m[n]["mean_ic"] - base[n]["mean_ic"] for n in with_hits if m[n]["any_hit"]]
                st = sign_test(d)
                row["sets_fewer_terms"], row["sets_more_terms"] = st["down"], st["up"]
                row["d_mean_ic"] = round(statistics.fmean(d_ic), 4) if d_ic else 0
                row["d_mean_ic_ci95"] = bootstrap_ci(d_ic)
                row.update(lost_profile[v])
            rows.append(row)
        nulls = [v for v in variants if v.startswith("null_")]
        if nulls:
            out[collection] = {
                "aigr_prune_sig_terms": sum(per_variant["aigr_prune"][n]["sig_terms"] for n in members),
                "null_sig_terms": [sum(per_variant[v][n]["sig_terms"] for n in members) for v in nulls],
                "aigr_prune_lost": lost_profile["aigr_prune"],
                "null_lost_mean_ic": [lost_profile[v]["lost_mean_ic"] for v in nulls],
                "null_lost_general_frac": [lost_profile[v]["lost_general_frac"] for v in nulls],
            }
        out.setdefault(collection, {})["exposure_median"] = round(statistics.median(exposure.values()), 3)
    write_tsv(report / "breadth_summary.tsv", rows)

    # Pooled over every breadth collection, split by how much of each set AIGR has reviewed.
    nulls = [v for v in variants if v.startswith("null_")]
    for tier, (lo, hi) in tiers.items():
        keys = [k for k, m in pooled["all"].items() if lo <= m["exposure"] < hi]
        sig = {v: {k: pooled[v][k]["sig"] for k in keys} for v in variants}
        for v in variants:
            if v.startswith("null_"):
                continue
            row = {"tier": tier, "variant": v, "sets": len(keys),
                   "sets_with_hits": sum(1 for k in keys if sig[v][k]),
                   "sig_terms": sum(len(sig[v][k]) for k in keys)}
            if v != "all":
                row.update(lost_profile_for(sig["all"], sig[v], keys))
            tier_rows.append(row)
        if nulls:
            profiles = [lost_profile_for(sig["all"], sig[v], keys) for v in nulls]
            sizes = [sum(len(sig[v][k]) for k in keys) for v in nulls]
            ics = [p["lost_mean_ic"] for p in profiles if p["lost_mean_ic"] is not None]
            gens = [p["lost_general_frac"] for p in profiles if p["lost_general_frac"] is not None]
            tier_rows.append({
                "tier": tier, "variant": f"null_prune (range of {len(nulls)})", "sets": len(keys),
                "sig_terms": f"{min(sizes)}-{max(sizes)}",
                "lost_mean_ic": f"{min(ics)}-{max(ics)}" if ics else None,
                "lost_general_frac": f"{min(gens)}-{max(gens)}" if gens else None,
            })
    write_tsv(report / "breadth_by_exposure.tsv", tier_rows)
    return {"summary": rows, "by_exposure": tier_rows, **out}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--eval-dir", required=True, type=Path)
    parser.add_argument("--genesets-dir", type=Path, default=ROOT / "curation" / "genesets")
    parser.add_argument("--binary", type=Path, default=ROOT / "target" / "release" / "genesets-rs")
    parser.add_argument("--max-p-adjust", type=float, default=0.05)
    args = parser.parse_args()

    eval_dir = args.eval_dir.resolve()
    meta = json.loads((eval_dir / "metadata.json").read_text())
    variants = [v for v in MAIN_VARIANTS if (eval_dir / v).exists()]
    variants += sorted(v for v in meta["variants"] if v.startswith("null_"))
    reviewed = read_ids(eval_dir / "reviewed_genes.txt")
    report = eval_dir / "report"
    report.mkdir(exist_ok=True)

    result = {
        "changeset": meta["changeset"],
        "matching": meta["matching"],
        "variants": meta["variants"],
        "curated": curated_arm(args, eval_dir, variants, report, reviewed),
        "breadth": breadth_arm(args, eval_dir, variants, report, reviewed),
    }
    with (report / "report.json").open("w") as handle:
        json.dump(result, handle, indent=2, default=str)
        handle.write("\n")
    print(f"report written to {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
