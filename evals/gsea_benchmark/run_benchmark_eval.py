#!/usr/bin/env python3
"""Score `genesets-rs gsea` methods against the curated GO-interpretation gold.

The gold in `curation/genesets/*.yaml` interprets gene sets that carry no
scores, so this eval gives each set synthetic ones: members draw from
N(shift, 1) and every other annotated gene from N(0, 1). Each method then has
to recover the curator's CORE terms (`category` core_process/core_component)
from the ranking alone. Fisher on the set itself, the existing benchmark
pipeline (`evals/iba_vs_benchmark`), is the reference; Fisher on the top 5% of
the synthetic ranking is the thresholded analogue.

What this measures: how the GSEA statistics and the graph-aware algorithms
trade recall of curator-level biology against the number of terms reported.
What it does not: performance on real expression data, where scores are not
a shift on set membership. The gold is only read, never refit (see the
guardrail in `evals/iba_vs_benchmark/README.md`).

Usage (after `cargo build --release`)::

    python3 scripts/prepare_go_eval.py --out-dir /tmp/go_eval --variants all
    uv run --project python/genesets-workflows --extra curation \\
        python evals/gsea_benchmark/run_benchmark_eval.py --go-dir /tmp/go_eval \\
        --work-dir /tmp/gsea_benchmark --out evals/gsea_benchmark/RESULTS.md
"""

from __future__ import annotations

import argparse
import csv
import random
import statistics
import subprocess
import sys
import time
import urllib.error
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "genesets-workflows" / "src"))

from genesets_workflows.curation import model  # noqa: E402
from genesets_workflows.curation.ranking import read_ranking  # noqa: E402
from genesets_workflows.sources.mygeneset import (  # noqa: E402
    normalize_gene_values,
    request_json,
)

CORE = {"core_process", "core_component"}
MIN_SIZE = 15
MAX_SIZE = 500
ALPHA = 0.05

GSEA_METHODS = [
    ("KS classic", ["--statistic", "ks", "--algorithm", "classic"]),
    ("KS elim", ["--statistic", "ks", "--algorithm", "elim"]),
    ("KS weight01", ["--statistic", "ks", "--algorithm", "weight01"]),
    ("Wilcoxon classic", ["--statistic", "wilcoxon", "--algorithm", "classic"]),
    ("Wilcoxon elim", ["--statistic", "wilcoxon", "--algorithm", "elim"]),
    ("Wilcoxon weight01", ["--statistic", "wilcoxon", "--algorithm", "weight01"]),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--go-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--binary", default=str(ROOT / "target" / "release" / "genesets-rs"))
    parser.add_argument("--genesets-dir", type=Path, default=ROOT / "curation" / "genesets")
    parser.add_argument(
        "--shifts",
        type=lambda value: [float(part) for part in value.split(",")],
        default=[1.5, 3.0],
        help="Comma-separated member score shifts, in SD; one table each.",
    )
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261017)
    return parser.parse_args()


def load_gold(genesets_dir: Path) -> dict[str, dict]:
    """gene_set_name -> source id, CORE term ids, and insight-tagged ids."""
    gold = {}
    for path in sorted(genesets_dir.glob("*.yaml")):
        interp = model.load_interpretation(path)
        core = {a.term.id for a in interp.associations if a.category in CORE}
        insight: dict[str, set[str]] = defaultdict(set)
        for association in interp.associations:
            if association.insight:
                insight[str(association.insight)].add(association.term.id)
        seeded_by = sorted(
            {
                a.enrichment_stats.method
                for a in interp.associations
                if a.enrichment_stats is not None and a.enrichment_stats.method
            }
        )
        gold[interp.gene_set_name] = {
            "id": interp.gene_set_id,
            "core": core,
            "insight": insight,
            "ranking": (path.parent / interp.ranking.file) if interp.ranking else None,
            "seeded_by": seeded_by,
        }
    return gold


def fetch_msigdb(name: str, attempts: int = 4) -> list[str] | None:
    params = {"q": name, "species": "9606", "fields": "_id,source,genes.symbol", "size": 10}
    for attempt in range(attempts):
        try:
            data = request_json(params)
            break
        except urllib.error.HTTPError as error:
            # mygeneset.info returns occasional 502s; retry those, raise the rest.
            if error.code < 500 or attempt == attempts - 1:
                raise
            time.sleep(2 ** (attempt + 1))
    for hit in data.get("hits", []):
        if hit.get("_id") == name and hit.get("source") == "msigdb":
            return normalize_gene_values(hit.get("genes"), "symbol")
    return None


def build_queries(gold: dict, genesets_dir: Path, work_dir: Path) -> dict[str, list[str]]:
    """Members of every evaluable set, cached in work_dir/queries.gmt."""
    path = work_dir / "queries.gmt"
    if path.exists():
        return {
            line.split("\t")[0]: line.rstrip("\n").split("\t")[2:]
            for line in path.read_text().splitlines()
            if line.strip()
        }
    members: dict[str, list[str]] = {}
    for name, entry in gold.items():
        if str(entry["id"]).startswith("MSIGDB:"):
            genes = fetch_msigdb(name)
            if genes:
                members[name] = genes
    for line in (genesets_dir / "lit_members.gmt").read_text().splitlines():
        fields = line.rstrip("\n").split("\t")
        if fields and fields[0] in gold:
            members[fields[0]] = fields[2:]
    with open(path, "w") as handle:
        for name, genes in members.items():
            handle.write(name + "\tna\t" + "\t".join(genes) + "\n")
    return members


def read_term_sizes(go_dir: Path, universe: set[str]) -> dict[str, int]:
    """Propagated gene counts per term within the universe."""
    ancestors: dict[str, list[str]] = defaultdict(list)
    with open(go_dir / "closure.tsv") as handle:
        next(handle)
        for line in handle:
            child, ancestor = line.rstrip("\n").split("\t")
            ancestors[child].append(ancestor)
    genes: dict[str, set[str]] = defaultdict(set)
    with open(go_dir / "all" / "gene_terms.tsv") as handle:
        next(handle)
        for line in handle:
            gene, term = line.rstrip("\n").split("\t")
            if gene in universe:
                for up in ancestors.get(term, [term]):
                    genes[up].add(gene)
    return {term: len(members) for term, members in genes.items()}


def significant(path: Path, column: str) -> dict[str, set[str]]:
    hits: dict[str, set[str]] = defaultdict(set)
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if float(row[column]) <= ALPHA:
                hits[row["query_id"]].add(row["target_id"])
    return hits


def run(command: list[str]) -> float:
    start = time.perf_counter()
    subprocess.run(command, check=True)
    return time.perf_counter() - start


def write_rankings(members, universe, shift, rng, ranks_path, top_path) -> None:
    universe_set = set(universe)
    with open(ranks_path, "w") as ranks, open(top_path, "w") as top:
        ranks.write("ranking_id\tgene\tscore\n")
        for name, genes in members.items():
            member_set = set(genes) & universe_set
            scores = {gene: rng.gauss(shift if gene in member_set else 0.0, 1.0) for gene in universe}
            ordered = sorted(scores.items(), key=lambda item: -item[1])
            for gene, score in ordered:
                ranks.write(f"{name}\t{gene}\t{score:.9f}\n")
            top_genes = [gene for gene, _ in ordered[: len(ordered) // 20]]
            top.write(name + "\tna\t" + "\t".join(top_genes) + "\n")


def main() -> None:
    args = parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    gold = load_gold(args.genesets_dir)
    members = build_queries(gold, args.genesets_dir, args.work_dir)

    background = args.go_dir / "background_all_goa_symbols.txt"
    universe = [line for line in background.read_text().splitlines()[1:] if line]
    target_args = [
        "--annotations",
        str(args.go_dir / "all" / "gene_terms.tsv"),
        "--terms",
        str(args.go_dir / "terms.tsv"),
        "--closure",
        str(args.go_dir / "closure.tsv"),
    ]
    fisher_args = ["--background", str(background), "--min-overlap", "2", "--max-p-adjust", "0.05"]

    fisher_out = args.work_dir / "fisher_set.tsv"
    fisher_seconds = run(
        [args.binary, "matrix", *target_args, "--queries", str(args.work_dir / "queries.gmt")]
        + ["--query-format", "gmt", *fisher_args, "-o", str(fisher_out)]
    )
    fisher_hits = significant(fisher_out, "p_adjust_bonferroni")
    # As in score_method_vs_benchmark.py: sets with any Fisher hit.
    scored = [name for name in members if fisher_hits.get(name)]
    term_sizes = read_term_sizes(args.go_dir, set(universe))

    rng = random.Random(args.seed)
    tables = []
    for shift in args.shifts:
        ranks_path = args.work_dir / "rankings.tsv"
        top_path = args.work_dir / "top_queries.gmt"
        write_rankings(members, universe, shift, rng, ranks_path, top_path)
        hits = {"Fisher on the set": fisher_hits}
        seconds = {"Fisher on the set": fisher_seconds}
        top_out = args.work_dir / "fisher_top5.tsv"
        seconds["Fisher on top 5%"] = run(
            [args.binary, "matrix", *target_args, "--queries", str(top_path)]
            + ["--query-format", "gmt", *fisher_args, "-o", str(top_out)]
        )
        hits["Fisher on top 5%"] = significant(top_out, "p_adjust_bonferroni")
        for name, method_args in GSEA_METHODS:
            output = args.work_dir / (name.replace(" ", "_").lower() + ".tsv")
            seconds[name] = run(
                [args.binary, "gsea", *target_args, *method_args, "--ranks", str(ranks_path)]
                + ["--permutations", str(args.permutations), "--max-p-adjust", str(ALPHA)]
                + ["-o", str(output)]
            )
            hits[name] = significant(output, "p_adjust")
        tables.append((shift, table(gold, scored, hits, seconds, term_sizes)))
    real_rows = real_rankings(args, gold, target_args, term_sizes)
    write_report(args, gold, members, scored, term_sizes, tables, real_rows)


def real_rankings(args, gold, target_args, term_sizes) -> list[str]:
    """Score the GSEA methods on curated entries that carry real rankings."""
    ranked = sorted(name for name, entry in gold.items() if entry["ranking"] is not None)
    if not ranked:
        return []
    ranks_path = args.work_dir / "real_rankings.tsv"
    with open(ranks_path, "w") as handle:
        handle.write("ranking_id\tgene\tscore\n")
        for name in ranked:
            for gene, score in read_ranking(gold[name]["ranking"]).items():
                handle.write(f"{name}\t{gene}\t{score!r}\n")
    hits: dict[str, dict[str, set[str]]] = {}
    seconds: dict[str, float] = {}
    for name, method_args in GSEA_METHODS:
        output = args.work_dir / ("real_" + name.replace(" ", "_").lower() + ".tsv")
        seconds[name] = run(
            [args.binary, "gsea", *target_args, *method_args, "--ranks", str(ranks_path)]
            + ["--permutations", str(args.permutations), "--max-p-adjust", str(ALPHA)]
            + ["-o", str(output)]
        )
        hits[name] = significant(output, "p_adjust")
    seeded = "; ".join(
        f"{name} seeded by {', '.join(gold[name]['seeded_by']) or 'curators only'}"
        for name in ranked
    )
    return [
        "",
        "## Curated real rankings",
        "",
        f"Entries with a `ranking` block, scored on their real scores: {', '.join(ranked)}.",
        f"Their CORE terms were seeded from enrichment runs ({seeded}), so the row for the",
        "seeding method is not an independent test; the other rows are.",
        "",
        *table(gold, ranked, hits, seconds, term_sizes),
    ]


def table(gold, scored, hits, seconds, term_sizes) -> list[str]:
    def testable(term: str) -> bool:
        return MIN_SIZE <= term_sizes.get(term, 0) <= MAX_SIZE

    core_total = sum(len(gold[n]["core"]) for n in scored)
    testable_total = sum(sum(map(testable, gold[n]["core"])) for n in scored)
    lines = [
        "| Method | recall | recall (testable) | confirm | mechan | terms/set | seconds |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for method, method_hits in hits.items():
        recovered = sum(len(gold[n]["core"] & method_hits.get(n, set())) for n in scored)
        recovered_testable = sum(
            sum(1 for t in gold[n]["core"] & method_hits.get(n, set()) if testable(t))
            for n in scored
        )
        cells = [
            f"{recovered / core_total:.3f}" if core_total else "-",
            f"{recovered_testable / testable_total:.3f}" if testable_total else "-",
        ]
        for kind in ("confirmatory", "mechanistic"):
            total = sum(len(gold[n]["insight"].get(kind, ())) for n in scored)
            got = sum(
                len(gold[n]["insight"].get(kind, set()) & method_hits.get(n, set()))
                for n in scored
            )
            cells.append(f"{got / total:.3f} (n={total})" if total else "-")
        counts = [len(method_hits.get(n, ())) for n in scored]
        cells.append(f"{statistics.median(counts):g} / {statistics.fmean(counts):.0f}")
        cells.append(f"{seconds[method]:.0f}")
        lines.append(f"| {method} | " + " | ".join(cells) + " |")
    return lines


def write_report(args, gold, members, scored, term_sizes, tables, real_rows) -> None:
    core_total = sum(len(gold[n]["core"]) for n in scored)
    testable_total = sum(
        sum(1 for t in gold[n]["core"] if MIN_SIZE <= term_sizes.get(t, 0) <= MAX_SIZE)
        for n in scored
    )
    lines = [
        "# GSEA methods against the curated benchmark (synthetic scores)",
        "",
        "Generated by `evals/gsea_benchmark/run_benchmark_eval.py`; see its docstring",
        "for the design and caveats. Do not edit by hand.",
        "",
        f"- Sets with fetched membership: {len(members)}; scored (any Fisher-on-set "
        f"hit): {len(scored)}",
        f"- CORE terms in scored sets: {core_total}; of those {testable_total} have "
        f"{MIN_SIZE}-{MAX_SIZE} annotated genes, the GSEA size window",
        f"- Scores: members N(shift, 1), other annotated genes N(0, 1); seed "
        f"{args.seed}; KS permutations {args.permutations}",
        "- Significance: Fisher Bonferroni <= 0.05 (the benchmark's setting); KS GSEA "
        "FDR q <= 0.05; Wilcoxon BH <= 0.05",
        "- Fisher on the set ignores the scores, so its row repeats across shifts",
        "",
        "`recall` = CORE terms recovered / all CORE terms; `recall (testable)` = the",
        "same over CORE terms inside the size window; `confirm` / `mechan` = recall",
        "of insight-tagged terms; `terms/set` = median / mean significant terms per",
        "set.",
    ]
    for shift, rows in tables:
        lines += ["", f"## Member shift {shift} SD", "", *rows]
    lines += real_rows
    args.out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
