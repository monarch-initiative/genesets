#!/usr/bin/env python3
"""Planted-signal evaluation of `genesets-rs gsea` on the real GO graph.

The curated benchmark in `curation/` holds gene sets without scores, so it
cannot score ranked-list methods. This simulation plants a known signal in a
known GO Biological Process term and asks how each method credits it:

- `specific`: half the genes of a 20-60 gene term are shifted up. A good
  method finds the term, ranks it first, and does not drag its ancestors
  along.
- `general`: half the genes of a 100-300 gene term with at least two
  children are shifted up. Here the signal belongs to the general term; a
  method that pushes credit down to children loses it.
- `null`: no signal; counts how many terms reach significance by chance.

Every ranking scores all genes with a BP annotation, N(0, 1) noise plus the
planted shift. Methods: GSEA KS and Wilcoxon rank-sum, each with the classic,
elim and weight01 algorithms, and Fisher over-representation (ORA) on the top
5% of the ranking as the thresholded baseline. Significance is p_adjust <=
0.05: the GSEA FDR q-value for KS, Benjamini-Hochberg for Wilcoxon and ORA.
For elim and weight01 that is a heuristic, since their p-values are not
independent tests.

Inputs come from `scripts/prepare_go_eval.py --variants all`. Usage::

    python3 scripts/prepare_go_eval.py --out-dir /tmp/go_eval --variants all
    cargo build --release
    python3 evals/gsea_planted/run_planted_eval.py --go-dir /tmp/go_eval \\
        --work-dir /tmp/gsea_planted --out evals/gsea_planted/RESULTS.md
"""

from __future__ import annotations

import argparse
import csv
import random
import statistics
import subprocess
import time
from collections import defaultdict
from pathlib import Path

BP_ROOT = "GO:0008150"
MIN_SIZE = 15
MAX_SIZE = 500
TOP_FRACTION = 0.05
ALPHA = 0.05

METHODS = [
    ("KS classic", ["gsea", "--statistic", "ks", "--algorithm", "classic"]),
    ("KS elim", ["gsea", "--statistic", "ks", "--algorithm", "elim"]),
    ("KS weight01", ["gsea", "--statistic", "ks", "--algorithm", "weight01"]),
    ("Wilcoxon classic", ["gsea", "--statistic", "wilcoxon", "--algorithm", "classic"]),
    ("Wilcoxon elim", ["gsea", "--statistic", "wilcoxon", "--algorithm", "elim"]),
    ("Wilcoxon weight01", ["gsea", "--statistic", "wilcoxon", "--algorithm", "weight01"]),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--go-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--binary", default="target/release/genesets-rs")
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--null-replicates", type=int, default=50)
    parser.add_argument("--fraction", type=float, default=0.5)
    parser.add_argument("--specific-shift", type=float, default=2.0)
    parser.add_argument("--general-shift", type=float, default=1.0)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261017)
    return parser.parse_args()


def read_pairs(path: Path) -> list[tuple[str, str]]:
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        return [(row[0], row[1]) for row in reader]


class Ontology:
    """BP terms, their propagated genes, and direct parents."""

    def __init__(self, go_dir: Path, work_dir: Path):
        ancestors: dict[str, set[str]] = defaultdict(set)
        for child, ancestor in read_pairs(go_dir / "closure.tsv"):
            ancestors[child].add(ancestor)
        self.terms = {term for term, up in ancestors.items() if BP_ROOT in up}
        self.ancestors = {term: ancestors[term] - {term} for term in self.terms}

        annotations = [
            (gene, term)
            for gene, term in read_pairs(go_dir / "all" / "gene_terms.tsv")
            if term in self.terms
        ]
        self.annotations_path = work_dir / "bp_gene_terms.tsv"
        with open(self.annotations_path, "w") as handle:
            handle.write("gene_id\tterm_id\n")
            for gene, term in annotations:
                handle.write(f"{gene}\t{term}\n")

        genes: dict[str, set[str]] = defaultdict(set)
        for gene, term in annotations:
            for up in ancestors[term] & self.terms:
                genes[up].add(gene)
        self.genes = dict(genes)
        self.universe = sorted(set().union(*self.genes.values()))
        self.tested = {
            term for term, members in self.genes.items() if MIN_SIZE <= len(members) <= MAX_SIZE
        }

        self.descendants: dict[str, set[str]] = defaultdict(set)
        for term, up in self.ancestors.items():
            for ancestor in up:
                self.descendants[ancestor].add(term)
        self.parents = {
            term: {
                a
                for a in up
                if not any(a in self.ancestors.get(b, ()) for b in up if b != a)
            }
            for term, up in self.ancestors.items()
        }
        self.children: dict[str, set[str]] = defaultdict(set)
        for term, parents in self.parents.items():
            for parent in parents:
                self.children[parent].add(term)

    def size(self, term: str) -> int:
        return len(self.genes.get(term, ()))

    def specific_candidates(self) -> list[str]:
        """20-60 gene terms with a tested, strictly larger ancestor."""
        return sorted(
            term
            for term in self.tested
            if 20 <= self.size(term) <= 60
            and any(
                a in self.tested and self.size(a) > self.size(term) for a in self.ancestors[term]
            )
        )

    def general_candidates(self) -> list[str]:
        """100-300 gene terms with at least two tested children."""
        return sorted(
            term
            for term in self.tested
            if 100 <= self.size(term) <= 300
            and sum(1 for child in self.children[term] if child in self.tested) >= 2
        )


def build_rankings(ontology: Ontology, args: argparse.Namespace) -> list[dict]:
    rng = random.Random(args.seed)
    specific = rng.sample(ontology.specific_candidates(), args.replicates)
    general = rng.sample(ontology.general_candidates(), args.replicates)
    plan = (
        [("specific", term, args.specific_shift) for term in specific]
        + [("general", term, args.general_shift) for term in general]
        + [("null", None, 0.0)] * args.null_replicates
    )
    rankings = []
    for index, (scenario, term, shift) in enumerate(plan):
        scores = {gene: rng.gauss(0.0, 1.0) for gene in ontology.universe}
        if term is not None:
            members = sorted(ontology.genes[term])
            for gene in rng.sample(members, round(args.fraction * len(members))):
                scores[gene] += shift
        rankings.append(
            {"id": f"{scenario}_{index:03d}", "scenario": scenario, "term": term, "scores": scores}
        )
    return rankings


def write_inputs(rankings: list[dict], work_dir: Path) -> tuple[Path, Path]:
    ranks_path = work_dir / "rankings.tsv"
    queries_path = work_dir / "top_queries.gmt"
    with open(ranks_path, "w") as ranks, open(queries_path, "w") as queries:
        ranks.write("ranking_id\tgene\tscore\n")
        for ranking in rankings:
            ordered = sorted(ranking["scores"].items(), key=lambda item: -item[1])
            for gene, score in ordered:
                ranks.write(f"{ranking['id']}\t{gene}\t{score:.9f}\n")
            top = [gene for gene, _ in ordered[: round(TOP_FRACTION * len(ordered))]]
            queries.write(ranking["id"] + "\tna\t" + "\t".join(top) + "\n")
    return ranks_path, queries_path


def run(command: list[str]) -> float:
    start = time.perf_counter()
    subprocess.run(command, check=True)
    return time.perf_counter() - start


def read_gsea(path: Path) -> dict[str, dict[str, tuple[float, float, float]]]:
    """ranking -> term -> (p_value, p_adjust, |NES|)."""
    results: dict[str, dict[str, tuple[float, float, float]]] = defaultdict(dict)
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            results[row["query_id"]][row["target_id"]] = (
                float(row["p_value"]),
                float(row["p_adjust"]),
                abs(float(row["nes"])),
            )
    return results


def read_ora(path: Path, ontology: Ontology) -> dict[str, dict[str, tuple[float, float, float]]]:
    """Fisher rows restricted to tested terms, with BH over all tested terms."""
    raw: dict[str, dict[str, float]] = defaultdict(dict)
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["target_id"] in ontology.tested:
                raw[row["query_id"]][row["target_id"]] = float(row["p_value"])
    tested = len(ontology.tested)
    results = {}
    for query, p_values in raw.items():
        # Terms with no overlap are absent from the output and have p = 1.
        ordered = sorted(p_values.items(), key=lambda item: item[1])
        adjusted = {}
        running = 1.0
        for position in range(len(ordered) - 1, -1, -1):
            term, p_value = ordered[position]
            running = min(running, p_value * tested / (position + 1))
            adjusted[term] = (p_value, running, 0.0)
        results[query] = adjusted
    return results


def order_key(row: tuple[float, float, float]) -> tuple[float, float]:
    return (row[0], -row[2])


def summarize(
    rankings: list[dict],
    results: dict[str, dict[str, tuple[float, float, float]]],
    ontology: Ontology,
) -> dict[str, dict[str, float]]:
    by_scenario: dict[str, list[dict[str, float]]] = defaultdict(list)
    for ranking in rankings:
        rows = results.get(ranking["id"], {})
        significant = {term for term, (_, p_adjust, _) in rows.items() if p_adjust <= ALPHA}
        term = ranking["term"]
        if term is None:
            by_scenario["null"].append({"significant": len(significant)})
            continue
        # Order as the output does: p-value, then |NES|; remaining ties share
        # the mid-rank.
        key_term = order_key(rows.get(term, (1.0, 1.0, 0.0)))
        keys = [order_key(row) for row in rows.values()]
        better = sum(1 for key in keys if key < key_term)
        tied = sum(1 for key in keys if key == key_term)
        rank = better + (max(tied, 1) + 1) / 2
        ancestors = ontology.ancestors[term]
        descendants = ontology.descendants[term]
        by_scenario[ranking["scenario"]].append(
            {
                "found": float(term in significant),
                "rank": rank,
                "top1": float(rank == 1),
                "ancestors": len(significant & ancestors),
                "descendants": len(significant & descendants),
                "other": len(significant - ancestors - descendants - {term}),
                "significant": len(significant),
            }
        )
    summary = {}
    for scenario, records in by_scenario.items():
        summary[scenario] = {
            key: (
                statistics.median(r[key] for r in records)
                if key == "rank"
                else statistics.fmean(r[key] for r in records)
            )
            for key in records[0]
        }
    return summary


def main() -> None:
    args = parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    ontology = Ontology(args.go_dir, args.work_dir)
    rankings = build_rankings(ontology, args)
    ranks_path, queries_path = write_inputs(rankings, args.work_dir)
    target_args = [
        "--annotations",
        str(ontology.annotations_path),
        "--terms",
        str(args.go_dir / "terms.tsv"),
        "--closure",
        str(args.go_dir / "closure.tsv"),
    ]

    summaries = {}
    seconds = {}
    for name, command in METHODS:
        output = args.work_dir / (name.replace(" ", "_").lower() + ".tsv")
        seconds[name] = run(
            [args.binary, *command, *target_args]
            + ["--ranks", str(ranks_path), "--permutations", str(args.permutations)]
            + ["--min-size", str(MIN_SIZE), "--max-size", str(MAX_SIZE), "-o", str(output)]
        )
        summaries[name] = summarize(rankings, read_gsea(output), ontology)

    background = args.work_dir / "background.txt"
    background.write_text("gene_id\n" + "\n".join(ontology.universe) + "\n")
    ora_output = args.work_dir / "ora_top5.tsv"
    seconds["ORA top 5%"] = run(
        [args.binary, "matrix", *target_args]
        + ["--queries", str(queries_path), "--query-format", "gmt"]
        + ["--background", str(background), "--correction", "none", "-o", str(ora_output)]
    )
    summaries["ORA top 5%"] = summarize(rankings, read_ora(ora_output, ontology), ontology)

    write_report(args, ontology, summaries, seconds)


def write_report(args, ontology, summaries, seconds) -> None:
    lines = [
        "# Planted-signal GSEA evaluation",
        "",
        "Generated by `evals/gsea_planted/run_planted_eval.py`; see its docstring for the",
        "design. Do not edit by hand.",
        "",
        f"- Ranked genes (any BP annotation): {len(ontology.universe):,}",
        f"- Tested BP terms ({MIN_SIZE}-{MAX_SIZE} genes): {len(ontology.tested):,}",
        f"- Replicates: {args.replicates} specific, {args.replicates} general, "
        f"{args.null_replicates} null; seed {args.seed}",
        f"- Signal: {args.fraction:.0%} of the planted term's genes shifted by "
        f"{args.specific_shift} SD (specific) or {args.general_shift} SD (general)",
        f"- KS permutations: {args.permutations}; significance: p_adjust <= {ALPHA} "
        "(GSEA FDR for KS, BH for Wilcoxon and ORA)",
        "",
        "Columns: `found` = planted term significant; `rank` = median rank of the",
        "planted term by p-value then |NES| (remaining ties share the mid-rank);",
        "`top1` = planted term ranked first; `anc` / `desc` / `other` = mean",
        "significant ancestors, descendants, and unrelated terms; `null sig` = mean",
        "significant terms with no signal planted.",
        "",
        "| Method | specific found | rank | top1 | anc | desc | other | general found | rank | "
        "top1 | anc | desc | other | null sig | seconds |",
        "| --- |" + " ---: |" * 14,
    ]
    for name, summary in summaries.items():
        s, g, n = summary["specific"], summary["general"], summary["null"]
        cells = [
            f"{s['found']:.2f}",
            f"{s['rank']:g}",
            f"{s['top1']:.2f}",
            f"{s['ancestors']:.2f}",
            f"{s['descendants']:.2f}",
            f"{s['other']:.2f}",
            f"{g['found']:.2f}",
            f"{g['rank']:g}",
            f"{g['top1']:.2f}",
            f"{g['ancestors']:.2f}",
            f"{g['descendants']:.2f}",
            f"{g['other']:.2f}",
            f"{n['significant']:.2f}",
            f"{seconds[name]:.0f}",
        ]
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    args.out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
