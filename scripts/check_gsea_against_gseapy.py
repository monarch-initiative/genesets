#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = ["gseapy>=1.1", "pandas"]
# ///
"""Cross-check ``genesets-rs gsea`` against gseapy's ``prerank``.

Builds a synthetic ranking (3,000 genes, Gaussian scores) and 60 gene sets, a
third biased toward the top of the ranking, a third toward the bottom, and a
third random. Runs both tools with weight 1 and gene-set permutation, then
reports how closely ES, NES and nominal p-values agree. ES is deterministic
and should match to output rounding; NES and p-values are Monte Carlo
estimates and should agree closely but not exactly.

Usage (from the repo root, after ``cargo build --release``)::

    uv run scripts/check_gsea_against_gseapy.py
"""

from __future__ import annotations

import argparse
import random
import subprocess
import tempfile
from pathlib import Path

import gseapy
import pandas as pd


def write_inputs(directory: Path, seed: int) -> None:
    rng = random.Random(seed)
    genes = [f"G{i:05d}" for i in range(3000)]
    scores = {gene: rng.gauss(0, 1) for gene in genes}
    ranked = sorted(genes, key=lambda gene: -scores[gene])
    with open(directory / "ranks.rnk", "w") as handle:
        for gene in genes:
            handle.write(f"{gene}\t{scores[gene]:.12f}\n")
    with open(directory / "sets.gmt", "w") as handle:
        for index in range(60):
            size = rng.randint(15, 300)
            pool = (ranked[:600], ranked[-600:], genes)[index % 3]
            members = rng.sample(pool, size // 2) + rng.sample(genes, size - size // 2)
            handle.write(f"SET{index:02d}\tna\t" + "\t".join(sorted(set(members))) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--binary", default="target/release/genesets-rs")
    parser.add_argument("--permutations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        write_inputs(directory, args.seed)
        ours_path = directory / "ours.tsv"
        subprocess.run(
            [
                args.binary,
                "gsea",
                "--target-sets",
                str(directory / "sets.gmt"),
                "--ranks",
                str(directory / "ranks.rnk"),
                "--permutations",
                str(args.permutations),
                "--output",
                str(ours_path),
            ],
            check=True,
        )
        ours = pd.read_csv(ours_path, sep="\t").set_index("target_id")
        ranks = pd.read_csv(directory / "ranks.rnk", sep="\t", header=None, index_col=0)[1]
        theirs = gseapy.prerank(
            rnk=ranks,
            gene_sets=str(directory / "sets.gmt"),
            permutation_num=args.permutations,
            weight=1,
            min_size=15,
            max_size=500,
            seed=1,
            outdir=None,
            no_plot=True,
            verbose=False,
        ).res2d
    theirs = theirs.set_index("Term")[["ES", "NES", "NOM p-val"]].astype(float)
    joined = theirs.join(ours, how="inner")

    significant_theirs = joined["NOM p-val"] < 0.05
    significant_ours = joined["p_value"] < 0.05
    print(f"gseapy {gseapy.__version__}")
    print(f"sets compared: {len(joined)} (gseapy {len(theirs)}, genesets-rs {len(ours)})")
    print(f"max |ES difference|: {(joined.ES - joined.es).abs().max():.3g}")
    print(f"NES Pearson r: {joined.NES.corr(joined.nes):.6f}")
    print(f"max |NES difference|: {(joined.NES - joined.nes).abs().max():.3g}")
    print(f"NES signs agree: {bool(((joined.NES > 0) == (joined.nes > 0)).all())}")
    print(
        "nominal p Spearman rho: "
        f"{joined['NOM p-val'].corr(joined.p_value, method='spearman'):.3f}"
    )
    print(
        "p < 0.05 (gseapy / genesets-rs / both): "
        f"{significant_theirs.sum()} / {significant_ours.sum()} / "
        f"{(significant_theirs & significant_ours).sum()}"
    )


if __name__ == "__main__":
    main()
