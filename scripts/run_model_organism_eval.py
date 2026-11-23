#!/usr/bin/env python3
"""Score enrichment against the curated gold for non-human (model-organism) sets.

Each species' sets must be scored against that species' own GO annotations: a
mouse gene list against goa_human would silently lose most genes to symbol
mismatch. This runner, per taxon present in ``curation/genesets/*.yaml``:

  1. selects that taxon's sets from ``curation/genesets/model_organism_members.gmt``
     into ``<out>/<db>/queries.gmt``;
  2. runs ``scripts/prepare_go_eval.py`` with the species GAF (all-evidence variant,
     symbol gene ids, species-wide GOA background);
  3. runs ``genesets-rs run`` on the generated config;
  4. scores with ``scripts/score_method_vs_benchmark.py --taxon <taxon>``.

See evals/model_organism_benchmark/README.md.

    >>> sorted(SPECIES_GAF)[:2]
    ['NCBITaxon:10090', 'NCBITaxon:10116']
    >>> SPECIES_GAF["NCBITaxon:6239"]
    'wb'
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
GAF_BASE = "https://current.geneontology.org/annotations"

# NCBITaxon -> GO Consortium GAF basename (current.geneontology.org/annotations/<db>.gaf.gz).
SPECIES_GAF = {
    "NCBITaxon:10090": "mgi",
    "NCBITaxon:10116": "rgd",
    "NCBITaxon:7955": "zfin",
    "NCBITaxon:7227": "fb",
    "NCBITaxon:6239": "wb",
    "NCBITaxon:559292": "sgd",
    "NCBITaxon:4896": "pombase",
    "NCBITaxon:3702": "tair",
    "NCBITaxon:44689": "dictybase",
}


def taxa_by_set(genesets_dir: Path) -> dict[str, str]:
    """gene_set_name -> taxon CURIE for every non-human interpretation."""
    out: dict[str, str] = {}
    for path in sorted(genesets_dir.glob("*.yaml")):
        data = yaml.safe_load(path.read_text())
        taxon = (data.get("taxon") or {}).get("id")
        if taxon and taxon != "NCBITaxon:9606":
            out[data["gene_set_name"]] = taxon
    return out


def run(cmd: list[str]) -> None:
    print("$", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--genesets-dir", type=Path, default=ROOT / "curation" / "genesets")
    parser.add_argument("--members", type=Path, default=ROOT / "curation" / "genesets" / "model_organism_members.gmt")
    parser.add_argument("--genesets-rs", default=str(ROOT / "target" / "release" / "genesets-rs"))
    args = parser.parse_args()

    taxa = taxa_by_set(args.genesets_dir)
    lines_by_taxon: dict[str, list[str]] = {}
    for line in args.members.read_text().splitlines():
        name = line.split("\t", 1)[0]
        if name not in taxa:
            raise SystemExit(f"{args.members}: {name} has no non-human interpretation YAML")
        lines_by_taxon.setdefault(taxa[name], []).append(line)
    missing = sorted(set(taxa) - {ln.split("\t", 1)[0] for lns in lines_by_taxon.values() for ln in lns})
    if missing:
        raise SystemExit(f"no captured membership for: {', '.join(missing)}")

    for taxon, lines in sorted(lines_by_taxon.items()):
        db = SPECIES_GAF.get(taxon)
        if db is None:
            raise SystemExit(f"no GAF configured for {taxon}; add it to SPECIES_GAF")
        eval_dir = args.out_dir / db
        eval_dir.mkdir(parents=True, exist_ok=True)
        (eval_dir / "queries.gmt").write_text("\n".join(lines) + "\n")
        run([sys.executable, str(ROOT / "scripts" / "prepare_go_eval.py"), "--out-dir", str(eval_dir),
             "--gaf-url", f"{GAF_BASE}/{db}.gaf.gz", "--variants", "all"])
        run([args.genesets_rs, "run", str(eval_dir / "all" / "config.yaml")])
        run([sys.executable, str(ROOT / "scripts" / "score_method_vs_benchmark.py"), "--eval-dir", str(eval_dir),
             "--genesets-dir", str(args.genesets_dir), "--variants", "all", "--taxon", taxon,
             "--out", str(eval_dir / "benchmark_scores.tsv")])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
