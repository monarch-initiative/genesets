#!/usr/bin/env python3
"""Adapter to run cGSA on this repo's gene sets and collect its output.

cGSA has **no API and no hosted service**. It is a local Python app
(``python PathDis.py``) that requires an Azure OpenAI key hard-coded into three
of its source files, plus live STRING / Enrichr / OLS / PubMed access. So this
adapter cannot *execute* cGSA in this session — it does the two halves we can
own deterministically and documents the manual step in between:

  prepare   our GMT + sets TSV        ->  Data/<name>.xlsx in cGSA's input schema
  (manual)  clone ncbi-nlp/cGSA, set Azure keys, point PathDis at the xlsx,
            `python PathDis.py`        ->  cGSA_Results/<id>/*.json
  collect   cGSA_Results/<id>/         ->  predictions TSV (result_id, pathway, score)

The collected TSV matches ``evals/cgsa/data/cgsa_predictions.tsv`` (minus the
expert category, which only the authors' blinded annotation supplies), so a
fresh cGSA run scores through ``scripts/ground_cgsa_to_go.py`` and the repo's
benchmark machinery. See ``evals/cgsa/runner/README.md``.

cGSA's PathDis.py expects an .xlsx sheet named "Sheet1" with columns:
  PMC, Topic, Database, Context, Gene List, Species
where Database is one of its DB keys (GOBP/GOMF/GOCC/KEGG/MSigDB/Reactome/
Enrichment/Pathway/Disease/Family) and Species is 9606 or 10090.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

# cGSA's Database keys (PathDis.py DB dict). We pass the benchmark's own value
# through when it is already one of these; otherwise fall back to Enrichment.
CGSA_DB_KEYS = {
    "GOBP", "GOMF", "GOCC", "KEGG", "MSigDB", "Reactome",
    "Enrichment", "Pathway", "Disease", "Family",
}


def _read_gmt(path: Path) -> dict[str, list[str]]:
    members: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        members[parts[0]] = parts[2:]
    return members


def _read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def _norm_db(value: str) -> str:
    for token in (value or "").split(";"):
        token = token.strip()
        if token in CGSA_DB_KEYS:
            return token
    return "Enrichment"


def _taxon_to_species(taxon: str) -> str:
    return (taxon or "").replace("NCBITaxon:", "").strip() or "9606"


def prepare(args: argparse.Namespace) -> int:
    try:
        import openpyxl
    except ImportError:
        print("openpyxl is required: pip install openpyxl", file=sys.stderr)
        return 1

    members = _read_gmt(args.gmt)
    _, sets = _read_tsv(args.sets)
    by_id = {row["set_id"]: row for row in sets}

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet1"
    ws.append(["PMC", "Topic", "Database", "Context", "Gene List", "Species"])
    n = 0
    for set_id, genes in members.items():
        meta = by_id.get(set_id, {})
        topic = meta.get("keywords") or meta.get("research_objective") or set_id
        context = (meta.get("experimental_conditions") or meta.get("curated_context")
                   or meta.get("research_objective") or topic)
        database = _norm_db(meta.get("database", ""))
        species = _taxon_to_species(meta.get("taxon", ""))
        ws.append([set_id, topic, database, context, " ".join(genes), int(species)])
        n += 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(args.out)
    print(f"wrote {args.out}  ({n} gene sets in cGSA input schema)")
    print("\nNext (manual — needs the cGSA repo + an Azure OpenAI key):")
    print("  1. git clone https://github.com/ncbi-nlp/cGSA && cd cGSA")
    print("  2. set openai.api_key/api_base/api_version in evaluation.py, exploration.py, confidence.py")
    print(f"  3. cp {args.out} Data/  and point PathDis.py's read_excel at it")
    print("  4. python PathDis.py   ->  cGSA_Results/<id>/*.json")
    print("  5. run:  python run_cgsa.py collect --results-dir <cGSA>/cGSA_Results --out predictions.tsv")
    return 0


def _extract_pairs(obj) -> list[tuple[str, str]]:
    """Best-effort (pathway, score) extraction from cGSA's scored-results JSON.

    Summarized_Results_Refine_with_Scores.json is ``{id: scored_pathways}`` and
    cGSA has revised the inner shape over time, so we walk it defensively: dicts
    of pathway->score, lists of {pathway, score}-ish records, or nested combos.
    """
    pairs: list[tuple[str, str]] = []
    score_keys = ("score", "relevance", "relevance_score", "Relevance Score")
    name_keys = ("pathway", "name", "term", "Pathway")

    def walk(node):
        if isinstance(node, dict):
            name = next((str(node[k]) for k in name_keys if k in node), None)
            score = next((str(node[k]) for k in score_keys if k in node), None)
            if name is not None:
                pairs.append((name, score or ""))
                return
            for key, val in node.items():
                if isinstance(val, (int, float)) or (isinstance(val, str) and val.replace(".", "", 1).isdigit()):
                    pairs.append((str(key), str(val)))
                else:
                    walk(val)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(obj)
    return pairs


def collect(args: argparse.Namespace) -> int:
    rows: list[list[str]] = []
    result_dirs = sorted(p for p in args.results_dir.iterdir() if p.is_dir())
    for d in result_dirs:
        scored = d / "Summarized_Results_Refine_with_Scores.json"
        if not scored.exists():
            print(f"skip {d.name}: no scored-results json", file=sys.stderr)
            continue
        data = json.loads(scored.read_text())
        for pathway, score in _extract_pairs(data):
            rows.append([d.name, pathway, score])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["result_id", "cgsa_pathway", "relevance_score"])
        writer.writerows(rows)
    print(f"wrote {args.out}  ({len(rows)} pathways from {len(result_dirs)} result dirs)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prepare", help="Write a cGSA-schema .xlsx from a GMT + sets TSV.")
    p.add_argument("--gmt", type=Path, required=True, help="e.g. evals/cgsa/data/case_study_members.gmt")
    p.add_argument("--sets", type=Path, required=True, help="e.g. evals/cgsa/data/case_study_sets.tsv")
    p.add_argument("--out", type=Path, required=True, help="Output .xlsx (cGSA Data/ input).")
    p.set_defaults(func=prepare)

    c = sub.add_parser("collect", help="Flatten cGSA_Results/<id>/ JSON to a predictions TSV.")
    c.add_argument("--results-dir", type=Path, required=True, help="cGSA's cGSA_Results/ directory.")
    c.add_argument("--out", type=Path, required=True)
    c.set_defaults(func=collect)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
