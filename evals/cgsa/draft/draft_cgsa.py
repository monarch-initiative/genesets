#!/usr/bin/env python3
"""Semi-automate cGSA-benchmark curation: do the deterministic gene-list work,
leave the judgment to a curator.

The split (per the repo's curation contract): everything about the *gene list* is
deterministic and done here; GO-term **category/insight** and **context** grounding
need a curator's judgment and are left blank/flagged.

Deterministic steps this tool performs:
  1. Pull the set's DEG membership + species + contexts from evals/cgsa/data/.
  2. Run the repo's real enrichment kernel (genesets-rs) against species-matched
     GOA (human goa_human, mouse MGI), so each candidate GO term arrives with its
     true overlap (carrier) genes and Bonferroni stats — no hand-asserted carriers.
  3. Assemble a draft interpretation: one association per significant term,
     seed_source=enrichment_recovered, recovery_status=annotation_supported,
     carriers in the curator_note, **category/insight/confidence left UNSET**.
  4. Deterministic screening: ground the benchmark's own ground-truth pathway
     names to GO (reusing the grounding cache) and report how many are recovered
     by enrichment — a GT-vs-membership *mismatch flag* (the PANC_IMMUNE failure
     mode) so unfit sets are caught before a curator spends time on them.
  5. Auto-suggest a disease context (top MONDO hit for the benchmark keywords) as
     a candidate for the curator to confirm — never silently accepted.

What the curator then does (judgment): set each association's category /
confidence / specificity / insight; confirm or replace the suggested context;
prune nonspecific/false terms; add any membership_gap terms; validate.

Usage:
  python evals/cgsa/draft/draft_cgsa.py --ids PMC9750880-1,PMC8983726-1 \
      --out-dir curation/cgsa/drafts --cache-dir /tmp/cgsa_annot
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "evals" / "cgsa" / "data"
PKG_SRC = ROOT / "python" / "genesets-workflows" / "src"
sys.path.insert(0, str(PKG_SRC))

from genesets_workflows.curation import draft as draft_mod  # noqa: E402
from genesets_workflows.curation import model  # noqa: E402
from genesets_workflows.curation.enrichment import parse_enrichment_tsv  # noqa: E402

SPECIES = {
    "9606": ("human", "https://current.geneontology.org/annotations/goa_human.gaf.gz"),
    "10090": ("mouse", "https://current.geneontology.org/annotations/mgi.gaf.gz"),
}


def load_sets() -> dict[str, dict]:
    with (DATA / "benchmark_sets.tsv").open() as h:
        return {r["set_id"]: r for r in csv.DictReader(h, delimiter="\t")}


def load_members() -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for line in (DATA / "benchmark_members.gmt").read_text().splitlines():
        if line.strip():
            p = line.split("\t")
            out[p[0]] = p[2:]
    return out


def load_gt_go() -> dict[str, set[str]]:
    """set_id -> set of GO ids its ground-truth names ground to (from the cache)."""
    cache_path = DATA / "go_grounding_cache.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    gt_names: dict[str, list[str]] = defaultdict(list)
    with (DATA / "benchmark_ground_truth.tsv").open() as h:
        for r in csv.DictReader(h, delimiter="\t"):
            gt_names[r["set_id"]].append(r["ground_truth_pathway"].strip())
    out: dict[str, set[str]] = {}
    for sid, names in gt_names.items():
        ids = set()
        for n in names:
            hit = cache.get(n)
            if hit and hit.get("id") and not hit.get("obsolete") and hit.get("exact"):
                ids.add(hit["id"])
        out[sid] = ids
    return out


def ensure_annotations(species_dir: Path, gaf_url: str, force: bool) -> None:
    if (species_dir / "all" / "gene_terms.tsv").exists() and not force:
        return
    cmd = [
        sys.executable, str(ROOT / "scripts" / "prepare_go_eval.py"),
        "--out-dir", str(species_dir), "--variants", "all", "--gaf-url", gaf_url,
    ]
    print(f"  preparing {species_dir.name} GOA annotations …", file=sys.stderr)
    subprocess.run(cmd, check=True)


def run_enrichment(species_dir: Path, queries: dict[str, list[str]], binary: str) -> str:
    gmt = species_dir / "queries.gmt"
    with gmt.open("w") as h:
        for sid, genes in queries.items():
            h.write("\t".join([sid, "cgsa_benchmark_degs", *genes]) + "\n")
    config = species_dir / "all" / "config.yaml"
    text = config.read_text()
    if "overlap_genes:" not in text:  # emit the carrier genes for grounding
        config.write_text(text.rstrip("\n") + "\noverlap_genes: true\n")
    subprocess.run([binary, "run", str(config)], check=True)
    return (species_dir / "all" / "results.tsv").read_text()


def aspect_lookup_from_terms(terms_tsv: Path):
    """Map GO id -> BP/CC/MF aspect from the prepared terms.tsv, if it carries one."""
    ns = {"biological_process": "biological_process",
          "cellular_component": "cellular_component",
          "molecular_function": "molecular_function",
          "P": "biological_process", "C": "cellular_component", "F": "molecular_function"}
    table: dict[str, str] = {}
    if terms_tsv.exists():
        with terms_tsv.open() as h:
            reader = csv.reader(h, delimiter="\t")
            header = next(reader, None)
            asp_i = None
            if header:
                for i, col in enumerate(header):
                    if col.lower() in ("aspect", "namespace", "go_aspect"):
                        asp_i = i
            for row in reader:
                if asp_i is not None and len(row) > asp_i:
                    table[row[0]] = ns.get(row[asp_i], "")
    return lambda gid: table.get(gid) or None


def suggest_disease(keywords: str) -> tuple[str, str] | None:
    kw = (keywords or "").strip()
    if not kw:
        return None
    q = urllib.parse.urlencode({"q": kw, "ontology": "mondo", "rows": "1",
                                "fieldList": "obo_id,label,is_obsolete"})
    try:
        with urllib.request.urlopen(f"https://www.ebi.ac.uk/ols4/api/search?{q}", timeout=30) as r:
            docs = json.load(r)["response"]["docs"]
    except Exception:  # noqa: BLE001
        return None
    if docs and not docs[0].get("is_obsolete"):
        return docs[0]["obo_id"], docs[0]["label"]
    return None


def build_draft(sid: str, meta: dict, genes: list[str], rows, aspect_lookup,
                suggest: bool) -> model.GeneSetInterpretation:
    name = "CGSA_" + sid.replace("-", "_")
    interp = draft_mod.assemble_draft(
        gene_set_id=f"CGSA:{sid.replace('-', '_')}",
        gene_set_name=name,
        collection="CGSA:BENCHMARK",
        rows=rows,
        aspect_lookup=aspect_lookup,
        n_genes=len(genes),
    )
    taxon = meta["taxon"]
    interp.taxon = model.Term(
        id=taxon, label={"NCBITaxon:9606": "Homo sapiens",
                         "NCBITaxon:10090": "Mus musculus"}.get(taxon, ""))
    interp.direction = "NA"
    contexts = []
    if suggest:
        hit = suggest_disease(meta.get("keywords", ""))
        if hit:
            contexts.append(model.BiologicalContext(
                term=model.Term(id=hit[0], label=hit[1]),
                context_type="disease",
                role_note="AUTO-SUGGESTED from benchmark keywords — verify/replace"))
    interp.contexts = contexts
    # default recovery_status for enrichment-recovered terms (curator may override)
    for a in interp.associations:
        a.recovery_status = "annotation_supported"
    return interp


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", help="Comma-separated benchmark set ids (e.g. PMC9750880-1).")
    ap.add_argument("--ids-file", type=Path, help="File of set ids, one per line.")
    ap.add_argument("--out-dir", type=Path, default=ROOT / "curation" / "cgsa" / "drafts")
    ap.add_argument("--cache-dir", type=Path, default=Path("/tmp/cgsa_annot"))
    ap.add_argument("--binary", default=str(ROOT / "target" / "release" / "genesets-rs"))
    ap.add_argument("--no-suggest-context", action="store_true")
    ap.add_argument("--force-annot", action="store_true")
    ap.add_argument("--top", type=int, default=30, help="Keep the top-N terms by p_adjust per set.")
    ap.add_argument("--max-target-size", type=int, default=1200,
                    help="Drop terms annotated to more than this many genes (generic roots).")
    args = ap.parse_args()

    ids: list[str] = []
    if args.ids:
        ids += [x.strip() for x in args.ids.split(",") if x.strip()]
    if args.ids_file:
        ids += [l.strip() for l in args.ids_file.read_text().splitlines() if l.strip()]
    if not ids:
        print("no --ids given", file=sys.stderr)
        return 1

    sets, members, gt_go = load_sets(), load_members(), load_gt_go()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    by_species: dict[str, list[str]] = defaultdict(list)
    for sid in ids:
        taxon = sets[sid]["taxon"].replace("NCBITaxon:", "")
        by_species[taxon].append(sid)

    report = []
    for taxon, sids in by_species.items():
        if taxon not in SPECIES:
            print(f"  ! skipping unsupported taxon {taxon}", file=sys.stderr)
            continue
        label, gaf = SPECIES[taxon]
        species_dir = args.cache_dir / label
        ensure_annotations(species_dir, gaf, args.force_annot)
        queries = {sid: members[sid] for sid in sids}
        text = run_enrichment(species_dir, queries, args.binary)
        aspect_lookup = aspect_lookup_from_terms(species_dir / "terms.tsv")
        rows_by_set: dict[str, list] = defaultdict(list)
        for row in parse_enrichment_tsv(text):
            rows_by_set[row.query_id].append(row)

        for sid in sids:
            # Drop root-like/generic terms (e.g. "protein binding", "cytoplasm")
            # by annotation breadth — a deterministic specificity pre-filter so the
            # curator adjudicates informative terms, not the ontology's roots.
            candidates = [r for r in rows_by_set.get(sid, [])
                          if r.target_size <= args.max_target_size]
            rows = sorted(candidates,
                          key=lambda r: (r.p_adjust if r.p_adjust is not None else 1.0))[:args.top]
            enriched_ids = {r.target_id for r in rows}
            gt_ids = gt_go.get(sid, set())
            recovered = gt_ids & enriched_ids
            mismatch = bool(gt_ids) and len(recovered) / len(gt_ids) < 0.2
            interp = build_draft(sid, sets[sid], members[sid], rows, aspect_lookup,
                                 not args.no_suggest_context)
            out = args.out_dir / f"{interp.gene_set_name}.draft.yaml"
            model.dump_interpretation(interp, out)
            report.append({
                "set_id": sid, "species": label, "n_genes": len(members[sid]),
                "n_sig_terms": len(rows), "gt_grounded": len(gt_ids),
                "gt_recovered": len(recovered),
                "mismatch_flag": "MISMATCH" if mismatch else "",
                "draft": str(out),
            })
            print(f"  {sid} [{label}] -> {out.name}  "
                  f"({len(rows)} terms, GT {len(recovered)}/{len(gt_ids)} recovered"
                  f"{'  ⚠ MISMATCH' if mismatch else ''})")

    rep_path = args.out_dir / "screening_report.tsv"
    with rep_path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["set_id", "species", "n_genes", "n_sig_terms",
                                          "gt_grounded", "gt_recovered", "mismatch_flag", "draft"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(report)
    print(f"\nwrote {len(report)} drafts + {rep_path}")
    print("Next: adjudicate each draft (set category/insight, confirm context, prune), then "
          "`curate validate` and move into curation/cgsa/genesets/.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
