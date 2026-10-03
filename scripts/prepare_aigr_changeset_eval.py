#!/usr/bin/env python3
"""Prepare GOA variants edited by an AI Gene Review (AIGR) change set.

Question: if the AIGR review decisions were applied to GOA human, would
enrichment of external gene sets get better?

This builds on ``prepare_go_eval.py`` (same GO/GOA downloads, terms/closure
tables, common all-GOA background and genesets-rs config shape) and adds
variants in which the AIGR change set is applied *row by row* to the GAF:

- ``aigr_remove``   drop rows reviewed REMOVE
- ``aigr_prune``    drop REMOVE and MARK_AS_OVER_ANNOTATED (over-annotation removal)
- ``aigr_full``     aigr_prune + MODIFY swaps the term for its GO replacement(s)
                    + NEW annotations added
- ``aigr_core``     reviewed rows kept only if ACCEPT (plus MODIFY replacements
                    and NEW); KEEP_AS_NON_CORE/UNDECIDED also dropped. The most
                    aggressive arm.
- ``null_prune_<k>`` negative control: among the same reviewed rows, drop a
                    random sample with the same per-evidence-code counts that
                    aigr_prune drops.

A GAF row is matched to a change-set row on (UniProt accession, GO term,
evidence code, reference), where the reference matches any of the GAF row's
pipe-separated references. GAF rows that no review covers (genes not reviewed,
or annotations added to GOA after the review) are left unchanged. A
(gene, term) pair survives if any of its supporting rows survives, so removing
one row does not remove a term that other rows still support.

The change set is the TSV written by ai-gene-review's
``scripts/export_goa_changeset.py``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_go_eval  # noqa: E402

POLICIES = ("aigr_remove", "aigr_prune", "aigr_full", "aigr_core")
DROP = {
    "aigr_remove": {"REMOVE"},
    "aigr_prune": {"REMOVE", "MARK_AS_OVER_ANNOTATED"},
    "aigr_full": {"REMOVE", "MARK_AS_OVER_ANNOTATED"},
    "aigr_core": {"REMOVE", "MARK_AS_OVER_ANNOTATED", "KEEP_AS_NON_CORE", "UNDECIDED", "PENDING", ""},
}
APPLIES_MODIFY_AND_NEW = {"aigr_full", "aigr_core"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--changeset", required=True, type=Path, help="AIGR change-set TSV")
    parser.add_argument("--ontology-url", default=prepare_go_eval.DEFAULT_ONTOLOGY_URL)
    parser.add_argument("--gaf-url", default=prepare_go_eval.DEFAULT_GAF_URL)
    parser.add_argument(
        "--drop-direct-terms",
        default=",".join(prepare_go_eval.DEFAULT_DROP_DIRECT_TERMS),
        help="GO ids whose direct annotations are dropped from every variant (default: protein binding).",
    )
    parser.add_argument("--null-seeds", type=int, default=20, help="number of null_prune replicates")
    parser.add_argument("--max-p-adjust", type=float, default=0.05)
    parser.add_argument("--force-download", action="store_true")
    return parser.parse_args()


def read_changeset(path: Path) -> tuple[dict, list[dict], dict]:
    """Index reviewed existing annotations by (acc, term, evidence) -> {ref: row}."""
    index: dict[tuple[str, str, str], dict[str, dict]] = defaultdict(dict)
    new_rows: list[dict] = []
    stats: Counter = Counter()
    with path.open() as handle:
        header = [line for line in handle if line.startswith("#")]
    with path.open(newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        for row in reader:
            stats["rows"] += 1
            if row["negated"] == "true":
                stats["negated_skipped"] += 1
                continue
            if row["action"] == "NEW":
                new_rows.append(row)
                continue
            key = (row["uniprot_id"], row["term_id"], row["evidence_type"])
            ref = row["original_reference_id"]
            if ref in index[key]:
                stats["duplicate_keys"] += 1
                if index[key][ref]["action"] != row["action"]:
                    stats["conflicting_duplicate_keys"] += 1
                continue
            index[key][ref] = row
    stats["indexed_existing"] = sum(len(v) for v in index.values())
    stats["new_rows"] = len(new_rows)
    return index, new_rows, {"header": [h.strip() for h in header], **stats}


def load_gaf_rows(
    path: Path, valid_terms: set[str], drop_terms: frozenset[str] = frozenset()
) -> list[tuple[str, str, str, str, tuple[str, ...]]]:
    """NOT-filtered GAF rows as (accession, symbol, term, evidence, refs)."""
    rows = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15 or fields[4] not in valid_terms or not fields[2].strip():
                continue
            if fields[4] in drop_terms:
                continue
            if "NOT" in prepare_go_eval.parse_qualifiers(fields[3]):
                continue
            rows.append((fields[1], fields[2].strip(), fields[4], fields[6], tuple(fields[5].split("|"))))
    return rows


def match_rows(gaf_rows, index) -> tuple[list[dict | None], dict]:
    matched: list[dict | None] = []
    used: set[tuple] = set()
    for acc, _symbol, term, evidence, refs in gaf_rows:
        hit = None
        candidates = index.get((acc, term, evidence))
        if candidates:
            for ref in refs:
                if ref in candidates:
                    hit = candidates[ref]
                    used.add((acc, term, evidence, ref))
                    break
        matched.append(hit)
    stats = {
        "gaf_rows": len(gaf_rows),
        "gaf_rows_matched": sum(1 for m in matched if m is not None),
        "changeset_existing_rows_matched": len(used),
    }
    return matched, stats


def write_pairs(out_dir: Path, name: str, pairs: set[tuple[str, str]]) -> None:
    (out_dir / name).mkdir(parents=True, exist_ok=True)
    with (out_dir / name / "gene_terms.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "term_id"])
        writer.writerows(sorted(pairs))


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    downloads = out_dir / "downloads"
    ontology_file = prepare_go_eval.download(args.ontology_url, downloads / "go-basic.obo", args.force_download)
    gaf_file = prepare_go_eval.download(args.gaf_url, downloads / "goa_human.gaf.gz", args.force_download)

    terms, parents, _ = prepare_go_eval.parse_obo(downloads / "go-basic.obo", {"is_a", "part_of"})
    prepare_go_eval.write_terms_and_closure(out_dir, terms, parents)
    valid_terms = set(terms)

    index, new_rows, cs_stats = read_changeset(args.changeset)
    drop_terms = prepare_go_eval.drop_direct_terms(args.drop_direct_terms)
    gaf_rows = load_gaf_rows(downloads / "goa_human.gaf.gz", valid_terms, drop_terms)
    matched, match_stats = match_rows(gaf_rows, index)

    acc_to_symbol: dict[str, str] = {}
    for acc, symbol, *_ in gaf_rows:
        acc_to_symbol.setdefault(acc, symbol)

    # Baseline: every NOT-filtered row (identical to prepare_go_eval's `all`).
    baseline = {(symbol, term) for _acc, symbol, term, _ev, _refs in gaf_rows}
    background = sorted({symbol for symbol, _ in baseline})
    with (out_dir / "background_all_goa_symbols.txt").open("w") as handle:
        handle.write("gene_id\n")
        handle.writelines(f"{g}\n" for g in background)
    write_pairs(out_dir, "all", baseline)

    reviewed_symbols = sorted({acc_to_symbol[r[0]] for r, m in zip(gaf_rows, matched) if m is not None})
    with (out_dir / "reviewed_genes.txt").open("w") as handle:
        handle.write("gene_id\n")
        handle.writelines(f"{g}\n" for g in reviewed_symbols)

    variant_stats: dict[str, dict] = {"all": {"pairs": len(baseline)}}
    for policy in POLICIES:
        drop = DROP[policy]
        pairs: set[tuple[str, str]] = set()
        counts: Counter = Counter()
        for (acc, symbol, term, _ev, _refs), hit in zip(gaf_rows, matched):
            if hit is None:
                pairs.add((symbol, term))
                continue
            action = hit["action"]
            if policy in APPLIES_MODIFY_AND_NEW and action == "MODIFY":
                replacements = [t for t in hit["replacement_term_ids"].split("|") if t in valid_terms]
                if replacements:
                    counts["modify_rows_replaced"] += 1
                    pairs.update((symbol, t) for t in replacements)
                    continue
                counts["modify_rows_without_valid_replacement_kept"] += 1
                pairs.add((symbol, term))
                continue
            if action in drop:
                counts[f"rows_dropped_{action or 'BLANK'}"] += 1
                continue
            pairs.add((symbol, term))
        if policy in APPLIES_MODIFY_AND_NEW:
            for row in new_rows:
                symbol = acc_to_symbol.get(row["uniprot_id"])
                if symbol and row["term_id"] in valid_terms:
                    counts["new_pairs_added"] += 1
                    pairs.add((symbol, row["term_id"]))
                else:
                    counts["new_rows_unmappable"] += 1
        write_pairs(out_dir, policy, pairs)
        variant_stats[policy] = {
            "pairs": len(pairs),
            "pairs_removed_vs_all": len(baseline - pairs),
            "pairs_added_vs_all": len(pairs - baseline),
            **counts,
        }

    # Null: same per-evidence-code number of dropped rows as aigr_prune, drawn
    # uniformly from all matched (reviewed) rows.
    prune_drop = DROP["aigr_prune"]
    target: Counter = Counter()
    pool: dict[str, list[int]] = defaultdict(list)
    for i, ((_acc, _symbol, _term, evidence, _refs), hit) in enumerate(zip(gaf_rows, matched)):
        if hit is None:
            continue
        pool[evidence].append(i)
        if hit["action"] in prune_drop:
            target[evidence] += 1
    null_names = []
    for seed in range(args.null_seeds):
        rng = random.Random(seed)
        dropped: set[int] = set()
        for evidence, n in target.items():
            dropped.update(rng.sample(pool[evidence], n))
        pairs = {(r[1], r[2]) for i, r in enumerate(gaf_rows) if i not in dropped}
        name = f"null_prune_{seed:02d}"
        write_pairs(out_dir, name, pairs)
        null_names.append(name)
        variant_stats[name] = {"pairs": len(pairs), "pairs_removed_vs_all": len(baseline - pairs)}

    queries_path = out_dir / "queries.gmt"
    variants = ["all", *POLICIES, *null_names]
    prepare_go_eval.write_run_configs(out_dir, variants, queries_path, max_p_adjust=args.max_p_adjust)

    metadata = {
        "generated_at_utc": prepare_go_eval.utc_now(),
        "ontology_file": ontology_file,
        "gaf_file": gaf_file,
        "gaf_header": prepare_go_eval.parse_gaf_header(downloads / "goa_human.gaf.gz"),
        "changeset": {"path": str(args.changeset), **cs_stats},
        "matching": match_stats,
        "dropped_direct_terms": sorted(drop_terms),
        "reviewed_gene_count": len(reviewed_symbols),
        "background_gene_count": len(background),
        "null_target_rows_by_evidence": dict(sorted(target.items())),
        "variants": variant_stats,
    }
    with (out_dir / "metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")

    print(json.dumps({"changeset": cs_stats, "matching": match_stats}, indent=2, default=str))
    for name in ["all", *POLICIES, null_names[0] if null_names else None]:
        if name:
            print(name, variant_stats[name])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
