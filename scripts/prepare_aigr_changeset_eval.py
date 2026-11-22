#!/usr/bin/env python3
"""Prepare GOA variants edited by an AI Gene Review (AIGR) change set.

Question: if the AIGR review decisions were applied to GOA human, would
enrichment of external gene sets get better?

This builds on ``prepare_go_eval.py`` (same GO/GOA downloads, terms/closure
tables, common all-GOA background and genesets-rs config shape) and adds
variants in which the AIGR change set is applied *row by row* to the GAF:

Every AIGR variant applies the uncontroversial part of the change set: rows
reviewed REMOVE are dropped, MODIFY swaps the term for its GO replacement(s)
(a MODIFY whose replacements are all proposed terms keeps the original term),
and NEW annotations are added. The open question is what to do with the two
"true but" actions, so they are crossed as a 2x2 on top of that base:

- ``aigr_base``          REMOVE dropped, MODIFY and NEW applied
- ``aigr_drop_moa``      aigr_base + MARK_AS_OVER_ANNOTATED rows dropped
- ``aigr_drop_noncore``  aigr_base + KEEP_AS_NON_CORE rows dropped
- ``aigr_drop_both``     aigr_base + both dropped
- ``null_moa_<k>`` / ``null_noncore_<k>``  negative controls: starting from
                         aigr_base, drop a random sample of the reviewed rows
                         that aigr_base keeps unchanged, with the same
                         per-evidence-code counts as the factor drops.

ACCEPT, TOO_GENERAL, UNDECIDED and PENDING rows are kept in every variant.

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

BASE_DROP = {"REMOVE"}
FACTORS = {"moa": {"MARK_AS_OVER_ANNOTATED"}, "noncore": {"KEEP_AS_NON_CORE"}}
POLICIES = ("aigr_base", "aigr_drop_moa", "aigr_drop_noncore", "aigr_drop_both")
DROP = {
    "aigr_base": BASE_DROP,
    "aigr_drop_moa": BASE_DROP | FACTORS["moa"],
    "aigr_drop_noncore": BASE_DROP | FACTORS["noncore"],
    "aigr_drop_both": BASE_DROP | FACTORS["moa"] | FACTORS["noncore"],
}

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
    parser.add_argument(
        "--gaf-db",
        default=",".join(prepare_go_eval.DEFAULT_GAF_DBS),
        help="GAF DB column values kept (default: UniProtKB proteins only).",
    )
    parser.add_argument("--null-seeds", type=int, default=10, help="null replicates per factor")
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
            if row.get("retired") == "true":  # no longer in GOA; nothing to match
                stats["retired_skipped"] += 1
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
    path: Path,
    valid_terms: set[str],
    drop_terms: frozenset[str] = frozenset(),
    gaf_dbs: frozenset[str] = frozenset(prepare_go_eval.DEFAULT_GAF_DBS),
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
            if gaf_dbs and fields[0] not in gaf_dbs:
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
    gaf_dbs = prepare_go_eval.drop_direct_terms(args.gaf_db)
    gaf_rows = load_gaf_rows(downloads / "goa_human.gaf.gz", valid_terms, drop_terms, gaf_dbs)
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
    def build_pairs(drop: set[str], skip_rows: set[int] = frozenset()) -> tuple[set, Counter]:
        """Apply REMOVE/MODIFY/NEW plus ``drop``; rows in ``skip_rows`` are also dropped."""
        pairs: set[tuple[str, str]] = set()
        counts: Counter = Counter()
        for i, ((_acc, symbol, term, _ev, _refs), hit) in enumerate(zip(gaf_rows, matched)):
            if i in skip_rows:
                counts["rows_dropped_null"] += 1
                continue
            if hit is None:
                pairs.add((symbol, term))
                continue
            action = hit["action"]
            if action == "MODIFY":
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
        for row in new_rows:
            symbol = acc_to_symbol.get(row["uniprot_id"])
            if symbol and row["term_id"] in valid_terms:
                counts["new_pairs_added"] += 1
                pairs.add((symbol, row["term_id"]))
            else:
                counts["new_rows_unmappable"] += 1
        return pairs, counts

    for policy in POLICIES:
        pairs, counts = build_pairs(DROP[policy])
        write_pairs(out_dir, policy, pairs)
        if policy == "aigr_base":
            base_pairs = set(pairs)
        variant_stats[policy] = {
            "pairs": len(pairs),
            "pairs_removed_vs_all": len(baseline - pairs),
            "pairs_added_vs_all": len(pairs - baseline),
            **counts,
        }

    # Nulls, one family per factor: start from aigr_base and drop a random
    # sample of the reviewed rows that aigr_base keeps unchanged (not REMOVE,
    # not MODIFY), matching the factor's per-evidence-code row counts. The
    # comparison is then "these rows" against "the same number of other rows".
    base_kept = [
        i for i, hit in enumerate(matched)
        if hit is not None and hit["action"] not in BASE_DROP and hit["action"] != "MODIFY"
    ]
    pool: dict[str, list[int]] = defaultdict(list)
    for i in base_kept:
        pool[gaf_rows[i][3]].append(i)
    null_names = []
    null_targets = {}
    for factor, actions in FACTORS.items():
        target: Counter = Counter(gaf_rows[i][3] for i in base_kept if matched[i]["action"] in actions)
        null_targets[factor] = dict(sorted(target.items()))
        for seed in range(args.null_seeds):
            rng = random.Random(f"{factor}:{seed}")
            dropped: set[int] = set()
            for evidence, n in target.items():
                dropped.update(rng.sample(pool[evidence], n))
            pairs, _ = build_pairs(BASE_DROP, dropped)
            name = f"null_{factor}_{seed:02d}"
            write_pairs(out_dir, name, pairs)
            null_names.append(name)
            variant_stats[name] = {
                "pairs": len(pairs),
                "pairs_removed_vs_all": len(baseline - pairs),
                "pairs_removed_vs_base": len(base_pairs - pairs),
            }

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
        "gaf_db_kept": sorted(gaf_dbs) or "all",
        "reviewed_gene_count": len(reviewed_symbols),
        "background_gene_count": len(background),
        "null_target_rows_by_evidence": null_targets,
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
