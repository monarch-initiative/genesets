#!/usr/bin/env python3
"""Prepare official GO ontology and human GAF tables for genesets-rs evals.

This helper intentionally lives outside the Rust core. It downloads official GO
files, parses the OBO graph for is_a/part_of closure, parses GAF annotations,
and writes lean TSV files consumed by the CLI.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path


DEFAULT_ONTOLOGY_URL = "https://current.geneontology.org/ontology/go-basic.obo"
DEFAULT_GAF_URL = "https://current.geneontology.org/annotations/goa_human.gaf.gz"
VARIANTS = {
    "all": {
        "description": "NOT-filtered GOA human annotations, all evidence codes, contributes_to retained.",
        "exclude_qualifiers": {"NOT"},
        "include_evidence": None,
    },
    "no_contributes_to": {
        "description": "NOT-filtered GOA human annotations, all evidence codes, contributes_to removed.",
        "exclude_qualifiers": {"NOT", "contributes_to"},
        "include_evidence": None,
    },
    "iba": {
        "description": "NOT-filtered GOA human annotations, IBA evidence only.",
        "exclude_qualifiers": {"NOT"},
        "include_evidence": {"IBA"},
    },
    "iba_iea": {
        "description": "NOT-filtered GOA human annotations, IBA or IEA evidence only.",
        "exclude_qualifiers": {"NOT"},
        "include_evidence": {"IBA", "IEA"},
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--ontology-url", default=DEFAULT_ONTOLOGY_URL)
    parser.add_argument("--gaf-url", default=DEFAULT_GAF_URL)
    parser.add_argument(
        "--gene-id-field",
        choices=["db_object_id", "db_object_symbol"],
        default="db_object_symbol",
        help="GAF column to use as the gene ID. MyGeneset query evals use symbols.",
    )
    parser.add_argument(
        "--relations",
        default="is_a,part_of",
        help="Comma-separated ontology relations used for closure.",
    )
    parser.add_argument(
        "--variants",
        default="all,no_contributes_to,iba,iba_iea",
        help=f"Comma-separated variants from: {', '.join(VARIANTS)}",
    )
    parser.add_argument(
        "--max-p-adjust",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff to write into generated genesets-rs configs.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Download inputs even if cached files already exist.",
    )
    return parser.parse_args()


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def download(url: str, path: Path, force: bool) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    if force or not path.exists():
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "genesets-rs-eval/0.1"},
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            path.write_bytes(response.read())
    return file_metadata(path, url)


def file_metadata(path: Path, url: str) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {
        "url": url,
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def parse_obo(path: Path, relation_names: set[str]) -> tuple[dict[str, dict], dict[str, set[str]], dict]:
    terms: dict[str, dict] = {}
    parents: dict[str, set[str]] = defaultdict(set)
    current: dict[str, object] | None = None
    stats = {
        "stanzas": 0,
        "non_obsolete_terms": 0,
        "obsolete_terms": 0,
        "is_a_edges": 0,
        "part_of_edges": 0,
        "relation_edges_used": 0,
    }

    def commit(term: dict[str, object] | None) -> None:
        if not term or "id" not in term:
            return
        stats["stanzas"] += 1
        term_id = str(term["id"])
        if term.get("is_obsolete") == "true":
            stats["obsolete_terms"] += 1
            return
        stats["non_obsolete_terms"] += 1
        terms[term_id] = {
            "id": term_id,
            "name": str(term.get("name", "")),
            "namespace": str(term.get("namespace", "")),
        }
        for parent in term.get("parents", set()):
            parents[term_id].add(parent)

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line == "[Term]":
                commit(current)
                current = {"parents": set()}
                continue
            if line.startswith("["):
                commit(current)
                current = None
                continue
            if current is None:
                continue

            key, _, value = line.partition(":")
            value = value.strip()
            if key == "id":
                current["id"] = value
            elif key == "name":
                current["name"] = value
            elif key == "namespace":
                current["namespace"] = value
            elif key == "is_obsolete":
                current["is_obsolete"] = value
            elif key == "is_a" and "is_a" in relation_names:
                parent = value.split()[0]
                current["parents"].add(parent)
                stats["is_a_edges"] += 1
                stats["relation_edges_used"] += 1
            elif key == "relationship":
                parts = value.split()
                if len(parts) >= 2 and parts[0] in relation_names:
                    current["parents"].add(parts[1])
                    if parts[0] == "part_of":
                        stats["part_of_edges"] += 1
                    stats["relation_edges_used"] += 1
        commit(current)

    # Drop parent references to obsolete or external terms not present in go-basic.
    for child in list(parents):
        parents[child] = {parent for parent in parents[child] if parent in terms}
    return terms, parents, stats


def ancestors_for(term_id: str, parents: dict[str, set[str]], memo: dict[str, set[str]]) -> set[str]:
    if term_id in memo:
        return memo[term_id]
    seen = {term_id}
    stack = list(parents.get(term_id, ()))
    while stack:
        parent = stack.pop()
        if parent in seen:
            continue
        seen.add(parent)
        stack.extend(parents.get(parent, ()))
    memo[term_id] = seen
    return seen


def write_terms_and_closure(
    out_dir: Path,
    terms: dict[str, dict],
    parents: dict[str, set[str]],
) -> dict:
    terms_path = out_dir / "terms.tsv"
    closure_path = out_dir / "closure.tsv"
    term_ids = sorted(terms)

    with terms_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["term_id", "name"])
        for term_id in term_ids:
            writer.writerow([term_id, terms[term_id]["name"]])

    closure_rows = 0
    memo: dict[str, set[str]] = {}
    with closure_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["child", "ancestor"])
        for term_id in term_ids:
            for ancestor in sorted(ancestors_for(term_id, parents, memo)):
                writer.writerow([term_id, ancestor])
                closure_rows += 1

    return {
        "terms_path": str(terms_path),
        "closure_path": str(closure_path),
        "term_count": len(term_ids),
        "closure_rows": closure_rows,
    }


def parse_qualifiers(value: str) -> set[str]:
    return {part for part in value.split("|") if part}


def parse_gaf_header(path: Path) -> dict:
    header: dict[str, object] = {
        "raw": [],
        "all_generated_by": [],
        "all_date_generated": [],
        "all_go_version": [],
    }
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("!"):
                break
            clean = line.strip()
            header["raw"].append(clean)
            if clean.startswith("!gaf-version:"):
                header["gaf_version"] = clean.partition(":")[2].strip()
            elif clean.startswith("!date-generated:"):
                value = clean.partition(":")[2].strip()
                header["all_date_generated"].append(value)
                header.setdefault("date_generated", value)
            elif clean.startswith("!generated-by:"):
                value = clean.partition(":")[2].strip()
                header["all_generated_by"].append(value)
                header.setdefault("generated_by", value)
            elif clean.startswith("!go-version:"):
                value = clean.partition(":")[2].strip()
                header["all_go_version"].append(value)
                header.setdefault("go_version", value)
    return header


def variant_allows(fields: list[str], variant: dict) -> bool:
    qualifiers = parse_qualifiers(fields[3])
    if qualifiers.intersection(variant["exclude_qualifiers"]):
        return False
    evidence = fields[6]
    include_evidence = variant["include_evidence"]
    return include_evidence is None or evidence in include_evidence


def write_annotations(
    path: Path,
    out_dir: Path,
    selected_variants: list[str],
    gene_id_field: str,
    valid_terms: set[str],
) -> tuple[dict, dict[str, set[str]]]:
    for name in selected_variants:
        (out_dir / name).mkdir(parents=True, exist_ok=True)
    variant_files = {
        name: (out_dir / name / "gene_terms.tsv").open("w", newline="")
        for name in selected_variants
    }
    writers = {
        name: csv.writer(handle, delimiter="\t", lineterminator="\n")
        for name, handle in variant_files.items()
    }
    for writer in writers.values():
        writer.writerow(["gene_id", "term_id"])

    variant_pairs: dict[str, set[tuple[str, str]]] = {name: set() for name in selected_variants}
    variant_genes: dict[str, set[str]] = {name: set() for name in selected_variants}
    evidence_counts: dict[str, int] = defaultdict(int)
    qualifier_counts: dict[str, int] = defaultdict(int)
    stats = {
        "raw_annotation_lines": 0,
        "malformed_lines": 0,
        "unknown_go_terms": 0,
        "not_qualified_lines": 0,
        "contributes_to_lines": 0,
        "variant_rows": {},
        "variant_unique_pairs": {},
        "variant_gene_counts": {},
        "evidence_counts": evidence_counts,
        "qualifier_counts": qualifier_counts,
    }

    gene_index = 1 if gene_id_field == "db_object_id" else 2
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("!"):
                continue
            stats["raw_annotation_lines"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15:
                stats["malformed_lines"] += 1
                continue
            term_id = fields[4]
            if term_id not in valid_terms:
                stats["unknown_go_terms"] += 1
                continue
            gene_id = fields[gene_index].strip()
            if not gene_id:
                continue
            evidence_counts[fields[6]] += 1
            qualifiers = parse_qualifiers(fields[3])
            if not qualifiers:
                qualifier_counts["<none>"] += 1
            for qualifier in qualifiers:
                qualifier_counts[qualifier] += 1
            if "NOT" in qualifiers:
                stats["not_qualified_lines"] += 1
            if "contributes_to" in qualifiers:
                stats["contributes_to_lines"] += 1

            for variant_name in selected_variants:
                if variant_allows(fields, VARIANTS[variant_name]):
                    variant_pairs[variant_name].add((gene_id, term_id))
                    variant_genes[variant_name].add(gene_id)

    for variant_name, pairs in variant_pairs.items():
        for gene_id, term_id in sorted(pairs):
            writers[variant_name].writerow([gene_id, term_id])
        stats["variant_unique_pairs"][variant_name] = len(pairs)
        stats["variant_gene_counts"][variant_name] = len(variant_genes[variant_name])
        stats["variant_rows"][variant_name] = len(pairs)

    for handle in variant_files.values():
        handle.close()

    # Use the broad all-annotations universe as the common background for every
    # variant so evidence-code comparisons are not confounded by a changing
    # population size.
    background_genes = variant_genes[selected_variants[0] if "all" not in variant_genes else "all"]
    background_path = out_dir / "background_all_goa_symbols.txt"
    with background_path.open("w", newline="") as handle:
        handle.write("gene_id\n")
        for gene_id in sorted(background_genes):
            handle.write(f"{gene_id}\n")

    stats["background_path"] = str(background_path)
    stats["background_gene_count"] = len(background_genes)
    stats["evidence_counts"] = dict(sorted(evidence_counts.items()))
    stats["qualifier_counts"] = dict(sorted(qualifier_counts.items()))
    return stats, variant_genes


def yaml_scalar(value: object) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value)
    if text == "" or any(ch in text for ch in ":#[]{}&,*!|>'\"%@`") or text.strip() != text:
        return json.dumps(text)
    return text


def write_yaml_value(handle, value: object, indent: int = 0) -> None:
    pad = " " * indent
    if isinstance(value, dict):
        for key, item in value.items():
            if isinstance(item, (dict, list)):
                handle.write(f"{pad}{key}:\n")
                write_yaml_value(handle, item, indent + 2)
            else:
                handle.write(f"{pad}{key}: {yaml_scalar(item)}\n")
    elif isinstance(value, list):
        for item in value:
            if isinstance(item, (dict, list)):
                handle.write(f"{pad}-\n")
                write_yaml_value(handle, item, indent + 2)
            else:
                handle.write(f"{pad}- {yaml_scalar(item)}\n")
    else:
        handle.write(f"{pad}{yaml_scalar(value)}\n")


def write_run_configs(
    out_dir: Path,
    selected_variants: list[str],
    queries_path: Path,
    min_overlap: int = 2,
    max_p_adjust: float | None = 0.05,
) -> None:
    for variant_name in selected_variants:
        config = {
            "mode": "matrix",
            "ontology": {
                "terms": str(out_dir / "terms.tsv"),
                "closure": str(out_dir / "closure.tsv"),
                "annotations": str(out_dir / variant_name / "gene_terms.tsv"),
            },
            "input": {
                "queries": str(queries_path),
                "query_format": "gmt",
            },
            "background": {"file": str(out_dir / "background_all_goa_symbols.txt")},
            "min_overlap": min_overlap,
            "correction": "bonferroni",
            "output": str(out_dir / variant_name / "results.tsv"),
        }
        if max_p_adjust is not None:
            config["max_p_adjust"] = max_p_adjust
        with (out_dir / variant_name / "config.yaml").open("w") as handle:
            write_yaml_value(handle, config)


def main() -> int:
    args = parse_args()
    selected_variants = [name.strip() for name in args.variants.split(",") if name.strip()]
    if not 0 <= args.max_p_adjust <= 1:
        raise SystemExit("--max-p-adjust must be between 0 and 1")
    unknown = [name for name in selected_variants if name not in VARIANTS]
    if unknown:
        raise SystemExit(f"unknown variants: {', '.join(unknown)}")
    if "all" not in selected_variants:
        raise SystemExit("the all variant is required to define the common background")

    started = time.perf_counter()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    downloads = args.out_dir / "downloads"

    ontology_path = downloads / "go-basic.obo"
    gaf_path = downloads / "goa_human.gaf.gz"
    ontology_file = download(args.ontology_url, ontology_path, args.force_download)
    gaf_file = download(args.gaf_url, gaf_path, args.force_download)

    relation_names = {name.strip() for name in args.relations.split(",") if name.strip()}
    terms, parents, ontology_stats = parse_obo(ontology_path, relation_names)
    table_stats = write_terms_and_closure(args.out_dir, terms, parents)
    gaf_header = parse_gaf_header(gaf_path)
    annotation_stats, _ = write_annotations(
        gaf_path,
        args.out_dir,
        selected_variants,
        args.gene_id_field,
        set(terms),
    )

    queries_path = args.out_dir / "queries.gmt"
    if not queries_path.exists():
        print(
            f"warning: {queries_path} does not exist; run fetch_mygeneset_eval.py first "
            "or copy queries.gmt before running genesets-rs configs",
            file=sys.stderr,
        )
    write_run_configs(
        args.out_dir,
        selected_variants,
        queries_path,
        max_p_adjust=args.max_p_adjust,
    )

    metadata = {
        "generated_at_utc": utc_now(),
        "runtime_seconds": round(time.perf_counter() - started, 3),
        "prep_strategy": "native_python_obo_gaf_parser",
        "ontology": {
            "source": "Gene Ontology official current release",
            "file": ontology_file,
            "relations_for_closure": sorted(relation_names),
            "stats": ontology_stats,
            "tables": table_stats,
        },
        "annotations": {
            "source": "GOA human official GAF from Gene Ontology current annotations",
            "file": gaf_file,
            "gaf_header": gaf_header,
            "gene_id_field": args.gene_id_field,
            "filters": {
                "always_exclude_qualifier": "NOT",
                "variant_definitions": {
                    name: {
                        "description": VARIANTS[name]["description"],
                        "exclude_qualifiers": sorted(VARIANTS[name]["exclude_qualifiers"]),
                        "include_evidence": sorted(VARIANTS[name]["include_evidence"])
                        if VARIANTS[name]["include_evidence"]
                        else "all",
                    }
                    for name in selected_variants
                },
            },
            "stats": annotation_stats,
        },
        "background": {
            "policy": "Common background from NOT-filtered all-evidence GOA human genes, using the configured GAF gene_id_field.",
            "path": annotation_stats["background_path"],
            "gene_count": annotation_stats["background_gene_count"],
        },
        "outputs": {
            "result_filter": {
                "max_p_adjust": args.max_p_adjust,
                "correction": "bonferroni",
            },
            "terms": str(args.out_dir / "terms.tsv"),
            "closure": str(args.out_dir / "closure.tsv"),
            "variants": {
                name: {
                    "annotations": str(args.out_dir / name / "gene_terms.tsv"),
                    "config": str(args.out_dir / name / "config.yaml"),
                    "intended_results": str(args.out_dir / name / "results.tsv"),
                }
                for name in selected_variants
            },
        },
    }

    with (args.out_dir / "metadata.yaml").open("w") as handle:
        write_yaml_value(handle, metadata)
    with (args.out_dir / "metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Prepared GO eval tables under {args.out_dir}")
    for name in selected_variants:
        print(
            f"{name}: {annotation_stats['variant_unique_pairs'][name]} pairs, "
            f"{annotation_stats['variant_gene_counts'][name]} genes"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
