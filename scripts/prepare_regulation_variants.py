#!/usr/bin/env python3
"""Prepare GO enrichment variants that differ only in how "regulation of X" is handled.

Variants (all share one terms table, one all-GOA background and NOT-filtered
GOA human annotations; only the closure and/or the gene-term table change):

``current``
    What ``prepare_go_eval.py`` does today: is_a + part_of closure. A gene
    annotated to "positive regulation of X" reaches "regulation of X" (is_a)
    but never X.
``drop_regulation``
    is_a + part_of closure, with every regulation term removed both as an
    annotation and as a closure ancestor, so only participant terms are tested.
``regulates_closure``
    Closure over is_a, part_of, regulates, positively_regulates and
    negatively_regulates (the AmiGO/GOlr ``regulates_closure``). A regulator of
    X is counted as a member of X, alongside the participants.
``collapse_to_target``
    Annotations to a regulation term are re-pointed at the process it regulates
    (followed through regulation-of-regulation chains), then regulation terms
    are dropped as in ``drop_regulation``. Same gene membership for X as
    ``regulates_closure`` but without the regulation terms competing for
    significance (and the multiple-testing burden they add).

The full ``go.obo`` is required (go-basic omits the regulates relations).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from prepare_go_eval import file_metadata, write_yaml_value  # noqa: E402

DEFAULT_ONTOLOGY_URL = "https://current.geneontology.org/ontology/go.obo"
DEFAULT_GAF_URL = "https://current.geneontology.org/annotations/goa_human.gaf.gz"
REG_RELS = ("regulates", "positively_regulates", "negatively_regulates")
VARIANTS = ("current", "drop_regulation", "regulates_closure", "collapse_to_target")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--obo", type=Path, help="Local full go.obo (downloaded if omitted).")
    p.add_argument("--gaf", type=Path, help="Local goa_human.gaf.gz (downloaded if omitted).")
    p.add_argument("--queries", type=Path, help="queries.gmt (default: <out-dir>/queries.gmt).")
    p.add_argument("--max-p-adjust", type=float, default=0.05)
    return p.parse_args()


def fetch(url: str, path: Path) -> Path:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(url, headers={"User-Agent": "genesets-rs-eval/0.1"})
        with urllib.request.urlopen(req, timeout=300) as resp:
            path.write_bytes(resp.read())
    return path


def parse_obo(path: Path) -> tuple[dict[str, str], dict[str, str], dict[str, dict[str, set[str]]]]:
    """Return (names, alt_id -> id, {relation: {child: parents}}) for non-obsolete terms.

    Regulates-family edges are kept only between biological_process terms: the
    MF "regulator activity" edges (e.g. enzyme inhibitor activity
    negatively_regulates catalytic activity) are out of scope here.
    """
    names: dict[str, str] = {}
    namespace: dict[str, str] = {}
    alts: dict[str, str] = {}
    edges: dict[str, dict[str, set[str]]] = {r: defaultdict(set) for r in ("is_a", "part_of", *REG_RELS)}
    cur: dict | None = None

    def commit(t: dict | None) -> None:
        if not t or "id" not in t or t.get("obsolete"):
            return
        names[t["id"]] = t.get("name", "")
        namespace[t["id"]] = t.get("namespace", "")
        for a in t.get("alt", ()):
            alts[a] = t["id"]
        for rel, tgt in t.get("rels", ()):
            if rel in edges:
                edges[rel][t["id"]].add(tgt)

    with path.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if line == "[Term]":
                commit(cur)
                cur = {"rels": []}
                continue
            if line.startswith("["):
                commit(cur)
                cur = None
                continue
            if cur is None or not line:
                continue
            key, _, val = line.partition(": ")
            if key == "id":
                cur["id"] = val
            elif key == "name":
                cur["name"] = val
            elif key == "namespace":
                cur["namespace"] = val
            elif key == "is_obsolete" and val == "true":
                cur["obsolete"] = True
            elif key == "alt_id":
                cur.setdefault("alt", []).append(val)
            elif key == "is_a":
                cur["rels"].append(("is_a", val.split()[0]))
            elif key == "relationship":
                parts = val.split()
                if len(parts) >= 2:
                    cur["rels"].append((parts[0], parts[1]))
        commit(cur)
    bp = {t for t, ns in namespace.items() if ns == "biological_process"}
    for rel in edges:
        for child in list(edges[rel]):
            keep = {p for p in edges[rel][child] if p in names}
            if rel in REG_RELS:
                keep = keep & bp if child in bp else set()
            edges[rel][child] = keep
    return names, alts, edges


def regulation_terms(names: dict[str, str], edges: dict) -> set[str]:
    """BP terms with a direct regulates-family edge, plus regulation-named is_a descendants.

    The descendants catch e.g. "regulation of X by Y", which often carries only
    an is_a to "regulation of X". Descendants are restricted to names containing
    "regulation" because GO places signal transduction (and so every signaling
    pathway) is_a "regulation of cellular process"; following is_a blindly would
    classify all signaling as regulation.
    """
    direct = {t for r in REG_RELS for t, ps in edges[r].items() if ps}
    children: dict[str, set[str]] = defaultdict(set)
    for c, ps in edges["is_a"].items():
        for p in ps:
            children[p].add(c)
    out, stack = set(direct), list(direct)
    while stack:
        t = stack.pop()
        for c in children.get(t, ()):
            if c not in out and "regulation" in names.get(c, ""):
                out.add(c)
                stack.append(c)
    return out


def closure(names: dict[str, str], parents: dict[str, set[str]]) -> dict[str, set[str]]:
    memo: dict[str, set[str]] = {}

    def anc(t: str) -> set[str]:
        if t in memo:
            return memo[t]
        seen, stack = {t}, list(parents.get(t, ()))
        while stack:
            p = stack.pop()
            if p not in seen:
                seen.add(p)
                stack.extend(parents.get(p, ()))
        memo[t] = seen
        return seen

    return {t: anc(t) for t in names}


def merge(*maps: dict[str, set[str]]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for m in maps:
        for k, v in m.items():
            out[k] |= v
    return out


def regulated_targets(term: str, edges: dict, reg: set[str], names: dict[str, str]) -> set[str]:
    """Non-regulation process(es) a regulation term ultimately points at.

    Follows regulates edges (and is_a for terms with no direct edge) until a
    non-regulation term is reached.
    """
    out: set[str] = set()
    seen: set[str] = set()
    stack = [term]
    while stack:
        t = stack.pop()
        if t in seen:
            continue
        seen.add(t)
        tg = set().union(*(edges[r].get(t, set()) for r in REG_RELS))
        if not tg:
            tg = {p for p in edges["is_a"].get(t, ()) if p in reg}
        for x in tg:
            if x in reg:
                stack.append(x)
            elif x in names:
                out.add(x)
    return out


def load_gaf(path: Path, names: dict[str, str], alts: dict[str, str]) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15 or "NOT" in f[3].split("|"):
                continue
            term = alts.get(f[4], f[4])
            if term in names and f[2].strip():
                pairs.add((f[2].strip(), term))
    return pairs


def write_closure(path: Path, clos: dict[str, set[str]], exclude: set[str]) -> int:
    n = 0
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["child", "ancestor"])
        for t in sorted(clos):
            for a in sorted(clos[t] - exclude):
                w.writerow([t, a])
                n += 1
    return n


def write_pairs(path: Path, pairs: set[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id", "term_id"])
        w.writerows(sorted(pairs))


def main() -> int:
    args = parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    obo = args.obo or fetch(DEFAULT_ONTOLOGY_URL, out / "downloads" / "go.obo")
    gaf = args.gaf or fetch(DEFAULT_GAF_URL, out / "downloads" / "goa_human.gaf.gz")
    queries = args.queries or out / "queries.gmt"

    names, alts, edges = parse_obo(obo)
    reg = regulation_terms(names, edges)
    pairs = load_gaf(gaf, names, alts)

    with (out / "terms.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["term_id", "name"])
        for t in sorted(names):
            w.writerow([t, names[t]])
    with (out / "regulation_terms.txt").open("w") as fh:
        fh.write("\n".join(sorted(reg)) + "\n")
    with (out / "background_all_goa_symbols.txt").open("w") as fh:
        fh.write("gene_id\n" + "\n".join(sorted({g for g, _ in pairs})) + "\n")

    isa_po = merge(edges["is_a"], edges["part_of"])
    clos_std = closure(names, isa_po)
    clos_reg = closure(names, merge(isa_po, *(edges[r] for r in REG_RELS)))

    collapsed: set[tuple[str, str]] = set()
    unmapped = 0
    for g, t in pairs:
        if t in reg:
            tg = regulated_targets(t, edges, reg, names)
            unmapped += not tg
            collapsed |= {(g, x) for x in tg}
        else:
            collapsed.add((g, t))

    plan = {
        "current": (clos_std, set(), pairs),
        "drop_regulation": (clos_std, reg, {(g, t) for g, t in pairs if t not in reg}),
        "regulates_closure": (clos_reg, set(), pairs),
        "collapse_to_target": (clos_std, reg, collapsed),
    }
    stats: dict = {
        "terms": len(names),
        "regulation_terms": len(reg),
        "annotation_pairs": len(pairs),
        "regulation_annotation_pairs": sum(t in reg for _, t in pairs),
        "collapse_unmapped_pairs": unmapped,
        "variants": {},
    }
    for v, (clos, exclude, vpairs) in plan.items():
        vdir = out / v
        vdir.mkdir(exist_ok=True)
        rows = write_closure(vdir / "closure.tsv", clos, exclude)
        write_pairs(vdir / "gene_terms.tsv", vpairs)
        stats["variants"][v] = {"closure_rows": rows, "gene_term_pairs": len(vpairs)}

        def rel(p: Path) -> str:
            return os.path.relpath(p, start=vdir)

        config = {
            "mode": "matrix",
            "ontology": {
                "terms": rel(out / "terms.tsv"),
                "closure": rel(vdir / "closure.tsv"),
                "annotations": rel(vdir / "gene_terms.tsv"),
            },
            "input": {"queries": rel(queries), "query_format": "gmt"},
            "background": {"file": rel(out / "background_all_goa_symbols.txt")},
            "min_overlap": 2,
            "correction": "bonferroni",
            "max_p_adjust": args.max_p_adjust,
            "output": rel(vdir / "results.tsv"),
        }
        with (vdir / "config.yaml").open("w") as fh:
            write_yaml_value(fh, config)
    stats["inputs"] = {"obo": file_metadata(obo, DEFAULT_ONTOLOGY_URL), "gaf": file_metadata(gaf, DEFAULT_GAF_URL)}
    (out / "metadata.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in stats.items() if k != "inputs"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
