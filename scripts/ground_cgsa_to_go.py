#!/usr/bin/env python3
"""Ground cGSA free-text pathway names to GO term IDs via OLS, for GO-ID scoring.

cGSA and its benchmark speak in free-text pathway names; ``genesets-rs`` and this
repo's ``score_method_vs_benchmark.py`` speak in GO term IDs. This bridge closes
the gap: it resolves each free-text name to a GO ID with an OLS4 exact-label
lookup (falling back to top fuzzy hit, flagged), so a cGSA-style output — or the
benchmark's own paper-derived ground truth — can be scored with the same exact-ID
join the rest of the repo uses.

It is deliberately conservative and transparent: only GO-flavored names ground
(the benchmark also carries KEGG / Reactome / MSigDB-signature names that have no
GO ID, and those are reported as unmapped, never forced). Results are cached to a
committed JSON so re-runs are deterministic and offline.

Usage:
  # ground the benchmark's ground-truth names (default GO-database subset only):
  python scripts/ground_cgsa_to_go.py \
      --in evals/cgsa/data/benchmark_ground_truth.tsv \
      --id-col set_id --name-col ground_truth_pathway \
      --map-out evals/cgsa/data/benchmark_ground_truth_go.tsv \
      --results-out evals/cgsa/data/benchmark_ground_truth_go.results.tsv

The ``results-out`` file has ``query_id<TAB>target_id`` columns — the same schema
``score_method_vs_benchmark.load_hits`` reads — so a grounded cGSA prediction
file drops straight into the existing scorer.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = ROOT / "evals" / "cgsa" / "data" / "go_grounding_cache.json"
OLS = "https://www.ebi.ac.uk/ols4/api/search"


def _norm(text: str) -> str:
    """Case/space/punctuation-insensitive key for label equality."""
    return "".join(ch for ch in text.lower() if ch.isalnum())


def _query_go(name: str, exact: str, rows: str) -> list[dict]:
    q = urllib.parse.urlencode({
        "q": name, "ontology": "go", "exact": exact,
        "fieldList": "obo_id,label,is_obsolete", "rows": rows,
    })
    try:
        with urllib.request.urlopen(f"{OLS}?{q}", timeout=30) as resp:
            return json.load(resp)["response"]["docs"]
    except Exception as exc:  # noqa: BLE001 - network best-effort
        print(f"  ! OLS error for {name!r}: {exc}", file=sys.stderr)
        return []


def _ols_lookup(name: str) -> dict | None:
    """Return {'id','label','exact','obsolete'} for the best GO hit, or None.

    ``exact`` is TRUE only when a returned GO label actually equals ``name``
    (normalized): OLS's exact= flag boosts but does not guarantee label
    equality, so we verify it ourselves and fall back to the top fuzzy hit
    (flagged exact=False) otherwise.
    """
    target = _norm(name)
    for doc in _query_go(name, "true", "10"):
        if _norm(doc.get("label", "")) == target:
            return {"id": doc.get("obo_id"), "label": doc.get("label"),
                    "exact": True, "obsolete": bool(doc.get("is_obsolete"))}
    fuzzy = _query_go(name, "false", "1")
    if fuzzy:
        d = fuzzy[0]
        return {"id": d.get("obo_id"), "label": d.get("label"),
                "exact": False, "obsolete": bool(d.get("is_obsolete"))}
    return None


def _load_cache(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--in", dest="infile", required=True, type=Path)
    parser.add_argument("--id-col", required=True)
    parser.add_argument("--name-col", required=True)
    parser.add_argument("--map-out", required=True, type=Path)
    parser.add_argument("--results-out", type=Path, default=None,
                        help="Optional query_id<TAB>target_id file for the repo scorer.")
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--fuzzy-ok", action="store_true",
                        help="Accept top non-exact hit as a match (flagged). Default: exact only.")
    parser.add_argument("--sleep", type=float, default=0.1, help="Delay between OLS calls (s).")
    args = parser.parse_args()

    with args.infile.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    names = sorted({r[args.name_col].strip() for r in rows if r[args.name_col].strip()})

    cache = _load_cache(args.cache)
    resolved = 0
    for name in names:
        if name in cache:
            continue
        hit = _ols_lookup(name)
        cache[name] = hit
        resolved += 1
        if resolved % 25 == 0:  # checkpoint so long runs survive interruption
            args.cache.parent.mkdir(parents=True, exist_ok=True)
            args.cache.write_text(json.dumps(cache, indent=2, sort_keys=True))
        time.sleep(args.sleep)
    if resolved:
        args.cache.parent.mkdir(parents=True, exist_ok=True)
        args.cache.write_text(json.dumps(cache, indent=2, sort_keys=True))
        print(f"resolved {resolved} new names (cache: {args.cache})", file=sys.stderr)

    def accept(hit) -> bool:
        return bool(hit and hit["id"] and not hit["obsolete"] and (hit["exact"] or args.fuzzy_ok))

    map_rows, results_rows = [], []
    mapped = set()
    for r in rows:
        name = r[args.name_col].strip()
        sid = r[args.id_col].strip()
        hit = cache.get(name)
        ok = accept(hit)
        map_rows.append([sid, name,
                         hit["id"] if ok else "",
                         hit["label"] if ok else "",
                         "exact" if (ok and hit["exact"]) else ("fuzzy" if ok else "unmapped")])
        if ok:
            results_rows.append([sid, hit["id"]])
            mapped.add(name)

    with args.map_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["set_id", "pathway_name", "go_id", "go_label", "match"])
        writer.writerows(map_rows)

    if args.results_out:
        with args.results_out.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["query_id", "target_id"])
            writer.writerows(results_rows)

    n_names = len(names)
    print(f"names: {n_names}   grounded to GO: {len(mapped)} "
          f"({len(mapped) / n_names:.0%})   rows emitted: {len(results_rows)}")
    if args.results_out:
        print(f"wrote {args.map_out} and {args.results_out}")
    else:
        print(f"wrote {args.map_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
