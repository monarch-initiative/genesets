"""Build blinded pooled-judgment packets for the regulation-closure recuration."""
import csv, glob, json, sys, collections, yaml
from pathlib import Path
E = Path(sys.argv[1]); OUT = Path(sys.argv[2]); NAGENTS = int(sys.argv[3])
K, MAXSIZE = 3, 1000
VS = ["current", "drop_regulation", "regulates_closure", "collapse_to_target"]
docs = {}
for p in sorted(glob.glob("curation/genesets/*.yaml")):
    d = yaml.safe_load(open(p)); docs[d["gene_set_name"]] = (p, d)
gold = {n: {a["term"]["id"] for a in d.get("associations", [])} for n, (p, d) in docs.items()}
members = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(E / "queries.gmt")}
res = {v: collections.defaultdict(dict) for v in VS}
for v in VS:
    for r in csv.DictReader(open(E / v / "results.tsv"), delimiter="\t"):
        res[v][r["query_id"]][r["target_id"]] = r
names = {r["term_id"]: r["name"] for r in csv.DictReader(open(E / "terms.tsv"), delimiter="\t")}
def term_genes(v):
    clos = collections.defaultdict(set)
    for r in csv.DictReader(open(E / v / "closure.tsv"), delimiter="\t"): clos[r["child"]].add(r["ancestor"])
    tg = collections.defaultdict(set)
    for r in csv.DictReader(open(E / v / "gene_terms.tsv"), delimiter="\t"):
        for a in clos.get(r["term_id"], ()): tg[a].add(r["gene_id"])
    return tg
std, regc = term_genes("current"), term_genes("regulates_closure")
pool = {}
for q in gold:
    if q not in members: continue
    sig = [set(res[v].get(q, {})) for v in VS]
    allv = set.intersection(*sig)
    cands = set()
    for v in VS:
        rows = [r for t, r in res[v].get(q, {}).items()
                if t not in allv and int(r["target_size"]) <= MAXSIZE and t not in gold[q]]
        rows.sort(key=lambda r: float(r["p_value"]))
        cands |= {r["target_id"] for r in rows[:K]}
    if cands: pool[q] = sorted(cands)
# per-candidate blinded evidence
packets = []
for q, cands in sorted(pool.items()):
    p, d = docs[q]; m = set(members[q])
    items = []
    for t in cands:
        direct = sorted(m & std.get(t, set()))
        via_reg = sorted((m & regc.get(t, set())) - set(direct))
        items.append({"term_id": t, "term_label": names[t], "genome_genes_annotated": len(std.get(t, ())),
                      "set_genes_annotated_isa_partof": direct, "set_genes_only_via_regulates": via_reg})
    packets.append({"gene_set_name": q, "file": p, "description": d.get("description", ""),
                    "context": d.get("context"), "existing_associations": [
                        {"id": a["term"]["id"], "label": a["term"].get("label"), "category": a.get("category")}
                        for a in d.get("associations", [])],
                    "membership_size": len(m), "membership": sorted(m), "candidates": items})
OUT.mkdir(parents=True, exist_ok=True)
batches = [[] for _ in range(NAGENTS)]
load = [0] * NAGENTS
for pk in sorted(packets, key=lambda x: -len(x["candidates"])):
    i = load.index(min(load)); batches[i].append(pk); load[i] += len(pk["candidates"])
for i, b in enumerate(batches):
    (OUT / f"batch_{i}.json").write_text(json.dumps(b, indent=1))
json.dump({q: c for q, c in pool.items()}, open(OUT / "pool.json", "w"), indent=1)
print(len(pool), sum(map(len, pool.values())), load)
