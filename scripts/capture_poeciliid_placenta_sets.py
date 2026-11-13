#!/usr/bin/env python3
"""Capture the Poeciliidae placenta gene sets (Safian, Kraunsoe et al. 2026) as human orthologs.

Source: bioRxiv doi:10.64898/2026.11.11.765921 and its data repository
https://github.com/cardosomoreiralab/fish-placentas-2026 (MIT licence), pinned at
commit 96536fe48c1b4fd9dce2690f97c874007308325c.

Four sets are written, one GMT row each, for ``curation/genesets/lit_members.gmt``:

* ``POECILIID_PLACENTAL_LFC_UP`` / ``_DN`` -- the union, over the four placental
  lineages (P. turneri, P. prolifica, P. bifurca, H. formosa), of genes up- or
  down-regulated in the late follicle cell (LFC) of the placental species
  relative to its non-placental sister (``data/reference/rna/DEGs/*_LFC.csv``,
  rows whose ``lineage`` is the placental species; the per-lineage counts
  reproduce Fig. 5A).
* ``POECILIID_PLACENTATION_LFC_GENES`` / ``_YOLK_SAC_GENES`` -- the genes the
  authors display in Fig. 6D as carrying several kinds of placentation-associated
  change (DEG, positive selection, parallel substitution, marker), transcribed
  from the figure.

Fish gene names are mapped to human symbols with ZFIN's curated human orthology
(``human_orthos.txt``; ``ensembl_1_to_1.txt`` for names that are bare Ensembl
zebrafish ids) and normalized to approved HGNC symbols (pinned monthly release).
A fish name with no ZFIN orthology is accepted only if its upper-cased form is
itself an approved HGNC symbol. Nothing else is inferred; unmapped names are
reported, not guessed.

The DEG tables do not share one gene-id namespace: each pair uses the ids of one
species' annotation (gracilis for gracilis-turneri, formosa for
brachi-formosa, bifurca for picta-bifurca, infans for infans-prolifica). Names
come from the table's ``full_gene_name`` column where present, else from that
species' GFF3 ``gene_name``. The namespaces were established by matching
``full_gene_name`` against each GFF3 (>95% agreement for the chosen species).

Usage::

    python3 scripts/capture_poeciliid_placenta_sets.py \
        --repo /path/to/fish-placentas-2026 --zfin-orthos human_orthos.txt \
        --zfin-ensembl ensembl_1_to_1.txt --hgnc hgnc_complete_set_2026-10-06.txt \
        --out poeciliid.gmt --report poeciliid_mapping.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

# (file stem, id-namespace species, placental lineage label)
DEG_TABLES = [
    ("gracilis_turneri_LFC", "gracilis", "turneri"),
    ("infans_prolifica_LFC", "infans", "prolifica"),
    ("picta_bifurca_LFC", "bifurca", "bifurca"),
    ("brachi_formosa_LFC", "formosa", "formosa"),
]

# Fig. 6D, maternal side, transcribed in reading order. ``g6532`` and ``g1066``
# are printed as bare P. bifurca gene ids; g6532 is named below through its
# orthogroup (ENSDARG00000076534 in all five named annotations); g1066 has no
# name in any annotation and cannot be mapped. ``ldldr`` is printed in the figure
# and is read as ldlr (the paralogs in these annotations are ldlra/ldlrb), and
# ``slc29a`` is read as slc29a1a, the only slc29a1 paralog in these annotations.
FIG6D_LFC = [
    "cldn3d", "diaph2", "gene:ENSDARG00000076534", "ssh1", "irfg3", "sema4ba", "rhoca",
    "hsd17b12a", "pdhb", "papss2b", "abat", "ptgs2b", "degs2", "pde1c", "fut9a",
    "cd63", "cd9a", "vps26bl", "golga7", "rgl1", "shisa4", "emp1", "anxa13", "cav2", "wbp1",
    "plekha1a", "arhgap36", "g1066",
    "esr2a", "pcnx4", "rorcb", "kdm1a", "zfand3", "dhx35", "a2ml",
    "rca2.1", "tnfrsf9", "vrk2", "thbs1",
    "slc29a1a", "slc11a2", "slc1a3a", "piezo2", "slc26a6", "xkrx",
    "paqr7a", "ldlra", "cxadr",
]
FIG6D_YOLK_SAC = [
    "tnfrsf9", "vrk2", "alcama", "rhoca", "irfg3", "fut9a", "nos1", "hsd11b2",
    "esyt2", "cd63", "cd9a", "tbx3a", "macrod2", "g1066",
    "slc29a1a", "slc3a2", "piezo2", "slc24a1", "slc1a3a", "xkrx", "slc26a6",
]


def load_hgnc(path: Path) -> tuple[set[str], dict[str, str]]:
    approved: set[str] = set()
    other: dict[str, set[str]] = defaultdict(set)
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["status"] != "Approved":
                continue
            sym = row["symbol"]
            approved.add(sym)
            for col in ("prev_symbol", "alias_symbol"):
                for s in filter(None, (row.get(col) or "").strip('"').split("|")):
                    other[s].add(sym)
    unambiguous = {s: next(iter(v)) for s, v in other.items() if len(v) == 1 and s not in approved}
    return approved, unambiguous


def load_zfin(orthos: Path, ensembl: Path) -> tuple[dict[str, set[str]], dict[str, str], set[str]]:
    sym2human: dict[str, set[str]] = defaultdict(set)
    with orthos.open() as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[3]:
                sym2human[p[1]].add(p[3])
    ens2sym: dict[str, str] = {}
    zfin_symbols: set[str] = set(sym2human)
    with ensembl.open() as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                ens2sym[p[3]] = p[2]
                zfin_symbols.add(p[2])
    return sym2human, ens2sym, zfin_symbols


class Mapper:
    def __init__(self, hgnc: Path, orthos: Path, ensembl: Path):
        self.approved, self.hgnc_other = load_hgnc(hgnc)
        self.sym2human, self.ens2sym, self.zfin_symbols = load_zfin(orthos, ensembl)

    def norm(self, sym: str) -> str | None:
        if sym in self.approved:
            return sym
        return self.hgnc_other.get(sym)

    def map(self, name: str) -> tuple[list[str], str]:
        """Return (human symbols, method)."""
        if not name or re.fullmatch(r"g\d+", name):
            return [], "unnamed"
        zsym = name
        if name.startswith("gene:ENSDARG"):
            zsym = self.ens2sym.get(name.split(":", 1)[1], "")
            if not zsym:
                return [], "ensdarg_unmapped"
        if zsym in self.sym2human:
            hs = sorted({h for h in (self.norm(s) for s in self.sym2human[zsym]) if h})
            if hs:
                return hs, "zfin_ortholog"
        up = self.norm(zsym) if zsym != zsym.lower() else None  # name already human-style
        if up:
            return [up], "hgnc_symbol"
        if zsym.upper() in self.approved:
            return [zsym.upper()], "uppercase_hgnc"
        return [], "no_ortholog"


def gff_names(path: Path) -> dict[str, str]:
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) > 8 and p[2] == "gene":
                m = re.match(r"ID=([^;\n]+);(?:gene_name=([^;\n]+))?", p[8])
                if m and m.group(2):
                    out[m.group(1)] = m.group(2).split(",")[0]
    return out


def deg_names(repo: Path, direction: str) -> tuple[list[str], dict[str, int]]:
    ref = repo / "data" / "reference"
    names: list[str] = []
    counts = {}
    for stem, ns, lineage in DEG_TABLES:
        gff = gff_names(ref / "annotation" / f"{ns}.gff3.gz")
        n = 0
        with (ref / "rna" / "DEGs" / f"{stem}.csv").open(newline="") as fh:
            for row in csv.DictReader(fh):
                row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
                if row["lineage"] != lineage or row["direction"] != direction:
                    continue
                n += 1
                fg = row.get("full_gene_name", "")
                name = fg.split("-", 1)[1] if "-" in fg else gff.get(row["gene_id"], row["gene_id"])
                names.append(name)
        counts[lineage] = n
    return names, counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--zfin-orthos", type=Path, required=True)
    ap.add_argument("--zfin-ensembl", type=Path, required=True)
    ap.add_argument("--hgnc", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    a = ap.parse_args()
    mapper = Mapper(a.hgnc, a.zfin_orthos, a.zfin_ensembl)

    sets: dict[str, list[str]] = {}
    for d, label in (("up", "UP"), ("down", "DN")):
        names, counts = deg_names(a.repo, d)
        print(f"LFC {label}: per-lineage DEG counts {counts}")
        sets[f"POECILIID_PLACENTAL_LFC_{label}"] = names
    sets["POECILIID_PLACENTATION_LFC_GENES"] = FIG6D_LFC
    sets["POECILIID_PLACENTATION_YOLK_SAC_GENES"] = FIG6D_YOLK_SAC

    with a.out.open("w") as gmt, a.report.open("w") as rep:
        rep.write("set\tfish_name\tmethod\thuman\n")
        for set_name, names in sets.items():
            human: list[str] = []
            methods: dict[str, int] = defaultdict(int)
            for name in names:
                hs, method = mapper.map(name)
                methods[method] += 1
                rep.write(f"{set_name}\t{name}\t{method}\t{','.join(hs)}\n")
                human.extend(h for h in hs if h not in human)
            human.sort()
            gmt.write("\t".join([set_name, "curated_literature_membership", *human]) + "\n")
            print(f"{set_name}: {len(names)} fish entries ({len(set(names))} distinct) -> "
                  f"{len(human)} human symbols; methods {dict(methods)}")


if __name__ == "__main__":
    main()
