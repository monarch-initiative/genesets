"""End-to-end check of prepare-chem2go and normalize-chemicals on a tiny world.

palmitic acid (CHEBI:20) and palmitate (CHEBI:21) are conjugates and collapse
to the pH 7.3 form palmitate. GO defines fatty acid metabolism over the class
fatty acid (CHEBI:10, collapsed into fatty acid anion CHEBI:11), so palmitate
reaches it by roll-up. One Rhea reaction (palmitate + water = glucose, a fake
equation) carries the ligase MF and the enzyme bridge.
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

from genesets_workflows.sources import chem2go, normalize_chemicals

OBO = "http://purl.obolibrary.org/obo/"

CHEBI_OBO = """format-version: 1.2
data-version: chebi/test

[Term]
id: CHEBI:1
name: chemical entity

[Term]
id: CHEBI:10
name: fatty acid
is_a: CHEBI:1
relationship: RO:0018034 CHEBI:11

[Term]
id: CHEBI:11
name: fatty acid anion
is_a: CHEBI:1

[Term]
id: CHEBI:20
name: palmitic acid
is_a: CHEBI:10
relationship: RO:0018034 CHEBI:21

[Term]
id: CHEBI:21
name: palmitate
alt_id: CHEBI:29
is_a: CHEBI:11

[Term]
id: CHEBI:15377
name: water
is_a: CHEBI:1

[Term]
id: CHEBI:40
name: glucose
is_a: CHEBI:1
"""


def go_node(curie: str, label: str) -> dict:
    return {"id": OBO + curie.replace(":", "_"), "lbl": label, "type": "CLASS"}


def is_a(child: str, parent: str) -> dict:
    return {"sub": OBO + child.replace(":", "_"), "pred": "is_a", "obj": OBO + parent.replace(":", "_")}


def definition(term: str, genus: str, prop: str, filler: str) -> dict:
    uri = lambda curie: OBO + curie.replace(":", "_")  # noqa: E731
    return {
        "definedClassId": uri(term),
        "genusIds": [uri(genus)],
        "restrictions": [{"propertyId": uri(prop), "fillerId": uri(filler)}],
    }


GO_PLUS = {
    "graphs": [
        {
            "id": "go-plus",
            "meta": {"version": "go/test"},
            "nodes": [
                go_node("GO:0008152", "metabolic process"),
                go_node("GO:0006631", "fatty acid metabolic process"),
                go_node("GO:0003824", "catalytic activity"),
                go_node("GO:0004467", "long-chain fatty acid-CoA ligase activity"),
                go_node("GO:0050896", "response to stimulus"),
                go_node("GO:0070542", "response to fatty acid"),
            ],
            "edges": [
                is_a("GO:0006631", "GO:0008152"),
                is_a("GO:0004467", "GO:0003824"),
                is_a("GO:0070542", "GO:0050896"),
            ],
            "logicalDefinitionAxioms": [
                definition("GO:0006631", "GO:0008152", "RO:0004007", "CHEBI:10"),
                definition("GO:0070542", "GO:0050896", "RO:0004007", "CHEBI:10"),
            ],
        }
    ]
}


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def build_sources(tmp_path: Path) -> list[str]:
    src = tmp_path / "src"
    src.mkdir()
    write(src / "go-plus.json", json.dumps(GO_PLUS))
    with gzip.open(src / "chebi_lite.obo.gz", "wt") as handle:
        handle.write(CHEBI_OBO)
    write(
        src / "database_accession.tsv",
        "id\tcompound_id\taccession_number\ttype\tstatus_id\tsource_id\n"
        "1\t20\tHMDB0000220\tMANUAL_X_REF\t1\t35\n"
        "2\t40\tC00031\tMANUAL_X_REF\t1\t45\n",
    )
    write(src / "source.tsv", "id\tname\turl\tprefix\tdescription\n35\tHMDB\t\thmdb\t\n45\tKEGG\t\tkegg.compound\t\n")
    write(
        src / "rhea-reactions.txt",
        "ENTRY       RHEA:100\nEQUATION    CHEBI:21 + CHEBI:15377 = CHEBI:40\n///\n",
    )
    write(src / "rhea2go.tsv", "RHEA_ID\tDIRECTION\tMASTER_ID\tID\n100\tUN\t100\tGO:0004467\n")
    write(src / "chebi_pH7_3_mapping.tsv", "CHEBI\tCHEBI_PH7_3\tORIGIN\n20\t21\tcuration\n21\t21\tcuration\n10\t11\tcuration\n")
    write(src / "rhea2uniprot.tsv", "RHEA_ID\tDIRECTION\tMASTER_ID\tID\n100\tUN\t100\tP1\n")
    write(src / "goa.gaf", "!gaf-version: 2.2\nUniProtKB\tP1\tLIG\tinvolved_in\tGO:0006631\tPMID:1\tIDA\t\tP\n")
    return [
        "--go-plus", str(src / "go-plus.json"),
        "--chebi", str(src / "chebi_lite.obo.gz"),
        "--chebi-xrefs", str(src / "database_accession.tsv"),
        "--chebi-xref-sources", str(src / "source.tsv"),
        "--rhea-reactions", str(src / "rhea-reactions.txt"),
        "--rhea2go", str(src / "rhea2go.tsv"),
        "--rhea-ph73", str(src / "chebi_pH7_3_mapping.tsv"),
        "--rhea2uniprot", str(src / "rhea2uniprot.tsv"),
        "--gaf", str(src / "goa.gaf"),
    ]


def read_pairs(path: Path) -> set[tuple[str, str]]:
    with path.open() as handle:
        return {(row["gene_id"], row["term_id"]) for row in csv.DictReader(handle, delimiter="\t")}


def test_prepare_chem2go_end_to_end(tmp_path):
    out = tmp_path / "out"
    args = build_sources(tmp_path) + ["--out-dir", str(out), "--max-rollup-class-size", "3"]
    assert chem2go.main(args) == 0

    with (out / "normalization.tsv").open() as handle:
        canonical = {row["chebi_id"]: row["canonical_id"] for row in csv.DictReader(handle, delimiter="\t")}
    assert canonical["CHEBI:20"] == "CHEBI:21"
    assert canonical["CHEBI:29"] == "CHEBI:21"
    assert canonical["CHEBI:10"] == "CHEBI:11"

    # Roll-up through fatty acid; "response to fatty acid" is filtered by link kind.
    assert read_pairs(out / "annotations.logical_definition.tsv") == {
        ("CHEBI:11", "GO:0006631"),
        ("CHEBI:21", "GO:0006631"),
    }
    # Water is currency, so it is not a reaction participant.
    assert read_pairs(out / "annotations.rhea.tsv") == {
        ("CHEBI:21", "GO:0004467"),
        ("CHEBI:40", "GO:0004467"),
    }
    assert read_pairs(out / "annotations.enzyme.tsv") == {
        ("CHEBI:21", "GO:0006631"),
        ("CHEBI:40", "GO:0006631"),
    }
    assert (out / "annotations.species_independent.tsv").exists()
    assert len(read_pairs(out / "annotations.all.tsv")) == 5

    with (out / "provenance.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rollup = [r for r in rows if r["chemical_id"] == "CHEBI:21" and r["source"] == "logical_definition"]
    assert rollup[0]["via_chemical_id"] == "CHEBI:11"
    assert rollup[0]["rollup"] == "true"

    metadata = json.loads((out / "metadata.json").read_text())
    assert metadata["versions"] == {"go_plus": "go/test", "chebi": "chebi/test"}

    closure = (out / "closure.tsv").read_text().splitlines()
    assert "GO:0006631\tGO:0008152" in closure
    assert "GO:0006631\tGO:0006631" in closure

    query = write(tmp_path / "query.txt", "id\nHMDB00220\nCHEBI:29\nkegg:C00031\nCHEBI:999\nnonsense\n")
    report = tmp_path / "report.tsv"
    canonical_out = tmp_path / "query.canonical.txt"
    assert normalize_chemicals.main([str(query), "--prep-dir", str(out), "--out", str(canonical_out), "--report", str(report)]) == 0
    assert canonical_out.read_text().split() == ["CHEBI:21", "CHEBI:40"]
    with report.open() as handle:
        statuses = {row["input_id"]: row["status"] for row in csv.DictReader(handle, delimiter="\t")}
    assert statuses == {
        "HMDB00220": "xref",
        "CHEBI:29": "secondary",
        "kegg:C00031": "xref",
        "CHEBI:999": "unknown_chebi",
        "nonsense": "unparsed",
    }


def test_rollup_limit_and_currency_flags(tmp_path):
    out = tmp_path / "out"
    args = build_sources(tmp_path) + ["--out-dir", str(out), "--max-rollup-class-size", "0", "--keep-currency"]
    assert chem2go.main(args) == 0
    assert ("CHEBI:15377", "GO:0004467") in read_pairs(out / "annotations.rhea.tsv")
