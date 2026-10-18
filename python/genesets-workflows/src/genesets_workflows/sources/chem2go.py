"""Prepare chemical-to-GO annotation tables for metabolite set enrichment.

GO annotates gene products, not chemicals, so a measured metabolite reaches GO
only through an intermediary. This module materializes up to three such
bridges as ordinary genesets-rs annotation tables (``gene_id`` = canonical
ChEBI id):

``logical_definition``
    GO-plus class definitions and restrictions that name a ChEBI class, e.g.
    ``fatty acid metabolic process`` has_primary_input_or_output
    ``fatty acid``. A measured chemical is annotated to the term when the
    ChEBI class is the chemical itself or one of its ``is_a`` ancestors, so
    palmitate reaches fatty acid metabolism. Species-independent.

``rhea``
    Rhea reaction participants joined to GO molecular functions through
    ``rhea2go``. The same ``is_a`` roll-up applies, so a specific fatty acid
    reaches reactions written with a generic ``long-chain fatty acid``.

``enzyme`` (optional, needs ``--gaf``)
    Rhea participants joined to the UniProt enzymes catalysing the reaction,
    then to those enzymes' GO annotations in one organism's GAF (biological
    process by default). This is the only bridge that reaches pathway-level
    processes such as the tricarboxylic acid cycle, which GO itself rarely
    links to chemicals. Species-specific.

Currency metabolites (water, H+, ATP, NAD(P)H, CoA, ...) are dropped wherever a
chemical is one participant among many (Rhea reactions, GO ``has_participant``),
because they take part in a large share of all reactions. Primary
inputs/outputs keep them, so ``ATP metabolic process`` still reaches ATP.

Both bridges first collapse ChEBI classes that metabolomics cannot tell apart:
conjugate acids/bases, tautomers, and Rhea's pH 7.3 major-microspecies
mapping. Each collapsed group is represented by its pH 7.3 form when one
exists, which is also the form GO and Rhea mostly use.

Outputs, under ``--out-dir``:

- ``terms.tsv`` / ``closure.tsv``: GO terms and reflexive is_a+part_of closure
  from the same GO-plus release.
- ``annotations.<variant>.tsv``: ``gene_id``/``term_id`` tables, one per
  bridge plus ``all`` (and ``species_independent`` = logical_definition + rhea
  when the enzyme bridge is on).
- ``provenance.tsv``: one row per chemical/term/bridge with the ChEBI class it
  matched through and whether that was a roll-up.
- ``chemicals.tsv``: canonical chemicals with names and collapsed members.
- ``normalization.tsv``: every ChEBI id (primary, alt_id, member) to its
  canonical id; read by ``normalize-chemicals``.
- ``xrefs.tsv``: HMDB/KEGG/LIPID MAPS/MetaCyc/DrugBank/CAS to ChEBI.
- ``currency.tsv``: the currency metabolites dropped as participants.
- ``metadata.yaml``/``metadata.json``: sources, policy and counts.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import IO, Any, Callable, Iterable, Iterator

from genesets_workflows.chemid import XREF_PREFIXES, chebi_sort_key, parse_chemical_id
from genesets_workflows.runtime import utc_now
from genesets_workflows.sources.reactome_flat import download, file_metadata
from genesets_workflows.yaml_io import write_yaml_value


DEFAULT_OUT_DIR = Path("evals/chem2go/generated/current")

SOURCES = {
    "go_plus": "https://purl.obolibrary.org/obo/go/extensions/go-plus.json",
    "chebi": "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo.gz",
    "chebi_xrefs": "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/database_accession.tsv.gz",
    "chebi_xref_sources": "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/source.tsv.gz",
    "rhea_reactions": "https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz",
    "rhea2go": "https://ftp.expasy.org/databases/rhea/tsv/rhea2go.tsv",
    "rhea_ph73": "https://ftp.expasy.org/databases/rhea/tsv/chebi_pH7_3_mapping.tsv",
}

#: Fetched only when the enzyme bridge is requested with --gaf/--gaf-url.
ENZYME_SOURCES = {
    "rhea2uniprot": "https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv",
}

OBO = "http://purl.obolibrary.org/obo/"

#: GO-plus properties that tie a GO class to a chemical it acts on.
RELATIONS = {
    "RO:0004009": "has_primary_input",
    "RO:0004008": "has_primary_output",
    "RO:0004007": "has_primary_input_or_output",
    "RO:0002332": "regulates_levels_of",
    "RO:0002505": "has_intermediate",
    "RO:0000057": "has_participant",
    "RO:0002608": "process_has_causal_agent",
    "RO:0012000": "has_small_molecule_regulator",
    "RO:0012001": "has_small_molecule_activator",
    "RO:0002588": "results_in_assembly_of",
}
DEFAULT_RELATIONS = ("RO:0004009", "RO:0004008", "RO:0004007", "RO:0002332", "RO:0002505", "RO:0000057")

#: Link kinds by GO ancestor, in priority order. A term takes the first kind
#: whose root is in its closure.
LINK_KIND_ROOTS = (
    ("metabolism", ("GO:0008152",)),
    ("transport", ("GO:0006810", "GO:0005215")),
    ("catalysis", ("GO:0003824",)),
    ("homeostasis", ("GO:0042592",)),
    ("signaling", ("GO:0023052", "GO:0038023")),
    ("binding", ("GO:0005488",)),
    ("response", ("GO:0050896",)),
)
DEFAULT_LINK_KINDS = ("metabolism", "transport", "catalysis", "homeostasis")

#: Currency metabolites as pH 7.3 ChEBI ids (names from Rhea's chebiId_name).
CURRENCY = {
    "CHEBI:15377": "H2O",
    "CHEBI:15378": "H(+)",
    "CHEBI:15379": "O2",
    "CHEBI:16526": "CO2",
    "CHEBI:17544": "hydrogencarbonate",
    "CHEBI:16240": "H2O2",
    "CHEBI:18421": "superoxide",
    "CHEBI:28938": "NH4(+)",
    "CHEBI:43474": "phosphate",
    "CHEBI:33019": "diphosphate",
    "CHEBI:30616": "ATP",
    "CHEBI:456216": "ADP",
    "CHEBI:456215": "AMP",
    "CHEBI:37565": "GTP",
    "CHEBI:58189": "GDP",
    "CHEBI:58115": "GMP",
    "CHEBI:46398": "UTP",
    "CHEBI:58223": "UDP",
    "CHEBI:37563": "CTP",
    "CHEBI:58069": "CDP",
    "CHEBI:61402": "ITP",
    "CHEBI:57540": "NAD(+)",
    "CHEBI:57945": "NADH",
    "CHEBI:58349": "NADP(+)",
    "CHEBI:57783": "NADPH",
    "CHEBI:57692": "FAD",
    "CHEBI:58307": "FADH2",
    "CHEBI:58210": "FMN",
    "CHEBI:57618": "FMNH2",
    "CHEBI:13193": "A (electron acceptor)",
    "CHEBI:17499": "AH2 (reduced acceptor)",
    "CHEBI:57287": "CoA",
    "CHEBI:57288": "acetyl-CoA",
    "CHEBI:59789": "S-adenosyl-L-methionine",
    "CHEBI:57856": "S-adenosyl-L-homocysteine",
    "CHEBI:29101": "Na(+)",
    "CHEBI:29103": "K(+)",
    "CHEBI:17996": "chloride",
    "CHEBI:18420": "Mg(2+)",
    "CHEBI:29108": "Ca(2+)",
    "CHEBI:29035": "Mn(2+)",
    "CHEBI:29033": "Fe(2+)",
    "CHEBI:29034": "Fe(3+)",
    "CHEBI:29105": "Zn(2+)",
}

CHARGED_NAME = re.compile(r"\(\d*[+-]\)|\b(?:anion|cation|zwitterion)\b")

#: Prefixes whose presence marks a ChEBI class as a measurable metabolite for
#: the default annotation universe.
METABOLITE_XREF_PREFIXES = ("hmdb", "kegg.compound", "lipidmaps", "metacyc.compound")

MERGE_RELATIONS = {
    "RO:0018033": "conjugate",  # is conjugate base of
    "RO:0018034": "conjugate",  # is conjugate acid of
    "RO:0018036": "tautomer",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="genesets-workflows prepare-chem2go",
        description="Prepare chemical-to-GO annotation tables for metabolite enrichment.",
    )
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR, type=Path)
    for key, url in {**SOURCES, **ENZYME_SOURCES}.items():
        flag = key.replace("_", "-")
        parser.add_argument(
            f"--{flag}",
            type=Path,
            help=f"Local copy of the {key} source; downloaded from {url} when omitted.",
        )
    parser.add_argument(
        "--gaf",
        type=Path,
        help=(
            "GO annotation file (GAF 2.x, UniProtKB ids) for one organism; enables the "
            "species-specific enzyme bridge."
        ),
    )
    parser.add_argument(
        "--gaf-url",
        help="Download the GAF for the enzyme bridge, e.g. http://current.geneontology.org/annotations/goa_human.gaf.gz",
    )
    parser.add_argument(
        "--enzyme-aspects",
        default="P",
        help="GAF aspects (P, F, C) carried through the enzyme bridge.",
    )
    parser.add_argument(
        "--relations",
        default=",".join(DEFAULT_RELATIONS),
        help="Comma-separated GO-plus properties (RO CURIEs) for the logical_definition bridge.",
    )
    parser.add_argument(
        "--link-kinds",
        default=",".join(DEFAULT_LINK_KINDS),
        help=(
            "Comma-separated GO link kinds kept in logical_definition: "
            + ", ".join(kind for kind, _ in LINK_KIND_ROOTS)
            + ", other."
        ),
    )
    parser.add_argument(
        "--universe",
        type=Path,
        help=(
            "Optional list of chemical ids (any supported scheme) to annotate, e.g. "
            "the platform's measured metabolites. By default: chemicals linked "
            "directly by a bridge, Rhea participants, and chemicals with an HMDB, "
            "KEGG, LIPID MAPS or MetaCyc cross-reference."
        ),
    )
    parser.add_argument(
        "--max-rollup-class-size",
        type=int,
        default=2500,
        help=(
            "Roll a chemical up only through ChEBI classes with at most this many "
            "descendant chemicals (0 = no limit). Broader GO terms are still reached "
            "through the GO closure."
        ),
    )
    parser.add_argument(
        "--keep-currency",
        action="store_true",
        help="Keep currency metabolites in the rhea bridge.",
    )
    parser.add_argument("--force-download", action="store_true")
    return parser.parse_args(argv)


# --------------------------------------------------------------------------
# Readers


def open_text(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8", newline="")


def obo_curie(uri: str) -> str:
    """Contract an OBO PURL to a CURIE.

    >>> obo_curie("http://purl.obolibrary.org/obo/GO_0008152")
    'GO:0008152'
    >>> obo_curie("http://purl.obolibrary.org/obo/go#foo")
    'http://purl.obolibrary.org/obo/go#foo'
    """
    if uri.startswith(OBO):
        local = uri[len(OBO):]
        prefix, sep, rest = local.partition("_")
        if sep and prefix.isalpha() and prefix.isupper():
            return f"{prefix}:{rest}"
    return uri


class Chebi:
    """The parts of ChEBI the bridges need."""

    def __init__(self) -> None:
        self.names: dict[str, str] = {}
        self.is_a: dict[str, set[str]] = defaultdict(set)
        self.merge_edges: list[tuple[str, str, str]] = []
        self.alt_ids: dict[str, str] = {}
        self.version = ""

    @classmethod
    def from_obo(cls, handle: Iterable[str]) -> "Chebi":
        """Parse ChEBI OBO, dropping obsolete stanzas.

        >>> text = '''data-version: chebi/1
        ... [Term]
        ... id: CHEBI:2
        ... name: acid
        ... alt_id: CHEBI:9
        ... is_a: CHEBI:1
        ... relationship: RO:0018034 CHEBI:3
        ...
        ... [Term]
        ... id: CHEBI:3
        ... name: anion
        ... [Term]
        ... id: CHEBI:4
        ... name: gone
        ... is_obsolete: true
        ... '''
        >>> chebi = Chebi.from_obo(text.splitlines())
        >>> sorted(chebi.names), chebi.alt_ids, chebi.merge_edges
        (['CHEBI:2', 'CHEBI:3'], {'CHEBI:9': 'CHEBI:2'}, [('CHEBI:2', 'CHEBI:3', 'conjugate')])
        """
        chebi = cls()
        stanza: dict[str, Any] | None = None

        def commit() -> None:
            if not stanza or "id" not in stanza or stanza.get("obsolete"):
                return
            term_id = stanza["id"]
            chebi.names[term_id] = stanza.get("name", "")
            chebi.is_a[term_id].update(stanza["is_a"])
            for alt in stanza["alt_ids"]:
                chebi.alt_ids[alt] = term_id
            for target, kind in stanza["merge"]:
                chebi.merge_edges.append((term_id, target, kind))

        for raw in handle:
            line = raw.strip()
            if line.startswith("["):
                commit()
                stanza = {"is_a": set(), "alt_ids": [], "merge": []} if line == "[Term]" else None
                continue
            if stanza is None:
                if line.startswith("data-version:"):
                    chebi.version = line.partition(":")[2].strip()
                continue
            key, _, value = line.partition(": ")
            if key == "id":
                stanza["id"] = value
            elif key == "name":
                stanza["name"] = value
            elif key == "is_a":
                stanza["is_a"].add(value.split()[0])
            elif key == "alt_id":
                stanza["alt_ids"].append(value.split()[0])
            elif key == "is_obsolete" and value == "true":
                stanza["obsolete"] = True
            elif key == "relationship":
                parts = value.split()
                if len(parts) >= 2 and parts[0] in MERGE_RELATIONS:
                    stanza["merge"].append((parts[1], MERGE_RELATIONS[parts[0]]))
        commit()
        for term_id in list(chebi.is_a):
            chebi.is_a[term_id] = {p for p in chebi.is_a[term_id] if p in chebi.names}
        return chebi


def read_ph73(path: Path) -> list[tuple[str, str]]:
    pairs = []
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pairs.append((f"CHEBI:{row['CHEBI']}", f"CHEBI:{row['CHEBI_PH7_3']}"))
    return pairs


def read_rhea_equations(path: Path) -> dict[str, list[str]]:
    """Parse the Rhea flat file into reaction id -> ChEBI participants.

    >>> text = '''ENTRY       RHEA:10000
    ... EQUATION    CHEBI:16459 + CHEBI:15377 = CHEBI:31011 + 2 CHEBI:28938
    ... ///
    ... ENTRY       RHEA:10010
    ... EQUATION    POLYMER:1 + CHEBI:15377 = GENERIC:2
    ... ///'''
    >>> parse_rhea_equations(text.splitlines())
    {'10000': ['CHEBI:15377', 'CHEBI:16459', 'CHEBI:28938', 'CHEBI:31011'], '10010': ['CHEBI:15377']}
    """
    with open_text(path) as handle:
        return parse_rhea_equations(handle)


def parse_rhea_equations(lines: Iterable[str]) -> dict[str, list[str]]:
    equations: dict[str, list[str]] = {}
    entry = None
    for line in lines:
        if line.startswith("ENTRY"):
            entry = line.split()[1].removeprefix("RHEA:")
        elif line.startswith("EQUATION") and entry:
            equations[entry] = sorted(set(re.findall(r"CHEBI:\d+", line)))
        elif line.startswith("///"):
            entry = None
    return equations


def read_rhea2go(path: Path) -> list[tuple[str, str]]:
    with open_text(path) as handle:
        return [(row["RHEA_ID"], row["ID"]) for row in csv.DictReader(handle, delimiter="\t")]


def read_gaf(path: Path, aspects: set[str]) -> Iterator[tuple[str, str, str]]:
    r"""Yield (UniProtKB accession, GO term, evidence) from a GAF, skipping NOT.

    >>> line = "UniProtKB\tP1\tG\tinvolved_in\tGO:1\tPMID:1\tIDA\t\tP\n"
    >>> list(parse_gaf([line, line.replace("involved_in", "NOT|involved_in")], {"P"}))
    [('P1', 'GO:1', 'IDA')]
    """
    with open_text(path) as handle:
        yield from parse_gaf(handle, aspects)


def parse_gaf(lines: Iterable[str], aspects: set[str]) -> Iterator[tuple[str, str, str]]:
    for line in lines:
        if line.startswith("!"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[0] != "UniProtKB" or fields[8] not in aspects:
            continue
        if "NOT" in fields[3].split("|"):
            continue
        yield fields[1], fields[4], fields[6]


def read_rhea2uniprot(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (Rhea master reaction id, UniProtKB accession)."""
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            yield row["MASTER_ID"], row["ID"]


def read_xrefs(accessions: Path, sources: Path) -> Iterator[tuple[str, str]]:
    """Yield (xref CURIE, ChEBI id) for the prefixes in XREF_PREFIXES."""
    with open_text(sources) as handle:
        prefix_for = {row["id"]: row["prefix"] for row in csv.DictReader(handle, delimiter="\t")}
    with open_text(accessions) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            prefix = "cas" if row["type"] == "CAS" else prefix_for.get(row["source_id"], "")
            if prefix not in XREF_PREFIXES:
                continue
            curie = parse_chemical_id(f"{prefix}:{row['accession_number']}")
            if curie and not curie.startswith("CHEBI:"):
                yield curie, f"CHEBI:{row['compound_id']}"


class GoPlus:
    def __init__(self) -> None:
        self.names: dict[str, str] = {}
        self.parents: dict[str, set[str]] = defaultdict(set)
        # (GO term, RO property, ChEBI class)
        self.chemical_links: set[tuple[str, str, str]] = set()
        self.version = ""

    @classmethod
    def from_obograph(cls, path: Path) -> "GoPlus":
        with open_text(path) as handle:
            graph = json.load(handle)["graphs"][0]
        go = cls()
        go.version = graph.get("meta", {}).get("version", "")
        for node in graph["nodes"]:
            curie = obo_curie(node["id"])
            if not curie.startswith("GO:") or node.get("type") != "CLASS":
                continue
            if node.get("meta", {}).get("deprecated"):
                continue
            go.names[curie] = node.get("lbl", "")
        for edge in graph["edges"]:
            sub, pred, obj = obo_curie(edge["sub"]), edge["pred"], obo_curie(edge["obj"])
            if sub not in go.names:
                continue
            if obj.startswith("CHEBI:"):
                go.chemical_links.add((sub, obo_curie(pred), obj))
            elif obj in go.names and (pred == "is_a" or obo_curie(pred) == "BFO:0000050"):
                go.parents[sub].add(obj)
        for axiom in graph.get("logicalDefinitionAxioms", []):
            defined = obo_curie(axiom["definedClassId"])
            if defined not in go.names:
                continue
            for restriction in axiom.get("restrictions", []):
                filler = obo_curie(restriction["fillerId"])
                if filler.startswith("CHEBI:"):
                    go.chemical_links.add((defined, obo_curie(restriction["propertyId"]), filler))
        return go


# --------------------------------------------------------------------------
# Graph helpers


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, item: str) -> str:
        self.parent.setdefault(item, item)
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != root:
            self.parent[item], item = root, self.parent[item]
        return root

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb

    def groups(self) -> dict[str, list[str]]:
        groups: dict[str, list[str]] = defaultdict(list)
        for item in self.parent:
            groups[self.find(item)].append(item)
        return groups


def ancestors(node: str, parents: dict[str, set[str]], memo: dict[str, frozenset[str]]) -> frozenset[str]:
    """Reflexive ancestors; tolerates cycles.

    >>> sorted(ancestors("a", {"a": {"b"}, "b": {"c", "a"}}, {}))
    ['a', 'b', 'c']
    """
    if node in memo:
        return memo[node]
    seen = {node}
    stack = list(parents.get(node, ()))
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        if current in memo:
            seen.update(memo[current])
            continue
        seen.add(current)
        stack.extend(parents.get(current, ()))
    result = frozenset(seen)
    memo[node] = result
    return result


def collapse_chemicals(chebi: Chebi, ph73: list[tuple[str, str]]) -> tuple[dict[str, str], dict[str, str]]:
    """Group indistinguishable ChEBI classes; return (id -> canonical, id -> method).

    The canonical member is the pH 7.3 major microspecies when the group has
    one (the most frequent mapping target), otherwise the lowest ChEBI id.

    >>> chebi = Chebi()
    >>> chebi.names = {"CHEBI:1": "acid", "CHEBI:2": "anion", "CHEBI:3": "other"}
    >>> chebi.merge_edges = [("CHEBI:1", "CHEBI:2", "conjugate")]
    >>> chebi.alt_ids = {"CHEBI:9": "CHEBI:1"}
    >>> canonical, method = collapse_chemicals(chebi, [("CHEBI:1", "CHEBI:2")])
    >>> canonical == {"CHEBI:1": "CHEBI:2", "CHEBI:2": "CHEBI:2", "CHEBI:3": "CHEBI:3", "CHEBI:9": "CHEBI:2"}
    True
    >>> method["CHEBI:1"], method["CHEBI:9"]
    ('conjugate', 'alt_id')

    Without pH 7.3 data the neutral-looking name wins over a lower charged id.

    >>> chebi.names = {"CHEBI:5": "amino-acid cation", "CHEBI:6": "amino acid"}
    >>> chebi.merge_edges, chebi.alt_ids = [("CHEBI:6", "CHEBI:5", "conjugate")], {}
    >>> collapse_chemicals(chebi, [])[0]["CHEBI:5"]
    'CHEBI:6'
    """
    uf = UnionFind()
    method: dict[str, str] = {}
    for term_id in chebi.names:
        uf.find(term_id)
    for a, b, kind in chebi.merge_edges:
        if a in chebi.names and b in chebi.names:
            uf.union(a, b)
            for member in (a, b):
                method.setdefault(member, kind)
    targets: Counter[str] = Counter()
    for a, b in ph73:
        if a in chebi.names and b in chebi.names:
            if a != b:
                uf.union(a, b)
                method.setdefault(a, "ph7_3")
            targets[b] += 1
    canonical: dict[str, str] = {}
    for members in uf.groups().values():
        preferred = [m for m in members if targets[m]]
        if preferred:
            rep = min(preferred, key=lambda m: (-targets[m], chebi_sort_key(m)))
        else:
            rep = min(members, key=lambda m: (bool(CHARGED_NAME.search(chebi.names.get(m, ""))), chebi_sort_key(m)))
        for member in members:
            canonical[member] = rep
    for term_id in chebi.names:
        if canonical[term_id] == term_id:
            method[term_id] = "canonical"
    for alt, primary in chebi.alt_ids.items():
        if primary in canonical and alt not in canonical:
            canonical[alt] = canonical[primary]
            method[alt] = "alt_id"
    return canonical, method


def quotient_parents(chebi: Chebi, canonical: dict[str, str]) -> dict[str, set[str]]:
    parents: dict[str, set[str]] = defaultdict(set)
    for child, direct in chebi.is_a.items():
        rep = canonical[child]
        parents[rep].update(canonical[p] for p in direct if canonical[p] != rep)
    return parents


def link_kind(term: str, go_parents: dict[str, set[str]], memo: dict[str, frozenset[str]]) -> str:
    closure = ancestors(term, go_parents, memo)
    for kind, roots in LINK_KIND_ROOTS:
        if any(root in closure for root in roots):
            return kind
    return "other"


# --------------------------------------------------------------------------
# Build


Index = dict[str, dict[str, str]]


def add_link(index: Index, chemical: str, term: str, detail: str) -> None:
    """Record chemical -> term, keeping the first detail seen as the example."""
    index.setdefault(chemical, {}).setdefault(term, detail)


def build_annotations(
    universe: Iterable[str],
    chem_parents: dict[str, set[str]],
    bridges: dict[str, Index],
    too_broad: frozenset[str] | set[str] = frozenset(),
) -> list[tuple[str, str, str, str, str, bool]]:
    """Join chemicals to GO through each bridge's ChEBI-class index.

    ``bridges[source][chebi_class]`` maps GO terms to an example detail. A
    chemical reaches every term indexed under itself or an ``is_a`` ancestor,
    except ancestors in ``too_broad`` (direct links are always kept). Rows are
    (chemical, term, source, detail, via, rollup), keeping the
    lowest-numbered ``via`` per chemical/term/source, with the chemical itself
    first.

    >>> parents = {"palmitate": {"fatty acid"}, "fatty acid": set()}
    >>> bridges = {"logical_definition": {"fatty acid": {"GO:fam": "has_primary_input_or_output"}}}
    >>> build_annotations(["palmitate"], parents, bridges)
    [('palmitate', 'GO:fam', 'logical_definition', 'has_primary_input_or_output', 'fatty acid', True)]
    >>> build_annotations(["palmitate"], parents, bridges, {"fatty acid"})
    []
    >>> build_annotations(["fatty acid"], parents, bridges, {"fatty acid"})[0][-1]
    False
    """
    memo: dict[str, frozenset[str]] = {}
    rows = []
    for chemical in universe:
        lineage = ancestors(chemical, chem_parents, memo)
        for source, index in bridges.items():
            best: dict[str, tuple[str, str]] = {}
            for via in sorted(lineage, key=lambda c: (c != chemical, chebi_sort_key(c))):
                if via != chemical and via in too_broad:
                    continue
                for term, detail in sorted(index.get(via, {}).items()):
                    best.setdefault(term, (detail, via))
            for term, (detail, via) in sorted(best.items()):
                rows.append((chemical, term, source, detail, via, via != chemical))
    return rows


def build_enzyme_index(
    rhea2uniprot: Path,
    gaf: Path,
    participants: Callable[[str], list[str]],
    go_terms: dict[str, str],
    aspects: set[str],
) -> tuple[Index, dict[str, Any]]:
    """Chemical -> GO terms annotated to an enzyme whose reaction involves it.

    Organism-specific: only enzymes present in the GAF contribute, so a
    human GAF gives human pathways. Detail is ``RHEA:<id>|UniProtKB:<acc>``.
    """
    enzyme_terms: dict[str, set[str]] = defaultdict(set)
    for accession, term, _evidence in read_gaf(gaf, aspects):
        if term in go_terms:
            enzyme_terms[accession].add(term)
    index: Index = {}
    enzymes_used: set[str] = set()
    reactions_used: set[str] = set()
    for reaction, accession in read_rhea2uniprot(rhea2uniprot):
        terms = enzyme_terms.get(accession)
        if not terms:
            continue
        enzymes_used.add(accession)
        reactions_used.add(reaction)
        detail = f"RHEA:{reaction}|UniProtKB:{accession}"
        for participant in participants(reaction):
            for term in terms:
                add_link(index, participant, term, detail)
    return index, {
        "gaf_annotated_proteins": len(enzyme_terms),
        "enzymes_with_rhea_reactions": len(enzymes_used),
        "rhea_reactions_with_enzymes": len(reactions_used),
        "aspects": sorted(aspects),
    }


def write_tsv(path: Path, header: list[str], rows: Iterable[Iterable[Any]]) -> int:
    count = 0
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def resolve_sources(args: argparse.Namespace) -> tuple[dict[str, Path], dict[str, Any]]:
    downloads = args.out_dir / "downloads"
    sources = dict(SOURCES)
    if args.gaf or args.gaf_url:
        sources.update(ENZYME_SOURCES)
        sources["gaf"] = args.gaf_url or str(args.gaf)
    paths: dict[str, Path] = {}
    provenance: dict[str, Any] = {}
    for key, url in sources.items():
        local = getattr(args, key, None)
        if local is not None:
            paths[key] = local
            provenance[key] = file_metadata(local, url)
        else:
            path = downloads / url.rsplit("/", 1)[-1]
            provenance[key] = download(url, path, args.force_download)
            paths[key] = path
    return paths, provenance


def read_universe(path: Path, canonical: dict[str, str], xrefs: dict[str, set[str]]) -> tuple[set[str], int]:
    chemicals: set[str] = set()
    unmapped = 0
    with path.open() as handle:
        for line in handle:
            raw = line.split("\t")[0].strip()
            if not raw or raw.startswith("#"):
                continue
            curie = parse_chemical_id(raw)
            ids = {curie} if curie and curie.startswith("CHEBI:") else xrefs.get(curie or "", set())
            mapped = {canonical[c] for c in ids if c in canonical}
            if mapped:
                chemicals.add(min(mapped, key=chebi_sort_key))
            else:
                unmapped += 1
    return chemicals, unmapped


def prepare(args: argparse.Namespace) -> dict[str, Any]:
    relations = [r.strip() for r in args.relations.split(",") if r.strip()]
    link_kinds = {k.strip() for k in args.link_kinds.split(",") if k.strip()}
    known_kinds = {kind for kind, _ in LINK_KIND_ROOTS} | {"other"}
    if not link_kinds <= known_kinds:
        raise SystemExit(f"unknown link kinds: {sorted(link_kinds - known_kinds)}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    paths, download_metadata = resolve_sources(args)

    with open_text(paths["chebi"]) as handle:
        chebi = Chebi.from_obo(handle)
    canonical, method = collapse_chemicals(chebi, read_ph73(paths["rhea_ph73"]))
    chem_parents = quotient_parents(chebi, canonical)
    groups: dict[str, list[str]] = defaultdict(list)
    for member, rep in canonical.items():
        if member in chebi.names:
            groups[rep].append(member)

    xrefs: dict[str, set[str]] = defaultdict(set)
    for xref, chebi_id in read_xrefs(paths["chebi_xrefs"], paths["chebi_xref_sources"]):
        if chebi_id in canonical:
            xrefs[xref].add(chebi_id)

    go = GoPlus.from_obograph(paths["go_plus"])
    go_memo: dict[str, frozenset[str]] = {}

    # Currency metabolites are dropped wherever a chemical is one participant
    # among many (has_participant, reactions); primary inputs/outputs keep
    # them, so "ATP metabolic process" still reaches ATP.
    currency = {canonical[c] for c in CURRENCY if c in canonical}
    skip_participant = set() if args.keep_currency else currency

    logical: Index = {}
    kind_counts: Counter[str] = Counter()
    relation_counts: Counter[str] = Counter()
    for term, prop, filler in sorted(go.chemical_links):
        if prop not in relations or filler not in canonical:
            continue
        if prop == "RO:0000057" and canonical[filler] in skip_participant:
            continue
        kind = link_kind(term, go.parents, go_memo)
        kind_counts[kind] += 1
        if kind not in link_kinds:
            continue
        relation_counts[prop] += 1
        add_link(logical, canonical[filler], term, f"{RELATIONS.get(prop, prop)}|{kind}")

    equations = read_rhea_equations(paths["rhea_reactions"])

    def participants(reaction: str) -> list[str]:
        found = {canonical[p] for p in equations.get(reaction, ()) if p in canonical}
        return sorted(found - skip_participant, key=chebi_sort_key)

    rhea: Index = {}
    participant_degree: Counter[str] = Counter()
    rhea_links = 0
    for reaction, term in read_rhea2go(paths["rhea2go"]):
        if term not in go.names:
            continue
        rhea_links += 1
        for participant in equations.get(reaction, ()):
            if participant in canonical:
                participant_degree[canonical[participant]] += 1
        for participant in participants(reaction):
            add_link(rhea, participant, term, f"RHEA:{reaction}")

    bridges: dict[str, Index] = {"logical_definition": logical, "rhea": rhea}
    enzyme_stats: dict[str, Any] = {}
    if "gaf" in paths:
        bridges["enzyme"], enzyme_stats = build_enzyme_index(
            paths["rhea2uniprot"], paths["gaf"], participants, go.names, set(args.enzyme_aspects)
        )

    if args.universe:
        universe, universe_unmapped = read_universe(args.universe, canonical, xrefs)
    else:
        universe_unmapped = 0
        universe = {c for index in bridges.values() for c in index} | {
            canonical[c]
            for xref, ids in xrefs.items()
            if xref.partition(":")[0] in METABOLITE_XREF_PREFIXES
            for c in ids
        }

    chem_memo: dict[str, frozenset[str]] = {}
    class_size: Counter[str] = Counter()
    for rep_id in groups:
        class_size.update(ancestors(rep_id, chem_parents, chem_memo))
    limit = args.max_rollup_class_size
    too_broad = {c for c, n in class_size.items() if limit and n > limit}

    rows = build_annotations(sorted(universe, key=chebi_sort_key), chem_parents, bridges, too_broad)

    out = args.out_dir
    write_tsv(out / "terms.tsv", ["term_id", "name"], sorted(go.names.items()))
    closure_rows = write_tsv(
        out / "closure.tsv",
        ["child", "ancestor"],
        ((t, a) for t in sorted(go.names) for a in sorted(ancestors(t, go.parents, go_memo))),
    )
    variants: dict[str, set[tuple[str, str]]] = {name: set() for name in bridges}
    for c, t, source, *_ in rows:
        variants[source].add((c, t))
    if "enzyme" in variants:
        variants["species_independent"] = variants["logical_definition"] | variants["rhea"]
    variants["all"] = set().union(*(variants[name] for name in bridges))
    variant_stats = {}
    for name, pairs in variants.items():
        ordered = sorted(pairs, key=lambda p: (chebi_sort_key(p[0]), p[1]))
        write_tsv(out / f"annotations.{name}.tsv", ["gene_id", "term_id"], ordered)
        variant_stats[name] = {
            "annotation_rows": len(pairs),
            "chemicals": len({c for c, _ in pairs}),
            "terms": len({t for _, t in pairs}),
        }
    write_tsv(
        out / "provenance.tsv",
        ["chemical_id", "chemical_name", "term_id", "term_name", "source", "detail", "via_chemical_id", "via_chemical_name", "rollup"],
        (
            (c, chebi.names.get(c, ""), t, go.names.get(t, ""), s, d, v, chebi.names.get(v, ""), str(r).lower())
            for c, t, s, d, v, r in rows
        ),
    )
    annotated = {c for c, _ in variants["all"]}
    write_tsv(
        out / "chemicals.tsv",
        ["chemical_id", "name", "members", "annotated"],
        (
            (rep, chebi.names.get(rep, ""), "|".join(sorted(members, key=chebi_sort_key)), str(rep in annotated).lower())
            for rep, members in sorted(groups.items(), key=lambda kv: chebi_sort_key(kv[0]))
        ),
    )
    write_tsv(
        out / "normalization.tsv",
        ["chebi_id", "canonical_id", "method"],
        ((c, canonical[c], method.get(c, "canonical")) for c in sorted(canonical, key=chebi_sort_key)),
    )
    xref_rows = write_tsv(
        out / "xrefs.tsv",
        ["xref", "chebi_id"],
        ((x, c) for x in sorted(xrefs) for c in sorted(xrefs[x], key=chebi_sort_key)),
    )
    write_tsv(
        out / "currency.tsv",
        ["chemical_id", "name", "rhea_reaction_links", "excluded"],
        (
            (c, chebi.names.get(c, ""), participant_degree[c], str(not args.keep_currency).lower())
            for c in sorted(currency, key=chebi_sort_key)
        ),
    )

    group_sizes = Counter(len(m) for m in groups.values())
    metadata = {
        "generated_at_utc": utc_now(),
        "source": "GO-plus chemical logical definitions and Rhea rhea2go, over collapsed ChEBI",
        "versions": {"go_plus": go.version, "chebi": chebi.version},
        "downloads": download_metadata,
        "policy": {
            "logical_definition_relations": {r: RELATIONS.get(r, r) for r in relations},
            "logical_definition_link_kinds": sorted(link_kinds),
            "chemical_rollup": "chemical annotated via itself and every is_a ancestor in collapsed ChEBI",
            "max_rollup_class_size": limit,
            "chemical_collapse": "conjugate acid/base, tautomer, Rhea pH 7.3 mapping; representative is the pH 7.3 form",
            "currency_dropped_as_participant": not args.keep_currency,
            "go_closure_relations": ["is_a", "part_of"],
            "go_closure_reflexive": True,
            "universe": str(args.universe) if args.universe else (
                "directly linked chemicals, Rhea participants, and chemicals with "
                + "/".join(METABOLITE_XREF_PREFIXES) + " cross-references"
            ),
            "background_warning": (
                "Use the chemicals your platform measured as background, normalized with "
                "normalize-chemicals; the annotated chemicals are not a valid metabolomics background."
            ),
        },
        "stats": {
            "chebi_classes": len(chebi.names),
            "canonical_chemicals": len(groups),
            "largest_collapsed_group": max(group_sizes) if group_sizes else 0,
            "go_terms": len(go.names),
            "go_closure_rows": closure_rows,
            "go_chemical_links": len(go.chemical_links),
            "go_chemical_links_by_kind": dict(sorted(kind_counts.items())),
            "go_chemical_links_used_by_relation": {RELATIONS.get(k, k): v for k, v in sorted(relation_counts.items())},
            "rhea2go_links": rhea_links,
            "enzyme_bridge": enzyme_stats or "disabled (no --gaf)",
            "xref_rows": xref_rows,
            "universe_size": len(universe),
            "rollup_classes_too_broad": len(too_broad),
            "universe_unmapped": universe_unmapped,
            "variants": variant_stats,
            "top_rhea_participants": [
                {"chemical_id": c, "name": chebi.names.get(c, ""), "links": n, "currency": c in currency}
                for c, n in participant_degree.most_common(25)
            ],
        },
    }
    with (out / "metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (out / "metadata.yaml").open("w") as handle:
        write_yaml_value(handle, metadata)
    return metadata


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    metadata = prepare(args)
    variants = metadata["stats"]["variants"]
    for name, stats in variants.items():
        print(
            f"{name}: {stats['annotation_rows']} annotations, "
            f"{stats['chemicals']} chemicals, {stats['terms']} GO terms"
        )
    print(f"Wrote chemical-to-GO tables under {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
