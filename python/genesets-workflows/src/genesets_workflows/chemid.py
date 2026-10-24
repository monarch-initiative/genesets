"""Chemical identifier parsing and normalization.

Metabolomics results arrive under many identifier schemes (ChEBI in several
spellings, HMDB with 5- or 7-digit accessions, KEGG COMPOUND, LIPID MAPS, CAS).
GO and Rhea refer to chemicals by ChEBI, usually the major microspecies at
pH 7.3. This module turns a raw identifier into a CURIE and, given the
``normalization.tsv`` written by ``prepare-chem2go``, into the canonical ChEBI
class used by the chemical-to-GO annotation tables.

>>> parse_chemical_id("chebi:17234")
'CHEBI:17234'
>>> parse_chemical_id("http://purl.obolibrary.org/obo/CHEBI_17234")
'CHEBI:17234'
>>> parse_chemical_id("HMDB00122")
'hmdb:HMDB0000122'
>>> parse_chemical_id("KEGG:C00031")
'kegg.compound:C00031'
>>> parse_chemical_id("LMFA01010001")
'lipidmaps:LMFA01010001'
>>> parse_chemical_id("50-99-7")
'cas:50-99-7'
>>> parse_chemical_id("glucose") is None
True
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path


CHEBI_URI = re.compile(r"^https?://purl\.obolibrary\.org/obo/CHEBI_(\d+)$", re.IGNORECASE)
CHEBI_CURIE = re.compile(r"^(?:CHEBI)[:_](\d+)$", re.IGNORECASE)
HMDB = re.compile(r"^(?:HMDB:)?HMDB(\d{5,7})$", re.IGNORECASE)
KEGG_COMPOUND = re.compile(r"^(?:KEGG(?:\.COMPOUND)?:|CPD:)?(C\d{5})$", re.IGNORECASE)
LIPIDMAPS = re.compile(r"^(?:LIPIDMAPS:)?(LM[A-Z]{2}[0-9A-Z]{8,10})$", re.IGNORECASE)
CAS = re.compile(r"^(?:CAS:)?(\d{2,7}-\d{2}-\d)$", re.IGNORECASE)
PREFIXED = re.compile(r"^([A-Za-z][A-Za-z0-9._-]*):(\S+)$")

#: Lower-case CURIE prefixes for the cross-reference sources we index.
XREF_PREFIXES = ("hmdb", "kegg.compound", "lipidmaps", "metacyc.compound", "drugbank", "cas")


def parse_chemical_id(raw: str) -> str | None:
    """Return a normalized CURIE for a chemical identifier, or None.

    ChEBI ids become ``CHEBI:<n>``; cross-references use the lower-case
    bioregistry prefixes in :data:`XREF_PREFIXES`. HMDB accessions are padded
    to the current 7-digit form. Bare integers are not accepted because they
    are ambiguous across databases.

    >>> parse_chemical_id(" CHEBI_15377 ")
    'CHEBI:15377'
    >>> parse_chemical_id("hmdb:HMDB0000122")
    'hmdb:HMDB0000122'
    >>> parse_chemical_id("metacyc.compound:CPD-8756")
    'metacyc.compound:CPD-8756'
    >>> parse_chemical_id("17234") is None
    True
    """
    value = raw.strip()
    if not value:
        return None
    for pattern in (CHEBI_URI, CHEBI_CURIE):
        match = pattern.match(value)
        if match:
            return f"CHEBI:{int(match.group(1))}"
    match = HMDB.match(value)
    if match:
        return f"hmdb:HMDB{int(match.group(1)):07d}"
    match = KEGG_COMPOUND.match(value)
    if match:
        return f"kegg.compound:{match.group(1).upper()}"
    match = LIPIDMAPS.match(value)
    if match:
        return f"lipidmaps:{match.group(1).upper()}"
    match = CAS.match(value)
    if match:
        return f"cas:{match.group(1)}"
    match = PREFIXED.match(value)
    if match and match.group(1).lower() in XREF_PREFIXES:
        return f"{match.group(1).lower()}:{match.group(2)}"
    return None


@dataclass
class Normalized:
    input_id: str
    curie: str | None
    canonical_id: str | None
    status: str
    alternatives: list[str] = field(default_factory=list)


@dataclass
class ChemicalNormalizer:
    """Map raw identifiers to canonical ChEBI classes.

    ``to_canonical`` maps every known ChEBI id (primary, secondary/alt, and
    every conjugate acid/base or tautomer) to its canonical representative.
    ``xrefs`` maps cross-reference CURIEs to the ChEBI ids that cite them.
    """

    to_canonical: dict[str, str]
    methods: dict[str, str]
    xrefs: dict[str, list[str]]
    names: dict[str, str]
    annotated: set[str] = field(default_factory=set)

    @classmethod
    def load(cls, prep_dir: Path, extra_xrefs: Path | None = None) -> "ChemicalNormalizer":
        """Load ``prepare-chem2go`` outputs.

        ``extra_xrefs`` is an optional two-column TSV (``xref``, ``chebi_id``)
        that supplements ChEBI's own cross-references, e.g. curated fixes or a
        UniChem/MetaNetX export; its rows are added to, not replacing, ChEBI's.
        """
        to_canonical: dict[str, str] = {}
        methods: dict[str, str] = {}
        with (prep_dir / "normalization.tsv").open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                to_canonical[row["chebi_id"]] = row["canonical_id"]
                methods[row["chebi_id"]] = row["method"]
        xrefs: dict[str, list[str]] = {}
        for xref_path in (prep_dir / "xrefs.tsv", extra_xrefs):
            if xref_path is None or not xref_path.exists():
                continue
            with xref_path.open(newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    xref = parse_chemical_id(row["xref"]) or row["xref"]
                    chebi_id = parse_chemical_id(row["chebi_id"]) or row["chebi_id"]
                    xrefs.setdefault(xref, []).append(chebi_id)
        names: dict[str, str] = {}
        annotated: set[str] = set()
        with (prep_dir / "chemicals.tsv").open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                names[row["chemical_id"]] = row["name"]
                if row.get("annotated") == "true":
                    annotated.add(row["chemical_id"])
        return cls(to_canonical, methods, xrefs, names, annotated)

    def normalize(self, raw: str) -> Normalized:
        """Normalize one identifier.

        ``status`` is one of ``canonical`` (already the representative),
        ``merged`` (a conjugate acid/base, tautomer or pH 7.3 form),
        ``secondary`` (an alt_id), ``xref``, ``xref_ambiguous``,
        ``unknown_chebi``, ``unmapped_xref`` or ``unparsed``.
        """
        curie = parse_chemical_id(raw)
        if curie is None:
            return Normalized(raw, None, None, "unparsed")
        if curie.startswith("CHEBI:"):
            canonical = self.to_canonical.get(curie)
            if canonical is None:
                return Normalized(raw, curie, None, "unknown_chebi")
            method = self.methods[curie]
            if canonical == curie:
                status = "canonical"
            elif method == "alt_id":
                status = "secondary"
            else:
                status = "merged"
            return Normalized(raw, curie, canonical, status)
        targets = sorted(
            {self.to_canonical[c] for c in self.xrefs.get(curie, []) if c in self.to_canonical},
            key=lambda c: (c not in self.annotated, chebi_sort_key(c)),
        )
        if not targets:
            return Normalized(raw, curie, None, "unmapped_xref")
        status = "xref" if len(targets) == 1 else "xref_ambiguous"
        # Prefer a target that carries annotations, then the lowest ChEBI
        # number, usually the established, least-specific entry (e.g.
        # D-glucose rather than one anomer).
        return Normalized(raw, curie, targets[0], status, targets[1:])


def chebi_sort_key(curie: str) -> tuple[int, str]:
    """Sort ChEBI CURIEs numerically.

    >>> sorted(["CHEBI:100", "CHEBI:9"], key=chebi_sort_key)
    ['CHEBI:9', 'CHEBI:100']
    """
    prefix, _, local = curie.partition(":")
    return (int(local), prefix) if local.isdigit() else (1 << 62, curie)
