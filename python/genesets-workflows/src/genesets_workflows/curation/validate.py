"""Build and run the three validators over an interpretation file."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from genesets_workflows.curation.model import GeneSetInterpretation

TARGET_CLASS = "GeneSetInterpretation"

# NCBITaxon ids are stable and its adapter is large; skip it in the obsolescence
# sweep. Every other configured ontology prefix is checked.
OBSOLESCENCE_SKIP_PREFIXES = frozenset({"NCBITaxon"})


def term_ids_by_prefix(
    interp: GeneSetInterpretation,
    skip_prefixes: frozenset[str] = OBSOLESCENCE_SKIP_PREFIXES,
) -> dict[str, set[str]]:
    """Map ontology prefix -> set of term ids asserted in contexts/associations.

    >>> from genesets_workflows.curation import model
    >>> interp = model.GeneSetInterpretation(
    ...     gene_set_id="LIT:X", gene_set_name="X",
    ...     taxon=model.Term(id="NCBITaxon:9606", label="Homo sapiens"),
    ...     contexts=[model.BiologicalContext(
    ...         term=model.Term(id="CL:0000129", label="microglial cell"),
    ...         context_type="cell_type")],
    ...     associations=[model.TermAssociation(
    ...         term=model.Term(id="GO:0006909", label="phagocytosis"),
    ...         seed_source="enrichment_recovered")])
    >>> term_ids_by_prefix(interp) == {"CL": {"CL:0000129"}, "GO": {"GO:0006909"}}
    True
    """
    by_prefix: dict[str, set[str]] = {}
    terms = [c.term for c in interp.contexts] + [a.term for a in interp.associations]
    for term in terms:
        tid = term.id
        if not tid or ":" not in tid:
            continue
        prefix = tid.split(":", 1)[0]
        if prefix in skip_prefixes:
            continue
        by_prefix.setdefault(prefix, set()).add(tid)
    return by_prefix


def obsolete_term_ids(path: Path, oak_config: Path) -> list[str]:
    """Return sorted obsolete ontology term ids asserted in the interpretation.

    The term-validator checks id+label match but NOT obsolescence (the OAK
    sqlite builds retain labels for obsolete classes, so an obsolete id with its
    old label still passes). This sweeps each term id against its ontology's
    OAK ``obsoletes()`` set, using the same adapters as the term-validator.
    Prefixes absent from ``oak_config`` are skipped (the term-validator skips
    them too).
    """
    import yaml
    from oaklib import get_adapter

    from genesets_workflows.curation import model

    adapters = (yaml.safe_load(Path(oak_config).read_text()) or {}).get(
        "ontology_adapters", {}
    )
    interp = model.load_interpretation(path)
    found: list[str] = []
    for prefix, ids in term_ids_by_prefix(interp).items():
        adapter_str = adapters.get(prefix)
        if not adapter_str:
            continue
        obsoletes = set(get_adapter(adapter_str).obsoletes())
        found.extend(tid for tid in ids if tid in obsoletes)
    return sorted(found)


def build_commands(
    path: Path,
    *,
    schema: Path,
    oak_config: Path,
    cache_dir: Path,
    skip_references: bool,
    ref_config: Path | None = None,
) -> list[list[str]]:
    commands = [
        [
            "linkml-validate",
            "--schema", str(schema),
            "--target-class", TARGET_CLASS,
            str(path),
        ],
        [
            "linkml-term-validator",
            "validate-data", str(path),
            "-s", str(schema),
            "-t", TARGET_CLASS,
            "--labels",
            "-c", str(oak_config),
        ],
    ]
    if not skip_references:
        reference_command = [
            "linkml-reference-validator",
            "validate", "data", str(path),
            "--schema", str(schema),
            "--target-class", TARGET_CLASS,
            "--cache-dir", str(cache_dir),
        ]
        if ref_config is not None:
            reference_command += ["--config", str(ref_config)]
        commands.append(reference_command)
    return commands


def validate_file(
    path: Path,
    *,
    schema: Path,
    oak_config: Path,
    cache_dir: Path,
    skip_references: bool = False,
    ref_config: Path | None = None,
) -> int:
    """Run each validator; return 0 if all succeed, else the first non-zero code."""
    for command in build_commands(
        path,
        schema=schema,
        oak_config=oak_config,
        cache_dir=cache_dir,
        skip_references=skip_references,
        ref_config=ref_config,
    ):
        print(f"$ {' '.join(command)}")
        completed = subprocess.run(command)
        if completed.returncode != 0:
            return completed.returncode

    print("$ obsolescence check (OAK obsoletes)")
    obsolete = obsolete_term_ids(path, oak_config)
    if obsolete:
        for tid in obsolete:
            print(f"  ERROR: obsolete term id in use: {tid}")
        return 1
    return 0
