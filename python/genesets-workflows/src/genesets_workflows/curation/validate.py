"""Build and run the validators over an interpretation file."""

from __future__ import annotations

import subprocess
from pathlib import Path

TARGET_CLASS = "GeneSetInterpretation"


def ontology_drift_errors(path: Path, oak_config: Path) -> tuple[list[str], list[str]]:
    """Return (errors, warnings) for ontology drift in one interpretation.

    The term-validator checks id+label but NOT obsolescence (the OAK sqlite
    builds keep labels for obsolete classes, so an obsolete id with its old
    label still passes), and it accepts synonyms, so a GO relabel passes
    silently. Every obsolete, merged or unknown id is an error, because GOA
    never annotates to it and the term could never be recovered. A relabelled
    live term is only a warning, so a routine GO rename does not break the
    build. Each message names the fix (`just curate-migrate`, or a curator
    `--map` decision).
    """
    from genesets_workflows.curation import ontology_drift as od

    prefixes = {tid.split(":", 1)[0] for tid, _ in od.asserted_terms(path)}
    drifts = od.scan([Path(path)], od.load_views(oak_config, prefixes))
    errors, warnings = [], []
    for d in drifts:
        if d.kind == "label_changed":
            warnings.append(f"{d.describe()} (fix: just curate-migrate)")
        elif d.auto_fixable:
            errors.append(f"{d.describe()} (fix: just curate-migrate)")
        else:
            hint = "curator decision: genesets-workflows curate drift --map OLD=NEW"
            errors.append(f"{d.describe()} ({hint}){' reason: ' + d.reason if d.reason else ''}")
    return errors, warnings


NOISE_CATEGORIES = {"nonspecific", "false_association"}


def noise_label_errors(interp) -> list[str]:
    """Consistency of the noise labels used for precision scoring.

    `artifact_cause` explains why a term is noise, so it is only meaningful on
    `nonspecific` / `false_association` terms. A set that declares
    `hit_adjudication` claims every pooled hit was labelled, so none of its
    associations may be left without a category.
    """
    errors = []
    for assoc in interp.associations:
        if assoc.artifact_cause and assoc.category not in NOISE_CATEGORIES:
            errors.append(
                f"{assoc.term.id}: artifact_cause set on category {assoc.category!r}; "
                "it is only for nonspecific / false_association"
            )
    if interp.hit_adjudication is not None:
        if not interp.hit_adjudication.pooled_variants:
            errors.append("hit_adjudication has no pooled_variants")
        for assoc in interp.associations:
            if not assoc.category:
                errors.append(
                    f"{assoc.term.id}: no category, but hit_adjudication says every hit was labelled"
                )
    return errors


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

    print("$ ontology drift check (obsolete / merged / relabelled terms)")
    drift_errors, drift_warnings = ontology_drift_errors(path, oak_config)
    for warning in drift_warnings:
        print(f"  WARNING: {warning}")
    if drift_errors:
        for error in drift_errors:
            print(f"  ERROR: {error}")
        return 1

    print("$ ranking check (file, checksum, leading edges, directions)")
    from genesets_workflows.curation import model
    from genesets_workflows.curation.ranking import ranking_errors

    interp = model.load_interpretation(path)
    errors = ranking_errors(interp, Path(path).parent)
    print("$ noise-label check (artifact_cause, hit_adjudication)")
    errors += noise_label_errors(interp)
    for error in errors:
        print(f"  ERROR: {error}")
    return 1 if errors else 0
