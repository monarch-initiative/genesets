"""Typer CLI for the `genesets-workflows curate` command group."""

from __future__ import annotations

from pathlib import Path

import typer

from genesets_workflows.curation import report as report_mod
from genesets_workflows.curation import validate as validate_mod

app = typer.Typer(no_args_is_help=True, add_completion=False, help="Curate GO interpretations of gene sets.")


@app.command()
def report(
    input: list[Path] = typer.Option(None, "--input", "-i", help="Interpretation YAML file(s)."),
    curation_dir: Path = typer.Option(None, "--dir", help="Directory of interpretation YAMLs."),
    output: Path = typer.Option(..., "--output", "-o", help="Output TSV path."),
) -> None:
    """Compute precision/recall/F1 for curated interpretations."""
    paths = list(input or [])
    if curation_dir:
        paths.extend(report_mod.discover(curation_dir))
    if not paths:
        raise typer.BadParameter("provide --input and/or --dir")
    report_mod.write_report(paths, output)
    typer.echo(f"wrote {output}")


@app.command()
def validate(
    path: Path = typer.Argument(..., help="Interpretation YAML file to validate."),
    schema: Path = typer.Option(Path("curation/schema/genesets_interpretation.yaml"), "--schema"),
    oak_config: Path = typer.Option(Path("curation/conf/oak_config.yaml"), "--oak-config"),
    cache_dir: Path = typer.Option(Path("curation/references_cache"), "--cache-dir"),
    ref_config: Path = typer.Option(Path("curation/conf/reference_validator_config.yaml"), "--ref-config"),
    skip_references: bool = typer.Option(False, "--skip-references", help="Skip the reference-validator step."),
) -> None:
    """Run structural, term, and reference validation on an interpretation file."""
    code = validate_mod.validate_file(
        path,
        schema=schema,
        oak_config=oak_config,
        cache_dir=cache_dir,
        skip_references=skip_references,
        ref_config=ref_config,
    )
    raise typer.Exit(code)


@app.command()
def draft(
    gene_set_id: str = typer.Argument(..., help="e.g. MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"),
    enrichment_tsv: Path = typer.Option(..., "--enrichment-tsv", help="genesets-rs enrich TSV for this set."),
    output: Path = typer.Option(..., "--output", "-o"),
    collection: str = typer.Option("C8", "--collection"),
    go_version: str = typer.Option("unspecified", "--go-version"),
) -> None:
    """Build a draft interpretation from a precomputed enrichment TSV."""
    from genesets_workflows.curation import draft as draft_mod
    from genesets_workflows.curation import model
    from genesets_workflows.curation.enrichment import parse_enrichment_tsv

    rows = parse_enrichment_tsv(Path(enrichment_tsv).read_text())
    name = gene_set_id.split(":", 1)[-1]
    interp = draft_mod.assemble_draft(
        gene_set_id=gene_set_id,
        gene_set_name=name,
        collection=collection,
        rows=rows,
        provenance=model.EnrichmentProvenance(tool="genesets-rs", go_version=go_version),
    )
    model.dump_interpretation(interp, output)
    typer.echo(f"wrote {output} ({len(interp.associations)} associations)")


def main(argv: list[str] | None = None) -> int:
    """argv bridge so the existing argparse dispatcher can invoke this Typer app."""
    from typer.main import get_command

    command = get_command(app)
    try:
        return command.main(
            args=list(argv) if argv is not None else None,
            prog_name="genesets-workflows curate",
            standalone_mode=False,
        ) or 0
    except SystemExit as exc:
        return int(exc.code or 0)
