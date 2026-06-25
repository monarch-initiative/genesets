from pathlib import Path

from typer.testing import CliRunner

from genesets_workflows.curation import cli as curation_cli
from genesets_workflows.cli import dispatch

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "curation" / "c8" / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"

runner = CliRunner()


def test_report_command_writes_tsv(tmp_path):
    out = tmp_path / "r.tsv"
    result = runner.invoke(curation_cli.app, ["report", "--input", str(EXAMPLE), "--output", str(out)])
    assert result.exit_code == 0, result.output
    assert out.exists()
    assert "precision" in out.read_text().splitlines()[0]


def test_help_lists_subcommands():
    result = runner.invoke(curation_cli.app, ["--help"])
    assert result.exit_code == 0
    for name in ("draft", "validate", "report"):
        assert name in result.output


def test_dispatcher_routes_curate(tmp_path):
    # curate is dispatched lazily (mirrors the `explore` extra pattern), so the
    # routed handler runs the curate CLI end-to-end rather than being identical
    # to curation_cli.main.
    out = tmp_path / "r.tsv"
    rc = dispatch("curate")(["report", "-i", str(EXAMPLE), "-o", str(out)])
    assert rc == 0
    assert out.exists()
