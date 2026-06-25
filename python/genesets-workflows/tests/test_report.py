from pathlib import Path

import pytest

from genesets_workflows.curation import report

ROOT = Path(__file__).resolve().parents[3]
C8_DIR = ROOT / "curation" / "c8"


def test_report_rows_for_example():
    rows = report.report_rows([C8_DIR / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"])
    assert len(rows) == 1
    row = rows[0]
    assert row["gene_set_id"] == "MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"
    # 5 recovered+adjudicated; 4 of them good (translation is nonspecific)
    assert row["recovered_adjudicated"] == 5
    assert row["precision"] == 0.8
    # core_total = antigen presentation (core_process, recovered) +
    # MHC class II protein complex (core_component, recovered) +
    # adaptive immune response (core_process, curator_added/missed)
    assert row["core_total"] == 3
    # 2 of the 3 core terms were recovered by enrichment
    assert row["recall"] == pytest.approx(2 / 3)


def test_write_report_tsv(tmp_path):
    out = tmp_path / "report.tsv"
    report.write_report([C8_DIR / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"], out)
    header = out.read_text().splitlines()[0]
    assert header.split("\t")[0] == "gene_set_id"
    assert "precision" in header
