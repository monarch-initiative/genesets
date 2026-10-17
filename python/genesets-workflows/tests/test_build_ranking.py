import gzip
import hashlib
from pathlib import Path

import pytest
import yaml

from genesets_workflows.curation.build_ranking import HgncMapper, build
from genesets_workflows.curation.ranking import read_ranking

HGNC = (
    "hgnc_id\tsymbol\tstatus\tprev_symbol\talias_symbol\tentrez_id\tensembl_gene_id\n"
    "HGNC:1\tELOB\tApproved\tTCEB2\tSIII\t6923\tENSG00000103363\n"
    "HGNC:2\tGRIN2A\tApproved\t\tNR2A\t2903\tENSG00000183454\n"
    "HGNC:3\tAAA1\tApproved\t\tSHARED\t1\tENSG00000000001\n"
    "HGNC:4\tAAA2\tApproved\t\tSHARED\t2\tENSG00000000002\n"
    "HGNC:5\tOLD1\tEntry Withdrawn\t\t\t\t\n"
)


def test_mapper_resolves_ids_and_refuses_ambiguous_aliases():
    hgnc = HgncMapper(HGNC)
    assert hgnc.map("ENSG00000183454.12", "ensembl") == "GRIN2A"
    assert hgnc.map("6923.0", "entrez") == "ELOB"
    assert hgnc.map("TCEB2", "symbol") == "ELOB"
    assert hgnc.map("NR2A", "symbol") == "GRIN2A"
    assert hgnc.map("SHARED", "symbol") is None
    assert hgnc.map("OLD1", "symbol") is None


def write(path: Path, data: bytes) -> tuple[str, str]:
    path.write_bytes(data)
    return path.as_uri(), hashlib.sha256(data).hexdigest()


def recipe(tmp_path: Path, table: bytes, **overrides) -> Path:
    hgnc_url, hgnc_sha = write(tmp_path / "hgnc.txt", HGNC.encode())
    source_url, source_sha = write(tmp_path / "source.tsv.gz", gzip.compress(table))
    body = {
        "name": "TEST",
        "source": {"url": source_url, "sha256": source_sha},
        "table": {"format": "tsv"},
        "gene_column": "gene",
        "gene_id_type": "symbol",
        "score_column": "p",
        "transform": "neg_log10",
        "hgnc": {"url": hgnc_url, "sha256": hgnc_sha},
    }
    body.update(overrides)
    path = tmp_path / "rankings" / "TEST.recipe.yaml"
    path.parent.mkdir(exist_ok=True)
    path.write_text(yaml.safe_dump(body))
    return path


def test_build_maps_transforms_and_collapses(tmp_path):
    table = b"gene\tp\nTCEB2\t0.01\nELOB\t0.1\nGRIN2A\t1e-6\nUNKNOWN\t0.5\nAAA1\tNA\n"
    path = recipe(tmp_path, table)
    report = build(path, tmp_path / "cache")
    assert (report.rows, report.missing, report.unmapped, report.collapsed, report.ranked) == (
        5, 1, 1, 1, 2,
    )
    scores = read_ranking(path.parent / "TEST.rnk.gz")
    assert scores == pytest.approx({"GRIN2A": 6.0, "ELOB": 2.0})
    # Rebuilding gives the same file.
    assert build(path, tmp_path / "cache").sha256 == report.sha256


def test_build_refuses_a_changed_download(tmp_path):
    path = recipe(tmp_path, b"gene\tp\nELOB\t0.1\n")
    body = yaml.safe_load(path.read_text())
    body["source"]["sha256"] = "0" * 64
    path.write_text(yaml.safe_dump(body))
    with pytest.raises(ValueError, match="sha256"):
        build(path, tmp_path / "cache")
