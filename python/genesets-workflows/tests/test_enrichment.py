from pathlib import Path

from genesets_workflows.curation import enrichment

DATA = Path(__file__).resolve().parent / "data" / "enrich_sample.tsv"


def test_parse_enrichment_tsv():
    rows = enrichment.parse_enrichment_tsv(DATA.read_text())
    assert len(rows) == 2
    first = rows[0]
    assert first.target_id == "GO:0019882"
    assert first.target_name == "antigen processing and presentation"
    assert first.overlap == 12
    assert first.p_adjust == 3.4e-07
    assert first.overlap_genes == ["HLA-DRA", "CD74"]


def test_parse_handles_bonferroni_column():
    text = "query_id\tquery_name\ttarget_id\ttarget_name\toverlap\tquery_size\ttarget_size\tbackground_size\tp_value\tp_adjust_bonferroni\nq\tq\tGO:1\tt\t1\t2\t3\t4\t0.1\t0.5\n"
    rows = enrichment.parse_enrichment_tsv(text)
    assert rows[0].p_adjust == 0.5
    assert rows[0].overlap_genes == []


import os
import stat


def test_run_enrichment_invokes_binary(tmp_path):
    stub = tmp_path / "genesets-rs"
    stub.write_text("#!/usr/bin/env bash\ncat " + str(DATA) + "\n")
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    for name in ("sample", "annotations", "terms", "closure", "background"):
        (tmp_path / name).write_text("x\n")
    out = enrichment.run_enrichment(
        sample=tmp_path / "sample",
        annotations=tmp_path / "annotations",
        terms=tmp_path / "terms",
        closure=tmp_path / "closure",
        background=tmp_path / "background",
        binary=str(stub),
    )
    rows = enrichment.parse_enrichment_tsv(out)
    assert rows[0].target_id == "GO:0019882"
