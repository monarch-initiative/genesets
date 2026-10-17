import gzip
from pathlib import Path

import pytest
from typer.testing import CliRunner

from genesets_workflows.curation import cli as curation_cli
from genesets_workflows.curation import model
from genesets_workflows.curation.draft import assemble_ranked_draft
from genesets_workflows.curation.enrichment import GseaRow
from genesets_workflows.curation.ranking import (
    ranking_errors,
    read_ranking,
    sha256_file,
    write_ranking,
)

SCORES = {"A": 3.0, "B": 2.5, "C": 0.5, "D": -1.0, "E": -2.0}


def ranked(tmp_path: Path, *, signed: bool = True, scores=SCORES) -> model.GeneSetInterpretation:
    path = tmp_path / "rankings" / "r.rnk.gz"
    digest = write_ranking(scores, path)
    return model.GeneSetInterpretation(
        gene_set_id="LIT:R",
        gene_set_name="R",
        ranking=model.RankedSignature(
            source="PMID:1",
            score_type="t_statistic" if signed else "neg_log10_p",
            signed=signed,
            gene_id_type="HGNC symbol",
            n_ranked=len(scores),
            file="rankings/r.rnk.gz",
            sha256=digest,
        ),
    )


def association(direction: str, nes: float, genes: list[str]) -> model.TermAssociation:
    return model.TermAssociation(
        term=model.Term(id="GO:0006909", label="phagocytosis"),
        seed_source="enrichment_recovered",
        direction=direction,
        enrichment_stats=model.EnrichmentStats(
            method="gsea_ks", nes=nes, leading_edge_genes=genes
        ),
    )


def test_write_ranking_is_reproducible_and_sorted(tmp_path):
    first = write_ranking(SCORES, tmp_path / "a.rnk.gz")
    second = write_ranking(dict(reversed(SCORES.items())), tmp_path / "b.rnk.gz")
    assert first == second
    lines = gzip.open(tmp_path / "a.rnk.gz", "rt").read().splitlines()
    assert lines[:3] == ["gene\tscore", "A\t3.0", "B\t2.5"]
    assert read_ranking(tmp_path / "a.rnk.gz") == SCORES


@pytest.mark.parametrize(
    "body, message",
    [
        ("gene\tvalue\nA\t1\n", "header"),
        ("gene\tscore\nA\t1\nA\t2\n", "twice"),
        ("gene\tscore\nA\tinf\n", "non-finite"),
    ],
)
def test_read_ranking_rejects_malformed_files(tmp_path, body, message):
    path = tmp_path / "bad.rnk.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(body)
    with pytest.raises(ValueError, match=message):
        read_ranking(path)


def test_consistent_ranked_interpretation_passes(tmp_path):
    interp = ranked(tmp_path)
    interp.associations = [association("UP", 1.8, ["A", "B"]), association("DN", -1.5, ["E"])]
    assert ranking_errors(interp, tmp_path) == []


def test_ranking_gate_catches_each_inconsistency(tmp_path):
    interp = ranked(tmp_path)
    interp.ranking.sha256 = "0" * 64
    interp.ranking.n_ranked = 4
    interp.direction = "UP"
    interp.associations = [
        association("DN", 1.8, ["A", "Z"]),
        model.TermAssociation(
            term=model.Term(id="GO:0006955", label="immune response"),
            seed_source="enrichment_recovered",
            enrichment_stats=model.EnrichmentStats(method="fisher"),
        ),
    ]
    errors = "\n".join(ranking_errors(interp, tmp_path))
    for expected in (
        "sha256 mismatch",
        "n_ranked is 4",
        "set direction on associations",
        "not in the ranking: Z",
        "implies direction UP",
        "method fisher is not a ranked test",
    ):
        assert expected in errors


def test_unsigned_ranking_has_no_low_pole(tmp_path):
    interp = ranked(tmp_path, signed=False, scores={"A": 5.0, "B": 1.0, "C": 0.1})
    interp.associations = [association("DN", -1.2, ["C"])]
    errors = ranking_errors(interp, tmp_path)
    assert any("no low pole" in error for error in errors)

    negative = ranked(tmp_path, signed=False)
    assert any("negative scores" in error for error in ranking_errors(negative, tmp_path))


def test_missing_ranking_file_is_reported(tmp_path):
    interp = ranked(tmp_path)
    (tmp_path / "rankings" / "r.rnk.gz").unlink()
    assert ranking_errors(interp, tmp_path)[0].startswith("ranking file not found")


def test_ranked_draft_sets_direction_and_drops_low_pole_when_unsigned(tmp_path):
    rows = [
        GseaRow("r", "GO:1", "up", 0.6, 2.0, 1e-4, 1e-3, 20, 25, 5, 30, ["A", "B"]),
        GseaRow("r", "GO:2", "down", -0.5, -1.7, 1e-3, 1e-2, 18, 20, 4, 900, ["E"]),
    ]
    signed = ranked(tmp_path).ranking
    draft = assemble_ranked_draft(
        gene_set_id="LIT:R", gene_set_name="R", ranking=signed, rows=rows, method="gsea_ks"
    )
    assert [a.direction for a in draft.associations] == ["UP", "DN"]
    stats = draft.associations[0].enrichment_stats
    assert (stats.method, stats.p_adjust_method, stats.leading_edge_genes) == (
        "gsea_ks",
        "gsea_fdr",
        ["A", "B"],
    )
    assert ranking_errors(draft, tmp_path) == []

    signed.signed = False
    unsigned = assemble_ranked_draft(
        gene_set_id="LIT:R", gene_set_name="R", ranking=signed, rows=rows, method="wilcoxon"
    )
    assert [a.term.id for a in unsigned.associations] == ["GO:1"]
    assert unsigned.associations[0].enrichment_stats.p_adjust_method == "benjamini_hochberg"


def test_ranked_interpretation_round_trips(tmp_path):
    interp = ranked(tmp_path)
    interp.associations = [association("UP", 1.8, ["A", "B"])]
    out = tmp_path / "r.yaml"
    model.dump_interpretation(interp, out)
    reloaded = model.load_interpretation(out)
    assert reloaded.ranking == interp.ranking
    assert reloaded.associations[0].direction == "UP"
    assert reloaded.associations[0].enrichment_stats.leading_edge_genes == ["A", "B"]


def test_draft_ranked_command(tmp_path):
    ranks = tmp_path / "rankings" / "r.rnk.gz"
    write_ranking(SCORES, ranks)
    gsea = tmp_path / "gsea.tsv"
    gsea.write_text(
        "query_id\tquery_name\ttarget_id\ttarget_name\tes\tnes\tp_value\tp_adjust\tset_size\t"
        "eliminated\ttarget_size\tranked_size\tleading_edge_size\tpeak_rank\t"
        "leading_edge_genes\tleading_edge_gene_names\n"
        "r\t\tGO:0006909\tphagocytosis\t0.7\t1.9\t1e-3\t2e-2\t2\t0\t2\t5\t2\t2\tA;B\t\n"
    )
    out = tmp_path / "genesets" / "LIT_R.yaml"
    out.parent.mkdir()
    result = CliRunner().invoke(
        curation_cli.app,
        [
            "draft-ranked", "LIT:R", "--gsea-tsv", str(gsea), "--ranks", str(ranks),
            "--source", "PMID:1", "--score-type", "t_statistic", "--signed", "-o", str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    interp = model.load_interpretation(out)
    assert interp.ranking.file == "../rankings/r.rnk.gz"
    assert interp.ranking.sha256 == sha256_file(ranks)
    assert ranking_errors(interp, out.parent) == []
