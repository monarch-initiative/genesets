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


def test_fallback_symbol_column_accepts_only_approved_symbols(tmp_path):
    table = (
        b"gene\tsym\tp\n"
        b"ENSG00000183454\tGRIN2A\t0.01\n"
        b"ENSG00000999999\tELOB\t0.1\n"
        b"ENSG00000999998\tTCEB2\t0.5\n"
    )
    path = recipe(
        tmp_path, table, gene_column="gene", gene_id_type="ensembl",
        fallback_symbol_column="sym",
    )
    report = build(path, tmp_path / "cache")
    assert (report.fallback, report.unmapped, report.ranked) == (1, 1, 2)


def test_drop_discards_every_row_of_a_duplicated_symbol(tmp_path):
    table = b"gene\tp\nTCEB2\t0.01\nELOB\t0.1\nGRIN2A\t1e-6\n"
    report = build(recipe(tmp_path, table, duplicates="drop"), tmp_path / "cache")
    assert (report.collapsed, report.ranked) == (1, 1)


ROOT = Path(__file__).resolve().parents[3]


def test_committed_rankings_match_store_flags():
    """Every ranking has a recipe; exactly the `store: true` ones are committed."""
    rankings = ROOT / "curation" / "rankings"
    recipes = {p.name.removesuffix(".recipe.yaml"): yaml.safe_load(p.read_text())
               for p in rankings.glob("*.recipe.yaml")}
    whitelisted = {
        line.strip().removeprefix("!curation/rankings/").removesuffix(".rnk.gz")
        for line in (ROOT / ".gitignore").read_text().splitlines()
        if line.strip().startswith("!curation/rankings/")
    }
    assert whitelisted == {name for name, body in recipes.items() if body.get("store")}
    for name in whitelisted:
        assert (rankings / f"{name}.rnk.gz").exists(), name


def test_every_ranked_interpretation_has_a_recipe():
    from genesets_workflows.curation import model

    for path in (ROOT / "curation" / "genesets").glob("*.yaml"):
        interp = model.load_interpretation(path)
        if interp.ranking is None:
            continue
        built = (path.parent / interp.ranking.file).resolve()
        assert built.with_name(built.name.removesuffix(".rnk.gz") + ".recipe.yaml").exists(), path


def test_prefer_exact_lets_the_approved_symbol_row_win(tmp_path):
    # TCEB2 is ELOB's previous symbol: the exact ELOB row wins though weaker.
    table = b"gene\tp\nTCEB2\t1e-9\nELOB\t0.1\nGRIN2A\t1e-3\nGRIN2A\t1e-6\n"
    path = recipe(tmp_path, table, duplicates="prefer_exact")
    build(path, tmp_path / "cache")
    scores = read_ranking(path.parent / "TEST.rnk.gz")
    assert scores == pytest.approx({"ELOB": 1.0, "GRIN2A": 6.0})
