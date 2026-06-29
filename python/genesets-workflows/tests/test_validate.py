from pathlib import Path

from genesets_workflows.curation import model, validate


def test_term_ids_by_prefix_groups_and_skips_taxon():
    interp = model.GeneSetInterpretation(
        gene_set_id="LIT:X",
        gene_set_name="X",
        taxon=model.Term(id="NCBITaxon:9606", label="Homo sapiens"),
        contexts=[
            model.BiologicalContext(
                term=model.Term(id="MONDO:0004976", label="amyotrophic lateral sclerosis"),
                context_type="disease",
            )
        ],
        associations=[
            model.TermAssociation(
                term=model.Term(id="GO:0006909", label="phagocytosis"),
                seed_source="enrichment_recovered",
            ),
            model.TermAssociation(
                term=model.Term(id="GO:0005764", label="lysosome"),
                seed_source="enrichment_recovered",
            ),
        ],
    )
    by = validate.term_ids_by_prefix(interp)
    assert by == {"MONDO": {"MONDO:0004976"}, "GO": {"GO:0006909", "GO:0005764"}}
    assert "NCBITaxon" not in by


def test_build_commands():
    cmds = validate.build_commands(
        Path("f.yaml"),
        schema=Path("s.yaml"),
        oak_config=Path("oak.yaml"),
        cache_dir=Path("rc"),
        skip_references=False,
    )
    names = [c[0] for c in cmds]
    assert names == ["linkml-validate", "linkml-term-validator", "linkml-reference-validator"]
    assert "--target-class" in cmds[0] and "GeneSetInterpretation" in cmds[0]
    assert "validate-data" in cmds[1] and "--labels" in cmds[1]
    assert cmds[2][:3] == ["linkml-reference-validator", "validate", "data"]


def test_ref_config_is_wired_into_reference_command():
    cmds = validate.build_commands(
        Path("f.yaml"),
        schema=Path("s.yaml"),
        oak_config=Path("oak.yaml"),
        cache_dir=Path("rc"),
        skip_references=False,
        ref_config=Path("ref.yaml"),
    )
    assert "--config" in cmds[2]
    assert cmds[2][cmds[2].index("--config") + 1] == "ref.yaml"


def test_skip_references_drops_third_command():
    cmds = validate.build_commands(
        Path("f.yaml"),
        schema=Path("s.yaml"),
        oak_config=Path("oak.yaml"),
        cache_dir=Path("rc"),
        skip_references=True,
    )
    assert len(cmds) == 2
