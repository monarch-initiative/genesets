from pathlib import Path

from genesets_workflows.curation import validate


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


def test_skip_references_drops_third_command():
    cmds = validate.build_commands(
        Path("f.yaml"),
        schema=Path("s.yaml"),
        oak_config=Path("oak.yaml"),
        cache_dir=Path("rc"),
        skip_references=True,
    )
    assert len(cmds) == 2
