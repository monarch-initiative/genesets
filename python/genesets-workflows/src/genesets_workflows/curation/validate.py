"""Build and run the three validators over an interpretation file."""

from __future__ import annotations

import subprocess
from pathlib import Path

TARGET_CLASS = "GeneSetInterpretation"


def build_commands(
    path: Path,
    *,
    schema: Path,
    oak_config: Path,
    cache_dir: Path,
    skip_references: bool,
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
        commands.append(
            [
                "linkml-reference-validator",
                "validate", "data", str(path),
                "--schema", str(schema),
                "--target-class", TARGET_CLASS,
                "--cache-dir", str(cache_dir),
            ]
        )
    return commands


def validate_file(
    path: Path,
    *,
    schema: Path,
    oak_config: Path,
    cache_dir: Path,
    skip_references: bool = False,
) -> int:
    """Run each validator; return 0 if all succeed, else the first non-zero code."""
    for command in build_commands(
        path, schema=schema, oak_config=oak_config, cache_dir=cache_dir, skip_references=skip_references
    ):
        print(f"$ {' '.join(command)}")
        completed = subprocess.run(command)
        if completed.returncode != 0:
            return completed.returncode
    return 0
