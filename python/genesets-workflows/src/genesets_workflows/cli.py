from __future__ import annotations

import argparse
import sys
from typing import Callable

from genesets_workflows import __version__
from genesets_workflows.reports import go_impact, reactome_flat
from genesets_workflows.runtime import version_output, which
from genesets_workflows.sources import mygeneset
from genesets_workflows.sources import mygeneset_stratified
from genesets_workflows.sources import reactome_flat as reactome_source
from genesets_workflows.curation import cli as curation_cli


Command = Callable[[list[str] | None], int]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="genesets-workflows",
        description="Convenience workflows and reports around the genesets-rs Rust CLI.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("doctor", help="Check local workflow dependencies.")
    subparsers.add_parser(
        "go-impact",
        help="Compare two prepared GO snapshots over a query collection.",
        add_help=False,
    )
    subparsers.add_parser(
        "reactome-flat",
        help="Run expression queries against official Reactome as a flat library.",
        add_help=False,
    )
    subparsers.add_parser(
        "prepare-reactome-flat",
        help="Prepare official Reactome GMT as a flat genesets-rs target library.",
        add_help=False,
    )
    subparsers.add_parser(
        "fetch-mygeneset",
        help="Fetch MyGeneSet query results into GMT.",
        add_help=False,
    )
    subparsers.add_parser(
        "fetch-mygeneset-stratified",
        help="Fetch a stratified MyGeneSet query collection into GMT.",
        add_help=False,
    )
    subparsers.add_parser(
        "explore",
        help="Launch a local browser over workflow report bundles.",
        add_help=False,
    )
    subparsers.add_parser(
        "curate",
        help="Curate and validate GO interpretations of gene sets.",
        add_help=False,
    )
    parsed, remainder = parser.parse_known_args(argv)
    parsed.remainder = remainder
    return parsed


def doctor(_argv: list[str] | None = None) -> int:
    checks = [
        ("genesets-rs", which("genesets-rs"), version_output(["genesets-rs", "--version"])),
        ("duckdb CLI", which("duckdb"), version_output(["duckdb", "--version"])),
        ("uv", which("uv"), version_output(["uv", "--version"])),
    ]
    try:
        import yaml  # noqa: F401

        yaml_status = "ok"
    except ImportError:
        yaml_status = "missing"
    try:
        import duckdb  # noqa: F401

        duckdb_python_status = "ok"
    except ImportError:
        duckdb_python_status = "missing; CLI fallback is supported"
    try:
        import fastapi  # noqa: F401
        import uvicorn  # noqa: F401

        explorer_status = "ok"
    except ImportError:
        explorer_status = "missing; install the explorer extra for web browsing"

    print(f"genesets-workflows {__version__}")
    for name, path, version in checks:
        if path:
            suffix = f" ({version})" if version else ""
            print(f"ok      {name}: {path}{suffix}")
        else:
            print(f"missing {name}")
    print(f"ok      PyYAML: {yaml_status}" if yaml_status == "ok" else f"missing PyYAML: {yaml_status}")
    print(f"ok      Python duckdb: {duckdb_python_status}" if duckdb_python_status == "ok" else f"warn    Python duckdb: {duckdb_python_status}")
    print(f"ok      Explorer web deps: {explorer_status}" if explorer_status == "ok" else f"warn    Explorer web deps: {explorer_status}")
    return 0


def dispatch(command: str) -> Command:
    def explore(argv: list[str] | None = None) -> int:
        from genesets_workflows.explorer import server

        return server.main(argv)

    commands: dict[str, Command] = {
        "doctor": doctor,
        "go-impact": go_impact.main,
        "reactome-flat": reactome_flat.main,
        "prepare-reactome-flat": reactome_source.main,
        "fetch-mygeneset": mygeneset.main,
        "fetch-mygeneset-stratified": mygeneset_stratified.main,
        "explore": explore,
        "curate": curation_cli.main,
    }
    return commands[command]


def main(argv: list[str] | None = None) -> int:
    parsed = parse_args(argv)
    return dispatch(parsed.command)(parsed.remainder)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
