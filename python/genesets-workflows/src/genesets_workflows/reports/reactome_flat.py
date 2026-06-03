from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml

from genesets_workflows.runtime import duckdb_json, run_command, utc_now
from genesets_workflows.sources import reactome_flat as reactome_source
from genesets_workflows.yaml_io import write_yaml_value


DEFAULT_CONFIG = Path("evals/reactome_flat/config.yaml")
DEFAULT_METADATA = Path("evals/reactome_flat/generated/current/run_metadata.yaml")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="genesets-workflows reactome-flat",
        description="Run expression queries against official Reactome as a flat library."
    )
    parser.add_argument("--config", default=DEFAULT_CONFIG, type=Path)
    parser.add_argument("--metadata-output", default=DEFAULT_METADATA, type=Path)
    parser.add_argument("--binary", default="genesets-rs")
    parser.add_argument("--skip-prepare", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--reactome-out-dir", default=reactome_source.DEFAULT_OUT_DIR, type=Path)
    parser.add_argument("--force-download", action="store_true")
    return parser.parse_args(argv)


def resolve_config_path(config: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else config.parent / path


def load_output_path(config: Path) -> Path:
    with config.open() as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict) or "output" not in data:
        raise SystemExit(f"{config} must contain an output path")
    return resolve_config_path(config, data["output"])


def summarize_output(path: Path) -> dict[str, Any] | None:
    rows = duckdb_json(
        f"""
        SELECT
          count(*) AS rows,
          count(DISTINCT query_id) AS queries_with_hits,
          count(DISTINCT target_id) AS pathways_with_hits
        FROM read_parquet('{path}')
        """
    )
    return rows[0] if rows else None


def run_report(args: argparse.Namespace) -> dict[str, Any]:
    output = load_output_path(args.config)
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)

    prepare_result: dict[str, Any] | None = None
    if not args.skip_prepare:
        source_args = reactome_source.parse_args(
            [
                "--out-dir",
                str(args.reactome_out_dir),
                *(["--force-download"] if args.force_download else []),
            ]
        )
        prepare_result = reactome_source.prepare(source_args)

    if args.skip_existing and output.exists():
        run_result = {
            "command": [args.binary, "run", str(args.config)],
            "runtime_seconds": 0.0,
            "stdout": "",
            "stderr": "",
            "status": "reused_existing",
        }
    else:
        run_result = run_command([args.binary, "run", str(args.config)])
        run_result["status"] = "generated"

    metadata = {
        "generated_at_utc": utc_now(),
        "eval": "expression500_vs_reactome_flat",
        "description": "Expression-derived MSigDB/GSE query sets enriched against official Reactome pathway GMT as a flat target library.",
        "config": str(args.config),
        "prepare": prepare_result,
        "run": run_result,
        "output": str(output),
        "summary": summarize_output(output),
    }
    with args.metadata_output.open("w") as handle:
        write_yaml_value(handle, metadata)
    return metadata


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    metadata = run_report(args)
    print(f"{metadata['run']['status']}: {metadata['output']}")
    if metadata["summary"]:
        print(
            "Rows: "
            f"{metadata['summary']['rows']}, "
            f"queries_with_hits={metadata['summary']['queries_with_hits']}, "
            f"pathways_with_hits={metadata['summary']['pathways_with_hits']}"
        )
    print(f"metadata: {args.metadata_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
