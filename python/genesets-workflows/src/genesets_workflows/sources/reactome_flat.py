from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

from genesets_workflows.runtime import utc_now
from genesets_workflows.yaml_io import write_yaml_value


DEFAULT_URL = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
DEFAULT_OUT_DIR = Path("evals/reactome_flat/generated/current")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="genesets-workflows prepare-reactome-flat",
        description="Prepare official Reactome pathway gene sets as a flat target library."
    )
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR, type=Path)
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--species-prefix", default="R-HSA-")
    parser.add_argument("--min-genes", type=int, default=2)
    parser.add_argument(
        "--max-genes",
        type=int,
        default=0,
        help="Maximum genes per pathway; 0 means no maximum.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Download the Reactome zip even if it already exists locally.",
    )
    return parser.parse_args(argv)


def download(url: str, path: Path, force: bool) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    headers: dict[str, str] = {}
    if force or not path.exists():
        request = urllib.request.Request(url, headers={"User-Agent": "genesets-workflows/0.1"})
        with urllib.request.urlopen(request, timeout=120) as response:
            headers = {key.lower(): value for key, value in response.headers.items()}
            path.write_bytes(response.read())
    return {
        **file_metadata(path, url),
        "last_modified": headers.get("last-modified"),
        "etag": headers.get("etag"),
    }


def file_metadata(path: Path, url: str) -> dict[str, Any]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {
        "url": url,
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def extract_gmt(zip_path: Path, out_dir: Path) -> Path:
    with zipfile.ZipFile(zip_path) as archive:
        names = [name for name in archive.namelist() if name.endswith(".gmt")]
        if len(names) != 1:
            raise SystemExit(f"expected exactly one GMT in {zip_path}, found {names}")
        raw_path = out_dir / names[0]
        with archive.open(names[0]) as source, raw_path.open("wb") as target:
            target.write(source.read())
    return raw_path


def normalize_gmt(
    raw_path: Path,
    targets_path: Path,
    terms_path: Path,
    background_path: Path,
    species_prefix: str,
    min_genes: int,
    max_genes: int,
) -> dict[str, Any]:
    total_rows = 0
    emitted_rows = 0
    skipped_species = 0
    skipped_small = 0
    skipped_large = 0
    background: set[str] = set()
    pathway_sizes: list[int] = []

    with raw_path.open(newline="") as raw, targets_path.open("w", newline="") as targets, terms_path.open(
        "w", newline=""
    ) as terms:
        target_writer = csv.writer(targets, delimiter="\t", lineterminator="\n")
        term_writer = csv.writer(terms, delimiter="\t", lineterminator="\n")
        term_writer.writerow(["term_id", "name"])

        for row in csv.reader(raw, delimiter="\t"):
            if len(row) < 3:
                continue
            total_rows += 1
            name = row[0].strip()
            pathway_id = row[1].strip()
            if species_prefix and not pathway_id.startswith(species_prefix):
                skipped_species += 1
                continue
            genes = sorted({gene.strip() for gene in row[2:] if gene.strip()})
            if len(genes) < min_genes:
                skipped_small += 1
                continue
            if max_genes and len(genes) > max_genes:
                skipped_large += 1
                continue

            target_writer.writerow([pathway_id, name, *genes])
            term_writer.writerow([pathway_id, name])
            background.update(genes)
            pathway_sizes.append(len(genes))
            emitted_rows += 1

    with background_path.open("w", newline="") as handle:
        handle.write("gene_id\n")
        for gene in sorted(background):
            handle.write(f"{gene}\n")

    pathway_sizes.sort()
    return {
        "raw_rows": total_rows,
        "emitted_pathways": emitted_rows,
        "skipped_species": skipped_species,
        "skipped_small": skipped_small,
        "skipped_large": skipped_large,
        "background_gene_count": len(background),
        "min_pathway_genes": pathway_sizes[0] if pathway_sizes else 0,
        "median_pathway_genes": pathway_sizes[len(pathway_sizes) // 2] if pathway_sizes else 0,
        "max_pathway_genes": pathway_sizes[-1] if pathway_sizes else 0,
    }


def prepare(args: argparse.Namespace) -> dict[str, Any]:
    if args.min_genes < 1:
        raise SystemExit("--min-genes must be positive")
    if args.max_genes < 0:
        raise SystemExit("--max-genes must be non-negative")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    downloads = args.out_dir / "downloads"
    zip_path = downloads / "ReactomePathways.gmt.zip"
    download_metadata = download(args.url, zip_path, args.force_download)
    raw_gmt = extract_gmt(zip_path, downloads)

    targets_path = args.out_dir / "targets.gmt"
    terms_path = args.out_dir / "terms.tsv"
    background_path = args.out_dir / "background.txt"
    stats = normalize_gmt(
        raw_gmt,
        targets_path,
        terms_path,
        background_path,
        args.species_prefix,
        args.min_genes,
        args.max_genes,
    )

    metadata = {
        "generated_at_utc": utc_now(),
        "source": "Reactome official pathway gene set GMT",
        "download": download_metadata,
        "raw_gmt": str(raw_gmt),
        "normalization": {
            "input_format": "Reactome GMT: pathway name, pathway stable id, gene symbols...",
            "output_format": "genesets-rs GMT: pathway stable id, pathway name, gene symbols...",
            "species_prefix": args.species_prefix,
            "min_genes": args.min_genes,
            "max_genes": args.max_genes,
        },
        "stats": stats,
        "outputs": {
            "target_sets": str(targets_path),
            "terms": str(terms_path),
            "background": str(background_path),
        },
    }
    with (args.out_dir / "metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (args.out_dir / "metadata.yaml").open("w") as handle:
        write_yaml_value(handle, metadata)
    return metadata


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    metadata = prepare(args)
    stats = metadata["stats"]
    print(
        f"Prepared {stats['emitted_pathways']} Reactome pathways "
        f"over {stats['background_gene_count']} genes under {args.out_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
