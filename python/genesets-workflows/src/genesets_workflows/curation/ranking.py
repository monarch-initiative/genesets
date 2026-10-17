"""Ranked-signature files and the checks that tie an interpretation to one.

A ranked interpretation points at a gzip-compressed two-column TSV (`gene`,
`score`) through its ``ranking`` block. The gate here checks that the file is
the one the curator recorded (checksum, gene count, well-formed scores) and
that every association's ranked statistics agree with it.
"""

from __future__ import annotations

import gzip
import hashlib
import math
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from genesets_workflows.curation.model import GeneSetInterpretation

RANKED_METHODS = frozenset({"gsea_ks", "wilcoxon"})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_ranking(path: Path) -> dict[str, float]:
    """Read a gzipped `gene<TAB>score` TSV with that header into a dict.

    Raises ValueError on a wrong header, a malformed or non-finite score, or a
    gene listed twice.
    """
    scores: dict[str, float] = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["gene", "score"]:
            raise ValueError(f"{path}: header must be 'gene<TAB>score', got {header}")
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[0]:
                raise ValueError(f"{path}:{line_number}: expected two columns")
            gene, raw = fields
            score = float(raw)
            if not math.isfinite(score):
                raise ValueError(f"{path}:{line_number}: non-finite score for {gene}")
            if gene in scores:
                raise ValueError(f"{path}:{line_number}: {gene} listed twice")
            scores[gene] = score
    return scores


def write_ranking(scores: dict[str, float], path: Path) -> str:
    """Write scores highest first as a reproducible gzipped TSV; return its SHA-256.

    The gzip header carries no timestamp or file name, so the same scores give
    the same checksum on every run.
    """
    ordered = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    text = "gene\tscore\n" + "".join(f"{gene}\t{score!r}\n" for gene, score in ordered)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, mtime=0
    ) as handle:
        handle.write(text.encode())
    return sha256_file(path)


def ranking_errors(interp: GeneSetInterpretation, base_dir: Path) -> list[str]:
    """Problems linking ``interp`` to its ranking file, or misusing ranked fields.

    ``base_dir`` is the interpretation file's directory; ``ranking.file`` is
    resolved against it.

    >>> from genesets_workflows.curation import model
    >>> plain = model.GeneSetInterpretation(
    ...     gene_set_id="LIT:X", gene_set_name="X",
    ...     associations=[model.TermAssociation(
    ...         term=model.Term(id="GO:0006909", label="phagocytosis"),
    ...         seed_source="enrichment_recovered", direction="UP",
    ...         enrichment_stats=model.EnrichmentStats(method="gsea_ks"))])
    >>> for error in ranking_errors(plain, Path(".")):
    ...     print(error)
    GO:0006909: method gsea_ks needs a ranking block
    GO:0006909: association direction is only for ranked signatures
    """
    errors: list[str] = []
    ranking = interp.ranking
    if ranking is None:
        for association in interp.associations:
            stats = association.enrichment_stats
            term = association.term.id
            if stats is not None and stats.method in RANKED_METHODS:
                errors.append(f"{term}: method {stats.method} needs a ranking block")
            if association.direction is not None:
                errors.append(f"{term}: association direction is only for ranked signatures")
        return errors

    if interp.direction not in (None, "NA"):
        errors.append(
            "a ranked signature carries both poles: set direction on associations, not the set"
        )
    path = base_dir / ranking.file
    if not path.exists():
        recipe = path.with_name(path.name.removesuffix(".rnk.gz") + ".recipe.yaml")
        hint = (
            f"; build it with `genesets-workflows curate build-ranking {recipe}`"
            if recipe.exists()
            else ""
        )
        return errors + [f"ranking file not found: {path}{hint}"]
    actual = sha256_file(path)
    if actual != ranking.sha256:
        errors.append(f"ranking sha256 mismatch: recorded {ranking.sha256}, file has {actual}")
    try:
        scores = read_ranking(path)
    except ValueError as error:
        return errors + [str(error)]
    if len(scores) != ranking.n_ranked:
        errors.append(f"n_ranked is {ranking.n_ranked} but the file ranks {len(scores)} genes")

    for association in interp.associations:
        term = association.term.id
        stats = association.enrichment_stats
        if association.direction == "DN" and not ranking.signed:
            errors.append(f"{term}: direction DN on an unsigned ranking, which has no low pole")
        if stats is None:
            continue
        if stats.method is not None and stats.method not in RANKED_METHODS:
            errors.append(f"{term}: method {stats.method} is not a ranked test")
        missing = sorted(set(stats.leading_edge_genes) - scores.keys())
        if missing:
            errors.append(f"{term}: leading-edge genes not in the ranking: {', '.join(missing)}")
        if stats.nes is not None and stats.nes != 0:
            expected = "UP" if stats.nes > 0 else "DN"
            if association.direction != expected:
                errors.append(
                    f"{term}: NES {stats.nes} implies direction {expected}, "
                    f"recorded {association.direction}"
                )
    return errors
