#!/usr/bin/env bash
# Smoke-test chemical-to-GO enrichment: does each bridge recover the pathway a
# small, unambiguous metabolite set comes from?
#
# The background is a proxy (annotated chemicals with an HMDB cross-reference),
# not a platform's measured metabolites, so p-values are illustrative only.
set -euo pipefail

prep=${1:-evals/chem2go/generated/current}
here=$(cd "$(dirname "$0")" && pwd)
out="$here/generated/smoke"
gw=(uv run --project python/genesets-workflows genesets-workflows)
bin=${GENESETS_RS:-target/release/genesets-rs}
mkdir -p "$out"

awk -F'\t' 'NR > 1 { print $1 "\t" $2 }' "$prep/chemicals.tsv" > "$out/names.tsv"
uv run --project python/genesets-workflows python - "$prep" "$out/background.txt" <<'PY'
import csv, sys
from pathlib import Path
prep, out = Path(sys.argv[1]), Path(sys.argv[2])
canon = {r["chebi_id"]: r["canonical_id"] for r in csv.DictReader((prep / "normalization.tsv").open(), delimiter="\t")}
annotated = {r["gene_id"] for r in csv.DictReader((prep / "annotations.all.tsv").open(), delimiter="\t")}
hmdb = {canon[r["chebi_id"]] for r in csv.DictReader((prep / "xrefs.tsv").open(), delimiter="\t")
        if r["xref"].startswith("hmdb:") and r["chebi_id"] in canon}
out.write_text("".join(f"{c}\n" for c in sorted(hmdb & annotated)))
PY

variants=(logical_definition rhea all)
[[ -f "$prep/annotations.enzyme.tsv" ]] && variants=(logical_definition rhea enzyme all)

for query in tca urea; do
  "${gw[@]}" normalize-chemicals "$here/smoke/$query.txt" --prep-dir "$prep" \
    --extra-xrefs "$here/smoke/extra_xrefs.tsv" \
    --out "$out/$query.canonical.txt" --report "$out/$query.report.tsv"
  for variant in "${variants[@]}"; do
    "$bin" enrich \
      --annotations "$prep/annotations.$variant.tsv" \
      --terms "$prep/terms.tsv" --closure "$prep/closure.tsv" \
      --gene-names "$out/names.tsv" \
      --sample "$out/$query.canonical.txt" --background "$out/background.txt" \
      --min-overlap 2 --max-p-adjust 0.05 --overlap-genes \
      -o "$out/$query.$variant.tsv"
    echo "== $query / $variant"
    awk -F'\t' 'NR > 1 && NR <= 6 { printf "  %s %s (%s/%s, padj %s)\n", $3, $4, $5, $7, $10 }' "$out/$query.$variant.tsv"
  done
done
