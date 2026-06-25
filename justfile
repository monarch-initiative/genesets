set shell := ["bash", "-euo", "pipefail", "-c"]

iba_bundle := "notebooks/generated/go_iba_impact_expression5000_diverse"
go5y_bundle := "notebooks/generated/go_impact_5y_expression5000_diverse"
contributes_bundle := "notebooks/generated/go_contributes_to_impact_expression4313_nongo"
default_bundles := iba_bundle + " " + go5y_bundle + " " + contributes_bundle
browser_port := "8765"

# List available recipes.
default:
    @just --list

# Stop a local explorer server on the default port.
browser-stop:
    @pids="$(lsof -tiTCP:{{browser_port}} -sTCP:LISTEN || true)"; \
      if [[ -n "$pids" ]]; then \
        echo "Stopping explorer server on port {{browser_port}}: $pids"; \
        kill $pids; \
      else \
        echo "No explorer server running on port {{browser_port}}"; \
      fi

# Open the web explorer for the default report bundles.
browser: browser-stop
    uv run --project python/genesets-workflows --extra explorer \
      genesets-workflows explore {{default_bundles}} --open --port {{browser_port}}

# Open the web explorer for the current GOA all-vs-IBA comparison.
browser-iba: browser-stop
    uv run --project python/genesets-workflows --extra explorer \
      genesets-workflows explore {{iba_bundle}} --open --port {{browser_port}}

# Open the web explorer for the 2021-vs-2026 GO/GOA comparison.
browser-go5y: browser-stop
    uv run --project python/genesets-workflows --extra explorer \
      genesets-workflows explore {{go5y_bundle}} --open --port {{browser_port}}

# Open the web explorer for the current GOA all-vs-no-contributes_to comparison.
browser-contributes: browser-stop
    uv run --project python/genesets-workflows --extra explorer \
      genesets-workflows explore {{contributes_bundle}} --open --port {{browser_port}}

# Open the web explorer for a specific report bundle.
browser-bundle bundle: browser-stop
    uv run --project python/genesets-workflows --extra explorer \
      genesets-workflows explore {{bundle}} --open --port {{browser_port}}

# --- Curation ---------------------------------------------------------------

curation_schema := "curation/schema/genesets_interpretation.yaml"
curation_oak := "curation/conf/oak_config.yaml"
gw := "uv run --project python/genesets-workflows --extra curation"

# Run the curation test suite (pytest + doctests).
curate-test:
    {{gw}} --with pytest pytest python/genesets-workflows/tests -v
    {{gw}} --with pytest python -m pytest --doctest-modules \
      python/genesets-workflows/src/genesets_workflows/curation -v

# Validate every interpretation YAML under curation/genesets.
curate-validate:
    for f in curation/genesets/*.yaml; do \
      echo "== $f =="; \
      {{gw}} genesets-workflows curate validate "$f" || exit 1; \
    done

# Validate the schema's own ontology-term meanings.
curate-validate-schema:
    {{gw}} linkml-term-validator validate-schema {{curation_schema}} -c {{curation_oak}}

# Build the precision/recall report over all curated C8 interpretations.
curate-report out="curation/report.tsv":
    {{gw}} genesets-workflows curate report --dir curation/genesets -o {{out}}
