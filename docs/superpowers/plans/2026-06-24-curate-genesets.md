# Curated GO Interpretations for Non-GO Gene Sets — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a LinkML-schema'd, validator-backed curation system in `genesets-rs` for recording the "correct/reasonable" GO-term interpretation of non-GO gene sets (MSigDB C8 cell-type signatures first), seeded from real enrichment output and adjudicated by curators, doubling as a precision/recall gold standard.

**Architecture:** A new top-level `curation/` tree holds a LinkML schema, OAK/validator config, and one YAML interpretation per gene set. A `curation` subpackage in `python/genesets-workflows` provides a Typer `curate` command group (`draft`, `validate`, `report`) wired into the existing argparse dispatcher. Drafting shells out to the `genesets-rs` enrichment kernel; validation shells out to `linkml-validate` + `linkml-term-validator` + `linkml-reference-validator` (dismech's anti-hallucination pattern); reporting computes precision/recall/F1 against the gold standard.

**Tech Stack:** Python 3.10+, Typer, PyYAML, LinkML (`linkml`, `linkml-term-validator`, `linkml-reference-validator`), oaklib (GO aspect lookup), `genesets-rs` Rust CLI (enrichment), pytest. Package + deps managed with `uv`; tasks exposed via `just`.

**Conventions:** Per the user's global instructions — `uv`, Typer (not argparse) for new CLI code, `just` targets, pytest (no unittest), doctests, CLI tests, no try/except workarounds. The existing `genesets-workflows` commands use argparse; we deliberately use Typer for the new `curate` group (honoring the global "never low-level argparse" rule) and bridge it into the existing dispatcher with a thin `main(argv)` wrapper. We do **not** rewrite existing commands.

**Key validator CLIs (verified against installed tools):**
- `linkml-validate --schema SCHEMA --target-class CLASS FILE`
- `linkml-term-validator validate-data FILE -s SCHEMA -t CLASS --labels -c OAK_CONFIG`
- `linkml-term-validator validate-schema SCHEMA -c OAK_CONFIG`
- `linkml-reference-validator validate data FILE --schema SCHEMA --target-class CLASS --cache-dir CACHE_DIR`
- `linkml-reference-validator lookup PMID:XXXX` (prints cached title/abstract — used to source exact snippets)

**Enrichment TSV header (from `src/output.rs`):** `query_id, query_name, target_id, target_name, overlap, query_size, target_size, background_size, p_value, p_adjust[ _bonferroni], [overlap_genes, overlap_gene_names]`. `overlap_genes`/`overlap_gene_names` are `;`-separated.

**Test command (used throughout):**
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests -v
```

---

## File Structure

Created/modified by this plan:

```
curation/
  README.md                                  # Task 13
  schema/genesets_interpretation.yaml        # Task 2  (LinkML schema, source of truth)
  conf/oak_config.yaml                       # Task 3  (OAK adapters per prefix)
  conf/reference_validator_config.yaml       # Task 3
  c8/manifest.tsv                            # Task 13 (seed list of C8 sets)
  c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml     # Task 4/12 (canonical example interpretation)
  cache/.gitkeep                             # Task 3  (committed ontology label cache lands here)
python/genesets-workflows/
  pyproject.toml                             # Task 1  (add `curation` extra)
  src/genesets_workflows/cli.py              # Task 11 (wire `curate` into dispatcher)
  src/genesets_workflows/curation/__init__.py# Task 1
  src/genesets_workflows/curation/model.py   # Task 5  (dataclasses + YAML round-trip)
  src/genesets_workflows/curation/enrichment.py # Task 6 (TSV parse + subprocess runner)
  src/genesets_workflows/curation/aspect.py  # Task 7  (GO aspect from ancestors)
  src/genesets_workflows/curation/draft.py   # Task 8  (assemble draft interpretation)
  src/genesets_workflows/curation/metrics.py # Task 9  (precision/recall/F1)
  src/genesets_workflows/curation/report.py  # Task 10 (join tree vs enrichment)
  src/genesets_workflows/curation/validate.py# Task 12 (validator command builders)
  src/genesets_workflows/curation/cli.py     # Task 11 (Typer app: draft/validate/report)
  tests/                                      # Tasks 5-12 (pytest)
justfile                                     # Task 12 (curate-* targets)
docs/src/curation.md + docs/src/SUMMARY.md   # Task 13
README.md                                    # Task 13
.gitignore                                   # Task 3 (ignore curation/references_cache)
```

---

## Task 1: Add the `curation` extra and package skeleton

**Files:**
- Modify: `python/genesets-workflows/pyproject.toml`
- Create: `python/genesets-workflows/src/genesets_workflows/curation/__init__.py`
- Create: `python/genesets-workflows/tests/test_curation_imports.py`

- [ ] **Step 1: Add the optional-dependency extra**

In `python/genesets-workflows/pyproject.toml`, under `[project.optional-dependencies]` (after the existing `notebooks = [...]` block), add:

```toml
curation = [
  "typer>=0.12",
  "linkml>=1.8",
  "linkml-term-validator>=0.4.0",
  "linkml-reference-validator>=0.2",
  "oaklib>=0.6",
]
```

- [ ] **Step 2: Create the subpackage**

Create `python/genesets-workflows/src/genesets_workflows/curation/__init__.py`:

```python
"""Curation of GO-term interpretations for non-GO gene sets."""

from __future__ import annotations

__all__ = []
```

- [ ] **Step 3: Write the import smoke test**

Create `python/genesets-workflows/tests/test_curation_imports.py`:

```python
def test_curation_package_imports():
    import genesets_workflows.curation as curation

    assert curation.__doc__
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_curation_imports.py -v
```
Expected: PASS (1 passed). This also confirms the `curation` extra resolves and installs.

- [ ] **Step 5: Commit**

```bash
git add python/genesets-workflows/pyproject.toml \
        python/genesets-workflows/src/genesets_workflows/curation/__init__.py \
        python/genesets-workflows/tests/test_curation_imports.py
git commit -m "feat(curation): add curation subpackage and extra deps"
```

---

## Task 2: Write the LinkML schema

**Files:**
- Create: `curation/schema/genesets_interpretation.yaml`
- Create: `python/genesets-workflows/tests/test_schema_loads.py`

- [ ] **Step 1: Write the failing schema-load test**

Create `python/genesets-workflows/tests/test_schema_loads.py`:

```python
from pathlib import Path

from linkml_runtime import SchemaView

SCHEMA = Path(__file__).resolve().parents[3] / "curation" / "schema" / "genesets_interpretation.yaml"


def test_schema_has_expected_classes_and_enum():
    view = SchemaView(str(SCHEMA))
    classes = set(view.all_classes())
    assert {"GeneSetInterpretation", "TermAssociation", "Term", "BiologicalContext"} <= classes
    category = view.get_enum("AssociationCategoryEnum")
    values = set(category.permissible_values)
    assert values == {
        "core_process",
        "core_component",
        "supporting_process",
        "marker_driven_plausible",
        "nonspecific",
        "false_association",
    }


def test_tree_root_is_interpretation():
    view = SchemaView(str(SCHEMA))
    roots = [c.name for c in view.all_classes().values() if c.tree_root]
    assert roots == ["GeneSetInterpretation"]
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_schema_loads.py -v
```
Expected: FAIL (schema file does not exist).

- [ ] **Step 3: Write the schema**

Create `curation/schema/genesets_interpretation.yaml`:

```yaml
id: https://w3id.org/genesets-rs/genesets-interpretation
name: genesets_interpretation
title: Gene Set GO Interpretation Schema
description: >-
  Curated, validated GO-term interpretations of non-GO gene sets (e.g. MSigDB
  C8 cell-type signatures). Each interpretation records, per GO term, a curator
  judgment (role category, confidence, specificity), optional literature
  evidence, and the enrichment provenance that seeded it.
license: MIT
prefixes:
  linkml: https://w3id.org/linkml/
  genesets: https://w3id.org/genesets-rs/genesets-interpretation/
  GO: http://purl.obolibrary.org/obo/GO_
  CL: http://purl.obolibrary.org/obo/CL_
  UBERON: http://purl.obolibrary.org/obo/UBERON_
  MONDO: http://purl.obolibrary.org/obo/MONDO_
  CHEBI: http://purl.obolibrary.org/obo/CHEBI_
  PR: http://purl.obolibrary.org/obo/PR_
  HP: http://purl.obolibrary.org/obo/HP_
  NCBITaxon: http://purl.obolibrary.org/obo/NCBITaxon_
  EFO: http://www.ebi.ac.uk/efo/EFO_
  MSIGDB: https://www.gsea-msigdb.org/gsea/msigdb/cards/
  PMID: http://identifiers.org/pubmed/
  DOI: https://doi.org/
default_prefix: genesets
default_range: string
imports:
  - linkml:types

classes:

  GeneSetInterpretation:
    tree_root: true
    description: A curated GO interpretation of a single non-GO gene set.
    attributes:
      gene_set_id:
        identifier: true
        range: uriorcurie
        required: true
        description: Stable identifier for the gene set (e.g. MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL).
      gene_set_name:
        required: true
        description: Human-readable name of the gene set.
      collection:
        description: Source collection, e.g. C8, C2:CP:KEGG_LEGACY, C7.
      msigdb_release:
        description: MSigDB release the gene set came from, e.g. v2025.1.Hs.
      taxon:
        range: Term
        inlined: true
        description: Species the gene set is defined in (an NCBITaxon term).
      direction:
        range: DirectionEnum
        description: Whether the signature is the up, down, or undirected component.
      n_genes:
        range: integer
        description: Number of genes in the source set.
      description:
        description: Free-text summary of what the signature represents.
      contexts:
        range: BiologicalContext
        multivalued: true
        inlined_as_list: true
        description: What the signature is about (cell type, disease, chemical, etc.).
      provenance:
        range: EnrichmentProvenance
        inlined: true
        description: How the draft interpretation was generated.
      curation_status:
        range: CurationStatusEnum
      curator:
      curation_date:
        range: date
      associations:
        range: TermAssociation
        multivalued: true
        inlined_as_list: true
        recommended: true
        description: Per-GO-term curator judgments.

  BiologicalContext:
    description: An ontology term describing what the gene set is about.
    attributes:
      term:
        range: Term
        inlined: true
        required: true
      context_type:
        range: ContextTypeEnum
        required: true
      role_note:
        description: Optional note on the term's role, e.g. "treatment", "vs comparator".

  TermAssociation:
    description: A curated association between the gene set and a single GO term.
    attributes:
      term:
        range: Term
        inlined: true
        required: true
        description: The GO term.
      aspect:
        range: GOAspectEnum
        description: The GO aspect (BP/CC/MF) of the term.
      category:
        range: AssociationCategoryEnum
        recommended: true
        description: >-
          The curator's role judgment. Unset on freshly drafted (un-adjudicated)
          associations.
      confidence:
        range: ConfidenceLevelEnum
      confidence_score:
        range: float
        minimum_value: 0.0
        maximum_value: 1.0
      specificity:
        range: SpecificityEnum
        description: How specific the term is to this biological context.
      specificity_score:
        range: float
        minimum_value: 0.0
        maximum_value: 1.0
      seed_source:
        range: SeedSourceEnum
        required: true
        description: Whether the term was recovered by enrichment or added by a curator.
      enrichment_stats:
        range: EnrichmentStats
        inlined: true
      evidence:
        range: EvidenceItem
        multivalued: true
        inlined_as_list: true
        recommended: true
      curator_note:

  Term:
    description: An ontology term reference (validated by linkml-term-validator).
    attributes:
      id:
        identifier: true
        range: uriorcurie
        required: true
      label:
        required: true
        comments:
          - This is automatically validated by the linkml-term-validator tool.

  EnrichmentStats:
    description: Enrichment numbers that seeded an association.
    attributes:
      p_value:
        range: float
      p_adjust:
        range: float
      overlap:
        range: integer
      query_size:
        range: integer
      target_size:
        range: integer
      overlap_genes:
        multivalued: true

  EnrichmentProvenance:
    description: How a draft interpretation was generated.
    attributes:
      tool:
        description: Tool used to seed the draft, e.g. genesets-rs.
      config_ref:
        description: Reference to the eval/run config (path or commit) used.
      go_version:
        description: GO snapshot/version identifier used for enrichment.
      evidence_codes:
        description: Which GO evidence set was used (e.g. all, no-IBA).
      background:
        description: Background gene universe identifier.
      run_date:
        range: date

  EvidenceItem:
    description: A literature evidence item supporting or refuting an association.
    attributes:
      reference:
        range: uriorcurie
        description: Authoritative reference (PMID/DOI).
        implements:
          - linkml:authoritative_reference
      reference_title:
        recommended: true
      supports:
        range: EvidenceSupportEnum
      snippet:
        description: Exact excerpt from the reference supporting/refuting the claim.
        comments:
          - This is automatically validated by the linkml-reference-validator tool.
        implements:
          - linkml:excerpt
      explanation:

enums:

  DirectionEnum:
    permissible_values:
      UP:
      DN:
      NA:

  GOAspectEnum:
    permissible_values:
      biological_process:
        meaning: GO:0008150
      cellular_component:
        meaning: GO:0005575
      molecular_function:
        meaning: GO:0003674

  ContextTypeEnum:
    permissible_values:
      cell_type:
      anatomical_structure:
      disease:
      phenotype:
      chemical_perturbagen:
      genetic_perturbation:
      infection:
      developmental_stage:
      cell_line:
      experimental_condition:

  AssociationCategoryEnum:
    permissible_values:
      core_process:
        description: Central biological process this context actually executes.
      core_component:
        description: Characteristic cellular component or structure.
      supporting_process:
        description: Real enabling/supporting process, not defining.
      marker_driven_plausible:
        description: Plausible given the marker genes but peripheral or generic.
      nonspecific:
        description: Housekeeping or broadly expressed; not informative here.
      false_association:
        description: Artifactual or biologically wrong for this set.

  ConfidenceLevelEnum:
    permissible_values:
      high:
      medium:
      low:

  SpecificityEnum:
    permissible_values:
      cell_type_specific:
      lineage_shared:
      broadly_expressed:
      ubiquitous:

  SeedSourceEnum:
    permissible_values:
      enrichment_recovered:
      curator_added:

  EvidenceSupportEnum:
    permissible_values:
      SUPPORT:
      REFUTE:
      NEUTRAL:

  CurationStatusEnum:
    permissible_values:
      draft:
      in_review:
      reviewed:
      final:
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_schema_loads.py -v
```
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add curation/schema/genesets_interpretation.yaml python/genesets-workflows/tests/test_schema_loads.py
git commit -m "feat(curation): add LinkML schema for gene set GO interpretations"
```

---

## Task 3: OAK + reference-validator config and cache scaffolding

**Files:**
- Create: `curation/conf/oak_config.yaml`
- Create: `curation/conf/reference_validator_config.yaml`
- Create: `curation/cache/.gitkeep`
- Modify: `.gitignore`

- [ ] **Step 1: Write the OAK config**

Create `curation/conf/oak_config.yaml` (maps CURIE prefixes to OAK adapters; only listed prefixes are validated):

```yaml
# OAK configuration for gene set interpretation term validation.
# Maps ontology prefixes to OAK adapter strings used by linkml-term-validator.
# Add prefixes here as new biological contexts require them.
ontology_adapters:
  GO: sqlite:obo:go
  CL: sqlite:obo:cl
  UBERON: sqlite:obo:uberon
  MONDO: sqlite:obo:mondo
  CHEBI: sqlite:obo:chebi
  PR: sqlite:obo:pr
  HP: sqlite:obo:hp
  NCBITaxon: sqlite:obo:ncbitaxon
```

- [ ] **Step 2: Write the reference-validator config**

Create `curation/conf/reference_validator_config.yaml`:

```yaml
# Configuration for linkml-reference-validator.
# Snippets are matched against cited references by deterministic substring match.
validation:
  cache_dir: references_cache
  unknown_prefix_severity: INFO
```

- [ ] **Step 3: Keep the committed label cache directory**

Create `curation/cache/.gitkeep` (empty file). The per-ontology label cache (e.g. `curation/cache/go/terms.csv`) is written here by `linkml-term-validator` and committed so reviewers see label drift.

- [ ] **Step 4: Ignore the references cache**

Append to `.gitignore`:

```
# Reference-validator paper cache (rebuilt on demand)
curation/references_cache/
```

- [ ] **Step 5: Verify configs parse**

Run:
```bash
uv run --project python/genesets-workflows --extra curation python -c "import yaml,pathlib; [yaml.safe_load(pathlib.Path(p).read_text()) for p in ['curation/conf/oak_config.yaml','curation/conf/reference_validator_config.yaml']]; print('configs ok')"
```
Expected: prints `configs ok`.

- [ ] **Step 6: Commit**

```bash
git add curation/conf/oak_config.yaml curation/conf/reference_validator_config.yaml curation/cache/.gitkeep .gitignore
git commit -m "feat(curation): add OAK + reference-validator config and cache scaffolding"
```

---

## Task 4: Hand-author the canonical example interpretation (no evidence yet)

**Files:**
- Create: `curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml`
- Create: `python/genesets-workflows/tests/test_example_valid.py`

- [ ] **Step 1: Write the failing structural test**

Create `python/genesets-workflows/tests/test_example_valid.py`:

```python
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "curation" / "c8" / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"


def test_example_loads_and_has_required_fields():
    data = yaml.safe_load(EXAMPLE.read_text())
    assert data["gene_set_id"] == "MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"
    assert data["contexts"], "must record at least one biological context"
    assert data["associations"], "must record at least one association"
    for assoc in data["associations"]:
        assert assoc["term"]["id"].startswith("GO:")
        assert assoc["term"]["label"]
        assert assoc["seed_source"] in {"enrichment_recovered", "curator_added"}


def test_example_categories_are_in_enum():
    data = yaml.safe_load(EXAMPLE.read_text())
    allowed = {
        "core_process",
        "core_component",
        "supporting_process",
        "marker_driven_plausible",
        "nonspecific",
        "false_association",
    }
    for assoc in data["associations"]:
        if "category" in assoc:
            assert assoc["category"] in allowed
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_example_valid.py -v
```
Expected: FAIL (example file does not exist).

- [ ] **Step 3: Write the example interpretation**

Create `curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml`. (GO IDs/labels are real and will be checked in Step 5; evidence is added in Task 12.)

```yaml
gene_set_id: MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL
gene_set_name: HAY_BONE_MARROW_DENDRITIC_CELL
collection: C8
msigdb_release: v2025.1.Hs
taxon:
  id: NCBITaxon:9606
  label: Homo sapiens
direction: NA
description: >-
  Single-cell marker signature of dendritic cells from human bone marrow
  (Hay et al.). Used here as a worked example of a curated GO interpretation.
contexts:
  - term:
      id: CL:0000451
      label: dendritic cell
    context_type: cell_type
  - term:
      id: UBERON:0002371
      label: bone marrow
    context_type: anatomical_structure
provenance:
  tool: genesets-rs
  config_ref: curation example (hand-authored)
  go_version: pinned-snapshot-TBD
  evidence_codes: all
  background: goa-human-symbols
curation_status: reviewed
curator: example
curation_date: 2026-06-24
associations:
  - term:
      id: GO:0019882
      label: antigen processing and presentation
    aspect: biological_process
    category: core_process
    confidence: high
    specificity: cell_type_specific
    seed_source: enrichment_recovered
  - term:
      id: GO:0042613
      label: MHC class II protein complex
    aspect: cellular_component
    category: core_component
    confidence: high
    specificity: cell_type_specific
    seed_source: enrichment_recovered
  - term:
      id: GO:0045087
      label: innate immune response
    aspect: biological_process
    category: supporting_process
    confidence: medium
    specificity: lineage_shared
    seed_source: enrichment_recovered
  - term:
      id: GO:0050900
      label: leukocyte migration
    aspect: biological_process
    category: supporting_process
    confidence: medium
    specificity: lineage_shared
    seed_source: enrichment_recovered
  - term:
      id: GO:0006412
      label: translation
    aspect: biological_process
    category: nonspecific
    confidence: high
    specificity: ubiquitous
    seed_source: enrichment_recovered
  - term:
      id: GO:0002250
      label: adaptive immune response
    aspect: biological_process
    category: core_process
    confidence: medium
    specificity: lineage_shared
    seed_source: curator_added
```

Note: `go_version: pinned-snapshot-TBD` is a free-text provenance note for a hand-authored example, not a schema placeholder; real drafts (Task 8) fill it from the actual run.

- [ ] **Step 4: Run the structural test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_example_valid.py -v
```
Expected: PASS (2 passed).

- [ ] **Step 5: Validate structure + terms with the real validators**

Run (this downloads the GO/CL/UBERON/NCBITaxon SQLite DBs on first run — may take minutes):
```bash
uv run --project python/genesets-workflows --extra curation \
  linkml-validate --schema curation/schema/genesets_interpretation.yaml \
  --target-class GeneSetInterpretation curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml

uv run --project python/genesets-workflows --extra curation \
  linkml-term-validator validate-data curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml \
  -s curation/schema/genesets_interpretation.yaml -t GeneSetInterpretation \
  --labels -c curation/conf/oak_config.yaml
```
Expected: `linkml-validate` reports "No issues found"; term-validator reports all terms valid with matching labels. If any GO/CL/UBERON label mismatch is reported, correct the `label:` in the YAML to the validator's reported label and re-run.

- [ ] **Step 6: Commit (including the populated label cache)**

```bash
git add curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml \
        python/genesets-workflows/tests/test_example_valid.py \
        curation/cache
git commit -m "feat(curation): add canonical dendritic-cell example interpretation"
```

---

## Task 5: Data model (dataclasses + YAML round-trip)

**Files:**
- Create: `python/genesets-workflows/src/genesets_workflows/curation/model.py`
- Create: `python/genesets-workflows/tests/test_model.py`

- [ ] **Step 1: Write the failing round-trip test**

Create `python/genesets-workflows/tests/test_model.py`:

```python
from pathlib import Path

from genesets_workflows.curation import model

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "curation" / "c8" / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"


def test_load_example():
    interp = model.load_interpretation(EXAMPLE)
    assert interp.gene_set_id == "MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"
    assert len(interp.associations) == 6
    first = interp.associations[0]
    assert first.term.id == "GO:0019882"
    assert first.category == "core_process"
    assert first.seed_source == "enrichment_recovered"


def test_roundtrip_drops_none(tmp_path):
    interp = model.load_interpretation(EXAMPLE)
    out = tmp_path / "rt.yaml"
    model.dump_interpretation(interp, out)
    text = out.read_text()
    assert "null" not in text  # None-valued optional fields are omitted
    reloaded = model.load_interpretation(out)
    assert reloaded.gene_set_id == interp.gene_set_id
    assert len(reloaded.associations) == len(interp.associations)


def test_build_minimal_interpretation():
    interp = model.GeneSetInterpretation(
        gene_set_id="MSIGDB:X",
        gene_set_name="X",
        associations=[
            model.TermAssociation(
                term=model.Term(id="GO:0006955", label="immune response"),
                seed_source="enrichment_recovered",
            )
        ],
    )
    assert interp.to_dict()["associations"][0]["term"]["id"] == "GO:0006955"
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_model.py -v
```
Expected: FAIL (`model` has no attribute `load_interpretation`).

- [ ] **Step 3: Write the model**

Create `python/genesets-workflows/src/genesets_workflows/curation/model.py`:

```python
"""Dataclasses mirroring the genesets_interpretation LinkML schema, with
YAML round-trip helpers. Optional (None / empty) fields are omitted on dump."""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass
class Term:
    id: str
    label: str


@dataclass
class EnrichmentStats:
    p_value: float | None = None
    p_adjust: float | None = None
    overlap: int | None = None
    query_size: int | None = None
    target_size: int | None = None
    overlap_genes: list[str] = field(default_factory=list)


@dataclass
class EvidenceItem:
    reference: str | None = None
    reference_title: str | None = None
    supports: str | None = None
    snippet: str | None = None
    explanation: str | None = None


@dataclass
class TermAssociation:
    term: Term
    seed_source: str
    aspect: str | None = None
    category: str | None = None
    confidence: str | None = None
    confidence_score: float | None = None
    specificity: str | None = None
    specificity_score: float | None = None
    enrichment_stats: EnrichmentStats | None = None
    evidence: list[EvidenceItem] = field(default_factory=list)
    curator_note: str | None = None


@dataclass
class BiologicalContext:
    term: Term
    context_type: str
    role_note: str | None = None


@dataclass
class EnrichmentProvenance:
    tool: str | None = None
    config_ref: str | None = None
    go_version: str | None = None
    evidence_codes: str | None = None
    background: str | None = None
    run_date: str | None = None


@dataclass
class GeneSetInterpretation:
    gene_set_id: str
    gene_set_name: str
    collection: str | None = None
    msigdb_release: str | None = None
    taxon: Term | None = None
    direction: str | None = None
    n_genes: int | None = None
    description: str | None = None
    contexts: list[BiologicalContext] = field(default_factory=list)
    provenance: EnrichmentProvenance | None = None
    curation_status: str | None = None
    curator: str | None = None
    curation_date: str | None = None
    associations: list[TermAssociation] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return _to_dict(self)


def _to_dict(obj: Any) -> Any:
    """Recursively convert dataclasses to plain dicts, dropping None and empties."""
    if is_dataclass(obj):
        result: dict[str, Any] = {}
        for f in fields(obj):
            value = _to_dict(getattr(obj, f.name))
            if value is None:
                continue
            if isinstance(value, (list, dict)) and len(value) == 0:
                continue
            result[f.name] = value
        return result
    if isinstance(obj, list):
        return [_to_dict(item) for item in obj]
    return obj


# Mapping of dataclass type -> nested dataclass fields needing reconstruction.
_TERM_FIELDS = {"term", "taxon"}


def _term(value: dict[str, Any] | None) -> Term | None:
    return None if value is None else Term(**value)


def _association(value: dict[str, Any]) -> TermAssociation:
    value = dict(value)
    value["term"] = Term(**value["term"])
    if "enrichment_stats" in value and value["enrichment_stats"] is not None:
        value["enrichment_stats"] = EnrichmentStats(**value["enrichment_stats"])
    if value.get("evidence"):
        value["evidence"] = [EvidenceItem(**e) for e in value["evidence"]]
    return TermAssociation(**value)


def _context(value: dict[str, Any]) -> BiologicalContext:
    value = dict(value)
    value["term"] = Term(**value["term"])
    return BiologicalContext(**value)


def from_dict(data: dict[str, Any]) -> GeneSetInterpretation:
    data = dict(data)
    data["taxon"] = _term(data.get("taxon"))
    if data.get("provenance") is not None:
        data["provenance"] = EnrichmentProvenance(**data["provenance"])
    data["contexts"] = [_context(c) for c in data.get("contexts", [])]
    data["associations"] = [_association(a) for a in data.get("associations", [])]
    return GeneSetInterpretation(**data)


def load_interpretation(path: str | Path) -> GeneSetInterpretation:
    data = yaml.safe_load(Path(path).read_text())
    return from_dict(data)


def dump_interpretation(interp: GeneSetInterpretation, path: str | Path) -> None:
    text = yaml.safe_dump(interp.to_dict(), sort_keys=False, allow_unicode=True)
    Path(path).write_text(text)
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_model.py -v
```
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add python/genesets-workflows/src/genesets_workflows/curation/model.py \
        python/genesets-workflows/tests/test_model.py
git commit -m "feat(curation): add dataclass model with YAML round-trip"
```

---

## Task 6: Enrichment TSV parser + subprocess runner

**Files:**
- Create: `python/genesets-workflows/src/genesets_workflows/curation/enrichment.py`
- Create: `python/genesets-workflows/tests/test_enrichment.py`
- Create: `python/genesets-workflows/tests/data/enrich_sample.tsv`

- [ ] **Step 1: Create the fixture TSV**

Create `python/genesets-workflows/tests/data/enrich_sample.tsv` (tab-separated; matches the real header):

```
query_id	query_name	target_id	target_name	overlap	query_size	target_size	background_size	p_value	p_adjust	overlap_genes	overlap_gene_names
dc_up	dc_up	GO:0019882	antigen processing and presentation	12	100	340	20000	1.2e-09	3.4e-07	HLA-DRA;CD74	HLA-DRA;CD74
dc_up	dc_up	GO:0006412	translation	8	100	1500	20000	2.0e-02	8.0e-01	RPL3;RPS6	RPL3;RPS6
```

(Ensure real tab characters separate the columns, not spaces.)

- [ ] **Step 2: Write the failing parser test**

Create `python/genesets-workflows/tests/test_enrichment.py`:

```python
from pathlib import Path

from genesets_workflows.curation import enrichment

DATA = Path(__file__).resolve().parent / "data" / "enrich_sample.tsv"


def test_parse_enrichment_tsv():
    rows = enrichment.parse_enrichment_tsv(DATA.read_text())
    assert len(rows) == 2
    first = rows[0]
    assert first.target_id == "GO:0019882"
    assert first.target_name == "antigen processing and presentation"
    assert first.overlap == 12
    assert first.p_adjust == 3.4e-07
    assert first.overlap_genes == ["HLA-DRA", "CD74"]


def test_parse_handles_bonferroni_column():
    text = "query_id\tquery_name\ttarget_id\ttarget_name\toverlap\tquery_size\ttarget_size\tbackground_size\tp_value\tp_adjust_bonferroni\nq\tq\tGO:1\tt\t1\t2\t3\t4\t0.1\t0.5\n"
    rows = enrichment.parse_enrichment_tsv(text)
    assert rows[0].p_adjust == 0.5
    assert rows[0].overlap_genes == []
```

- [ ] **Step 3: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_enrichment.py -v
```
Expected: FAIL (`enrichment` module missing).

- [ ] **Step 4: Write the enrichment module**

Create `python/genesets-workflows/src/genesets_workflows/curation/enrichment.py`:

```python
"""Parse genesets-rs enrichment TSV and run the Rust enrichment kernel."""

from __future__ import annotations

import csv
import io
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class EnrichmentRow:
    query_id: str
    target_id: str
    target_name: str
    overlap: int
    query_size: int
    target_size: int
    p_value: float | None
    p_adjust: float | None
    overlap_genes: list[str]


def _to_float(value: str) -> float | None:
    value = value.strip()
    if value == "":
        return None
    return float(value)


def parse_enrichment_tsv(text: str) -> list[EnrichmentRow]:
    """Parse genesets-rs enrich/matrix TSV output into EnrichmentRow records.

    The adjusted-p column may be named ``p_adjust`` or ``p_adjust_bonferroni``.
    ``overlap_genes`` is optional and ``;``-separated.
    """
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    p_adjust_key = "p_adjust_bonferroni" if "p_adjust_bonferroni" in (reader.fieldnames or []) else "p_adjust"
    rows: list[EnrichmentRow] = []
    for record in reader:
        genes_raw = record.get("overlap_genes") or ""
        rows.append(
            EnrichmentRow(
                query_id=record["query_id"],
                target_id=record["target_id"],
                target_name=record.get("target_name", ""),
                overlap=int(record["overlap"]),
                query_size=int(record["query_size"]),
                target_size=int(record["target_size"]),
                p_value=_to_float(record.get("p_value", "")),
                p_adjust=_to_float(record.get(p_adjust_key, "")),
                overlap_genes=[g for g in genes_raw.split(";") if g],
            )
        )
    return rows


def run_enrichment(
    *,
    sample: Path,
    annotations: Path,
    terms: Path,
    closure: Path,
    background: Path,
    min_overlap: int = 3,
    max_p_adjust: float = 0.05,
    binary: str = "genesets-rs",
) -> str:
    """Run `genesets-rs enrich` and return its TSV stdout."""
    cmd = [
        binary,
        "enrich",
        "--annotations", str(annotations),
        "--terms", str(terms),
        "--closure", str(closure),
        "--sample", str(sample),
        "--background", str(background),
        "--min-overlap", str(min_overlap),
        "--max-p-adjust", str(max_p_adjust),
        "--overlap-genes",
    ]
    completed = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return completed.stdout
```

- [ ] **Step 5: Run the parser test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_enrichment.py -v
```
Expected: PASS (2 passed).

- [ ] **Step 6: Write a subprocess-runner test with a stub binary**

Append to `python/genesets-workflows/tests/test_enrichment.py`:

```python
import os
import stat


def test_run_enrichment_invokes_binary(tmp_path):
    stub = tmp_path / "genesets-rs"
    stub.write_text("#!/usr/bin/env bash\ncat " + str(DATA) + "\n")
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    for name in ("sample", "annotations", "terms", "closure", "background"):
        (tmp_path / name).write_text("x\n")
    out = enrichment.run_enrichment(
        sample=tmp_path / "sample",
        annotations=tmp_path / "annotations",
        terms=tmp_path / "terms",
        closure=tmp_path / "closure",
        background=tmp_path / "background",
        binary=str(stub),
    )
    rows = enrichment.parse_enrichment_tsv(out)
    assert rows[0].target_id == "GO:0019882"
```

- [ ] **Step 7: Run the full enrichment test file to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_enrichment.py -v
```
Expected: PASS (3 passed).

- [ ] **Step 8: Commit**

```bash
git add python/genesets-workflows/src/genesets_workflows/curation/enrichment.py \
        python/genesets-workflows/tests/test_enrichment.py \
        python/genesets-workflows/tests/data/enrich_sample.tsv
git commit -m "feat(curation): add enrichment TSV parser and runner"
```

---

## Task 7: GO aspect derivation

**Files:**
- Create: `python/genesets-workflows/src/genesets_workflows/curation/aspect.py`
- Create: `python/genesets-workflows/tests/test_aspect.py`

- [ ] **Step 1: Write the failing test**

Create `python/genesets-workflows/tests/test_aspect.py`:

```python
from genesets_workflows.curation import aspect


def test_aspect_from_ancestors_bp():
    assert aspect.aspect_from_ancestors("GO:0019882", {"GO:0008150", "GO:0019882"}) == "biological_process"


def test_aspect_from_ancestors_cc():
    assert aspect.aspect_from_ancestors("GO:0042613", {"GO:0005575", "GO:0032991"}) == "cellular_component"


def test_aspect_from_ancestors_unknown():
    assert aspect.aspect_from_ancestors("GO:9999999", set()) is None


def test_aspect_for_term_uses_lookup():
    def fake_ancestors(go_id):
        return {"GO:0003674", go_id}

    assert aspect.aspect_for_term("GO:0004672", fake_ancestors) == "molecular_function"
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_aspect.py -v
```
Expected: FAIL (`aspect` module missing).

- [ ] **Step 3: Write the aspect module**

Create `python/genesets-workflows/src/genesets_workflows/curation/aspect.py`:

```python
"""Derive the GO aspect (BP/CC/MF) of a term from its ancestors."""

from __future__ import annotations

from typing import Callable, Iterable

ASPECT_ROOTS = {
    "GO:0008150": "biological_process",
    "GO:0005575": "cellular_component",
    "GO:0003674": "molecular_function",
}


def aspect_from_ancestors(go_id: str, ancestors: Iterable[str]) -> str | None:
    """Return the GO aspect implied by an ancestor set, or None if undetermined.

    >>> aspect_from_ancestors("GO:0019882", {"GO:0008150"})
    'biological_process'
    """
    ancestor_set = set(ancestors)
    for root, name in ASPECT_ROOTS.items():
        if root in ancestor_set:
            return name
    return None


def aspect_for_term(go_id: str, ancestors_lookup: Callable[[str], Iterable[str]]) -> str | None:
    """Derive aspect using a caller-supplied ancestors lookup (e.g. an OAK adapter)."""
    return aspect_from_ancestors(go_id, ancestors_lookup(go_id))


def oak_ancestors_lookup(adapter):  # pragma: no cover - exercised in integration runs
    """Build an ancestors lookup backed by an oaklib adapter (IS_A only)."""
    from oaklib.datamodels.vocabulary import IS_A

    def lookup(go_id: str):
        return set(adapter.ancestors(go_id, predicates=[IS_A]))

    return lookup
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_aspect.py -v
```
Expected: PASS (4 passed).

- [ ] **Step 5: Run the doctest to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest python -m pytest --doctest-modules python/genesets-workflows/src/genesets_workflows/curation/aspect.py -v
```
Expected: PASS (1 doctest passed).

- [ ] **Step 6: Commit**

```bash
git add python/genesets-workflows/src/genesets_workflows/curation/aspect.py \
        python/genesets-workflows/tests/test_aspect.py
git commit -m "feat(curation): add GO aspect derivation"
```

---

## Task 8: Draft assembly

**Files:**
- Create: `python/genesets-workflows/src/genesets_workflows/curation/draft.py`
- Create: `python/genesets-workflows/tests/test_draft.py`

- [ ] **Step 1: Write the failing test**

Create `python/genesets-workflows/tests/test_draft.py`:

```python
from genesets_workflows.curation import draft, model
from genesets_workflows.curation.enrichment import EnrichmentRow


def _rows():
    return [
        EnrichmentRow("dc", "GO:0019882", "antigen processing and presentation", 12, 100, 340, 1e-9, 3e-7, ["HLA-DRA"]),
        EnrichmentRow("dc", "GO:0006412", "translation", 8, 100, 1500, 2e-2, 8e-1, ["RPL3"]),
    ]


def test_assemble_draft_basic():
    interp = draft.assemble_draft(
        gene_set_id="MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL",
        gene_set_name="HAY_BONE_MARROW_DENDRITIC_CELL",
        collection="C8",
        rows=_rows(),
        provenance=model.EnrichmentProvenance(tool="genesets-rs", go_version="snap-X"),
    )
    assert interp.gene_set_id == "MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"
    assert interp.curation_status == "draft"
    assert len(interp.associations) == 2
    a0 = interp.associations[0]
    assert a0.term.id == "GO:0019882"
    assert a0.seed_source == "enrichment_recovered"
    assert a0.category is None  # un-adjudicated
    assert a0.enrichment_stats.overlap == 12
    assert a0.enrichment_stats.p_adjust == 3e-7


def test_assemble_draft_applies_aspect_lookup():
    interp = draft.assemble_draft(
        gene_set_id="MSIGDB:X",
        gene_set_name="X",
        collection="C8",
        rows=_rows()[:1],
        provenance=model.EnrichmentProvenance(),
        aspect_lookup=lambda go_id: "biological_process",
    )
    assert interp.associations[0].aspect == "biological_process"
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_draft.py -v
```
Expected: FAIL (`draft` module missing).

- [ ] **Step 3: Write the draft module**

Create `python/genesets-workflows/src/genesets_workflows/curation/draft.py`:

```python
"""Assemble a draft GeneSetInterpretation from enrichment output."""

from __future__ import annotations

from typing import Callable, Iterable

from genesets_workflows.curation import model
from genesets_workflows.curation.enrichment import EnrichmentRow


def assemble_draft(
    *,
    gene_set_id: str,
    gene_set_name: str,
    collection: str,
    rows: Iterable[EnrichmentRow],
    provenance: model.EnrichmentProvenance,
    contexts: list[model.BiologicalContext] | None = None,
    aspect_lookup: Callable[[str], str | None] | None = None,
    msigdb_release: str | None = None,
    n_genes: int | None = None,
) -> model.GeneSetInterpretation:
    """Build an un-adjudicated draft: every enrichment hit becomes a
    TermAssociation with seed_source=enrichment_recovered and no category."""
    associations: list[model.TermAssociation] = []
    for row in rows:
        aspect = aspect_lookup(row.target_id) if aspect_lookup else None
        associations.append(
            model.TermAssociation(
                term=model.Term(id=row.target_id, label=row.target_name),
                seed_source="enrichment_recovered",
                aspect=aspect,
                enrichment_stats=model.EnrichmentStats(
                    p_value=row.p_value,
                    p_adjust=row.p_adjust,
                    overlap=row.overlap,
                    query_size=row.query_size,
                    target_size=row.target_size,
                    overlap_genes=list(row.overlap_genes),
                ),
            )
        )
    return model.GeneSetInterpretation(
        gene_set_id=gene_set_id,
        gene_set_name=gene_set_name,
        collection=collection,
        msigdb_release=msigdb_release,
        n_genes=n_genes,
        contexts=contexts or [],
        provenance=provenance,
        curation_status="draft",
        associations=associations,
    )
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_draft.py -v
```
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add python/genesets-workflows/src/genesets_workflows/curation/draft.py \
        python/genesets-workflows/tests/test_draft.py
git commit -m "feat(curation): add draft interpretation assembly"
```

---

## Task 9: Precision/recall/F1 metrics

**Files:**
- Create: `python/genesets-workflows/src/genesets_workflows/curation/metrics.py`
- Create: `python/genesets-workflows/tests/test_metrics.py`

- [ ] **Step 1: Write the failing test**

Create `python/genesets-workflows/tests/test_metrics.py`:

```python
from genesets_workflows.curation import metrics, model


def _assoc(go_id, category, seed):
    return model.TermAssociation(
        term=model.Term(id=go_id, label=go_id),
        seed_source=seed,
        category=category,
    )


def test_precision_recall_f1():
    assocs = [
        _assoc("GO:1", "core_process", "enrichment_recovered"),       # good, recovered, core
        _assoc("GO:2", "supporting_process", "enrichment_recovered"), # good, recovered
        _assoc("GO:3", "false_association", "enrichment_recovered"),  # bad, recovered
        _assoc("GO:4", "nonspecific", "enrichment_recovered"),        # bad, recovered
        _assoc("GO:5", "core_component", "curator_added"),            # core, missed
        _assoc("GO:6", None, "enrichment_recovered"),                 # un-adjudicated -> ignored
    ]
    m = metrics.score(assocs)
    assert m["recovered_adjudicated"] == 4
    assert m["precision"] == 0.5      # 2 good of 4 recovered+adjudicated
    assert m["core_total"] == 2       # GO:1 (recovered) + GO:5 (missed)
    assert m["recall"] == 0.5         # 1 of 2 core recovered
    assert round(m["f1"], 3) == 0.5


def test_score_empty_is_none():
    m = metrics.score([])
    assert m["precision"] is None
    assert m["recall"] is None
    assert m["f1"] is None
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_metrics.py -v
```
Expected: FAIL (`metrics` module missing).

- [ ] **Step 3: Write the metrics module**

Create `python/genesets-workflows/src/genesets_workflows/curation/metrics.py`:

```python
"""Precision/recall/F1 of enrichment against a curated gold standard.

Precision: of the recovered terms that a curator adjudicated, how many are
biologically good (core/supporting)?
Recall: of the curator-marked core terms, how many did enrichment recover?
"""

from __future__ import annotations

from typing import Iterable

GOOD_CATEGORIES = {"core_process", "core_component", "supporting_process"}
CORE_CATEGORIES = {"core_process", "core_component"}
RECOVERED = "enrichment_recovered"


def precision(associations: Iterable) -> float | None:
    """Fraction of recovered, adjudicated terms that are good.

    >>> from types import SimpleNamespace as S
    >>> a = S(seed_source="enrichment_recovered", category="core_process")
    >>> b = S(seed_source="enrichment_recovered", category="false_association")
    >>> precision([a, b])
    0.5
    """
    recovered = [a for a in associations if a.seed_source == RECOVERED and a.category]
    if not recovered:
        return None
    good = [a for a in recovered if a.category in GOOD_CATEGORIES]
    return len(good) / len(recovered)


def recall(associations: Iterable) -> float | None:
    """Fraction of core terms that enrichment recovered.

    >>> from types import SimpleNamespace as S
    >>> a = S(seed_source="enrichment_recovered", category="core_process")
    >>> b = S(seed_source="curator_added", category="core_component")
    >>> recall([a, b])
    0.5
    """
    core = [a for a in associations if a.category in CORE_CATEGORIES]
    if not core:
        return None
    recovered = [a for a in core if a.seed_source == RECOVERED]
    return len(recovered) / len(core)


def f1(p: float | None, r: float | None) -> float | None:
    if p is None or r is None or (p + r) == 0:
        return None
    return 2 * p * r / (p + r)


def score(associations: Iterable) -> dict:
    associations = list(associations)
    p = precision(associations)
    r = recall(associations)
    recovered_adjudicated = sum(1 for a in associations if a.seed_source == RECOVERED and a.category)
    core_total = sum(1 for a in associations if a.category in CORE_CATEGORIES)
    return {
        "precision": p,
        "recall": r,
        "f1": f1(p, r),
        "recovered_adjudicated": recovered_adjudicated,
        "core_total": core_total,
    }
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_metrics.py -v
```
Expected: PASS (2 passed).

- [ ] **Step 5: Run the doctests to verify they pass**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest python -m pytest --doctest-modules python/genesets-workflows/src/genesets_workflows/curation/metrics.py -v
```
Expected: PASS (2 doctests passed).

- [ ] **Step 6: Commit**

```bash
git add python/genesets-workflows/src/genesets_workflows/curation/metrics.py \
        python/genesets-workflows/tests/test_metrics.py
git commit -m "feat(curation): add precision/recall/F1 metrics"
```

---

## Task 10: Report (join curation tree vs metrics, write TSV)

**Files:**
- Create: `python/genesets-workflows/src/genesets_workflows/curation/report.py`
- Create: `python/genesets-workflows/tests/test_report.py`

- [ ] **Step 1: Write the failing test**

Create `python/genesets-workflows/tests/test_report.py`:

```python
from pathlib import Path

from genesets_workflows.curation import report

ROOT = Path(__file__).resolve().parents[3]
C8_DIR = ROOT / "curation" / "c8"


def test_report_rows_for_example():
    rows = report.report_rows([C8_DIR / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"])
    assert len(rows) == 1
    row = rows[0]
    assert row["gene_set_id"] == "MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"
    # 5 recovered+adjudicated; 4 of them good (translation is nonspecific)
    assert row["recovered_adjudicated"] == 5
    assert row["precision"] == 0.8
    # core_total = antigen presentation (recovered) + adaptive immune (curator_added)
    assert row["core_total"] == 2
    assert row["recall"] == 0.5


def test_write_report_tsv(tmp_path):
    out = tmp_path / "report.tsv"
    report.write_report([C8_DIR / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"], out)
    header = out.read_text().splitlines()[0]
    assert header.split("\t")[0] == "gene_set_id"
    assert "precision" in header
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_report.py -v
```
Expected: FAIL (`report` module missing).

- [ ] **Step 3: Write the report module**

Create `python/genesets-workflows/src/genesets_workflows/curation/report.py`:

```python
"""Build precision/recall/F1 report rows from curated interpretation files."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from genesets_workflows.curation import metrics, model

COLUMNS = [
    "gene_set_id",
    "gene_set_name",
    "collection",
    "curation_status",
    "n_associations",
    "recovered_adjudicated",
    "core_total",
    "precision",
    "recall",
    "f1",
]


def report_rows(paths: Iterable[Path]) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        interp = model.load_interpretation(path)
        scored = metrics.score(interp.associations)
        rows.append(
            {
                "gene_set_id": interp.gene_set_id,
                "gene_set_name": interp.gene_set_name,
                "collection": interp.collection or "",
                "curation_status": interp.curation_status or "",
                "n_associations": len(interp.associations),
                **scored,
            }
        )
    return rows


def write_report(paths: Iterable[Path], out: Path) -> Path:
    rows = report_rows(paths)
    with Path(out).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return Path(out)


def discover(curation_dir: Path) -> list[Path]:
    return sorted(Path(curation_dir).glob("*.yaml"))
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_report.py -v
```
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add python/genesets-workflows/src/genesets_workflows/curation/report.py \
        python/genesets-workflows/tests/test_report.py
git commit -m "feat(curation): add precision/recall report builder"
```

---

## Task 11: Typer CLI + dispatcher wiring

**Files:**
- Create: `python/genesets-workflows/src/genesets_workflows/curation/cli.py`
- Modify: `python/genesets-workflows/src/genesets_workflows/cli.py`
- Create: `python/genesets-workflows/tests/test_curation_cli.py`

- [ ] **Step 1: Write the failing CLI test**

Create `python/genesets-workflows/tests/test_curation_cli.py`:

```python
from pathlib import Path

from typer.testing import CliRunner

from genesets_workflows.curation import cli as curation_cli
from genesets_workflows.cli import dispatch

ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "curation" / "c8" / "HAY_BONE_MARROW_DENDRITIC_CELL.yaml"

runner = CliRunner()


def test_report_command_writes_tsv(tmp_path):
    out = tmp_path / "r.tsv"
    result = runner.invoke(curation_cli.app, ["report", "--input", str(EXAMPLE), "--output", str(out)])
    assert result.exit_code == 0, result.output
    assert out.exists()
    assert "precision" in out.read_text().splitlines()[0]


def test_help_lists_subcommands():
    result = runner.invoke(curation_cli.app, ["--help"])
    assert result.exit_code == 0
    for name in ("draft", "validate", "report"):
        assert name in result.output


def test_dispatcher_routes_curate():
    assert dispatch("curate") is curation_cli.main
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_curation_cli.py -v
```
Expected: FAIL (`curation.cli` missing; dispatcher has no "curate").

- [ ] **Step 3: Write the Typer CLI**

Create `python/genesets-workflows/src/genesets_workflows/curation/cli.py`:

```python
"""Typer CLI for the `genesets-workflows curate` command group."""

from __future__ import annotations

from pathlib import Path

import typer

from genesets_workflows.curation import report as report_mod
from genesets_workflows.curation import validate as validate_mod

app = typer.Typer(no_args_is_help=True, add_completion=False, help="Curate GO interpretations of gene sets.")


@app.command()
def report(
    input: list[Path] = typer.Option(None, "--input", "-i", help="Interpretation YAML file(s)."),
    curation_dir: Path = typer.Option(None, "--dir", help="Directory of interpretation YAMLs."),
    output: Path = typer.Option(..., "--output", "-o", help="Output TSV path."),
) -> None:
    """Compute precision/recall/F1 for curated interpretations."""
    paths = list(input or [])
    if curation_dir:
        paths.extend(report_mod.discover(curation_dir))
    if not paths:
        raise typer.BadParameter("provide --input and/or --dir")
    report_mod.write_report(paths, output)
    typer.echo(f"wrote {output}")


@app.command()
def validate(
    path: Path = typer.Argument(..., help="Interpretation YAML file to validate."),
    schema: Path = typer.Option(Path("curation/schema/genesets_interpretation.yaml"), "--schema"),
    oak_config: Path = typer.Option(Path("curation/conf/oak_config.yaml"), "--oak-config"),
    cache_dir: Path = typer.Option(Path("curation/references_cache"), "--cache-dir"),
    skip_references: bool = typer.Option(False, "--skip-references", help="Skip the reference-validator step."),
) -> None:
    """Run structural, term, and reference validation on an interpretation file."""
    code = validate_mod.validate_file(
        path, schema=schema, oak_config=oak_config, cache_dir=cache_dir, skip_references=skip_references
    )
    raise typer.Exit(code)


@app.command()
def draft(
    gene_set_id: str = typer.Argument(..., help="e.g. MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL"),
    enrichment_tsv: Path = typer.Option(..., "--enrichment-tsv", help="genesets-rs enrich TSV for this set."),
    output: Path = typer.Option(..., "--output", "-o"),
    collection: str = typer.Option("C8", "--collection"),
    go_version: str = typer.Option("unspecified", "--go-version"),
) -> None:
    """Build a draft interpretation from a precomputed enrichment TSV."""
    from genesets_workflows.curation import draft as draft_mod
    from genesets_workflows.curation import model
    from genesets_workflows.curation.enrichment import parse_enrichment_tsv

    rows = parse_enrichment_tsv(Path(enrichment_tsv).read_text())
    name = gene_set_id.split(":", 1)[-1]
    interp = draft_mod.assemble_draft(
        gene_set_id=gene_set_id,
        gene_set_name=name,
        collection=collection,
        rows=rows,
        provenance=model.EnrichmentProvenance(tool="genesets-rs", go_version=go_version),
    )
    model.dump_interpretation(interp, output)
    typer.echo(f"wrote {output} ({len(interp.associations)} associations)")


def main(argv: list[str] | None = None) -> int:
    """argv bridge so the existing argparse dispatcher can invoke this Typer app."""
    from typer.main import get_command

    command = get_command(app)
    try:
        return command.main(
            args=list(argv) if argv is not None else None,
            prog_name="genesets-workflows curate",
            standalone_mode=False,
        ) or 0
    except SystemExit as exc:
        return int(exc.code or 0)
```

- [ ] **Step 4: Wire `curate` into the dispatcher**

In `python/genesets-workflows/src/genesets_workflows/cli.py`:

Add to the imports block (after the existing `from genesets_workflows.sources import ...` lines):

```python
from genesets_workflows.curation import cli as curation_cli
```

In `parse_args`, after the `explore` subparser block, add:

```python
    subparsers.add_parser(
        "curate",
        help="Curate and validate GO interpretations of gene sets.",
        add_help=False,
    )
```

In `dispatch`, add `"curate": curation_cli.main,` to the `commands` dict:

```python
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
```

- [ ] **Step 5: Run the test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_curation_cli.py -v
```
Expected: PASS (3 passed). Note: `validate` imports `validate.py` (Task 12) — if running this task before Task 12, the `validate` subcommand import will fail. Implement Task 12's `validate.py` **before** running Step 5, or temporarily run only `test_report_command_writes_tsv`. Recommended order: do Task 12 Step 3 (write `validate.py`) first, then this step.

- [ ] **Step 6: Smoke-test the integrated CLI end to end**

Run:
```bash
uv run --project python/genesets-workflows --extra curation genesets-workflows curate --help
uv run --project python/genesets-workflows --extra curation genesets-workflows curate report -i curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml -o /tmp/curate_report.tsv
```
Expected: help lists draft/validate/report; report writes `/tmp/curate_report.tsv` and echoes the path.

- [ ] **Step 7: Commit**

```bash
git add python/genesets-workflows/src/genesets_workflows/curation/cli.py \
        python/genesets-workflows/src/genesets_workflows/cli.py \
        python/genesets-workflows/tests/test_curation_cli.py
git commit -m "feat(curation): add Typer curate CLI and wire into dispatcher"
```

---

## Task 12: Validation command builders + add evidence to the example

**Files:**
- Create: `python/genesets-workflows/src/genesets_workflows/curation/validate.py`
- Create: `python/genesets-workflows/tests/test_validate.py`
- Modify: `curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml`

- [ ] **Step 1: Write the failing command-builder test**

Create `python/genesets-workflows/tests/test_validate.py`:

```python
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_validate.py -v
```
Expected: FAIL (`validate` module missing).

- [ ] **Step 3: Write the validate module**

Create `python/genesets-workflows/src/genesets_workflows/curation/validate.py`:

```python
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests/test_validate.py -v
```
Expected: PASS (2 passed).

- [ ] **Step 5: Source a real snippet for the example's evidence**

Look up a real reference for dendritic-cell antigen presentation and read its cached text (this populates `curation/references_cache/`):
```bash
uv run --project python/genesets-workflows --extra curation \
  linkml-reference-validator lookup PMID:9521319
```
This prints the title/abstract of Banchereau & Steinman (1998), "Dendritic cells and the control of immunity." Copy one **exact, contiguous** sentence or phrase from the printed abstract — you will paste it verbatim as `snippet`. If `PMID:9521319` has no usable abstract text, pick another DC review from the printed candidates and use a phrase from it.

- [ ] **Step 6: Add the evidence item to the example**

In `curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml`, add an `evidence` block to the `antigen processing and presentation` association (replace `<EXACT PHRASE FROM STEP 5>` with the verbatim text you copied; do not paraphrase):

```yaml
  - term:
      id: GO:0019882
      label: antigen processing and presentation
    aspect: biological_process
    category: core_process
    confidence: high
    specificity: cell_type_specific
    seed_source: enrichment_recovered
    evidence:
      - reference: PMID:9521319
        reference_title: Dendritic cells and the control of immunity
        supports: SUPPORT
        snippet: <EXACT PHRASE FROM STEP 5>
        explanation: >-
          Dendritic cells are professional antigen-presenting cells, so antigen
          processing and presentation is a core process for this signature.
```

- [ ] **Step 7: Run all three validators on the example**

Run:
```bash
uv run --project python/genesets-workflows --extra curation \
  genesets-workflows curate validate curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml
```
Expected: all three validators pass. If the reference-validator reports the snippet was not found, the text was not verbatim — re-copy the exact phrase from the cached paper (`curation/references_cache/`) and retry.

- [ ] **Step 8: Run the whole test suite**

Run:
```bash
uv run --project python/genesets-workflows --extra curation --with pytest pytest python/genesets-workflows/tests -v
```
Expected: all tests pass.

- [ ] **Step 9: Commit**

```bash
git add python/genesets-workflows/src/genesets_workflows/curation/validate.py \
        python/genesets-workflows/tests/test_validate.py \
        curation/c8/HAY_BONE_MARROW_DENDRITIC_CELL.yaml \
        curation/cache
git commit -m "feat(curation): add validator runner and literature evidence to example"
```

---

## Task 13: just targets, docs, and the C8 manifest

**Files:**
- Modify: `justfile`
- Create: `curation/README.md`
- Create: `curation/c8/manifest.tsv`
- Create: `docs/src/curation.md`
- Modify: `docs/src/SUMMARY.md`
- Modify: `README.md`

- [ ] **Step 1: Add just targets**

Append to `justfile`:

```just
# --- Curation ---------------------------------------------------------------

curation_schema := "curation/schema/genesets_interpretation.yaml"
curation_oak := "curation/conf/oak_config.yaml"
gw := "uv run --project python/genesets-workflows --extra curation"

# Run the curation test suite (pytest + doctests).
curate-test:
    {{gw}} --with pytest pytest python/genesets-workflows/tests -v
    {{gw}} --with pytest python -m pytest --doctest-modules \
      python/genesets-workflows/src/genesets_workflows/curation -v

# Validate every interpretation YAML under curation/.
curate-validate:
    for f in curation/**/*.yaml; do \
      echo "== $f =="; \
      {{gw}} genesets-workflows curate validate "$f" || exit 1; \
    done

# Validate the schema's own ontology-term meanings.
curate-validate-schema:
    {{gw}} linkml-term-validator validate-schema {{curation_schema}} -c {{curation_oak}}

# Build the precision/recall report over all curated C8 interpretations.
curate-report out="curation/report.tsv":
    {{gw}} genesets-workflows curate report --dir curation/c8 -o {{out}}
```

- [ ] **Step 2: Verify the targets are listed and the report runs**

Run:
```bash
just --list | grep curate
just curate-report /tmp/curate_report.tsv
```
Expected: the four `curate-*` targets appear; `curate-report` writes `/tmp/curate_report.tsv`.

- [ ] **Step 3: Write the C8 manifest**

Create `curation/c8/manifest.tsv` (seed list; columns: `gene_set_id`, `context_term`, `context_type`, `context_label`). Verify exact set IDs against MSigDB before drafting each.

```
gene_set_id	context_term	context_type	context_label
MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL	CL:0000451	cell_type	dendritic cell
MSIGDB:HAY_BONE_MARROW_NK_CELLS	CL:0000623	cell_type	natural killer cell
MSIGDB:TRAVAGLINI_LUNG_ALVEOLAR_EPITHELIAL_TYPE_2_CELL	CL:0002063	cell_type	type II pneumocyte
MSIGDB:DESCARTES_FETAL_LIVER_HEPATOBLASTS	CL:0000182	cell_type	hepatocyte
MSIGDB:AIZARANI_LIVER_C9_KUPFFER_CELLS	CL:0000091	cell_type	Kupffer cell
```

- [ ] **Step 4: Write the curation README**

Create `curation/README.md`:

```markdown
# Curated GO Interpretations of Non-GO Gene Sets

Curated, validated GO-term interpretations of non-GO MSigDB gene sets
(C8 cell-type signatures first). Each `c8/<SET>.yaml` records, per GO term, a
curator judgment (role category, confidence, specificity), optional literature
evidence, and the enrichment provenance that seeded it. The data doubles as a
precision/recall gold standard for `genesets-rs`.

## Layout
- `schema/genesets_interpretation.yaml` — LinkML schema (source of truth).
- `conf/` — OAK adapters and reference-validator config.
- `c8/manifest.tsv` — seed list of C8 sets and their biological context terms.
- `c8/<SET>.yaml` — one interpretation per gene set.
- `cache/` — committed ontology label cache (written by linkml-term-validator).
- `references_cache/` — cited-paper cache (gitignored; rebuilt on demand).

## Workflow
1. `just curate-validate-schema` — sanity-check the schema's enum meanings.
2. Draft from enrichment output:
   `genesets-workflows curate draft MSIGDB:<SET> --enrichment-tsv <tsv> -o c8/<SET>.yaml`
3. Adjudicate: set each association's `category`, `confidence`, `specificity`;
   add `curator_added` core terms the tool missed; add `evidence`.
4. `just curate-validate` — structural + term + reference validation.
5. `just curate-report` — precision/recall/F1.

## Categories
`core_process`, `core_component`, `supporting_process`,
`marker_driven_plausible`, `nonspecific`, `false_association`.
```

- [ ] **Step 5: Write the docs book page**

Create `docs/src/curation.md`:

```markdown
# Curated GO Interpretations

`genesets-rs` curates "correct/reasonable" GO interpretations for non-GO gene
sets (MSigDB C8 cell-type signatures first) as a precision/recall gold standard.

Interpretations live in `curation/c8/*.yaml`, conform to the LinkML schema at
`curation/schema/genesets_interpretation.yaml`, and are validated by
`linkml-validate`, `linkml-term-validator` (GO/CL/UBERON IDs + labels), and
`linkml-reference-validator` (literature snippets).

Each GO term carries a curator role category (`core_process`, `core_component`,
`supporting_process`, `marker_driven_plausible`, `nonspecific`,
`false_association`), a confidence and a specificity, the seeding enrichment
stats, and optional cited evidence. See `curation/README.md` for the workflow.

Build the report:

```bash
just curate-report
```
```

- [ ] **Step 6: Link the page into the book**

In `docs/src/SUMMARY.md`, add a line under the appropriate section (after the `evals.md` entry):

```markdown
- [Curated GO Interpretations](./curation.md)
```

Confirm the entry references an existing file:
```bash
grep -n "curation.md" docs/src/SUMMARY.md
```
Expected: the new line is present.

- [ ] **Step 7: Add a README section**

In the top-level `README.md`, after the `## Evals` section, add:

```markdown
## Curated GO Interpretations

`curation/` holds expert GO-term interpretations of non-GO gene sets (MSigDB C8
cell-type signatures first), validated with LinkML term/reference validators and
used as a precision/recall gold standard. See `curation/README.md`.

```bash
just curate-validate
just curate-report
```
```

- [ ] **Step 8: Final full verification**

Run:
```bash
just curate-test
just curate-validate
just curate-report /tmp/curate_report.tsv
```
Expected: tests pass; the example validates against all three validators; the report writes with `precision`/`recall`/`f1` columns.

- [ ] **Step 9: Commit**

```bash
git add justfile curation/README.md curation/c8/manifest.tsv \
        docs/src/curation.md docs/src/SUMMARY.md README.md
git commit -m "feat(curation): add just targets, docs, and C8 manifest"
```

---

## Self-Review

**Spec coverage:**
- §4 layout → Tasks 1–4, 13. §5 schema (all classes/enums) → Task 2; `BiologicalContext` ontology-agnostic context → Task 2 + example Task 4. §6 category enum → Task 2 + metrics Task 9. §7 enrichment-seeded drafting → Tasks 6–8, 11 (`curate draft`); the optional LLM pre-fill (§7 step 2) is intentionally **not** automated here — the schema/validators support it and the README workflow leaves adjudication to a curator; add as a follow-up. §8 validation (all three validators) → Tasks 4, 12. §9 reporting (precision/recall/F1) → Tasks 9–10, 13 (`curate report`); stratification by GO version/evidence code reuses existing eval bundles and is a reporting follow-up once interpretations exist. §10 testing conventions → tests in every task + doctests (Tasks 7, 9) + `just curate-test`. §11 phasing → task order.
- Gap noted intentionally: GO-version/evidence-code **stratified** report and LLM pre-fill are deferred (need a corpus of curated sets first); the schema + provenance fields already capture what they need.

**Placeholder scan:** No "TBD"/"TODO" in code. The two free-text strings (`go_version: pinned-snapshot-TBD` in the hand-authored example, `<EXACT PHRASE FROM STEP 5>`) are explicitly called out as human-fill provenance/snippet values, with concrete commands to source them — not schema placeholders.

**Type consistency:** `assemble_draft` (Task 8) returns `model.GeneSetInterpretation`; `model.dump_interpretation`/`load_interpretation` used consistently (Tasks 5, 8, 10, 11). `EnrichmentRow` fields produced in Task 6 are consumed unchanged in Task 8. `metrics.score` keys (`precision`, `recall`, `f1`, `recovered_adjudicated`, `core_total`) produced in Task 9 are consumed in Task 10 `report_rows` and asserted in Task 10 tests. `validate.build_commands`/`validate_file` signature (Task 12) matches the call in `cli.validate` (Task 11). Category strings identical across schema enum (Task 2), example (Task 4), and `GOOD_/CORE_CATEGORIES` (Task 9). Note flagged in Task 11 Step 5: implement Task 12's `validate.py` before running the Task 11 CLI test, since `cli.py` imports it.
