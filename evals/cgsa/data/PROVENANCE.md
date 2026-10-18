# Provenance

The TSV/GMT files in this directory are a faithful re-serialization of data
released by the cGSA authors. They are **not** authored by this repository.

- **Source:** https://github.com/ncbi-nlp/cGSA (`Supplementary/` workbooks).
- **Paper:** Wang Z, Day C-P, Wei C-H, et al. *Knowledge-guided contextual gene
  set analysis with large language models.* Bioinformatics 2026;42:btag214.
  https://doi.org/10.1093/bioinformatics/btag214
- **License:** the paper is Open Access (CC BY 4.0); the cGSA repository is
  public. Consult the upstream repository for the code/data license terms before
  redistribution.
- **How built:** `python evals/cgsa/data/build_cgsa_data.py` downloads
  `Curated DEGs.xlsx` and `Results of cGSA.xlsx` from
  `raw.githubusercontent.com/ncbi-nlp/cGSA/main/Supplementary/` (cached under
  `data/_workbooks/`, gitignored) and flattens the sheets. Re-run to refresh.

The source `.xlsx` are intentionally **not** committed (they are binary and
re-fetchable); only the derived text artifacts are tracked.
