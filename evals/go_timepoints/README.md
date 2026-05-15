# GO Timepoints Eval

This eval runs the Disease20 query set against official GO archive releases at three timepoints:

- `2026-03-25`: current release used for the eval;
- `2021-05-01`: roughly five years before the May 2026 analysis date;
- `2016-05-01`: roughly ten years before the May 2026 analysis date.

Each timepoint uses official release URLs:

```text
https://release.geneontology.org/<date>/ontology/go-basic.obo
https://release.geneontology.org/<date>/annotations/goa_human.gaf.gz
```

Run:

```bash
python3 scripts/run_go_timepoints_eval.py
```

The runner prepares each GO release, runs the `all`, `no_contributes_to`, `iba`, and `iba_iea` variants, and writes:

- `generated/<timepoint>/metadata.yaml`: prep metadata for that timepoint;
- `generated/<timepoint>/run_metadata.yaml`: run metadata for that timepoint;
- `generated/<timepoint>/<variant>/results.tsv`: significant result rows;
- `generated/summary.yaml`: cross-timepoint counts and pairwise diff summaries.

The default result filter is `max_p_adjust: 0.05`.
