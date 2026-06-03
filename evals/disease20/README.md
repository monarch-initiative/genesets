# Disease20 Eval

Disease20 is a starter eval manifest of 20 human Disease Ontology gene sets from MyGeneset.info.

The manifest is hand-selected for breadth:

- neurodegenerative: Alzheimer's disease, Parkinson's disease, amyotrophic lateral sclerosis, Huntington's disease;
- psychiatric/neurodevelopmental: schizophrenia, bipolar disorder, autism spectrum disorder;
- metabolic/cardiovascular: type 2 diabetes mellitus, obesity, coronary artery disease, myocardial infarction;
- immune/inflammatory/respiratory: asthma, rheumatoid arthritis, systemic lupus erythematosus, Crohn's disease, ulcerative colitis;
- cancer: breast, colorectal, prostate, and lung cancer.

Fetch a local GMT snapshot:

```bash
python3 scripts/fetch_mygeneset_eval.py \
  --manifest evals/disease20/sets.tsv \
  --out-dir evals/disease20/generated
```

Run disease-vs-disease overlap enrichment:

```bash
genesets-rs run evals/disease20/config.yaml
```

The generated files are ignored by git. Commit a generated snapshot only when we intentionally want a fully frozen benchmark fixture.
