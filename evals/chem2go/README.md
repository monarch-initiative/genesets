# chem2go smoke checks

Small, unambiguous metabolite sets that each bridge of `prepare-chem2go` should
resolve to a known pathway. They are smoke checks, not a benchmark: the
background is a proxy (annotated chemicals with an HMDB cross-reference,
~12.4k) rather than a platform's measured metabolites, so p-values are only
illustrative. See `docs/src/metabolite-enrichment.md` for the method.

```bash
just chem2go-prep     # writes generated/current (gitignored)
cargo build --release
just chem2go-smoke    # writes generated/smoke
```

Queries (`smoke/`), given on purpose as mixed id schemes and neutral forms:

- `tca.txt`: citrate, isocitrate, 2-oxoglutarate, succinate, fumarate, malate,
  oxaloacetate, succinyl-CoA.
- `urea.txt`: ornithine, citrulline, argininosuccinate, arginine, urea,
  carbamoyl phosphate, aspartate, glutamate.
- `extra_xrefs.tsv`: two HMDB → ChEBI fixes ChEBI itself lacks.

## Results (GO-plus 2026-07-26, ChEBI 256, human GAF, 2026-10-18)

Top terms by adjusted p (overlap / term size):

| Query | Bridge | Top terms |
|---|---|---|
| TCA | logical_definition | C4-dicarboxylate transport (6/17), dicarboxylic acid transport (6/100) |
| TCA | rhea | C4-dicarboxylate transport (5/10), antiporter activity (5/23) |
| TCA | enzyme | malate metabolic process (5/9), **tricarboxylic acid cycle** (6/34) |
| TCA | all | C4-dicarboxylate transport (6/17), tricarboxylic acid transport (6/32), **tricarboxylic acid cycle** (6/34) |
| Urea | logical_definition | linear-amidine C–N hydrolase activity (6/23), arginase activity (3/3) |
| Urea | rhea | linear-amidine C–N hydrolase activity (6/25), arginase activity (3/3) |
| Urea | enzyme | arginine metabolic process (6/14), **urea cycle** (6/21) |
| Urea | all | arginine metabolic process (6/16), **urea cycle** (6/21) |

The species-independent bridges recover the chemistry (transporters,
catalytic activities). Only the enzyme bridge reaches the pathway processes,
because GO's own axioms rarely link those to chemicals.

Prep stats for that run (`metadata.json`): 219,011 ChEBI classes collapsed to
208,730; 49,190 chemicals in the default universe; 32,501 annotated in `all`
(1.99M rows over 13,858 GO terms), of which `logical_definition` 339k rows,
`rhea` 162k, `enzyme` 1.67M.
