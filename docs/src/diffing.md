# Diffing

Diffing should operate on canonical, unpostprocessed or minimally filtered result tables. Post-processing can make a term disappear for reasons unrelated to enrichment change, which is useful for presentation but bad for version comparison.

## Use Cases

Common comparisons:

- GO release A vs GO release B;
- `contributes_to` retained vs filtered;
- IBA-only vs all evidence;
- IBA plus IEA vs all evidence;
- alternative backgrounds;
- alternative correction methods.

## Stable Join Keys

The default diff key should be:

```text
query_id, target_id
```

For ontology version comparisons, classify missing terms explicitly:

- `present_both`;
- `missing_left`;
- `missing_right`;
- `obsolete_or_removed`;
- `merged_or_replaced`, when ontology metadata can prove it;
- `filtered_by_p_cutoff`, when a term exists in the full run but not in the cutoff TSV.

This means serious diffs need either full result tables or at least a top-k/threshold run that is generous enough not to hide important deltas.

## Delta Metrics

Useful columns:

| Column | Meaning |
| --- | --- |
| `delta_log10_p_adjust` | `-log10(p_adjust_right) - -log10(p_adjust_left)` |
| `delta_log10_p_value` | raw p-value version of the above |
| `delta_overlap` | overlap count difference |
| `delta_target_size` | target annotation size difference |
| `delta_query_size` | query size difference, usually should be zero |
| `left_rank`, `right_rank` | within-query rank by adjusted p-value |
| `delta_rank` | rank movement |

Rank deltas are often more interpretable than tiny p-value deltas at extreme significance.

For thresholded result TSVs, the primary endpoint should be significance crossing:

| Class | Meaning |
| --- | --- |
| `shared_significant` | pair is significant in both versions |
| `lost_significant` | pair is significant on the left but not on the right |
| `gained_significant` | pair is significant on the right but not on the left |

A term moving from `1e-10` to `1e-20` is usually less important than a term crossing the cutoff. Large log-p deltas are still useful as a secondary diagnostic, especially when target size or background size changed dramatically.

## Ordering

Recommended pipeline:

1. Prepare each ontology/annotation variant.
2. Run enrichment.
3. Diff canonical results.
4. Optionally post-process each side for display.
5. Optionally diff reduced display tables, but label that as a presentation diff.

This ordering keeps version comparisons honest. A term that disappears due to redundancy pruning should not be confused with a term that disappeared because the ontology or annotations changed.

## Metadata

Every diff output should include companion metadata:

- left/right result file paths and digests;
- left/right ontology and annotation metadata paths;
- cutoff policy;
- whether missing rows were looked up in a full table;
- diff metrics and ranking policy;
- post-processing status.
