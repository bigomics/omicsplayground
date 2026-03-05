EXPERIMENT: {{experiment}}
ANALYSIS TYPE: {{analysis_type}}
ANALYSIS TYPE DESCRIPTION: {{analysis_type_description}}

## Contrast Signal Classification
(Signal strength ranking; score = weighted combination of significant drug count, max|NES|, MOA convergence, target enrichment; use tiers to decide which contrasts get individual subsections vs. Minor Signals)

{{rank_table}}

## Cross-Contrast MOA Convergence
(NES per pharmacological class across all contrasts; ** q<0.05, * p<0.05, — = not in top enriched results for that contrast; negative = opposes signature, positive = mimics)

{{moa_matrix}}

## Contrast Detail
(Each contrast block starts with a deterministic interpretation evidence summary. Treat that summary as the primary interpretation anchor. Raw drug, MOA, and target tables that follow are supporting evidence only. Negative NES = drug opposes the experimental signature; positive NES = drug mimics it; use as evidence for the narrative, do not reproduce as tables)

{{contrast_detail}}
