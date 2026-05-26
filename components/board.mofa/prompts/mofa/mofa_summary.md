## EXPERIMENT CONTEXT

This analysis uses Multi-Omics Factor Analysis (MOFA) to decompose a
multi-omics dataset into latent factors capturing shared or view-specific
variation.

Experiment: {{experiment}}

## Task

Interpret the biological function and significance of MOFA factor
**{{factor}}** based on the per-view loadings, pathway enrichment, and
trait coordination provided below.

## Factor data

{{factor_detail}}

## Output Instructions

Always begin your summary with this heading: **Factor {{factor}}**

Then synthesise the data block above to explain:
- The dominant biological programme this factor captures
- Which omics views contribute and whether the factor is cross-view
- How the enriched pathways collectively support that interpretation
- How the top-weighted features (especially cross-view ones) anchor the
  programme
- The biological relevance to the experimental traits

Treat every verbal label in the data block (`dominant` / `strong` /
`major` / `hub` / `correlated` / `highly significant` / `up-strong` …)
as authoritative — do not re-derive from raw numbers; raw numbers are
omitted by design.
