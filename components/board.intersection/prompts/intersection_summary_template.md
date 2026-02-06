## Task

Interpret the biological significance of the contrast intersection analysis for the selected contrasts based on the data provided below.

## Analysis Overview

**Selected contrasts:** {contrasts}

**Analysis level:** {level}

This intersection analysis compares differential expression profiles across the selected contrasts to identify shared and unique biological signals.

## Venn Diagram Summary

The following summarizes the overlap of significantly differentially expressed features across the selected contrasts:

{venn_summary}

Focus on identifying:
1. Whether the contrasts share a substantial common signature or are largely distinct
2. The biological meaning of shared versus unique genes
3. Whether the overlap size is expected given the experimental design

## Fold-Change Heatmap: Top Genes

The following genes show the strongest differential expression across the selected contrasts:

{heatmap_genes}

## Contrast Correlation

Pairwise Pearson correlation of fold-change profiles between contrasts:

{correlation_info}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**Intersection analysis** | {n_contrasts} contrasts compared

Then synthesize the above information to explain:
- The degree of similarity or divergence between the selected contrasts
- Key biological processes represented by the shared genes
- Condition-specific signals revealed by unique gene sets
- How the pairwise correlations support or refine the overlap interpretation
- Notable genes that are concordantly or discordantly regulated across contrasts
- Use quantitative metrics to support your interpretation (cite specific overlap counts, correlation values, and fold changes when available)
