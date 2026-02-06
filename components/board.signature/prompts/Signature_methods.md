## EXPERIMENT CONTEXT

This analysis tests a gene signature against known gene set databases and contrast profiles to identify biological overlap and functional enrichment.

Experiment: {experiment}

## GENE SIGNATURE ANALYSIS GUIDANCE

- A gene signature is a defined set of genes hypothesized to represent a biological process, disease state, or experimental condition
- Overlap analysis compares the signature against curated gene set databases using Fisher's exact test
- GSEA enrichment tests whether the signature genes are enriched at the top or bottom of ranked contrast profiles
- High overlap with known pathways suggests the signature captures recognized biological processes
- Novel signatures with low overlap may represent previously uncharacterized biology

## QUANTITATIVE INTERPRETATION GUIDELINES

Overlap metrics:
- Score combines odds ratio and q-value significance; higher is better
- Odds ratio > 10: very strong overlap; 5-10 strong; 2-5 moderate; < 2 weak
- Q-value (Fisher) < 0.01: high confidence; 0.01-0.05 standard; 0.05-0.10 marginal; > 0.10 not significant
- k/K ratio: higher ratio with smaller K indicates more specific overlap

Gene-level metrics:
- Log2 fold change |LogFC| > 2: strong differential expression; 1-2 moderate; 0.5-1 mild; < 0.5 minimal
- Q-value < 0.01: highly significant; 0.01-0.05 significant; 0.05-0.10 marginal

Integration rules:
- Prioritize pathways with both high scores and low q-values
- Cross-reference overlapping genes with their fold changes to assess direction of regulation
- Multiple related pathways enriched together strengthen biological interpretation
- Treat marginal q-values or low overlap ratios as tentative and avoid over-interpretation
