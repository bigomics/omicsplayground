## EXPERIMENT CONTEXT

This analysis compares differential expression signatures across multiple contrasts to identify shared and unique biological signals.

Experiment: {experiment}

## INTERSECTION-SPECIFIC GUIDANCE

- Contrast intersection analysis identifies genes that are differentially expressed across multiple experimental comparisons
- Venn diagrams show the overlap of significant genes (passing FDR and logFC thresholds) between selected contrasts
- Genes appearing in the intersection of all selected contrasts represent a shared biological response across conditions
- Genes unique to a single contrast represent condition-specific responses
- Fold-change heatmaps display the top differentially expressed genes across all selected contrasts, revealing concordant and discordant regulation patterns
- Pairwise correlation of fold-change profiles quantifies how similar two contrasts are at the genome-wide level

## QUANTITATIVE INTERPRETATION GUIDELINES

Venn diagram metrics:
- Overlap count: number of genes passing significance thresholds in multiple contrasts simultaneously
- Large overlaps suggest a shared biological mechanism or common upstream regulator
- Small or empty overlaps suggest distinct, condition-specific responses

Fold-change correlation:
- Pearson r > 0.7: strong similarity between contrast profiles; likely shared biology
- Pearson r 0.4-0.7: moderate similarity; partially overlapping mechanisms
- Pearson r < 0.4: weak similarity; largely distinct responses
- Negative correlation: opposing regulation patterns between contrasts

Heatmap interpretation:
- Concordant rows (same direction across contrasts) indicate conserved regulation
- Discordant rows (mixed direction) indicate condition-specific regulation
- Genes with large absolute fold change in multiple contrasts are high-confidence shared drivers

Integration rules:
- Prioritize genes that are significant and concordantly regulated across multiple contrasts
- Genes with discordant regulation (up in one contrast, down in another) may reveal context-dependent biology
- Consider the biological relationship between the contrasts when interpreting overlaps
- Large shared gene sets between related conditions strengthen mechanistic interpretations
- Treat genes unique to a single contrast as candidates for condition-specific follow-up
