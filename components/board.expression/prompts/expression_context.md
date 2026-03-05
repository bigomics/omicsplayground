## EXPERIMENT CONTEXT

This analysis uses differential gene expression analysis to identify genes with significant changes between experimental conditions.

Experiment: {{experiment}}

## DIFFERENTIAL EXPRESSION-SPECIFIC GUIDANCE

- Differential expression analysis compares gene expression levels between two conditions (e.g., treatment vs. control) to identify significantly upregulated or downregulated genes
- Multiple statistical methods are used (limma-trend, edgeR QLF/LRT, DESeq2 Wald/LRT) and results are combined into a meta q-value
- The meta q-value corresponds to the maximum q-value across methods, ensuring conservative significance calls
- The number of stars indicates how many methods detected a significant differential expression
- logFC (log2 fold change) represents the magnitude and direction of expression change between conditions
- Volcano plots visualize fold change (x-axis) against statistical significance (y-axis) for all tested genes

## QUANTITATIVE INTERPRETATION GUIDELINES

Fold change metrics:
- |logFC| > 2: strong differential expression; 1-2 moderate; 0.5-1 mild; < 0.5 minimal
- Positive logFC: upregulated in the first group relative to the second
- Negative logFC: downregulated in the first group relative to the second

Significance metrics:
- meta.q < 0.01: high confidence; 0.01-0.05 standard significance; 0.05-0.10 marginal; > 0.10 not significant
- Stars: more stars indicate higher cross-method agreement and reliability
- FDR (false discovery rate) corrects for multiple testing across thousands of genes

Integration rules:
- Prioritize genes with both strong logFC and low meta.q values
- Look for functional coherence among top DE genes (e.g., co-regulated pathway members, known interactors)
- Consider the biological plausibility of expression changes in the context of the experimental comparison
- Genes with concordant results across multiple statistical methods are most reliable
- Treat marginal q-values or minimal fold changes as tentative and avoid over-interpretation
- Very large fold changes with high q-values may reflect noisy or lowly expressed genes
