## EXPERIMENT CONTEXT

This analysis uses gene set enrichment analysis to identify differentially enriched pathways and biological processes.

Experiment: {{experiment}}

## ENRICHMENT-SPECIFIC GUIDANCE

- Gene set enrichment analysis tests whether predefined sets of genes show statistically significant differences between experimental conditions
- Multiple methods are used (CAMERA, GSEA, fGSEA, GSVA, ssGSEA, Fisher) and results are combined into a meta q-value
- The meta q-value corresponds to the maximum q-value across methods, ensuring conservative significance calls
- The number of stars indicates how many methods detected a significant enrichment
- logFC represents the average fold change of genes within each gene set

## QUANTITATIVE INTERPRETATION GUIDELINES

Enrichment metrics:
- logFC > 1: strong upregulation of the pathway; 0.5-1 moderate; 0.2-0.5 mild; < 0.2 minimal
- logFC < -1: strong downregulation; -0.5 to -1 moderate; -0.2 to -0.5 mild; > -0.2 minimal
- meta.q < 0.01: high confidence; 0.01-0.05 standard significance; 0.05-0.10 marginal; > 0.10 not significant
- Stars: more stars indicate higher cross-method agreement and reliability

Gene-level metrics within enriched gene sets:
- Gene fold change (fc) reflects individual gene contribution to the pathway signal
- Genes with high absolute fold change are likely drivers of the enrichment signal
- Concordant direction of gene-level changes strengthens the pathway interpretation

Integration rules:
- Prioritize pathways with both strong logFC and low meta.q values
- Look for thematic coherence among top enriched pathways (e.g., multiple immune-related or metabolic pathways)
- Consider gene set size: very large sets may be less specific, very small sets may be less reliable
- Cross-reference leading edge genes across related pathways to identify key driver genes
- Treat marginal q-values or minimal fold changes as tentative and avoid over-interpretation
