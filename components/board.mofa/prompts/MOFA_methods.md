## EXPERIMENT CONTEXT

This analysis uses Multi-Omics Factor Analysis (MOFA) to identify latent factors that capture the major sources of variation across multiple omics data types.

Experiment: {experiment}

## MOFA-SPECIFIC GUIDANCE

- MOFA factors are latent variables that capture shared and unique sources of variation across multiple data modalities (e.g., transcriptomics, proteomics, metabolomics)
- Factor loadings (weights) indicate the contribution of each feature to a given factor; high absolute weight means the feature is strongly associated with that factor
- Factor scores represent sample-level values for each factor and can be correlated with phenotypes/traits
- Enriched pathways for a factor reveal the biological processes captured by that factor
- Variance explained per factor and per view indicates which factors and data types drive the main heterogeneity

## QUANTITATIVE INTERPRETATION GUIDELINES

Definitions and relationships:
- Factor loadings (weights) quantify the relationship between each feature and a latent factor; features with high absolute loadings are the most important contributors to that factor
- Factor scores represent each sample's position along a factor axis; differences in scores between conditions suggest the factor captures condition-specific variation
- Variance explained (R2) measures how much of the total variance in each view is captured by each factor
- Factor-trait correlation quantifies the association between a factor and a phenotype/trait variable

Important note on thresholds:
- MOFA does not prescribe universal cutoffs; interpretation depends on the dataset and number of factors
- Treat any numeric threshold as a guideline, not a rule

Enrichment metrics:
- NES (Normalized Enrichment Score): |NES| > 2.0 very strong; 1.5-2.0 strong; 1.0-1.5 moderate; < 1.0 weak
- padj < 0.01: high confidence; 0.01-0.05 standard; 0.05-0.10 marginal; > 0.10 not significant
- Pathway size: larger gene sets provide more robust enrichment estimates; very small sets (< 10) should be interpreted cautiously

Gene/feature metrics:
- Loading weight: features with the highest absolute weights are key drivers of the factor
- Centrality: high centrality features are well-connected in the factor's feature network
- Consider both the magnitude and sign of loadings when interpreting factor biology

Integration rules:
- Prioritize factors with high variance explained and significant trait correlations
- Favor pathways with strong NES and low padj that are consistent with the top loading features
- When multiple omics layers contribute to a factor, look for convergent biological themes across data types
- Treat marginal enrichments or low-weight features as tentative and avoid over-interpretation
