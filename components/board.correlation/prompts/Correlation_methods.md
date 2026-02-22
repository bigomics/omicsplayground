## EXPERIMENT CONTEXT

This analysis uses gene-gene and sample-gene correlation analysis to identify co-expressed features and potential regulatory relationships.

Experiment: {{experiment}}

## CORRELATION-SPECIFIC GUIDANCE

- Correlation analysis computes Pearson correlation coefficients between a selected gene and all other features across samples
- Partial correlation (via the Glasso algorithm) adjusts for indirect effects, estimating direct gene-gene interactions
- The correlation table reports both canonical Pearson correlation (cor) and partial correlation (pcor) for each feature
- Partial correlation removes confounding effects of other variables, providing a more specific measure of direct association
- Gene filtering can restrict the analysis to specific gene families or custom gene lists

## QUANTITATIVE INTERPRETATION GUIDELINES

Correlation metrics:
- |cor| > 0.8: very strong correlation; 0.6-0.8 strong; 0.4-0.6 moderate; 0.2-0.4 weak; < 0.2 negligible
- |pcor| > 0.5: strong direct association; 0.3-0.5 moderate; 0.1-0.3 weak; < 0.1 negligible
- Positive correlation: genes are co-expressed (increase/decrease together across samples)
- Negative correlation: genes are inversely expressed (one increases as the other decreases)
- Partial correlation values are typically smaller in magnitude than canonical correlation values

Network and module context:
- Highly correlated gene clusters may indicate co-regulation or shared biological function
- Genes with strong partial correlation are more likely to have direct regulatory relationships
- The partial correlation network reveals the local neighborhood of direct interactions around the selected gene
- Hub genes (many strong partial correlations) may serve as key regulators in the network

Integration rules:
- Prioritize genes with both strong canonical and partial correlation for the most robust associations
- Genes with high canonical but low partial correlation may be indirectly associated through a common regulator
- Genes with relatively high partial correlation but moderate canonical correlation suggest direct but subtle interactions
- Consider the biological function of top correlated genes to identify coherent functional modules
- Cross-reference correlated genes with known pathways to infer biological relevance
