## EXPERIMENT CONTEXT

This analysis uses time series clustering to identify groups of genes with similar temporal expression patterns.

Experiment: {{experiment}}

## TIME SERIES-SPECIFIC GUIDANCE

- Modules are clusters of genes with similar temporal expression profiles, identified by k-means clustering
- Each module represents a distinct temporal trend (e.g., early upregulation, late downregulation, oscillatory)
- Enriched gene sets within a module reveal the biological processes that follow that temporal pattern
- Genes with high standard deviation across time points are the most dynamically regulated
- The correlation between gene set expression and module expression indicates functional relevance

## QUANTITATIVE INTERPRETATION GUIDELINES

Temporal clustering metrics:
- Genes are clustered based on scaled, centered expression averaged across time points
- Module assignment is determined by k-means clustering of SVD-reduced expression profiles
- Gene-module association strength can be inferred from how closely a gene's profile matches the module average

Enrichment metrics (correlation-based):
- rho > 0.8: very strong correlation with module; 0.6-0.8 strong; 0.4-0.6 moderate; < 0.4 weak
- rho < -0.6: strong anti-correlation (opposite temporal pattern)
- p-value < 0.01: high confidence; 0.01-0.05 standard; 0.05-0.10 marginal; > 0.10 not significant

Differential expression metrics:
- Log2 fold change |logFC| > 2: strong; 1-2 moderate; 0.5-1 mild; < 0.5 minimal
- Q-value < 0.01: high confidence; 0.01-0.05 standard; > 0.05 not significant
- Interaction p-value tests whether the treatment effect changes over time

Integration rules:
- Prioritize gene sets with strong positive correlation (rho > 0.6) and low p-values
- Look for coherent biological themes across multiple enriched gene sets within a module
- Consider temporal context: early-response modules may reflect signaling, late-response modules may reflect adaptation
- Anti-correlated gene sets (negative rho) suggest opposing biological processes
- Interaction effects indicate time-dependent treatment responses, which are often the most biologically interesting
