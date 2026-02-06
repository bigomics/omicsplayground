## EXPERIMENT CONTEXT

This analysis uses Weighted Gene Co-expression Network Analysis (WGCNA) to identify modules of co-expressed genes.

Experiment: {experiment}

## WGCNA-SPECIFIC GUIDANCE

- Modules are clusters of genes with highly correlated expression patterns
- Module eigengenes summarize the expression profile of each module
- Module-trait correlations indicate potential biological relevance
- Hub genes have high intramodular connectivity and may be key regulators
- Enriched pathways suggest the biological function of the module

## QUANTITATIVE INTERPRETATION GUIDELINES

Definitions and relationships:
- Module eigengene (ME) is the first principal component of the module and represents its expression profile
- Module membership (MM, kME) is the correlation between a gene and the module eigengene; high absolute MM indicates strong module membership
- Gene significance (GS) is commonly defined as the absolute correlation between a gene and a trait
- Module significance is the average GS across genes in a module
- Hub genes are highly connected intramodular genes and typically have high absolute MM

Important note on thresholds:
- Core WGCNA sources define metrics but do not prescribe universal numeric cutoffs
- Application studies use explicit, dataset-specific heuristics; treat any numeric threshold as a guideline, not a rule

Enrichment metrics:
- Score > 0.8: very strong enrichment; 0.6-0.8 strong; 0.4-0.6 moderate; < 0.4 weak
- Q-value < 0.01: high confidence; 0.01-0.05 standard; 0.05-0.10 marginal; > 0.10 not significant
- Overlap: higher overlap with smaller denominators is more specific; low overlap suggests weaker specificity

Gene metrics:
- MM (kME) > 0.8: core hub gene; 0.6-0.8 strong member; 0.4-0.6 moderate; < 0.4 peripheral
- GS |r| > 0.7 strong; 0.5-0.7 moderate; 0.3-0.5 weak; < 0.3 minimal
- Log2 fold change |LogFC| > 2 strong; 1-2 moderate; 0.5-1 mild; < 0.5 minimal
- Centrality > 0.8 highly connected; 0.6-0.8 well connected; 0.4-0.6 moderate; < 0.4 peripheral

Integration rules:
- Prioritize genes with both high MM and high GS in trait-related modules
- Favor pathways with strong scores and low q-values that match hub gene functions
- Treat marginal q-values or low MM/GS values as tentative and avoid over-interpretation
