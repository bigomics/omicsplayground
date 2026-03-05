## EXPERIMENT CONTEXT

This analysis uses pathway enrichment analysis to identify biological pathways and gene sets with coordinated expression changes.

Experiment: {{experiment}}

## PATHWAY ANALYSIS GUIDANCE

- Pathway enrichment analysis tests whether predefined sets of genes show statistically significant, concordant differences between biological states
- Multiple pathway databases are used: Gene Ontology (GO), WikiPathways, and Reactome
- Each pathway is scored using multiple enrichment methods (e.g., fgsea, camera, fisher) and results are combined via meta-analysis
- The meta log fold change (logFC) indicates the direction and magnitude of pathway activation
- The meta q-value represents the combined statistical significance across methods

## QUANTITATIVE INTERPRETATION GUIDELINES

Enrichment metrics:
- |logFC| > 1.5: very strong activation/repression; 1.0-1.5 strong; 0.5-1.0 moderate; < 0.5 weak
- meta.q < 0.001: very high confidence; 0.001-0.01 high; 0.01-0.05 standard; 0.05-0.10 marginal; > 0.10 not significant
- Positive logFC indicates pathway upregulation; negative logFC indicates pathway downregulation

Pathway database specifics:
- GO terms are organized hierarchically (Biological Process, Molecular Function, Cellular Component); child terms are more specific than parent terms
- Reactome pathways represent curated reaction networks with detailed molecular mechanisms
- WikiPathways are community-curated pathway maps often capturing disease-specific or emerging biology

Integration rules:
- Prioritize pathways with both strong effect sizes (high |logFC|) and strong statistical evidence (low q-value)
- Look for convergence across pathway databases (same biological theme enriched in GO, Reactome, and WikiPathways)
- Group related pathways to identify dominant biological themes rather than interpreting each pathway independently
- Treat marginal q-values (0.05-0.10) or weak effect sizes (|logFC| < 0.5) as tentative and avoid over-interpretation
