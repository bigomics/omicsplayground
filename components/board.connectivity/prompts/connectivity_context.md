## EXPERIMENT CONTEXT

This analysis uses gene signature connectivity mapping to identify experiments with similar or opposite expression profiles to the query contrast.

Experiment: {{experiment}}

## CONNECTIVITY-SPECIFIC GUIDANCE

- Connectivity analysis compares a query gene expression signature (fold-change profile) against a reference database of signatures from public datasets
- Similarity is measured using multiple metrics: Pearson correlation (rho), normalized enrichment score (NES), odds ratio, and cosine distance (tau)
- The composite connectivity score combines these metrics as abs(rho) * NES * odd.ratio * abs(tau)
- Positive rho indicates the reference signature changes in the same direction as the query; negative rho indicates opposite direction
- Leading edge genes are those shared between the query and reference signatures among the top differentially expressed genes

## QUANTITATIVE INTERPRETATION GUIDELINES

Connectivity metrics:
- score > 5: very strong similarity; 2-5 strong; 1-2 moderate; < 1 weak
- |rho| > 0.5: strong correlation; 0.3-0.5 moderate; 0.1-0.3 weak; < 0.1 negligible
- NES > 2: strong enrichment; 1.5-2 moderate; 1-1.5 mild; < 1 not enriched
- Positive rho: concordant regulation (similar biological response)
- Negative rho: discordant regulation (opposite biological response, potential drug reversal)

Biological interpretation:
- Signatures with high positive connectivity suggest shared biological mechanisms or disease processes
- Signatures with high negative connectivity suggest opposing mechanisms, relevant for drug repurposing (reversal signatures)
- Clusters of similar signatures from the same tissue or disease context strengthen biological relevance
- Leading edge genes shared across multiple similar signatures are likely key drivers of the shared biology
- Consider the source of reference signatures (cell line, tissue, treatment) when interpreting biological relevance

Integration rules:
- Prioritize signatures with both high composite score and high absolute correlation
- Look for thematic coherence among top similar experiments (e.g., same disease, same pathway perturbation)
- Cross-reference leading edge genes across top similar signatures to identify core driver genes
- Consider cumulative fold-change and cumulative enrichment patterns to identify robust cross-study signals
- Treat low-score matches or single-metric hits as tentative and avoid over-interpretation
