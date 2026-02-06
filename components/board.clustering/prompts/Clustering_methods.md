## EXPERIMENT CONTEXT

This analysis uses unsupervised clustering to discover patterns, subgroups, and gene modules in high-dimensional omics data.

Experiment: {experiment}

## CLUSTERING-SPECIFIC GUIDANCE

- Hierarchical clustering simultaneously groups genes (rows) and samples (columns) based on expression similarity
- Gene modules are groups of co-expressed genes identified by cutting the hierarchical dendrogram at a specified number of clusters (K)
- Dimensionality reduction methods (PCA, t-SNE, UMAP) project high-dimensional data into 2D for visual exploration of sample relationships
- Functional annotation of gene modules uses correlation with reference gene sets (e.g., Hallmark, GO, KEGG) to assign biological meaning
- Marker gene selection identifies genes that are specifically overexpressed in one group compared to all others

## QUANTITATIVE INTERPRETATION GUIDELINES

Cluster annotation metrics:
- Correlation (R) > 0.5: strong association between gene module and annotation term
- Correlation (R) 0.3-0.5: moderate association
- Correlation (R) 0.1-0.3: weak association
- Correlation (R) < 0.1: negligible association

Gene module interpretation:
- Modules with coherent functional annotations indicate coordinated biological processes
- Large modules (many genes) may reflect broad transcriptional programs; small modules may capture specific pathways
- Modules that separate sample groups by phenotype suggest phenotype-associated expression programs

Sample clustering interpretation:
- Samples clustering by known phenotype confirms expected biological variation
- Unexpected sample groupings may reveal batch effects, outliers, or novel biological subtypes
- Dimensionality reduction preserving group separation indicates strong transcriptomic differences between conditions

Integration rules:
- Prioritize gene modules with strong, coherent functional annotations
- Look for correspondence between gene modules and sample phenotype groupings
- Consider whether clustering patterns align with known biology of the experimental system
- Identify marker genes that drive cluster separation as candidate biomarkers or regulators
- Treat modules with weak or incoherent annotations as potentially noise-driven
