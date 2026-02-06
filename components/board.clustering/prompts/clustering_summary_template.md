## Task

Interpret the biological significance of the unsupervised clustering analysis based on the gene modules, cluster annotations, and sample groupings provided below.

## Analysis Overview

**Experiment:** {experiment}

**Clustering method:** {cluster_method}

**Number of gene modules (K):** {n_modules}

The analysis uses hierarchical clustering to group co-expressed genes into modules and identify patterns across samples.

## Gene Module Summary

{cluster_info}

## Functional Annotation of Gene Modules

The following annotation terms show the strongest correlation with each gene module:

{enriched_terms}

## Top Marker Genes

{top_genes}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**Clustering analysis** | {experiment}

Then synthesize the above information to explain:
- The major gene expression programs revealed by the clustering analysis
- How functional annotations characterize each gene module biologically
- Whether the gene modules correspond to known biological processes or phenotype differences
- Key marker genes that define each module and their potential biological roles
- Use quantitative metrics to support your interpretation (cite specific correlation values and gene counts when available)
