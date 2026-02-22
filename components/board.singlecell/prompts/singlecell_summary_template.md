## Task

Interpret the cell type composition and marker gene analysis results for the experiment based on the single-cell / cell profiling data provided below.

## Analysis Overview

**Reference dataset:** {{refset}}

**Deconvolution method:** {{dcmethod}}

The cell profiling analysis estimates the proportion of different cell types in each sample using computational deconvolution against the selected reference dataset. Marker genes highlight the most variable genes across samples or clusters.

## Top Predicted Cell Types

The following cell types were identified as the dominant populations across all samples:

{{cell_types}}

Focus on identifying:
1. The dominant cell populations and their relative abundance
2. How cell type composition varies across experimental groups
3. Whether the predicted cell types are consistent with the experimental context

## Marker Genes

{{marker_genes}}

## Phenotype Associations

{{phenotype_info}}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{{refset}}** | Cell profiling summary

Then synthesize the above information to explain:
- The primary cell type composition of the samples and how it relates to the experimental context
- Which cell populations show the most variation across groups or conditions
- How marker gene expression patterns support or refine the cell type predictions
- The biological significance of the observed cell type proportions and any notable shifts between conditions
- Use quantitative metrics to support your interpretation (cite specific proportions, marker genes, and group differences when available)
