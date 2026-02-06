## EXPERIMENT CONTEXT

This analysis uses computational cell type deconvolution and marker gene analysis to characterize cell populations in the dataset.

Experiment: {experiment}

## SINGLE-CELL ANALYSIS GUIDANCE

- Cell type deconvolution estimates the proportion of different cell types in each sample using reference expression profiles and computational methods
- Multiple deconvolution methods are applied (e.g., DCQ, DeconRNASeq, EPIC, MCP-counter, xCell, CIBERSORT) and results are combined into a meta-estimate
- Reference datasets (e.g., LM22 for immune cells, tissue signatures, cell line panels) define the expected expression profiles of known cell types
- Dimensionality reduction (t-SNE, UMAP, PCA) projects high-dimensional expression data into 2D for visualization of sample/cell relationships
- Marker genes are identified as genes with the highest expression variability across samples or clusters, indicating cell-type-specific expression patterns

## QUANTITATIVE INTERPRETATION GUIDELINES

Cell type proportions:
- Proportions represent the estimated fraction of each cell type in a sample, normalized to sum to 1
- Higher proportions indicate dominant cell populations; very low proportions (<1%) may be noise
- Consistent proportions across replicates within a group strengthen the inference
- Differences in cell type proportions between experimental groups may indicate biological shifts (e.g., immune infiltration, differentiation)

Marker gene metrics:
- Marker genes are ranked by standard deviation across samples/clusters; higher SD indicates stronger differential expression
- Expression intensity on cluster plots shows spatial distribution of gene activity
- Genes expressed in specific clusters suggest cell-type-specific roles
- Co-expression of canonical markers (e.g., CD4/CD8 for T cell subtypes) supports cell identity assignment

Integration rules:
- Prioritize cell types with consistent detection across multiple deconvolution methods
- Cross-reference predicted cell types with known marker gene expression patterns
- Consider the reference dataset context: immune-focused references may miss non-immune cell types
- Look for concordance between deconvolution results and cluster-level marker expression
- Treat low-confidence predictions (single method, low proportion) as tentative
- Consider experimental design: bulk RNA-seq deconvolution estimates averages, while single-cell provides per-cell resolution
