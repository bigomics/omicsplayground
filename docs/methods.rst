.. _Methods:


Methods
================================================================================


Below are snippets that you can use to describe the methods when using
the Omics Playground. These are just examples and you need to extract
and modify the parts you used and need.


Batch correction 
-----------------

Batch effects, or contamination by unwanted variables, was identified
by an F-test for the first three principal components. Continuous
variables were dichotomized into high/low before testing. Highly
confounding variables would appear as having high relative
contribution in the first or second principal component, often higher
than the variable of interest. Batch effects were also visually
assessed (before and after correction) using annotated heatmaps and
t-SNE plots colored by variables.

[Setting: 'low'] Batch correction was performed for explicit batch
variables or unwanted covariates. Parameters with a correlation r>0.3
with any of variables of interest (i.e. the model parameters) were
omitted from the regression. Correction was performed by regressing
out the covariate using the 'removeBatchEffect' function in the limma
R/Bioconductor package.

[Setting: 'medium'] Additional batch correction was performed for for
intrinsic technical parameters such as library size (i.e. total
counts), mitochondrial and ribosomal proportions, cell cycle and
gender. These parameters were estimated from the data. The cell cycle
was estimated using the Seurat R/Bioconductor package. Gender (if not
given) was estimated by checking the absence/presence of expression of
gender specific genes on the X/Y chromosome. Parameters with a
correlation r>0.3 with any of the model parameters were omitted from
the regression. Correction was performed by regressing out the
covariate using the 'removeBatchEffect' function in the limma
R/Bioconductor package.

[Setting: 'high'] Unsupervised batch correction was performed using
'surrogate variable analysis' (SVA) (Leek 2007) by estimating the
latent surrogate variables and regressing out using the
'removeBatchEffect' function in the limma R/Bioconductor package.


Clustering
---------------------------

Heatmaps were generated using the ComplexHeatmap R/Bioconductor
package (Gu 2016) on scaled log-expression values (z-score) using
euclidean distance and Ward linkage. The standard deviation was used
to rank the genes for the reduced heatmaps.

T-distributed stochastic neighbor embedding (t-SNE) was computed using
the top 1000 most varying genes, then reduced to 50 PCA dimensions
before computing the t-SNE embedding. The perplexity heuristically set
to 25% of the sample size or 30 at maximum, and 2 at
minimum. Calculation was performed using the `Rtsne` R package.

Uniform Manifold Approximation and Projection (UMAP) was computed
using the top 1000 most varying genes, then reduced to 50 PCA
dimensions before computing the UMAP embedding. The number of
neighbours was heuristically set to 25% of the sample size or 30 at
maximum, and 2 at minimum. Calculation was performed using the `uwot`
R package.

Principal component analysis (PCA) was computed using the `irlba` R
package.


Statistical testing
---------------------------

Multi-method statistical testing. For gene-level testing, statistical
significance was assessed using three independent statistical methods:
DESeq2 (Wald test), edgeR (QLF test) and limma-trend (Love 2014;
Robinson 2010; Ritchie 2015). The maximum q-value of the three methods
was taken as aggregate q-value, which corresponds to taking the
intersection of significant genes from all three tests.


Statistical testing of differential enrichment of genesets was
performed using an aggregation of multiple statistical methods:
Fisher's exact test, fGSEA (Korotkevich 2019), Camera (Wu 2012) and
GSVA/limma (Hanzelmann 2013, Ritchie 2015). The maximum q-value of the
selected methods was taken as aggregate meta.q value, which
corresponds to taking the intersection of significant genes from all
tests. As each method uses different estimation parameters (NES for
GSEA, odd-ratio for fisher, etc.) for the effect size, for
consistency, we took the average log fold-change of the genes in the
geneset as sentinel value. We used more than 50000 genesets from
various public databases including: MSigDB (Subramanian 2005; Liberzon
2015), Gene Ontology (Ashburner 2000), and Kyoto Encyclopedia of Genes
and Genomes (KEGG) (Kanehisa 2000).


Functional analysis
---------------------------

Graph-weighted GO analysis. The enrichment score of a GO term was
defined as the sum of q-weighted average fold-changes, (1-q)*logFC, of
the GO term and all its higher order terms along the shortest path to
the root in the GO graph. The fold-change of a gene set was defined as
the average of the fold-change values of its member genes. This
graph-weighted enrichment score thus reflects the enrichment of a GO
term with evidence that is corroborated by its parents in the GO graph
and therefore provides a more robust estimate of enrichment. The
activation map visualizes the scores of the top-ranked GO terms for
multiple contrasts as a heatmap.

KEGG pathway visualization was performed using the Pathview
R/Bioconductor package using the foldchange as node color.



Scripting and visualization
---------------------------

Data preprocessing was performed using bespoke scripts using R (R Core
Team 2013) and packages from Bioconductor (Huber 2015). Statistical
computation and visualization have been performed using the Omics
Playground version vX.X.X (Akhmedov 2020).



REFERENCES 
---------------------------

Akhmedov M, Martinelli A, Geiger R and Kwee I. Omics Playground: A
comprehensive self-service platform forvisualization, analytics and
exploration of Big Omics Data. NAR Genomics and Bioinformatics, Volume
2, Issue 1, March 2020,

Ashburner et al. Gene ontology: tool for the unification of
biology. Nat Genet. May 2000;25(1):25-9.


Huber W, et al. (2015) Orchestrating high-throughput genomic analysis
with Bioconductor. Nature Methods 12:115-121; doi:10.1038/nmeth.3252

Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and
Genomes. Nucleic Acids Res. 28, 27-30 (2000).

Leek J., Storey J. Capturing heterogeneity in gene expression studies
by ‘surrogate variable analysis’ PLoS Genet. 2007

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold
change and dispersion for RNA-seq data with DESeq2.” Genome Biology,
15, 550.

R Core Team (2013). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
URL http://www.R-project.org/.

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK
(2015). “limma powers differential expression analyses for
RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7)


Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor
package for differential expression analysis of digital gene
expression data.” Bioinformatics, 26(1), 139-140.
