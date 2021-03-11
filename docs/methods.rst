.. _Methods:


Methods
================================================================================

Below are snippets that you can use to describe the methods when using
the Omics Playground. These are just examples and you need to extract
and modify the parts you used and need.


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


Functional analysis
---------------------------

Statistical testing of differential enrichment of genesets was
performed using an aggregation of three statistical methods: fGSEA
(Korotkevich 2019), Camera (Wu 2012) and GSVA/limma (Hanzelmann 2013,
Ritchie 2015). The maximum q-value of the selected methods was taken
as aggregate meta.q value, which corresponds to taking the
intersection of significant genes from all tests. We used more than
50000 genesets from various public databases including: MSigDB
(Subramanian 2005; Liberzon 2015), Gene Ontology (Ashburner 2000), and
Kyoto Encyclopedia of Genes and Genomes (KEGG) (Kanehisa 2000).


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


Scripting and visualization
---------------------------

Data preprocessing was performed using bespoke scripts using R (R Core
Team 2013) and packages from Bioconductor (Huber 2015). Statistical
computation and visualization have been performed using the Omics
Playground version vX.X.X (Akhmedov 2020).

