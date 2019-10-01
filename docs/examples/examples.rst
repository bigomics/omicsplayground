.. _examples:

Reanalyzing Public Datasets
================================================================================
To illustrate the use case of the Omics Playground, we reanalyzed different types
of publics datasets, including microarray, bulk RNA-seq, single-cell RNA-seq and
proteomic datasets to recapitulate the results.


Single-cell RNA-seq data
--------------------------------------------------------------------------------
For single-cell RNA-seq data, we downloaded the melanoma data set 
`GSE72056 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056>`__ of
`Tirosh et al. <https://www.ncbi.nlm.nih.gov/pubmed/27124452>`__.
Our platform recapitulates well the original findings of the paper. 
The t-SNE clustering (`Figure 1`_) separates the different cell types. 
`Figure 2`_ and `Figure 3`_ show the volcano plot, MA plot and most differentially
expressed genes between malignant and non-malignant cells. 
The CNV map (`Figure 4`_) confirms the major chromosomal copy number 
variations found in the malignant cells. `Figure 5`_ shows high enrichment
of a immune checkpoint signature, particularly concentrated in the T cells.
The biomarker heatmap (`Figure 6`_) highlights the marker genes for
each cell type. Each gene cluster is furthermore automatically
annotated with the most correlated gene sets (`Figure 7`_).


tSNE plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`Figure 1`:

.. figure:: figures/fig2_a.png
    :align: center
    :width: 100%

**Figure 1**. The t-SNE clustering with cell type annotation. 
To reproduce the same figure on the platform, select and load ``GSE72056-scmelanoma`` dataset, 
and go to the **PCA/tSNE** panel of the **Clustering** module. From the plot *Settings*, 
set the ``color``: group and ``layout``: tsne.


Volcano and MA plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`Figure 2`:

.. figure:: figures/fig2_b.png
    :align: center
    :width: 100%    

**Figure 2**. Volcano and MA plot for the malignant versus non-malignant contrast.
To reproduce the same figure on the platform, select and load ``GSE72056-scmelanoma`` dataset, 
and go to the **Plots** panel of the **Expression** module. From the input slider, 
set the ``Contrast``: yes_vs_no, ``Gene family``: all, ``FDR``: 0.2, and 
``logFC threshold``: 0.5.    


Differentially expressed genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`Figure 3`:

.. figure:: figures/fig2_c.png
    :align: center
    :width: 100%     

**Figure 3**. Barplot of corresponding differentially expressed genes.
To reproduce the figure on the platform, select and load ``GSE72056-scmelanoma`` dataset, 
and go to the **Top genes** panel of the **Expression** module. From the input slider, 
set the ``Contrast``: yes_vs_no, ``Gene family``: all, ``FDR``: 0.2, and 
``logFC threshold``: 0.5.
    

Inferred copy number
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`Figure 4`:

.. figure:: figures/fig2_d.png
    :align: center
    :width: 100% 

**Figure 4**. Inferred copy number for sample Cy80.
To reproduce the figure on the platform, select and load ``GSE72056-scmelanoma`` dataset, 
and go to the **CNV** panel of the **scProfiling** module. From the plot *Settings*, 
set the ``Annotate with``: malignant and ``Order samples by``: clust.
    
    
Immune checkpoint signature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`Figure 5`:

.. figure:: figures/fig2_e.png
    :align: center
    :width: 100%   

**Figure 5**. Enrichment distribution for an immune checkpoint signature showing high
enrichment in T and B cells .
To reproduce the figure on the platform, select and load ``GSE72056-scmelanoma`` dataset, 
and go to the **Marker** panel in the **Signature** module. From the input slider, 
set the ``Contrast``: custom and ``Signature``: immune_chkpt from the provided sample list.
    

Biomarker heatmap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`Figure 6`:

.. figure:: figures/fig2_f.png
    :align: center
    :width: 100% 

**Figure 6**. Biomarker heatmap for non-malignant cells.
To reproduce the figure on the platform, select and load ``GSE72056-scmelanoma`` dataset, 
and go to the **Heatmap** panel in the **Clustering** module. From the input slider, 
set the ``Level``: gene, ``Features``: all, and ``Filter samples``: {cell.type=Bcell,
cell.type=CAF, cell.type=endothelial, cell.type=Macrophage, cell.type=NK, cell.type=Tcell}.
In the plot *Settings*, set ``Plot type``: ComplexHeatmap, ``split by``: cell.type,
``top mode``: specific, ``top N``: 50 and ``scale``: relative.
    

Annotate heatmap clusters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`Figure 7`:

.. figure:: figures/fig2_g.png
    :align: center
    :width: 100%     

**Figure 7**. Enrichment annotation of corresponding heatmap clusters from the `Figure 6`_.
To reproduce the figure on the platform, generate the heatmap in `Figure 6`_ first, 
then go to the **Annotate clusters** panel. From the plot *Settings*, 
set the ``Reference level``: geneset and ``Reference set``: GOBP.



Bulk RNA-seq Data
--------------------------------------------------------------------------------
To elucidate the mechanism of action of a new drug, or for the intention of drug
repurposing, it is often useful to find other drugs that have similar or opposing
signatures compared to some given fold change profile. 
This example is illustriated using the RNA-sequencing dataset 
`GSE114716 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114716>`__ from
`Goswami et al. <https://www.ncbi.nlm.nih.gov/pubmed/29905573>`__, 
which contains CD4 T cells following ipilimumab therapy.
`Figure 8`_ shows the top ranked drugs with most similar or most opposing signatures
to ipilimumab, a novel monoclonal antibody targeting CTLA-4 used in tumour therapy.
The list contains several known anti-tumoral drugs, such as bortezomib and 
palbociclib, but also highlights relationships with other compounds not normally
used in tumor therapy, such as emetine, an anti-protozoal drug with anti-tumoral
properties (41, 43).

Drug enrichment profiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`Figure 8`: 

.. figure:: figures/fig2_h.png
    :align: center
    :width: 100% 

**Figure 8**. Drug enrichment profiles for most similar and opposing drugs
compared to ipilimumab treatment.



    
Microarray Data
--------------------------------------------------------------------------------
In this section, we perform the biomarker selection and survival analysis using
the `GSE10846 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846>`__
from `Lenz et al. <https://www.ncbi.nlm.nih.gov/pubmed/19038878>`__,
which is the microarray gene expression dataset of diffuse large B-cell 
lymphoma (DLBCL) patients.

Fig. 3b shows a microarray gene expression data set, GSE10846 (44), of diffuse
large B-cell lymphoma (DLBCL). Fig. 3c and 3d show the variable importance 
plot and a survival tree on the overall survival of the DLBCL patients, respectively.

Biomarker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Survival analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






Proteomic Data
--------------------------------------------------------------------------------

Proteome profiles of activated vs resting human naive T cells at different 
times (Geiger et al., Cell 2016).

Fig. 3a shows the volcano plots corresponding to eight different statistical tests
comparing time-dependent activation of T cells at 48h vs. 12h (42). We see that 
both standard t-test and the Welch t-test show much less power to detect 
significant genes compared to the other methods. The result from edgeR-QLF 
is close to those of the two limma based methods, while edgeR-LRT is very 
similar to the results of DESeq2-Wald. 


With larger data sets, often the number of contrasts increases and complicates 
the overall analysis. For example, the proteomics data set of 
`Rieckmann et al. 2017 <https://www.ncbi.nlm.nih.gov/pubmed/28263321>`__
comprises 26 populations of seven major immune cell types, measured during resting
and activated states. There are more than 300 possible comparisons to make.
To gain a better overview, gene set connectivity heatmaps (Fig. 3e) help 
visualize the similarities between multiple contrasts on a functional level. 
Alternatively, similarities can be visualized as a connectivity graph (Fig. 3f). 
For the same data set, Fig. 3g shows a computed partition tree that classifies 
the major cell types.


Connectivity heatmap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Connectivity graph
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Classification tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


