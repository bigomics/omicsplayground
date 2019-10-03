.. _Signature:

Signature Analysis
================================================================================
In the **Signature Analysis** module, users can test their gene signature by
calculating an enrichment score. 

After uploading a gene list, the **Markers** panel produces a t-SNE plot of 
samples for each gene. The **Enrichment** panel performs the enrichment analysis
of the gene list against all contrasts by running the GSEA algorithm and plots 
enrichment outputs. Under the **Overlap/similarity** panel, users can find the 
similarity of their gene list with all the gene sets and pathways in the platform.

.. note::

    This module is supported in the EXPERT MODE ONLY.


Input panel
--------------------------------------------------------------------------------
Users need to specify the contrasts of their interest to start the analysis in 
the ``Contrast`` settings. They can use a sample list provided on the platform
or upload their own gene list. Instead of a short list, a contrast profile can 
also be selected, which is a complete gene list resulted from one of the contrasts
in the analysis.
Under the main *Options*, users can set the ``Enrichment method``, where
``rcor`` performs a correlation based computation while ``fgsea`` runs the
fast implementation of `GSEA <https://www.biorxiv.org/content/10.1101/060012v1.full>`__ 
algorithm.

.. figure:: figures/psc8.0.png
    :align: center
    :width: 30%


Markers
--------------------------------------------------------------------------------
After uploading a gene list, the **Markers** panel produces a t-SNE plot of 
samples for each gene, where the samples are colored with respect to the 
upregulation (in red) or downregulation (in blue) of that particular gene.

.. figure:: figures/psc8.1.png
    :align: center
    :width: 100%


Enrichment
--------------------------------------------------------------------------------
The **Enrichment** panel performs the enrichment analysis of the gene list 
against all contrasts by computing a correlation based enrichment or running the
`GSEA <https://www.biorxiv.org/content/10.1101/060012v1.full>`__ 
algorithm and plots enrichment outputs. Under the plot *Settings*, users can
quickly check the enrichment of their gene list in other contrasts from 
the relevant public datasets by setting the ``Test dataset``.

.. figure:: figures/psc8.2.0.png
    :align: center
    :width: 30%

The enrichment plots are shown below. They show the enrichment of the query 
signature across all constrasts. Positive enrichment means that this particular
contrast shows similar expression changes as the query signature.
Furthermore, the enrichment statistics can be found in the right tables, where

:**a**: Reports the summary of correlation/enrichment of the query signature 
        in all contrasts. 
:**b**: Reports the summary of fold-changes of genes in the query signature.

.. figure:: figures/psc8.2.png
    :align: center
    :width: 100%


Overlap/similarity
--------------------------------------------------------------------------------
Under the **Overlap/similarity** panel, users can compare their gene
list with all the gene sets and pathways in the platform through
overlap analysis, or also known as over-representation analysis. The
significance of overlap is computed by the Fisher's exact test. A
score is computed as the geometric mean of the absolute logarithm of
the odds ratio and q-value of the Fisher's test.

The table reports the :option:`score`, total number of genes in the
gene set (:option:`K`), the number of intersecting genes between the
list and the gene set (:option:`k`), the overlapping ratio of
:option:`k/K`, as well as the :option:`odds.ratio` and
:option:`q.fisher` values by the Fisher's test for the overlap test.

.. figure:: figures/psc8.3.png
    :align: center
    :width: 100%

Under the plot settings, users can specify the number to top features
to show, or users can select to hide/show the feature names in the plot.
	    
.. figure:: figures/psc8.3.0.png
    :align: center
    :width: 100%
	   

