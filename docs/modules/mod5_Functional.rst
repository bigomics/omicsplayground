.. _Functional:

Functional Analysis
================================================================================
This module performs specialized pathway and drug enrichment analysis. 
It contains three panels where it provides higher level functional and 
visual analysis of the contrast space using the 
`KEGG <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102409/>`__ graph structure
in the **KEGG pathways** panel. Under the **GO** panel, very similar functional
analysis is done using the Gene Ontology (`GO <http://geneontology.org/>`__) 
graph structure. 
Given a particular contrast profile, it also searches for the closest 
drug profiles from the `L1000 <https://www.ncbi.nlm.nih.gov/pubmed/29195078>`__
drug expression database under the **Drug Connectivity Map** panel.

.. note::

    This module is supported in the EXPERT MODE ONLY.


Input slider
--------------------------------------------------------------------------------
It is possible to more information about the module in the ``Info`` from the 
input slider. Users can specify the contrast of their interest in 
the ``Contrast`` settings. Under the main *Options*, users can select
``normalize activation matrix`` to fine-tune the coloring of an activation 
matrices and ``filter significant (tables)`` to filter the significant entries
in the tables.

.. figure:: figures/psc6.0.png
    :align: center
    :width: 30%


KEGG pathways
--------------------------------------------------------------------------------
`KEGG <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102409/>`__ is a collection
of manually curated pathways representing the current knowledge of molecular 
interactions, reactions and relation networks as pathway maps. In the 
**KEGG pathway** panel, each pathway is scored for the selected contrast profile
and reported in the table. A unique feature of the platform is that it provides 
an activation-heatmap comparing the activation levels of pathways across multiple
contrast profiles. This facilitates to quickly see and detect the similarities 
between profiles in certain pathways. More detailed explaination of each output
is provided below.

:**a**: In the pathway map, genes are colored according to their upregulation 
        (red) or downregulation (blue) in the contrast profile. Each pathway 
        is scored for the selected contrast profile and reported in the table 
        below.

:**b**: Enrichment table. The table is interactive; enabling user to sort on 
        different variables and select a pathway by clicking on the row in the 
        table. The scoring is performed by considering the total number of genes
        in the pathway (:option:`n`), the number of genes in the pathway supported by the 
        contrast profile (:option:`k`), the ratio of :option:`k/n`, and the ratio of 
        :option:`|upregulated or downregulated genes|/k`. Additionally, the table 
        contains the list of the upregulated and downregulated genes for each
        pathway and a :option:`q` value from the Fisher's test for the overlap.

:**c**: The KEGG activation matrix visualizes the activation levels of pathways
        (or pathway keywords) across multiple contrast profiles. This facilitates
        to quickly see and detect the similarities of certain pathways between
        contrasts. The size of the circles correspond to their relative activation,
        and are colored according to their upregulation (red) or downregulation
        (blue) in the contrast profile.

.. figure:: figures/psc6.1.png
    :align: center
    :width: 100%


GO graph
--------------------------------------------------------------------------------
In the **GO** panel, users can perform `GO <http://geneontology.org/>`__ analysis.
GO defines functional concepts/classes and their relationships as a hierarchical
graph. 
The GO database provides a computational representation of the current knowledge 
about roles of genes for many organisms in terms of molecular functions, cellular
components and biological processes. All the features described under the 
**KEGG pathway** panel, such as scoring the gene sets and drawing an 
activation-heatmap,
can be performed for the GO database under the GO graph tab. Instead of pathway
maps, an annotated graph structure provided by the GO database is potted for
every selected gene set. 
Each output chart/table of the panel is describer below in detail.

:**a**: The structure of GO can be described in terms of a graph, where each
        GO term is a node, and the relationships between the terms are edges 
        between the nodes. GO is loosely hierarchical, with 'child' terms being
        more specialized than their 'parent' terms. The graph is interactive. 
        You can move the graph and zoom in using the mouse.
        Under the graph *Settings*, users can select ``Prune tree`` to prune
        the tree only with significant branches and ``color custers`` to 
        highlight clusters with different colors

        .. figure:: figures/psc6.2.a.png
            :align: center
            :width: 35%

:**b**: GO score table. The scoring of a GO term is performed by considering
        the cumulative score of all terms from that term to the root node. 
        That means that GO terms that are supported by higher level terms
        levels are preferentially scored.

:**c**: The GO activation matrix visualizes the activation of GO terms
        across conditions. From this figure, you can easily detect GO terms
        that are consistently up/down across conditions. The size of the circles
        correspond to their relative activation, and are colored according to 
        their upregulation (red) or downregulation (blue) in the contrast
        profile.

.. figure:: figures/psc6.2.png
    :align: center
    :width: 100%

    
Drug C-Map
--------------------------------------------------------------------------------
In the **Drug Connectivity Map** panel, users can correlate their signature with
more than 5000 known drug profiles from the 
`L1000 <https://www.ncbi.nlm.nih.gov/pubmed/29195078>`__ database. 
An activation-heatmap compares drug activation profiles across multiple contrasts. 
This facilitates to quickly see and detect the similarities between contrasts
for certain drugs.

:**a**: The Drug Connectivity Map correlates your signature with more than 
        5000 known drug profiles from the L1000 database, and shows the top
        N=10 similar and opposite profiles by running the GSEA algorithm on 
        the contrast-drug profile correlation space. Under the plots *Settings*,
        users can select the type of drug enrichment analysis: ``mono`` or 
        ``combo`` (if available).

        .. figure:: figures/psc6.3.a.png
            :align: center
            :width: 35%

:**b**: Drug profile enrichment table. Enrichment is calculated by correlating
        your signature with more than 5000 known drug profiles from the L1000
        database. Because the L1000 has multiple perturbation experiment for a
        single drug, drugs are scored by running the GSEA algorithm on the 
        contrast-drug profile correlation space. In this way, we obtain a 
        single score for multiple profiles of a single drug.

:**c**: This plot visualizes the mechanism of action (MOA) across the enriched
        drug profiles. On the vertical axis, the number of drugs with the same
        MOA are plotted. You can switch to visualize between MOA or target gene.
        Under the plots *Settings*, users can select the plot type of MOA
        analysis: by class description (``drug class``) or by target gene 
        (``target gene``).

        .. figure:: figures/psc6.3.c.png
            :align: center
            :width: 35%

.. figure:: figures/psc6.3.png
    :align: center
    :width: 100%
    
    
    