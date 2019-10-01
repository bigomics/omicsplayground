.. _Dataprep:

Data cleaning and preprocessing
================================================================================

The data cleaning and preprocessing includes preparing the input data, filtering, 
normalising, and precomputing statistics for some analyses. The data cleaning and 
preprocessing is performed offline using scripts in order to support real-time 
interaction and minimize user interface latency.

    
Input data
--------------------------------------------------------------------------------
The platform requires the transcriptomics and proteomics data to be in a 
structured format as an input. Users can prepare an input data from
their own FASTQ files, gene counts tables, or from a dataset of interest stored 
in public repositories such as `GEO <https://www.ncbi.nlm.nih.gov/geo/>`__.
Similarly, they can also prepare an input from LC-MS/MS proteomics data.
For all of these cases, the platform comes with the necessary scripts for data 
cleaning and preprocessing under the ``/scripts`` folder.


Filtering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The data preprocessing includes some filtering criteria, such as filtering of 
genes based on variance, the expression across the samples, and the number of 
missing values. Similarly, samples can also be filtered based on the read quality, 
total abundance, unrelated phenotype, or an outlier criterion.


Normalisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The raw counts are converted into counts per million (CPM) and log2. Depending on 
the data set, a quantile normalization can be applied. Known batches in the data 
can be corrected with limma or ComBat. Other unknown batch 
effects and unwanted variation can be further removed using surrogate variable 
analysis in the sva package.


Offline computation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Statistics for the differentially expressed genes analysis and gene set enrichment
analysis are precomputed to accelerate the visualisation on the interface.


.. seealso::

    See :ref:`data preparation examples <Dataprep_example>` how
    to prepare an input data for the platform. You can find more detailed 
    information regarding the filtering and normalisation methods for preparing
    an input from different sources of experiments.



