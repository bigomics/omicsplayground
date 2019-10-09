.. _Dataprep:

Data import and precomputation
================================================================================

The data import and precomputation involves preparing the input data through 
filtering, normalising and precomputing statistics for some analyses and 
importing it into the platform. The data cleaning and precomputation is 
performed offline to support real-time interaction by minimizing user interface
latency.

    
Data import
--------------------------------------------------------------------------------
The platform requires the transcriptomics and proteomics data to be in a 
structured format as an input. Users can prepare an input data from
their own FASTQ files, gene counts tables, or from a dataset of interest stored 
in public repositories such as `GEO <https://www.ncbi.nlm.nih.gov/geo/>`__.
Similarly, they can also prepare an input from LC-MS/MS proteomics data.
For all of these cases, the platform comes with the necessary scripts for data 
cleaning and preprocessing under the ``/scripts`` folder.

Users can import their data to the platform by either uploading under the 
:ref:`Home` module of the interface or preparing an input object using scripts.
To upload, the platform requires the tables of counts, samples info, genes info
and contrasts in CSV format. Users can provide their own counts or download the
relevant data from repositories such as `GEO <https://www.ncbi.nlm.nih.gov/geo/>`__. 
On the other hand, an input object can be prepared with more detailed data 
cleaning, filtering, normalisation and preprocessing using scripts. 
The platform contains the required scripts and examples.

.. seealso::

    See :ref:`data preparation examples <Dataprep_example>` how
    to prepare an input data for the platform. You can find more detailed 
    information regarding the filtering and normalisation methods for preparing
    an input from different sources of experiments.
    

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





