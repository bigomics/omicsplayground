.. _Dataprep:

Data import and precomputation
================================================================================

The data import and precomputation involve preparing the input data through 
filtering, normalising and precomputing statistics for some analyses and 
importing it into the platform. The data cleaning and precomputation is 
performed offline to support real-time interaction by minimizing user interface
latency.

    
Data import
--------------------------------------------------------------------------------
Users can import their transcriptomics or proteomics data to the platform by 
either uploading the data through
the interface or preparing an input object using scripts.
For uploading, the platform requires the counts, samples information, genes 
information and contrasts tables in CSV format. 
On the other hand, an input object can be prepared using scripts from different 
types and formats of data, including counts and FASTQ.
With scripts it is also possible to do more detailed data 
cleaning, filtering, normalisation and preprocessing. 
The platform contains the required example cases for the preparation of input 
objects under the ``scripts/`` folder.


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





