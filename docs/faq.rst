.. _FAQ:


Frequently Asked Questions (FAQ)
================================================================================


Can I use LIMMA, EdgeR and DESeq2 for my proteomics data?
----------------------------------------------------------

EdgeR and DESeq2 are statistical methods based on the negative
binomial model (an overdispersed Poisson model). Poisson and negative
binomial models, naturally account for heteroscedasticity and zero
values in the data. Also LIMMA is widely used in RNA-seq analysis, and
is an emperical Bayesian method based on moderated t-test
statistics. While these methods were originally conceived for the for
differential expression analysis in RNA-seq data, there is
increasingly more acceptance that these methods can also be used for
proteomics data.

Langley and Mayer (2015) assessed seven methods for differential
expression analysis in proteomics, and showed that DESeq to outperform
the more commonly used t-test. Kammers et al. (2015) show that LIMMA
shows better results than the t-test. Branson and Fretais (2016) also
highlighted the statistical analogy of proteomics label-free spectral
count quantification with count data from RNA sequencing and propose
to use edgeR, DESEq and baySeq for proteomics data. Gregori et
al. (2013) uses EdgeR and Poisson based methods for their LC/MS-MS
data. They also created the 'msmsTests' R package. Medo et al (2019)
use EdgeR for their study on missing values in the differential
analysis of proteomic and phosphoproteomics data. Gatto proposes both
edgeR and LIMMA based methods for proteomics statistical
analysis. Chen et al. (2020) reviewed bioinformatics methods for
proteomics data analysis and reaffirmed that LIMMA can achieve more
robust and accurate results than the traditional t-test.

Notice that proper scaling/normalization of quantitative proteomics
data is important before using RNA-seq based methods as sometimes the
proteomics intensity values may be far greater than those found in
RNA-seq. Therefore Omics Playground scales proteomics data
automatically to 'counts per million' (CPM). We have also seen that
batch correction (e.g. using ComBat) may improve the downstream
statistical analysis in proteomics.


- Langley SR, Mayr M. "Comparative analysis of statistical methods
used for detecting differential expression in label-free mass
spectrometry proteomics." J Proteomics. 2015 Nov 3;129:83-92.
- Kammers et al. "Detecting significant changes in protein abundance."
EuPA Open Proteomics Volume 7, June 2015.
- Branson OE, Freitas MA. A multi-model statistical approach for
proteomic spectral count quantitation. J Proteomics. 2016.
- Gregori et al. "msmsTests-package: LC-MS/MS Differential Expression
Tests". R package.
- Gregori et al. "An Effect Size Filter Improves the Reproducibility
in Spectral Counting-based Comparative Proteomics." Journal of
Proteomics, 2013.
- Medo, M, Aebersold, DM and Medová, M "ProtRank: bypassing the
imputation of missing values in differential expression analysis of
proteomic data." BMC Bioinformatics 20, 563 (2019).
- Gatto L. "Bioconductor tools for mass spectrometry and proteomics",
https://lgatto.github.io/bioc-ms-prot/lab.html#8_statistical_analysis
- Chen et al. "Bioinformatics Methods for Mass Spectrometry-Based
Proteomics Data Analysis." Int J Mol Sci. 2020;21(8):2873, 2020.


How are duplicated gene/protein names handled in the counts file?
----------------------------------------------------------

Duplicated row identifiers (genes/proteins with same name) are handled
by summing up their linear intensities/counts. If the data was in
logarithm, it will be (hopefully) automatically detected and
exponentiated. The rational of summing up the counts (or linear
intensities in proteomics) is that we don't differentiate between
possible gene/protein isoforms and sum them up as a group. If you want
to retain the isoforms, you may keep the names as GENE.1 and GENE.2
but you must turn off any gene filter. However as currently such
gene/protein variants are not recognized in the gene sets, this will
result in wrong enrichment test.


How are "missing values" handled in the counts file?
----------------------------------------------------------

At the moment missing values are imputed to 0 (zero). The Omics
Playground uses the function `pgx.createPGX()` in the file
pgx-compute.R. Note that only real NA (or empty) values in the counts
file are supposed to be "missing". Zero values are assumed to be real
zeros. If zero-value imputaton is not what you want, and you want to
impute your "missing" values differently, (currently) you must do that
manually before uploading the CSV file.
