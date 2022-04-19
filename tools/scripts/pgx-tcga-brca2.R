##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
## Author: BigOmics Analytics (IK)
## Date:   2020
## 

RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")
source("../R/pgx-tcga.R")
##source("options.R")

##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

##rda.file = "../data/tcga-brca_pub2.pgx"
rda.file = "../data/tcga-brca_pub.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "RNA-seq"
##ngs$datatype = "multi-omics"
ngs$description = "TCGA breast cancer data set. Gene expression from 526 patients annotated with PAM50 classification. Data from cBioPortal."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

load("../../omxp5/omxdata/tcga/brca_tcga_pub-omx.rda",verbose=TRUE)
gx <- omx$level[[1]]$mat$gx
px <- omx$level[[1]]$mat$px
samples <- omx$pheno

## Create contrasts 
ct <- makeDirectContrasts(samples, ref=rep(0,ncol(samples)))
contr.matrix <- ct$exp.matrix

dtype='gx'
dtype='px'
##for(dtype in c("gx","px")) {
for(dtype in c("px")) {
    
    counts <- 2**(gx - min(gx,na.rm=TRUE)) - 1
    if(dtype == 'px') counts <- 2**(px - min(px,na.rm=TRUE)) - 1
    kk <- sort(intersect(rownames(samples),colnames(counts)))
    
    ngs <- pgx.createPGX(
        counts = counts[,kk],
        samples = samples[kk,],
        contrasts = contr.matrix[kk,],
        is.logx = FALSE
    )
    names(ngs)
    
    gx.methods    = c("ttest","ttest.welch","trend.limma","edger.qlf","deseq2.wald")
    gset.methods  = c("fisher","gsva","fgsea","spearman","camera","fry")
    extra.methods  = c("meta.go","infer","drugs","deconv","wordcloud","connectivity")
    
    gx.methods    = c("trend.limma")
    gset.methods  = c("spearman")
    ##gx.methods    = c("trend.limma","edger.qlf","deseq2.wald")
    ## gset.methods  = c("fisher","gsva","fgsea","spearman")
    extra.methods  = c("meta.go","infer","drugs","wordcloud","connectivity")
    
    ngs <- pgx.computePGX(
        ngs,
        max.genes = 20000,
        max.genesets = 5000, 
        gx.methods = gx.methods,
        gset.methods = gset.methods,
        extra.methods = extra.methods,
        use.design = TRUE,      ## no.design+prune are combined 
        prune.samples = FALSE,  ##
        do.cluster = TRUE,                
        progress = NULL,
        lib.dir = FILES 
    )
    
    ##-------------------------------------------------------------------
    ## save PGX object
    ##-------------------------------------------------------------------
    
    rda.file = paste0("../data/tcga-brca2-",dtype,".pgx")
    ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
    ngs$date = date()
    ngs$datatype = "proteomics"
    ngs$datatype = "mRNA"
    ngs$description = "TCGA BRCA gene expression (from OMX object)"
    if(dtype=='px') ngs$description = "TCGA BRCA proteomics (from OMX object)"    
    ngs$organism = 'human'
    
    rda.file
    ngs.save(ngs, file=rda.file)

}

    
##===================================================================
##========================= END OF FILE =============================
##===================================================================


