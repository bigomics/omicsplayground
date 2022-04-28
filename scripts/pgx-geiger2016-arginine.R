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
FILESX = "../libx"
PGX.DIR = "../data"
source("../R/pgx-include.R")
##source("options.R")
FILES

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
##devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)    
library(org.Hs.eg.db)
##library(proteus)    
meta=metadataFile="../ext/arginine/meta.txt"
file=proteinGroupsFile="../ext/arginine/proteinGroups.txt"
## Read the proteinGroups file
prot <- prot.readProteinGroups(
    proteinGroupsFile, meta=metadataFile,
    is.log=TRUE, collapse.gene=TRUE, use.LFQ=FALSE)

if(0) {
    prot <- proteus.readProteinGroups(
        proteinGroupsFile, meta=metadataFile,
        collapse.gene=TRUE, use.LFQ=FALSE)
}

##-------------------------------------------------------------------
## scale/normalize counts for proteomics
##-------------------------------------------------------------------    
## impute missing values
norm.counts <- prot.imputeMissing(
    prot$counts, groups=prot$samples$condition,
    method="group.median", zero.na=TRUE)

## scale/normalize
norm.counts <- prot.normalizeCounts(
    norm.counts, scale=1e6, scaling="global", ## prior.count=0, 
    qnormalize=TRUE)

##hist(log2(1+norm.counts), breaks=200)

##-------------------------------------------------------------------
## create ngs object
##-------------------------------------------------------------------

samples = prot$samples
counts = norm.counts
short.name <- sub(".*_tcell_","",colnames(counts))
rownames(samples)=colnames(counts)=short.name

## Create contrasts 
group.levels <- unique(samples$condition)
group.levels
## 10 contrasts in total
contrasts <- limma::makeContrasts(
    act_vs_notact = (act96h + act72h + act48h + act24h + act12h)/5 - notact,
    act12h_vs_notact = act12h - notact,
    act24h_vs_notact = act24h - notact,
    act48h_vs_notact = act48h - notact,
    act72h_vs_notact = act72h - notact,
    act96h_vs_notact = act96h - notact,
    levels = group.levels)
contrasts
##contrasts = contrasts[,1:3]

ngs <- pgx.createPGX(
    counts = counts,
    samples = samples,
    contrasts = contrasts
)
names(ngs)

gx.methods    = c("trend.limma")
gset.methods  = c("fisher")
extra.methods  = c("meta.go","infer","connectivity")

gx.methods    = c("ttest","ttest.welch","trend.limma","edger.qlf","deseq2.wald")
gset.methods  = c("fisher","gsva","fgsea","spearman","camera","fry")
gset.methods  = c("fisher","ssgsea","fgsea","spearman","camera","fry")
extra.methods  = c("meta.go","infer","drugs","deconv","wordcloud","connectivity")

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

if(0) {
    extra <- c("drugs")
    extra <- c("drugs","connectivity")    
    extra <- c("connectivity")    
    ngs <- pgx.computeExtra(ngs, extra=extra, lib.dir=FILES)
    names(ngs)
    names(ngs$connectivity)
    names(ngs$drugs[[1]])
}


##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------
names(ngs)

rda.file="../data/geiger2016-arginine2.pgx"
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "LC/MS proteomics"
ngs$description = "Proteome profiles of activated  vs resting human naive T cells at different times (Geiger et al., Cell 2016)."
ngs$organism = 'human'

rda.file
ngs.save(ngs, file=rda.file)


if(0) {
    
    source("../R/pgx-include.R")    
    load("../data/geiger2016-arginine-test.pgx")
    names(ngs)

    names(ngs$gx.meta$meta)
    dim(ngs$gx.meta$meta[[1]])
    m1 <- ngs$gx.meta$meta[[1]]
    head(m1[order(-m1$meta.fx),])
    head(sort(rowMeans(ngs$X)))
    tail(sort(rowMeans(ngs$X)))
    
}



##===================================================================
##========================= END OF FILE =============================
##===================================================================
