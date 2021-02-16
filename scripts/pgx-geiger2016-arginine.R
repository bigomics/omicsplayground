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
##source("options.R")
FILES

##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/geiger2016-arginine-test.pgx"
##rda.file="../data/geiger2016-arginine.pgx"
rda.file

## load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "LC/MS proteomics"
ngs$description = "Proteome profiles of activated  vs resting human naive T cells at different times (Geiger et al., Cell 2016)."

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
##ngs$samples = prot$metadata
ngs$samples = prot$samples
##colnames(ngs$samples) <- sub("condition","group",colnames(ngs$samples))
##ngs$counts = data$counts
ngs$counts = norm.counts
colnames(ngs$counts)==ngs$samples$sample
short.name <- sub(".*_tcell_","",colnames(ngs$counts))
rownames(ngs$samples)=colnames(ngs$counts)=short.name

ngs$X <- log2(ngs$counts + 1)

## relevel factors??
##ngs$samples$group <- relevelFactorFirst(ngs$samples$group)
ngs$samples$condition <- relevelFactorFirst(ngs$samples$condition)
ngs$samples$activated <- relevelFactorFirst(ngs$samples$activated)
ngs$samples$time <- relevelFactorFirst(ngs$samples$time)   

require(org.Hs.eg.db)
GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
## ngs$genes = data.frame( gene_name = prot$gene,
##                        gene_alias = prot$gene.names,
##                        gene_title = GENE.TITLE[prot$gene] )
ngs$genes <- prot$genes

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=3)
ngs <- pgx.clusterSamples2(ngs, perplexity=3)
head(ngs$samples)
table(ngs$samples$cluster)

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
group.levels <- levels(ngs$samples$condition)
group.levels
## 10 contrasts in total
contr.matrix <- makeContrasts(
    act_vs_notact = (act96h + act72h + act48h + act24h + act12h)/5 - notact,
    act12h_vs_notact = act12h - notact,
    act24h_vs_notact = act24h - notact,
    act48h_vs_notact = act48h - notact,
    act72h_vs_notact = act72h - notact,
    act96h_vs_notact = act96h - notact,
    act48h_vs_act12h = act48h - act12h,
    act72h_vs_act12h = act72h - act12h,
    act72h_vs_act48h = act72h - act48h,
    act96h_vs_act48h = act96h - act48h,
    act96h_vs_act72h = act96h - act72h,
    levels = group.levels)
contr.matrix
##contr.matrix = contr.matrix[,1:3]

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------
source("../R/pgx-include.R")

ngs$timings <- c()

GENE.METHODS=c("trend.limma")
GENESET.METHODS = c("fisher","fgsea")

GENE.METHODS=c("ttest","ttest.welch", ## "ttest.rank",
                   "voom.limma","trend.limma","notrend.limma",
                   "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
GENESET.METHODS = c("fisher","gsva","ssgsea","spearman",
                    "camera", "fry","fgsea") ## no GSEA, too slow...

MAX.GENES = 8000
MAX.GENESETS = 8000

## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features = MAX.GENES,
    test.methods = GENE.METHODS)
names(ngs)
head(ngs$gx.meta$meta[[1]])

ngs <- compute.testGenesets (
    ngs, max.features=MAX.GENES,
    test.methods = GENESET.METHODS,
    remove.outputs = FALSE,
    lib.dir=FILES)

extra <- c("drugs-combo")
extra <- c("connectivity")
extra <- c("meta.go","infer","deconv","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

if(0) {
    sigdb = NULL
    sigdb = c("../libx/sigdb-lincs.h5","../libx/sigdb-virome.h5")
    ngs <- compute.extra(ngs, extra, lib.dir=FILES, sigdb=sigdb) 
    names(ngs$connectivity)
}

names(ngs)
names(ngs$gset.meta)
ngs$timings

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------
rda.file
ngs.save(ngs, file=rda.file)

if(0) {

    load("../data/geiger2016-arginine-test.pgx")
    names(ngs$gx.meta$meta)
    dim(ngs$gx.meta$meta[[1]])
    m1 <- ngs$gx.meta$meta[[1]]
    head(m1[order(-m1$meta.fx),])
    head(sort(ngs$gx.meta$meta[[1]]$meta.avg),50)
    tail(sort(ngs$gx.meta$meta[[1]]$meta.avg),50)
    head(sort(rowMeans(ngs$X)))
    tail(sort(rowMeans(ngs$X)))

    s0 <- names(which(ngs$model.parameters$exp.matrix[,1]==-1))
    s1 <- names(which(ngs$model.parameters$exp.matrix[,1]==1))
    mean(ngs$X["TNIP3",s0])
    mean(ngs$X["TNIP3",s1])
    
    
    
}



##===================================================================
##========================= END OF FILE =============================
##===================================================================
