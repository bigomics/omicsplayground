##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
##
##
##

RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")
##source("options.R")


##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/GSE5972-sars.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE5972. Gene expression profiling of patients with severe acute respiratory syndrome (SARS) [...] from onset of symptoms to discharge from hospital or death. We modeled gene expression in 60 datasets from 40 unique patients of varied clinical evolution with emphasis on correlating innate and adaptive immune responses with discrete clinical phases of the disease."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

## ##############################################################
##   Differential expression analysis with limma
library(Biobase)
library(GEOquery)

## load series and platform data from GEO
geo <- getGEO("GSE5972", GSEMatrix=TRUE, getGPL=TRUE)
attr(geo, "names")
X <- exprs(geo[[1]])
dim(X)
head(X)[,1:4]

## extract GENE symbol from featureData
colnames(featureData(geo[[1]])@data)
gene.symbol <- as.character(featureData(geo[[1]])@data$"GENE_SYMBOL")
gene.symbol <- gsub("[ ]","",sapply(strsplit(gene.symbol,split="///"),"[",1))
gene.symbol <- alias2hugo(gene.symbol)
gene.symbol[100 + 1:10]    
jj <- which(!gene.symbol %in% c(NA,"-",""))
X <- X[jj,]
rownames(X) <- gene.symbol[jj]

## Get sample info
pdata = pData(geo[[1]])
head(pdata)
ch2 <- as.character(pdata$characteristics_ch2)
sampleTable <- eset.parseCharacteristicsInfo(ch2, split=",")
rownames(sampleTable) <- rownames(pdata)
head(sampleTable)

## simplify
sampleTable$Status <- gsub(" nadir| control","",sampleTable$Status)
sampleTable$Status <- gsub("[-]","_",sampleTable$Status)

##-------------------------------------------------------------------
## gene annotation
##-------------------------------------------------------------------
require(org.Hs.eg.db)
GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
head(GENE.TITLE)
gene_title <- GENE.TITLE[rownames(X)]

## get chromosome locations
chrloc = sapply(as.list(org.Hs.egMAP),"[",1)
names(chrloc) = gene.symbol
chrloc <- chrloc[rownames(X)]

genes = data.frame( gene_name=rownames(X),
                   gene_title=gene_title,
                   chr=chrloc)
##genes = apply(genes,2,as.character)
head(genes)

## take out duplicated
jj <- order(-apply(X,1,sd))
X <- X[jj,]
genes <- genes[jj,]    
jj <- which(!duplicated(genes$gene_name) & !is.na(genes$gene_name))
X <- X[jj,]
genes <- genes[jj,]
rownames(X) <- rownames(genes) <- genes$gene_name

##-------------------------------------------------------------------
## Now create an DGEList object  (see tximport Vignette)
##-------------------------------------------------------------------
library(limma)
ngs$counts <- 2**(8 + X)  ## treat as counts
ngs$samples <- data.frame(sampleTable)
ngs$genes = genes

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, perplexity=NULL, skipifexists=FALSE, prefix="C")
head(ngs$samples)


##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------

##load(file=rda.file, verbose=1)    
head(ngs$samples)
ngs$samples$group <- ngs$samples$Status
levels = unique(ngs$samples$group)
levels

contr.matrix <- limma::makeContrasts(
    prePO2_vs_healthy = pre_pO2 - healthy,
    postPO2_vs_healthy = post_pO2 - healthy,
    postPO2_vs_prePO2 = post_pO2 - pre_pO2,
    convalescent_vs_healthy = convalescent - healthy,
    convalescent_vs_postPO2 = convalescent - post_pO2,
    convalescent_vs_prePO2 = convalescent - pre_pO2,        
    levels = levels)
contr.matrix

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------

rda.file
ngs$timings <- c()

GENE.METHODS = c("trend.limma","edger.qlf","deseq2.wald")
GENESET.METHODS = c("fisher","gsva","fgsea") ## no GSEA, too slow...

MAX.GENES = 20000
MAX.GENESETS = 5000

## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features = MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features=MAX.GENESETS,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("connectivity")
extra <- c("meta.go","infer","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

names(ngs)
ngs$timings

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------

rda.file
ngs.save(ngs, file=rda.file)

##===================================================================
##========================= END OF FILE =============================
##===================================================================












