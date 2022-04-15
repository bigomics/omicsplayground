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

rda.file="../data/GSE100509-merscov.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE100509. Primary human microvascular endothelial cell transcriptome response to wild type MERS-CoV (icMERS-CoV EMC2012), MMVE003. The purpose is to obtain samples for mRNA analysis in primary human microvascular endothelial cells infected with wild type MERS-coronavirus (MERS-CoV) (icMERS)."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
## READ/PARSE DATA

library(Biobase)
library(GEOquery)
library(limma)
library(hgu133plus2.db)

## load series and platform data from GEO
geo <- getGEO("GSE100509", GSEMatrix=TRUE, getGPL=TRUE)
attr(geo, "names")
X <- exprs(geo[[1]])
head(X)[,1:4]

## extract GENE symbol from featureData
colnames(featureData(geo[[1]])@data)
gene.symbol <- as.character(featureData(geo[[1]])@data$GENE_SYMBOL)
##gene.symbol <- gsub("[ ]","",sapply(strsplit(gene.annot,split="//"),"[",2))
gene.symbol[10000 + 1:10]    
jj <- which( !gene.symbol %in% c(NA,"-",""))
X <- X[jj,]
rownames(X) <- gene.symbol[jj]

## Get sample info
pdata = pData(geo[[1]])
head(pdata)
tt <- as.character(pdata$title)
head(tt)
tt <- sub("_RNA","",tt)
tt <- sub("HRS","H",tt)
tt <- sub("MMVE003_","",tt)
sampleTable <- do.call(rbind,strsplit(tt,split="_"))
head(sampleTable)
colnames(sampleTable) <- c("infected","time","replicate")
colnames(X) <- rownames(sampleTable) <- tt
head(sampleTable)

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
X <- limma::normalizeQuantiles(X)
ngs$counts <- 2**X  ## treat as counts
ngs$samples <- data.frame(sampleTable)
ngs$genes = genes

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, perplexity=30, skipifexists=FALSE)
head(ngs$samples)

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
head(ngs$samples)
grp <- paste(ngs$samples$infected,ngs$samples$time,sep="_")
ngs$samples$group <- grp
levels = unique(ngs$samples$group)
levels

contr.matrix <- limma::makeContrasts(
    icMERS_0H_vs_MOCK_0H = icMERS_0H - MOCK_0H,
    icMERS_12H_vs_MOCK_12H = icMERS_12H - MOCK_12H,
    icMERS_24H_vs_MOCK_24H = icMERS_24H - MOCK_24H,
    icMERS_48H_vs_MOCK_48H = icMERS_48H - MOCK_48H,        
    levels = levels)
contr.matrix

rda.file
ngs$timings <- c()

GENE.METHODS=c("ttest","ttest.welch","ttest.rank",
                   "voom.limma","trend.limma","notrend.limma",
                   "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
GENESET.METHODS = c("fisher","gsva","ssgsea","spearman",
                    "camera", "fry","fgsea") ## no GSEA, too slow...
GENE.METHODS=c("trend.limma","edger.qlf","deseq2.wald")
GENESET.METHODS = c("fisher","gsva","fgsea") ## no GSEA, too slow...

MAX.GENES = 20000
MAX.GENESETS = 5000

## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features = MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features = MAX.GENESETS,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("connectivity")
extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

names(ngs)
ngs$timings

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------
rda.file
ngs.save(ngs, file=rda.file)












