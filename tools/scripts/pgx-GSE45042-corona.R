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

rda.file="../data/GSE45042-corona.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "RNA-seq"
ngs$description = "Cell host-response to infection with novel human coronavirus EMC predict potential antivirals and important differences with SARS-coronavirus."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

library(GEOquery)
library(limma)
library(hgu133plus2.db)

## load series and platform data from GEO
geo <- getGEO("GSE45042", GSEMatrix=TRUE, getGPL=TRUE)
attr(geo, "names")
X <- exprs(geo[[1]])
head(X)

## extract GENE symbol from featureData
gene.symbol <- as.character(featureData(geo[[1]])@data$GENE_SYMBOL)
dim(X)
length(gene.symbol)
## gene.symbol <- gsub("[ ]","",sapply(strsplit(gene.annot,split="//"),"[",2))
## gene.symbol <- alias2hugo(gene.symbol)
gene.symbol[10000 + 1:10]    
jj <- which( !gene.symbol %in% c(NA,"-",""))
X <- X[jj,]
rownames(X) <- gene.symbol[jj]

## Get sample info
pdata = pData(geo[[1]])
head(pdata)
tt <- as.character(pdata$title)
sampleTable <- do.call(rbind,strsplit(tt,split="_"))
head(sampleTable)
colnames(sampleTable) <- c("code","treatment","time","replicate")
colnames(X) <- rownames(sampleTable) <- tt

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
ngs <- pgx.clusterSamples(ngs, perplexity=NULL, skipifexists=FALSE, prefix="C")
head(ngs$samples)


##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------

##load(file=rda.file, verbose=1)
head(ngs$samples)
grp <- paste(ngs$samples$treatment,ngs$samples$time,sep="_")
ngs$samples$group <- grp
levels = unique(ngs$samples$group)
levels
    
contr.matrix <- limma::makeContrasts(
    EMC_0h_vs_Mock_0h = EMC_0hr - Mock_0hr,
    EMC_3h_vs_Mock_3h = EMC_3hr - Mock_3hr,
    EMC_7h_vs_Mock_7h = EMC_7hr - Mock_7hr,
    EMC_12h_vs_Mock_12h = EMC_12hr - Mock_12hr,
    EMC_18h_vs_Mock_24h = EMC_18hr - Mock_24hr,
    EMC_24h_vs_Mock_24h = EMC_24hr - Mock_24hr,
    levels = levels)
contr.matrix
    
contr.matrix
##contr.matrix = contr.matrix[,1:3]

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------
rda.file
ngs$timings <- c()

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
