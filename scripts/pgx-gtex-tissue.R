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


##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/gtex-tissue.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GTEx. Samples were collected from 54 non-diseased tissue sites across nearly 1000 individuals."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
library(Biobase)
library(limma)
library(rhdf5)

## load series and platform data from GEO
h5.file = "/data/PublicData/archs4data/gtex_matrix.h5"
##h5.file = "../../data/archs4data/gtex_matrix.h5"
h5ls(h5.file)

X0 <- h5read(h5.file, "data/expression")
rownames(X0) <- h5read(h5.file, "meta/genes")
colnames(X0) <- h5read(h5.file, "meta/sample")

tissue <- h5read(h5.file, "meta/tissue")
table(tissue)
X <- c()
s = tissue[1]
tissues <- setdiff(unique(tissue),c(NA,""))
for(s in tissues) {
    ii <- which(tissue == s)
    xx <- X0[,ii]
    if(length(ii) > 10) {
        gx <- log2(100 + xx)
        gx <- head(gx[order(-apply(gx,1,sd)),],2000)
        gx <- gx - rowMeans(gx)
        idx <- cutree(hclust(as.dist(1 - cor(gx))),10)
        xx1 <- tapply(1:ncol(xx),idx,function(j) rowMeans(xx[,j,drop=FALSE]))
        xx <- do.call(cbind, xx1)
    }
    colnames(xx) <- paste0(s,".",1:ncol(xx))
    X <- cbind(X, xx)
}
dim(X)

## Get sample info
xtissue <- sub("[.].*","",colnames(X))
xtissue <- gsub("[ ]","",xtissue)
xtissue[xtissue==""] <- "NA"    
sampleTable <- data.frame( tissue = xtissue)
rownames(sampleTable) <- colnames(X)
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
ngs$counts <- X  ## treat as counts
ngs$samples <- data.frame(sampleTable)
ngs$genes = genes

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters 
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, perplexity=30, skipifexists=FALSE, prefix="C")
head(ngs$samples)
table(ngs$samples$cluster)

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------    
head(ngs$samples)
ngs$samples$group <- ngs$samples$tissue
levels = unique(ngs$samples$group)
levels

ct1 <- makeClusterContrasts(
    ngs$samples$tissue, min.freq=0.01, full=FALSE, by.sample=FALSE ) 

ct2 <- makeClusterContrasts(
    ngs$samples$tissue, min.freq=0.01, full=TRUE, by.sample=FALSE ) 

contr.matrix <- cbind( ct1, ct2)
dim(contr.matrix)


##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------
rda.file
ngs$timings <- c()

GENE.METHODS=c("ttest.welch","trend.limma","edger.qlf")
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
extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
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












