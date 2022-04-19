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

##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/GSE53784-wnvjev.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "RNA-seq"
ngs$description = "GSE53784 (Clarke et al., MBio 2014). Gene expression in the brain following WNV or JEV infection. WNV- or JEV-infected (N=3) vs. mock-infected (N=3) mouse brain."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
library(Biobase)
library(GEOquery)
library(limma)
library(hgu133plus2.db)

## load series and platform data from GEO
geo <- getGEO("GSE53784", GSEMatrix=TRUE, getGPL=TRUE)
attr(geo, "names")
X <- exprs(geo[[1]])
head(X)

## extract GENE symbol from featureData
colnames(featureData(geo[[1]])@data)
gene.annot <- featureData(geo[[1]])@data$gene
gene.symbol <- gsub("[ ]","",sapply(strsplit(gene.annot,split="//"),"[",2))
gene.symbol[10000 + 1:10]    
jj <- which( !gene.symbol %in% c(NA,"-",""))
X <- X[jj,]
rownames(X) <- gene.symbol[jj]

## Get sample info
pdata = pData(geo[[1]])
head(pdata)
tt <- as.character(pdata$title)    
treatment <- sub("_.*","",tt)
replicate <- sub(".*_","",tt)
sampleTable <- data.frame(sample=tt, treatment=treatment)
colnames(X) <- rownames(sampleTable) <- tt

## conform tables
sample.names <- as.character(sampleTable$sample)
rownames(sampleTable) = colnames(X) = sample.names

##-------------------------------------------------------------------
## gene annotation
##-------------------------------------------------------------------
require(org.Mm.eg.db)
GENE.TITLE = unlist(as.list(org.Mm.egGENENAME))
gene.symbol = unlist(as.list(org.Mm.egSYMBOL))
names(GENE.TITLE) = gene.symbol
head(GENE.TITLE)
gene_title <- GENE.TITLE[rownames(X)]

## get chromosome locations
chrloc = as.list(org.Mm.egCHRLOC)
names(chrloc) = gene.symbol
chrloc <- chrloc[rownames(X)]
loc <- sapply(chrloc, "[", 1)
chrom <- sapply(chrloc, function(s) names(s)[1])
loc[sapply(loc,is.null)] <- NA
chrom[sapply(chrom,is.null)] <- NA
chrom <- as.vector(unlist(chrom))
loc   <- as.vector(unlist(loc))

genes = data.frame( gene_name=rownames(X),
                   gene_title=gene_title,
                   chr=chrom, pos=loc)
##genes = apply(genes,2,as.character)
head(genes)

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
ngs$samples <- sampleTable
ngs$genes = genes

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, perplexity=2, skipifexists=FALSE, prefix="C")
head(ngs$samples)


##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
head(ngs$samples)
ngs$samples$group <- ngs$samples$treatment
levels = unique(ngs$samples$group)
levels

contr.matrix <- limma::makeContrasts(
    JEV_vs_MOCK = JEV - MOCK,
    WNV_vs_Mock = WNV - Mock,
    levels = levels)
contr.matrix

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------

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
