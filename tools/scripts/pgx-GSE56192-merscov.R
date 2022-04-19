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

rda.file="../data/GSE56192-mers.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE56192. Transcriptomic Analysis Of The Novel Middle East Respiratory Syndrome Coronavirus (MERS-CoV)."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
library(Biobase)
library(GEOquery)
library(limma)
library(hgu133plus2.db)

## load series and platform data from GEO
geo <- getGEO("GSE56192", GSEMatrix=TRUE, getGPL=TRUE)
attr(geo, "names")

df <- fread("../data/GSE56192_GeneLevel_Raw_data.csv")
head(df)[,1:4]
X <- as.matrix(df[,3:ncol(df)])
rownames(X) <- df$gene_symbol
head(X)[,1:4]
max(X)    

## take out duplicated
X1 <- tapply(1:nrow(X), rownames(X), function(i) colMeans(X[i,,drop=FALSE],na.rm=TRUE))
X <- do.call(rbind, X1)
head(X)[,1:4]
remove(X1)

## Get sample info from title
pdata = pData(geo[[1]])
head(pdata)
title <- as.character(unlist(sapply(geo, function(d) pData(d)$title)))
gsm <- unlist(sapply(geo, function(d) pData(d)$geo_accession))
names(title) <- gsm

## clean up titles.... :(
head(title)
tt <- sub("lowMOI","_lowMOI",title)
tt <- sub("HighMOI","_HighMOI",tt)
tt <- sub("hrs$","hr_0",tt)
tt <- sub("MOCK_MRC5_([1-3])","MOCK_MRC5_mock_0hr_\\1",tt)
tt <- sub("MRC5_","",tt)
head(tt)

sampleTable <- do.call(rbind,strsplit(tt,split="_"))
head(sampleTable)
colnames(sampleTable) <- c("infected","treatment","time","replicate")
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
rownames(genes) <- rownames(X)
head(genes)

##-------------------------------------------------------------------
## Put in NGS/PGX object
##-------------------------------------------------------------------
library(limma)
##X <- limma::normalizeQuantiles(X)
ngs$counts <- X  ## treat as counts
ngs$samples <- data.frame(sampleTable)
ngs$genes = genes

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, perplexity=30, skipifexists=FALSE, prefix="C")
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
                
    MERS_24hr_vs_MOCK_24hr = MERS_24hr - MOCK_24hr,
    MERS_48hr_vs_MOCK_48hr = MERS_48hr - MOCK_48hr,
    
    SARS_24hr_vs_MOCK_24hr = SARS_24hr - MOCK_24hr,
    SARS_48hr_vs_MOCK_48hr = SARS_48hr - MOCK_48hr,
    
    MOCK_24hr_vs_MOCK_0hr = MOCK_24hr - MOCK_0hr,
    MOCK_48hr_vs_MOCK_0hr = MOCK_48hr - MOCK_0hr,
    
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












