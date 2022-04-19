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
FILES

##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/GSE88808-prostate.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "mRNA array"
ngs$description = "GSE88808 data set. Gleason-score matched tumor and adjacent normal samples were collected to compare gene expression differences in early-onset versus late-onset prostate cancer patients (Ding, PLOS Genet 2016)."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

library(Biobase)
library(GEOquery)
library(data.table)

##--------------------------------------------------------------
## Read SC counts
##--------------------------------------------------------------
## load series and platform data from GEO
geo <- getGEO("GSE88808", GSEMatrix =TRUE, AnnotGPL=FALSE)
length(geo)
attr(geo, "names")
counts = 2**exprs(geo[[1]])
head(counts)[,1:10]
summary(colSums(counts))

##--------------------------------------------------------------
## gene annotation
##--------------------------------------------------------------
require(org.Hs.eg.db)
GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
head(GENE.TITLE)
gene_title <- GENE.TITLE[rownames(counts)]
genes = data.frame( gene_name=rownames(counts), gene_title=gene_title)
rownames(genes) <- rownames(counts)

##--------------------------------------------------------------
## Prepare sample table
##--------------------------------------------------------------
##geo <- getGEO("GSE99795", GSEMatrix =TRUE, AnnotGPL=TRUE)
pdata = pData(geo[[1]])
head(pdata)
sampleTable <- pdata[,grep(":ch1$",colnames(pdata))]
colnames(sampleTable) <- sub(":ch1","",colnames(sampleTable))
head(sampleTable)

age <- as.numeric(sampleTable$age)
summary(age)
sampleTable$age.group <- c("young","old")[ 1 + 1*(age >= median(age)) ]
head(sampleTable)

##-------------------------------------------------------------------
## Now create an PGX object
##-------------------------------------------------------------------
ngs$counts <- round(counts)
ngs$samples <- sampleTable
ngs$genes = data.frame(genes)
rownames(ngs$genes) <- genes$gene_name

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE,
                          perplexity=30, kclust=2)
head(ngs$samples)

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
##sampleTable$group <- paste0(sampleTable$Stage,"_",sampleTable$tissue)
ngs$samples$group <- paste0(ngs$samples$age.group,"_",ngs$samples$tissue)
table(ngs$samples$group)
levels = unique(ngs$samples$group)
levels

contr.matrix <- limma::makeContrasts(
    old_tumor_vs_old_normal = old_tumor - old_normal,
    young_tumor_vs_young_normal = young_tumor - young_normal,
    old_normal_vs_young_normal = old_normal - young_normal,
    old_tumor_vs_young_tumor = old_tumor - young_tumor,
    levels = levels)

dim(contr.matrix)
contr.matrix

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------

GENE.METHODS=c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
GENESET.METHODS = c("fisher","gsva","camera","fgsea")

MAX.GENES = 8000
MAX.GENESETS = 8000

## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features = MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features = MAX.GENESETS,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("drugs-combo")
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
