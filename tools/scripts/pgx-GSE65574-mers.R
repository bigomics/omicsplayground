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

rda.file="../data/GSE65574-mers.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE65574. MERS-CoV Accessory ORFs Play Key Role for Infection and Pathogenesis (Menachery 2017). Human Calu-3 cell transcriptome response to a wild type infectious clone of Middle Eastern Respiratory Syndrome (icMERS) coronavirus and icMERS mutant viruses."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

library(Biobase)
library(GEOquery)
library(limma)
library(hgu133plus2.db)

## load series and platform data from GEO
geo <- getGEO("GSE65574", GSEMatrix=TRUE, getGPL=TRUE)
attr(geo, "names")
X <- exprs(geo[[1]])
head(X)[,1:4]
max(X)

## extract GENE symbol from featureData
colnames(featureData(geo[[1]])@data)
gene.symbol <- as.character(featureData(geo[[1]])@data$GENE_SYMBOL)
##gene.symbol <- gsub("[ ]","",sapply(strsplit(gene.annot,split="//"),"[",2))
gene.symbol[10000 + 1:10]    
jj <- which( !gene.symbol %in% c(NA,"-",""))
X <- X[jj,]
rownames(X) <- gene.symbol[jj]

## Get sample info from title
pdata = pData(geo[[1]])
head(pdata)
tt <- as.character(pdata$title)
head(tt)
tt <- sub("MCL001_","",tt)
tt <- sub("_mRNA","",tt)
tt <- sub("RFP-MERS","RFPMERS",tt)
tt <- sub("d4B-MERS","d4BMERS",tt)
tt <- sub("d3-5-MERS","d35MERS",tt)
tt <- sub("dNSP16-MERS","dNSP16MERS",tt)

head(tt)

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
ngs <- pgx.clusterSamples(ngs, perplexity=20, skipifexists=FALSE, prefix="C")
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
    icMERS_0hr_vs_mock_0hr = icMERS_0hr - mock_0hr,
    icMERS_7hr_vs_mock_7hr = icMERS_7hr - mock_7hr,
    icMERS_12hr_vs_mock_12hr = icMERS_12hr - mock_12hr,
    icMERS_24hr_vs_mock_24hr = icMERS_24hr - mock_24hr,
    
    RFPMERS_0hr_vs_mock_0hr = RFPMERS_0hr - mock_0hr,
    RFPMERS_7hr_vs_mock_7hr = RFPMERS_7hr - mock_7hr,
    RFPMERS_12hr_vs_mock_12hr = RFPMERS_12hr - mock_12hr,
    RFPMERS_24hr_vs_mock_24hr = RFPMERS_24hr - mock_24hr,
    
    dNSP16MERS_0hr_vs_mock_0hr = dNSP16MERS_0hr - mock_0hr,
    dNSP16MERS_7hr_vs_mock_7hr = dNSP16MERS_7hr - mock_7hr,
    dNSP16MERS_12hr_vs_mock_12hr = dNSP16MERS_12hr - mock_12hr,
    dNSP16MERS_24hr_vs_mock_24hr = dNSP16MERS_24hr - mock_24hr,
    
    d4BMERS_0hr_vs_mock_0hr = d4BMERS_0hr - mock_0hr,
    d4BMERS_7hr_vs_mock_7hr = d4BMERS_7hr - mock_7hr,
    d4BMERS_12hr_vs_mock_12hr = d4BMERS_12hr - mock_12hr,
    d4BMERS_24hr_vs_mock_24hr = d4BMERS_24hr - mock_24hr,
    
    d35MERS_0hr_vs_mock_0hr = d35MERS_0hr - mock_0hr,
    d35MERS_7hr_vs_mock_7hr = d35MERS_7hr - mock_7hr,
    d35MERS_12hr_vs_mock_12hr = d35MERS_12hr - mock_12hr,
    d35MERS_24hr_vs_mock_24hr = d35MERS_24hr - mock_24hr,
    
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












