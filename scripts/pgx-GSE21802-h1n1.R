##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
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

rda.file="../data/GSE21802-h1n1.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE21802. Host adaptive immunity deficiency in severe pandemic influenza (Bermejo-Martin, Crit Care 2010)."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
library(GEOquery)
library(hgu133plus2.db)

## load series and platform data from GEO
geo <- getGEO("GSE21802", GSEMatrix=TRUE, getGPL=TRUE)
attr(geo, "names")
X <- as.matrix(exprs(geo[[1]]))
head(X)[,1:4]
max(X)
min(X)
X <- log2(X)

## extract GENE symbol from featureData
colnames(featureData(geo[[1]])@data)
gene.symbol <- as.character(featureData(geo[[1]])@data$Symbol)
##gene.symbol <- gsub("[ ]","",sapply(strsplit(gene.annot,split="//"),"[",2))
gene.symbol[10000 + 1:10]    
jj <- which( !gene.symbol %in% c(NA,"-",""))
X <- X[jj,]
rownames(X) <- gene.symbol[jj]

X <- X[order(-apply(X,1,sd)),]
X <- X[!duplicated(rownames(X)),]

## Get sample info from CH1
pdata = pData(geo[[1]])
head(pdata)
sampleTable <- pdata[,grep(":ch1",colnames(pdata))]
head(sampleTable)
##colnames(sampleTable) <- c("code","infected","time","replicate")
colnames(sampleTable) <- sub(":ch1$","",colnames(sampleTable))
head(sampleTable)

colnames(sampleTable) <- gsub("[ ]","_",colnames(sampleTable))
colnames(sampleTable) <- sub("mechanical_ventilation","ventilation",colnames(sampleTable))
head(sampleTable)

##sampleTable$time <- sub(" hours post infection","H",sampleTable$time)
##ngs$samples$induce_IFNB <- sub("YES","IFNB",ngs$samples$induce_IFNB)
sampleTable$virus_strain  <- sub("Pandemic influenza \\(A/H1N1 new subtype\\)","H1N1",sampleTable$virus_strain)
sampleTable$disease_phase <- sub(" period","",sampleTable$disease_phase)
sampleTable$tissue <- NULL
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

## take out duplicated
jj <- order(-apply(X,1,sd))
X <- X[jj,]
genes <- genes[jj,]    
jj <- which(!duplicated(genes$gene_name) & !is.na(genes$gene_name) & !is.na(genes$chr))
X <- X[jj,]
genes <- genes[jj,]
rownames(X) <- rownames(genes) <- genes$gene_name


##-------------------------------------------------------------------
## Now create an DGEList object  (see tximport Vignette)
##-------------------------------------------------------------------
library(limma)
##X <- limma::normalizeQuantiles(X)
ngs$counts <- 2**X  ## treat as counts
ngs$samples <- data.frame(sampleTable)
ngs$genes = genes

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
dim(ngs$counts)
ngs <- pgx.clusterSamples(ngs, perplexity=10, skipifexists=FALSE, prefix="C")
head(ngs$samples)


##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
##load(file=rda.file, verbose=1)

head(ngs$samples)
grp <- paste(ngs$samples$virus_strain, ngs$samples$disease_phase, ngs$samples$ventilation, sep="_")
grp <- sub("none_control_NA","control",grp)
ngs$samples$group <- grp
levels = unique(ngs$samples$group)
levels

contr.matrix <- limma::makeContrasts(                
    H1N1earlyVTLno_vs_control = H1N1_early_no - control,
    H1N1lateVTLno_vs_control = H1N1_late_no - control,
    H1N1earlyVTLyes_vs_control = H1N1_early_yes - control,
    H1N1lateVTLyes_vs_control = H1N1_late_yes - control,
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












