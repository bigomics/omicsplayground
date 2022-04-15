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

rda.file="../data/GSE17400-sarscov.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE17400. Dynamic Innate Immune Responses of Human Bronchial Epithelial Cells against SARS-CoV and DOHV infection (Yoshikawa, PLoS One 2010)."


## load series and platform data from GEO
library(GEOquery)
library(hgu133plus2.db)

geo <- getGEO("GSE17400", GSEMatrix=TRUE, getGPL=TRUE)
attr(geo, "names")
X <- exprs(geo[[1]])
head(X)[,1:4]
max(X)
min(X)
##X <- log2(X)

## extract GENE symbol from featureData
colnames(featureData(geo[[1]])@data)
gene.symbol <- as.character(featureData(geo[[1]])@data$"Gene Symbol")
gene.symbol <- gsub("[ ]","",sapply(strsplit(gene.symbol,split="///"),"[",1))
gene.symbol[10000 + 1:10]    
jj <- which( !gene.symbol %in% c(NA,"-",""))
X <- X[jj,]
rownames(X) <- gene.symbol[jj]

## Get sample info from CH1
pdata = pData(geo[[1]])
head(pdata)
sampleTable <- pdata[,grep(":ch1",colnames(pdata))]
head(sampleTable)
##colnames(sampleTable) <- c("code","infected","time","replicate")
colnames(sampleTable) <- sub(":ch1$","",colnames(sampleTable))
head(sampleTable)
colnames(sampleTable) <- c("cell_line","infection","time")
sampleTable$time <- sub(" hours post infection","H",sampleTable$time)
##ngs$samples$induce_IFNB <- sub("YES","IFNB",ngs$samples$induce_IFNB)
sampleTable$infection <- sub("-infected","",sampleTable$infection)
sampleTable$infection <- sub("SARS-CoV","SARSCoV",sampleTable$infection)
sampleTable$infection <- sub(" \\(MOI=0.1\\)","",sampleTable$infection)
sampleTable$cell_line <- NULL
head(sampleTable)

##colnames(X) <- rownames(sampleTable) <- tt


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
ngs <- pgx.clusterSamples(ngs, perplexity=8, skipifexists=FALSE, prefix="C")
head(ngs$samples)

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------

head(ngs$samples)
grp <- paste(ngs$samples$infection,ngs$samples$time,sep="_")
grp <- gsub("SARS-CoV","SARSCoV",grp)
ngs$samples$group <- grp
levels = unique(ngs$samples$group)
levels

contr.matrix <- limma::makeContrasts(
    
    SARSCoV_12H_vs_Mock_12H = SARSCoV_12H - Mock_12H,
    SARSCoV_24H_vs_Mock_24H = SARSCoV_24H - Mock_24H,        
    SARSCoV_48H_vs_Mock_48H = SARSCoV_48H - Mock_48H,
    
    DOHV_12H_vs_Mock_12H = DOHV_12H - Mock_12H,
    DOHV_24H_vs_Mock_24H = DOHV_24H - Mock_24H,        
    DOHV_48H_vs_Mock_48H = DOHV_48H - Mock_48H,
    
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
    ngs, max.features = MAX.GENESETS,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("connectivity")
extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")

##FILESX <- "../libx"
ngs <- compute.extra(ngs, extra, lib.dir=FILES)

names(ngs)
ngs$timings

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------
rda.file
ngs.save(ngs, file=rda.file)












