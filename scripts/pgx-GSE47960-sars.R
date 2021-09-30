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

rda.file="../data/GSE47960-sars.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE47960. A network integration approach to predict conserved regulators related to pathogenicity of influenza and SARS-CoV respiratory viruses (Mitchell 2013). HAE cultures were infected with SARS-CoV, SARS-dORF6 or SARS-BatSRBD and were directly compared to A/CA/04/2009 H1N1 influenza-infected cultures. Cell samples were collected at various hours post-infection for analysis."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
library(Biobase)
library(GEOquery)
library(limma)
library(hgu133plus2.db)

## load series and platform data from GEO
geo <- getGEO("GSE47960", GSEMatrix=TRUE, getGPL=TRUE)
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
tt <- sub("_B$","B",tt)
sampleTable <- do.call(rbind,strsplit(tt,split="_"))
colnames(sampleTable) <- c("code","infected","time","replicate")
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
    BAT_0h_vs_mock_0h = BAT_0h - mock_0h,
    BAT_12h_vs_mock_12h = BAT_12h - mock_12h,
    BAT_24h_vs_mock_24h = BAT_24h - mock_24h,
    BAT_36h_vs_mock_36h = BAT_36h - mock_36h,
    BAT_48h_vs_mock_48h = BAT_48h - mock_48h,
    BAT_60h_vs_mock_60h = BAT_60h - mock_60h,
    BAT_72h_vs_mock_72h = BAT_72h - mock_72h,
    BAT_84h_vs_mock_84h = BAT_84h - mock_84h,
    BAT_96h_vs_mock_96h = BAT_96h - mock_96h,
    
    H1N1_0h_vs_mock_0h = H1N1_0h - mock_0h,
    H1N1_6h_vs_mock_6h = H1N1_6h - mock_6h,
    H1N1_12h_vs_mock_12h = H1N1_12h - mock_12h,
    H1N1_18h_vs_mock_18h = H1N1_18h - mock_18h,
    H1N1_24h_vs_mock_24h = H1N1_24h - mock_24h,
    H1N1_36h_vs_mock_36h = H1N1_36h - mock_36h,
    H1N1_48h_vs_mock_48h = H1N1_48h - mock_48h,
    
    dORF6_0h_vs_mock_0h = dORF6_0h - mock_0h,
    dORF6_12h_vs_mock_12h = dORF6_12h - mock_12h,
    dORF6_24h_vs_mock_24h = dORF6_24h - mock_24h,
    dORF6_36h_vs_mock_36h = dORF6_36h - mock_36h,
    dORF6_48h_vs_mock_48h = dORF6_48h - mock_48h,
    dORF6_60h_vs_mock_60h = dORF6_60h - mock_60h,
    dORF6_72h_vs_mock_72h = dORF6_72h - mock_72h,
    dORF6_84h_vs_mock_84h = dORF6_84h - mock_84h,
    dORF6_96h_vs_mock_96h = dORF6_96h - mock_96h,
    
    icSARS_0h_vs_mock_0h = icSARS_0h - mock_0h,
    icSARS_12h_vs_mock_12h = icSARS_12h - mock_12h,
    icSARS_24h_vs_mock_24h = icSARS_24h - mock_24h,
    icSARS_36h_vs_mock_36h = icSARS_36h - mock_36h,
    icSARS_48h_vs_mock_48h = icSARS_48h - mock_48h,
    icSARS_60h_vs_mock_60h = icSARS_60h - mock_60h,
    icSARS_72h_vs_mock_72h = icSARS_72h - mock_72h,
    icSARS_84h_vs_mock_84h = icSARS_84h - mock_84h,
    icSARS_96h_vs_mock_96h = icSARS_96h - mock_96h,
    
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












