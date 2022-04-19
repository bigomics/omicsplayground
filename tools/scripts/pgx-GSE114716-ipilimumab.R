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
FILES
MAX.GENES = 8000
BATCH.CORRECT=1


##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/GSE114716-ipilimumab.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "RNA-seq"
ngs$description = "GSE114716 data set. CD4 t cell from patients with metastatic melanoma who received Ipilimumab at baseline and after 3 doses of therapy with ipilimumab."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

library(Biobase)
library(GEOquery)
library(limma)
library(hgu133plus2.db)

## load series and platform data from GEO
gset <- getGEO("GSE114716", GSEMatrix=TRUE, getGPL=FALSE)
attr(gset, "names")
gset <- gset[[1]]
dim(exprs(gset))

## still need to get the matrix...
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114716/suppl/GSE114716_raw.counts.hs.xlsx -P /tmp")
if(!require(xlsx)) {
    ## need to do: >> sudo R CMD javareconf
    install.packages("rJava")
    install.packages("xlsx")
}
library(xlsx)
counts <- read.xlsx2("/tmp/GSE114716_raw.counts.hs.xlsx", 1)  ## yikes.. xlsx...
genes  <- as.character(counts[,1])
counts <- apply(counts[,-1],2,function(x) as.integer(as.character(x)))
X = apply(counts[,], 2, function(x) tapply(x, genes, sum))
dim(X)
head(X)[,1:4]
colSums(X)    

pdata = pData(gset)
head(pdata)
tt <- as.character(pdata$title)
tt <- gsub("CD4 T .* at ","",tt)
tt <- gsub("CD4 T .* of ","",tt)
tt <- gsub("patient|from","",tt)
tt <- gsub("  "," ",tt)
tt <- gsub("  "," ",tt)
sampleTable <- data.frame(do.call(rbind,strsplit(tt, split=" ")))
colnames(sampleTable) <- c("treatment","patient")
sampleTable$patient <- paste0("pt",sampleTable$patient)

## convert affymetrix ID to GENE symbol
symbol <- unlist(as.list(org.Hs.egSYMBOL))
jj <- which(rownames(X) %in% symbol)
X = X[jj,]
dim(X)
head(X)[,1:4]

## conform tables
sample.names <- apply(sampleTable,1,paste,collapse="_")
rownames(sampleTable) = colnames(X) = sample.names

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
chrloc = as.list(org.Hs.egCHRLOC)
names(chrloc) = gene.symbol
chrloc <- chrloc[rownames(X)]
loc <- abs(sapply(chrloc, "[", 1))
chrom <- sapply(chrloc, function(s) names(s)[1])
chrom[sapply(chrom,is.null)] <- NA
chrom <- as.vector(unlist(chrom))

genes = data.frame( gene_name=rownames(X),
                   gene_title=gene_title,
                   chr=chrom, pos=loc)
##genes = apply(genes,2,as.character)
head(genes)
table(rownames(genes) == rownames(X))

if(0) {
    ## check txlen vs average counts
    genes <- ngs.getGeneAnnotation(rownames(X))
    head(genes)
    par(mfrow=c(2,2), mar=c(4,4,2,2))
    plot( genes$tx_len, rowMeans(X),pch=".")
    plot( log(1+genes$tx_len), log(1+rowMeans(X)),pch=".")
    plot( log(1+genes$tx_len), log(1+rowMeans(X)),pch=".")
    xlen <- imputeMedian(genes$tx_len)
    plot( log(1+genes$tx_len), log(1+rowMeans(X)/xlen),pch=".")

}

##--------------------------------------------------------------------
## check if batch correction is needed
##--------------------------------------------------------------------
BATCH.CORRECT
if(BATCH.CORRECT) {
    require(sva)
    batch <- sampleTable$patient
    design = model.matrix( ~ as.character(sampleTable$treatment))
    ##bX = ComBat(X, batch=as.character(batch), mod=design)
    bX = removeBatchEffect(log(0.0001+X), batch=as.character(batch),
                           ##batch2=as.character(sampleTable$gender),
                           design=design)
    ##X = exp(normalizeQuantiles(bX))
    X = exp(bX)
}

##-------------------------------------------------------------------
## Now create an DGEList object  (see tximport Vignette)
##-------------------------------------------------------------------
ngs$counts <- X  ## treat as counts
ngs$samples <- sampleTable
ngs$genes = genes

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters 
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
    Ipi_vs_baseline = Ipi - baseline,
    baseline_vs_Ipi = baseline - Ipi,
    levels = levels)
contr.matrix

GENE.METHODS=c("ttest.welch", "trend.limma", "edger.qlf")
GENESET.METHODS = c("fisher","gsva","fgsea") ## no GSEA, too slow...

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------
## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features=MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features=MAX.GENES,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("drugs")
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

## load(rda.file, verbose=TRUE)












