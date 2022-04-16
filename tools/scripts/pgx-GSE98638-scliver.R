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
FILES
MAX.GENES = 8000

COMPARE.CLUSTERS=FALSE
##COMPARE.CLUSTERS=TRUE
DOWNSAMPLE=75

##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/GSE98638-scliver-TEST.pgx"
rda.file="../data/GSE98638-scliver.pgx"
##if(COMPARE.CLUSTERS) rda.file <- sub(".pgx$",paste0("-vsCLUST.pgx"),rda.file)
##if(DOWNSAMPLE>0) rda.file <- sub(".pgx$",paste0("-s",DOWNSAMPLE,".pgx"),rda.file)
##if(SMALL>0) rda.file <- sub(".pgx$",paste0("-",EXT,".pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "scRNA-seq"
ngs$description = "GSE98638 data set (Zheng et al., 2017). Single cell RNA sequencing of different subtypes from HCC patients, including CD8+ T cells (CD3+ and CD8+), T helper cells (CD3+, CD4+ and CD25-), and regulatory T cells (CD3+, CD4+ and CD25high). Ref: Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing. Cell 2017."

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
##geo <- getGEO("GSE98638", GSEMatrix =TRUE, AnnotGPL=TRUE)
datafile <- "/tmp/GSE98638_HCC.TCell.S5063.count.txt.gz"
if(!file.exists(datafile)) {
    system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98638/suppl/GSE98638_HCC.TCell.S5063.count.txt.gz -P /tmp")
}
counts = fread(datafile,nrow=-1000)
dim(counts)
head(counts)[,1:10]
counts = counts[!is.na(counts[["symbol"]]),]
gene = counts[["symbol"]]
counts = as.matrix(counts[,3:ncol(counts)])
rownames(counts) = gene
head(counts)[,1:10]
summary(colSums(counts))

##--------------------------------------------------------------
## DOWNSAMPLE
##--------------------------------------------------------------

DOWNSAMPLE
if(DOWNSAMPLE>0) {
    ## sample each category
    code <- substring(colnames(counts),1,3)
    table(code)
    jj <- order(-colSums(counts))
    jj <- tapply(jj, code[jj], function(x) head(x,DOWNSAMPLE))
    ##jj <- tapply(jj, code[jj], function(x) head(x,20))
    jj <- unlist(jj)
    table(code[jj])
    counts <- counts[,jj]
    dim(counts)
    summary(colSums(counts))
}

##--------------------------------------------------------------
## gene annotation
##--------------------------------------------------------------
require(org.Hs.eg.db)
GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
head(GENE.TITLE)
genes = data.frame( gene_name=rownames(counts),
                   gene_title=GENE.TITLE[rownames(counts)] )
rownames(genes) = rownames(counts)
head(genes)

##--------------------------------------------------------------
## Prepare sample table
##--------------------------------------------------------------
sample = colnames(counts)
group = substring(colnames(counts),1,3)
cell.type = c("C"="CD8.T","H"="CD4.Th","R"="CD4.Treg",
              "S"="unkown")[substring(group,3,3)]
tissue.type = c("P"="blood","N"="normal.liver","T"="tumor",
                "J"="joint")[substring(group,1,1)]
sampleTable = data.frame(
    ##sample = sample,
    group = group,
    cell.type = cell.type,
    tissue.type = tissue.type
)
rownames(sampleTable) = colnames(counts)
table(group)

##-------------------------------------------------------------------
## Now create an PGX object
##-------------------------------------------------------------------
if(is.null(sampleTable$group)) stop("samples need group")
table(sampleTable$group)
ngs$counts <- round(counts)
ngs$samples <- sampleTable
ngs$genes = genes
##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
ngs$samples$batch <- NULL
##ngs$samples$batch <- as.integer(lib.size2)

## tagged rownames???
row.id = paste0("tag",1:nrow(ngs$genes),":",ngs$genes[,"gene_name"])
rownames(ngs$genes) = rownames(ngs$counts) = row.id
names(ngs)

##-------------------------------------------------------------------
## sample QC filtering
##-------------------------------------------------------------------
##

##-------------------------------------------------------------------
## collapse multiple row for genes by summing up counts
##-------------------------------------------------------------------
sum(duplicated(ngs$genes$gene_name))
x1 = apply(ngs$counts, 2, function(x) tapply(x, ngs$genes$gene_name, sum))
ngs$genes = ngs$genes[match(rownames(x1),ngs$genes$gene_name),]
ngs$counts = x1
rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
remove(x1)
dim(ngs$counts)

##-------------------------------------------------------------------
## gene filtering
##-------------------------------------------------------------------
##keep <- rep(TRUE,nrow(ngs$counts))
##keep <- filterByExpr(ngs)  ## default edgeR filter
##keep <- (rowSums(cpm(ngs$counts, log=TRUE) > 1) >= 0.01)
if(0) {
    keep <- (rowMeans(ngs$counts >= 3) > 0.01)
    table(keep)
    ngs$counts <- ngs$counts[keep,]
    ngs$genes  <- ngs$genes[keep,]
}

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=30 )
head(ngs$samples)

##-------------------------------------------------------------------
## take top varying
##-------------------------------------------------------------------
MAX.GENES
if(TRUE && MAX.GENES>0) {
    cat("shrinking data matrices: n=",MAX.GENES,"\n")
    logcpm = edgeR::cpm(ngs$counts, log=TRUE)
    jj <- head( order(-apply(logcpm,1,sd)), MAX.GENES )  ## how many genes?
    head(jj)
    ##bX <- bX[jj,]
    ngs$counts <- ngs$counts[jj,]
    ngs$genes  <- ngs$genes[jj,]
}
dim(ngs$counts)

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
COMPARE.CLUSTERS
if(COMPARE.CLUSTERS) {
    ## make model matrix for group vs. rest
    clusters <- ngs$samples$cluster
    table(clusters)
    contr.matrix <- makeClusterContrasts(clusters)
    contr.matrix
    
} else {
    table(ngs$samples$group)
    levels = levels(ngs$samples$group)
    levels
    contr.matrix <- limma::makeContrasts(
        ## tumor_vs_blood = (TTC + TTH + TTR + TTS)/4 - (PTC + PTH + PTR + PTS)/4,
        ## tumor_vs_normal = (TTC + TTH + TTR + TTS)/4 - (NTC + NTH + NTR)/3,
        ## tumor_vs_joint  = (TTC + TTH + TTR + TTS)/4 - (JTC + JTH + JTS)/3,
        ## joint_vs_normal = (JTC + JTH + JTS)/3 - (NTC + NTH + NTR)/3,
        ## normal_vs_blood = (NTC + NTH + NTR)/3 - (PTC + PTH + PTR + PTS)/4,
        CD8_tumor_vs_normal = TTC - NTC,
        Th_tumor_vs_normal = TTH - NTH,
        Treg_tumor_vs_normal = TTR - NTR,
        CD8_tumor_vs_joint = TTC - JTC,
        Th_tumor_vs_joint = TTH - JTH,
        ##Treg_tumor_vs_joint = TTR - JTR,
        ##Ts_tumor_vs_joint = TTS - JTS,
        CD8_tumor_vs_blood = TTC - PTC,
        Th_tumor_vs_blood = TTH - PTH,
        Treg_tumor_vs_blood = TTR - PTR,
        CD8_joint_vs_normal = JTC - NTC,
        Th_joint_vs_normal = JTH - NTH,
        ##Treg_joint_vs_normal = JTR - NTR,
        levels = levels)
    dim(contr.matrix)
    contr.matrix
}

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------
rda.file
ngs$timings <- c()

GENE.METHODS=c("ttest","ttest.welch", ## "ttest.rank",
               "voom.limma","trend.limma","notrend.limma",
               "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
GENESET.METHODS = c("fisher","gsva","ssgsea","spearman",
                    "camera", "fry","fgsea") ## no GSEA, too slow...

GENE.METHODS=c("ttest.welch", "trend.limma", "edger.qlf")
GENESET.METHODS = c("fisher","gsva","fgsea") ## no GSEA, too slow...

## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features=MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features=MAX.GENES,
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



