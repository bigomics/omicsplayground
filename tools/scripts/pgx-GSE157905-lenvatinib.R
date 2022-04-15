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
FILESX = "../libx"
PGX.DIR = "../data"
source("../R/pgx-include.R")
source("../R/pgx-getgeo.R")
##source("options.R")

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
library(Biobase)
library(GEOquery)
library(limma)

## load series and platform data from GEO
##archs.h5 = file.path(FILESX,"human_matrix.h5")
##ngs <- pgx.getGEOseries("GSE157905", archs.h5=archs.h5)
##names(ngs)

gsedata <- fread.csv("~/Playground/data/GSE/GSE157905_readcounts.txt.gz")

genes  <- gsedata[,17:22]
counts <- as.matrix(gsedata[,1:16])
samples <- do.call(rbind,strsplit(colnames(counts),split="[_-]"))
samples <- data.frame(samples)
colnames(samples) <- c("batch","nr","cell.line","treatment","barcode")
samples <- samples[,c("nr","cell.line","treatment")]

short.names <- paste(samples$nr,samples$cell.line,samples$treatment,sep='-')
rownames(samples) <- short.names
colnames(counts)  <- short.names

Y <- samples[,'treatment',drop=FALSE]
ct <- makeDirectContrasts( Y, ref='CT')
head(ct$exp.matrix)

if(1) {
    ## Batch correct on cell-line
    X  <- edgeR::cpm(counts, log=TRUE, prior.count=2)
    b <- samples$cell.line
    y <- samples$treatment
    
    ##y <- pgx.getConditions(ct$exp.matrix)
    ##table(y)
    
    bX <- sva::ComBat(X, batch=b)
    bX <- gx.nnmcorrect(bX, y)$X
    bX <- pgx.svaCorrect(bX, y)
    
    counts <- pmax(2**bX, 0)
    max(counts)
    min(counts)

    if(0) {

        clust1 <- pgx.clusterBigMatrix(X, dims=2)
        clust2 <- pgx.clusterBigMatrix(bX, dims=2)
        names(clust1)
        
        par(mfrow=c(3,3))
        pgx.scatterPlotXY(clust1[[1]], var=samples$treatment)
        pgx.scatterPlotXY(clust2[[1]], var=samples$treatment)    
        pgx.scatterPlotXY(clust1[[2]], var=samples$treatment)
        pgx.scatterPlotXY(clust2[[2]], var=samples$treatment)    
        pgx.scatterPlotXY(clust1[[3]], var=samples$treatment)
        pgx.scatterPlotXY(clust2[[3]], var=samples$treatment)    
    }
}

ngs <- pgx.createPGX(
    counts = counts,
    samples = samples,
    contrasts = ct$exp.matrix
)
names(ngs)

gx.methods    = c("trend.limma")
gset.methods  = c("fisher")
extra.methods  = c("meta.go","infer")

gx.methods    = c("trend.limma","edger.qlf","deseq2.wald")
gset.methods  = c("fisher","gsva","fgsea","spearman")
extra.methods  = c("meta.go","infer","drugs","deconv","wordcloud","connectivity")

ngs <- pgx.computePGX(
    ngs,
    max.genes = 10000,
    max.genesets = 5000, 
    gx.methods = gx.methods,
    gset.methods = gset.methods,
    extra.methods = extra.methods,
    use.design = TRUE,      ## no.design+prune are combined 
    prune.samples = FALSE,  ##
    do.cluster = TRUE,                
    progress = NULL,
    lib.dir = FILES
)

##ngs <- pgx.computeExtra(ngs, extra="drugs", lib.dir=FILES)
names(ngs$drugs)
names(ngs$drugs[[1]])

##-------------------------------------------------------------------
## save object
##-------------------------------------------------------------------

rda.file="../data/GSE157905-lenvatinib-bc.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$organism = 'human'
ngs$datatype = "RNA-seq"
ngs$description = "GSE157905. RNA Sequencing of HCC cells after lenvatinib, gefitinib, and combination treatment (Jin et al., Nature 2021). Batchcorrected: combat+nnm+sva"

rda.file
ngs.save(ngs, file=rda.file)

