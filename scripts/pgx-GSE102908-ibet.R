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

rda.file="../data/GSE102908-ibet.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
library(Biobase)
library(GEOquery)
library(limma)

## load series and platform data from GEO
archs.h5 = file.path(FILESX,"human_matrix.h5")
ngs <- pgx.getGEOseries("GSE102908", archs.h5=archs.h5)
names(ngs)

ngs <- pgx.createPGX(
    counts = ngs$counts,
    samples = ngs$samples,
    contrasts = ngs$contrasts
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
    max.genes = 20000,
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

ngs <- pgx.computeExtra(ngs, extra="drugs", lib.dir=FILES)
names(ngs$drugs)
names(ngs$drugs[[1]])

##-------------------------------------------------------------------
## save object
##-------------------------------------------------------------------

ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$organism = 'human'
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE102908. RNAseq data from SCCOHT1 and OVCAR8 ovarian cancer cells treated with BET inhibitors. Analysis of mRNA profile of 2 cell lines exposed to DMSO, OTX015 for 4 and 24 hours in duplicate."

rda.file
ngs.save(ngs, file=rda.file)


if(0) {
    
    source("../R/pgx-include.R")
    load(rda.file,verbose=1)

    extra = "drugs"
    ngs <- compute.extra(ngs, extra, lib.dir=FILES, sigdb=sigdb) 
    
    names(ngs$drugs[[1]])
    pos <- ngs$drugs[[1]]$clust
    var <- ngs$drugs[[1]]$R[,1]
    x1  <- ngs$drugs[[1]]$X[,1]
    aa <- ngs$drugs[[1]][['annot']]
    tail(sort(x1))

    xdrugs <- sub("_.*","",rownames(pos))
    xpos <- apply(pos,2,function(x) tapply(x,xdrugs,median))
    dim(xpos)
    
    pgx.scatterPlotXY.BASE(xpos, var=x1, col='grey70', ## hilight=h1,
                           cex=1.7, softmax=0, opacity=0.9,
                           zsym=TRUE, cex.lab=1.3)

    pos1 <- rbind(pos, xpos)
    xvar <- tapply(var,xdrugs,median)
    var1 <- c(var, xvar)
    xdrugs1 <- c(xdrugs, rownames(xpos))
    h1 <- tail(names(sort(x1)),20)
    h1 <- rownames(xpos)
    h1 <- rownames(pos1)[which(xdrugs1 == 'warfarin')]


    d = 'wortmannin'
    d = 'digoxin'
    d = 'digitoxin'
    top.drugs <- names(sort(x1,decreasing=TRUE))
    par(mfrow=c(5,7))
    for(d in head(top.drugs,35)) {
        h1 <- rownames(pos1)[which(xdrugs1 == d)]
        pgx.scatterPlotXY.BASE(pos1, var=var1,
                               hilight=h1, hilight2=NULL, 
                               hilight.col = 'red', hilight.cex = 1.2,
                               cex=0.7, opacity=0.05, softmax=0, 
                               rstep = 0.1, repel=1, title=d,
                               zsym=TRUE, cex.lab=1.3, dlim=0.01)
    }
    
    pgx.scatterPlotXY.GGPLOT(pos1, var=var1, hilight=h1, legend=TRUE,
                             cex=0.7, softmax=1, opacity=0.3,
                             zsym=TRUE, cex.lab=1.3)


    rda.file
    ngs.save(ngs, file=rda.file)
    
    cpal <- colorRampPalette(c("blue4","lightskyblue1","lightyellow","rosybrown1","red4"))(64)
    cpal <- colorRampPalette(c("blue4","lightskyblue1","grey95","rosybrown1","red4"))(64)
    cpal <- colorRampPalette(c("#313695","#FFFFDF","#A50026"))(64)
    cpal <- colorspace::diverge_hcl(64, c=60, l=c(30,100), power=1)
    plot( 1:64, 1:64, col=cpal, pch=20, cex=5)

}













