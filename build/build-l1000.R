##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##BiocManager::install("crossmeta")
##BiocManager::install("ccmap")

library(ccmap)
library(ccdata)
library(fgsea)

##data(cmap_es)
##dim(cmap_es)
system.time( data(l1000_es) )
dim(l1000_es)

D <- read.csv("../lib/L1000_repurposing_drugs.txt",sep="\t",
              skip=13, header=FALSE,row.names=1)

X <- l1000_es
sum(grepl("_oe_",colnames(X)))
sum(grepl("_sh_",colnames(X)))
sum(grepl("_lig_",colnames(X)))
colnames(X) <- gsub("_oe_","-oe_",colnames(X))
colnames(X) <- gsub("_sh_","-sh_",colnames(X))
colnames(X) <- gsub("_lig_","-lig_",colnames(X))

X <- X[,order(-colMeans(X**2))]
xdrugs <- gsub("_.*$","",colnames(X))
names(xdrugs) <- colnames(X)
length(table(xdrugs))
table(xdrugs %in% rownames(D))

## only with annotation...
if(0) {
    sel <- xdrugs %in% rownames(D)
    X <- X[,sel]
}

xdrugs <- gsub("_.*$","",colnames(X))
##sum(table(xdrugs)>=10)
sum(table(xdrugs)>=15)
sum(table(xdrugs)>=20)

gmt <- tapply(colnames(X), xdrugs, list)
gmt.size <- sapply(gmt,length)  ## how many profiles per drug
table(gmt.size)

head(grep("-sh",names(gmt),value=TRUE))
table(gmt.size[grep("-sh",names(gmt))])

nmin = 15
nmin = 20
gmt1 <- gmt[which(gmt.size >= nmin)]
length(gmt1)
gmt1 <- lapply(gmt1, function(g) head(g,nmin))

X1 <- X[,which(colnames(X) %in% unlist(gmt1))]
dim(X)
dim(X1)
length(gmt1)

saveRDS( X1[,], file="../lib/l1000_es_n20d1043.rds")
saveRDS( X1[,], file="../lib/l1000_es_n15d3756.rds")

if(0) {
    library(data.table)
    X <- readRDS(file="../lib/l1000_es_n15d3756.rds")
    X <- readRDS(file="../lib/l1000_es_n20d1043.rds")
    dim(X)
    
    system.time( X <- readRDS(file="../lib/l1000_es.rds") )
    system.time( X <- data.table::fread(file="../lib/l1000_es.csv.gz") )
    dim(X)
    X1 <- as.matrix(fread("../lib/l1000_es.csv.gz"), rownames=1)
    dim(X1)
    head(X1)[,1:4]
    
    source("../R/pgx-drugs.R")
    source("../R/pgx-functions.R")
    load("../pgx/guarda2019-myc-12k-LT.pgx")
    FILES="../files"
    X <- readRDS(file=file.path(FILES,"l1000_es_5685drugs.rds"))
    xdrugs <- gsub("_.*$","",colnames(X))
    length(table(xdrugs))
    avgX <- readRDS(file=file.path(FILES,"l1000_es_5685drugsAVG.rds"))
    
    contrast <- c("IL15_Torin_vs_IL15_nd_1","IL15_Torin_vs_NT_nd_1")
    
    out <- pgx.computeDrugEnrichment(
        ngs, X, xdrugs, methods=c("cor","GSEA"),
        avgX=avgX, contrast=contrast)
    names(out)
    
    F <- pgx.getMetaFoldChangeMatrix(ngs)$fc
    dim(F)
    rownames(F) <- toupper(rownames(F))
    out <- pgx.computeDrugEnrichment(
        F, X, xdrugs, methods=c("cor","GSEA"),
        avgX=avgX, contrast=contrast)
    
    names(out)
    head(out[[1]]$X)
    head(out[[2]]$X)
    
    dim(avgX)
    dx.rho <- pgx.computeDrugEnrichment(
        obj=avgX, X, xdrugs, methods="cor", contrast=NULL)
    
    dim(dx.rho[[1]]$X)
    xdist <- 1- dx.rho[[1]]$X * (1 - dx.rho[[1]]$Q)
    pos <- Rtsne( xdist, is_distance=TRUE)$Y
    rownames(pos) <- rownames(xdist)
    write.csv(pos, file="l1000-drugs-tsne.csv")
    
    require(scatterD3)
    scatterD3( pos[,1], pos[,2], lab=rownames(xdist))

}
