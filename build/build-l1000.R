##BiocManager::install("crossmeta")
##BiocManager::install("ccmap")

library(ccmap)
library(ccdata)
##data(cmap_es)
##dim(cmap_es)
system.time( data(l1000_es) )
dim(l1000_es)

D <- read.csv("../lib/L1000_repurposing_drugs.txt",sep="\t",
              skip=13, header=FALSE,row.names=1)

require(fgsea)
X <- l1000_es
X <- X[,order(-colMeans(X**2))]
x.drugs <- gsub("_.*$","",colnames(X))
names(x.drugs) <- colnames(X)
length(table(x.drugs))
table(x.drugs %in% rownames(D))

## only with annotation...
sel <- x.drugs %in% rownames(D)
X <- X[,sel]
x.drugs <- gsub("_.*$","",colnames(X))

sum(table(x.drugs)>=10)
sum(table(x.drugs)>=15)
sum(table(x.drugs)>=20)

gmt <- tapply(colnames(X), x.drugs, list)
gmt.size <- sapply(gmt,length)  ## how many profiles per drug
table(gmt.size)

nmin = 20
gmt1 <- gmt[which(gmt.size >= nmin)]
length(gmt1)
gmt1 <- lapply(gmt1, function(g) head(g,nmin))

X1 <- X[,which(colnames(X) %in% unlist(gmt1))]
dim(X)
dim(X1)
length(gmt1)

saveRDS( X1[,], file="../lib/l1000_es_5685drugs.rds")
saveRDS( X1[,], file="../lib/l1000_es_520drugs.rds")


if(0) {
    system.time( X <- readRDS(file="../files/l1000_es_5685drugs.rds"))
    source("../R/pgx-drugs.R")
    source("../R/pgx-functions.R")
    load("../pgx/guarda2019-myc-12k-LT.pgx")
    FILES="../files"
    X <- readRDS(file=file.path(FILES,"l1000_es_5685drugs.rds"))
    x.drugs <- gsub("_.*$","",colnames(X))
    length(table(x.drugs))
    avgX <- readRDS(file=file.path(FILES,"l1000_es_5685drugsAVG.rds"))
    
    contrast <- c("IL15_Torin_vs_IL15_nd_1","IL15_Torin_vs_NT_nd_1")
    
    out <- pgx.computeDrugEnrichment(
        ngs, X, x.drugs, methods=c("cor","GSEA"),
        avgX=avgX, contrast=contrast)
    names(out)
    
    F <- pgx.getMetaFoldChangeMatrix(ngs)$fc
    dim(F)
    rownames(F) <- toupper(rownames(F))
    out <- pgx.computeDrugEnrichment(
        F, X, x.drugs, methods=c("cor","GSEA"),
        avgX=avgX, contrast=contrast)
    
    names(out)
    head(out[[1]]$X)
    head(out[[2]]$X)
    
    dim(avgX)
    dx.rho <- pgx.computeDrugEnrichment(
        obj=avgX, X, x.drugs, methods="cor", contrast=NULL)
    
    dim(dx.rho[[1]]$X)
    xdist <- 1- dx.rho[[1]]$X * (1 - dx.rho[[1]]$Q)
    pos <- Rtsne( xdist, is_distance=TRUE)$Y
    rownames(pos) <- rownames(xdist)
    write.csv(pos, file="l1000-drugs-tsne.csv")
    
    require(scatterD3)
    scatterD3( pos[,1], pos[,2], lab=rownames(xdist))

}
