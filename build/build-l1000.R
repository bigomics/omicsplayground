##BiocManager::install("crossmeta")
##BiocManager::install("ccmap")

library(ccmap)
library(ccdata)
##data(cmap_es)
##dim(cmap_es)
system.time( data(l1000_es) )
dim(l1000_es)


require(fgsea)
X <- l1000_es
x.drugs <- gsub("_.*$","",colnames(X))
length(table(x.drugs))

sum(table(x.drugs)>=10)
sum(table(x.drugs)>=20)

gmt <- tapply(colnames(X), x.drugs, list)
gmt.size <- sapply(gmt,length)  ## how many profiles per drug
table(gmt.size)
##gmt <- gmt[which(gmt.size>=10)]
gmt <- gmt[which(gmt.size>=15)]
##gmt <- gmt[which(gmt.size>=20)]
length(gmt)

gg <- intersect(rownames(l1000_es),rownames(cmap_es))
sel <- which(x.drugs %in% names(gmt))
length(sel)
##saveRDS( X[,sel], file="../files/l1000_es_1763drugs.rds")
saveRDS( X[,sel], file="../files/l1000_es_5685drugs.rds")
##saveRDS( X[,sel], file="../files/l1000_es_8221drugs.rds")






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
