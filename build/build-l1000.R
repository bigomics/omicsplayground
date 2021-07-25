##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. A                                        ll rights reserved.
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
dim(X)
sum(grepl("_oe_",colnames(X)))
head(grep("_oe_",colnames(X),value=TRUE))
sum(grepl("_sh_",colnames(X)))
sum(grepl("_lig_",colnames(X)))
colnames(X) <- gsub("_oe_","-oe_",colnames(X))
colnames(X) <- gsub("_sh_","-sh_",colnames(X))
colnames(X) <- gsub("_lig_","-lig_",colnames(X))
dim(X)

X <- X[,order(-colMeans(X**2))]  ## order on intensity
xdrugs <- gsub("_.*$","",colnames(X))
names(xdrugs) <- colnames(X)
length(table(xdrugs))
hist(table(xdrugs),breaks=100,xlim=c(0,100))
table(xdrugs %in% rownames(D))
table(unique(xdrugs) %in% rownames(D))

if(0) {
    ## only with drug annotation??? Gene perturbation will be removed!!
    length(xdrugs)
    has.annot <- xdrugs %in% rownames(D)
    table(has.annot)
    X <- X[,has.annot]
    xdrugs <- xdrugs[has.annot]
}

##sum(table(xdrugs)>=10)
sum(table(xdrugs)>=10)
sum(table(xdrugs)>=15)
sum(table(xdrugs)>=20)

gmt <- tapply(colnames(X), xdrugs, list)
gmt.size <- sapply(gmt,length)  ## how many profiles per drug
table(gmt.size)

head(grep("-sh",names(gmt),value=TRUE))
head(grep("-oe",names(gmt),value=TRUE))
head(grep("-lig",names(gmt),value=TRUE))
table(gmt.size[grep("-sh",names(gmt))])
table(gmt.size[grep("-oe",names(gmt))])
table(gmt.size[grep("-lig",names(gmt))])

if(0) {
    ## -----------------------------------------------------------
    ## DRUG SIGNATURES
    ## -----------------------------------------------------------

    nmin = 10
    ##nmin = 15
    ##nmin = 20
    sel1 <- grepl("-sh|-oe|-lig",names(gmt))  ## gene target
    has.annot <- (names(gmt) %in% rownames(D))    
    ##gmt1 <- gmt[which(gmt.size >= nmin & !sel1 & has.annot)]
    gmt1 <- gmt[which(gmt.size >= nmin & !sel1)]    
    length(gmt1)
    ##gmt1 <- lapply(gmt1, function(g) head(g,nmin))
    gmt1 <- lapply(gmt1, function(g) head(g,2*nmin))    
    length(gmt1)

    X1 <- X[,which(colnames(X) %in% unlist(gmt1))]  ## filtered on intensity
    ## X1 <- X[,which(xdrugs %in% names(gmt1))]    
    dim(X)
    dim(X1)    
    X1 <- round(X1, digits=3)
    
    ## file must be smaller than 100Mb for GitHub
    ##saveRDS(X1, file="../lib/l1000_es_n20d1011.rds")
    ##write.csv(X1, file=gzfile("../lib/l1000_es_n20d1011.csv.gz")) 

    saveRDS(X1, file="../lib/l1000_es_n10d5263.rds")  ## a???
    write.csv(X1, file=gzfile("../lib/l1000_es_n10d5263.csv.gz")) 
    ##readr::write_csv(data.frame(X1), file="../lib/l1000_es_n15d3479.csv2.gz")      

    rds <- readRDS(file="../lib/l1000_es_n10d5263.rds")  ## a???
    fz <- fread(file="../lib/l1000_es_n10d5263.csv.gz")
    
}


if(0) {
    ## -----------------------------------------------------------
    ## GENE PERTURBATIONS. Those OE/SH/LIG profiles that target a
    ## single gene.    
    ## -----------------------------------------------------------

    nmin = 8
    ##nmin=10;nmax=10
    nmin=8;nmax=20
    sel1 <- grepl("-sh|-oe|-lig",names(gmt))
    gmt1 <- gmt[which(gmt.size >= nmin & sel1)]
    length(gmt1)
    gmt1 <- lapply(gmt1, function(g) head(g,nmax))    
    X1 <- X[,which(colnames(X) %in% unlist(gmt1))]
    X1 <- round(X1, digits=3)
    dim(X)
    dim(X1)
    length(gmt1)
    ##saveRDS( X1, file="../lib/l1000_gpert_n10g1766.rds")
    ##write.csv(X1, file=gzfile("../lib/l1000_gpert_n10g1766.csv.gz"))
    saveRDS(X1, file="../lib/l1000_gpert_n8m20g5812.rds")
    write.csv(X1, file=gzfile("../lib/l1000_gpert_n8m20g5812.csv.gz"))
    ##saveRDS( X1[,], file="../lib/l1000_es_n8g5812.rds")
    ##write.csv(X1, file=gzfile("../lib/l1000_es_n8g5812.csv.gz"))
        
}

if(0) {
    ##
    ## TESTING
    ##
    source("../R/pgx-include.R")
    source("../R/pgx-drugs.R")

    load("../data/GSE102908-ibet.pgx")
    F <- pgx.getMetaMatrix(ngs)$fc
    names(ngs$drugs)
    names(ngs$drugs[[1]])
    lapply(ngs$drugs[[1]],dim)    
    dr <- ngs$drugs[[1]]
    
    FILES="../files"
    library(data.table)
    X <- readRDS(file="../lib/l1000_es.rds")    
    
    dim(X)
    xdrugs <- gsub("_.*$","",colnames(X))
    length(table(xdrugs))

    gg <- intersect(rownames(X),rownames(F))
    jj <- sample(1:ncol(X),3)
    X1 <- scale(X[gg,-jj])
    X1 <- scale(X[gg,])    
    pos <- uwot::umap(t(X1), metric='euclidean', fast_sgd=TRUE)    
    rownames(pos) <- colnames(X1)
    ##pos <- pos.compact(pos)

    dsea <- ngs$drugs[[1]]
    names(dsea)
    gmt <- grep("panobino",rownames(dsea$stats),value=TRUE)
    gsea.enplot(dsea$stats[,1], gmt, main="panobinostat")

    gg <- intersect(rownames(X),rownames(F))
    rho <- cor(X[gg,], F[gg,], use='pairwise')[,1]
    names(rho) <- colnames(X1)
    tail(sort(rho),20)
    nnb <- sort(table(sub("_.*","",names(tail(sort(rho),100)))),decreasing=TRUE)
    
    xdrugs <- sub("_.*","",colnames(X))
    out <- pgx.computeDrugEnrichment(ngs, X, xdrugs, methods="GSEA")[[1]]
    head(out)
    head(out$X)        
    
    ## Train UMAP
    obj_umap <- uwot::umap(t(X1), metric='correlation', ret_model=TRUE)        
    pos <- obj_umap$embedding
    rownames(pos) <- colnames(X1)

    jj <- grep("scriptaid",colnames(X))
    X2 <- X1[,jj]
    X2 <- X1[,head(jj,2)]
    X2 <- F[match(rownames(X1),rownames(F)),]
    rownames(X2) <- rownames(X1)

    ## Project query onto previous UMAP
    qpos <- uwot::umap_transform(t(X2), obj_umap, verbose=TRUE)    
    rownames(qpos) <- colnames(X2)

    head(X2)
    rho <- cor(X1, X2, use='pairwise')[,1]
    ##rho <- cor(apply(X1,2,rank), apply(X2,2,rank), use='pairwise')[,2]    
    names(rho) <- colnames(X1)
    tail(sort(rho),20)

    pgx.scatterPlotXY.BASE(pos, var=rho, cex=0.8)
    ##pgx.scatterPlotXY.PLOTLY(pos, var=varx)
    points(qpos, pch=20, cex=1.5)
    text(qpos, rownames(qpos), cex=0.9, pos=3, offset=0.3)
    
    qpos
    dpos <- sqrt(colMeans((t(pos) - qpos[1,])**2))
    head(sort(dpos),10)
    topdr.tab <- sort(table(sub("_.*","",names(head(sort(dpos),100)))),decreasing=TRUE)
    head(topdr.tab)
    aa <- ngs$drugs[[1]]$annot
    topdr <- intersect(names(topdr.tab),rownames(aa))
    aa[head(topdr,10),]
    xdrugs <- sub("_.*","",rownames(pos))
    hi <- which(xdrugs == 'scriptaid')
    points(pos[hi,], pch=20, cex=1, col='red2')
    text(pos[hi,], rownames(pos)[hi], cex=0.7, pos=3, offset=0.3)
    
    out <- pgx.computeDrugEnrichment(
        ngs, X, xdrugs, methods="GSEA",
        contrast=NULL)[[1]]
    head(out)

    jj <- grep("crizotinib",colnames(X))
    pairs(X[,jj], pch='.')
    cor(X[,jj])    

       
    pgx.scatterPlotXY.BASE(
        pos, var=rho, opacity=0.2, cex=0.6,
        hilight=h1, hilight.col='red', hilight.cex=1.2)

    

}
