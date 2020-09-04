##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

source("../R/pgx-functions.R")
source("../R/pgx-graph.R")
source("../R/xcr-math.r")

n=100
##G <- head(R1,100)
KNNreduceCorMatrix <- function(G, n=100) {

    ##knn.list <- apply(G,1,function(x) head(order(-x),n))
    if(0) {
        knn.list <- list()
        for(i in 1:nrow(G)) {
            ##knn.list[[i]] <- head(order(-G[i,]),n)
            knn.list[[i]] <- head(order(-G[i,]),n)
        }
    } else {
        require(parallel)
        ##knn.list <- mclapply(1:nrow(G), function(i) head(order(-abs(G[i,])),n))
        knn.list <- lapply(1:2, function(i) head(order(-G[i,]),n))
        ##knn.list <- mclapply(1:2, function(i) head(order(-G[i,]),n), mc.cores=16)

    }

    ii <- unlist(lapply(1:length(knn.list),function(i) rep(i,n)))
    jj <- as.vector(unlist(knn.list))
    xx <- G[cbind(ii,jj)]
    ##ii2 <- c(ii,jj)  ## symmetric matrix
    ##jj2 <- c(jj,ii)  ## symmetrix matrix
    ##xx <- G[cbind(ii2,jj2)]
    S = sparseMatrix( i=ii, j=jj, x=xx, use.last.ij=TRUE,
                     dims=c(nrow(G), ncol(G)) )
    colnames(S) <- colnames(G)
    rownames(S) <- rownames(G)
    return(S)
}

## Calculate correlation matrix. cosSparse is 'better' as it retains sparseness.
##

if(0) {
    require(qlcMatrix)
    load("../files/gset-sparseG-XL.rda",verbose=1)

    cat("computing GENES cosSparse...\n")
    R1 <- cosSparse( G!=0 )
    saveRDS( R1, file="GENExGENE-cosSparse-XL.rds")

    R1 <- readRDS(file="GENExGENE-cosSparse-XL.rds")
    dim(R1)

    cat("computing GENES cosSparse knn=100...\n")
    R2 <- KNNreduceCorMatrix(R1[,], 100)
    saveRDS( R2, file="GENExGENE-cosSparseKNN100-XL.rds")


    cat("computing GENES cosSparse knn=500...\n")
    R2x <- KNNreduceCorMatrix(R1, 500)
    saveRDS( R2x, file="GENExGENE-cosSparseKNN500-XL.rds")

    remove(R1)
    remove(R2)

    ##corGSETxGSET <- corSparse( t(G!=0) )
    cat("computing GSET cosSparse...\n")
    R3 <- cosSparse( t(G!=0) )
    dim(R3)
    saveRDS( R3, file="GSETxGSET-cosSparse-XL.rds")
}

if(0) {
    require(qlcMatrix)
    load("../files/gset-sparseG-XL.rda",verbose=1)

    ##X=G[1:1000,]
    ##system.time( R3  <- tcosine.sparse(X, k=100, th=0.01, block=1000, ties.method="") )
    ##system.time( R3g <- tcosine.sparse(X, k=100, th=0.01, block=1000, gpu=TRUE, ties.method="") )
    ##system.time( R3p <- tcosine.sparse.paral(X, k=100, th=0.01, block=5000, mc.cores=24 ))
    ##dim(R3p)

    cat("computing GSET cosSparse knn=100...\n")
    ##R4 <- KNNreduceCorMatrix(R3, 100)
    R4 <- tcosine.sparse.paral(G[,], k=100, th=0.01, block=5000, mc.cores=48)
    dim(R4)
    saveRDS( R4, file="GSETxGSET-cosSparseKNN100-XL.rds")

    cat("computing GSET cosSparse knn=500...\n")
    ##R4x <- KNNreduceCorMatrix(R3, 500)
    R4x <- tcosine.sparse.paral(G, k=500, th=0.01, block=5000, mc.cores=48)
    saveRDS( R4x, file="GSETxGSET-cosSparseKNN500-XL.rds")

    remove(R3)
    remove(R4)
}

require(Rtsne.multicore)
require(Rtsne)

##------------------- genes -----------------------

R <- readRDS(file="GENExGENE-cosSparseKNN500-XL.rds")
##R <- R[1:1000,1:1000]

require(sparsesvd)
require(irlba)
V = sparsesvd(1-R, rank=20)$v
V = irlba(1-R, rank=100)$v
dim(V)

tsne2d <- Rtsne( as.matrix(V), is_distance=FALSE,
                check_duplicates=FALSE, perplexity=30)$Y
rownames(tsne2d) <- rownames(R)
saveRDS(tsne2d, file="GENExGENE-tSNE2d-cosKNN500XL.rds")

tsne3d <- Rtsne( as.matrix(1-R), is_distance=TRUE,
              dim=3, check_duplicates=FALSE, perplexity=30)$Y
rownames(tsne3d) <- rownames(R)
saveRDS(tsne3d, file="GENExGENE-tSNE3d-cosKNN500XL.rds")


##------------------- gene sets -----------------------
R <- readRDS(file="GSETxGSET-cosSparseKNN500-XL.rds")
tsne2d <- Rtsne( as.matrix(1-R), is_distance=TRUE,
                  check_duplicates=FALSE, perplexity=30)$Y
rownames(tsne2d) <- rownames(R)
saveRDS(tsne2d, file="GSETxGSET-tSNE2d-cosKNN500XL.rds")

tsne3d <- Rtsne( as.matrix(1-R), is_distance=TRUE,
                dims=3, check_duplicates=FALSE, perplexity=30)$Y
rownames(tsne3d) <- rownames(R)
saveRDS(tsne3d, file="GSETxGSET-tSNE3d-cosKNN500XL.rds")


