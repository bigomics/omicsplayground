##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

source("../R/pgx-functions.R")
source("../R/pgx-graph.R")
source("../R/xcr-math.r")

n=100
##G <- head(R1,100)
KNNreduceCorMatrix <- function(G, n=100) {

    require(parallel)
    ##knn.list <- mclapply(1:nrow(G), function(i) head(order(-abs(G[i,])),n))
    knn.list <- lapply(1:2, function(i) head(order(-G[i,]),n))
    ##knn.list <- mclapply(1:2, function(i) head(order(-G[i,]),n), mc.cores=16)

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