##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

########################################################################
##
## Functions for batch correction
##
########################################################################

##use.design=TRUE;dist.method="cor";center.x=FALSE;replace=FALSE
library(sva)
library(limma)

gx.nnmcorrect2 <- function(X, y, use.design=TRUE, dist.method="cor",
                           center.x=FALSE, replace=FALSE) 
{
    ##-----------------------------------------------------
    ## nearest-neighbour matching for batch correction
    ##-----------------------------------------------------

    ## compute full NNM-pairing
    grp <- unique(y)
    ny <- length(grp)
    dX <- scale(X)
    if(center.x) dX <- dX - rowMeans(dX,na.rm=TRUE)
    if(dist.method=="cor") {
        D <- 1 - cor(dX)
    } else {
        D <- as.matrix(dist(t(dX)))
    }
    remove(dX)

    ##replace=TRUE
    ##replace=FALSE
    maxD <- max(D)
    P <- matrix(NA,nrow=ncol(X),ncol=ny)
    rownames(P) <- colnames(X)
    colnames(P) <- grp
    i=1;k=2
    for(i in 1:nrow(P)) {
        for(k in 1:ny) {
            jj <- which(y == grp[k])
            ii <- which(y == y[i])
            dd <- D[i,jj]
            if(!replace) dd <- dd + maxD*table(factor(P[ii,k],levels=jj)) ## penalty factor
            P[i,k] <- jj[which.min(dd)]
        }
    }
    P
    
    ## imputing full paired data set
    kk <- as.vector(P)
    full.y <- y[kk]
    full.pairs <- rep(rownames(P),ncol(P))
    full.X <- X[,kk]

    ## remove pairing effect
    if(use.design) {
        design <- model.matrix( ~ full.y )
        full.X <- limma::removeBatchEffect(full.X, batch=full.pairs,
                                           design=design)
    } else {
        full.X <- limma::removeBatchEffect(full.X, batch=full.pairs)
    }

    ## now contract to original samples
    full.idx <- rownames(P)[kk]
    cx <- do.call(cbind, tapply(1:ncol(full.X), full.idx,
                                function(i) rowMeans(full.X[,i,drop=FALSE])))
    cx <- cx[,colnames(X)]
    
    res <- list(X=cx, pairings=P)
    return(res)
}

gx.nnmcorrect <- function(x, y, k=3, dist.method="cor") {
    ##-----------------------------------------------------
    ## nearest-neighbour matching for batch correction
    ##-----------------------------------------------------
    ## distance metric for matching
    if(dist.method=="cor") {
        D <- 1 - cor(x)
    } else {
        D <- as.matrix(dist(t(x)))
    }

    ## masking for other groups
    diag(D) <- NA
    i=1
    for(i in 1:nrow(D)) {
        jj <- which(y == y[i])
        D[i,jj] <- NA
    }
    
    dx <- x
    j=1
    pairings <- matrix(NA,ncol(x),k)
    for(j in 1:ncol(x)) {
        nj <- which(y!=y[j])
        nn <- intersect(order(D[j,]),nj)
        nn <- head(nn,k)
        dx[,j] <- x[,j] - rowMeans(x[,nn,drop=FALSE])
        pairings[j,] <- head(c(nn,rep(NA,k)),k)
    }
    cx <- rowMeans(x) + dx / 2    
    res <- list(X=cx, pairings=pairings)
    return(res)
}

gx.nnmcorrect.SAVE <- function(x, y, k=3) {
    ##-----------------------------------------------------
    ## nearest-neighbour matching for batch correction
    ##-----------------------------------------------------
    xcor <- cor(x)
    diag(xcor) <- 0
    nx <- x
    j=1
    for(j in 1:ncol(x)) {
        nj <- which(y!=y[j])
        nn <- intersect(order(-xcor[j,]),nj)
        nn <- head(nn,k)
        nx[,j] <- x[,j] - rowMeans(x[,nn,drop=FALSE])
    }
    nx <- nx + rowMeans(x)
    return(nx)
}

gx.qnormalize <- function(X) {
    ##-----------------------------------------------------
    ## quantile normalization
    ##-----------------------------------------------------
    require(affy)
    if( max(X, na.rm=TRUE) < 40 && min(X, na.rm=TRUE) > 0) {
        tmp <- normalize.qspline( 2**X, na.rm=TRUE, verbose=FALSE)
        tmp <- log2(tmp)
    } else {
        tmp <- normalize.qspline( X, na.rm=TRUE, verbose=FALSE)
    }
    rownames(tmp) <- rownames(X)
    colnames(tmp) <- colnames(X)
    X <- tmp
}


