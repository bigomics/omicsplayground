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

##replace=FALSE;center.x=FALSE
gx.nnmcorrect <- function(X, y, use.design=TRUE, dist.method="cor",
                          center.x=TRUE, center.m=TRUE,
                          b.method="new", replace=FALSE) 
{
    ## Nearest-neighbour matching for batch correction. This
    ## implementation creates a fully paired dataset with nearest
    ## matching neighbours when pairs are missing.
    
    ## compute distance matrix for NNM-pairing
    ##y <- factor(as.character(y))
    y <- paste0("y=",y)
    dX <- scale(X)
    if(center.x) dX <- dX - rowMeans(dX,na.rm=TRUE)
    if(center.m) {
        y1 <- paste0("y.",y)
        mX <- tapply(1:ncol(dX),y1,function(i) rowMeans(dX[,i,drop=FALSE]))
        mX <- do.call(cbind, mX)
        dX <- dX - mX[,y1]
    }

    if(dist.method=="cor") {
        D <- 1 - cor(dX)
    } else {
        D <- as.matrix(dist(t(dX)))
    }
    remove(dX)

    if(b.method=="old") {
        ##replace=TRUE
        ##replace=FALSE
        maxD <- max(D)
        grp <- sort(unique(y))
        ny <- length(grp)
        B <- matrix(NA,nrow=ncol(X),ncol=ny)
        rownames(B) <- colnames(X)
        colnames(B) <- grp
        i=1;k=2
        for(i in 1:nrow(B)) {
            for(k in 1:ny) {
                jj <- which(y == grp[k])
                ii <- which(y == y[i])
                dd <- D[i,jj]
                if(!replace) {
                    dd <- dd + maxD*table(factor(B[ii,k],levels=jj)) ## penalty factor???
                }
                B[i,k] <- jj[which.min(dd)]
            }
        }
        B <- apply(B,2,function(i)rownames(B)[i])
    } else {
        ## should be faster...
        B <- t(apply(D,1,function(r) tapply(r,y1,function(s)names(which.min(s)))))
    }
    rownames(B) <- colnames(X)
    head(B)
    
    ## imputing full paired data set
    kk <- match(as.vector(B), rownames(B))
    full.y <- y[kk]
    full.pairs <- rep(rownames(B),ncol(B))
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
    full.idx <- rownames(B)[kk]
    cX <- do.call(cbind, tapply(1:ncol(full.X), full.idx,
                                function(i) rowMeans(full.X[,i,drop=FALSE])))
    cX <- cX[,colnames(X)]

    dim(cX)
    res <- list(X=cX, pairings=B)
    return(res)
}


pairs=NULL;dist.method="cor"
gx.nnmcorrect2.NOTWORKING <- function(X, y, pairs=NULL, use.design=TRUE,
                                      center.x=TRUE, center.m=TRUE,
                                      dist.method="cor")
{
    ## Nearest-neighbour matching for batch correction. This
    ## implementation creates a fully paired dataset with nearest
    ## matching neighbours when pairs are missing.
    require(sva)
    require(limma)
    y1  <- paste0("y=",y)
    dX  <- scale(X)
    if(center.x)
        dX <- dX - rowMeans(dX,na.rm=TRUE)    
    if(center.m) {
        mX <- tapply(1:ncol(dX), y1, function(i) rowMeans(dX[,i,drop=FALSE],na.rm=TRUE))
        mX <- do.call(cbind,mX)
        dX <- dX - mX[,y1]
    }
    if(dist.method=="cor") {
        D <- 1 - cor(dX)
    } else {
        D <- as.matrix(dist(t(dX)))
    }
    B <- t(apply(D,1,function(r) tapply(r,y1,function(s)names(which.min(s)))))

    if(is.null(pairs)) {
        pairs <- rep(NA,length(y))
        pairs[which(y==y[1])] <- paste0("s",1:sum(y==y[1])) ## at least one group
    }
    pairs <- as.character(pairs)
    names(pairs) <- colnames(X)
    nnb.idx <- matrix(pairs[as.vector(B)],ncol=ncol(B))
    nnb.best <- apply(nnb.idx,1,function(s) table(s))
    nnb.best <- apply(nnb.idx,1,function(s) names(which.max(table(s))))
    nnb.best <- sapply(apply(nnb.idx,1,function(s) {
        names(which(table(s)==max(table(s)))) }),sample,1)
    names(nnb.best) <- colnames(X)
    nnb.best

    ## !!! here not working: some are not fully paired up and correct to zero...
    full.pairs <- pairs
    full.pairs[is.na(pairs)] <- nnb.best[is.na(pairs)]
    full.pairs
    
    ## batch correction on full-pair info
    mod = NULL
    if(use.design) {
        mod <- model.matrix(~y)
        ##cX <- sva::ComBat(X, batch=full.pairs, mod=mod)        
        cX <- limma::removeBatchEffect(X, batch=full.pairs, design=mod)
    } else {
        ##cX <- sva::ComBat(X, batch=full.pairs)        
        cX <- limma::removeBatchEffect(X, batch=full.pairs)
    }
    cX
    res <- list(X=cX, pairings=B)
    return(res)
}

gx.nnmcorrect.OLD <- function(x, y, k=3, dist.method="cor")
{
    ## Nearest-neighbour matching for batch correction. This
    ## implementation substracts explicitly the nearest matching
    ## neighbour. It does not work well with multiple groups.

    
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


