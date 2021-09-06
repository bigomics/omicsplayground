##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##========================================================================
##================ rank correlation based geneset testing ================
##========================================================================

gset.cor <- function(rnk, gset, compute.p=FALSE) {
    gset.rankcor(rnk=rnk, gset=gset, compute.p=compute.p, use.rank=FALSE)
}

gset.rankcor <- function(rnk, gset, compute.p=FALSE, use.rank=TRUE)
{
    if(!any(class(gset) %in% c("Matrix","dgCMatrix","matrix","array")) ) {
        stop("gset must be a matrix")
        ##gset = gmt2mat(gset, bg=names(rnk))
    }
    is.vec <- (NCOL(rnk)==1 && !class(rnk) %in% c('matrix','array','Matrix'))
    if(is.vec && is.null(names(rnk))) stop("rank vector must be named")
    if(!is.vec && is.null(rownames(rnk))) stop("rank matrix must have rownames")
    if(is.vec) rnk <- matrix(rnk,ncol=1,dimnames=list(names(rnk),'rnk'))
    dim(rnk)
    Matrix::head(rnk)
    n1 <- sum(rownames(rnk) %in% colnames(gset), na.rm=TRUE)
    n2 <- sum(rownames(rnk) %in% rownames(gset), na.rm=TRUE)
    if( n1 > n2 ) gset <- t(gset)
    
    gg <- intersect(rownames(gset),rownames(rnk))
    rnk1 <- rnk[gg,,drop=FALSE]
    if(use.rank) {
        rnk1 <- apply(rnk1, 2, rank )
    }
    rho1 <- qlcMatrix::corSparse(gset[gg,], rnk1)
    rownames(rho1) <- colnames(gset)
    colnames(rho1) <- colnames(rnk1)
    
    ## compute p-value by permutation
    pv=qv=NULL
    if(compute.p) {
        ## Permutation of the GSET matrix
        idx = which(gset!=0, arr.ind=TRUE)
        dim(idx)
        S = Matrix::sparseMatrix(sample(idx[,1]), sample(idx[,2]), x=rep(1,nrow(idx)),
                         dims = dim(gset), dimnames=dimnames(gset) )
        rho2 = qlcMatrix::corSparse(S[gg,], rnk1)
        rho1[is.na(rho1)] <- 0
        z1 = abs(rho1)/sd(rho2,na.rm=TRUE)
        pv = 2*pnorm(-z1)
        qv = apply(pv, 2, p.adjust, method='fdr')
        df = list(rho=rho1, p.value=pv, q.value=qv)
    } else {
        df = list(rho=rho1, p.value=NA, q.value=NA)
    }    
    df
}
