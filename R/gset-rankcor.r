##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##========================================================================
##================ rank correlation based geneset testing ================
##========================================================================


gset.rankcor <- function(rnk, gset , compute.p=FALSE)
{
    if(!any(class(gset) %in% c("Matrix","dgCMatrix")) ) {
        stop("gset must be a matrix")
        ##gset = gmt2mat(gset, bg=names(rnk))
    }
    gg <- intersect(rownames(gset),names(rnk))
    rnk1 <- matrix(rank(rnk[gg]),ncol=1)
    rho1 <- qlcMatrix::corSparse( gset[gg,], rnk1)[,1]
    names(rho1) <- colnames(gset)
    
    ## compute p-value by permutation
    pv=qv=NULL
    if(compute.p) {

        if(0) {
            ## Permutation of the ranking vector
            R <- sapply(1:1000, function(i) sample(rnk1))
            rho2 <- qlcMatrix::corSparse( gset[gg,], R)
            pv <- rowMeans( abs(rho2) >  abs(rho1))
            pv = pmax(pv,1/ncol(R))

        } else {

            ## Permutation of the GSET matrix
            idx = which(gset!=0, arr.ind=TRUE)
            dim(idx)
            S = sparseMatrix(sample(idx[,1]), sample(idx[,2]), x=rep(1,nrow(idx)),
                             dims = dim(gset), dimnames=dimnames(gset) )
            rho2 = qlcMatrix::corSparse(S[gg,], rnk1)[,1]
            ## hist(rho2,breaks=200)        
            pv = sapply( rho1, function(r) mean( abs(rho2) > abs(r),na.rm=TRUE))
            pv = pmax(pv,1/ncol(S))
            names(pv) = names(rho1)
        }
        qv = p.adjust(pv, method='fdr')
        df = data.frame( rho=rho1, p.value=pv, q.value=qv)
        df = df[order(df$p.value,-abs(df$rho)),]
    } else {
        df = data.frame( rho=rho1, p.value=NA, q.value=NA)
        df = df[order(-abs(df$rho)),]
    }    
    df
}
