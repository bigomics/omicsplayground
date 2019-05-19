##-------------------------------------------------------------------
## Pre-calculate geneset expression with different methods
##-------------------------------------------------------------------
pgx.computeGeneSetExpression <- function(X, gmt, method=NULL, center=TRUE) {
    
    library(GSVA)
    ALL.METHODS <- c("gsva","spearman","average")
    ALL.METHODS <- c("gsva","ssgsea","spearman","average")
    if(is.null(method))
        method <- ALL.METHODS

    ##X=ngs$X;gmt=GSETS[grep("HALLMARK",names(GSETS))]
    ## this is important!!! centering on genes (GSVA does)
    if(center) {
        X <- X - rowMeans(X,na.rm=TRUE)
    }
    dim(X)
    
    gmt.size <- sapply(gmt, function(x) sum(x %in% rownames(X)))
    gmt <- gmt[ gmt.size>=1 ]
    length(gmt)
    
    S <- list()
    if("gsva" %in% method) {
        S[["gsva"]] <- gsva(X, gmt, method="gsva")
    }
    if("ssgsea" %in% method) {
        S[["ssgsea"]] <- gsva(X, gmt, method="ssgsea", min.sz=1)
    }
    if(any(method %in% c("spearman","average"))) {
        gg <- rownames(X)
        G <- gmt2mat(gmt, bg=gg)
        if("spearman" %in% method) {
            ##rho <- cor(as.matrix(G[gg,]), apply(X[gg,],2,rank))
            rho <- t(G[gg,]) %*% scale(apply(X[gg,],2,rank)) / sqrt(nrow(X)-1)
            rho[is.na(rho)] <- 0
            S[["spearman"]] <- rho
        }
        if("average" %in% method) {
            ##rho <- cor(as.matrix(G[gg,]), apply(G[gg,],2,rank))
            avg.X <- t(G[gg,]) %*% X[gg,] / Matrix::colSums(G[gg,])
            avg.X[is.na(avg.X)] <- 0
            S[["average"]] <- avg.X
        }        
    }

    ## compute meta score
    S1 <- lapply(S,function(x) apply(x,2,rank))
    S[["meta"]] <- scale(Reduce('+',S1)/length(S1))   
    gs <- Reduce(intersect, lapply(S,rownames)) 
    S <- lapply(S, function(x) x[gs,])
    
    if(0) {
        names(S)
        pairs(sapply(S,function(x) x[,1]))
    }
            
    return(S)
}


