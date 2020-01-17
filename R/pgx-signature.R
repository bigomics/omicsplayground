
ntop=1000;nsig=100;nperm=10000
pgx.correlateSignatureH5 <- function(fc, h5.file, nsig=100, ntop=1000, nperm=10000)
{
    ##
    ##
    ##
    ##
    
    if(is.null(names(fc))) stop("fc must have names")
    
    ## or instead compute correlation on top100 fc genes (read from file)
    rn <- h5read(h5.file,"data/rownames")
    cn <- h5read(h5.file,"data/colnames")

    ## ---------------------------------------------------------------
    ## Compute simple correlation between query profile and signatures
    ## ---------------------------------------------------------------
    gg <- names(c(head(sort(fc),nsig), head(tail(fc),nsig)))
    ##gg <- intersect(names(fc),rn)
    gg <- intersect(gg,rn)
    row.idx <- match(gg,rn)
    G <- h5read(h5.file, "data/matrix", index=list(row.idx,1:length(cn)))
    dim(G)
    dimnames(G) <- list(rn[row.idx],cn)
    rho <- cor( G[gg,], fc[gg], use="pairwise")[,1]

    ## --------------------------------------------------
    ## test all signature on query profile using fGSEA
    ## --------------------------------------------------
    
    require(fgsea)
    sel <- head(names(sort(-abs(rho))), ntop)
    sel.idx <- match(sel, cn)
    sig100.up <- h5read(h5.file, "signature/sig100.up",
                        index = list(1:100, sel.idx) )
    sig100.dn <- h5read(h5.file, "signature/sig100.dn",
                        index = list(1:100, sel.idx) )                        
    ##head(sig100.up,2)    

    ## combine up/down into one
    gmt <- rbind(sig100.up, sig100.dn)
    gmt <- unlist(apply(gmt, 2, list),recursive=FALSE)
    names(gmt) <- cn[sel.idx]
    length(gmt)
    
    ##system.time( res <- fgsea(gmt, fc, nperm=10000))
    system.time( res <- fgsea(gmt, abs(fc), nperm=nperm))
    dim(res)
            
    ## ---------------------------------------------------------------
    ## Combine correlation+GSEA by combined score (NES*rho)
    ## ---------------------------------------------------------------
    jj <- match( res$pathway, names(rho))
    res$rho <- rho[jj]
    res$R2 <- rho[jj]**2
    res$R2.NES <- res$R2*res$NES
    res <- res[order(res$R2.NES, decreasing=TRUE),]

    if(0) {
        res$rho.p <- cor.pvalue(res$rho, n=length(gg))
        res$meta.p  <- apply( res[,c("pval","rho.p")], 1, function(p) sumz(p)$p)    
        res <- res[order(res$meta.p),]
    }
    
    head(res)
    return(res)
}

pgx.createSignatureH5 <- function(pgx.files, h5.file, chunk=100, pgx.names=NULL)
{
    require(rhdf5)

    ##pgx0 <- dir(pgx.dir, "[.]pgx$", full.names=FALSE)
    ##pgx <- dir(pgx.dir, "[.]pgx$", full.names=TRUE)
    if(is.null(pgx.names) && !is.null(names(pgx.files))) {
        pgx.names <- names(pgx.files)
    }
    if(is.null(pgx.names)) {
        pgx.names <- pgx.files
    }

    ##--------------------------------------------------
    ## make big FC signature matrix
    ##--------------------------------------------------
    F <- list()
    cat("reading FC from",length(pgx.files),"pgx files ")
    i=1
    for(i in 1:length(pgx.files)) {
        if(!file.exists(pgx.files[i])) next()
        cat(".")
        load(pgx.files[i], verbose=0)
        meta <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        rownames(meta$fc) <- toupper(rownames(meta$fc))
        pgx <- gsub(".*[/]|[.]pgx$","",pgx.names[i])
        colnames(meta$fc) <- paste0("[",pgx,"] ",colnames(meta$fc))
        F[[ pgx ]] <- meta$fc    
    }
    cat("\n")
    
    genes <- sort(unique(as.vector(sapply(F,rownames))))
    length(genes)    
    F <- lapply(F, function(x) x[match(genes,rownames(x)),,drop=FALSE])
    X <- do.call(cbind, F)
    rownames(X) <- genes    
    pgx.saveMatrixH5(X, h5.file, chunk=chunk)

    ##--------------------------------------------------
    ## Calculate top100 gene signatures
    ##--------------------------------------------------
    cat("Creating signatures...\n")

    ## X  <- h5read(h5.file, "data/matrix")
    rn <- h5read(h5.file,"data/rownames")
    cn <- h5read(h5.file,"data/colnames")
    h5ls(h5.file)

    dim(X)
    ##X <- X[,1:100]
    X[is.na(X)] <- 0
    orderx <- apply(X,2,function(x) {
        idx=order(x);
        list(DN=head(idx,100),UP=rev(tail(idx,100)))
    })    
    sig100.dn <- sapply(orderx,"[[","DN")
    sig100.dn <- apply(sig100.dn, 2, function(i) rn[i])
    sig100.up <- sapply(orderx,"[[","UP")
    sig100.up <- apply(sig100.up, 2, function(i) rn[i])

    h5createGroup(h5.file,"signature")    
    h5write( sig100.dn, h5.file, "signature/sig100.dn")  ## can write list???    
    h5write( sig100.up, h5.file, "signature/sig100.up")  ## can write list??

    ##--------------------------------------------------
    ## Precalculate t-SNE/UMAP
    ##--------------------------------------------------
    dim(X)
    h5createGroup(h5.file,"clustering")    

    cat("calculating t-SNE 2D...\n")
    require(Rtsne)
    system.time(
        pos <- Rtsne(
            t(X[,]),
            dims = 2,
            check_duplicates=FALSE )$Y
    )
    ## rownames(pos) <- colnames(X)[1:nrow(pos)]  ## not recorded anyway...
    h5write( pos, h5.file, "clustering/tsne2d")  ## can write list??    

    if(1) {
        cat("calculating t-SNE 3D...\n")
        require(Rtsne)
        system.time(
            pos <- Rtsne(
                t(X[,]),
                dims = 3,
                check_duplicates=FALSE )$Y
        )
        ## rownames(pos) <- colnames(X)[1:nrow(pos)]
        h5write( pos, h5.file, "clustering/tsne3d")  ## can write list??    
    }
    
    if(1) {
        cat("calculating UMAP...\n")
        require(umap)
        system.time(
            pos <- umap( t(X[,]) )$layout
        )
        ## rownames(pos) <- colnames(X)[1:nrow(pos)]
        h5write( pos, h5.file, "clustering/umap")  ## can write list??    
    }

    h5closeAll()
    ## return(X)
}

##-------------------------------------------------------------------
## Pre-calculate geneset expression with different methods
##-------------------------------------------------------------------

pgx.computeMultiOmicsGSE <- function(X, gmt, omx.type, 
                                     method=NULL, center=TRUE)
{
    if(0) {
        omx.type <- c("MRNA","MIR")[1+grepl("^MIR",rownames(X))]
        table(omx.type)
        omx.type <- sample(c("MRNA","CNV"),nrow(X),replace=TRUE)
    }
    if(is.null(omx.type))
        omx.type <- gsub("[:=].*","",rownames(X))
    omx.types <- setdiff(unique(omx.type),c("MIR",""))
    omx.types

    sx <- list()
    for(tp in omx.types) {
        x1 <- X[which(omx.type==tp),]
        rownames(x1) <- sub(":.*","",rownames(x1))
        sx[[tp]] <- pgx.computeGeneSetExpression(x1, gmt, method=method, center=center)
        sx[[tp]] <- lapply(sx[[tp]],function(x) {
            rownames(x)=paste0(tp,"=",rownames(x))
            x
        })
    }

    ## concatenate all omx-types
    cx <- sx[[1]]
    for(j in 1:length(sx[[1]])) {
        cx[[j]] <- do.call(rbind, lapply(sx,"[[",j))
    }
    
    return(cx)
}

pgx.computeGeneSetExpression <- function(X, gmt, method=NULL,
                                         min.size=10, center=TRUE)
{    
    library(GSVA)
    ALL.METHODS <- c("gsva","spearman","average")
    ALL.METHODS <- c("gsva","ssgsea","spearman","average")
    if(is.null(method))
        method <- ALL.METHODS
    if(0){
        X=ngs$X;gmt=GSETS[grep("HALLMARK",names(GSETS))]
    }
    ## this is important!!! centering on genes (GSVA does)
    if(center) {
        X <- X - rowMeans(X,na.rm=TRUE)
    }
    dim(X)
    
    gmt.size <- sapply(gmt, function(x) sum(x %in% rownames(X)))
    gmt <- gmt[ gmt.size >= min.size ]
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
    S1 <- lapply(S,function(x) apply(x,2,rank)) ## rank by sample
    S[["meta"]] <- scale(Reduce('+',S1)/length(S1))   
    gs <- Reduce(intersect, lapply(S,rownames)) 
    S <- lapply(S, function(x) x[gs,])
    
    if(0) {
        ## show pairs
        names(S)
        dim(S[[1]])
        pairs(sapply(S,function(x) x[,1])) ## corr by genesets
        pairs(sapply(S,function(x) x[1,])) ## corr by sample       
    }
            
    return(S)
}


