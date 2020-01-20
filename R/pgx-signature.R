if(0) {

    sigdb = "../data/datasets-allFC.csv"
    h5.file = "../data/sigdb-gse25k.h5"
    FILES="../lib"
    RDIR="../R"
    source("../R/pgx-include.R")
    source("../R/pgx-files.R")
}

chunk=100
pgx.createSignatureDatabaseH5 <- function(pgx.files, h5.file, chunk=100, update.only=FALSE)
{
    require(rhdf5)

    h5exists <- function(h5.file, obj) {
        xobjs <- apply(h5ls(h5.file)[,1:2],1,paste,collapse="/")
        obj %in% gsub("^/|^//","",xobjs)
    }

    if(update.only && h5exists(h5.file, "data/matrix")) {
        X  <- h5read(h5.file, "data/matrix")
        rn <- h5read(h5.file,"data/rownames")
        cn <- h5read(h5.file,"data/colnames")
        rownames(X) <- rn
        colnames(X) <- cn
    } else {
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
            rownames(meta$fc) <- toupper(rownames(meta$fc))  ## mouse-friendly
            pgx <- gsub(".*[/]|[.]pgx$","",pgx.files[i])
            colnames(meta$fc) <- paste0("[",pgx,"] ",colnames(meta$fc))
            F[[ pgx ]] <- meta$fc    
        }
        cat("\n")
        
        genes <- as.vector(unlist(sapply(F,rownames)))
        genes <- sort(unique(toupper(genes)))
        length(genes)    
        F <- lapply(F, function(x) x[match(genes,rownames(x)),,drop=FALSE])
        X <- do.call(cbind, F)
        rownames(X) <- genes    

        ## Filter out genes (not on known chromosomes...)
        genes <- rownames(X)
        gannot <- ngs.getGeneAnnotation(genes)
        table(is.na(gannot$chr))
        sel <- which(!is.na(gannot$chr))
        X <- X[sel,,drop=FALSE]
        dim(X)

        pgx.saveMatrixH5(X, h5.file, chunk=chunk)

        if(0) {
            h5ls(h5.file)
            h5write( X, h5.file, "data/matrix")  ## can write list??
            h5write( colnames(X), h5.file,"data/colnames")
            h5write( rownames(X), h5.file,"data/rownames")
        }        
        remove(F)
    }
    dim(X)
    
    ##--------------------------------------------------
    ## Calculate top100 gene signatures
    ##--------------------------------------------------
    cat("Creating signatures...\n")
    
    if(!update.only || !h5exists(h5.file, "signature")) {
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
        
        if(!h5exists(h5.file, "signature")) h5createGroup(h5.file,"signature")    
        h5write( sig100.dn, h5.file, "signature/sig100.dn")  ## can write list???    
        h5write( sig100.up, h5.file, "signature/sig100.up")  ## can write list??
        
        remove(orderx)
        remove(sig100.dn)
        remove(sig100.up)
    }
    
    ##--------------------------------------------------
    ## Precalculate t-SNE/UMAP
    ##--------------------------------------------------
    dim(X)

    if(!update.only || !h5exists(h5.file, "clustering")) {
        
        if(!h5exists(h5.file, "clustering")) h5createGroup(h5.file,"clustering")    
        h5ls(h5.file)
        
        pos <- pgx.clusterBigMatrix(
            abs(X),  ## on absolute foldchange!!
            methods=c("pca","tsne","umap"),
            dims=c(2,3),
            reduce.sd = 2000,
            reduce.pca = 200 )
        names(pos)
        
        h5write( pos[["pca2d"]], h5.file, "clustering/pca2d")  ## can write list??    
        h5write( pos[["pca3d"]], h5.file, "clustering/pca3d")  ## can write list??    
        h5write( pos[["tsne2d"]], h5.file, "clustering/tsne2d")  ## can write list??    
        h5write( pos[["tsne3d"]], h5.file, "clustering/tsne3d")  ## can write list??    
        h5write( pos[["umap2d"]], h5.file, "clustering/umap2d")  ## can write list??    
        h5write( pos[["umap3d"]], h5.file, "clustering/umap3d")  ## can write list??            

    }

    h5closeAll()
    ## return(X)
}

mc.cores=8
pgx.addEnrichmentSignaturesH5 <- function(h5.file, X=NULL, mc.cores=4, lib.dir) 
{

    require(rhdf5)
    
    h5exists <- function(h5.file, obj) {
        xobjs <- apply(h5ls(h5.file)[,1:2],1,paste,collapse="/")
        obj %in% gsub("^/|^//","",xobjs)
    }

    if(is.null(X)) {
        X  <- h5read(h5.file, "data/matrix")
        rn <- h5read(h5.file,"data/rownames")
        cn <- h5read(h5.file,"data/colnames")
        rownames(X) <- rn
        colnames(X) <- cn
    }

    ##sig100.dn <- h5read(h5.file, "signature/sig100.dn")  
    ##sig100.up <- h5read(h5.file, "signature/sig100.up")  
    
    G <- readRDS(file.path(lib.dir,"gset-sparseG-XL.rds"))
    dim(G)    
    sel <- grep("HALLMARK|C[1-9]|^GO", rownames(G))
    sel <- grep("HALLMARK", rownames(G))
    length(sel)
    G <- G[sel,,drop=FALSE]
    gmt <- apply( G, 1, function(x) colnames(G)[which(x!=0)])

    ##X <- X[,1:20]
    
    require(fgsea)
    ##F1 <- apply(X, 2, function(x) {fgsea( gmt, x, nperm=1000)$NES })
    F1 <- apply(X, 2, function(x) {fgsea( gmt, x, nperm=10)$NES })  ## FDR not important, small nperm
    rownames(F1) <- rownames(G)

    require(GSVA)
    ## mc.cores = 4
    F2 <- gsva(X, gmt, method="gsva", parallel.sz=mc.cores)
    dim(F2)
    
    if(!h5exists(h5.file, "enrichment")) h5createGroup(h5.file,"enrichment")
    h5write(rownames(F1), h5.file, "enrichment/genesets")
    h5write(F1, h5.file, "enrichment/GSEA")
    h5write(F2, h5.file, "enrichment/GSVA")

    h5ls(h5.file)
    h5closeAll()

}

pgx.computeConnectivityScores <- function(ngs, sigdb, ntop=-1, contrasts=NULL)
{

    meta = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
    colnames(meta$fc)
    
    is.h5ref <- grepl("h5$",sigdb)       
    ##cat("[calcConnectivityScores] sigdb =",sigdb,"\n")
    ##cat("[calcConnectivityScores] ntop =",ntop,"\n")
    h5.file <- NULL
    refmat <- NULL
    if(grepl("csv$",sigdb)) {
        refmat <- read.csv(sigdb,row.names=1,check.names=FALSE)
        dim(refmat)
    }
    if(grepl("h5$",sigdb)) {
        if(file.exists(sigdb)) h5.file <- sigdb
    }
    
    if(is.null(contrasts))
        contrasts <- colnames(meta$fc)
    contrasts <- intersect(contrasts, colnames(meta$fc))

    scores <- list()
    ct <- contrasts[1]
    for(ct in contrasts) {
        
        fc <- meta$fc[,ct]
        names(fc) <- rownames(meta$fc)
        names(fc) <- toupper(names(fc)) ## for mouse
        
        h5.file
        if(!is.null(h5.file))  {
            res <- pgx.correlateSignatureH5(
                fc, h5.file = h5.file,
                nsig=100, ntop=ntop, nperm=10000)            

        } else if(!is.null(refmat)) {                
            res <- pgx.correlateSignature(
                fc, refmat = refmat,
                nsig=100, ntop=ntop, nperm=10000)
            
        } else {
            stop("FATAL:: could not determine reference type")
        }
        dim(res)
        scores[[ct]] <- res
    }

    names(scores)
    return(scores)
}


## ntop=1000;nsig=100;nperm=10000
pgx.correlateSignature <- function(fc, refmat, nsig=100, ntop=1000, nperm=10000)
{
    ##
    ##
    ##
    ##
    
    if(is.null(names(fc))) stop("fc must have names")

    ## mouse... mouse...
    names(fc) <- toupper(names(fc))
    
    ## or instead compute correlation on top100 fc genes (read from file)
    ##refmat = PROFILES$FC
    rn <- rownames(refmat)
    cn <- colnames(refmat)
    
    ## ---------------------------------------------------------------
    ## Compute simple correlation between query profile and signatures
    ## ---------------------------------------------------------------
    fc <- sort(fc)
    gg <- unique(names(c(head(fc,nsig), tail(fc,nsig))))
    ##gg <- intersect(names(fc),rn)
    gg <- intersect(gg,rn)
    G  <- refmat[gg,,drop=FALSE]
    dim(G)
    rho <- cor( G[gg,], fc[gg], use="pairwise")[,1]
    names(rho) <- colnames(G)
    
    ## --------------------------------------------------
    ## test all signature on query profile using fGSEA
    ## --------------------------------------------------
    
    require(fgsea)
    sel <- head(names(sort(-abs(rho))),ntop)

    notx <- setdiff(sel,colnames(refmat))
    if(length(notx)>0) {
        ## should not happen...   
        cat("[pgx.correlateSignature] length(sel)=",length(sel),"\n")
        cat("[pgx.correlateSignature] head(sel)=",head(sel),"\n")
        cat("[pgx.correlateSignature] head.notx=",head(notx),"\n")
    }

    sel <- intersect(sel, colnames(refmat))  
    X <- refmat[,sel,drop=FALSE]
    dim(X)
    X[is.na(X)] <- 0
    orderx <- apply(X,2,function(x) {
        idx=order(x);
        list(DN=head(idx,100),UP=rev(tail(idx,100)))
    })    
    sig100.dn <- sapply(orderx,"[[","DN")
    sig100.dn <- apply(sig100.dn, 2, function(i) rn[i])
    sig100.up <- sapply(orderx,"[[","UP")
    sig100.up <- apply(sig100.up, 2, function(i) rn[i])
    dim(sig100.dn)
    
    ## combine up/down into one (unsigned GSEA test)
    gmt <- rbind(sig100.up, sig100.dn)
    gmt <- unlist(apply(gmt, 2, list),recursive=FALSE)
    names(gmt) <- colnames(X)
    length(gmt)
    
    ##system.time( res <- fgsea(gmt, fc, nperm=10000))
    suppressMessages( suppressWarnings(
        res <- fgsea(gmt, abs(fc), nperm=nperm)
    ))
    dim(res)
            
    ## ---------------------------------------------------------------
    ## Combine correlation+GSEA by combined score (NES*rho)
    ## ---------------------------------------------------------------
    jj <- match(res$pathway, names(rho))
    res$rho <- rho[jj]
    res$R2 <- rho[jj]**2
    res$score <- res$R2*res$NES
    res <- res[order(res$score, decreasing=TRUE),]

    if(0) {
        res$rho.p <- cor.pvalue(res$rho, n=length(gg))
        res$meta.p  <- apply( res[,c("pval","rho.p")], 1, function(p) sumz(p)$p)    
        res <- res[order(res$meta.p),]
    }
    
    head(res)
    return(res)
}

ntop=1000;nsig=100;nperm=10000
pgx.correlateSignatureH5 <- function(fc, h5.file, nsig=100, ntop=1000, nperm=10000)
{
    ##
    ##
    ##
    ##
    
    if(is.null(names(fc))) stop("fc must have names")
    
    ## mouse... mouse...
    names(fc) <- toupper(names(fc))

    ## or instead compute correlation on top100 fc genes (read from file)
    rn <- h5read(h5.file,"data/rownames")
    cn <- h5read(h5.file,"data/colnames")

    ## ---------------------------------------------------------------
    ## Compute simple correlation between query profile and signatures
    ## ---------------------------------------------------------------
    fc <- sort(fc)
    gg <- unique(names(c(head(fc,nsig), tail(fc,nsig))))
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

    ## combine up/down into one (unsigned GSEA test)
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
    res$score <- res$R2*res$NES
    res <- res[order(res$score, decreasing=TRUE),]

    if(0) {
        res$rho.p <- cor.pvalue(res$rho, n=length(gg))
        res$meta.p  <- apply( res[,c("pval","rho.p")], 1, function(p) sumz(p)$p)    
        res <- res[order(res$meta.p),]
    }
    
    head(res)
    return(res)
}

pgx.ReclusterSignatureDatabase <- function(h5.file, reduce.sd=1000, reduce.pca=100)
{
    require(rhdf5)

    h5exists <- function(h5.file, obj) {
        xobjs <- apply(h5ls(h5.file)[,1:2],1,paste,collapse="/")
        obj %in% gsub("^/|^//","",xobjs)
    }
    
    X  <- h5read(h5.file, "data/matrix")
    rn <- h5read(h5.file,"data/rownames")
    cn <- h5read(h5.file,"data/colnames")
    rownames(X) <- rn
    colnames(X) <- cn

    ##--------------------------------------------------
    ## Precalculate t-SNE/UMAP
    ##--------------------------------------------------
    dim(X)
    
    if(!h5exists(h5.file, "clustering")) h5createGroup(h5.file,"clustering")    
    
    pos <- pgx.clusterBigMatrix(
        abs(X),  ## on absolute foldchange!!
        methods = c("pca","tsne","umap"),
        dims = c(2,3),
        reduce.sd = reduce.sd,
        reduce.pca = reduce.pca )
    names(pos)
    
    h5write( pos[["pca2d"]], h5.file, "clustering/pca2d")  ## can write list??    
    h5write( pos[["pca3d"]], h5.file, "clustering/pca3d")  ## can write list??    
    h5write( pos[["tsne2d"]], h5.file, "clustering/tsne2d")  ## can write list??    
    h5write( pos[["tsne3d"]], h5.file, "clustering/tsne3d")  ## can write list??    
    h5write( pos[["umap2d"]], h5.file, "clustering/umap2d")  ## can write list??    
    h5write( pos[["umap3d"]], h5.file, "clustering/umap3d")  ## can write list??            
    h5closeAll()    
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


