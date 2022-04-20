##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##X=comboX;x.drugs=combo.drugs
##X=comboX;x.drugs=combo.drugs
##nprune=250;res.mono=NULL;ntop=10;nsample=20

pgx.createComboDrugAnnot <- function(combo, annot0) {
    drugs <- strsplit(combo,split="[+]")
    drug1 <- sapply(drugs,"[",1)
    drug2 <- sapply(drugs,"[",2)
    j1 <- match( drug1, rownames(annot0))
    j2 <- match( drug2, rownames(annot0))
    cmoa <- paste( annot0[j1,"moa"],"+",annot0[j2,"moa"])
    ctarget <- paste( annot0[j1,"target"],"+",annot0[j2,"target"])
    cmoa <- gsub("NA [+] NA","",cmoa)
    ctarget <- gsub("NA [+] NA","",ctarget)
    annot1 <- data.frame( drug=combo, moa=cmoa, target=ctarget)
    rownames(annot1) <- combo
    annot1
}

##obj=ngs;methods="GSEA";contrast=NULL;nprune=250;nmin=5
pgx.computeDrugEnrichment <- function(obj, X, xdrugs, methods=c("GSEA","cor"),
                                      nmin=15, nprune=250, contrast=NULL )
{
    ## 'obj'   : can be ngs object or fold-change matrix
    ## X       : drugs profiles (may have multiple for one drug)
    ## xdrugs : drug associated with profile
    
    names(obj)
    if("gx.meta" %in% names(obj)) {
        F <- pgx.getMetaMatrix(obj)$fc
        ## check if multi-omics
        is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]",rownames(F)))
        is.multiomics
        if(is.multiomics) {
            jj <- grep("\\[gx\\]|\\[mrna\\]",rownames(F))
            F <- F[jj,,drop=FALSE]
        }
        rownames(F) <- toupper(sub(".*:|.*\\]","",rownames(F)))
        ##colnames(F) <- names(obj$gx.meta$meta)
        F <- F[order(-rowMeans(F**2)),,drop=FALSE]
        F <- F[!duplicated(rownames(F)),,drop=FALSE]
        dim(F)
    } else {
        ## it is a matrix
        F <- obj
    }

    if(is.null(contrast))
        contrast <- colnames(F)
    contrast <- intersect(contrast, colnames(F))
    contrast
    F <- F[,contrast,drop=FALSE]

    ## create drug meta sets
    meta.gmt <- tapply(colnames(X), xdrugs, list)
    meta.gmt <- meta.gmt[which(sapply(meta.gmt,length)>=nmin)]
    length(meta.gmt)
    if(length(meta.gmt)==0) {
        message("WARNING::: pgx.computeDrugEnrichment : no valid genesets!!")
        return(NULL)
    }
    
    ## first level (rank) correlation
    message("Calculating first level rank correlation ...")
    gg <- intersect(rownames(X), rownames(F))
    length(gg)
    if(length(gg) < 20) {
        message("WARNING::: pgx.computeDrugEnrichment : not enough common genes!!")
        return(NULL)
    }
    rnk1 <- apply(X[gg,,drop=FALSE],2,rank,na.last="keep")
    rnk2 <- apply(F[gg,,drop=FALSE],2,rank,na.last="keep")
    system.time(R1 <- stats::cor(rnk1, rnk2, use="pairwise"))
    dim(R1)
    ##
    ##system.time(R1 <- gpuCor(rnk1, rnk2, use="pairwise")$coefficients)
    R1 <- R1 + 1e-8*matrix(rnorm(length(R1)),nrow(R1),ncol(R1))
    colnames(R1) <- colnames(F)
    rownames(R1) <- colnames(X)
    dim(R1)

    ## experiment to drug
    

    results <- list()
    if("cor" %in% methods) {
        message("Calculating drug enrichment using rank correlation ...")
        
        ##D <- model.matrix( ~ 0 + xdrugs)
        D <- Matrix::sparse.model.matrix( ~ 0 + xdrugs)
        dim(D)
        colnames(D) <- sub("^xdrugs","",colnames(D))
        rownames(D) <- colnames(X)          ## not necessary..
        rho2 <- qlcMatrix::corSparse(D, R1)
        rownames(rho2) <- colnames(D)
        colnames(rho2) <- colnames(R1)
        rho2 <- rho2[order(-rowMeans(rho2**2)),,drop=FALSE]
        cor.pvalue <- function(x,n) pnorm(-abs(x/((1-x**2)/(n-2))**0.5))
        P <- apply(rho2, 2, cor.pvalue, n=nrow(D))
        Q <- apply(P, 2, p.adjust, method="fdr")
        results[["cor"]] <- list( X=rho2, Q=Q, P=P)
    }

    if("GSEA" %in% methods) {
        message("Calculating drug enrichment using GSEA ...")
        res0 <- list()
        i=1
        for(i in 1:ncol(R1)) {
            suppressWarnings(res0[[i]] <- fgsea::fgseaSimple( meta.gmt, stats=R1[,i], nperm=1000))
        }
        names(res0) <- colnames(R1)
        length(res0)
        
        mNES <- sapply(res0, function(x) x$NES)
        mQ   <- sapply(res0, function(x) x$padj)
        mP   <- sapply(res0, function(x) x$pval)
        if(length(res0)==1) {
            mNES <- cbind(mNES)
            mP <- cbind(mP)
            mQ <- cbind(mQ)
        }
        
        pw <- res0[[1]]$pathway
        rownames(mNES) <- rownames(mQ) <- rownames(mP) <- pw
        colnames(mNES) <- colnames(mQ) <- colnames(mP) <- colnames(F)
        msize <- res0[[1]]$size
        dim(R1)        
        results[["GSEA"]] <- list( X=mNES, Q=mQ, P=mP, size=msize)
    }
    names(results)

    ## this takes only the top matching drugs for each comparison to
    ## reduce the size of the matrices
    nprune    
    if(nprune > 0) {
        message("[pgx.computeDrugEnrichment] pruning : nprune = ", nprune)        
        k=1
        for(k in 1:length(results)) {
            res <- results[[k]]
            ## reduce solution set with top-N of each comparison??
            mtop <- apply(abs(res$X), 2, function(x) Matrix::head(order(-x),nprune)) ## absolute NES!!
            mtop <- unique(as.vector(mtop))
            length(mtop)
            top.idx <- unique(unlist(meta.gmt[mtop]))
            length(top.idx)
            results[[k]]$X <- res$X[mtop,,drop=FALSE]
            results[[k]]$P <- res$P[mtop,,drop=FALSE]
            results[[k]]$Q <- res$Q[mtop,,drop=FALSE]
            results[[k]]$size <- res$size[mtop]
        }
    }

    ## UMAP clustering
    message("[pgx.computeDrugEnrichment] UMAP clustering...")
    top.drugs <- unique(unlist(sapply(results,function(res) rownames(res$X))))    
    sel <- which(xdrugs %in% top.drugs)
    ##sel <- unique(unlist(sapply(results,function(res) rownames(res$stats))))
    cX <- scale(X[,sel], center=FALSE)
    cX <- cX - rowMeans(cX)
    cX <- scale(cX, center=TRUE)
    pos <- uwot::umap(t(cX), fast_sgd=TRUE)
    dim(pos)
    rownames(pos) <- colnames(cX)
    results$clust <- pos
    results$stats <- R1[rownames(pos),,drop=FALSE]
    
    if(0) {
        pgx.scatterPlotXY.BASE(pos, var=2**R1[rownames(pos),1])
        pgx.scatterPlotXY.BASE(pos, var=log(colMeans(X**2)))
    }

    message("[pgx.computeDrugEnrichment] done!")
    
    return(results)
}


pgx.computeComboEnrichment <- function(obj, X, xdrugs,
                                       ntop=10, nsample=20, nprune=250,
                                       contrasts=NULL, res.mono=NULL )
{
    
    if(0) {
        X <- readRDS(file=file.path(FILES,"l1000_es.rds"))
        xdrugs <- gsub("_.*$","",colnames(X))
        length(table(xdrugs))
        dim(X)
        ntop=10;nsample=20
    }

    if("gx.meta" %in% names(obj)) {
        F <- sapply(obj$gx.meta$meta,function(x) x$meta.fx)
        rownames(F) <- rownames(obj$gx.meta$meta[[1]])
        ## check if multi-omics
        is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]",rownames(F)))
        is.multiomics
        if(is.multiomics) {
            jj <- grep("\\[gx\\]|\\[mrna\\]",rownames(F))
            F <- F[jj,,drop=FALSE]
        }
        rownames(F) <- toupper(sub(".*:|.*\\]","",rownames(F)))
        F <- F[order(-rowMeans(F**2)),,drop=FALSE]
        F <- F[!duplicated(rownames(F)),,drop=FALSE]
        dim(F)
    }

    if(is.null(contrasts)) contrasts <- colnames(F)
    contrasts <- intersect(contrasts, colnames(F))
    contrasts

    ## calculate average drug profile
    if(is.null(res.mono)) {
        cat("Calculating single drug enrichment using GSEA ...\n")
        er.mono <- pgx.computeDrugEnrichment(
            obj, X, xdrugs, methods="GSEA",
            nprune=nprune, contrasts=NULL )
        er.mono <- er.mono[["GSEA"]]
        names(er.mono)
    } else {
        cat("Using passed single drug enrichment results...\n")
        if("GSEA" %in% names(res.mono)) {
            er.mono <- res.mono[["GSEA"]]
        }
    }

    ## determine top-combinations
    ##ntop=10
    top.mono.up  <- apply(er.mono$X, 2, function(x) Matrix::head(order(-x),ntop))
    top.mono.dn  <- apply(er.mono$X, 2, function(x) Matrix::head(order(x),ntop))
    top.combo.up <- apply(top.mono.up,2,function(idx) list(combn(idx,2)))
    top.combo.dn <- apply(top.mono.dn,2,function(idx) list(combn(idx,2)))
    top.combo.up <- unlist(top.combo.up, recursive=FALSE)
    top.combo.dn <- unlist(top.combo.dn, recursive=FALSE)
    top.combo.up <- lapply(top.combo.up, function(d) apply(d,2,function(j) paste(sort(j),collapse="-")))
    top.combo.dn <- lapply(top.combo.dn, function(d) apply(d,2,function(j) paste(sort(j),collapse="-")))
    top.combo <- unique(c(unlist(top.combo.up),unlist(top.combo.dn)))
    length(top.combo)

    ##-------------- sample pairs from original mono-matrix
    sample.pairs <- list()
    k=1
    for(k in 1:length(top.combo)) {
        cmbn.idx <- as.integer(strsplit(top.combo[k],split="-")[[1]])
        cmbn.idx
        cmbn <- sort(rownames(er.mono$X)[cmbn.idx])
        Matrix::head(cmbn)
        p1 <- sample(which(xdrugs==cmbn[1]),nsample,replace=TRUE)
        p2 <- sample(which(xdrugs==cmbn[2]),nsample,replace=TRUE)
        pp <- cbind(p1, p2)
        sample.pairs[[k]] <- pp
    }
    sample.pairs <- do.call(rbind, sample.pairs)
    dim(sample.pairs)

    ##--------------- now create combination matrix X
    comboX <- apply(sample.pairs, 1, function(ii) rowMeans(X[,ii],na.rm=TRUE))
    dim(comboX)
    ##colnames(comboX) <- apply(combo.idx, 1, function(ii) paste(colnames(X)[ii],collapse="+"))
    combo.drugs <- apply(sample.pairs, 1, function(ii) paste(sort(xdrugs[ii]),collapse="+"))

    Matrix::tail(sort(table(combo.drugs)))
    sum(table(combo.drugs)>=15)
    sel.combo <- names(which(table(combo.drugs)>=15))
    jj <- which( combo.drugs %in% sel.combo)
    length(jj)
    comboX <- comboX[,jj,drop=FALSE]
    combo.drugs <- combo.drugs[jj]
    colnames(comboX) <- paste0(combo.drugs,"_combo",1:ncol(comboX))

    cat("Calculating drug-combo enrichment using GSEA ...\n")
    dim(comboX)
    res.combo <- pgx.computeDrugEnrichment(
        obj, X=comboX, xdrugs=combo.drugs, methods="GSEA", nprune=nprune)
    res.combo <- res.combo[["GSEA"]]
    names(res.combo)

    return(res.combo)
}
