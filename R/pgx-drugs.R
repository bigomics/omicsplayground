##X=comboX;x.drugs=combo.drugs
##X=comboX;x.drugs=combo.drugs
##nprune=250;res.mono=NULL;ntop=10;nsample=20

##obj=ngs;methods="GSEA";contrast=NULL;nprune=250
pgx.computeDrugEnrichment <- function(obj, X, x.drugs, methods=c("GSEA","cor"),
                                      nprune=250, contrast=NULL )
{
    require(fgsea)
    names(obj)
    if("gx.meta" %in% names(obj)) {
        F <- sapply(obj$gx.meta$meta,function(x) x$meta.fx)
        rownames(F) <- rownames(obj$gx.meta$meta[[1]])
        rownames(F) <- toupper(sub(".*:","",rownames(F)))
        ##colnames(F) <- names(obj$gx.meta$meta)
        F <- F[order(-rowMeans(F**2)),,drop=FALSE]
        F <- F[!duplicated(rownames(F)),,drop=FALSE]
        dim(F)
    } else {
        F <- obj
    }

    if(is.null(contrast))
        contrast <- colnames(F)
    contrast <- intersect(contrast, colnames(F))
    contrast
    F <- F[,contrast,drop=FALSE]

    ## create drug meta sets
    meta.gmt <- tapply(colnames(X), x.drugs, list)
    meta.gmt <- meta.gmt[which(sapply(meta.gmt,length)>=15)]
    length(meta.gmt)

    ## first level (rank) correlation
    cat("Calculating first level correlation ...\n")
    gg <- intersect(rownames(X), rownames(F))
    length(gg)
    if(length(gg) < 20) {
        cat("WARNING::: pgx.computeDrugEnrichment : not enough common genes!!\n")
        return(NULL)
    }
    rnk1 <- apply(X[gg,,drop=FALSE],2,rank,na.last="keep")
    rnk2 <- apply(F[gg,,drop=FALSE],2,rank,na.last="keep")
    system.time(R1 <- cor(rnk1, rnk2, use="pairwise"))
    dim(R1)
    ##require(gputools)
    ##system.time(R1 <- gpuCor(rnk1, rnk2, use="pairwise")$coefficients)
    R1 <- R1 + 1e-8*matrix(rnorm(length(R1)),nrow(R1),ncol(R1))
    colnames(R1) <- colnames(F)
    rownames(R1) <- colnames(X)
    dim(R1)

    ## experiment to drug
    require(Matrix)
    ##D <- model.matrix( ~ 0 + x.drugs)
    D <- sparse.model.matrix( ~ 0 + x.drugs)
    dim(X)
    dim(D)
    colnames(D) <- sub("^x.drugs","",colnames(D))
    rownames(D) <- colnames(X)          ## not necessary..

    results <- list()
    if("cor" %in% methods) {
        cat("Calculating drug enrichment using rank correlation ...\n")
        require(qlcMatrix)
        dim(D)
        dim(R1)
        rho2 <- qlcMatrix::corSparse(D, R1)
        rownames(rho2) <- colnames(D)
        colnames(rho2) <- colnames(R1)
        rho2 <- rho2[order(-rowMeans(rho2**2)),,drop=FALSE]
        cor.pvalue <- function(x,n) pnorm(-abs(x/((1-x**2)/(n-2))**0.5))
        P <- apply(rho2, 2, cor.pvalue, n=nrow(D))
        Q <- apply(P, 2, p.adjust, method="fdr")
        results[["cor"]] <- list( X=rho2, Q=Q, P=P, stats=R1)
    }

    if("GSEA" %in% methods) {
        cat("Calculating drug enrichment using GSEA ...\n")
        res0 <- list()
        i=1
        for(i in 1:ncol(R1)) {
            suppressWarnings(res0[[i]] <- fgsea( meta.gmt, stats=R1[,i], nperm=1000))
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
        results[["GSEA"]] <- list( X=mNES, Q=mQ, P=mP, size=msize, stats=R1 )
    }

    ## this takes only the top matching drugs for each comparison to
    ## reduce the size of the matrices
    if(nprune > 0) {
        names(results)
        k=1
        for(k in 1:length(results)) {
            res <- results[[k]]
        ## reduce solution set with top-N of each comparison??
            mtop <- apply( abs(res$X), 2, function(x) head(order(-x),nprune))
            mtop <- unique(as.vector(mtop))
            length(mtop)
            top.idx <- unique(unlist(meta.gmt[mtop]))
            length(top.idx)
            results[[k]]$stats <- res$stats[top.idx,,drop=FALSE]
            results[[k]]$X <- res$X[mtop,,drop=FALSE]
            results[[k]]$P <- res$P[mtop,,drop=FALSE]
            results[[k]]$Q <- res$Q[mtop,,drop=FALSE]
            results[[k]]$size <- res$size[mtop]
        }
    }

    ## now pool by comparison???
    meta.X <- list()
    meta.Q <- list()
    meta.P <- list()
    for(i in 1:ncol(F)) {
        meta.X[[i]]<- sapply(results, function(res) res$X[,i,drop=FALSE])
        meta.Q[[i]]<- sapply(results, function(res) res$Q[,i,drop=FALSE])
        meta.P[[i]]<- sapply(results, function(res) res$P[,i,drop=FALSE])
    }

    return(results)
}



pgx.computeComboEnrichment <- function(obj, X, x.drugs,
                                       ntop=10, nsample=20, nprune=250,
                                       contrasts=NULL, res.mono=NULL )
{
    require(fgsea)
    if(0) {
        X <- readRDS(file=file.path(FILES,"l1000_es_5685drugs.rds"))
        x.drugs <- gsub("_.*$","",colnames(X))
        length(table(x.drugs))
        dim(X)
        ntop=10;nsample=20
    }

    if("gx.meta" %in% names(obj)) {
        F <- sapply(obj$gx.meta$meta,function(x) x$meta.fx)
        rownames(F) <- rownames(obj$gx.meta$meta[[1]])
        rownames(F) <- toupper(sub(".*:","",rownames(F)))
        F <- F[order(-rowMeans(F**2)),,drop=FALSE]
        F <- F[!duplicated(rownames(F)),,drop=FALSE]
        dim(F)
    }

    if(is.null(contrasts))
        contrasts <- colnames(F)
    contrasts <- intersect(contrasts, colnames(F))
    contrasts

    ## calculate average drug profile
    if(is.null(res.mono)) {
        cat("Calculating single drug enrichment using GSEA ...\n")
        er.mono <- pgx.computeDrugEnrichment(
            obj, X, x.drugs, methods="GSEA",
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
    top.mono.up  <- apply(er.mono$X, 2, function(x) head(order(-x),ntop))
    top.mono.dn  <- apply(er.mono$X, 2, function(x) head(order(x),ntop))
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
        head(cmbn)
        p1 <- sample(which(x.drugs==cmbn[1]),nsample,replace=TRUE)
        p2 <- sample(which(x.drugs==cmbn[2]),nsample,replace=TRUE)
        pp <- cbind(p1, p2)
        sample.pairs[[k]] <- pp
    }
    sample.pairs <- do.call(rbind, sample.pairs)
    dim(sample.pairs)

    ##--------------- now create combination matrix X
    comboX <- apply(sample.pairs, 1, function(ii) rowMeans(X[,ii],na.rm=TRUE))
    dim(comboX)
    ##colnames(comboX) <- apply(combo.idx, 1, function(ii) paste(colnames(X)[ii],collapse="+"))
    combo.drugs <- apply(sample.pairs, 1, function(ii) paste(sort(x.drugs[ii]),collapse="+"))

    tail(sort(table(combo.drugs)))
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
        obj, X=comboX, x.drugs=combo.drugs, methods="GSEA", nprune=nprune)
    res.combo <- res.combo[["GSEA"]]
    names(res.combo)

    return(res.combo)
}
