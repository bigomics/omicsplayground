##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Batch correction of counts data
##
##

NORMALIZATION.METHODS <- c("no_normalization","cpm","TMM","RLE","quantile","SVA")
##nparam=NULL;niter=1;resample=1;normalization=NORMALIZATION.METHODS[1:3];show.progress=1

pgx.performNormalization <- function(zx, methods)
{
    if(length(methods)==0 && methods[1] %in% c("","no","none")) {
        return(zx)
    }

    for(mtd in methods) {
        if(mtd=="mean") {
            ## normalization on individual mean
            zx <- t(t(zx) - Matrix::colMeans(zx)) + mean(zx,na.rm=TRUE)
        } else if(mtd=="scale") {
            ## normalization on individual mean
            zx <- sd(zx)*scale(zx) + mean(zx)
        } else if(mtd=="CPM") {
            ## normalization on total counts (linear scale)
            ##nx <- t(t(nx) / colSums(nx)) * mean(nx,na.rm=TRUE)
            zx <- edgeR::cpm(2**zx, log=TRUE)
        } else if(mtd=="TMM") {
            ## normalization on total counts (linear scale)
            zx <- normalizeTMM(2**zx, log=TRUE) ## does TMM on counts
        } else if(mtd=="RLE") {
            ## normalization on total counts (linear scale)
            zx <- normalizeRLE(2**zx, log=TRUE) ## does TMM on counts
        } else if(mtd %in% c("upperquartile")) {
            ## normalization on total counts (linear scale)
            zx <- normalizeTMM(2**zx, log=TRUE, method=mtd) ## does TMM on counts
        } else if(mtd=="quantile") {
            require(preprocessCore)
            new.zx <- normalize.quantiles(as.matrix(2**zx))
            rownames(new.zx) <- rownames(zx)
            colnames(new.zx) <- colnames(zx)
            zx <- log2(new.zx + 1)
        } else if(mtd=="log-quantile") {
            require(preprocessCore)
            new.zx <- normalize.quantiles(as.matrix(zx))
            rownames(new.zx) <- rownames(zx)
            colnames(new.zx) <- colnames(zx)
            zx <- new.zx
        }
    } ## end of for method
    return(zx)
}

randomImputeMissing <- function(x) {
    i=1
    for(i in 1:ncol(x)) {
        jj <- which(is.na(x[,i]) | x[,i]=="NA")
        if(length(jj)) {
            rr <- sample( x[-jj,i], length(jj), replace=TRUE)
            x[jj,i] <- rr
        }
    }
    return(x)
}

pgx.performBatchCorrection <- function(ngs, zx, batchparams, method=c("ComBat","BMC","limma","MNN","fastMNN"))
{
    require(limma)
    require(irlba)
    require(qlcMatrix)

    ## precompute PCA
    svd <- irlba(zx - Matrix::rowMeans(zx), nv=3)
    Y <- ngs$samples[colnames(zx),]

    batchparams0 <- setdiff(batchparams, colnames(Y))
    batchparams1 <- intersect(batchparams, colnames(Y))

    ##---------------------------------------------------------------------
    ## Correct for conceptual parameters
    ##---------------------------------------------------------------------
    if(length(batchparams0)>0) {
        for(batchpar in batchparams0) {
            if(batchpar=="<PC1>") {
                sv1 <- svd$v[,1]
                zx <- limma::removeBatchEffect(zx, covariates=sv1)
            } else if(batchpar=="<PC2>") {
                sv2 <- svd$v[,2]
                zx <- limma::removeBatchEffect(zx, covariates=sv2)
            } else if(batchpar=="<PC3>") {
                sv3 <- svd$v[,3]
                zx <- limma::removeBatchEffect(zx, covariates=sv3)
            } else if(1 && batchpar=="<XY>") {
                ##svd <- svd(zx - rowMeans(zx))
                xgenes <- ngs$genes[rownames(X),]
                gx <- which(xgenes$chr %in% c("X",23))
                gy <- which(xgenes$chr %in% c("Y",24))
                xy <- NULL
                if(length(gx)) {
                    xy <- cbind(xy, Matrix::colMeans(ngs$X[gx,,drop=FALSE],na.rm=TRUE))
                }
                if(length(gy)) {
                    xy <- cbind(xy, Matrix::colMeans(ngs$X[gy,,drop=FALSE],na.rm=TRUE))
                }
                if(!is.null(xy)) {
                    zx <- limma::removeBatchEffect(zx, covariates=xy)
                }
            } else if(batchpar %in% colnames(Y)) {
                batch <- Y[,batchpar]
                class(batch)
                if(class(batch)=="numeric") {
                    ## treat as numeric
                    batch0 <- scale(as.numeric(batch))
                    batch0[is.na(batch0)] <- mean(batch0,na.rm=TRUE) ## impute NA at mean
                    zx <- limma::removeBatchEffect(zx, covariates=batch0)
                } else {
                    ## treat as factor variable
                    batch0 <- as.character(batch)
                    batch0[is.na(batch0)] <- "NA" ## NA as separate group??
                    zx <- pgx.removeBatchFactor(zx, batch0, method)
                }  ## end of iter

            } else if(batchpar=="<SVA>") {
                require(sva)
                ##group <- ngs$samples$group
                mod1 = model.matrix( ~ group, data=ngs$samples)
                mod0 = cbind(mod1[,1])
                ##mod0 = model.matrix( ~ 1, data=ngs$samples)
                sv <- sva( 0.0001+zx, mod1, mod0, n.sv=NULL)$sv
                ##sv <- svaseq( 2**zx, mod1, mod0, n.sv=NULL)$sv
                zx <- removeBatchEffect(zx, covariates=sv, design=mod1)
            } else if(batchpar=="<NNM>") {
                y <- ngs$samples$group
                table(y)
                zx <- gx.nnmcorrect( zx, y, k=5)
            } else if(batchpar=="<NNM2>") {
                y <- ngs$samples$group
                zx <- gx.nnmcorrect2( zx, y)
            } else {
                cat("warning:: unknown batch parameter\n")
            }
        }
    }

    ##---------------------------------------------------------------------
    ## Correct for parameters in phenotype: iterative or at once using
    ## model matrix.
    ## ---------------------------------------------------------------------
    ITERATIVE.CORRECT=TRUE
    ITERATIVE.CORRECT=FALSE
    if(length(batchparams1)>0) {
        if(ITERATIVE.CORRECT) {
            batchpar=batchparams1[1]
            for(batchpar in batchparams1) {
                batch <- Y[,batchpar]
                class(batch)
                if(class(batch)=="numeric") {
                    ## treat as numeric
                    batch0 <- scale(as.numeric(batch))
                    batch0[is.na(batch0)] <- mean(batch0,na.rm=TRUE) ## impute NA at mean
                    zx <- limma::removeBatchEffect(zx, covariates=batch0)
                } else {
                    ## treat as factor variable
                    batch0 <- as.character(batch)
                    batch0[is.na(batch0)] <- "NA" ## NA as separate group??
                    zx <- pgx.removeBatchFactor(zx, batch0, method)
                }
            }  ## end of iter
        } else {
            Y <- tidy.dataframe(ngs$samples[colnames(zx),])
            ny <- apply(Y,2,function(x) length(setdiff(unique(x),NA)))
            Y <- Y[,which(ny>1),drop=FALSE]
            pp <- intersect(batchparams1, colnames(Y))
            Y <- Y[,pp,drop=FALSE]
            Y <- randomImputeMissing(Y)  ## NEED RETHINK!!!
            ##batch.formula <- formula(paste("~ 0 + ",paste(pp,collapse=" + ")))
            batch.formula <- formula(paste("~ ",paste(pp,collapse=" + ")))
            B <- model.matrix(batch.formula, data=Y)
            dim(B)
            ##B1 <- expandAnnotationMatrix(Y[,pp])
            B <- B[match(colnames(zx),rownames(B)),,drop=FALSE]
            rownames(B) <- colnames(zx)
            group <- ngs$samples[colnames(zx),"group"]
            design <- model.matrix( ~ group)
            log <- capture.output({
                suppressWarnings(
                    bx <- limma::removeBatchEffect(zx, covariates=B, design=design)
                )
            })
            zx <- bx
        }

    }

    return(zx)
}


pgx.removeBatchFactor <- function(zx, batch,
                                  method=c("ComBat","BMC","limma","MNN","fastMNN"))
{
    ## treat as factor variable
    method <- method[1]
    batch0 <- as.character(batch)
    batch0[is.na(batch0)] <- "NA" ## NA as separate group??
    if(method=="MNN") {
        require(scran)
        matlist <- tapply(1:ncol(zx), batch0, function(i) zx[,i,drop=FALSE])
        ##out <- mnnCorrect( matlist[[1]], matlist[[2]])
        suppressWarnings( out <- do.call( scran::mnnCorrect,
                                         c(matlist, pc.approx=TRUE)))
        new.zx <- do.call(cbind, out$corrected)
        colnames(new.zx) <- unlist(lapply(matlist,colnames))
        zx <- new.zx[,colnames(zx)]
    } else if(method=="fastMNN") {
        require(scran)
        d = min(50,ncol(zx)/2)
        matlist <- tapply(1:ncol(zx), batch0, function(i) zx[,i,drop=FALSE])
        out <- do.call(fastMNN, c(matlist,d=d))
        cor.exp <- tcrossprod(out$rotation[,], out$corrected)
        rownames(cor.exp) <- rownames(zx)
        colnames(cor.exp) <- unlist(lapply(matlist,colnames))
        zx <- cor.exp[,colnames(zx)]
    } else if(method=="limma") {
        require(limma)
        zx <- limma::removeBatchEffect(zx, batch=batch0)
    } else if(method=="ComBat") {
        require(sva)
        zx <- ComBat(zx, batch = batch0)
    } else if(method=="BMC") {
        ## batch mean center
        matlist <- tapply(1:ncol(zx), batch0, function(i) zx[,i,drop=FALSE])
        matlist <- lapply(matlist, function(x) (x - Matrix::rowMeans(x,na.rm=TRUE)))
        new.zx <- do.call(cbind, matlist)
        new.zx <- new.zx + Matrix::rowMeans(zx,na.rm=TRUE)
        zx <- new.zx[,colnames(zx)]
    } else {
        cat("ERROR! uknown method\n")
    }
    return(zx)
}


##================================================================================
## Automagic batch correction by trying all combinations of batch parameters and
## optimizing the number of significant genes.
##================================================================================

pgx.optimizeBatchCorrection <- function(ngs, batch, contrast, nparam=NULL,
                                        normalization=NORMALIZATION.METHODS,
                                        niter=1, resample=0.9, show.progress=1)
{
    require(parallel)

    ct <- ngs$model.parameters$contr.matrix
    if(is.null(contrast))  contrast <- colnames(ct)

    makeParameterCombinations <- function( params, n) {
        if(n==0) return(c("no_correction"))
        parcomb <- sapply(params, list)
        for(k in 2:n) {
            pp <- unlist(apply(combn(params,k),2,function(x) list(x)),recursive=FALSE)
            parcomb <- c(parcomb, pp)
        }
        parcomb <- c("no_correction",parcomb)
        parcomb
        names(parcomb) <- sapply(parcomb, paste, collapse="+")
        return(parcomb)
    }

    ##params <- c("gender","age","LDH.ratio","Chemotherapy")
    if(is.null(nparam))  nparam <- length(batch)
    nparam <- min(nparam, length(batch))
    nparam

    cat(">> correcting for",length(batch),"batch parameters: ",
        paste(batch,collapse=" "),"\n")

    parcomb <- makeParameterCombinations(batch, n=nparam)
    cat(">> optimizing for",length(parcomb),"parameter combinations\n")
    ##parcomb <- sample(parcomb,10)

    ##contrast=NULL;resample=0.9;niter=1
    if(niter==1) resample=-1
    if(niter!=1) show.progress=0
    niter
    out1 <- mclapply(1:niter, function(i)
        pgx._runComputeNumSig(ngs, parcomb, contrast=contrast,
                              normalization=normalization,
                              resample=resample,
                              show.progress=show.progress)
        )
    numsig <- unlist(out1)
    numsig <- tapply( numsig, names(numsig), mean, na.rm=TRUE)
    tail(sort(numsig),200)

    ## make nice data.frame
    numsig.details <- strsplit(names(numsig),split=":")
    norm.method <- sapply(numsig.details,"[",1)
    bc.method <- sapply(numsig.details,"[",2)
    S <- data.frame(normalization=norm.method, bc.method=bc.method, num.sig=numsig)
    bc <- strsplit(as.character(S$bc.method),split="[+]")
    bc.params <- sort(unique(unlist(bc)))
    B <- t(sapply(bc, function(bb) 1*(bc.params %in% bb)))
    if(length(bc.params)==1) B <- t(B)
    colnames(B) <- bc.params
    rownames(B) <- rownames(S)
    dim(B)
    S <- cbind(S, B)
    S <- S[order(-S$num.sig),]

    ## compute best normalized/corrected combination
    best.bcm <- names(numsig)[which.max(numsig)]
    best.bcm

    res <- list( results=S, corrected=NULL )
    return(res)
}

##contr=NULL;resample=0.9
pgx._runComputeNumSig <- function(ngs, parcomb, contrast, resample=-1,
                                  normalization=NORMALIZATION.METHODS,
                                  show.progress=1)
{
    require(limma)
    k="cpm"
    numsig <- c()
    ##NORMALIZATION.METHODS <- c("no_normalization","cpm","TMM","RLE","quantile","SVA")
    for(k in normalization) {
        aX <- NULL
        if(k=="no_normalization") aX <- log(1+ngs$counts)
        if(k=="cpm") aX <- edgeR::cpm(ngs$counts, log=TRUE)
        if(k=="TMM") aX <- log2(1+normalizeTMM(ngs$counts))
        if(k=="RLE") aX <- log2(1+normalizeRLE(ngs$counts))
        if(k=="quantile") aX <- normalizeQuantiles(log2(1+ngs$counts))
        if(k=="SVA") {
            require(sva)
            mod1 = model.matrix( ~ group, data=ngs$samples)
            mod0 = cbind(mod1[,1])
            logcpm <- edgeR::cpm(ngs$counts, log=TRUE) ## perform SVA on logCPM
            log <- capture.output({
                suppressWarnings(sv <- sva(logcpm, mod1, mod0, n.sv=NULL)$sv)
                suppressWarnings(aX <- removeBatchEffect( logcpm, covariates=sv, design=mod1))
            })
            dim(aX)
        }
        if(resample>0) {
            jj <- sample(colnames(aX),ncol(aX)*resample)
            aX <- aX[,jj]
        }
        dim(aX)

        ## solve for all combinations
        pp <- parcomb[[1]]
        if(show.progress) cat(paste0("[",k,"]")) ## show progress
        for(pp in parcomb) {
            if(show.progress) cat(".") ## show progress
            Y <- ngs$samples[colnames(aX),]
            bX <- NULL
            pp
            if(pp[1]!="no_correction") {
                model.formula <- formula(paste("~ 0 + ",paste(pp,collapse=" + ")))
                bvar <- model.matrix(model.formula, data=Y)
                group <- ngs$samples[colnames(aX),"group"]
                design <- model.matrix( ~ group)
                log <- capture.output({
                    suppressWarnings(bX <- removeBatchEffect(aX, covariates=bvar, design=design))
                })
            } else {
                bX <- aX
            }
            dim(bX)
            pp1 <- paste0(k,":",paste(pp,collapse="+"))
            pp1
            numsig[pp1] <- pgx._computeNumSig(ngs, X=bX, contrast=contrast)
        }

    }
    if(show.progress) cat("\n") ## show progress
    return(numsig)
}


##X=bX;fc=0;qv=0.05
pgx._computeNumSig <- function(ngs, X, contrast=NULL, fc=0, qv=0.05) {
    require(preprocessCore)
    ##group <- ngs$samples$group
    ##design <- model.matrix(~ 0+group)
    samples <- colnames(X)
    design <- ngs$model.parameters$design[samples,]
    contr.matrix <- ngs$model.parameters$contr.matrix
    if(is.null(contrast)) contrast <- colnames(contr.matrix)
    contr.matrix <- contr.matrix[,contrast,drop=FALSE]
    res <- ngs.fitContrastsWithLIMMA(
        X, is.logcpm=TRUE, contr.matrix, design,
        quantile.normalize=FALSE, method="limma", trend=TRUE,
        conform.output=FALSE, plot=FALSE)
    names(res)
    names(res$tables)
    fc0 <- sapply(res$tables,function(x) x$logFC)
    qv0 <- sapply(res$tables,function(x) x$adj.P.Val)
    numsig <- mean(colSums( abs(fc0) >= fc & qv0 <= qv, na.rm=TRUE))
    ##numsig <- mean(colSums(qv0 < qv))
    numsig
    return(numsig)
}
