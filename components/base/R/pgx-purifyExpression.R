##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

if(0) {
    
    source(file.path(RDIR,"pgx-deconv.R"))
    source(file.path(RDIR,"gx-heatmap.R"))

    load("../pgx/geiger2018b-liver2-asis-8k.pgx", verbose=1)
    load("../pgx/geiger2018b-liver2-bc-8k.pgx", verbose=1)
    load("../pgx/geiger2018b-liver2-sc-8k.pgx", verbose=1)
    load("../pgx/geiger2018b-liver2-bcsc-8k.pgx", verbose=1)
    load("../pgx/geiger2018b-liver2-qc-8k.pgx", verbose=1)
    load("../pgx/geiger2018b-liver2-qcbc-8k.pgx", verbose=1)

    ##X <- 2**ngs$X
    X <- ngs$counts
    ref="hepatocyte"
    X <- t( t(X) / Matrix::colSums(X)) * 1e6
    X <- Matrix::head(X[order(-apply(X[,j1],1,sd)),],1000)
    dim(X)

    ##jj <- which(!(ngs$samples$group %in% ref)) ## samples of interest
    j0 <- which(ngs$samples$cell.type==ref)
    j1 <- setdiff( 1:ncol(X), j0)

    normalX <- X[,j0]
    tumorX  <- X[,j1]

    if(0) {
        normalX <- X[,sample(j0,5)]
        tumorX  <- X[,sample(j1,19)]
        dim(tumorX)
        i=1
        for(i in 1:ncol(tumorX)) {
            a <- (i-1)*0.05
            rr <- runif(1:ncol(normalX))
            rr
            contamination <- rowMeans(normalX %*% diag(rr/sum(rr)))
            tumorX[,i]  <- (1-a) * tumorX[,i] + a * contamination
        }
        Matrix::colSums(tumorX)
        tumorX <- t( t(tumorX) / Matrix::colSums(tumorX)) * 1e6
    }


    method <- c("nnlm","nnmf","isopurer")
    method <- c("nnlm","nnmf")
    method <- c("nnlm")
    res <- pgx.purifyExpression( tumorX, normalX[,], method=method)
    summary(res$alpha[["nnlm"]])
    res$alpha[["nnlm"]]
    names(res$xhat)

    ## Compare all alphas
    A <- do.call( cbind, res$alpha)
    dim(A)
    pairs( A, xlim=c(0,1), ylim=c(0,1) )

    ##grp <- sampleTable$group
    logx <- log2(1+X[,])
    krt <- grep("^KRT",rownames(logx),value=TRUE)
    fabp <- grep("^FABP",rownames(logx),value=TRUE)
    fc <- rowMeans(logx[,j0]) - rowMeans(logx[,j1])
    Matrix::head(sort(fc,decreasing=TRUE),20)
    topHEPA <- names(Matrix::head(sort(-fc),200))
    markers <- unique(c(topHEPA, krt, fabp))

    xhat <- cbind( res$xhat[["nnlm"]], X[,j0])
    xhat <- t( t(xhat) / Matrix::colSums(xhat)) * 1e6
    purified <- log2(1+xhat)
    ##purified <- purified[,colnames(X)]
    dim(purified)
    prex <- logx[,colnames(purified)]


    pdf("../plots/pgx-geiger2018b-liver-topHEPA3.pdf", w=16,h=20)
    gx.heatmap(prex[markers,], col.annot=aa, nmax=-1, scale="none",
               keysize=0.5, key=FALSE, mar=c(14,8),
               main="pre-purified (top HEPA markers)")
    gx.heatmap(purified[markers,], col.annot=aa, nmax=-1, scale="none",
               keysize=0.5, key=FALSE, mar=c(14,8),
               main="post-purified (top HEPA markers)")

    topsd <- Matrix::head(order(-apply( prex,1,sd)),200)
    gx.heatmap(prex[topsd,], col.annot=aa, nmax=-1, scale="none",
               keysize=0.5, key=FALSE, mar=c(14,8),
               main="pre-purified (top all SD)")
    gx.heatmap(purified[topsd,], col.annot=aa, nmax=-1, scale="none",
               keysize=0.5, key=FALSE, mar=c(14,8),
               main="post-purified (top all SD)")

    kk <- colnames(tumorX)
    topsd <- Matrix::head(order(-apply( prex[,kk],1,sd)),200)
    gx.heatmap(prex[topsd,], col.annot=aa, nmax=-1, scale="none",
               keysize=0.5, key=FALSE, mar=c(14,8),
               main="pre-purified (top pre-pure SD)")
    gx.heatmap(purified[topsd,], col.annot=aa, nmax=-1, scale="none",
               keysize=0.5, key=FALSE, mar=c(14,8),
               main="post-purified (top pre-pure SD)")

    topsd <- Matrix::head(order(-apply( purified[,kk],1,sd)),200)
    gx.heatmap(prex[topsd,], col.annot=aa, nmax=-1, scale="none",
               keysize=0.5, key=FALSE, mar=c(14,8),
               main="pre-purified (top pure SD)")
    gx.heatmap(purified[topsd,], col.annot=aa, nmax=-1, scale="none",
               keysize=0.5, key=FALSE, mar=c(14,8),
               main="post-purified (top post-pure SD)")

    dev.off()



}

PURIFY.METHODS=c("nnlm","nnmf","isopurer","demixt","undo")
method=PURIFY.METHODS
pgx.purifyExpression <- function( tumorX, normalX, method=PURIFY.METHODS)
{
    alpha <- list()
    xhat <- list()

    method <- tolower(method)
    if("nnlm" %in% method) {
        ##----------------------------------------------------------------------
        ## NNLM (BigOmics's own method...)
        ##----------------------------------------------------------------------
        
        pen <- rep(0,3)
        ##pen <- c(1,1,1)*0.01
        res <- NNLM::nnlm( normalX, tumorX, alpha=pen)
        ##res <- NNLM::nnlm( cbind(rowMeans(tumorX), normalX), tumorX, alpha=pen)
        cf <- res$coefficients
        cf
        normal.frac <- (normalX %*% cf)
        alpha0 = (1 - Matrix::colSums(normal.frac) / Matrix::colSums(tumorX) )
        alpha0

        xhat[["nnlm"]] <- pmax(tumorX - normal.frac,0)
        dim(x.hat)
        alpha[["nnlm"]] <- alpha0
    }

    if("nnmf" %in% method) {
        ##----------------------------------------------------------------------
        ## NNMF (as in vignette)
        ##----------------------------------------------------------------------
        

        ## compute proportion of contaminant content using NNMF
        k=10
        res.nmf <- NNLM::nnmf(tumorX, k=k, init = list(W0 = normalX), check.k=FALSE)

        x.hat <- res.nmf$W[,1:k,drop=FALSE] %*% res.nmf$H[1:k,,drop=FALSE]
        nnlm.alpha <- with(res.nmf, Matrix::colSums(x.hat) / Matrix::colSums(W %*% H))
        round(nnlm.alpha, 2)

        xhat[["nnmf"]] <- x.hat
        alpha[["nnmf"]] <- nnlm.alpha
    }

    if("isopurer" %in% method) {
        ##----------------------------------------------------------------------
        ## IsoPureR (quite slow...)
        ## https://cran.r-project.org/web/packages/ISOpureR/vignettes/ISOpureRGuide.pdf
        ##----------------------------------------------------------------------

        ##install.packages('ISOpureR', repos = "http://cran.stat.sfu.ca/");
        
        ##path.to.data <- file.path(system.file(package = 'ISOpureR'), 'extdata/Beer');
        ##load(file.path(path.to.data, 'beer.normaldata.250.transcripts.RData'));
        ##load(file.path(path.to.data, 'beer.tumordata.250.transcripts.30.patients.RData'));

        ISOpureS1model <- ISOpureR::ISOpure.step1.CPE(tumorX, normalX)
        ISOpureS2model <- ISOpureR::ISOpure.step2.PPE(tumorX, normalX, ISOpureS1model)
        isopurer.alpha <- ISOpureS2model$alphapurities
        isopurer.alpha

        x.hat <- ISOpureS2model$cc_cancerprofiles
        dim(x.hat)
        alpha[["isopurer"]] <- isopurer.alpha
        xhat[["isopurer"]] <- x.hat
    }

    if("demixt" %in% method) {
        ##----------------------------------------------------------------------
        ## DeMixT (crashes often...)
        ## https://bioinformatics.mdanderson.org/main/DeMixT
        ##----------------------------------------------------------------------

        ##devtools::install_github("wwylab/DeMixT")
        
        ##data(test.data1)
        ##?DeMixT
        ##head(test.data1$y)

        res <- DeMixT::DeMixT(data.Y=tumorX, data.comp1=normalX, if.filter = FALSE)
        res$pi

        Matrix::head(res$decovExprT, 3)  ## purified tumor data
        Matrix::head(res$decovExprN1, 3)  ## normal contiminant profile
        Matrix::head(res$decovMu, 3)
        Matrix::head(res$decovSigma, 3)

        x.hat <- res$decovExprT

        demixt.alpha <- (1 - res$pi[1,])
        alpha[["demixt"]] <- demixt.alpha
        xhat[["demixt"]] <- x.hat
    }


    if("undo" %in% method) {
        ##----------------------------------------------------------------------
        ## UNDO
        ##----------------------------------------------------------------------
        ##source("http://bioconductor.org/biocLite.R")
        ##BiocManager::install("UNDO", version = "3.8")
        
        ##load tumor stroma mixing tissue samples
        ##data(PureMCF7HS27)
        ##S <- exprs(PureMCF7HS27)
        ##head(S)
        ##two_source_deconv(
        ##    X, lowper=0.4, highper=0.1, epsilon1=0.01,
        ##    epsilon2=0.01, A, S[,1], S[,2], return=0)

        res <- UNDO::two_source_deconv(
            tumorX, lowper=0.8, highper=0.1, epsilon1=0.4,
            epsilon2=0, A=NULL, S1=normalX, S2=NULL, return=1)
        str(res)
        res
        res$Estimated_Mixing_Matrix
        undo.alpha <- res$Estimated_Mixing_Matrix[,2]
        x.hat <- NULL

        alpha[["undo"]] <- undo.alpha
        xhat[["undo"]] <- NULL
    }

    return(list(alpha=alpha, xhat=xhat))
}
