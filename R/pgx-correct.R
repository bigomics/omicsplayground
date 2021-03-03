##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Batch correction methods
##
##

## max.rho=0.3;batch.par="*";max.iter=10;hc.top=50;partype=NULL;pca.correct=hc.correct=1;sva.correct=mnn.correct=nnm.correct=0;bio.correct=c("mito","ribo","cc.phase","cc.score","gender")
pgx.superBatchCorrect <- function(X, pheno, model.par, partype=NULL,
                                  batch.par="*", ## batch.cov="*",
                                  lib.correct = TRUE,
                                  bio.correct=c("mito","ribo","cell_cycle","gender"),
                                  sva.correct=TRUE, pca.correct=TRUE, hc.correct=TRUE,
                                  mnn.correct=NULL, nnm.correct=FALSE,
                                  max.rho=0.3, max.iter=10, hc.top=50)
{    

    getModelMatrix <- function(v) {
        y <- as.character(pheno[,v])
        y[is.na(y)] <- "NA"  ## or impute???
        m1 <- model.matrix( ~ y)[,-1,drop=FALSE]
        colnames(m1) <- sub("^y",paste0(v,"="),colnames(m1))
        m1
    }
    if(is.null(model.par) && is.null(batch.par)) {
        stop("ERROR:: model and batch cannot be both NULL")
    }

    ## tidy up pheno matrix?? get correct parameter types
    ##pheno <- tidy.dataframe(pheno)
    pheno <- type.convert(pheno)
    
    ## setup model matrix
    mod1 <- NULL
    if(!is.null(model.par) && length(model.par)>0 ) {
        model.par <- intersect(model.par, colnames(pheno))
        mod1 <- do.call(cbind,lapply(model.par, getModelMatrix))
        rownames(mod1) <- rownames(pheno)
    }
    model.par

    ## get technical/biological effects
    Y <- pgx.computeBiologicalEffects(X)
    colnames(Y) <- paste0(".",colnames(Y))
    head(Y)

    ## add to phenotype matrix
    pheno <- cbind(pheno, Y)
    not.na <- colSums(is.na(pheno))<1
    nlev <- apply(pheno,2,function(x) length(unique(x[!is.na(x)])))
    nlev
    pheno <- pheno[,which(nlev>1 & not.na),drop=FALSE]
    partype <- sapply(pheno,class)
    partype
    
    ##--------------------------------------------------------------------
    ## select parameters
    ##--------------------------------------------------------------------    
    ##batch.par <- NULL   

    ## select all non-model variables
    if(!is.null(batch.par) && batch.par[1]=="*") {
        batch.par <- setdiff(colnames(pheno),model.par)
    }
        
    bio.correct
    if("mito" %in% bio.correct) {
        b1 <- grep("mito", colnames(pheno), value=TRUE)
        batch.par <- c(batch.par, b1)
    }
    if("ribo" %in% bio.correct) {
        b1 <- grep("ribo", colnames(pheno), value=TRUE)
        batch.par <- c(batch.par, b1)
    }
    if("cell_cycle" %in% bio.correct) {
        ##b1 <- grep("cc[.]", colnames(pheno), value=TRUE)  ## both cc.phase and cc.score
        ##b1 <- grep("cc[.].*phase$", colnames(pheno), value=TRUE) ## only cc.phase
        b1 <- grep("cc[.].*score$", colnames(pheno), value=TRUE) ## only s.score and g2m.score
        batch.par <- c(batch.par, b1)
    }
    if("gender" %in% bio.correct) {
        b1 <- grep("gender", colnames(pheno), value=TRUE)
        batch.par <- c(batch.par, b1)
    }        
    
    model.par
    batch.par <- intersect(batch.par,colnames(pheno))
    batch.par <- setdiff(batch.par, model.par)
    batch.par <- setdiff(batch.par, c("group","cluster","condition"))  ## never???    
    batch.par
    
    ##--------------------------------------------------------------------
    ## guess parameter type
    ##--------------------------------------------------------------------    

    ## select which are (continuous) covariates or (discrete) factors 
    sel1 <- which(partype %in% c("factor","character","discrete","logical"))
    batch.prm <- intersect(batch.par, names(partype[sel1])) 
    sel2 <- which(partype %in% c("integer","numeric"))    
    batch.cov <- intersect(batch.par, names(partype[sel2]))
    
    batch.prm
    batch.cov

    model.par <- intersect(model.par, colnames(pheno))
    batch.prm <- intersect(batch.prm, colnames(pheno))
    batch.cov <- intersect(batch.cov, colnames(pheno))    
    if(length(model.par)==0)  model.par <- NULL
    if(length(batch.prm)==0)  batch.prm <- NULL
    if(length(batch.cov)==0)  batch.cov <- NULL
    
    ##--------------------------------------------------------------------
    ## Check confounding
    ##--------------------------------------------------------------------
    if(!is.null(batch.prm) && !is.null(mod1)) {
        mod0 <- do.call(cbind,lapply(batch.prm, getModelMatrix))        
        rho <- cor(mod0,mod1)
        rho
        rho[is.na(rho)] <- 0
        if(max(abs(rho),na.rm=TRUE) > max.rho) {
            idx <- which(abs(rho) > max.rho, arr.ind=TRUE)
            idx
            for(i in 1:nrow(idx)) {
                v0 <- colnames(mod0)[idx[i,1]]
                v1 <- colnames(mod1)[idx[i,2]]
                cat(paste0("WARNING:: '",v0,"' is confounded with '",v1,"' ",
                           ": rho= ",round(rho[idx[i,1],idx[i,2]],3),"\n"))
            }
            confounding.pars <- colnames(mod0)[idx[,1]]
            confounding.pars <- unique(gsub("=.*","",confounding.pars))
            cat("WARNING:: removing confounding batch factors:",confounding.pars,"\n")
            batch.prm <- setdiff(batch.prm, confounding.pars)
        }
    }

    if(!is.null(batch.cov) && !is.null(mod1)) {
        cvar <- data.matrix(pheno[,batch.cov])
        rho1 <- cor(cvar,mod1,use="pairwise")
        rho1
        rho1[is.na(rho1)] <- 0        
        if(max(abs(rho1),na.rm=TRUE) > max.rho) {
            idx <- which(abs(rho1) > max.rho, arr.ind=TRUE)    
            idx
            for(i in 1:nrow(idx)) {
                v0 <- colnames(cvar)[idx[i,1]]
                v1 <- colnames(mod1)[idx[i,2]]
                cat(paste0("WARNING:: '",v0,"' is confounded with '",v1,"' ",
                           ": rho= ",round(rho1[idx[i,1],idx[i,2]],3),"\n"))
            }
            confounding.cov <- colnames(cvar)[idx[,1]]
            confounding.cov <- unique(gsub("=.*","",confounding.cov))
            confounding.cov
            cat("WARNING:: removing confounding batch covariates:",confounding.cov,"\n")
            batch.cov <- setdiff(batch.cov, confounding.cov)
            batch.cov
        }
    }

    cat("[pgx.superBatchCorrect] model.par=",model.par,"\n")
    cat("[pgx.superBatchCorrect] batch.prm=",batch.prm,"\n")
    cat("[pgx.superBatchCorrect] batch.cov=",batch.cov,"\n")
    
    cX <- X    
    mod1x <- matrix(1,ncol(cX),1)
    if(!is.null(mod1)) mod1x <- cbind(1, mod1)

    B <- mod1x[,0]  ## accumulate batch-correction matrix
    
    ##--------------------------------------------------------------------
    ## Remove (unwanted) technical experiment effects (libsize, nfeature, etc.)
    ##--------------------------------------------------------------------
    if(lib.correct) {
        sel <- grep("libsize|nfeature",colnames(pheno),value=TRUE)
        sel
        if(length(sel)) {
            cat("[pgx.superBatchCorrect] Correcting for unwanted library effects:",sel,"\n")
            exp.pheno <- pheno[,sel,drop=FALSE]
            cX <- removeBatchEffect(cX, covariates=exp.pheno, design=mod1x)
            B <- cbind(B, exp.pheno)
        }
    }

    ##--------------------------------------------------------------------
    ## Remove (unwanted) biological effects
    ##--------------------------------------------------------------------
    if(!is.null(bio.correct) && length(bio.correct)>0 && bio.correct[1]!=FALSE ) {

        p1 <- intersect(batch.prm,colnames(Y))
        cat("[pgx.superBatchCorrect] Correcting for unwanted biological factors:",p1,"\n")        
        if(length(p1)) {
            i=1
            for(i in 1:length(p1)) {
                b1 <- pheno[,p1[i]]
                cX <- removeBatchEffect(cX, batch=b1, design=mod1x)
                b1x <- model.matrix( ~b1)[,-1,drop=FALSE]
                colnames(b1x) <- sub("^b1",paste0(p1[i],"."),colnames(b1x))
                B <- cbind(B, b1x)
            }
        }

        p2 <- intersect(batch.cov,colnames(Y))
        if(length(p2)) {
            cat("[pgx.superBatchCorrect] Correcting for unwanted biological covariates:",p2,"\n")
            b2 <- pheno[,p2,drop=FALSE]
            cX <- removeBatchEffect(cX, covariates=b2, design=mod1x)
            B <- cbind(B, b2)
        }
        
        ##out <- pgx.removeBiologicalEffect(cX, pheno, model.par=model.par,
        ##                                  correct=bio.correct, force=force)
        ##cY <- out$Y  ## extended phenotypes
        ##cX <- out$X
    }
    
    ##--------------------------------------------------------------------
    ## batch correct other parameters with limma
    ##--------------------------------------------------------------------
    if(!is.null(batch.prm) && length(batch.prm)>0) {
        batch.prm
        batch.prm1 <- setdiff(batch.prm, colnames(Y))
        cat("[pgx.superBatchCorrect] Batch correction for factors:",batch.prm1,"\n")
        b <- batch.prm1[1]
        for(b in batch.prm1) {
            ##cat("Performing batch correction for factor:",b,"\n")            
            batch <- as.character(pheno[,b])
            nna <- sum(is.na(batch))
            if(nna>0) {
                ## impute missing values
                batch[is.na(batch)] <- sample(batch[!is.na(batch)],nna,replace=TRUE)
            }
            mod1x <- matrix(1,ncol(cX),1)
            if(!is.null(mod1)) mod1x <- cbind(1, mod1)
            cX <- removeBatchEffect(cX, batch=batch, design=mod1x)
            
            b1x <- model.matrix( ~batch)[,-1,drop=FALSE]
            colnames(b1x) <- sub("^batch",paste0(b,"."),colnames(b1x))
            B <- cbind(B, b1x)
        }
    }

    batch.cov
    if(!is.null(batch.cov) && length(batch.cov)>0) {        
        batch.cov
        batch.cov1 <- setdiff(batch.cov, colnames(Y))
        cat("[pgx.superBatchCorrect] Batch correction for covariates:",batch.cov1,"\n")
        for(b in batch.cov1) {
            ##cat("Performing batch correction for covariate:",b,"\n")                        
            batch <- as.numeric(pheno[,b])
            ##batch[is.na(batch)] <- "NA"
            nna <- sum(is.na(batch))
            if(nna>0) {
                batch[is.na(batch)] <- sample(batch[!is.na(batch)],nna,replace=TRUE)
            }
            mod1x <- matrix(1,ncol(cX),1)
            if(!is.null(mod1)) mod1x <- cbind(1, mod1)
            cX <- removeBatchEffect(cX, covariates=batch, design=mod1x)
            B <- cbind(B, batch)
        }
    }
    
    ##--------------------------------------------------------------------
    ## MNN correction (e.g. for single-cell)
    ##--------------------------------------------------------------------
    if(!is.null(mnn.correct)) {
        mnn.correct <- intersect(mnn.correct, colnames(pheno))
        if(length(mnn.correct)==0) mnn.correct <- NULL
    }
    if(!is.null(mnn.correct)) {
        require(batchelor)
        cat("[pgx.superBatchCorrect] Mutual Nearest Neighbour (MNN) correction on",mnn.correct,"\n")
        b <- pheno[,mnn.correct]
        out <- batchelor::mnnCorrect(cX, batch=b, cos.norm.out=FALSE)
        cX <- out@assays@data[["corrected"]]
    }

    ##--------------------------------------------------------------------
    ## Nearest-neighbour matching (NNM)
    ##--------------------------------------------------------------------
    if(nnm.correct) {
        cat("[pgx.superBatchCorrect] Correcting with nearest-neighbour matching (NNM)\n")        
        for(i in 1:length(model.par)) {
            y1 <- pheno[,model.par[i]]
            cX <- gx.nnmcorrect(cX, y1, center.x=TRUE, center.m=TRUE)$X
        }
    }

    ##--------------------------------------------------------------------
    ## SVA correction (removing unwanted variation)
    ##--------------------------------------------------------------------
    if(sva.correct && !is.null(mod1)) {
        message("[pgx.superBatchCorrect] Calculating SVA...")
        require(SmartSVA)
        ##
        ## This is a combination of methods from SVA and SmartSVA
        ## because of speed. 
        ##

        ##cX=X
        ##y <- pheno[,"dlbcl.type"]
        ##df <- data.frame(var = y)    
        ## mod1x = model.matrix( ~var, df)
        ## mod0x = model.matrix( ~1, df)        
        mod1x <- cbind(1, mod1)
        mod0x <- mod1x[,1,drop=FALSE] ## just ones...        
        if(0) {
            ## original method using SVA
            n.sv = num.sv(cX, mod1x, method="be")
            n.sv            
        } else {
            ## fast method using SmartSVA
            pp <- paste0(model.par,collapse="+")
            pp
            lm.expr <- paste0("lm(t(cX) ~ ",pp,", data=pheno)")
            X.r <- t(resid(eval(parse(text=lm.expr))))
            n.sv <- EstDimRMT(X.r, FALSE)$dim + 1
            n.sv
        }
        sv <- try( sva(cX, mod1x, mod0=mod0x, n.sv=n.sv)$sv )
        ##sv <- SmartSVA::smartsva.cpp(cX, mod1x, mod0=mod0x, n.sv=n.sv)$sv
        if(any(class(sv)=="try-error")) {
            ## try again with little bit of noise...
            a <- 0.01*mean(apply(cX,1,sd))
            cX1 <- cX + a*matrix(rnorm(length(cX)),nrow(cX),ncol(cX))
            sv <- try( sva(cX1, mod1x, mod0=mod0x, n.sv=pmax(n.sv-1,1))$sv )
        }
        if(!any(class(sv)=="try-error")) {
            message("[pgx.superBatchCorrect] Performing SVA correction...")
            ##sv <- svaseq( 2**X, mod1, mod0, n.sv=NULL)$sv
            rownames(sv) <- colnames(cX)
            colnames(sv) <- paste0("SV.",1:ncol(sv))
            cX <- removeBatchEffect(cX, covariates=sv, design=mod1x)
            ##cX <- removeBatchEffect(X, covariates=sv)
            B <- cbind(B, sv)
        }
    }
    ##gx.heatmap(cX, nmax=100, col.annot=phenox, keysize=0.9)

    ##--------------------------------------------------------------------
    ## PCA correction: remove remaining batch effect using PCA
    ## (iteratively, only SV larger than max correlated SV)
    ## --------------------------------------------------------------------
    if(pca.correct && !is.null(mod1)) {
        ii <- 1:99
        niter=0
        nremove=0
        pX <- NULL
        while(length(ii)>0 && niter<max.iter) {
            nv <- min(10,ncol(cX)-1)
            suppressWarnings(suppressMessages(
                pc <- irlba::irlba(cX, nv=nv)$v
            ))
            pc.rho <- cor(pc,mod1)
            pc.rho
            pc.rho <- apply(abs(pc.rho),1,max)
            ii <- which(pc.rho < max.rho)
            ii <- ii[ ii < which.max(pc.rho) ]
            ii
            if(length(ii)>0) {
                mod1x <- cbind(1, mod1)
                cX <- removeBatchEffect(cX, covariates=pc[,ii], design=mod1x)
                pX <- cbind( pX, pc[,ii,drop=FALSE] )
                nremove = nremove +1
            }
            niter <- niter+1
        }
        niter
        if(niter==max.iter) {
            cat("WARNING:: PCA correction did not converge after",nremove,"iterations\n")
        } else {
            cat("Performed",nremove,"iterations of PCA batch correction\n")      
        }
        if(!is.null(pX)) {
            colnames(pX) <- paste0("PC.",1:ncol(pX))
            B <- cbind(B, pX)  ## update batch correction matrix
        }
    }

    ##--------------------------------------------------------------------
    ## HC correction: remove remaining batch effect iteratively using
    ## hclust
    ## --------------------------------------------------------------------
    if(hc.correct && !is.null(mod1)) {
        ii <- 1:99
        niter=0
        nremove=0
        pX <- NULL
        while(length(ii)>0 && niter<max.iter) {
            xx <- head(cX[order(-apply(cX,1,sd)),], hc.top)
            hc <- cutree(hclust(dist(t(xx)),method="ward.D2"),2)
            table(hc)
            hc.rho <- cor(hc,mod1)
            hc.rho
            hc.rho <- apply(abs(hc.rho),1,max)
            ii <- which(hc.rho < max.rho)
            ii
            if(length(ii)>0) {
                mod1x <- cbind(1,mod1)
                hc <- scale(hc)
                cX <- removeBatchEffect(cX, covariates=hc, design=mod1x)
                pX <- cbind(pX, hc) 
                nremove = nremove + 1
            }
            niter <- niter+1
        }
        niter
        if(niter==max.iter) {
            cat("WARNING:: HC correction did not converge after",nremove,"iterations\n")
        } else {
            cat("Performed",nremove,"iterations of HC batch correction\n")
        }
        if(!is.null(pX)) B <- cbind(B, pX)  ## update batch correction matrix
    }
    
    ##--------------------------------------------------------------------
    ## important: means seems to be affected!!! regressed out??
    ##--------------------------------------------------------------------
    cX <- cX - rowMeans(cX) + rowMeans(X) 

    cat("[pgx.superBatchCorrect] almost done!\n")

    ## matrix B contains the active batch correction vectors
    ## B <- type.convert(B)    
    res <- list(X=cX, Y=pheno, B=B)   

    cat("[pgx.superBatchCorrect] done!\n")
    
    return(res)
}

##nv=3;stat="F";plot=TRUE;main=NULL
pgx.PC_correlation <- function(X, pheno, nv=3, stat="F", plot=TRUE, main=NULL) {

    getF <- function(x,y) {
        x <- t(scale(t(x)))  ## rowscale
        ii <- which(!is.na(y))
        y1 <- y[ii]
        class(y1)
        if(class(y1) %in% c("factor","character","logical")) {
            y1 <- factor(as.character(y1))
        } else {
            y1 <- y1 + 1e-8*rnorm(length(y1))
            y1 <- (y1 > median(y1))
        }
        design <- model.matrix(~ 1 + y1)
        fit <- lmFit( x[,ii], design)
        suppressWarnings( fit <- try( eBayes(fit, trend=FALSE) ) )
        class(fit)
        if(class(fit)[1]=="try-error") {
            return(NULL)
        }
        suppressMessages( top <- topTableF(fit, number=nrow(x)) )
        ##top <- topTable(fit, number=nrow(x), coef=NULL)
        return(top$F)
    }
    getCor <- function(x,y) {
        ii <- which(!is.na(y))
        y1 <- y[ii]
        if(class(y1)=="factor") y1 <- factor(as.character(y1))
        design <- model.matrix(~ 0 + y1)
        ##r1 <- cor(t(x[,ii]), design[,-1,drop=FALSE])
        r1 <- cor(t(x[,ii]), design)
        rowMeans(abs(r1))
    }

    ##nv=5
    X <- X - rowMeans(X) ## center features
    V <- irlba::irlba(X, nv=nv)$v
    rho <- list()
    ##px <- tidy.dataframe(pheno)  ## get variable types correct
    px <- pheno
    p="Chemotherapy"
    for(p in c("<random>",colnames(px))) {
        p
        if(p=="<random>") {
            y <- sample(c("a","b"), ncol(X), replace=TRUE)
        } else {
            y <- px[,p]
        }
        nlevels <- length(unique(y[!is.na(y)]))
        nlevels
        nrep <- max(table(y))
        nrep
        if(nlevels>1 && nrep >= 2) {
            if(stat=="cor") {
                rho[[p]] <- getCor(x=t(V),y)
            }
            if(stat=="F") {
                rho[[p]] <- getF(x=t(V),y)
            }
        }
    }
    
    R <- do.call(rbind, rho)
    colnames(R) <- paste0("PC",1:ncol(R))
    if(stat=="F") R <- t(t(R) / colMeans(R)) 
    R

    if(plot) {
        stat0 <- c("correlation","F-statistic")[1 + 1*(stat=="F")]
        tt0   <- c("PC correlation","PC variation")[1 + 1*(stat=="F")]
        if(is.null(main)) main <- tt0
        R <- R[,ncol(R):1]
        plt <- ggbarplot(t(R), ylab=stat0, srt=45, group.name="") +
            ## theme(
            ##     legend.key.size = unit(0.65,"lines"),
            ##     legend.key.height = unit(0.35,"lines"),
            ##     legend.text = element_text(size=9),
            ##     legend.justification = c(1,1),
            ##     legend.position = c(0.98,0.98)) +
            theme(plot.margin = ggplot2::margin(2,2,0,2,"mm"),
                  plot.title = element_text(size=12)) +
            xlab("") + ggtitle(main)
        ## plt
        return(plt)
    }
    R
}


NORMALIZATION.METHODS <- c("none","mean","scale","NC","CPM","TMM","RLE","quantile")
##nparam=NULL;niter=1;resample=1;normalization=NORMALIZATION.METHODS[1:3];show.progress=1

pgx.performNormalization.CHECK <- function(zx, methods)
{
    ## Column-wise normalization (along samples).
    ##
    ## zx:      log-expression
    ## method:   single method

    methods <- methods[1]
    
    for(mtd in methods) {
        if(mtd=="none") {
            ## normalization on individual mean
            zx <- zx
        } else if(mtd=="mean") {
            ## normalization on individual mean
            zx <- t(t(zx) - Matrix::colMeans(zx)) + mean(zx,na.rm=TRUE)
        } else if(mtd=="scale") {
            ## normalization on individual mean
            zx <- sd(zx)*scale(zx) + mean(zx)
        } else if(mtd=="CPM") {
            ## normalization on total counts (linear scale)
            ##nx <- t(t(nx) / colSums(nx)) * mean(nx,na.rm=TRUE)
            zx <- logCPM(2**zx, total=1e6)
        } else if(mtd=="qCPM") {
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
            new.zx <- normalize.quantiles(as.matrix(zx))  ## shift to avoid clipping
            rownames(new.zx) <- rownames(zx)
            colnames(new.zx) <- colnames(zx)
            zx <- new.zx
        }
    } ## end of for method
    return(zx)
}

pgx.countNormalization <- function(x, methods, keep.zero=TRUE)
{
    ## Column-wise normalization (along samples).
    ##
    ## x:        counts (linear)
    ## method:   single method

    methods <- methods[1]
    which.zero <- which(x==0, arr.ind=TRUE)
    
    for(m in methods) {
        if(m=="none") {
            ## normalization on individual mean
            x <- x
        } else if(m=="scale") {
            ## normalization on individual mean
            mx <- mean(x,na.rm=TRUE)
            x <- t(t(x) / colMeans(x,na.rm=TRUE)) * mx
        } else if(m=="CPM") {
            ##x <- edgeR::cpm(2**x, log=TRUE)
            x <- t(t(x) / colSums(x,na.rm=TRUE)) * 1e6            
        } else if(m=="TMM") {
            ## normalization on total counts (linear scale)
            x <- normalizeTMM(x, log=FALSE) ## does TMM on counts (edgeR)
        } else if(m=="RLE") {
            ## normalization on total counts (linear scale)
            x <- normalizeRLE(x, log=FALSE) ## does RLE on counts (Deseq2)
##        } else if(m %in% c("upperquartile")) {
##            ## normalization on total counts (linear scale)
##            x <- normalizeTMM(x, log=FALSE, method=m) ## does upperquartile on counts
        } else if(m=="quantile") {
            require(preprocessCore)
            new.x <- 0.01 * normalize.quantiles(as.matrix(100*x)) ## shift to avoid clipping
            rownames(new.x) <- rownames(x)
            colnames(new.x) <- colnames(x)
            x <- new.x  
        }
    } ## end of for method

    x <- pmax(x, 0)  ## prevent negative values
    ## put back zeros as zeros
    if(keep.zero && nrow(which.zero)>0) {
        x[which.zero] <- 0
    }
    
    return(x)
}

pgx.performBatchCorrection <- function(ngs, zx, batchparams,
                                       method=c("ComBat","BMC","limma","MNN","fastMNN"))
{
    require(limma)
    require(irlba)
    require(qlcMatrix)

    ## precompute PCA
    suppressWarnings(suppressMessages(
        svd <- irlba(zx - Matrix::rowMeans(zx), nv=3)
    ))
    Y <- ngs$samples[colnames(zx),]

    batchparams0 <- setdiff(batchparams, colnames(Y))
    batchparams1 <- intersect(batchparams, colnames(Y))

    ## get group from design matrix
    dd <- ngs$model.parameters$design
    group <- colnames(dd)[max.col(dd)]

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
                    zx <- pgx.removeBatchEffect(zx, batch0, method)
                }  ## end of iter

            } else if(batchpar=="<SVA>") {
                require(sva)
                ##group <- ngs$samples$group                
                ##mod1 = model.matrix( ~ group, data=ngs$samples)
                mod1 = model.matrix( ~group)                
                ## mod1 <- ngs$model.parameters$design
                mod0 <- cbind(mod1[,1])
                ##mod0 = model.matrix( ~ 1, data=ngs$samples)
                sv <- sva( 0.0001+zx, mod1, mod0, n.sv=NULL)$sv
                ##sv <- svaseq( 2**zx, mod1, mod0, n.sv=NULL)$sv
                zx <- removeBatchEffect(zx, covariates=sv, design=mod1)
            } else if(batchpar=="<NNM>") {
                ##y <- ngs$samples$group
                y <- group
                ##zx <- gx.nnmcorrect( zx, y)
                zx <- gx.nnmcorrect( zx, y, center.x=TRUE, center.m=TRUE)$X                
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
                    zx <- pgx.removeBatchEffect(zx, batch0, method)
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

pgx.removeBatchEffect <- function(X, batch, model.vars=NULL,
                                  method=c("ComBat","BMC","limma","MNN","fastMNN"))
{
    ## treat as factor variable
    method <- method[1]
    batch0 <- as.character(batch)
    batch0[is.na(batch0)] <- "NA" ## NA as separate group??
    if(method=="MNN") {
        require(scran)
        matlist <- tapply(1:ncol(X), batch0, function(i) X[,i,drop=FALSE])
        ##out <- mnnCorrect( matlist[[1]], matlist[[2]])
        suppressWarnings( out <- do.call( scran::mnnCorrect,
                                         c(matlist, pc.approx=TRUE)))
        new.X <- do.call(cbind, out$corrected)
        colnames(new.X) <- unlist(lapply(matlist,colnames))
        X <- new.X[,colnames(X)]
    } else if(method=="fastMNN") {
        require(scran)
        d = min(50,ncol(X)/2)
        matlist <- tapply(1:ncol(X), batch0, function(i) X[,i,drop=FALSE])
        out <- do.call(fastMNN, c(matlist,d=d))
        cor.exp <- tcrossprod(out$rotation[,], out$corrected)
        rownames(cor.exp) <- rownames(X)
        colnames(cor.exp) <- unlist(lapply(matlist,colnames))
        X <- cor.exp[,colnames(X)]
    } else if(method=="limma") {
        require(limma)
        X <- limma::removeBatchEffect(X, batch=batch0)
    } else if(method=="ComBat") {
        require(sva)
        X <- ComBat(X, batch = batch0)
    } else if(method=="BMC") {
        ## batch mean center
        matlist <- tapply(1:ncol(X), batch0, function(i) X[,i,drop=FALSE])
        matlist <- lapply(matlist, function(x) (x - Matrix::rowMeans(x,na.rm=TRUE)))
        new.X <- do.call(cbind, matlist)
        new.X <- new.X + Matrix::rowMeans(X,na.rm=TRUE)
        X <- new.X[,colnames(X)]
    } else {
        cat("ERROR! uknown method\n")
    }
    return(X)
}


pgx.plotMitoRibo <- function(counts, percentage=TRUE) {

    tot.counts <- colSums(counts, na.rm=TRUE)
    sel.mt <- grep("^mt-",rownames(counts),ignore.case=TRUE)
    sel.rb <- grep("^rp[ls]",rownames(counts),ignore.case=TRUE)
    mito.counts <- colSums(counts[sel.mt,,drop=FALSE], na.rm=TRUE)
    ribo.counts <- colSums(counts[sel.rb,,drop=FALSE], na.rm=TRUE)
    other.counts <- tot.counts - mito.counts - ribo.counts
    ##df <- cbind( ribo=ribo.counts, mito=mito.counts, other=other.counts )
    df <- cbind( ribo=ribo.counts, mito=mito.counts )
    if(percentage) df <- round((df / tot.counts) * 100, digits=2)
    head(df)
    barplot( t(df), beside=FALSE, las=3 )
    
}

##max.rho=0.3;force.remove=TRUE;correct.mito=TRUE;correct.ribo=TRUE;correct.cc=TRUE
pgx.computeBiologicalEffects <- function(X, is.count=FALSE)
{    
    ## estimate biological variation
    ##
    ## X:     log-expression matrix
    ##

    message("[pgx.computeBiologicalEffects] estimating biological effects...")

    ## shift zero to 1% percentile
    if(!is.count) {
        q0 <- quantile(X[X>0], probs=0.01)
        q0
        tx <- pmax(X - q0,0)
        cx <- pmax(2**tx - 1, 0)  ## counts
    } else {
        cx <- X
        tx <- log2(cx + 1)
    }
    nfeature <- Matrix::colSums(cx>0)+1
    libsize  <- Matrix::colSums(cx)

    mt.genes <- grep("^MT-",rownames(X),ignore.case=TRUE,value=TRUE)
    rb.genes <- grep("^RP[SL]",rownames(X),ignore.case=TRUE,value=TRUE)
    mito = ribo = 0
    mt.genes
    rb.genes
    if(length(mt.genes)>0) {
        mito <- Matrix::colMeans(tx[mt.genes,,drop=FALSE])
        pct.mito <- Matrix::colSums(cx[mt.genes,,drop=FALSE],na.rm=TRUE) / libsize
    }
    if(length(rb.genes)>0) {
        ii <- rb.genes[order(-apply(tx[rb.genes,],1,sd,na.rm=TRUE))]
        sel20 <- head(ii,20)
        ribo <- Matrix::colMeans(tx[rb.genes,,drop=FALSE])
        ribo20 <- Matrix::colMeans(tx[sel20,,drop=FALSE])
        pct.ribo <- Matrix::colSums(cx[rb.genes,,drop=FALSE],na.rm=TRUE) / libsize
    }
    pheno <- data.frame(
        mito = mito,
        ribo = ribo,
        ##ribo20 = ribo20,
        ##pct.mito = pct.mito,
        ##pct.ribo = pct.ribo,
        libsize = log2(libsize+1),
        ##nfeature = log2(nfeature+1),
        check.names=FALSE
    )
    
    cc.score <- try(pgx.scoreCellCycle(cx))
    head(cc.score)
    if(!any(class(cc.score)=="try-error")) {
        ##cc.score <- cc.score[,c("s_score","g2m_score","diff_score")]
        cc.score <- cc.score[,c("s_score","g2m_score")]
        colnames(cc.score) <- paste0("cc.",colnames(cc.score))
        pheno <- cbind(pheno, cc.score)
    }
    pheno$gender <- pgx.inferGender(cx)

    head(pheno)    
    return(pheno)
}

##max.rho=0.3;force.remove=TRUE;correct.mito=TRUE;correct.ribo=TRUE;correct.cc=TRUE
pgx.removeBiologicalEffect.DEPRECATED <- function(X, pheno, model.par, 
                                                  correct = c("mito","ribo","cc.score","gender"),
                                                  max.rho=0.3, force=FALSE)
{    
    ## estimate biological variation
    message("[pgx.removeBiologicalEffect] Estimating mito/ribo content...")
    q0 <- quantile(X[X>0], probs=0.01)
    q0
    tX <- pmax(X - q0,0)
    cx <- 2**tX-1  ## counts
    
    mt.genes <- grep("^MT-",rownames(X),ignore.case=TRUE,value=TRUE)
    rb.genes <- grep("^RP[SL]",rownames(X),ignore.case=TRUE,value=TRUE)
    mito = ribo = 0
    mt.genes
    rb.genes
    if(length(mt.genes)>0) {
        mito <- Matrix::colSums(tX[mt.genes,,drop=FALSE]) / Matrix::colSums(tX)*100
    }
    if(length(rb.genes)>0) {
        ribo <- Matrix::colSums(tX[rb.genes,,drop=FALSE]) / Matrix::colSums(tX)*100
    }
    nfeature <- Matrix::colSums(tX>0)
    ncounts  <- Matrix::colSums(tX)
    pheno1 <- data.frame(
        pheno[,model.par,drop=FALSE],
        "<mito>" = mito,
        "<ribo>" = ribo,
        "<lib.size>" = ncounts,
        check.names=FALSE)
    
    batch.cov <- NULL
    batch.prm <- NULL
    if("mito" %in% correct) {
        batch.cov <- c(batch.cov, "<mito>")
    }
    if("ribo" %in% correct) {        
        batch.cov <- c(batch.cov, "<ribo>")
    }
    ## if("cc.phase" %in% correct) {
    ##     cc.phase <- try(pgx.inferCellCyclePhase(cx))   ## from pgx-deconv.R        
    ##     if(!any(class(cc.phase)=="try-error")) {
    ##         batch.prm <- c(batch.prm, "<phase>")
    ##         pheno1 <- cbind(pheno1, "<phase>"=cc.phase)
    ##     }    
    ## }
    if("cc.score" %in% correct) {
        message("[pgx.removeBiologicalEffect] Inferring cell cycle...")
        cc.score <- try(pgx.scoreCellCycle(cx))
        head(cc.score)
        if(!any(class(cc.score)=="try-error")) {
            batch.cov <- c(batch.cov,"<s.score>","<g2m.score>")
            cc.score$phase <- NULL
            colnames(cc.score) <- paste0("<",colnames(cc.score),">")
            pheno1 <- cbind(pheno1, cc.score)
        }    
    }
    if("gender" %in% correct) {
        if("gender" %in% colnames(pheno)) {
            gender <- pheno$gender
        } else {
            gender <- pgx.inferGender(cx)
        }
        if(!all(is.na(gender))) {
            ## pheno1$gender <- gender
            batch.prm <- c(batch.prm, "<gender>")
            pheno1 <- cbind(pheno1, "<gender>"=gender)
        }
    }
    batch.prm
    batch.cov
    
    batch.par2 = unique(c(batch.prm, batch.cov))

    pp <- unique(c(model.par,batch.par2))
    pheno1 <- pheno1[,pp,drop=FALSE]

    message("[pgx.removeBiologicalEffect] Correcting for unwanted biological variation...")    
    out <- pgx.superBatchCorrect(
        X=X, pheno = pheno1,
        model.par = model.par,
        batch.par = batch.par2,
        max.rho=max.rho, bio.correct=NULL,
        sva.correct=FALSE, pca.correct=FALSE, hc.correct=FALSE,
        max.iter=10, hc.top=50,
        force=force)

    cX <- out$X
    
    message("[pgx.removeBiologicalEffect] done")    
    Y1 <- pheno1[,2:ncol(pheno1),drop=FALSE]
    res <- list(X=cX, Y=Y1)
    return(res)
}

pgx.svaCorrect <- function(X, pheno) {
    ## 
    ## IK: not sure about this SVA correction stuff... 
    require(sva)
    require(SmartSVA)

    if(NCOL(pheno)==1) {
        pheno <- data.frame(pheno=pheno)
    }
    X <- as.matrix(X)

    ## setup model matrix
    mod1 <- c()
    pheno1 <- pheno
    colnames(pheno1) <- paste0(colnames(pheno1),"_IS_")
    for(v in colnames(pheno1)) {        
        expr <- paste0("model.matrix(~",v,",data=pheno1)")
        m1 <- eval(parse(text=expr))[,-1,drop=FALSE]
        mod1 <- cbind(mod1, m1)
    }
    colnames(mod1) <- sub("_IS_","=",colnames(mod1))
       
    ##df <- data.frame(var=y)    
    ##mod1x = model.matrix( ~var, data=df)
    ##mod0x = model.matrix( ~1, data=df)

    mod1x <- cbind(1, mod1)
    mod0x = mod1x[,1,drop=FALSE]
    ##mod0 = NULL

    message("Estimating number of surrogate variables...")
    if(0) {
        ## original method using SVA
        n.sv = num.sv(X, mod1x, method="be")
        n.sv            
    } else {
        ## fast method using SmartSVA
        ##X.r <- t(resid(lm(t(X) ~ var, data=df)))
        pp <- paste0(colnames(pheno),collapse="+")
        pp
        lm.expr <- paste0("lm(t(X) ~ ",pp,", data=pheno)")
        X.r <- t(resid(eval(parse(text=lm.expr))))        
        n.sv <- EstDimRMT(X.r, FALSE)$dim + 1
        n.sv
    }
    
    message("Calculating SVA...")
    sv <- sva(X, mod1x, mod0x, n.sv=n.sv)$sv
    ##sv <- SmartSVA::smartsva.cpp(X, mod1x, mod0=mod0x, n.sv=n.sv)$sv
    
    message("Perform batch correction...")
    ##sv <- svaseq( 2**X, mod1, mod0, n.sv=NULL)$sv
    cX <- removeBatchEffect(X, covariates=sv, design=mod1x)
    ##cX <- removeBatchEffect(X, covariates=sv)
    
    cX
}


##================================================================================
## Automagic batch correction by trying all combinations of batch parameters and
## optimizing the number of significant genes.
##================================================================================

pgx.optimizeBatchCorrection.NOTREADY <- function(ngs, batch, contrast, nparam=NULL,
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
    ##NORMALIZATION.METHODS <- c("none","cpm","TMM","RLE","quantile","SVA")
    for(k in normalization) {
        aX <- NULL
        if(k=="nono") aX <- log(1+ngs$counts)
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

##=====================================================================================
##=========================== END OF FILE =============================================
##=====================================================================================
