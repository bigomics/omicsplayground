##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## add me
if(0) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("tximport")
    biocLite("edgeR")
    biocLite("limma")
    biocLite("DESeq2")
    biocLite("ensembldb")
    biocLite("EnsDb.Hsapiens.v86")

}

ALL.GENETEST.METHODS=c("ttest","ttest.welch","voom.limma","trend.limma","notrend.limma",
                       "deseq2.wald","deseq2.lrt","edger.qlf","edger.lrt")
##ALL.GENETEST.METHODS=c("ttest","ttest.welch","voom.limma","trend.limma",
##                       "edger.qlf","edger.lrt","deseq2.wald")
methods=ALL.GENETEST.METHODS
methods=c("ttest.welch","trend.limma","deseq2.wald","edger.qlf")
methods

##-----------------------------------------------------------------------------
##-------------------- FIT ALL CONTRASTS --------------------------------------
##-----------------------------------------------------------------------------

ngs.fitContrastsWithAllMethods <- function(counts, X=NULL, samples, design, contr.matrix, genes=NULL, 
                                           prior.cpm=1, prune.samples=FALSE,
                                           conform.output=TRUE, do.filter=TRUE, cpm.scale=1e6,
                                           remove.batch=TRUE, methods=ALL.GENETEST.METHODS,
                                           correct.AveExpr=TRUE,custom=NULL, custom.name=NULL )
{
    ##--------------------------------------------------------------
    ## Run all tests on raw counts
    ##--------------------------------------------------------------
    
    
    if(0) {
        do.filter=FALSE;conform.output=TRUE;remove.batch=FALSE;prior.cpm=1;
        custom=NULL;custom.name=NULL;cpm.scale=1e6;prune.samples=FALSE
        counts=ngs$counts;X=ngs$X;samples=ngs$samples;genes=NULL;correct.AveExpr=TRUE
        design=ngs$model.parameters$design;contr.matrix=ngs$model.parameters$contr.matrix
    }

    if(methods[1]=="*") {
        methods <- ALL.GENETEST.METHODS
    }
    methods <- intersect(methods, ALL.GENETEST.METHODS)

    message("[ngs.fitContrastsWithAllMethods] calculating methods : ", methods)
    message("[ngs.fitContrastsWithAllMethods] prune.samples = ", prune.samples)
    
    ##cat("dim(X) = ",dim(X),"\n")
    
    ## If degenerate set design to NULL
    if(!is.null(design) && ncol(design)>=ncol(X) ) {
        ## "no-replicate" design!!!
        cat("WARNING: degenerate design. setting design to NULL\n")
        contr.matrix <- design %*% contr.matrix
        design <- NULL
    }

    message("[ngs.fitContrastsWithAllMethods] 2 : ")
    ##------------------------------------------------------------------
    ## define transformation methods: log2CPM for counts
    ##------------------------------------------------------------------        
    if(is.null(X)) {
        message("[ngs.fitContrastsWithAllMethods] prior CPM counts =",prior.cpm)
        message("[ngs.fitContrastsWithAllMethods] CPM scale =",cpm.scale)
        X <- log2(t(t(counts) / Matrix::colSums(counts)) * cpm.scale + prior.cpm)  ## CPM
        X <- limma::normalizeQuantiles(X)  ## in linear space
    } else {
        message("[ngs.fitContrastsWithAllMethods] using input log-expression matrix X...")
    }
    
    ##------------------------------------------------------------------    
    ## get main grouping variable for modeling
    ##------------------------------------------------------------------
    group <- NULL
    ##if(all(rownames(contr.matrix) %in% samples$group)) {
    ##    group <- samples$group
    ##}
    if(!is.null(design)) {
        group <- colnames(design)[max.col(design)]
        if(nrow(design) == ncol(design) &&
           all(rownames(design)==colnames(design))) {
            group <- NULL
        }
    }
    message("[ngs.fitContrastsWithAllMethods] [ngs.fitContrastsWithAllMethods] 3 : ")
    
    timings <- list()
    outputs = list()

    ##---------------- t-Test methods -------------------
    if("ttest" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting using Student's t-test")
        timings[["ttest"]] <- system.time(
            outputs[["ttest"]] <- ngs.fitContrastsWithTTEST(
                X, contr.matrix, design, method="equalvar",
                conform.output=conform.output)
        )
    }

    if("ttest.rank" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting using t-test on ranks")
        rX <- scale(apply(X, 2, rank, na.last="keep"))
        timings[["ttest.rank"]] <- system.time(
            outputs[["ttest.rank"]] <- ngs.fitContrastsWithTTEST(
                rX, contr.matrix, design, method="equalvar",
                conform.output=conform.output)
        )
    }

    if("ttest.welch" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting using Welch t-test")
        timings[["ttest.welch"]] <- system.time(
            outputs[["ttest.welch"]] <- ngs.fitContrastsWithTTEST(
                X, contr.matrix, design, method="welch",
                conform.output=conform.output)
        )
    }

    ##---------------- LIMMA methods -------------------
    if("trend.limma" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting using LIMMA trend")
        tt <- system.time(
            outputs[["trend.limma"]] <- ngs.fitContrastsWithLIMMA(
                X, contr.matrix, design, method="limma", trend=TRUE,
                prune.samples = prune.samples,
                conform.output=conform.output, plot=FALSE)
        )
        timings[["trend.limma"]] <- round(as.numeric(tt),digits=4)
    }
    if(TRUE && "notrend.limma" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting using LIMMA no-trend")
        timings[["notrend.limma"]] <- system.time(
            outputs[["notrend.limma"]] <- ngs.fitContrastsWithLIMMA(
                X, contr.matrix, design, method="limma", trend=FALSE,
                prune.samples = prune.samples,
                conform.output=conform.output, plot=FALSE)
        )
    }
    if("voom.limma" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting using voom/LIMMA ")
        timings[["voom.limma"]] <- system.time(
            outputs[["voom.limma"]] <- ngs.fitContrastsWithLIMMA(
                X, contr.matrix, design, method="voom",
                prune.samples = prune.samples,
                conform.output=conform.output, plot=FALSE)
        )
    }

    ##---------------- EdgeR methods -------------------
    if("edger.qlf" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting edgeR using QL F-test ")
        timings[["edger.qlf"]] <- system.time(
            outputs[["edger.qlf"]] <- ngs.fitContrastsWithEDGER(
                counts, group, contr.matrix, design, method="qlf", X=X,
                prune.samples = prune.samples,
                conform.output=conform.output, plot=FALSE)
        )
    }
    if("edger.lrt" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting edgeR using LRT")
        timings[["edger.lrt"]] <- system.time(
            outputs[["edger.lrt"]] <- ngs.fitContrastsWithEDGER(
                counts, group, contr.matrix, design, method="lrt", X=X,
                prune.samples = prune.samples,
                conform.output=conform.output, plot=FALSE)
        )
    }

    ##---------------- DESEQ2 methods -------------------
    if("deseq2.wald" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting using DESeq2 (Wald test)")
        timings[["deseq2.wald"]] <- system.time(
            outputs[["deseq2.wald"]] <- ngs.fitConstrastsWithDESEQ2(
                counts, group, contr.matrix, design, X=X, genes=genes,
                test="Wald", prune.samples = prune.samples,
                conform.output=conform.output )
        )
    }
    if("deseq2.lrt" %in% methods) {
        message("[ngs.fitContrastsWithAllMethods] fitting using DESeq2 (LRT test)")
        timings[["deseq2.lrt"]] <- system.time(
            outputs[["deseq2.lrt"]] <- ngs.fitConstrastsWithDESEQ2(
                counts, group, contr.matrix, design, X=X, genes=genes,
                test="LRT", prune.samples = prune.samples,
                conform.output=conform.output )
        )
    }

    if(!is.null(custom)) {
        message("[ngs.fitContrastsWithAllMethods] adding custom results table")
        if(is.null(custom.name)) custom.name="custom"
        if(!all( c("tables","expr") %in% names(custom)))
            stop("custom must have 'tables' and 'expr'")
        need.tests = names(outputs[[1]]$tables)
        need.tests
        if(!all(need.tests %in% names(custom$tables)))
            stop("custom must include tables: ",paste(need.tests,collapse=" "))
        need.cols = c("gene_name","AveExpr","adj.P.Val","P.Value","logFC")
        if(!all(need.cols %in% names(custom$tables[[1]])))
            stop("custom tables must include columns: ",paste(need.cols,collapse=" "))
        outputs[[custom.name]] <- custom
    }

    ##----------------------------------------------------------------------
    ## "corrections" ...
    ##----------------------------------------------------------------------
    if(correct.AveExpr) {
        message("[ngs.fitContrastsWithAllMethods] correcting AveExpr values...")
        ## Some methods like edgeR and Deseq2 compute some weird
        ## normalized expression matrix. We need to "correct" for
        ## those.
        ## 
        exp.matrix = contr.matrix
        if(!is.null(design)) exp.matrix <- (design %*% contr.matrix)
        samplesX <- lapply(apply(exp.matrix!=0,2,which),function(i) rownames(exp.matrix)[i])
        samples1 <- lapply(apply(exp.matrix>0,2,which),function(i) rownames(exp.matrix)[i])
        samples0 <- lapply(apply(exp.matrix<0,2,which),function(i) rownames(exp.matrix)[i])

        avgX <- sapply(samplesX, function(s) rowMeans(X[,s,drop=FALSE]))
        avg.1 <- sapply(samples1, function(s) rowMeans(X[,s,drop=FALSE]))
        avg.0 <- sapply(samples0, function(s) rowMeans(X[,s,drop=FALSE]))
        
        dim(avgX)
        i=j=1
        for(i in 1:length(outputs)) {
            for(j in 1:length(outputs[[i]]$tables)) {
                outputs[[i]]$tables[[j]]$AveExpr <- avgX[,j]
                outputs[[i]]$tables[[j]]$AveExpr1 <- avg.1[,j]
                outputs[[i]]$tables[[j]]$AveExpr0 <- avg.0[,j]
            }
        }

    }
    
    ##----------------------------------------------------------------------
    ## add some statistics
    ##----------------------------------------------------------------------
    message("[ngs.fitContrastsWithAllMethods] calculating statistics...")
    i=1
    for(i in 1:length(outputs)) {

        res = outputs[[i]]
        M  = sapply( res$tables, function(x) x[,"AveExpr"])  ## !!!! edgeR and Deseq2 are weird!!!
        M0  = sapply( res$tables, function(x) x[,"AveExpr0"])
        M1  = sapply( res$tables, function(x) x[,"AveExpr1"])
        Q  = sapply( res$tables, function(x) x[,"adj.P.Val"] )
        P  = sapply( res$tables, function(x) x[,"P.Value"] )
        logFC = sapply( res$tables, function(x) x[,"logFC"] )
        colnames(M) = colnames(logFC) = colnames(P) = colnames(Q) = colnames(contr.matrix)
        rownames(M) = rownames(logFC) = rownames(P) = rownames(Q) = rownames(res$tables[[1]])
        rownames(M0) = rownames(M1) = rownames(res$tables[[1]])

        ## count significant terms
        qvalues = c(1e-16,10**seq(-8,-2,2),0.05, 0.1, 0.2, 0.5,1)
        lfc=1
        sig.both = sapply(qvalues, function(q) Matrix::colSums( (Q<=q ) * (abs(logFC)>lfc), na.rm=TRUE))
        sig.up = sapply(qvalues, function(q) Matrix::colSums( (Q<=q ) * (logFC>lfc), na.rm=TRUE))
        sig.down = sapply(qvalues, function(q) Matrix::colSums( (Q<=q ) * (logFC < -lfc), na.rm=TRUE))
        sig.notsig = sapply(qvalues, function(q) Matrix::colSums( Q>q | (abs(logFC) < lfc), na.rm=TRUE))
        if(NCOL(Q)==1) {
            sig.both <- matrix(sig.both, nrow=1)
            sig.up <- matrix(sig.up, nrow=1)
            sig.down <- matrix(sig.down, nrow=1)
            sig.notsig <- matrix(sig.notsig, nrow=1)
            rownames(sig.both) = rownames(sig.up) = colnames(Q)
            rownames(sig.down) = rownames(sig.notsig) = colnames(Q)
        }
        colnames(sig.both) = colnames(sig.notsig) = qvalues
        colnames(sig.up) = colnames(sig.down) = qvalues

        res$sig.counts = list(both=sig.both, up=sig.up, down=sig.down, notsig=sig.notsig)

        ## need this? takes space!!!
        res$p.value = P
        res$q.value = Q
        res$logFC = logFC
        res$aveExpr = M
        res$aveExpr0 = M0
        res$aveExpr1 = M1
        
        outputs[[i]] <- res
    }

    ##--------------------------------------------------------------
    ## Reshape matrices by comparison
    ##--------------------------------------------------------------
    message("[ngs.fitContrastsWithAllMethods] reshape matrices...")

    ##fdr = 0.25
    tests = colnames(outputs[[1]]$p.value)
    ntest = length(tests)
    P = lapply(1:ntest, function(i) sapply( outputs, function(x) x$p.value[,i]))
    Q = lapply(1:ntest, function(i) sapply( outputs, function(x) x$q.value[,i]))
    logFC = lapply(1:ntest, function(i) sapply( outputs, function(x) x$logFC[,i]))
    aveExpr = lapply(1:ntest, function(i) sapply( outputs, function(x) x$aveExpr[,i]))
    aveExpr0 = lapply(1:ntest, function(i) sapply( outputs, function(x) x$aveExpr0[,i]))
    aveExpr1 = lapply(1:ntest, function(i) sapply( outputs, function(x) x$aveExpr1[,i]))

    sig.up = lapply(1:ntest, function(i) {
        do.call(rbind,lapply(outputs, function(x) x$sig.counts[["up"]][i,]))
    })
    sig.down = lapply(1:ntest, function(i) {
        do.call(rbind,lapply(outputs, function(x) x$sig.counts[["down"]][i,]))
    })
    sig.notsig = lapply(1:ntest, function(i) {
        do.call(rbind,lapply(outputs, function(x) x$sig.counts[["notsig"]][i,]))
    })
    sig.both = lapply(1:ntest, function(i) {
        do.call(rbind,lapply(outputs, function(x) x$sig.counts[["both"]][i,]))
    })
    names(P) = names(Q) = names(logFC) = names(aveExpr) = tests
    names(sig.up) = names(sig.down) = names(sig.both) = names(sig.notsig) = tests
    sig.counts = list(up=sig.up, down=sig.down, both=sig.both, notsig=sig.notsig)

    ##--------------------------------------------------
    ## meta analysis, aggregate p-values
    ##--------------------------------------------------
    message("[ngs.fitContrastsWithAllMethods] aggregating p-values...")
    
    all.meta <- list()
    i=1
    for(i in 1:ntest) {

        pv = P[[i]]
        qv = Q[[i]]
        fc = logFC[[i]]
        mx = aveExpr[[i]]
        mx0 = aveExpr0[[i]]
        mx1 = aveExpr1[[i]]
        
        ## avoid strange values
        fc[is.infinite(fc) | is.nan(fc)] <- NA
        pv <- pmax(pv, 1e-99)
        pv[is.na(pv)] <- 1
        qv[is.na(qv)] <- 1


        ##!!!!!!!!!!!!!!!!!!!!!!!! NEED RETHINK !!!!!!!!!!!!!!!!!!!!!!!!
        ##meta.p = apply(pv, 1, max, na.rm=TRUE ) ## maximum statistic
        ##meta.p = apply(pv, 1, function(p) metap::allmetap(p, method="sumlog")$p[[1]]) ## Fisher's method
        ##meta.p = apply(pv, 1, function(p) metap::sumlog(p)$p) 
        meta.p = apply(pv, 1, function(p) exp(mean(log(p)))) ## geometric mean
        meta.q = p.adjust(meta.p, method="BH")
        meta.fx = rowMeans(fc, na.rm=TRUE)
        ##p.meta   = apply(pv, 1, function(p) metap::allmetap(p, method="sumlog")$p[[1]])
        ## q.meta = p.adjust(p.meta, method="BH")

        meta.avg  <- rowMeans(mx, na.rm=TRUE)
        meta.avg0 <- rowMeans(mx0, na.rm=TRUE)
        meta.avg1 <- rowMeans(mx1, na.rm=TRUE)

        meta = data.frame(fx=meta.fx, p=meta.p, q=meta.q)
        ##avg <- data.frame(avg=meta.avg, avg0=meta.avg0, avg1=meta.avg1)
        avg <- data.frame(avg.0=meta.avg0, avg.1=meta.avg1)
        rownames(meta) <- rownames(logFC[[i]])
        rownames(avg)  <- rownames(logFC[[i]])
        rownames(fc) <- NULL  ## saves memory
        rownames(pv) <- NULL
        rownames(qv) <- NULL
        all.meta[[i]] = data.frame(meta=meta, avg, fc=I(fc), p=I(pv), q=I(qv))
        if(!is.null(genes)) all.meta[[i]] <- cbind(genes, all.meta[[i]])
    }
    names(all.meta) = tests
    
    timings0 <- do.call(rbind, timings)
    colnames(timings0) <- names(timings[[1]])
    colnames(timings0) <- c("user.self","sys.self","elapsed","user.child","sys.child")
    
    res = list( outputs=outputs,  meta=all.meta, sig.counts=sig.counts,
               timings = timings0, X=X )
    return(res)
}

##--------------------------------------------------------------------------------------------
##----------------------------------- FIT ALL CONTRASTS --------------------------------------
##--------------------------------------------------------------------------------------------

##dge=fish1$cooked;trend=TRUE
ngs.fitContrastsWithTTEST <- function( X, contr.matrix, design, method="welch",
                                      conform.output=0)
{
    ##if(class(dge)!="DGEList") stop("dge must be a DGEList object")
    ##install.packages("gmodels")
    ##install.packages("lsmeans")
    
    

    tables <- list()
    i=1
    exp.matrix = contr.matrix
    if(!is.null(design)) exp.matrix <- (design %*% contr.matrix)
    for(i in 1:ncol(exp.matrix)) {
        j1 = which(exp.matrix[,i] > 0)
        j0 = which(exp.matrix[,i] < 0)
        if(method=="welch") {
            suppressWarnings( rt <- matrixTests::row_t_welch( X[,j1,drop=FALSE], X[,j0,drop=FALSE] ) )
        } else {
            suppressWarnings( rt <- matrixTests::row_t_equalvar( X[,j1,drop=FALSE], X[,j0,drop=FALSE] ) )
        }
        kk <- c("mean.x","mean.y","mean.diff","stderr","df","statistic","pvalue")
        rt = rt[,kk]
        rt$qvalue = p.adjust(rt[,"pvalue"], method="BH")
        rt$mean.value = (rt[,"mean.x"] + rt[,"mean.y"])/2
        tables[[i]] <- rt
    }
    names(tables) <- colnames(contr.matrix)
    if(conform.output==TRUE) {
        i=1
        for(i in 1:length(tables)) {
            k1 = c("mean.diff","mean.value","statistic","pvalue","qvalue","mean.y","mean.x")
            k2 = c("logFC","AveExpr","statistic","P.Value","adj.P.Val","AveExpr0","AveExpr1")
            tables[[i]] = tables[[i]][,k1]
            colnames(tables[[i]]) = k2
            ##tables[[i]] = cbind(dge$genes, tables[[i]])
        }
    }
    res = list(tables=tables)
    return(res)
}


##trend=TRUE;robust=TRUE
ngs.fitContrastsWithLIMMA <- function( X, contr.matrix, design, method=c("voom","limma"),
                                      trend=TRUE, robust=TRUE, prune.samples=FALSE,
                                      conform.output=FALSE, plot=FALSE)
{
    design
    method <- method[1]
    
    if(!is.null(design)) {
        ## With no design (grouping) we perform LIMMA not on the
        ## entire contrast matrix but per contrast one-by-one.
        ##
        message("[ngs.fitContrastsWithLIMMA] fitting LIMMA contrasts using design matrix")
        exp0 <- design %*% contr.matrix
        kk <- rownames(exp0)
        if(prune.samples) {
            kk <- rownames(exp0)[which(rowSums(abs(exp0))>0)]
        }
        design1 <- design[kk,,drop=FALSE]
        X1 <- X[,kk,drop=FALSE]
        contr1 <- contr.matrix        
        if(method=="voom") {
            v <- limma::voom(2**X1, design1, plot=plot, normalize.method="none")
            vfit <- limma::lmFit(v, design1)
            trend = FALSE  ## no need
        } else  {
            vfit <- limma::lmFit(X1, design1)
        }
        vfit <- limma::contrasts.fit(vfit, contrasts=contr1)
        efit <- limma::eBayes(vfit, trend=trend, robust=robust)  ## robust YES
        if(plot==TRUE) limma::plotSA(efit)
        tables <- list()
        i=1
        exp1 = (design1 %*% contr1)
        for(i in 1:ncol(contr1)) {
            ##coef = colnames(contr.matrix)[i]
            top = limma::topTable(efit, coef=i, sort.by="none",number=Inf, adjust.method="BH")
            j1 = which( exp1[,i] > 0 )
            j2 = which( exp1[,i] < 0 )
            ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
            mean1 = rowMeans(X1[,j1,drop=FALSE], na.rm=TRUE)
            mean0 = rowMeans(X1[,j2,drop=FALSE], na.rm=TRUE)
            top = top[rownames(X1),]
            top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1)
            Matrix::head(top,10)
            tables[[i]] = top
        }
        names(tables) <- colnames(contr1)

    } else {
        message("[ngs.fitContrastsWithLIMMA] fitting LIMMA contrasts *without* design")
        tables <- list()
        exp0 = contr.matrix ## sample-wise contrasts...
        i=1
        for(i in 1:ncol(exp0)) {
            kk <- 1:nrow(exp0)
            if(prune.samples) {
                kk <- which(!is.na(exp0[,i]) & exp0[,i]!=0 )
            }            
            ct <- exp0[kk,i]
            y <- factor(c("neg","o","pos")[2+sign(ct)] )
            design1 <- model.matrix( ~ 0 + y)
            X1 <- X[,kk,drop=FALSE]
            if(method=="voom") {
                ##message("lmFit using voom")
                v <- limma::voom(2**X1, design1, plot=FALSE)
                suppressMessages( vfit <- limma::lmFit(v, design1) )
                trend=FALSE  ## no need
            } else {
                suppressMessages( vfit <- limma::lmFit(X1, design1) )
            }
            contr1 <- matrix(c(-1,0,1),nrow=3)
            rownames(contr1) <- c("yneg","yo","ypos")
            colnames(contr1) <- "pos_vs_neg"
            contr1 <- contr1[colnames(vfit),,drop=FALSE]
            vfit <- limma::contrasts.fit(vfit, contrasts=contr1)
            efit <- limma::eBayes(vfit, trend=trend, robust=robust)
            top = limma::topTable(efit, coef=1, sort.by="none",number=Inf, adjust.method="BH")
            Matrix::head(top)
            j1 = which( ct > 0 )
            j0 = which( ct < 0 )
            mean1 = rowMeans(X1[,j1,drop=FALSE], na.rm=TRUE)
            mean0 = rowMeans(X1[,j0,drop=FALSE], na.rm=TRUE)
            top = top[rownames(X1),]
            top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1 )
            Matrix::head(top,10)
            tables[[i]] = top
        }
        names(tables) <- colnames(exp0)
    }

    if(conform.output==TRUE) {
        for(i in 1:length(tables)) {
            k1 = c("logFC","AveExpr","t","P.Value","adj.P.Val","AveExpr0","AveExpr1")
            k2 = c("logFC","AveExpr","statistic","P.Value","adj.P.Val","AveExpr0","AveExpr1")
            tables[[i]] = tables[[i]][,k1]
            colnames(tables[[i]]) = k2
        }
    }
    res = list( tables=tables)
    return(res)
}

##method="qlf";robust=TRUE;plot=FALSE;conform.output=TRUE
ngs.fitContrastsWithEDGER <- function( counts, group, contr.matrix, design,
                                      method=c("qlf","lrt"), prune.samples=FALSE, X=NULL, 
                                      conform.output=FALSE, robust=TRUE, plot=TRUE)
{
    
    

    method=method[1]

    exp0 = contr.matrix
    if(!is.null(design)) exp0 = design %*% contr.matrix
    
    if(prune.samples) {
        kk <- which(rowSums(exp0!=0,na.rm=TRUE)>0)
        exp0 <- exp0[kk,,drop=FALSE]
        counts <- counts[,kk,drop=FALSE]
        if(!is.null(design)) design <- design[kk,]
        if(is.null(design))  contr.matrix <- contr.matrix[kk,,drop=FALSE]
        if(!is.null(group))  group <- group[kk]
    }
    
    dge <- edgeR::DGEList( round(counts), group=NULL)  ## we like integer counts...
    dge$samples$group <- group
    dge <- edgeR::calcNormFactors(dge, method="TMM")

    if(is.null(design) && !prune.samples) {
        message("[ngs.fitContrastsWithEDGER] fitting EDGER contrasts *without* design, no pruning ")
        res <- .ngs.fitContrastsWithEDGER.nodesign(
            dge=dge, contr.matrix=contr.matrix, method=method, 
            conform.output=conform.output, robust=robust, plot=plot)
        return(res)
    }
    if(is.null(design) && prune.samples) {
        message("[ngs.fitContrastsWithEDGER] fitting EDGER contrasts *without* design, with pruning")
        res <- .ngs.fitContrastsWithEDGER.nodesign.pruned(
            counts=counts, contr.matrix=contr.matrix, method=method, group=group,
            conform.output=conform.output, robust=robust, plot=plot)
        return(res)
    }

    message("[ngs.fitContrastsWithEDGER] fitting EDGER contrasts using design matrix")    
    dge <- edgeR::estimateDisp(dge, design=design, robust=robust)
    if(is.null(X)) X <- edgeR::cpm(counts, log=TRUE)

    ##method="qlf";robust=FALSE;plot=FALSE
    if(method=="qlf") {
        fit <- edgeR::glmQLFit(dge, design, robust=robust)
    } else if(method=="lrt") {
        fit <- edgeR::glmFit(dge, design, robust=robust)
        ## lrt <- edgeR::glmLRT(fit, design)
    } else {
        stop("unknown method")
    }

    ## get top table and calculate means
    exp.matrix <- (design %*% contr.matrix)
    tables <- list()
    i=1
    for(i in 1:ncol(contr.matrix)) {
        ##coef = colnames(contr.matrix)[i]
        cntr = contr.matrix[,i]
        if(method=="qlf") {
            ct = edgeR::glmQLFTest(fit, contrast=cntr)
        } else if(method=="lrt") {
            ct = edgeR::glmLRT(fit, contrast=cntr)
        } else {
            stop("unknown method")
        }
        ##summary(limma::decideTests(ct))
        ##if(plot) limma::plotMD(ct)
        top = edgeR::topTags(ct, n=Inf, sort.by="none")$table
        top = data.frame(top[rownames(X),])

        ## calculate means
        j1 = which( exp.matrix[,i] > 0 )
        j2 = which( exp.matrix[,i] < 0 )
        ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
        mean1 = rowMeans(X[,j1,drop=FALSE], na.rm=TRUE)
        mean0 = rowMeans(X[,j2,drop=FALSE], na.rm=TRUE)
        ## logFC of edgeR is not really reliable..
        if(conform.output) top$logFC <- (mean1 - mean0)
        
        top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1)
        tables[[i]] = top
    }
    names(tables) <- colnames(contr.matrix)

    if(conform.output==TRUE) {
        i=1
        for(i in 1:length(tables)) {
            if(method=="qlf") {
                k1 = c("logFC","logCPM","F","PValue","FDR","AveExpr0","AveExpr1")
            } else if(method=="lrt") {
                k1 = c("logFC","logCPM","LR","PValue","FDR","AveExpr0","AveExpr1")
            } else {
                stop("switch method error")
            }
            k2 = c("logFC","AveExpr","statistic","P.Value","adj.P.Val","AveExpr0","AveExpr1")
            tables[[i]] = tables[[i]][,k1]
            colnames(tables[[i]]) = k2
            ##tables[[i]] = cbind(dge$genes, tables[[i]])
        }
    }
    res = list(tables=tables)
    return(res)
}


##method="qlf";plot=TRUE;robust=FALSE;plot=TRUE;conform.output=TRUE
.ngs.fitContrastsWithEDGER.nodesign <- function( dge, contr.matrix, method=c("qlf","lrt"), X=NULL,
                                                conform.output=FALSE, robust=TRUE, plot=TRUE)
{
    ## With no design matrix, we must do EdgeR per contrast
    ## one-by-one. Warning this can become very slow.
    ##
    
    if(class(dge)!="DGEList") stop("dge must be a DGEList object")
    method=method[1]

    ##X = log2(0.0001+dge$counts)  ## assuming pseudocount already added
    if(is.null(X)) X <- edgeR::cpm(dge$counts, log=TRUE)
    dge <- edgeR::estimateDisp(dge, design=NULL, robust=robust)  ## fails...
    dge.disp <- edgeR::estimateDisp(dge$counts, design=NULL, robust=robust) 
    ##dge@.Data <- c(dge@.Data, dge.disp)
    dge$common.dispersion  <- dge.disp$common.dispersion
    dge$trended.dispersion <- dge.disp$trended.dispersion
    dge$tagwise.dispersion <- dge.disp$tagwise.dispersion
    
    contr0 <- matrix(c(-1,0,1),nrow=3)
    rownames(contr0) <- c("yneg","yo","ypos")
    colnames(contr0) <- "pos_vs_neg"

    tables <- list()
    i=1
    for(i in 1:ncol(contr.matrix)) {
        ct <- contr.matrix[,i]
        y <- factor(c("neg","o","pos")[2+sign(ct)] )
        design1 <- model.matrix( ~ 0 + y)
        ##method="qlf";robust=FALSE;plot=FALSE
        if(method=="qlf") {
            ##message("[ngs.fitContrastsWithEDGER] fitting using QL F-test")
            fit <- edgeR::glmQLFit(dge, design1, robust=robust)
            ctx <- contr0[colnames(coef(fit)),]
            res <- edgeR::glmQLFTest(fit, contrast=ctx)
        } else if(method=="lrt") {
            ##message("[ngs.fitContrastsWithEDGER] fitting using GLM/LRT")
            fit <- edgeR::glmFit(dge, design1, robust=robust)
            ctx <- contr0[colnames(coef(fit)),]
            res = edgeR::glmLRT(fit, contrast=ctx)
        } else {
            stop("unknown method: ",method)
        }
        ##summary(limma::decideTests(ct))
        top = edgeR::topTags(res, n=1e9)$table
        top = data.frame(top[rownames(X),])
        j1 = which( contr.matrix[,i] > 0 )
        j0 = which( contr.matrix[,i] < 0 )
        ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
        mean1 = rowMeans(X[,j1,drop=FALSE], na.rm=TRUE)
        mean0 = rowMeans(X[,j0,drop=FALSE], na.rm=TRUE)
        ## logFC of edgeR is not really reliable..
        if(conform.output) top$logFC <- (mean1 - mean0)

        top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1 )
        Matrix::head(top)
        tables[[i]] = top
    }
    names(tables) <- colnames(contr.matrix)

    if(conform.output==TRUE) {
        i=1
        for(i in 1:length(tables)) {
            if(method=="qlf") {
                k1 = c("logFC","logCPM","F","PValue","FDR","AveExpr0","AveExpr1")
            } else if(method=="lrt") {
                k1 = c("logFC","logCPM","LR","PValue","FDR","AveExpr0","AveExpr1")
            } else {
                stop("switch method error")
            }
            k2 = c("logFC","AveExpr","statistic","P.Value","adj.P.Val","AveExpr0","AveExpr1")
            tables[[i]] = tables[[i]][,k1]
            colnames(tables[[i]]) = k2
            ##tables[[i]] = cbind(dge$genes, tables[[i]])
        }
    }
    res = list( tables=tables)
    return(res)
}

.ngs.fitContrastsWithEDGER.nodesign.pruned <- function( counts, contr.matrix, group=NULL,
                                                       method=c("qlf","lrt"), X=NULL,
                                                       conform.output=FALSE, robust=TRUE, plot=TRUE)
{
    ## With no design matrix, we must do EdgeR per contrast
    ## one-by-one. Warning this can become very slow.
    ##


    message("[.ngs.fitContrastsWithEDGER.nodesign.pruned] called")
    message("[.ngs.fitContrastsWithEDGER.nodesign.pruned] dim(ct) = ",
            paste(dim(contr.matrix),collapse="x"))
    
    ##if(class(dge)!="DGEList") stop("dge must be a DGEList object")
    method=method[1]

    tables <- list()
    i=1
    for(i in 1:NCOL(contr.matrix)) {
        
        kk <- which(!is.na(contr.matrix[,i]) & contr.matrix[,i]!=0)
        ## X = log2(0.0001+dge$counts)  ## assuming pseudocount already added
        ## counts1 = counts[,kk,drop=FALSE]
        counts1 = counts[,kk,drop=FALSE]
        X1 <- NULL
        if(!is.null(X)) {
            X1 <- X[,kk,drop=FALSE]
        } else {
            X1 <- edgeR::cpm(counts1,log=TRUE)
        }        
        group1 <- group
        if(!is.null(group)) group1 <- group[kk]
        dge1 <- edgeR::DGEList(round(counts1), group=NULL)  ## we like integer counts...
        dge1$samples$group <- group1
        dge1 <- edgeR::calcNormFactors(dge1, method="TMM")
        ##dge1 <- edgeR::estimateDisp(dge1, design=NULL, robust=robust)  ## fails...
        dge.disp <- edgeR::estimateDisp(dge1$counts, design=NULL, robust=robust) 
        ##dge@.Data <- c(dge@.Data, dge.disp)
        dge1$common.dispersion  <- dge.disp$common.dispersion
        dge1$trended.dispersion <- dge.disp$trended.dispersion
        dge1$tagwise.dispersion <- dge.disp$tagwise.dispersion
        
        M <- matrix(c(-1,0,1),nrow=3)
        rownames(M) <- c("yneg","yo","ypos")
        colnames(M) <- "pos_vs_neg"
        
        ct <- contr.matrix[kk,i]
        y <- factor(c("neg","o","pos")[2+sign(ct)] )
        design1 <- model.matrix( ~ 0 + y)
        ##method="qlf";robust=FALSE;plot=FALSE
        if(method=="qlf") {
            ##message("[ngs.fitContrastsWithEDGER] fitting using QL F-test")
            fit <- edgeR::glmQLFit(dge1, design1, robust=robust)
            ctx <- M[colnames(coef(fit)),]
            res <- edgeR::glmQLFTest(fit, contrast=ctx)
        } else if(method=="lrt") {
            ##message("[ngs.fitContrastsWithEDGER] fitting using GLM/LRT")
            fit <- edgeR::glmFit(dge1, design1, robust=robust)
            ctx <- M[colnames(coef(fit)),]
            res <- edgeR::glmLRT(fit, contrast=ctx)
        } else {
            stop("unknown method: ",method)
        }
        ##summary(limma::decideTests(ct))
        top = edgeR::topTags(res, n=1e9)$table
        top = data.frame(top[rownames(X1),])
        contr1 <- contr.matrix[kk,i]        
        j1 = which( contr1 > 0 )
        j0 = which( contr1 < 0 )
        ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
        mean1 = rowMeans(X1[,j1,drop=FALSE], na.rm=TRUE)
        mean0 = rowMeans(X1[,j0,drop=FALSE], na.rm=TRUE)
        ## logFC of edgeR is not really reliable..
        if(conform.output) top$logFC <- (mean1 - mean0)
        top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1 )
        Matrix::head(top)
        tables[[i]] = top
    }
    names(tables) <- colnames(contr.matrix)

    if(conform.output==TRUE) {
        i=1
        for(i in 1:length(tables)) {
            if(method=="qlf") {
                k1 = c("logFC","logCPM","F","PValue","FDR","AveExpr0","AveExpr1")
            } else if(method=="lrt") {
                k1 = c("logFC","logCPM","LR","PValue","FDR","AveExpr0","AveExpr1")
            } else {
                stop("switch method error")
            }
            k2 = c("logFC","AveExpr","statistic","P.Value","adj.P.Val","AveExpr0","AveExpr1")
            tables[[i]] = tables[[i]][,k1]
            colnames(tables[[i]]) = k2
            ##tables[[i]] = cbind(dge$genes, tables[[i]])
        }
    }
    res = list( tables=tables)
    return(res)
}


ngs.fitConstrastsWithDESEQ2 <- function(counts, group, contr.matrix, design, 
                                        X=NULL, genes=NULL, test="Wald", prune.samples=FALSE,
                                        conform.output=FALSE)
{
    
    exp0 = contr.matrix
    if(!is.null(design)) exp0 = design %*% contr.matrix

    if(prune.samples) {
        kk <- which(rowSums(exp0!=0,na.rm=TRUE)>0)
        exp0 <- exp0[kk,,drop=FALSE]
        counts <- counts[,kk,drop=FALSE]
        if(!is.null(design)) design <- design[kk,]
        if(is.null(design))  contr.matrix <- contr.matrix[kk,,drop=FALSE]
        if(!is.null(group))  group <- group[kk]
    }

    if(is.null(design)) {
        message("[ngs.fitContrastsWithDESEQ2] fitting DESEQ2  *without* design")
        out <- .ngs.fitConstrastsWithDESEQ2.nodesign(
            counts=counts, contr.matrix=contr.matrix, test=test,
            prune.samples = prune.samples, conform.output=conform.output)
        return(out)
    }
    message("[ngs.fitContrastsWithDESEQ2] fitting DESEQ2 using design matrix")
    
    design.formula = formula(" ~ 0 + group")
    message("[ngs.fitContrastsWithDESEQ2] using model design: ",as.character(design.formula))

    rownames.counts <- rownames(counts)
    ##rownames(counts) <- NULL
    counts <- round(counts)  ## WARNING!!!
    if(all(rowSums(counts==0)>0)) {
        ## ERROR: 'every gene contains at least one zero, cannot compute log
        ## geometric means' so we fix it villager-style
        jmax <- which.max(rowSums(counts))
        counts[jmax,] <- pmax(counts[jmax,],1)
    }
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData=counts, design=design.formula, colData=data.frame(group))
    rownames(counts) <- rownames.counts
    ##Run DESeq : Modeling counts with generic 'group'
    fitType = 'parametric' ## sometime errors
    ##fitType = 'local'
    fitType = 'mean'
    if(test=="LRT") {
        dds <- try( DESeq2::DESeq(dds, fitType=fitType, test="LRT", reduced= ~ 1) )
        ##dds <- DESeq2::DESeq(dds, test="LRT", reduced= ~ 1)
    } else {
        dds <- try( DESeq2::DESeq(dds, fitType=fitType, test="Wald") )
        ##dds <- DESeq2::DESeq(dds, test="Wald")
    }

    ## sometime DESEQ2 fails and we resort to gene-wise estimates
    if('try-error' %in% class(dds)) {
        message("[.ngs.fitConstrastsWithDESEQ2.nodesign] retrying DESEQ2 with gene-wise estimates...")
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData=counts, design=design.formula, colData=data.frame(group))
        dds <- DESeq2::estimateSizeFactors(dds)
        dds <- DESeq2::estimateDispersionsGeneEst(dds)
        DESeq2::dispersions(dds) <- mcols(dds)$dispGeneEst
        suppressWarnings({
            if(test=="LRT") {
                dds <- try( DESeq2::nbinomLRT(dds) )
            } else {
                dds <- try( DESeq2::nbinomWaldTest(dds) )
            }
        })
    }

    ## we add the gene annotation here (not standard...)
    ##colnames(SummarizedExperiment::rowData(dds))
    if(!is.null(genes)) SummarizedExperiment::rowData(dds)$genes = genes  ## does this work??

    ## logCPM for calculating means
    ##X <- edgeR::cpm(DESeq2::counts(dds),log=TRUE,prior.count=0.000001)
    if(is.null(X)) X <- edgeR::cpm(counts, log=TRUE)
    dim(X)
    exp.matrix = contr.matrix
    if(!is.null(design)) exp.matrix = (design %*% contr.matrix)

    tables <- list()
    i=1
    for(i in 1:ncol(contr.matrix)) {
        ## manually specify contrast vector. See also https://support.bioconductor.org/p/69104/
        contr = contr.matrix[,i]
        contr[is.na(contr)] <- 0
        contr <- (contr>0)/sum(contr>0) - (contr<0)/sum(contr<0) ## mean zero, signed sum to one.
        DESeq2::resultsNames(dds)
        if(any(grepl("group", DESeq2::resultsNames(dds))) ) {
            grp.contr = contr
            names(grp.contr) = paste0("group",names(contr))
            contr = rep(0,length(DESeq2::resultsNames(dds)))
            names(contr) <- DESeq2::resultsNames(dds)
            contr[names(grp.contr)] <- grp.contr
        }
        ##resx <- DESeq2::results(dds, contrast=contr )
        ## do no set p values to NA
        resx <- DESeq2::results(dds, contrast=contr, cooksCutoff=FALSE, independentFiltering=FALSE)
        ## resx = dds.result[rownames(dds),]  ## seems not
        pos.samples = which(exp.matrix[,i] > 0)
        neg.samples = which(exp.matrix[,i] < 0)
        resx$AveExpr1 = rowMeans(X[,pos.samples,drop=FALSE], na.rm=TRUE )
        resx$AveExpr0 = rowMeans(X[,neg.samples,drop=FALSE], na.rm=TRUE )
        resx$log2BaseMean = log2(0.0001+resx$baseMean)
        if(conform.output) resx$log2FoldChange <- (resx$AveExpr1 - resx$AveExpr0)  ## recompute
        tables[[i]] <-  data.frame(resx)
        names(tables)[i] = colnames(contr.matrix)[i]
    }
    names(tables) <- colnames(contr.matrix)

    if(conform.output==TRUE) {
        for(i in 1:length(tables)) {
            k1 = c("log2FoldChange","log2BaseMean","stat","pvalue","padj","AveExpr0","AveExpr1")
            k2 = c("logFC","AveExpr","statistic","P.Value","adj.P.Val","AveExpr0","AveExpr1")
            tables[[i]] = tables[[i]][,k1]
            colnames(tables[[i]]) = k2
            ##genes = SummarizedExperiment::rowData(dds)$genes  ## BE SURE TO HAVE THIS!!
            ##tables[[i]] = cbind(genes, tables[[i]])
        }
    }
    res = list(tables=tables)
    return(res)
}

##dds=fish2$dds.object
.ngs.fitConstrastsWithDESEQ2.nodesign <- function(counts, contr.matrix, test="Wald",
                                                  prune.samples = FALSE,
                                                  conform.output=FALSE, X=NULL)
{
    
    
    ##if(class(dds)!="DESeqDataSet") stop("argument must be a DESeqDataSet object")
    ##if(!is.null(design)) stop("design must be NULL")
    counts <- round(counts)
    ##X <- edgeR::cpm(counts,log=TRUE)
    if(is.null(X)) {
        X <- edgeR::cpm(counts, log=TRUE)
    }
    if(nrow(contr.matrix) != ncol(X)) {
        stop("ngs.fitConstrastsWithDESEQ2.nodesign:: contrast matrix must be by sample")
    }

    dim(X)
    exp.matrix = contr.matrix
    
    tables <- list()
    i=1
    for(i in 1:ncol(exp.matrix)) {
        ## manual design matrix (CHECK THIS!!!)
        kk <- 1:nrow(exp.matrix)
        if(prune.samples) {
            ## prune unused samples if requested
            kk <- which(!is.na(exp.matrix[,i]) & exp.matrix[,i]!=0)
        }        
        ct <- exp.matrix[kk,i]
        y <- factor(c("neg","zero","pos")[2+sign(ct)],levels=c("neg","zero","pos"))
        ##design1 <- model.matrix( ~ factor(y[jj]==1))
        ## sample-wise model matrix (does this work???)
        ##design.formula = formula(" ~ 0 + group")
        design.formula = formula("~ 0+y")
        counts1 <- counts[,kk,drop=FALSE]
        ##rownames.counts <- rownames(counts)
        colnames(counts1) <- NULL
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData=counts1, design=design.formula, colData=data.frame(y))
        ##fitType = 'parametric'
        fitType ='mean'
        suppressWarnings({
            if(test=="LRT") {
                dds <- try(DESeq2::DESeq(dds, fitType=fitType, test="LRT", reduced=~1))
            } else {
                dds <- try(DESeq2::DESeq(dds, fitType=fitType, test="Wald")) 
            }
        })

        ## sometime DESEQ2 fails
        if('try-error' %in% class(dds)) {
            message("[.ngs.fitConstrastsWithDESEQ2.nodesign] retrying DESEQ2 with gene-wise estimates...")
            dds <- DESeq2::DESeqDataSetFromMatrix(
                countData=counts1, design=design.formula, colData=data.frame(y))
            dds <- DESeq2::estimateSizeFactors(dds)
            dds <- DESeq2::estimateDispersionsGeneEst(dds)
            DESeq2::dispersions(dds) <- mcols(dds)$dispGeneEst
            suppressWarnings({
                if(test=="LRT") {
                    dds <- try( DESeq2::nbinomLRT(dds) )
                } else {
                    dds <- try( DESeq2::nbinomWaldTest(dds) )
                }
            })
        }
                
        DESeq2::resultsNames(dds)
        ctx <- c("yneg"=-1, "yzero"=0, "ypos"=1)[DESeq2::resultsNames(dds)]
        resx <- DESeq2::results(dds, contrast=ctx, cooksCutoff=FALSE, independentFiltering=FALSE )

        ## we add the gene annotation here (not standard...)
        rownames(resx) <- rownames(SummarizedExperiment::rowData(dds))
        
        ## resx = dds.result[rownames(dds),]  ## seems not
        X1 = X[,kk,drop=FALSE]
        pos.samples = which(exp.matrix[kk,i] > 0)
        neg.samples = which(exp.matrix[kk,i] < 0)
        resx$AveExpr1 = rowMeans(X1[,pos.samples,drop=FALSE], na.rm=TRUE)
        resx$AveExpr0 = rowMeans(X1[,neg.samples,drop=FALSE], na.rm=TRUE)
        resx$log2BaseMean = log2(0.0001+resx$baseMean)
        if(conform.output) resx$log2FoldChange <- (resx$AveExpr1 - resx$AveExpr0)  ## recompute
        tables[[i]] <-  data.frame(resx)
    }
    names(tables) <- colnames(contr.matrix)

    if(conform.output==TRUE) {
        i=1
        for(i in 1:length(tables)) {
            k1 = c("log2FoldChange","log2BaseMean","stat","pvalue","padj","AveExpr0","AveExpr1")
            k2 = c("logFC","AveExpr","statistic","P.Value","adj.P.Val","AveExpr0","AveExpr1")
            tables[[i]] = tables[[i]][,k1]
            colnames(tables[[i]]) = k2
            ##genes = SummarizedExperiment::rowData(dds)$genes  ## BE SURE TO HAVE THIS!!
            ##tables[[i]] = cbind(genes, tables[[i]])
        }
    }

    res = list( tables=tables)
    return(res)
}


##--------------------------------------------------------------------------------------------
##------------------------------------ END OF FILE -------------------------------------------
##--------------------------------------------------------------------------------------------
