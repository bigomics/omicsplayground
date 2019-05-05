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

TEST.METHODS=c("ttest","voom.limma","trend.limma","notrend.limma","deseq2",
               "deseq2.wald","deseq2.lrt","edger.qlf","edger.lrt")
methods=TEST.METHODS
methods=c("ttest","trend.limma","deseq2","edger.qlf","edger.lrt")
methods

##-----------------------------------------------------------------------------
##--------------------- FIND ALL MARKERS --------------------------------------
##-----------------------------------------------------------------------------

## ngs.findAllMarkers <- function(counts, clusters) {
##     ## make model matrix for group vs. rest
##     ##clusters <- ngs$samples$cluster
##     table(clusters)
##     contr.matrix <- makeClusterContrasts(clusters)
##     contr.matrix
##     out <- ngs.fitContrastsWithAllMethods(
##         counts, samples, genes, design, contr.matrix,
##         prior.cpm=1, conform.output=TRUE, do.filter=TRUE,
##         quantile.normalize=TRUE,
##         remove.batch=TRUE, methods=TEST.METHODS,
##         custom=NULL, custom.name=NULL )
##     return(out)
## }


##-----------------------------------------------------------------------------
##-------------------- FIT ALL CONTRASTS --------------------------------------
##-----------------------------------------------------------------------------

ngs.fitContrastsWithAllMethods <- function(counts, samples, genes, design, contr.matrix,
                                           prior.cpm=1, conform.output=TRUE, do.filter=TRUE,
                                           quantile.normalize=TRUE,
                                           remove.batch=TRUE, methods=TEST.METHODS,
                                           custom=NULL, custom.name=NULL )
{
    ##--------------------------------------------------------------
    ## Run all tests on raw counts
    ##--------------------------------------------------------------
    require(edgeR)
    require(limma)
    if(0) {
        prior.cpm=1;do.filter=FALSE;conform.output=TRUE;remove.batch=FALSE;custom=NULL
        custom.name=NULL;quantile.normalize=FALSE

        counts=ngs$count;samples=ngs$samples;genes=ngs$genes
        design=ngs$model.parameters$design;contr.matrix=ngs$model.parameters$contr.matrix
    }

    ## filter counts?
    if(do.filter) {
        cat("filtering non-expressed counts\n")
        ##keep = (rowSums(edgeR::cpm(counts)>1) >= 3)
        cpm.counts <- t(t(counts) / Matrix::colSums(counts)) * 1e6
        keep = (rowSums(cpm.counts>1) >= max(2,0.01*ncol(counts)))
        table(keep)
        counts = counts[keep,]
        genes  = genes[keep,]
        dim(counts)
    }

    ##------------------------------------------------------------------
    ## prior count regularization. Stabilized logCPM with some prior
    ## count
    ##------------------------------------------------------------------
    ## based on the above, we add a prior count
    ##PRIOR.CPM = 3  ## MAGIC NUMBER .... FIDDLE....
    if(prior.cpm > 0) {
        cat("adding prior counts at PRIOR.CPM=",prior.cpm,"\n")
        CPM.FACTORS = Matrix::colSums(counts)  / 1e6
        prior.counts = (prior.cpm * CPM.FACTORS)
        counts = t(t(counts) + prior.counts)
        counts <- round(counts) ## meh...
    }

    ## Do we quantile normalize??? EdgeR and DESeq2 have their own normalizations...
    if(quantile.normalize)  X <- normalizeQuantiles(X)

    ## reference log2CPM matrix
    X <- edgeR::cpm(counts, log=TRUE, prior.count=0.0001)  ## assumed to have priorcount...

    ## main grouping variable for modeling
    group <- NULL
    if(all(rownames(contr.matrix) %in% samples$group)) {
        group <- samples$group
    }
    
    timings <- list()


    ## ------------------- prepare objects ----------------------------
    dge=dds=NULL
    timings[["ngs.cookForEDGER"]] <- system.time(
        dge <- ngs.cookForEDGER(counts, samples, genes, prior.cpm=0,
                                normalization="none",  ## really?
                                remove.batch=remove.batch, filter=FALSE)
    )
    cat("cookForEDGER took ",timings[["ngs.cookForEDGER"]][1],"s\n")

    if("deseq2" %in% methods) {
        require(DESeq2)
        timings[["ngs.cookForDESEQ2"]] <- system.time(
            ##dds <- ngs.cookForDESEQ2(counts, samples, genes, prior.cpm=prior.cpm,
            ##                         remove.batch=remove.batch, filter=FALSE)
            dds <- ngs.cookForDESEQ2(counts, samples, genes, prior.cpm=0,
                                     remove.batch=remove.batch, filter=FALSE)
        )
        cat("CookForDEseq2 took ",timings[["ngs.cookForDESEQ2"]][1],"s\n")
    }
    
    ##---------------- t-Test methods -------------------
    outputs = list()
    if("ttest" %in% methods) {
        cat("fitting using t-test")
        timings[["ttest"]] <- system.time(
            outputs[["ttest"]] <- ngs.fitContrastsWithTTEST(
                X, contr.matrix, design, conform.output=conform.output)
        )
        cat(paste0(" (",timings[["ttest"]][1],"s)\n"))
    }

    ##---------------- LIMMA methods -------------------
    if("trend.limma" %in% methods) {
        cat("fitting using LIMMA trend")
        tt <- system.time(
            outputs[["trend.limma"]] <- ngs.fitContrastsWithLIMMA(
                counts, contr.matrix, design, method="limma", trend=TRUE,
                conform.output=conform.output, plot=FALSE)
        )
        timings[["trend.limma"]] <- round(as.numeric(tt),digits=4)
        cat(paste0(" (",timings[["trend.limma"]][1],"s)\n"))
    }
    if("notrend.limma" %in% methods) {
        cat("fitting using LIMMA no-trend")
        timings[["notrend.limma"]] <- system.time(
            outputs[["notrend.limma"]] <- ngs.fitContrastsWithLIMMA(
                counts, contr.matrix, design, method="limma", trend=FALSE,
                conform.output=conform.output, plot=FALSE)
        )
        cat(paste0(" (",timings[["notrend.limma"]][1],"s)\n"))
    }
    if("voom.limma" %in% methods) {
        cat("fitting using voom/LIMMA ")
        timings[["voom.limma"]] <- system.time(
            outputs[["voom.limma"]] <- ngs.fitContrastsWithLIMMA(
                counts, contr.matrix, design, method="voom",
                conform.output=conform.output, plot=FALSE)
        )
        cat(paste0(" (",timings[["voom.limma"]][1],"s)\n"))
    }

    ##---------------- EdgeR methods -------------------
    if("edger.qlf" %in% methods) {
        cat("fitting edgeR using QL F-test")
        timings[["edger.qlf"]] <- system.time(
            outputs[["edger.qlf"]] <- ngs.fitContrastsWithEDGER(
                counts, group, contr.matrix, design, method="qlf",
                conform.output=conform.output, plot=FALSE)
        )
        cat(paste0(" (",timings[["edger.qlf"]][1],"s)\n"))
    }
    if("edger.lrt" %in% methods) {
        cat("fitting edgeR using LRT")
        timings[["edger.lrt"]] <- system.time(
            outputs[["edger.lrt"]] <- ngs.fitContrastsWithEDGER(
                counts, group, contr.matrix, design, method="lrt",
                conform.output=conform.output, plot=FALSE)
        )
        ##cat("fitting edgeR using QL F-test (",timings[["edger.lrt"]][1],"s)\n")
        cat(paste0(" (",timings[["edger.lrt"]][1],"s)\n"))
    }

    ##---------------- DESEQ2 methods -------------------
    if(!is.null(design) && "deseq2" %in% methods) {
        cat("fitting using DESeq2 (old)")
        timings[["deseq2"]] <- system.time(
            outputs[["deseq2"]] <- ngs.fitConstrastsWithDESEQ2.fromDDS(
                dds, contr.matrix, design, ## test="Wald",
                conform.output=conform.output )
        )
        cat(paste0(" (",timings[["deseq2"]][1],"s)\n"))
    }
    if("deseq2.wald" %in% methods) {
        cat("fitting using DESeq2 (Wald test)")
        timings[["deseq2.wald"]] <- system.time(
            outputs[["deseq2.wald"]] <- ngs.fitConstrastsWithDESEQ2(
                counts, group, contr.matrix, design, test="Wald",
                conform.output=conform.output )
        )
        cat(paste0(" (",timings[["deseq2.wald"]][1],"s)\n"))
    }
    if("deseq2.lrt" %in% methods) {
        cat("fitting using DESeq2 (LRT test)")
        timings[["deseq2.lrt"]] <- system.time(
            outputs[["deseq2.lrt"]] <- ngs.fitConstrastsWithDESEQ2(
                counts, group, contr.matrix, design, test="LRT",
                conform.output=conform.output )
        )
        cat(paste0(" (",timings[["deseq2.lrt"]][1],"s)\n"))
    }

    if(!is.null(custom)) {
        cat("adding custom results table\n")
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

    ## add some statistics
    i=1
    for(i in 1:length(outputs)) {

        res = outputs[[i]]
        M = sapply( res$tables, function(x) x[,"AveExpr"])
        Q = sapply( res$tables, function(x) x[,"adj.P.Val"] )
        P = sapply( res$tables, function(x) x[,"P.Value"] )
        logFC = sapply( res$tables, function(x) x[,"logFC"] )
        colnames(M) = colnames(logFC) = colnames(P) = colnames(Q) = colnames(contr.matrix)
        rownames(M) = rownames(logFC) = rownames(P) = rownames(Q) = rownames(res$tables[[1]])

        ## count significant terms
        qvalues = c(1e-16,10**seq(-8,-2,2),0.05, 0.1, 0.2, 0.5,1)
        lfc=1
        sig.both = sapply(qvalues, function(q) colSums( (Q<=q ) * (abs(logFC)>lfc), na.rm=TRUE))
        sig.up = sapply(qvalues, function(q) colSums( (Q<=q ) * (logFC>lfc), na.rm=TRUE))
        sig.down = sapply(qvalues, function(q) colSums( (Q<=q ) * (logFC < -lfc), na.rm=TRUE))
        sig.notsig = sapply(qvalues, function(q) colSums( Q>q | (abs(logFC) < lfc), na.rm=TRUE))
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

        outputs[[i]] <- res
    }

    ##--------------------------------------------------------------
    ## Reshape matrices by comparison
    ##--------------------------------------------------------------
    ##fdr = 0.25
    tests = colnames(outputs[[1]]$p.value)
    ntest = length(tests)
    P = lapply(1:ntest, function(i) sapply( outputs, function(x) x$p.value[,i]))
    Q = lapply(1:ntest, function(i) sapply( outputs, function(x) x$q.value[,i]))
    logFC = lapply(1:ntest, function(i) sapply( outputs, function(x) x$logFC[,i]))
    aveExpr = lapply(1:ntest, function(i) sapply( outputs, function(x) x$aveExpr[,i]))
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
    require(metap)
    all.meta <- list()
    i=1
    for(i in 1:ntest) {
        pv = P[[i]]
        qv = Q[[i]]
        fc = logFC[[i]]
        ##meta.p = apply(pv, 1, max, na.rm=TRUE ) ## maximum statistic
        ##meta.q = apply(qv, 1, max, na.rm=TRUE ) ## maximum statistic
        meta.p = apply(pv, 1, function(p) metap::allmetap(p, method="sumlog")$p[[1]])
        meta.q = p.adjust(meta.p, method="fdr")
        meta.fx = rowMeans(logFC[[i]], na.rm=TRUE)
        ##p.meta   = apply(pv, 1, function(p) metap::allmetap(p, method="sumlog")$p[[1]])
        ## q.meta = p.adjust(p.meta, method="fdr")
        meta = data.frame(fx=meta.fx, p=meta.p, q=meta.q)
        all.meta[[i]] = data.frame( genes, meta=meta, fc=I(fc), p=I(pv), q=I(qv))
    }
    names(all.meta) = tests

    timings0 <- do.call(rbind, timings)
    res = list( outputs=outputs,  meta=all.meta, sig.counts=sig.counts,
               timings = timings0, X=X )
    return(res)
}

##--------------------------------------------------------------------------------------------
##----------------------------------- FIT ALL CONTRASTS --------------------------------------
##--------------------------------------------------------------------------------------------

##dge=fish1$cooked;trend=TRUE
ngs.fitContrastsWithTTEST <- function( X, contr.matrix, design,
                                      conform.output=0)
{
    ##if(class(dge)!="DGEList") stop("dge must be a DGEList object")
    ##install.packages("gmodels")
    ##install.packages("lsmeans")
    require(gmodels)
    require(matrixTests)

    tables <- list()
    i=1
    exp.matrix = contr.matrix
    if(!is.null(design)) exp.matrix <- (design %*% contr.matrix)
    for(i in 1:ncol(exp.matrix)) {
        j1 = which(exp.matrix[,i] > 0)
        j0 = which(exp.matrix[,i] < 0)
        suppressWarnings( rt <- row_t_welch( X[,j1,drop=FALSE], X[,j0,drop=FALSE] ) )
        rt = rt[,c(4:6,9:12)]
        rt$qvalue = p.adjust(rt[,"pvalue"])
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
ngs.fitContrastsWithLIMMA <- function( counts, contr.matrix, design,
                                      method=c("voom","limma"),
                                      trend=TRUE, prior.count=1,
                                      conform.output=FALSE, plot=TRUE)
{
    X <- edgeR::cpm(counts, log=TRUE, prior.count=prior.count)
    design    
    if(!is.null(design)) {
        cat("fitting contrasts using LIMMA with design matrix....\n")
        if(method=="voom") {
            ##cat("lmFit using voom\n")
            v <- voom(counts, design, plot=plot, normalize.method="none")
            vfit <- lmFit(v, design)
            trend=FALSE  ## no need
        } else  {
            vfit <- lmFit(X, design)
        }
        vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
        efit <- eBayes(vfit, trend=trend, robust=TRUE)
        ##efit <- eBayes(vfit, trend=trend, robust=FALSE)
        if(plot==TRUE) plotSA(efit)

        tables <- list()
        i=1
        exp.matrix = (design %*% contr.matrix)
        for(i in 1:ncol(contr.matrix)) {
            ##coef = colnames(contr.matrix)[i]
            top = topTable(efit, coef=i, sort.by="none",number=Inf, adjust.method="BH")
            j1 = which( exp.matrix[,i] > 0 )
            j2 = which( exp.matrix[,i] < 0 )
            ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
            mean1 = rowMeans(X[,j1,drop=FALSE], na.rm=TRUE)
            mean0 = rowMeans(X[,j2,drop=FALSE], na.rm=TRUE)
            top = top[rownames(X),]
            top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1)
            head(top,10)
            tables[[i]] = top
        }
        names(tables) <- colnames(contr.matrix)

    } else {
        cat("fitting contrasts using LIMMA without design (explicit contrasts)....\n")
        tables <- list()
        i=1
        for(i in 1:ncol(contr.matrix)) {            
            ct <- contr.matrix[,i]
            y <- factor(c("neg","o","pos")[2+sign(ct)],levels=c("neg","o","pos"))
            y <- factor(c("neg","o","pos")[2+sign(ct)] )
            design0 <- model.matrix( ~ 0+y)            
            if(method=="voom") {
                ##cat("lmFit using voom\n")
                v <- voom(counts[,], design0, plot=FALSE)
                suppressMessages( vfit <- lmFit(v, design0) )
                trend=FALSE  ## no need
            } else {
                suppressMessages( vfit <- lmFit(X[,], design0) )
            }
            contr.matrix0 <- matrix(c(-1,0,1),nrow=3)
            rownames(contr.matrix0) <- c("yneg","yo","ypos")
            colnames(contr.matrix0) <- "pos_vs_neg"
            contr.matrix0 <- contr.matrix0[colnames(vfit),,drop=FALSE]
            vfit <- contrasts.fit(vfit, contrasts=contr.matrix0)
            efit <- eBayes(vfit, trend=trend, robust=TRUE)
            top = topTable(efit, coef=1, sort.by="none",number=Inf, adjust.method="BH")
            head(top)
            j1 = which( contr.matrix[,i] > 0 )
            j0 = which( contr.matrix[,i] < 0 )
            mean1 = rowMeans(X[,j1,drop=FALSE], na.rm=TRUE)
            mean0 = rowMeans(X[,j0,drop=FALSE], na.rm=TRUE)
            top = top[rownames(X),]
            top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1 )
            head(top,10)
            tables[[i]] = top
        }
        names(tables) <- colnames(contr.matrix)
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


##method="qlf";robust=FALSE;plot=FALSE;conform.output=TRUE
ngs.fitContrastsWithEDGER <- function( counts, group, contr.matrix, design,
                                      method=c("qlf","lrt"), 
                                      conform.output=FALSE, robust=FALSE, plot=TRUE)
{
    require(limma)
    require(edgeR)
    dge <- DGEList( round(counts), group=NULL)  ## we like integer counts...
    dge$samples$group <- group
    ##dge$samples = cbind(dge$samples, samples)    
    dge <- calcNormFactors(dge, method="TMM")
    if(is.null(design)) {
        res <- .ngs.fitContrastsWithEDGER.nodesign(
            dge=dge, contr.matrix=contr.matrix, method=method,
            conform.output=conform.output, robust=robust, plot=plot)
        return(res)
    }
    dge <- estimateDisp(dge, design=design, robust=robust)
    
    method=method[1]
    ##X = log2(0.0001+dge$counts)  ## hmmm....
    X = edgeR::cpm(counts, log=TRUE) 

    ##method="qlf";robust=FALSE;plot=FALSE
    if(method=="qlf") {
        ##cat("fitting using QL F-test\n")
        fit <- glmQLFit(dge, design, robust=robust)
        ##if(plot) plotQLDisp(fit)
    } else if(method=="lrt") {
        ##cat("fitting using GLM/LRT\n")
        fit <- glmFit(dge, design, robust=robust)
        lrt <- glmLRT(fit, design)
    } else {
        stop("unknown method")
    }

    ## get top table and calculate means
    exp.matrix <- design %*% contr.matrix
    tables <- list()
    i=1
    for(i in 1:ncol(contr.matrix)) {
        ##coef = colnames(contr.matrix)[i]
        cntr = contr.matrix[,i]
        if(method=="qlf") {
            ct = glmQLFTest(fit, contrast=cntr)
        } else if(method=="lrt") {
            ct = glmLRT(fit, contrast=cntr)
        } else {
            stop("unknown method")
        }
        ##summary(decideTests(ct))
        ##if(plot) plotMD(ct)
        top = topTags(ct, n=1e9)$table
        top = data.frame(top[rownames(X),])
        j1 = which( exp.matrix[,i] > 0 )
        j2 = which( exp.matrix[,i] < 0 )
        ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
        mean1 = rowMeans(X[,j1,drop=FALSE], na.rm=TRUE)
        mean0 = rowMeans(X[,j2,drop=FALSE], na.rm=TRUE)
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
.ngs.fitContrastsWithEDGER.nodesign <- function( dge, contr.matrix, method=c("qlf","lrt"),
                                                conform.output=FALSE, robust=FALSE, plot=TRUE)
{
    if(class(dge)!="DGEList") stop("dge must be a DGEList object")
    method=method[1]
    ##X = log2(0.0001+dge$counts)  ## assuming pseudocount already added
    X = edgeR::cpm(counts, log=TRUE)   
    dge <- estimateDisp(dge, design=NULL, robust=robust)
    
    contr.matrix0 <- matrix(c(-1,0,1),nrow=3)
    rownames(contr.matrix0) <- c("yneg","yo","ypos")
    colnames(contr.matrix0) <- "pos_vs_neg"
    
    tables <- list()
    i=1
    for(i in 1:ncol(contr.matrix)) {
        ct <- contr.matrix[,i]
        y <- factor(c("neg","o","pos")[2+sign(ct)] )
        design1 <- model.matrix( ~ 0+y)            
        ##method="qlf";robust=FALSE;plot=FALSE        
        if(method=="qlf") {
            ##cat("fitting using QL F-test\n")
            fit <- glmQLFit(dge[,], design1, robust=robust)
            ctx <- contr.matrix0[colnames(coef(fit)),]
            res <- glmQLFTest(fit, contrast=ctx)
        } else if(method=="lrt") {
            ##cat("fitting using GLM/LRT\n")
            fit <- glmFit(dge[,], design1, robust=robust)
            ctx <- contr.matrix0[colnames(coef(fit)),]
            res = glmLRT(fit, contrast=ctx)
        } else {
            stop("unknown method: ",method)
        }
        ##summary(decideTests(ct))
        ##if(plot) plotMD(ct)        
        top = topTags(res, n=1e9)$table
        top = data.frame(top[rownames(X),])
        j1 = which( contr.matrix[,i] > 0 )
        j0 = which( contr.matrix[,i] < 0 )
        ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
        mean1 = rowMeans(X[,j1,drop=FALSE], na.rm=TRUE)
        mean0 = rowMeans(X[,j0,drop=FALSE], na.rm=TRUE)
        top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1 )
        head(top)
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

ngs.fitConstrastsWithDESEQ2 <- function(counts, group,
                                        contr.matrix, design, test="Wald",
                                        conform.output=FALSE)
{
    if(is.null(design)) {
        out <- ngs.fitConstrastsWithDESEQ2.nodesign(
            counts=counts, contr.matrix=contr.matrix, test=test,
            conform.output=conform.output)
        return(out)
    }
    design.formula = formula(" ~ 0 + group")
    cat("using model design: ",as.character(design.formula),"\n")
    rownames.counts <- rownames(counts)
    rownames(counts) <- NULL
    dds <- DESeqDataSetFromMatrix(
        countData=counts, design=design.formula, colData=data.frame(group))
    rownames(counts) <- rownames.counts
    ##Run DESeq : Modeling counts with generic 'group'
    fitType = 'parametric'
    fitType ='local'
    if(test=="LRT") {
        dds <- DESeq(dds, fitType=fitType, test="LRT", reduced= ~ 1)
    } else {
        dds <- DESeq(dds, fitType=fitType, test="Wald")
    }
    ## we add the gene annotation here (not standard...)
    colnames(rowData(dds))
    rowData(dds)$genes = genes

    out <- ngs.fitConstrastsWithDESEQ2.fromDDS(
        dds, contr.matrix=contr.matrix, design=design,
        conform.output=conform.output)
    return(out)
}

##dds=fish2$dds.object
ngs.fitConstrastsWithDESEQ2.fromDDS <- function(dds, contr.matrix, design,
                                                conform.output=FALSE)
{
    if(class(dds)!="DESeqDataSet") stop("argument must be a DESeqDataSet object")

    dim(design)
    if(is.null(design)) {
        stop("FATAL:: ngs.fitConstrastsWithDESEQ2: not yet implemented for design==NULL!")
    }
    X <- edgeR::cpm(counts(dds),log=TRUE)
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
        resultsNames(dds)        
        if(any(grepl("group",resultsNames(dds))) ) {
            grp.contr = contr
            names(grp.contr) = paste0("group",names(contr))
            contr = rep(0,length(resultsNames(dds)))
            names(contr) <- resultsNames(dds)
            contr[names(grp.contr)] <- grp.contr
        }
        ##resx <- results(dds, contrast=contr )
        ## do no set p values to NA
        resx <- results(dds, contrast=contr, cooksCutoff=FALSE, independentFiltering=FALSE )

        ## resx = dds.result[rownames(dds),]  ## seems not
        pos.samples = which(exp.matrix[,i] > 0)
        neg.samples = which(exp.matrix[,i] < 0)
        resx$AveExpr1 = rowMeans( X[,pos.samples,drop=FALSE], na.rm=TRUE )
        resx$AveExpr0 = rowMeans( X[,neg.samples,drop=FALSE], na.rm=TRUE )
        resx$log2BaseMean = log2(0.0001+resx$baseMean)
        head(resx)

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
            ##genes = rowData(dds)$genes  ## BE SURE TO HAVE THIS!!
            ##tables[[i]] = cbind(genes, tables[[i]])
        }
    }
    res = list( tables=tables)
    return(res)
}

##dds=fish2$dds.object
ngs.fitConstrastsWithDESEQ2.nodesign <- function(counts, contr.matrix, test="Wald",
                                                 conform.output=FALSE)
{
    ##if(class(dds)!="DESeqDataSet") stop("argument must be a DESeqDataSet object")
    ##if(!is.null(design)) stop("design must be NULL")
    X <- edgeR::cpm(counts,log=TRUE)
    dim(X)
    if(nrow(contr.matrix) != ncol(X)) {
        stop("ngs.fitConstrastsWithDESEQ2.nodesign:: contrast matrix must be by sample")
    }

    dim(X)
    exp.matrix = contr.matrix
    counts1 <- counts
    ##rownames.counts <- rownames(counts)
    colnames(counts1) <- NULL

    i=1
    for(i in 1:ncol(exp.matrix)) {
        ## manual design matrix (CHECK THIS!!!)
        ct <- contr.matrix[,i]
        y <- factor(c("neg","zero","pos")[2+sign(ct)],levels=c("neg","zero","pos"))
        ##design1 <- model.matrix( ~ factor(y[jj]==1))
        ## sample-wise model matrix (does this work???)
        ##design.formula = formula(" ~ 0 + group")
        design.formula = formula("~ 0+y")
        dds <- DESeqDataSetFromMatrix(
            countData=counts1[,], design=design.formula, colData=data.frame(y))
        ##fitType = 'parametric'
        fitType ='local'
        if(test=="LRT") {
            suppressWarnings( dds <- DESeq(dds, fitType=fitType, test="LRT", reduced=~1) )
        } else {
            suppressWarnings( dds <- DESeq(dds, fitType=fitType, test="Wald") )
        }
        resultsNames(dds)
        ctx <- c("yneg"=-1, "yzero"=0, "ypos"=1)[resultsNames(dds)]
        resx <- results(dds, contrast=ctx, cooksCutoff=FALSE, independentFiltering=FALSE )
        
        ## we add the gene annotation here (not standard...)
        rownames(resx) <- rownames(rowData(dds))
        dim(contr.matrix)
        
        ## resx = dds.result[rownames(dds),]  ## seems not
        pos.samples = which(exp.matrix[,i] > 0)
        neg.samples = which(exp.matrix[,i] < 0)
        resx$AveExpr1 = rowMeans(X[,pos.samples,drop=FALSE], na.rm=TRUE)
        resx$AveExpr0 = rowMeans(X[,neg.samples,drop=FALSE], na.rm=TRUE)
        resx$log2BaseMean = log2(0.0001+resx$baseMean)
        head(resx)        

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
            ##genes = rowData(dds)$genes  ## BE SURE TO HAVE THIS!!
            ##tables[[i]] = cbind(genes, tables[[i]])
        }
    }

    res = list( tables=tables)
    return(res)
}


##--------------------------------------------------------------------------------------------
##------------------------------------ END OF FILE -------------------------------------------
##--------------------------------------------------------------------------------------------
