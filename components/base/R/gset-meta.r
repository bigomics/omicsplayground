##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

ALL.GENESET.METHODS = c("fisher","ssgsea","gsva", "spearman", "camera", "fry",
                        ##"plage", "cuda.gsea", "enricher",
                        "gsea.permPH","gsea.permGS","gseaPR","fgsea")
methods=ALL.GENESET.METHODS
methods=c("fisher","ssgsea","gsva","fgsea","gseaPR")
methods=c("gsva","camera")
methods=c("fisher","gsva","fgsea")
use.multicore=TRUE

mc.threads=1
if(0) {
    X=ngs$X;Y=ngs$samples;design=ngs$model.parameters$design;G=ngs$GMT
    gmt <- 
    contr.matrix=ngs$model.parameters$contr.matrix;
    mc.cores=1;mc.threads=1;batch.correct=TRUE
    methods=c("fisher","gsva","fgsea")    
}

gset.fitContrastsWithAllMethods <- function(gmt, X, Y, G, design, contr.matrix, methods, 
                                            mc.threads=1, mc.cores=NULL, batch.correct=TRUE)
{
    timings <- c()

    if(is.null(mc.cores)) {
        mc.cores <- round( 0.5 * parallel::detectCores(all.tests = TRUE, logical = FALSE) )
        mc.cores <- pmax(mc.cores,1)
        mc.cores <- pmin(mc.cores,16)
        mc.cores
    }
    cat("using",mc.cores,"number of cores\n")
    cat("using",mc.threads,"number of threads\n")

    if(methods[1]=="*") {
        methods <- ALL.GENESET.METHODS
    }
    methods <- intersect(methods, ALL.GENESET.METHODS)
    cat("calculating methods:",methods,"\n")

    ##  filter small gene sets
    cat("filtering gene sets...\n")
    gmt <- lapply( gmt, function(s) intersect(s, rownames(X)))
    gmt.size <- sapply(gmt, length)
    summary(gmt.size)
    keep = (gmt.size >= 15 & gmt.size < 99999)
    table(keep)
    gmt <- gmt[which(keep)]
    length(gmt)
    
    ## If degenerate set design to NULL
    if(!is.null(design) && ncol(design)>=ncol(X) ) {
        ## "no-replicate" design!!!
        cat("WARNING: degenerate design. setting design to NULL\n")
        contr.matrix <- design %*% contr.matrix
        design <- NULL
    }

    ## experiment matrix
    if(!is.null(design)) {
        exp.matrix = (design %*% contr.matrix)[colnames(X),,drop=FALSE]
    } else {
        exp.matrix <- contr.matrix[colnames(X),,drop=FALSE]
    }
    Y <- Y[colnames(X),,drop=FALSE]
    dim(exp.matrix)

    ## some "normalization"...
    ## remove batch-effects with LIMMA. Be sure to include batch in the
    ## model to avoid overfitting.
    if("batch" %in% colnames(Y) && batch.correct) cat("correcting for batch effects\n")
    if("nnm" %in% colnames(Y) && batch.correct) cat("correcting for NNM\n")
    ##zx=zx.gsva
    my.normalize <- function(zx, Y) {
        if("batch" %in% colnames(Y) && batch.correct) {
            nbatch <- length(unique(Y$batch))
            if(!is.null(design) && nbatch>1) {
                zx <- limma::removeBatchEffect( zx, batch=Y$batch, design=design)
            } else if(nbatch>1) {
                zx <- limma::removeBatchEffect( zx, batch=Y$batch)
            }
        }
        if("nnm" %in% colnames(Y) && batch.correct) {
            yy = Y$nnm
            ny <- length(unique(yy))
            if(ny>1) zx <- gx.nnmcorrect( zx, yy, k=3)
        }
        zx <- scale(limma::normalizeQuantiles(zx))
        return(zx)
    }
    
    all.results <- list()
    ## pre-compute matrices
    zx.gsva = zx.ssgsea = zx.rnkcorr = NULL
    res.gsva = res.ssgsea = res.rnkcorr = NULL
    methods
    
    dim(G)
    table(rownames(X) %in% rownames(G))
    table(colnames(G) %in% names(gmt))
    G <- G[rownames(X),names(gmt)]

    if("spearman" %in% methods) {

        cat("fitting contrasts using spearman/limma... \n")
        
        ## single-sample gene set enrichment using (fast) rank correlation
        xx1 <-  X - rowMeans(X,na.rm=TRUE)  ## center it...
        xx1 <- apply(xx1,2,rank,na.last="keep")  ## rank correlation (like spearman)
        jj = intersect(rownames(G),rownames(xx1))
        tt <- system.time({

            zx.rnkcorr <- qlcMatrix::corSparse(G[jj,], xx1[jj,])  ## superfast
            rownames(zx.rnkcorr) <- colnames(G)
            colnames(zx.rnkcorr) <- colnames(X)

            ## row-wise (per feature) scaling is 'good practice', see
            ## tests comparing rankcor and ssGSEA/gsva
            ##zx.rnkcorr <- t(scale(t(zx.rnkcorr))) ## ??? 2020.10.26 IK.. not sure

            ## additional batch correction and NNM???
            zx.rnkcorr <- my.normalize(zx.rnkcorr, Y)
            zx.rnkcorr <- zx.rnkcorr[names(gmt),colnames(X)] ## make sure..

            ## compute LIMMA
            all.results[["spearman"]] <- gset.fitContrastsWithLIMMA(
                zx.rnkcorr, contr.matrix,  design=design, trend=TRUE, conform.output=TRUE)

        })
        timings <- rbind(timings, c("spearman", tt))
        sum(is.na(zx.rnkcorr))
    }

    methods
    if("gsva" %in% methods) {
        cat("fitting contrasts using GSVA/limma... \n")
        tt <- system.time({
            zx.gsva <- NULL
            zx.gsva <- try( GSVA::gsva(as.matrix(X), gmt[], method="gsva",
                                 parallel.sz=mc.cores, verbose=FALSE))
            dim(zx.gsva)
            if(is.null(zx.gsva) || "try-error" %in% class(zx.gsva) ) {
                ## switch to single core...
                cat("WARNING:: GSVA ERROR : retrying single core ... \n")
                zx.gsva <- try(GSVA::gsva(as.matrix(X), gmt[], method="gsva",
                                    parallel.sz=1, verbose=FALSE))
            }
            class(zx.gsva)
            if("try-error" %in% class(zx.gsva)) {
                stop("FATAL ERROR in GSVA\n")
            }
            zx.gsva <- my.normalize(zx.gsva, Y)
            jj <- match(names(gmt), rownames(zx.gsva))
            zx.gsva <- zx.gsva[jj,colnames(X)] ## make sure..
            zx.gsva[is.na(zx.gsva)] <- 0
            all.results[["gsva"]] <- gset.fitContrastsWithLIMMA(
                zx.gsva, contr.matrix,  design=design, trend=TRUE, conform.output=TRUE)
        })
        timings <- rbind(timings, c("gsva", tt))
        sum(is.na(zx.gsva))
    }

    if("ssgsea" %in% methods) {
        cat("fitting contrasts using ssGSEA/limma... \n")
        tt <- system.time({
            zx.ssgsea <- GSVA::gsva(as.matrix(X), gmt[], method="ssgsea",
                              parallel.sz=mc.cores, verbose=FALSE)
            dim(zx.ssgsea)
            zx.ssgsea <- my.normalize(zx.ssgsea, Y)
            jj <- match(names(gmt), rownames(zx.ssgsea))
            zx.ssgsea <- zx.ssgsea[jj,colnames(X)] ## make sure..
            zx.ssgsea[is.na(zx.ssgsea)] <- 0
            all.results[["ssgsea"]] <- gset.fitContrastsWithLIMMA(
                zx.ssgsea, contr.matrix, design, trend=TRUE,conform.output=TRUE)
        })
        timings <- rbind(timings, c("ssgsea", tt))
    }

    k=1
    fitThisContrastWithMethod <- function(method, k) {

        ##cat("fitting contrast",colnames(contr.matrix)[k],"\n")
        jj = which(exp.matrix[,k]!=0)
        yy <- 1 * (exp.matrix[jj,k] > 0)
        xx <- X[,jj]
        dim(xx)
        table(yy)
        ##ref = names(which(contr.matrix[,k] < 0))
        ##ref = rownames(contr.matrix)[which(contr.matrix[,k] < 0)]
        ref = 0
        ref

        timings <- c()
        res = list()

        ## Standard Fisher exact test
        if("fisher" %in% method) {

            ## calculate significant genes with LIMMA (we need all genes for GSEA-PR)
            lfc = 0
            lfc05=0.2; fdr=0.25   ## OLD thresholds 
            lfc05=0.0; fdr=0.05  ## NEW thresholds (since oct2021)
            suppressWarnings( suppressMessages(
                limma0 <- gx.limma( xx, yy, fdr=1.0, lfc=0,
                                   ref=ref, trend=TRUE, verbose=0)  ## trend true for NGS
            ))
            which.up = which(limma0[,"adj.P.Val"] <= fdr & limma0[,"logFC"] > lfc05)
            which.dn = which(limma0[,"adj.P.Val"] <= fdr & limma0[,"logFC"] < -lfc05)
            ## which.up = which(limma0[,"P.Value"] < 0.05 & limma0[,"logFC"] > lfc05)
            ## which.dn = which(limma0[,"P.Value"] < 0.05 & limma0[,"logFC"] < -lfc05)
            genes.up = rownames(limma0)[which.up]
            genes.dn = rownames(limma0)[which.dn]
            ##cat("siggenes.down=",length(genes.dn),"\n")
            ##cat("siggenes.up=",length(genes.up),"\n")

            ## Always take at least first 100.. (HACK??!!!)
            if(length(genes.dn) < 100) {
                genes.dn0 <-  rownames(limma0)[order(limma0[,"logFC"])]
                genes.dn <- head(unique(c(genes.dn,genes.dn0)),100)
            }
            if(length(genes.up) < 100) {
                genes.up0 <-  rownames(limma0)[order(-limma0[,"logFC"])]
                genes.up <- head(unique(c(genes.up,genes.up0)),100)
            }
            
            ##cat("fisher: testing...\n")
            tt <- system.time({
                output <- gset.fisher2(genes.up, genes.dn, genesets=gmt, fdr=1.0,
                                       background = rownames(X), check.background=FALSE,
                                       min.genes=0, max.genes=99999,
                                       common.genes=FALSE, verbose=0)
            })

            timings <- rbind(timings, c("fisher", tt))
            Matrix::head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <-  names(gmt)
            output <- output[,c("sign","p.value","q.value","odd.ratio","overlap")]
            colnames(output) <- c("score","p.value","q.value","odd.ratio","overlap")
            res[["fisher"]] <- output
        }

        ##----------------------------------------------------
        ## GSVA-limma
        ##----------------------------------------------------
        LIMMA.TREND=FALSE
        LIMMA.TREND=TRUE
        if("ssgsea" %in% method) {
            
            zx <- zx.ssgsea[,colnames(xx)]
            gs <- intersect(names(gmt),rownames(zx))
            tt <- system.time(
                output <-  gx.limma( zx[gs,], yy, fdr=1, lfc=0, ref=ref, trend=LIMMA.TREND, verbose=0)  ## ssgsea
            )
            timings <- rbind(timings, c("ssgsea", tt))
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <-  names(gmt)
            output <- output[,c("logFC","P.Value","adj.P.Val","0","1")]
            colnames(output) <- c("score","p.value","q.value","AveExpr0","AveExpr1")
            res[["ssgsea"]] = output
        }

        if("gsva" %in% method) {
            zx <- zx.gsva[,colnames(xx)]
            gs <- intersect(names(gmt),rownames(zx))
            tt <- system.time({
                output <-  gx.limma( zx[gs,], yy, fdr=1, lfc=0, ref=ref,
                                    trend=LIMMA.TREND, verbose=0)  ## ssgsea                
            })
            timings <- rbind(timings, c("gsva", tt))
            Matrix::head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <-  names(gmt)
            output <- output[,c("logFC","P.Value","adj.P.Val","0","1")]
            colnames(output) <- c("score","p.value","q.value","AveExpr0","AveExpr1")
            res[["gsva"]] = output
        }

        if("spearman" %in% method && !is.null(zx.rnkcorr) ) {
            tt <- system.time(
                output <-  gx.limma( zx.rnkcorr[,], yy, fdr=1, lfc=0, ref=ref, trend=LIMMA.TREND, verbose=0)  ## ssgsea
            )
            timings <- rbind(timings, c("spearman", tt))
            Matrix::head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <-  names(gmt)
            output <- output[,c("logFC","P.Value","adj.P.Val","0","1")]
            colnames(output) <- c("score","p.value","q.value","AveExpr0","AveExpr1")
            res[["spearman"]] = output
        }

        ##----------------------------------------------------
        ## LIMMA methods
        ##----------------------------------------------------
        if("camera" %in% method) {
            
            
            cdesign <- cbind(Intercept=1,Group=yy)
            ##design <- model.matrix( ~ 0 + as.factor(yy))
            ##colnames(design) <- levels(yy)
            tt <- system.time({
                suppressWarnings( suppressMessages(
                    output <- limma::camera(xx, gmt, cdesign, contrast=2)
                ))
            })
            timings <- rbind(timings, c("camera", tt))
            ## note: camera does not provide any score!!!!
            output$score = c(-1,1)[ 1 + 1*(output$Direction=="Up")] * -log10(output$PValue)
            Matrix::head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <-  names(gmt)
            output <- output[,c("score","PValue","FDR","NGenes","Direction")]
            colnames(output) <- c("score","p.value","q.value","NGenes","Direction")
            res[["camera"]] = output
        }
        if("fry" %in% method) {
            
            cdesign <- cbind(Intercept=1,Group=yy)
            ##design <- model.matrix( ~ 0 + as.factor(yy))
            ##colnames(design) <- levels(yy)
            tt <- system.time(
                output <- limma::fry(xx, gmt, cdesign, contrast=2)
            )
            timings <- rbind(timings, c("fry", tt))
            ## note: camera does not provide any logFC!!!!
            output$score = c(-1,1)[ 1 + 1*(output$Direction=="Up")] * -log10(output$PValue)
            Matrix::head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <-  names(gmt)
            output <- output[,c("score","PValue","FDR","NGenes","Direction")]
            colnames(output) <- c("score","p.value","q.value","NGenes","Direction")
            res[["fry"]] = output
        }

        ##----------------------------------------------------
        ## GSEA methods
        ##----------------------------------------------------
        if("gsea.permPH" %in% method) {
            tt <- system.time(
                output <- run.GSEA( xx, yy, gmt, fdr=1.0, do.leading.edge = FALSE,
                                  set.min=0, set.max=99999, ref.type=ref, permute="phenotype")
            )
            timings <- rbind(timings, c("gsea.permPH", tt))
            dim(output)
            rownames(output) = output$GS
            Matrix::head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <- names(gmt)
            output <- output[,c("NES","NOM p-val","FDR q-val","NES","SIZE","LEADING EDGE","LEADING GENES")]
            colnames(output) <- c("score","p.value","q.value","NES","SIZE","LEADING EDGE","LEADING GENES")
            res[["gsea.permPH"]] = output
        }

        if("gsea.permGS" %in% method) {
            tt <- system.time(
                output <- run.GSEA( xx, yy, gmt, fdr=1.0, do.leading.edge = FALSE,
                                   set.min=0, set.max=99999, ref.type=ref, permute="gene_set")
            )
            timings <- rbind(timings, c("gsea.permGS", tt))
            dim(output)
            rownames(output) = output$GS
            Matrix::head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <- names(gmt)
            output <- output[,c("NES","NOM p-val","FDR q-val","NES","SIZE","LEADING EDGE","LEADING GENES")]
            colnames(output) <- c("score","p.value","q.value","NES","SIZE","LEADING EDGE","LEADING GENES")
            res[["gsea.permGS"]] = output
        }

        ## GSEA preranked
        if("gseaPR" %in% method) {
            ##rnk = limma0[,"t"]  ## or FC???
            ##rnk = limma0[,"logFC"]  ## or FC???
            rnk = rowMeans(xx[,which(yy==1),drop=FALSE]) - rowMeans(xx[,which(yy==0),drop=FALSE])
            tt <- system.time(
                output <- run.GSEA.preranked( rnk, gmt, fdr=1.0, do.leading.edge = FALSE,
                                             set.min=0, set.max=99999, output.dir=NULL)
            )
            timings <- rbind(timings, c("gseaPR", tt))
            dim(output)
            rownames(output) = output$GS
            Matrix::head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <- names(gmt)
            output <- output[,c("NES","NOM p-val","FDR q-val","NES","SIZE","LEADING EDGE","LEADING GENES")]
            colnames(output) <- c("score","p.value","q.value","NES","SIZE","LEADING EDGE","LEADING GENES")
            res[["gseaPR"]] = output
        }

        ## fast GSEA
        if("fgsea" %in% method) {
            
            ##rnk = limma0[,"t"]  ## or FC???
            ##rnk = limma0[,"logFC"]  ## or FC???
            rnk = rowMeans(xx[,which(yy==1),drop=FALSE]) - rowMeans(xx[,which(yy==0),drop=FALSE])
            rnk <- rnk + 1e-8*rnorm(length(rnk))
            ##output <- fgsea::fgsea(gmt[1:2], rnk, nperm=1001, nproc=1)
            tt <- system.time(                
                output <- fgsea::fgseaSimple(gmt, rnk, nperm=10000,
                                minSize=1, maxSize=9999, nproc=1) ## nproc0 fails!!!
            )
            timings <- rbind(timings, c("fgsea", tt))
            names(output)
            output = as.data.frame(output)
            output = output[match(names(gmt),output$pathway),]
            rownames(output) = names(gmt)
            Matrix::head(output)
            dim(output)
            output <- output[,c("NES","pval","padj","size","leadingEdge")]
            colnames(output) <- c("score","p.value","q.value","SIZE","LEADING GENES")
            res[["fgsea"]] = output
        }

        res2 <- list(results=res, timings=timings)
        return(res2)
    }

    method="fisher"
    method="fgsea"
    method="gsva"
    fitContrastsWithMethod <- function(method) {
        cat("fitting contrasts using",method,"... \n")
        results <- list()
        timings <- c()
        k=1
        ncontrasts <- ncol(contr.matrix)
        for(k in 1:ncontrasts) {
            res <- fitThisContrastWithMethod(method=method, k)
            results[[k]] <- res$results[[1]]
            timings <- rbind( timings, res$timings)
            names(results)[k] <- colnames(contr.matrix)[k]
        }
        return( list(results=results, timings=timings) )
    }

    ##--------------------------------------------------------------
    ## Fit remaining methods
    ##--------------------------------------------------------------
    names(all.results)
    methods2 = setdiff(methods, names(all.results))
    methods2
    m="camera"
    m="fgsea"
    m="fisher"
    for(m in methods2) {
        res <- fitContrastsWithMethod(method=m)
        all.results[[m]] <- res$results
        timings <- rbind( timings, res$timings)
    }
    
    ##--------------------------------------------------------------
    ## Reshape matrices by comparison
    ##--------------------------------------------------------------
    names(all.results)
    cat("[gset.fitContrastsWithAllMethods] length(all.results)=",length(all.results),"\n")
    tests = names(all.results[[1]])
    ntest = length(tests)
    P = lapply(tests, function(k) sapply( all.results, function(x) x[[k]][,"p.value"]))
    Q = lapply(tests, function(k) sapply( all.results, function(x) x[[k]][,"q.value"]))
    S = lapply(tests, function(k) sapply( all.results, function(x) x[[k]][,"score"]))

    for(i in 1:ntest) {
        rownames(P[[i]]) <- names(gmt)
        rownames(Q[[i]]) <- names(gmt)
        rownames(S[[i]]) <- names(gmt)
    }
    names(P) = names(Q) = names(S) = tests

    ##--------------------------------------------------------------
    ## Compute sig counts (by method)
    ##--------------------------------------------------------------
    cat("computing sigcounts... \n")
    
    methods = colnames(Q[[1]])
    nmethod = length(methods)
    nmethod
    pv.list <-  sort( c(1e-16, 10**seq(-8,-2,2), 0.05, 0.1, 0.2, 0.5,1))
    sig.counts <- list()
    i=1;fdr=0.05;
    lfc = 1e-3
    for(i in 1:nmethod) {
        ##p0 <- sapply( P, function(x) x[,1])
        q0 <- sapply( Q, function(x) x[,i])
        s0 <- sapply( S, function(x) x[,i])
        q0[is.na(q0)] <- 1
        s0[is.na(s0)] <- 0
        up0 <- sapply( pv.list, function(p) Matrix::colSums( q0 <= p & s0 > 0, na.rm=TRUE))
        dn0 <- sapply( pv.list, function(p) Matrix::colSums( q0 <= p & s0 < 0, na.rm=TRUE))
        ns0 <- sapply( pv.list, function(p) Matrix::colSums( q0 > p, na.rm=TRUE))
        both0 <- sapply( pv.list, function(p) Matrix::colSums( q0 <= p  & abs(s0) > 0, na.rm=TRUE))
        if(ncol(q0)==1) {
            up0 = matrix(up0, nrow=1)
            dn0 = matrix(dn0, nrow=1)
            ns0 = matrix(ns0, nrow=1)
            both0 = matrix(both0, nrow=1)
            rownames(up0) <- rownames(dn0) <- rownames(ns0) <- rownames(both0) <- colnames(q0)[1]
        }
        colnames(up0) <- colnames(dn0) <- colnames(ns0) <- colnames(both0) <- pv.list
        m = methods[i]
        sig.counts[[m]] = list(both=both0, up=up0, down=dn0, notsig=ns0)
    }

    ##--------------------------------------------------
    ## meta analysis, aggregate p-values
    ##--------------------------------------------------
    cat("computing meta-p values... \n")

    all.meta <- list()
    i=1
    for(i in 1:ntest) {
        pv = P[[i]]
        qv = Q[[i]]
        fc = S[[i]]
        meta.p = apply(pv, 1, max, na.rm=TRUE ) ## maximum p-statistic (simple & fast)
        meta.q = apply(qv, 1, max, na.rm=TRUE ) ## maximum q-statistic (simple & fast)
        ## meta.p = apply(pv, 1, function(p) metap::allmetap(p, method="sumlog")$p[[1]])
        ## meta.q = p.adjust(meta.p, method="fdr")
        ss.rank <- function(x) scale(sign(x)*rank(abs(x),na.last="keep"),center=FALSE)
        meta.fx = rowMeans( apply(S[[i]], 2, ss.rank), na.rm=TRUE)
        meta = data.frame(fx=meta.fx, p=meta.p, q=meta.q)
        rownames(fc) <- NULL  ## saves memory...
        rownames(pv) <- NULL
        rownames(qv) <- NULL
        all.meta[[i]] = data.frame(meta=meta, fc=I(fc), p=I(pv), q=I(qv))
        rownames(all.meta[[i]]) <- rownames(S[[i]])
    }
    names(all.meta) = tests

    ##--------------------------------------------------
    ## Add meta matrices (this becomes quite large...)
    ##--------------------------------------------------
    cat("computing meta-matrix... \n")
    
    m <- list(gsva=zx.gsva, ssgsea=zx.ssgsea, rnkcorr=zx.rnkcorr)        
    m = m[which(!sapply(m,is.null))]
    names(m)

    if(0) {
        if(length(m)>1) {
            ## Average normalized single-sample values
            ##m <- lapply(m, function(x) apply(x,2,rank,na.last="keep"))
            m = lapply(m, function(x) scale(x,center=FALSE))
            avg.m  <- Reduce('+',m) / length(m)  ## matrix mean
            meta.matrix <- scale(avg.m,center=FALSE)
            meta.matrix <- limma::normalizeQuantiles(meta.matrix)
        } else {
            meta.matrix <- m[[1]]
        }
        m[["meta"]] <- meta.matrix
    } else {
        ## average expression of geneset members
        ng <- Matrix::colSums(G!=0)
        meta.matrix <- as.matrix(Matrix::t(G!=0) %*% X) / ng
    }
    ## meta.matrix <- meta.matrix - rowMeans(meta.matrix,na.rm=TRUE)  ## center??
    m[["meta"]] <- meta.matrix
    
    ##timings0 <- do.call(rbind, timings)
    timings <- as.matrix(timings)
    rownames(timings) <- timings[,1]
    timings0 <- matrix(timings[,-1],nrow=nrow(timings))
    timings0 <- matrix(as.numeric(timings0),nrow=nrow(timings0))
    rownames(timings0) <- rownames(timings)
    if(nrow(timings0)>1 && sum(duplicated(rownames(timings0))>0) ) {
        ## timings0 <- apply(timings0, 2, function(x) tapply(x,rownames(timings0),sum))
        timings0 <- do.call(rbind,tapply(1:nrow(timings0),rownames(timings0),function(i) colSums(timings0[i,,drop=FALSE])))        
    }

    res = list( meta = all.meta, sig.counts = sig.counts,  outputs = all.results,
               matrices = m, timings = timings0)

    return(res)
}

##trend=TRUE;gsetX=zx.gsva;conform.output=TRUE
##trend=TRUE;gsetX=zx.rnkcorr;conform.output=TRUE
gset.fitContrastsWithLIMMA <- function( gsetX, contr.matrix, design,
                                       trend=TRUE,conform.output=FALSE)
{
    if(!is.null(design)) {
        cat("fitting gset.LIMMA contrasts with design matrix....\n")

        ##xfit = limma::normalizeQuantiles(xfit)
        vfit <- limma::lmFit(gsetX, design)
        vfit <- limma::contrasts.fit(vfit, contrasts=contr.matrix)
        efit <- limma::eBayes(vfit, trend=trend, robust=TRUE)
        ##efit <- limma::eBayes(vfit, trend=trend, robust=FALSE)

        tables <- list()
        i=1
        exp.matrix = (design %*% contr.matrix)
        for(i in 1:ncol(contr.matrix)) {
            ##coef = colnames(contr.matrix)[i]
            top = limma::topTable(efit, coef=i, sort.by="none",number=Inf, adjust.method="BH")
            j1 = which( exp.matrix[,i] > 0 )
            j0 = which( exp.matrix[,i] < 0 )
            ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
            mean1 = rowMeans(gsetX[,j1,drop=FALSE], na.rm=TRUE)
            mean0 = rowMeans(gsetX[,j0,drop=FALSE], na.rm=TRUE)
            top = top[rownames(gsetX),]
            top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1)
            Matrix::head(top,10)
            tables[[i]] = top
        }
        names(tables) <- colnames(contr.matrix)
    } else {
        cat("fitting gset.LIMMA contrasts without design....\n")

        ##trend=TRUE
        tables <- list()
        i=1
        for(i in 1:ncol(contr.matrix)) {
            design0 <- cbind( 1, contr.matrix[,i])
            colnames(design0) <- c( "ref", colnames(contr.matrix)[i])
            vfit <- limma::lmFit(gsetX, design0)
            efit <- limma::eBayes(vfit, trend=trend, robust=TRUE)
            top = limma::topTable(efit, coef=2, sort.by="none",number=Inf, adjust.method="BH")
            j1 = which( contr.matrix[,i] > 0 )
            j0 = which( contr.matrix[,i] < 0 )
            ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
            mean1 = rowMeans(gsetX[,j1,drop=FALSE], na.rm=TRUE)
            mean0 = rowMeans(gsetX[,j0,drop=FALSE], na.rm=TRUE)
            top = top[rownames(gsetX),]
            top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1)
            Matrix::head(top,10)
            tables[[i]] = top
        }
        names(tables) <- colnames(contr.matrix)
    }

    if(conform.output==TRUE) {
        for(i in 1:length(tables)) {
            ##k1 = c("logFC","AveExpr","t","P.Value","adj.P.Val","AveExpr1","AveExpr2")
            ##k2 = c("logFC","AveExpr","statistic","P.Value","adj.P.Val","AveExpr1","AveExpr2")
            jj <- match(rownames(gsetX),rownames(tables[[i]]))
            k1 <- c("logFC","P.Value","adj.P.Val","AveExpr0","AveExpr1")
            k2 <- c("score","p.value","q.value","AveExpr0","AveExpr1")
            tables[[i]] = tables[[i]][jj,k1]
            colnames(tables[[i]]) = k2
            ## tables[[i]] = cbind( rownames(gsetX), tables[[i]])
        }
    }

    return(tables)
}



##======================================================================
##======================= GSET METHODS =================================
##======================================================================

##path="output_GSEA"
getGeneSetTables <- function(path) {
    dd = dir(path,full.names=TRUE)
    dd
    tables = vector("list",length(dd))
    d=dd[1]
    i=1
    for(i in 1:length(dd)) {
        d = dd[i]
        table_dir = file.path(d,"tables")
        if(!file.exists(table_dir)) next
        ff = dir(table_dir,full.names=TRUE)
        ff.short = dir(table_dir,full.names=FALSE)
        ff.name = gsub("output-|-results.csv|-test.csv","",ff.short)
        tables[[i]] = lapply(ff, read.csv, row.names=1,check.names=FALSE)
        names(tables[[i]]) = ff.name
        names(tables)[i] = gsub(".*/","",d)
    }

    gsets = sort(unique(unlist(lapply(tables[[1]], rownames))))
    length(gsets)

    meta <- list()
    i=1
    for(i in 1:length(dd)) {
        pv = c()
        fx = c()
        j=1
        for(j in 1:length(tables[[i]])) {
            x = tables[[i]][[j]]
            pv.col = grep("p.value|^p$|p-val|pval|nom.p.val|nom p-val|p.val",tolower(colnames(x)))[1]
            pv = cbind(pv, x[match(gsets,rownames(x)),pv.col] )
            fx.col = grep("sign|nes|logfc|fc",tolower(colnames(x)))[1]
            fx = cbind(fx, (x[match(gsets,rownames(x)),fx.col] ))
            colnames(pv)[ncol(pv)] = names(tables[[i]])[j]
            colnames(fx)[ncol(fx)] = names(tables[[i]])[j]
        }
        Matrix::head(pv)
        rownames(pv) = gsets
        rownames(fx) = gsets
        meta[[i]] = data.frame(geneset=gsets, fx=fx, p=pv)
    }
    names(meta) = names(tables)
    p

    res = list(tables=tables, meta=meta)
    return(res)
}


##======================================================================
##======================= GSEA METHODS =================================
##======================================================================

gmt2mat.nocheck <- function(gmt, bg=NULL, use.multicore=TRUE)
{
    ##max.genes=-1;ntop=-1;sparse=TRUE;bg=NULL;normalize=FALSE;r=0.01;use.multicore=TRUE
    
    
    ##gmt <- gmt[!duplicated(names(gmt))]
    if(is.null(bg)) {
        bg <- names(sort(table(unlist(gmt)),decreasing=TRUE))
    }
    gmt <- lapply(gmt, function(s) intersect(bg,s))
    ##kk <- unique(names(gmt))
    ## D <- Matrix::Matrix(0, nrow=length(bg),ncol=length(kk), sparse=TRUE)
    ##D <- Matrix::sparseMatrix(1, 1, x=0, dims=c(length(bg),ncol=length(kk)))
    ##rownames(D) <- bg
    ##colnames(D) <- kk
    j=1
    idx <- c()
    if(use.multicore) {
        idx <- parallel::mclapply(gmt, function(s) match(s,bg))
        idx[sapply(idx,length)==0] <- 0
        idx <- sapply(1:length(idx), function(i) rbind(idx[[i]],i))
        idx <- matrix(unlist(idx[]),byrow=TRUE,ncol=2)
        idx <- idx[!is.na(idx[,1]),]
        idx <- idx[idx[,1]>0,]
        ##D[idx] <- 1
    } else {
        idx <- c()
        for(j in 1:length(gmt)) {
            ii0 <- which(bg %in% gmt[[j]])
            if(length(ii0)>0) {
                ##D[ii0,j] <- +1
                idx <- rbind(idx, cbind(ii0,j))
            }
        }
    }
    D <- Matrix::sparseMatrix(idx[,1], idx[,2], x=rep(1,nrow(idx)),
                      dims=c(length(bg),ncol=length(gmt)))
    dim(D)
    rownames(D) <- bg
    colnames(D) <- names(gmt)
    ##D <- Matrix::Matrix(D, sparse=TRUE)
    D
}

shortstring = function(s,n) {
    s=as.character(s);
    ifelse(nchar(s)<=n,s,paste0(substring(s,1,n),"..."))
}

getGseaOutputDir <- function(path) {
    ## untangle Gsea subfolder
    gsea_dir = dir(path)[grep("\\.Gsea\\.",dir(path))]
    gsea_dir = file.path(path,gsea_dir)
    gsea_dir
    if(length(gsea_dir)==0) {
        cat("Missing results")
        return(NULL)
    }
    return(gsea_dir)
}

getGseaTable <- function(path) {
    ff = dir(path,full.names=TRUE)
    report_name = ff[grep("gsea_report.txt$",ff)]
    if(is.null(report_name) || length(report_name)==0) return(NULL)
    R = read.csv(report_name,sep="\t",check.names=FALSE)
    R = R[order(-abs(R$NES)),]
    colnames(R)[1] = "NAME"
    R$NES = round(R$NES, digits=3)
    kk = grep("NAME|SIZE|NES|NOM|FDR|LEADING",colnames(R),ignore.case=TRUE)
    R = R[,kk]
    ##colnames(R)[4:5] = c("p","q")
    return(R)
}

gseaSnapshot <- function(gsets, gsea_dir) {
    
    
    
    enplots = dir(gsea_dir, pattern="enplot_")
    enplots0 = dir(gsea_dir, pattern="enplot_", full.names=TRUE)
    kk = match(gsets, gsub("enplot_|_[0-9]*.png$","",enplots))
    imgs = lapply( enplots0[kk], function(p) grid::rasterGrob(as.raster(png::readPNG(p)),
                                                        interpolate=TRUE))
    grid.arrange( grobs=imgs, ncol=5)
}


##======================================================================
##======================================================================
##======================================================================
