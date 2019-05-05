

ALL.GENESET.METHODS=c("fisher","ssgsea","gsva", "spearman", "camera", "fry",
                  ##"plage", "cuda.gsea", "enricher",
                  "gsea.permPH","gsea.permGS","gseaPR","fgsea")
methods=ALL.GENESET.METHODS
methods=c("fisher","ssgsea","gsva","fgsea","gseaPR")
methods=c("gsva","camera")
methods=c("fisher","gsva","camera")

gmt2mat.nocheck <- function(gmt, bg=NULL, use.multicore=TRUE)
{
    ##max.genes=-1;ntop=-1;sparse=TRUE;bg=NULL;normalize=FALSE;r=0.01;use.multicore=TRUE
    require(Matrix)
    require(parallel)
    if(is.null(bg)) {
        bg <- names(sort(table(unlist(gmt)),decreasing=TRUE))
    }
    gmt <- lapply(gmt, function(s) intersect(bg,s))
    kk <- unique(names(gmt))
    D <- Matrix(0, nrow=length(bg),ncol=length(kk), sparse=TRUE)
    dim(D)
    rownames(D) <- bg
    colnames(D) <- kk
    j=1
    if(use.multicore) {
        idx <- mclapply(gmt, function(s) match(s,bg))
        idx[sapply(idx,length)==0] <- 0
        idx <- sapply(1:length(idx), function(i) rbind(idx[[i]],i))
        idx <- matrix(unlist(idx[]),byrow=TRUE,ncol=2)
        idx <- idx[!is.na(idx[,1]),]
        idx <- idx[idx[,1]>0,]
        D[idx] <- 1
    } else {
        for(j in 1:ncol(D)) {
            k0 <- match(kk[j], names(gmt))
            ii0 <- which(bg %in% gmt[[k0]])
            if(length(ii0)>0) D[ii0,j] <- +1
        }
    }
    D
}

##X=ngs$X;Y=ngs$Y;design=ngs$model.parameters$design;contr.matrix=ngs$model.parameters$contr.matrix
##mc.cores=4;mc.threads=1;batch.correct=TRUE;gmt=ngs$gmt.all
gset.fitContrastsWithAllMethods <- function(gmt, X, Y, design, contr.matrix, methods,
                                            mc.threads=1, mc.cores=NULL, batch.correct=TRUE)
{
    require(GSVA)
    require(Matrix)
    require(qlcMatrix)
    timings <- c()

    if(is.null(mc.cores)) {
        mc.cores <- round( 0.5*detectCores(all.tests = TRUE, logical = FALSE) )
        mc.cores <- pmax(mc.cores,1)
        mc.cores <- pmin(mc.cores,16)
        mc.cores
    }
    cat("using",mc.cores,"number of cores\n")
    cat("using",mc.threads,"number of threads\n")

    if(methods[1]=="*") {
        methods <- ALL.GENESET.METHODS
    }
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

    ## pre-compute the big GMT matrix (again??)
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
    normalize <- function(zx, Y) {
        if("batch" %in% colnames(Y) && batch.correct) {
            ##design0 = model.matrix( ~ Y$group )
            design0 = design
            zx <- removeBatchEffect( zx, batch=Y$batch, design=design0)
        }
        if("nnm" %in% colnames(Y) && batch.correct) {
            yy = Y$nnm
            zx = gx.nnmcorrect( zx, yy, k=3)
        }
        zx <- scale(normalizeQuantiles(zx))
        return(zx)
    }

    all.results <- list()
    ## pre-compute matrices
    zx.gsva = zx.ssgsea = zx.rnkcorr = NULL
    res.gsva = res.ssgsea = res.rnkcorr = NULL
    methods

    if("spearman" %in% methods) {

        cat("fitting contrasts using spearman/limma... \n")
        require(qlcMatrix)
        ##GMT = Matrix(sapply( gmt[], function(s) 1*(rownames(X) %in% s)),sparse=TRUE)
        ##rownames(GMT) = rownames(X)
        GMT <- gmt2mat.nocheck(gmt[], bg=rownames(X))  ## in gset-gsea.r
        dim(GMT)
        table(rownames(X) %in% rownames(GMT))
        table(colnames(GMT) %in% names(gmt))
        table(colnames(GMT)==names(gmt))
        GMT <- GMT[rownames(X),names(gmt)]

        ## single-sample gene set enrichment using (fast) rank correlation
        xx1 <-  X - rowMeans(X,na.rm=TRUE)
        xx1 <- apply(xx1,2,rank,na.last="keep")
        ##G = sapply( gmt[], function(s) 1*(gg %in% s))
        jj = intersect(rownames(GMT),rownames(xx1))
        tt <- system.time({

            zx.rnkcorr <- qlcMatrix::corSparse(GMT[jj,], xx1[jj,])
            rownames(zx.rnkcorr) <- colnames(GMT)
            colnames(zx.rnkcorr) <- colnames(X)
            ## row-wise (per feature) scaling is 'good practice', see
            ## tests comparing rankcor and ssGSEA/gsva
            zx.rnkcorr <- t(scale(t(zx.rnkcorr)))

            ## additional batch correction and NNM
            zx.rnkcorr <- normalize(zx.rnkcorr, Y)
            zx.rnkcorr <- zx.rnkcorr[names(gmt),colnames(X)] ## make sure..

    ## compute LIMMA
            all.results[["spearman"]] <- gset.fitContrastsWithLIMMA(
                zx.rnkcorr, contr.matrix,  design=design, trend=TRUE, conform.output=TRUE)
        })
        timings <- rbind(timings, c("spearman", tt))
        sum(is.na(zx.rnkcorr))
    }

    if("gsva" %in% methods) {
        cat("fitting contrasts using GSVA/limma... \n")
        tt <- system.time({
            zx.gsva <- NULL
            zx.gsva <- try(gsva(as.matrix(X), gmt[], method="gsva", parallel.sz=mc.cores))
            if(is.null(zx.gsva) || class(zx.gsva)=="try-error") {
                ## switch to single core...
                cat("WARNING:: GSVA ERROR : retrying single core ... \n")
                zx.gsva <- try(gsva(as.matrix(X), gmt[], method="gsva", parallel.sz=1))
            }
            class(zx.gsva)
            if(class(zx.gsva)=="try-error") {
                stop("FATAL ERROR in GSVA\n")
            }
            zx.gsva <- normalize(zx.gsva, Y)
            zx.gsva <- zx.gsva[names(gmt),colnames(X)] ## make sure..

            all.results[["gsva"]] <- gset.fitContrastsWithLIMMA(
                zx.gsva, contr.matrix,  design=design, trend=TRUE, conform.output=TRUE)
        })
        timings <- rbind(timings, c("gsva", tt))
        sum(is.na(zx.gsva))
    }

    if("ssgsea" %in% methods) {
        cat("fitting contrasts using ssGSEA/limma... \n")
        tt <- system.time({
            zx.ssgsea <- gsva(as.matrix(X), gmt[], method="ssgsea", parallel.sz=mc.cores )
            dim(zx.ssgsea)
            zx.ssgsea <- normalize(zx.ssgsea, Y)
            zx.ssgsea <- zx.ssgsea[names(gmt),colnames(X)] ## make sure..
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
        ref = names(which(contr.matrix[,k] < 0))
        ref = 0
        ref

        timings <- c()
        res = list()

        ## Standard Fisher exact test
        if("fisher" %in% method) {
            ## calculate significant genes with LIMMA (we need all gene for GSEA-PR)
            lfc = 0
            lfc05 = 0.2  ## for genes
            fdr = 0.25
            limma0 = gx.limma( xx, yy, fdr=1.0, lfc=0, ref=ref, trend=TRUE, verbose=0)  ## trend true for NGS
            which.up = which(limma0[,"adj.P.Val"] <= fdr & limma0[,"logFC"] > lfc05)
            which.dn = which(limma0[,"adj.P.Val"] <= fdr & limma0[,"logFC"] < -lfc05)
            ## which.up = which(limma0[,"P.Value"] < 0.05 & limma0[,"logFC"] > lfc05)
            ## which.dn = which(limma0[,"P.Value"] < 0.05 & limma0[,"logFC"] < -lfc05)
            genes.up = rownames(limma0)[which.up]
            genes.dn = rownames(limma0)[which.dn]
            ##cat("siggenes.down=",length(genes.dn),"\n")
            ##cat("siggenes.up=",length(genes.up),"\n")

            ## at least take first 20
            if(length(genes.dn) < 20) {
                genes.dn <-  head(rownames(limma0)[order(limma0[,"logFC"])],20)
            }
            if(length(genes.up) < 20) {
                genes.up <-  head(rownames(limma0)[order(-limma0[,"logFC"])],20)
            }
            ##cat("siggenes.down=",length(genes.dn),"\n")
            ##cat("siggenes.up=",length(genes.up),"\n")
            tt <- system.time({
                output <- gset.fisher2(genes.up, genes.dn, genesets=gmt, fdr=1.0,
                                       background = rownames(X), check.background=FALSE,
                                       min.genes=0, max.genes=99999,
                                       common.genes=FALSE, verbose=0)
            })

            timings <- rbind(timings, c("fisher", tt))
            head(output)
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
            require(GSVA)
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
            tt <- system.time(
                output <-  gx.limma( zx[gs,], yy, fdr=1, lfc=0, ref=ref,
                                    trend=LIMMA.TREND, verbose=0)  ## ssgsea
            )
            timings <- rbind(timings, c("gsva", tt))
            head(output)
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
            head(output)
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
            require(limma)
            require(limma)
            cdesign <- cbind(Intercept=1,Group=yy)
            ##design <- model.matrix( ~ 0 + as.factor(yy))
            ##colnames(design) <- levels(yy)
            tt <- system.time(
                output <- camera(xx, gmt, cdesign, contrast=2)
            )
            timings <- rbind(timings, c("camera", tt))
            ## note: camera does not provide any score!!!!
            output$score = c(-1,1)[ 1 + 1*(output$Direction=="Up")] * -log10(output$PValue)
            head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <-  names(gmt)
            output <- output[,c("score","PValue","FDR","NGenes","Direction")]
            colnames(output) <- c("score","p.value","q.value","NGenes","Direction")
            res[["camera"]] = output
        }
        if("fry" %in% method) {
            require(limma)
            cdesign <- cbind(Intercept=1,Group=yy)
            ##design <- model.matrix( ~ 0 + as.factor(yy))
            ##colnames(design) <- levels(yy)
            tt <- system.time(
                output <- fry(xx, gmt, cdesign, contrast=2)
            )
            timings <- rbind(timings, c("fry", tt))
            ## note: camera does not provide any logFC!!!!
            output$score = c(-1,1)[ 1 + 1*(output$Direction=="Up")] * -log10(output$PValue)
            head(output)
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
            head(output)
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
            head(output)
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
            head(output)
            dim(output)
            output = output[match(names(gmt),rownames(output)),]
            rownames(output) <- names(gmt)
            output <- output[,c("NES","NOM p-val","FDR q-val","NES","SIZE","LEADING EDGE","LEADING GENES")]
            colnames(output) <- c("score","p.value","q.value","NES","SIZE","LEADING EDGE","LEADING GENES")
            res[["gseaPR"]] = output
        }

        ## fast GSEA
        if("fgsea" %in% method) {
            require(fgsea)
            ##rnk = limma0[,"t"]  ## or FC???
            ##rnk = limma0[,"logFC"]  ## or FC???
            rnk = rowMeans(xx[,which(yy==1),drop=FALSE]) - rowMeans(xx[,which(yy==0),drop=FALSE])
            rnk <- rnk + 1e-8*rnorm(length(rnk))
            ##output <- fgsea(gmt[1:2], rnk, nperm=1001, nproc=1)
            tt <- system.time(                
                output <- fgsea(gmt, rnk, nperm=10000,
                                minSize=1, maxSize=9999, nproc=1) ## nproc0 fails!!!
            )
            timings <- rbind(timings, c("fgsea", tt))
            names(output)
            output = as.data.frame(output)
            output = output[match(names(gmt),output$pathway),]
            rownames(output) = names(gmt)
            head(output)
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
    require(metap)
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
        up0 <- sapply( pv.list, function(p) colSums( q0 <= p & s0 > 0, na.rm=TRUE))
        dn0 <- sapply( pv.list, function(p) colSums( q0 <= p & s0 < 0, na.rm=TRUE))
        ns0 <- sapply( pv.list, function(p) colSums( q0 > p, na.rm=TRUE))
        both0 <- sapply( pv.list, function(p) colSums( q0 <= p  & abs(s0) > 0, na.rm=TRUE))
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
    require(metap)
    all.meta <- list()
    i=1
    for(i in 1:ntest) {
        pv = P[[i]]
        qv = Q[[i]]
        fc = S[[i]]
        ##meta.p = apply(pv, 1, max, na.rm=TRUE ) ## maximum statistic
        ##meta.q = apply(qv, 1, max, na.rm=TRUE ) ## maximum statistic
        meta.p = apply(pv, 1, function(p) metap::allmetap(p, method="sumlog")$p[[1]])
        meta.q = p.adjust(meta.p, method="fdr")
        ss.rank <- function(x) scale(sign(x)*rank(abs(x)),center=FALSE)
        meta.fx = rowMeans( apply(S[[i]], 2, ss.rank), na.rm=TRUE)
        meta = data.frame(fx=meta.fx, p=meta.p, q=meta.q)
        ##all.meta[[i]] = data.frame( gene.set=names(gmt), meta=meta, fc=I(fc), p=I(pv), q=I(qv))
        all.meta[[i]] = data.frame( meta=meta, fc=I(fc), p=I(pv), q=I(qv))
    }
    names(all.meta) = tests

    ##--------------------------------------------------
    ## Add meta matrices (this becomes quite large...)
    ##--------------------------------------------------
    cat("computing meta-matrix... \n")
    m <- list(gsva=zx.gsva, ssgsea=zx.ssgsea, rnkcorr=zx.rnkcorr)        
    m = m[which(!sapply(m,is.null))]
    names(m)
    ## conform
    if(length(m)>1) {
        ##m <- lapply(m, function(x) apply(x,2,rank,na.last="keep"))
        m = lapply(m, function(x) scale(x,center=FALSE))
        avg.m  <- Reduce('+',m) / length(m)  ## matrix mean
        meta.matrix <- scale(avg.m,center=FALSE)
        meta.matrix <- normalizeQuantiles(meta.matrix)
    } else {
        meta.matrix <- m[[1]]
    }
    m[["meta"]] <- meta.matrix
    
    ##timings0 <- do.call(rbind, timings)
    timings <- as.matrix(timings)
    rownames(timings) <- timings[,1]
    timings0 <- apply(as.matrix(timings[,-1]),2,as.numeric)
    rownames(timings0) <- rownames(timings)
    timings0 <- apply( timings0, 2, function(x) tapply(x,rownames(timings0),sum))
    res = list( meta = all.meta, sig.counts = sig.counts,  outputs = all.results,
               matrices = m, timings = timings0)

    return(res)
}


##X=ngs$X;Y=ngs$Y;design=ngs$model.parameters$design;contr.matrix=ngs$model.parameters$contr.matrix
##mc.cores=4;batch.correct=TRUE;mc.cores=4
##trend=TRUE;gsetX=zx.gsva
##trend=TRUE;gsetX=zx.rnkcorr
gset.fitContrastsWithLIMMA <- function( gsetX, contr.matrix, design,
                                       trend=TRUE,conform.output=FALSE)
{

    if(!is.null(design)) {
        ##xfit = normalizeQuantiles(xfit)
        vfit <- lmFit(gsetX, design)
        vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
        efit <- eBayes(vfit, trend=trend, robust=TRUE)
        ##efit <- eBayes(vfit, trend=trend, robust=FALSE)

        tables <- list()
        i=1
        exp.matrix = (design %*% contr.matrix)
        for(i in 1:ncol(contr.matrix)) {
            ##coef = colnames(contr.matrix)[i]
            top = topTable(efit, coef=i, sort.by="none",number=Inf, adjust.method="BH")
            j1 = which( exp.matrix[,i] > 0 )
            j0 = which( exp.matrix[,i] < 0 )
            ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
            mean1 = rowMeans(gsetX[,j1,drop=FALSE], na.rm=TRUE)
            mean0 = rowMeans(gsetX[,j0,drop=FALSE], na.rm=TRUE)
            top = top[rownames(gsetX),]
            top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1)
            head(top,10)
            tables[[i]] = top
        }
        names(tables) <- colnames(contr.matrix)
    } else {

        ##trend=TRUE
        tables <- list()
        i=1
        for(i in 1:ncol(contr.matrix)) {
            design0 <- cbind( 1, contr.matrix[,i])
            colnames(design0) <- c( "ref", colnames(contr.matrix)[i])
            vfit <- lmFit(gsetX, design0)
            efit <- eBayes(vfit, trend=trend, robust=TRUE)
            top = topTable(efit, coef=2, sort.by="none",number=Inf, adjust.method="BH")
            j1 = which( contr.matrix[,i] > 0 )
            j0 = which( contr.matrix[,i] < 0 )
            ## if(!( length(cf)==6 || length(cf)==7)) stop("wrong coef format")
            mean1 = rowMeans(gsetX[,j1,drop=FALSE], na.rm=TRUE)
            mean0 = rowMeans(gsetX[,j0,drop=FALSE], na.rm=TRUE)
            top = top[rownames(gsetX),]
            top = cbind(top, "AveExpr0"=mean0, "AveExpr1" = mean1)
            head(top,10)
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
        head(pv)
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
    require(png, quietly=TRUE)
    require(grid, quietly=TRUE)
    require(gridExtra, quietly=TRUE)
    enplots = dir(gsea_dir, pattern="enplot_")
    enplots0 = dir(gsea_dir, pattern="enplot_", full.names=TRUE)
    kk = match(gsets, gsub("enplot_|_[0-9]*.png$","",enplots))
    imgs = lapply( enplots0[kk], function(p) rasterGrob(as.raster(readPNG(p)),
                                                        interpolate=TRUE))
    grid.arrange( grobs=imgs, ncol=5)
}


##======================================================================
##======================================================================
##======================================================================
