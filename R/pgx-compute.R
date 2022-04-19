##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


if(0) {
    X=NULL;is.logx=NULL;do.cluster=TRUE;auto.scale=TRUE;only.known=TRUE
    max.genes=max.genesets=25000;lib.dir=FILES;progress=NULL;only.hugo=1;
    extra.methods=c("meta.go","deconv","infer","drugs","wordcloud")
    gx.methods=c("trend.limma");gset.methods=c("fisher","gsva");
    only.hugo=TRUE;only.proteincoding=TRUE;rik.orf=FALSE;prune.samples=FALSE    
    do.clustergenes=cluster.contrasts=batch.correct=TRUE
    
    counts <- ngs$counts
    samples <- ngs$samples
    contrasts <- ngs$contrasts
}

pgx.createPGX <- function(counts, samples, contrasts, X=NULL, ## genes,
                          is.logx=NULL, batch.correct=TRUE,
                          auto.scale=TRUE, filter.genes=TRUE, prune.samples=FALSE,
                          only.known=TRUE, only.hugo=TRUE, convert.hugo=TRUE,
                          do.cluster=TRUE,  cluster.contrasts = FALSE, do.clustergenes=TRUE,
                          only.proteincoding=TRUE)
{
    
    ##if(!is.null(progress)) progress$inc(0.01, detail = "creating PGX object")
    if(0 && !"group" %in% colnames(samples)) {
        stop("samples information must have 'group' column\n")
        return(NULL)
    }

    if(!is.null(X) && !all(dim(counts)==dim(X))) {
        stop("dimension of counts and X do not match\n")
    }
    
    ##-------------------------------------------------------------------
    ## clean up input files
    ##-------------------------------------------------------------------
    samples <- data.frame(samples)
    counts <- as.matrix(counts)
    if(is.null(contrasts))  contrasts <- samples[,0]
    ## contrasts[is.na(contrasts)] <- 0
    
    ## contrast matrix
    colnames(contrasts)
    is.numbered <- all(unique(as.vector(contrasts)) %in% c(-1,0,1))
    is.numbered <- all(sapply(type.convert(data.frame(contrasts)),class) %in% c("numeric","integer"))
    is.numbered
    ct.type <- c("labeled (new style)","numbered (old style)")[1 + 1*is.numbered]
    message("[createPGX] contrast type :",ct.type)
    if(is.numbered) {
        message("[createPGX] converting numbered contrasts to LABELED")
        ##contrasts <- makeContrastsFromLabelMatrix(contrasts)
        contrasts <- contrastAsLabels(contrasts) 
    }

    ## convert group-wise contrast to sample-wise
    is.group1=is.group2=FALSE
    if("group" %in% colnames(samples)) is.group1 <- all(rownames(contrasts) %in% samples$group)
    if("condition" %in% colnames(samples)) is.group2 <- all(rownames(contrasts) %in% samples$condition)
    is.group.contrast <- (is.group1 || is.group2)
    if(is.group.contrast && nrow(contrasts)<nrow(samples)) {
        ## group
        message("[createPGX] converting group contrast to sample-wise contrasts...")
        if(is.group1) grp <- as.character(samples$group)
        if(is.group2) grp <- as.character(samples$condition)            
        contrasts.new <- contrasts[grp,,drop=FALSE]
        rownames(contrasts.new) <- rownames(samples)
        contrasts <- contrasts.new
    }
    
    ## sanity check...
    if( !all(rownames(contrasts)==rownames(samples)) &&
        !all(rownames(contrasts)==colnames(counts)) ) {
        stop("[createPGX] FATAL :: matrices do not match")
    }

    ## prune.samples=FALSE    
    ##used.samples <- names(which(rowSums(contrasts!=0)>0))
    contrasts[contrasts==""] <- NA
    used.samples <- names(which(rowSums(!is.na(contrasts))>0))        
    if(prune.samples && length(used.samples) < ncol(counts) ) {

        message("[createPGX] pruning unused samples...")
        counts    <- counts[,used.samples,drop=FALSE]
        samples   <- samples[used.samples,,drop=FALSE]
        contrasts <- contrasts[used.samples,,drop=FALSE] ## sample-based!!! 

        message("[createPGX] dim(counts) = ",dim(counts))
        message("[createPGX] dim(samples) = ",dim(samples))
        message("[createPGX] dim(contrasts) = ",dim(contrasts))
        
    }

    ##-------------------------------------------------------------------
    ## conform
    ##-------------------------------------------------------------------
    message("[createPGX] conforming matrices...")
    kk <- intersect(colnames(counts),rownames(samples))
    counts  <- counts[,kk,drop=FALSE]
    samples <- samples[kk,,drop=FALSE]
    samples <- type.convert(samples) ## automatic type conversion
    if(!is.null(X)) X <- X[,kk,drop=FALSE]
    if(all(kk %in% rownames(contrasts))) {
        contrasts <- contrasts[kk,,drop=FALSE]
    }

    ##-------------------------------------------------------------------
    ## check counts
    ##-------------------------------------------------------------------
    message("[createPGX] check logarithm/linear...")
    guess.log <- (min(counts,na.rm=TRUE) < 0 || max(counts,na.rm=TRUE) < 100)
    guess.log <- guess.log && is.null(X) && (is.null(is.logx) || is.logx==TRUE)
    guess.log
    if(is.null(is.logx)) is.logx <- guess.log
    is.logx
    if(is.logx) {
        cat("[createPGX] input assumed log-expression (logarithm)\n")
        cat("[createPGX] ...undo-ing logarithm\n")
        counts <- pmax(2**counts-1,0)  ## undo logarithm
    } else {
        cat("[createPGX] input assumed counts (not logarithm)\n")
    }
    
    ##-------------------------------------------------------------------
    ## How to deal with missing values??
    ##-------------------------------------------------------------------
    if( any(is.na(counts)) || any(is.infinite(counts))) {
        message("[createPGX] setting missing values to zero")
        counts[is.na(counts)|is.infinite(counts)] <- 0
    }

    ##-------------------------------------------------------------------
    ## global scaling (no need for CPM yet)
    ##-------------------------------------------------------------------
    message("[createPGX] scaling counts...")
    counts_multiplier = 1
    totcounts = Matrix::colSums(counts, na.rm=TRUE)
    totcounts
    if(auto.scale) {        

        ## If the difference in total counts is too large, we need to
        ## euqalize them because the thresholds can become
        ## strange. Here we decide if normalizing is necessary (WARNING
        ## changes total counts!!!)
        totratio <- log10(max(totcounts,na.rm=TRUE) / min(totcounts,na.rm=TRUE))
        totratio
        if(totratio > 6) {
            cat("[createPGX:autoscale] WARNING: too large differences in total counts. forcing normalization.")
            meancounts <- exp(mean(log(totcounts)))
            meancounts
            counts <- t( t(counts) / totcounts) * meancounts
        }
        
        ## check if too big (more than billion reads)
        mean.counts <- mean(Matrix::colSums(counts,na.rm=TRUE))
        mean.counts
        is.toobig <- log10(mean.counts) > 9
        is.toobig
        if(is.toobig) {
            ## scale to about 10 million reads
            ##progress$inc(0.01, detail = "scaling down counts")
            cat("[createPGX:autoscale] WARNING: too large total counts. Scaling down to 10e6 reads.\n")
            unit <- 10**(round(log10(mean.counts)) - 7)
            unit
            counts <- counts / unit
            counts_multiplier = unit
        }
        counts_multiplier
        cat("[createPGX:autoscale] count_multiplier= ",counts_multiplier,"\n")
    }

    if(0 && auto.scale) {
        ## auto-scale down billions of counts like sometimes for proteomics
        q10 <- quantile(counts[counts>0.25],probs=0.10)
        q10
        if(q10 > 100) {
            counts <- counts / q10
            counts_multiplier = q10            
        }
        cat("[createPGX:autoscale] count_multiplier= ",counts_multiplier,"\n")
    }

    ##-------------------------------------------------------------------
    ## convert probe-IDs to gene symbol (do not translate yet to HUGO)
    ##-------------------------------------------------------------------
    message("[createPGX] converting probes to symbol...")
    symbol <- probe2symbol(rownames(counts), type=NULL)  ## auto-convert function
    if(mean(rownames(counts) == symbol,na.rm=TRUE) < 0.5) {  ## why??
        jj <- which(!is.na(symbol))
        counts <- as.matrix(counts[jj,])
        rownames(counts) <- symbol[jj]
        if(!is.null(X)) {
            rownames(X) <- rownames(counts)
        }
    }
    dim(counts)
    
    ##-------------------------------------------------------------------
    ## create ngs object
    ##-------------------------------------------------------------------
    message("[createPGX] creating pgx object...")

    ##load(file=rda.file, verbose=1)
    ngs <- list()  ## empty object
    ngs$name = "data set"
    this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ##ngs$date = date()
    ngs$date = this.date
    ngs$datatype = "unknown"
    ngs$description = "data set"

    ngs$samples = data.frame(samples, check.names=FALSE)
    ngs$counts  = as.matrix(counts)
    ngs$contrasts = contrasts
    ngs$X <- X  ## input normalized log-expression (can be NULL)
    
    ngs$total_counts = totcounts
    ngs$counts_multiplier = counts_multiplier
    
    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    ## take only first gene as rowname, retain others as alias
    gene0 <- rownames(ngs$counts)
    gene1 <- gene0
    gene1 <- sapply(gene0, function(s) strsplit(s,split="[;,\\|]")[[1]][1])

    if(convert.hugo) {
        message("[createPGX] converting to HUGO symbols...")
        gene1 <- alias2hugo(gene1)  ## convert to latest HUGO
    } else {
        message("[createPGX] skip conversion to HUGO symbols")
    }
    ndup <- sum(duplicated(gene1))
    ndup
    if(ndup>0) {
        message("[createPGX:autoscale] duplicated rownames detected: summing up rows (counts).")
        x1 = tapply(1:nrow(ngs$counts), gene1, function(i)
            Matrix::colSums(ngs$counts[i,,drop=FALSE]))
        x1 <- do.call(rbind, x1)
        ngs$counts = x1
        remove(x1)
    }
    if(ndup>0 && !is.null(ngs$X)) {
        ## cat("[createPGX:autoscale] duplicated rownames detected: collapsing rows (X).\n")
        x1 = tapply(1:nrow(ngs$X), gene1, function(i)
            log2(Matrix::colSums(2**ngs$X[i,,drop=FALSE])) )
        x1 <- do.call(rbind, x1)
        ngs$X = x1
        remove(x1)
    }

    ##-------------------------------------------------------------------
    ## create gene annotation if not given (no HUGO conversion)
    ##-------------------------------------------------------------------
    message("[createPGX] annotating genes...")
    ngs$genes <- ngs.getGeneAnnotation(genes=rownames(ngs$counts))  
    rownames(ngs$genes) <- rownames(ngs$counts)
    ngs$genes[is.na(ngs$genes)] <- ""
    
    ##-------------------------------------------------------------------
    ## Filter out not-expressed
    ##-------------------------------------------------------------------
    if(filter.genes) {
        ## There is second filter in the statistics computation. This
        ## first filter is primarily to reduce the counts table.
        message("[createPGX] filtering out not-expressed genes...")
        keep <- (Matrix::rowMeans(ngs$counts > 0) > 0) ## at least in one...
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,,drop=FALSE]
        if(!is.null(ngs$X)) {
            ngs$X <- ngs$X[keep,]
        }
    }
    
    ##-------------------------------------------------------------------
    ## Filter genes?
    ##-------------------------------------------------------------------
    is.mouse <- (mean(grepl("[a-z]",rownames(ngs$counts))) > 0.9)
    org = ifelse(is.mouse, "mouse", "human")
    org
    message("[createPGX] detected organism: ",org,"")
    do.filter <- (only.hugo | only.known | only.proteincoding )        
    if(do.filter && org == "mouse") {
        message("[createPGX] filtering genes...")
        require(org.Mm.eg.db)
        SYMBOL = unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))        
        is.hugo = is.known = is.protcoding = TRUE
        if(only.hugo) is.hugo <- (ngs$genes$gene_name %in% SYMBOL)
        ##imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
        if(only.known) {
            is.known <- !grepl("Rik|^Orf|^Loc",ngs$genes$gene_name) ## ???
        }
        if(only.proteincoding) {
            is.protcoding <- ngs$genes$gene_biotype %in% c("protein_coding")
        }
        keep <- (is.known & is.hugo & is.protcoding)
        table(keep)
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]
        if(!is.null(ngs$X)) ngs$X <- ngs$X[keep,]
    }
    if(do.filter && org == "human") {
        message("[createPGX] filtering genes...")        
        require(org.Hs.eg.db)
        SYMBOL = unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
        is.hugo = is.protcoding = is.known = TRUE
        if(only.hugo) is.hugo <- (ngs$genes$gene_name %in% SYMBOL)
        if(only.known) {
            is.known <- !grepl("^ORF|^LOC",ngs$genes$gene_name) ## ???
        }
        if(only.proteincoding) {
            is.protcoding <- ngs$genes$gene_biotype %in% c("protein_coding")
        }
        keep <- (is.known & is.hugo & is.protcoding)
        table(keep)
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]
        if(!is.null(ngs$X)) ngs$X <- ngs$X[keep,]
    }
    
    ##-------------------------------------------------------------------
    ## Infer cell cycle/gender here (before any batchcorrection)
    ##-------------------------------------------------------------------
    ngs <- compute.cellcycle.gender(ngs)
    Matrix::head(ngs$samples)
    
    ##-------------------------------------------------------------------
    ## Batch-correction (if requested. WARNING: changes counts )
    ##-------------------------------------------------------------------
    batch.par <- c("batch","batch2")
    has.batchpar <- any(batch.par %in% colnames(ngs$samples))
    if(batch.correct && has.batchpar) {
        b="batch"
        bb <- intersect(colnames(ngs$samples),batch.par)
        bb
        for(b in bb) {
            message("[createPGX] batch correcting for parameter '",b,"'\n")
            batch <- ngs$samples$batch
            zz <- which(ngs$counts==0, arr.ind=TRUE)
            cX <- log2(1 + ngs$counts)
            bx <- ngs$sample[,b]

            cX <- limma::removeBatchEffect(cX, batch=bx) ## in log-space
            cX <- pmax(2**cX - 1,0)
            cX[zz] <- 0
            ngs$counts <- cX   ## batch corrected counts...

            if(!is.null(ngs$X)) {
                message("[createPGX] batch correcting for logX\n")
                ngs$X <- limma::removeBatchEffect(ngs$X, batch=bx) ## in log-space
                ngs$X[zz] <- 0
            }
            
        }
        remove(cX)
    }

    
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    if(do.cluster) {
        message("[createPGX] clustering samples...")        
        ## if(!is.null(progress)) progress$inc(0.01, detail = "clustering")
        ## ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=NULL)
        ngs <- pgx.clusterSamples2(
            ngs, dims = c(2,3),
            ## replace.orig = FALSE,
            perplexity = NULL,
            methods = c("pca","tsne","umap")
        )

        ## NEED RETHINK: for the moment we use combination of t-SNE/UMAP
        ## posx <- scale(do.call(cbind,ngs$cluster$pos))
        posx <- scale(cbind(ngs$cluster$pos[["umap2d"]],ngs$cluster$pos[["tsne2d"]]))
        ##posx <- scale(cbind(ngs$cluster$pos[["umap3d"]],ngs$cluster$pos[["tsne3d"]]))
        idx <- pgx.findLouvainClusters(posx, level=1, prefix='c', small.zero=0.0)
        table(idx)
        if(length(unique(idx))==1) {
            ## try again with finer settings if single cluster...
            idx <- pgx.findLouvainClusters(posx, level=2, prefix='c', small.zero=0.01)
        }
        ngs$samples$cluster <- idx        
        Matrix::head(ngs$samples)
        table(ngs$samples$cluster)
    }    

    if(cluster.contrasts) {
        ## Add cluster contrasts
        message("[createPGX] adding cluster contrasts...")
        Y = ngs$samples[,"cluster",drop=FALSE]
        if(length(unique(Y)) < 2) {
            message("[createPGX] warning: only one cluster.")
        } else {
            ct <- makeDirectContrasts(Y, ref="others")
            ctx <- contrastAsLabels(ct$exp.matrix)
            if(ncol(ngs$contrasts)==0) {
                ngs$contrasts <- ctx
            } else {
                ngs$contrasts <- cbind(ngs$contrasts, ctx)
            }
        }
    }

    ##-------------------------------------------------------------------
    ## Add normalized log-expression
    ##-------------------------------------------------------------------
    if(is.null(ngs$X)) {
        message("[createPGX] calculating log-expression matrix X...")        
        ##ngs$X <- logCPM(ngs$counts, total=NULL)
        ngs$X <- logCPM(ngs$counts, total=1e6, prior=1) 
        ## ngs$X <- limma::normalizeQuantiles(ngs$X)  ## Sure ???
        dim(ngs$X)
    } else {
        message("[createPGX] using passed log-expression X...")        
    }

    if(!all(dim(ngs$X) == dim(ngs$counts))) {
        stop("[createPGX] dimensions of X and counts do not match\n")        
    }
    
    if(do.clustergenes) {
        message("[createPGX] clustering genes...")
        ngs <- pgx.clusterGenes(ngs, methods='umap', dims=c(2,3), level='gene')
        ##ngs <- pgx.clusterGenes(ngs, methods='umap', dims=c(2,3), level='geneset')  ## gsetX not ready!!
    }
    
    return(ngs)
}

if(0) {
    max.genes=19999;max.genesets=9999
    gx.methods = c("trend.limma")
    gset.methods = c("fisher")
    extra.methods = c()

    gx.methods = c("ttest.welch","trend.limma","edger.qlf")
    gset.methods = c("fisher","gsva","fgsea")
    extra.methods = c("meta.go","deconv","infer","drugs","wordcloud")

    do.cluster=TRUE;use.design=TRUE;prune.samples=FALSE
}

.EXTRA.METHODS = c("meta.go","deconv","infer","drugs","wordcloud")

pgx.computePGX <- function(ngs, 
                           max.genes = 19999, max.genesets = 9999, 
                           gx.methods = c("ttest.welch","trend.limma","edger.qlf"),
                           gset.methods = c("fisher","gsva","fgsea"),
                           do.cluster = TRUE, use.design = TRUE, prune.samples = FALSE,
                           extra.methods = .EXTRA.METHODS,
                           lib.dir = "../lib", progress=NULL)
{
    
    ##======================================================================
    ##======================================================================
    ##======================================================================

    if(!"contrasts" %in% names(ngs)) {
        stop("[pgx.computePGX] FATAL:: no contrasts in object")
    }
    message("[pgx.computePGX] called.")
    
    ## make proper contrast matrix
    contr.matrix <- ngs$contrasts
    contr.values <- unique(as.vector(contr.matrix))
    is.numcontrast <- all(contr.values %in% c(NA,-1,0,1))    
    is.numcontrast <- is.numcontrast && (-1 %in% contr.values) && (1 %in% contr.values)
    is.numcontrast
    if(!is.numcontrast) {
        message("[pgx.computePGX] converting label to numeric contrast...")        
        contr.matrix <- makeContrastsFromLabelMatrix(contr.matrix)
        contr.matrix <- sign(contr.matrix) ## sign is fine
    }

    ## select valid contrasts
    sel <- Matrix::colSums(contr.matrix == -1)>0 & Matrix::colSums(contr.matrix == 1)>0
    contr.matrix <- contr.matrix[,sel,drop=FALSE]
    
    ##======================================================================
    ##======================================================================
    ##======================================================================

    ngs$timings <- c()    
    GENETEST.METHODS = c("ttest","ttest.welch","ttest.rank",
                         "voom.limma","trend.limma","notrend.limma",
                         "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
    GENESETTEST.METHODS = c("fisher","gsva","ssgsea","spearman",
                            "camera", "fry","fgsea") ## no GSEA, too slow...
    
    ## ------------------ gene level tests ---------------------
    if(!is.null(progress)) progress$inc(0.1, detail = "testing genes")
    message("[pgx.computePGX] testing genes...")

    ngs <- compute.testGenes(
        ngs, contr.matrix,
        max.features = max.genes,
        test.methods = gx.methods,
        use.design = use.design,
        prune.samples = prune.samples
    )
    Matrix::head(ngs$gx.meta$meta[[1]])        
       
    ## ------------------ gene set tests -----------------------
    if(!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

    message("[pgx.computePGX] testing genesets...")    
    ##max.features=max.genes;test.methods=gset.methods
    ngs <- compute.testGenesets(
        ngs, max.features = max.genesets,
        test.methods = gset.methods,
        lib.dir = lib.dir )
    Matrix::head(ngs$gset.meta$meta[[1]])

    
    if(do.cluster) {
        message("[pgx.computePGX] clustering genes...")
        ##ngs <- pgx.clusterGenes(ngs, methods='umap', dims=c(2,3), level='gene') 
        ngs <- pgx.clusterGenes(ngs, methods='umap', dims=c(2,3), level='geneset')  ## gsetX not ready!!
    }

    
    ## ------------------ extra analyses ---------------------
    if(!is.null(progress)) progress$inc(0.3, detail = "extra modules")
    message("[pgx.computePGX] computing extra modules...")

    ##extra <- c("meta.go","deconv","infer","drugs")
    ##extra <- c("meta.go","infer","drugs")
    ngs <- compute.extra(ngs, extra=extra.methods, lib.dir=lib.dir)

    message("[pgx.computePGX] done!")    
    ngs$timings
    return(ngs)
}



##=====================================================================================
##========================== END OF FILE ==============================================
##=====================================================================================






