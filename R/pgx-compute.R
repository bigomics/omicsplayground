##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


source(file.path(RDIR,"ngs-functions.R"))

if(0) {
    is.logx=FALSE;do.cluster=TRUE;auto.scale=TRUE;only.chrom=TRUE
    max.genes=max.genesets=25000;lib.dir=FILES;progress=NULL;only.hugo=1;
    extra.methods=c("meta.go","deconv","infer","drugs","wordcloud")
    gx.methods=c("ttest.welch","trend.limma");gset.methods=c("fisher","gsva");
    only.hugo=TRUE;only.proteincoding=TRUE;rik.orf=FALSE
    batch.correct=TRUE

    counts <- ngs$counts
    samples <- ngs$samples
    contrasts <- ngs$model.parameters$contr.matrix

}

pgx.createPGX <- function(counts, samples, contrasts, X=NULL, ## genes,
                          is.logx=NULL, do.cluster=TRUE, batch.correct=TRUE,
                          auto.scale=TRUE, filter.genes=TRUE, prune.samples=FALSE,
                          only.chrom=TRUE, rik.orf=FALSE,
                          only.hugo=TRUE, convert.hugo=TRUE,
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
    counts <- as.matrix(counts)
    contrasts <- as.matrix(contrasts)
    contrasts[is.na(contrasts)] <- 0

    ## convert group-wise contrast to sample-wise
    if("group" %in% names(samples) && nrow(contrasts)!=nrow(samples)) {
        is.group.contrast <- all(rownames(contrasts) %in% samples$group)
        is.group.contrast
        if(is.group.contrast) {
            ## group
            message("[pgx.createPGX] converting group contrast to sample-wise contrasts...")
            contrasts.new <- contrasts[samples$group,,drop=FALSE]
            rownames(contrasts.new) <- rownames(samples)
            contrasts <- contrasts.new
        }
    }
    
    ## prune.samples=FALSE
    used.samples <- names(which(rowSums(contrasts!=0)>0))
    if(prune.samples && length(used.samples) < ncol(counts) ) {

        message("[pgx.createPGX] pruning not-used samples...")
        counts  <- counts[,used.samples,drop=FALSE]
        samples <- samples[used.samples,,drop=FALSE]
        contrasts <- contrasts[used.samples,,drop=FALSE] ## sample-based!!! 

        message(cat("[pgx.createPGX] dim(counts) = ",dim(counts)))
        message(cat("[pgx.createPGX] dim(samples) = ",dim(samples)))
        message(cat("[pgx.createPGX] dim(contrasts) = ",dim(contrasts)))
        
    }
    
    ##-------------------------------------------------------------------
    ## conform
    ##-------------------------------------------------------------------
    kk <- intersect(colnames(counts),rownames(samples))
    counts  <- counts[,kk,drop=FALSE]
    samples <- samples[kk,,drop=FALSE]
    samples <- type.convert(samples) ## automatic type conversion
    if(!is.null(X)) X <- X[,kk,drop=FALSE]
        
    ##-------------------------------------------------------------------
    ## check counts
    ##-------------------------------------------------------------------
    guess.log <- (min(counts,na.rm=TRUE) < 0 || max(counts,na.rm=TRUE) < 100)
    guess.log <- guess.log && is.null(X) && (is.null(is.logx) || is.logx==TRUE)
    if(is.null(is.logx)) is.logx <- guess.log
    if(is.logx) {
        cat("[pgx.createPGX] input assumed log-expression (logarithm)\n")
        cat("[pgx.createPGX] ...undo-ing logarithm\n")
        counts <- pmax(2**counts-1,0)  ## undo logarithm
    } else {
        cat("[pgx.createPGX] input assumed counts (not logarithm)\n")
    }
    
    ##-------------------------------------------------------------------
    ## global scaling (no need for CPM yet)
    ##-------------------------------------------------------------------
    counts_multiplier = 1
    totcounts = colSums(counts)
    if(auto.scale) {        

        ## decide if normalizing is necessary (WARNING changes total counts!!!)
        totratio <- log10(max(totcounts) / min(totcounts))
        totratio
        if(totratio > 6) {
            cat("[pgx.createPGX:autoscale] normalizing necessary!\n")
            meancounts <- exp(mean(log(totcounts)))
            meancounts
            counts <- t( t(counts) / totcounts) * meancounts
        }
        
        ## check if too big (more than billion reads)
        mean.counts <- mean(colSums(counts,na.rm=TRUE))
        mean.counts
        is.toobig <- log10(mean.counts) > 9
        is.toobig
        if(is.toobig) {
            ## scale to about 10 million reads
            ##progress$inc(0.01, detail = "scaling down counts")
            unit <- 10**(round(log10(mean.counts)) - 7)
            unit
            counts <- counts / unit
            counts_multiplier = unit
        }
        counts_multiplier
        cat("[pgx.createPGX:autoscale] count_multiplier= ",counts_multiplier,"\n")
    }

    if(0 && auto.scale) {
        ## auto-scale down billions of counts like sometimes for proteomics
        q10 <- quantile(counts[counts>0.25],probs=0.10)
        q10
        if(q10 > 100) {
            counts <- counts / q10
            counts_multiplier = q10            
        }
        cat("[pgx.createPGX:autoscale] count_multiplier= ",counts_multiplier,"\n")
    }

    ##-------------------------------------------------------------------
    ## convert probe-IDs to gene symbol (do not translate yet to HUGO)
    ##-------------------------------------------------------------------
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
    ##ngs$genes   = data.frame(genes)
    ngs$contrasts = as.matrix(contrasts)
    ngs$X <- X  ## input normalized log-expression (can be NULL)
    
    ngs$total_counts = totcounts
    ngs$counts_multiplier = counts_multiplier
    
    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    ## take only first gene as rowname, retain others as alias
    gene0 <- rownames(ngs$counts)
    gene1 <- sapply(gene0, function(s) strsplit(s,split="[;,\\|]")[[1]][1])
    if(convert.hugo) {
        cat("[pgx.createPGX] converting to HUGO symbols...\n")
        gene1 <- alias2hugo(gene1)  ## convert to latest HUGO
    } else {
        cat("[pgx.createPGX] skip conversion to HUGO symbols\n")
    }
    ndup <- sum(duplicated(gene1))
    ndup
    if(ndup>0) {        
        x1 = tapply(1:nrow(ngs$counts), gene1, function(i)
            colSums(ngs$counts[i,,drop=FALSE]))
        x1 <- do.call(rbind, x1)
        ngs$counts = x1
        remove(x1)
    }
    if(ndup>0 && !is.null(ngs$X)) {        
        x1 = tapply(1:nrow(ngs$X), gene1, function(i)
            log2(colSums(2**ngs$X[i,,drop=FALSE])) )
        x1 <- do.call(rbind, x1)
        ngs$X = x1
        remove(x1)
    }

    ##-------------------------------------------------------------------
    ## create gene annotation if not given (no HUGO conversion)
    ##-------------------------------------------------------------------
    cat("[pgx.createPGX] annotating genes...\n")
    ngs$genes <- ngs.getGeneAnnotation(genes=rownames(ngs$counts))  
    rownames(ngs$genes) <- rownames(ngs$counts)
    ngs$genes[is.na(ngs$genes)] <- ""
    
    ##-------------------------------------------------------------------
    ## Filter out not-expressed
    ##-------------------------------------------------------------------
    if(filter.genes) {
        cat("[pgx.createPGX] filtering out not-expressed genes...\n")
        keep <- (Matrix::rowMeans(ngs$counts > 0) > 0) ## at least in 
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,,drop=FALSE]
        if(!is.null(ngs$X)) ngs$X <- ngs$X[keep,]
    }
    
    ##-------------------------------------------------------------------
    ## Filter genes?
    ##-------------------------------------------------------------------
    is.mouse <- (mean(grepl("[a-z]",rownames(ngs$counts))) > 0.9)
    org = ifelse(is.mouse, "mouse", "human")
    org
    cat("[pgx.createPGX] detected organism: ",org,"\n")
    do.filter <- (only.hugo | only.chrom | only.proteincoding | !rik.orf)        
    if(do.filter && org == "mouse") {
        SYMBOL = unlist(as.list(org.Mm.egSYMBOL))        
        has.chrloc = is.official = not.rik = is.protcoding = TRUE
        if(only.hugo) is.official <- (ngs$genes$gene_name %in% SYMBOL)
        if(!rik.orf) not.rik <- !grepl("Rik",ngs$genes$gene_name) ## ???
        ##imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
        if(only.chrom) has.chrloc <- !is.na(ngs$genes$chr)
        if(only.proteincoding) is.protcoding <- ngs$genes$gene_biotype %in% c("protein_coding")
        keep <- (not.rik & is.official & has.chrloc & is.protcoding)
        table(keep)
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]
        if(!is.null(ngs$X)) ngs$X <- ngs$X[keep,]
    }
    if(do.filter && org == "human") {
        SYMBOL = unlist(as.list(org.Hs.egSYMBOL))
        has.chrloc = is.official = is.protcoding = not.orf = TRUE
        if(only.hugo) is.official <- (ngs$genes$gene_name %in% SYMBOL)
        if(!rik.orf) not.orf <- !grepl("ORF",ngs$genes$gene_name)
        ##imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
        if(only.chrom) has.chrloc <- !is.na(ngs$genes$chr)
        if(only.proteincoding) is.protcoding <- ngs$genes$gene_biotype %in% c("protein_coding")
        keep <- (not.orf & is.official & has.chrloc & is.protcoding)
        table(keep)
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]
        if(!is.null(ngs$X)) ngs$X <- ngs$X[keep,]
    }
    
    ##-------------------------------------------------------------------
    ## Do infer cell cycle/gender here (before any batchcorrection)
    ##-------------------------------------------------------------------
    ngs <- compute.cellcycle.gender(ngs)
    head(ngs$samples)
    
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
            message("[pgx.createPGX] batch correcting for parameter '",b,"'\n")
            batch <- ngs$samples$batch
            zz <- which(ngs$counts==0, arr.ind=TRUE)
            cX <- log2(1 + ngs$counts)
            bx <- ngs$sample[,b]

            cX <- limma::removeBatchEffect(cX, batch=bx) ## in log-space
            cX <- pmax(2**cX - 1,0)
            cX[zz] <- 0
            ngs$counts <- cX   ## batch corrected counts...

            if(!is.null(ngs$X)) {
                message("[pgx.createPGX] batch correcting for logX\n")
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
        cat("[pgx.createPGX] clustering samples...\n")        
        ##if(!is.null(progress)) progress$inc(0.01, detail = "clustering")
        perplexity=30
        perplexity <- max(1,min(30,round(ncol(ngs$counts)/4)))
        perplexity=NULL
        ## ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=perplexity)
        ngs <- pgx.clusterSamples2(ngs, dims=c(2,3), perplexity=perplexity,
                                   methods=c("pca","tsne","umap"))
        head(ngs$samples)
        table(ngs$samples$cluster)

        ## Extra clustering methods: PCA, t-SNE, UMAP.
        ##ngs <- pgx.clusterSamples2(ngs, dims=c(2,3))
        ##names(ngs$cluster)
        ##ngs$tsne2d <- ngs$cluster$pos[["tsne2d"]] ## old style
        ##ngs$tsne3d <- ngs$cluster$pos[["tsne3d"]] ## old style        
    }

    ##-------------------------------------------------------------------
    ## Add normalized log-expression
    ##-------------------------------------------------------------------
    if(is.null(ngs$X)) {
        cat("[pgx.createPGX] calculating log-expression matrix...\n")        
        ngs$X <- logCPM(ngs$counts, total=NULL) 
        dim(ngs$X)
    }
    if(!all(dim(ngs$X) == dim(ngs$counts))) {
        stop("[pgx.createPGX] dimensions of X and counts do not match\n")        
    }

    
    return(ngs)
}

if(0) {
    max.genes=19999;max.genesets=9999
    gx.methods = c("ttest.welch","trend.limma","edger.qlf")
    gset.methods = c("fisher","gsva","fgsea")
    extra.methods = c("meta.go","deconv","infer","drugs","wordcloud")
    extra.methods = c("wordcloud")    
}

pgx.computePGX <- function(ngs, 
                           max.genes = 19999, max.genesets = 9999, 
                           gx.methods = c("ttest.welch","trend.limma","edger.qlf"),
                           gset.methods = c("fisher","gsva","fgsea"),
                           do.cluster = TRUE,
                           extra.methods = c("meta.go","deconv","infer","drugs","wordcloud"),
                           lib.dir = "../lib", progress=NULL)
{
    
    ##======================================================================
    ##======================================================================
    ##======================================================================

    if(!"contrasts" %in% names(ngs)) {
        stop("[pgx.computePGX] FATAL:: no contrasts in object")
    }
    
    ## contrast matrix
    colnames(ngs$contrasts)
    contr.matrix <- ngs$contrasts
    contr.matrix <- as.matrix(contr.matrix) ## must be numeric matrix...

    ## create groups if not exists. DEPRECATED: maybe not necessary anymore (IK)
    if(FALSE  && !"group" %in% colnames(ngs$samples)) {
        ct = apply(ngs$samples,2,function(px) mean(px %in% rownames(ngs$contrasts),na.rm=TRUE))
        group.col = which(ct > 0.95)
        group.col
        if(length(group.col)>0) {
            grp.var = colnames(ngs$samples)[group.col[1]]
            cat(paste0("INFO [pgx-upload] assigning '",grp.var,"' variable as 'group' column\n"))
            ## colnames(ngs$samples)[group.col[1]] <- "group"
            ngs$samples$group <- ngs$samples[,group.col[1]] 
        } else {
            ## stop("sample annotation file must have 'group' column\n")
            cat("creating statistical groups...\n")
            ngs$samples$group <- pgx.getConditions(contr.matrix, nmax=1)             
        }
    }
    
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
    
    ngs <- compute.testGenes(
        ngs, contr.matrix,
        max.features = max.genes,
        test.methods = gx.methods)
    head(ngs$gx.meta$meta[[1]])        
       
    ## ------------------ gene set tests -----------------------
    if(!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

    max.features=max.genes;test.methods=gset.methods
    ngs <- compute.testGenesets(
        ngs, max.features = max.genesets,
        test.methods = gset.methods,
        lib.dir = lib.dir )
    head(ngs$gset.meta$meta[[1]])
    
    ## ------------------ extra analyses ---------------------
    if(!is.null(progress)) progress$inc(0.3, detail = "extra modules")
    ##extra <- c("meta.go","deconv","infer","drugs")
    ##extra <- c("meta.go","infer","drugs")
    ngs <- compute.extra(ngs, extra=extra.methods, lib.dir=lib.dir)
    
    ngs$timings
    return(ngs)
}


pgx.computeObjectPGX <- function(counts, samples, contrasts, ## genes, 
                                 ##gx.methods = c("trend.limma","edger.qlf","deseq2.wald"),
                                 max.genes = 9999, max.genesets = 9999, only.hugo=TRUE,
                                 gx.methods = c("ttest.welch","trend.limma","edger.qlf"),
                                 gset.methods = c("fisher","gsva","fgsea"),
                                 extra.methods = c("meta.go","deconv","infer","drugs","wordcloud"),
                                 lib.dir = "../lib", do.cluster=TRUE,
                                 progress=NULL)
{

    library(org.Hs.eg.db)
    cat("***************************************************************\n")    
    cat("****  DEPRECATED: please use createPGX() and computePGX() *****\n")
    cat("***************************************************************\n")
    
    if(!is.null(progress)) progress$inc(0.01, detail = "creating object")
    
    ##-------------------------------------------------------------------
    ## convert to gene symbol
    ##-------------------------------------------------------------------
    if(1) {
        symbol <- probe2symbol(rownames(counts), type=NULL)
        jj <- which(!is.na(symbol))
        counts <- as.matrix(counts[jj,])
        rownames(counts) <- symbol[jj]
    }

    ##-------------------------------------------------------------------
    ## create ngs object
    ##-------------------------------------------------------------------
    ##load(file=rda.file, verbose=1)
    ngs <- list()  ## empty object
    ngs$name = "data set"
    ngs$date = date()
    ngs$datatype = "unknown"
    ngs$description = "data set"

    ngs$samples = data.frame(samples, check.names=FALSE)
    ngs$counts  = as.matrix(counts)
    ##ngs$genes   = data.frame(genes)
    ngs$contrasts = contrasts
    
    ##-------------------------------------------------------------------
    ## check counts
    ##-------------------------------------------------------------------
    is.log <- (min(ngs$counts,na.rm=TRUE) < 0 || max(ngs$counts,na.rm=TRUE) < 100)
    is.log
    if(is.log) {
        ngs$counts <- 2**ngs$counts  ## undo logarithm
        cat("[pgx.computeObjectPGX] undo logarithm\n")
    }

    mean.counts <- mean(colSums(ngs$counts,na.rm=TRUE))
    is.toobig <- log10(mean.counts) > 10
    is.toobig
    ngs$counts_multiplier = 1
    if(is.toobig) {
        ## scale to about 10 million reads
        ##progress$inc(0.01, detail = "scaling down counts")
        unit <- 10**(round(log10(mean.counts)) - 7)  
        ngs$counts <- ngs$counts / unit
        ngs$counts_multiplier = unit
    }
    ngs$counts_multiplier
    cat("[pgx.computeObjectPGX] count_multiplier= ",ngs$counts_multiplier,"\n")
    
    ##-------------------------------------------------------------------
    ## create gene annotation if not given
    ##-------------------------------------------------------------------
    is.mouse <- (mean(grepl("[a-z]",rownames(counts))) > 0.9)
    org = ifelse(is.mouse, "mouse", "human")
    org
    gene.symbol <- NULL
    if(org == "human") {
        require(org.Hs.eg.db)
        GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
        gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
        gene.map <- sapply(as.list(org.Hs.egMAP),"[",1)
        names(GENE.TITLE) = gene.symbol
        names(gene.map) = gene.symbol

        ## get gene biotype
        require("EnsDb.Hsapiens.v86")
        daf <- transcripts(EnsDb.Hsapiens.v86,
                           columns = c("gene_name", "gene_biotype"),
                           return.type="DataFrame")
        GENE.BIOTYPE = daf$gene_biotype
        names(GENE.BIOTYPE) = daf$gene_name
        
    }
    if(org == "mouse") {
        require(org.Mm.eg.db)
        GENE.TITLE = unlist(as.list(org.Mm.egGENENAME))
        gene.symbol = unlist(as.list(org.Mm.egSYMBOL))
        gene.map <- sapply(as.list(org.Mm.egCHR),"[",1)
        names(GENE.TITLE) = gene.symbol
        names(gene.map) = gene.symbol

        ## get gene biotype
        require("EnsDb.Mmusculus.v79")
        daf <- transcripts(EnsDb.Mmusculus.v79,
                           columns = c("gene_name", "gene_biotype"),
                           return.type="DataFrame")
        GENE.BIOTYPE = daf$gene_biotype
        names(GENE.BIOTYPE) = daf$gene_name
    }
    if(!is.null(progress)) {
        aa <- paste("detected organism: ",org)
        progress$inc(0.01, detail = aa)
    }
    cat("[pgx.computeObjectPGX] detected organism: ",org,"\n")
    
    gene <- rownames(ngs$counts)
    gene1 <- sapply(gene, function(s) strsplit(s,split="[;,]")[[1]][1])
    gene1 <- alias2hugo(gene1)  ## convert to HUGO
    ngs$genes = data.frame( gene_name = gene1,
                           gene_alias = gene,
                           chr = gene.map[gene1],
                           gene_title = GENE.TITLE[gene1],
                           gene_biotype = GENE.BIOTYPE[gene1] )
    ##rownames(ngs$genes) <- gene1
    
    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    ndup <- sum(duplicated(ngs$genes$gene_name))
    ndup
    if(ndup>0) {        
        gene = as.character(ngs$genes$gene_name)        
        x1 = tapply(1:nrow(ngs$counts), gene, function(i) colSums(ngs$counts[i,,drop=FALSE]))
        x1 <- do.call(rbind, x1)
        ngs$genes = ngs$genes[match(rownames(x1), ngs$genes$gene_name),]
        ngs$counts = x1
        rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
        remove(x1)
    }

    ##-------------------------------------------------------------------
    ## Filter genes?
    ##-------------------------------------------------------------------
    if(org == "mouse") {
        has.name <- !is.na(ngs$genes$gene_name)
        is.hugo = TRUE
        if(only.hugo) is.hugo <- (ngs$genes$gene_name %in% gene.symbol)
        rik.genes <- grepl("Rik",ngs$genes$gene_name)
        ##imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
        has.chrloc <- !is.na(ngs$genes$chr)        
        keep <- (has.name & !rik.genes & is.hugo & has.chrloc)
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]        
    }
    if(org == "human") {
        has.name <- !is.na(ngs$genes$gene_name)
        is.hugo = TRUE
        if(only.hugo) is.hugo <- (ngs$genes$gene_name %in% gene.symbol)
        ##orf.genes <- grepl("ORF",ngs$genes$gene_name)
        ##imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
        has.chrloc <- !is.na(ngs$genes$chr)
        keep <- (has.name & is.hugo & has.chrloc)
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]        
    }
    
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    if(do.cluster) {
        if(!is.null(progress)) progress$inc(0.01, detail = "clustering")
        perplexity <- max(2,min(30,round(ncol(ngs$counts)/4)))
        perplexity
        ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=perplexity)
        head(ngs$samples)
        table(ngs$samples$cluster)
    }
    
    ##======================================================================
    ##======================================================================
    ##======================================================================

    ## contrast matrix
    contr.matrix <- ngs$contrasts

    ## create groups if not exists
    if(!"group" %in% colnames(ngs$samples)) {
        ct = apply(ngs$samples,2,function(px) mean(px %in% rownames(ngs$contrasts),na.rm=TRUE))
        group.col = which(ct > 0.95)
        if(length(group.col)>0) {
            grp.var = colnames(ngs$samples)[group.col[1]]
            cat(paste0("INFO [pgx-upload] assigning '",grp.var,"' variable as 'group' column\n"))
            ## colnames(ngs$samples)[group.col[1]] <- "group"
            ngs$samples$group <- ngs$samples[,group.col[1]] 
        } else {
            stop("sample annotation file must have 'group' column\n")
        }
    }
    
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
    
    ngs <- compute.testGenes(
        ngs, contr.matrix,
        max.features = max.genes,
        test.methods = gx.methods)
    head(ngs$gx.meta$meta[[1]])        
    
    ## ------------------ gene set tests -----------------------
    if(!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

    max.features=max.genes;test.methods=gset.methods
    system.time(
        ngs <- compute.testGenesets(
            ngs, max.features = max.genesets,
            test.methods = gset.methods,
            lib.dir = lib.dir )
    )
    head(ngs$gset.meta$meta[[1]])
    
    ## ------------------ extra analyses ---------------------
    if(!is.null(progress)) progress$inc(0.3, detail = "extra modules")
    ##extra <- c("meta.go","deconv","infer","drugs")
    ##extra <- c("meta.go","infer","drugs")
    ngs <- compute.extra(ngs, extra=extra.methods, lib.dir=lib.dir)
    
    ngs$timings
    return(ngs)
}



##=====================================================================================
##========================== END OF FILE ==============================================
##=====================================================================================






