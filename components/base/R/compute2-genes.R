##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Workflow for differential expression analysis (at gene level)
##
## Input:  contr.matrix to be defined using group as levels
##
##
##
##
##

##SAVE.PARAMS <- ls()
if(0) {
    max.features=20000;type="counts"
    test.methods=c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
    test.methods=c("trend.limma","edger.qlf","deseq2.wald")
}

compute.testGenes <- function(pgx, contr.matrix, max.features=1000, 
                              test.methods=c("trend.limma","deseq2.wald","edger.qlf"),
                              use.design = TRUE, prune.samples=FALSE,
                              remove.outputs=TRUE )
{
    single.omics <- mean(grepl("\\[",rownames(pgx$counts))) < 0.1
    single.omics
    data.types <- unique(gsub("\\[|\\].*","",rownames(pgx$counts)))
    ##data.types
    if(single.omics || length(data.types)==1) {
        ## single-omics, no missing values
        cat(">>> computing gene tests for SINGLE-OMICS\n")
        pgx <- compute.testGenesSingleOmics(
            pgx = pgx, 
            contr.matrix = contr.matrix,
            max.features = max.features,
            test.methods = test.methods,
            use.design = use.design,
            prune.samples = prune.samples,
            remove.outputs = remove.outputs
        )
    } else {
        ## multi-omics, missing values allowed
        cat(">>> computing gene tests for MULTI-OMICS\n")
        pgx <- compute.testGenesMultiOmics(
            pgx=pgx,  ## type is inferred
            contr.matrix = contr.matrix,
            max.features = max.features,
            test.methods = test.methods,
            use.design = use.design,            
            remove.outputs = remove.outputs
        )
    }
    return(pgx)
}

if(0) {
    test.methods=c("trend.limma")
    test.methods=c("ttest.welch","trend.limma","edger.qlf")
    max.features=25000;type="counts";filter.low=TRUE
    use.design = TRUE; prune.samples=FALSE;
    contr.matrix = pgx$contrast
}

##contr.matrix=pgx$contrasts
compute.testGenesSingleOmics <- function(pgx, contr.matrix, max.features=1000,
                                         filter.low = TRUE, remove.outputs=TRUE,
                                         use.design = TRUE, prune.samples=FALSE,
                                         test.methods = c("trend.limma","deseq2.wald","edger.qlf"))
{

    message("[compute.testGenesSingleOmics] called")
    contr.matrix0 <- contr.matrix  ## SAVE
    
    ##-----------------------------------------------------------------------------
    ## Check parameters, decide group level
    ##-----------------------------------------------------------------------------    
    if(!("counts" %in% names(pgx))) {
        stop("[compute.testGenesSingleOmics] FATAL: cannot find counts in pgx object")
    }
    if(!("X" %in% names(pgx))) {
        stop("[compute.testGenesSingleOmics] FATAL: cannot find normalized expression X in pgx object")
    }
    
    ##is.expmatrix <- all(rownames(contr.matrix)==rownames(pgx$samples))
    is.expmatrix <- all(rownames(contr.matrix) %in% rownames(pgx$samples))
    is.expmatrix
    if(!is.expmatrix) {
        stop("[compute.testGenesSingleOmics] FATAL: contrast must be sample-wise")
    }

    stat.group = NULL
    if(use.design) {
        message("[compute.testGenesSingleOmics] detecting stat groups...")        
        stat.group <- pgx.getConditions(contr.matrix, nmax=0)  ## !!!
        names(stat.group) <- rownames(contr.matrix)
        nlev <- length(unique(stat.group))
        nlev        
        if (nlev >= nrow(contr.matrix)) {
            message("[compute.testGenesSingleOmics] cannot use groups, switching to no design")
            use.design = FALSE
        }
    }
    
    if(use.design) {

        message("[compute.testGenesSingleOmics] contrasts on groups (use design)")
        ##stat.group = rownames(pgx$samples)
        ## convert sample-wise contrasts to group-wise contrasts
        message("replacing contrast matrix...")
        stat0 <- sort(unique(stat.group))
        contr.matrix <- contr.matrix[match(stat0,stat.group),,drop=FALSE]
        rownames(contr.matrix) <- stat0

    } else if(!use.design) {

        message("[compute.testGenesSingleOmics] contrasts on samples (no design)")
        stat.group = rownames(contr.matrix)
        names(stat.group) <- rownames(contr.matrix)
    }

    ## table(stat.group)
    dim(contr.matrix)
    
    if(1) {
        message("[compute.testGenesSingleOmics] pruning unused contrasts")
        ## take out any empty comparisons
        sel <- which(Matrix::colSums(contr.matrix>0) & Matrix::colSums(contr.matrix<0))
        contr.matrix <- contr.matrix[,sel,drop=FALSE]
        contr.matrix[is.na(contr.matrix)] <- 0
    }
    dim(contr.matrix)
    
    ##-----------------------------------------------------------------------------
    ## normalize contrast matrix to zero mean and signed sums to one
    ##-----------------------------------------------------------------------------
    ## normalize?? why??
    message("[compute.testGenesSingleOmics] normalizing contrasts")
    for(i in 1:ncol(contr.matrix)) {
        m <- contr.matrix[,i]
        m[is.na(m)] <- 0
        contr.matrix[,i] <- 1*(m>0)/sum(m>0) - 1*(m<0)/sum(m<0)
    }
    dim(contr.matrix)
    
    ##-----------------------------------------------------------------------------
    ## create design matrix from defined contrasts (group or clusters)
    ##-----------------------------------------------------------------------------
    
    no.design <- all(stat.group %in% rownames(pgx$samples))  ## sample-wise design
    ## no.design <- is.expmatrix && !use.design
    design=NULL
    no.design
    
    if(no.design || !use.design) {
        
        message("[compute.testGenesSingleOmics] 6 : no design matrix ")
        ## SAMPLE-WISE DESIGN
        design = NULL
        exp.matrix <- contr.matrix

    } else {

        message("[compute.testGenesSingleOmics] 6 : creating model design matrix ")

        ## GROUP DESIGN
        ##stat.group[is.na(stat.group)] <- "_"
        notk <- which(!stat.group %in% rownames(contr.matrix))        
        length(notk)
        if(length(notk)) {
            stat.group[notk] <- "_"
        }
        design <- model.matrix(~ 0 + stat.group )  ## clean design no batch effects...
        colnames(design) <- sub("^stat.group", "", colnames(design))
        ## rownames(design) <- rownames(pgx$samples)
        if(is.null(names(stat.group))) {
            stop("[compute.testGenesSingleOmics] FATAL:: stat.group must have names")
        }
        rownames(design) <- names(stat.group)
        Matrix::head(design)
        
        ## make sure matrix align and compute experiment matrix
        design <- design[,match(rownames(contr.matrix),colnames(design)),drop=FALSE]
        colnames(design) <- rownames(contr.matrix)
        ##design = design[,rownames(contr.matrix),drop=FALSE]
        exp.matrix = (design %*% contr.matrix)

        ## check contrasts for sample sizes (at least 2 in each group) and
        ## remove otherwise        
        keep <- rep(TRUE,ncol(contr.matrix))
        keep = (Matrix::colSums(exp.matrix > 0) >= 1 & Matrix::colSums(exp.matrix < 0) >= 1)
        ##keep = ( Matrix::colSums(exp.matrix > 0) >= 2 & Matrix::colSums(exp.matrix < 0) >= 2 )
        table(keep)
        contr.matrix <- contr.matrix[,keep,drop=FALSE]
        exp.matrix   <- exp.matrix[,keep,drop=FALSE]
    }
    
    model.parameters <- list(
        design = design,
        contr.matrix = contr.matrix, 
        exp.matrix = exp.matrix,
        group = stat.group
    )
    pgx$model.parameters <- model.parameters
    
    ##-----------------------------------------------------------------------------
    ## Filter genes
    ##-----------------------------------------------------------------------------    
    if(is.null(names(stat.group))) {
        stop("[compute.testGenesSingleOmics] FATAL2:: stat.group must have names")
    }
    
    ## notice original counts will not be affected
    ss <- names(stat.group)
    gg <- rownames(pgx$counts)
    if(!is.null(pgx$X)) gg <- intersect(gg,rownames(pgx$X))
    counts  = pgx$counts[gg,ss,drop=FALSE]  
    genes   = pgx$genes[gg,]
    samples = pgx$samples[ss,]
    
    ## Rescale if too low. Often EdgeR/DeSeq can give errors of total counts
    ## are too low. Happens often with single-cell (10x?). We rescale
    ## to a minimum of 1 million counts (CPM)
    if(1) {
        mean.counts <- mean(Matrix::colSums(counts,na.rm=TRUE))
        mean.counts
        if( mean.counts < 1e6) {
            cat("[compute.testGenesSingleOmics] WARNING:: low total counts = ",mean.counts,"\n")        
            cat("[compute.testGenesSingleOmics] applying global mean scaling to 1e6...\n")        
            counts = counts * 1e6 / mean.counts
        }
        mean(Matrix::colSums(counts,na.rm=TRUE))
    }
    
    ## prefiltering for low-expressed genes (recommended for edgeR and
    ## DEseq2). Require at least in 2 or 1% of total. Specify the
    ## PRIOR CPM amount to regularize the counts and filter genes
    PRIOR.CPM = 1
    if(filter.low) {
        PRIOR.CPM = 0.25
        PRIOR.CPM = 1
        PRIOR.CPM
        AT.LEAST = ceiling(pmax(2,0.01*ncol(counts)))    
        cat("filtering for low-expressed genes: >",PRIOR.CPM,"CPM in >=",AT.LEAST,"samples\n")
        keep <- (rowSums( edgeR::cpm(counts) > PRIOR.CPM, na.rm=TRUE) >= AT.LEAST)
        ##keep <- edgeR::filterByExpr(counts)  ## default edgeR filter
        pgx$filtered <- NULL
        pgx$filtered[["low.expressed"]] <-
            paste(rownames(counts)[which(!keep)],collapse=";")
        table(keep)
        counts <- counts[which(keep),,drop=FALSE]
        genes <- genes[which(keep),,drop=FALSE]
        cat("filtering out",sum(!keep),"low-expressed genes\n")
        cat("keeping",sum(keep),"expressed genes\n")
    }
    
    ##-----------------------------------------------------------------------------
    ## Shrink number of genes before testing (highest SD/var)
    ##-----------------------------------------------------------------------------
    if(is.null(max.features)) max.features <- -1
    if(max.features > 0 && nrow(counts) > max.features) {
        cat("shrinking data matrices: n=",max.features,"\n")
        ##avg.prior.count <- mean(PRIOR.CPM * Matrix::colSums(counts) / 1e6)  ##
        ##logcpm = edgeR::cpm(counts, log=TRUE, prior.count=avg.prior.count)
        ##logcpm <- log2(PRIOR.CPM + edgeR::cpm(counts, log=FALSE))
        logcpm <- logCPM(counts, total=NULL)
        sdx <- apply(logcpm,1,sd)
        jj <- Matrix::head( order(-sdx), max.features )  ## how many genes?
        ## always add immune genes??
        if(FALSE && "gene_biotype" %in% colnames(genes)) {
            imm.gene <- grep("^TR_|^IG_",genes$gene_biotype)
            imm.gene <- imm.gene[which(sdx[imm.gene] > 0.001)]
            jj <- unique(c(jj,imm.gene))
        }
        jj0 <- setdiff(1:nrow(counts),jj)
        ##pgx$filtered[["low.variance"]] <- NULL
        pgx$filtered[["low.variance"]] <- paste(rownames(counts)[jj0],collapse=";")
        counts <- counts[jj,]
        genes <- genes[jj,]
    }
    Matrix::head(genes)
    genes  = genes[,c("gene_name","gene_title")]
    dim(counts)
    
    ##-----------------------------------------------------------------------------
    ## Do the fitting
    ##-----------------------------------------------------------------------------
    methods <- test.methods
    methods
    cat(">>> Testing differential expressed genes (DEG) with methods:",methods,"\n")

    ## Run all test methods
    ##
    ##X=counts;design=design,
    X <- pgx$X[rownames(counts),colnames(counts)]
    ##X <- X[rownames(counts),colnames(counts)]
    dim(X)

    message("[compute.testGenesSingleOmics] 12 : start fitting... ")
    
    ## quantile.normalize=TRUE;remove.batch=FALSE;conform.output=TRUE;do.filter=FALSE;custom=NULL;custom.name=NULL

    gx.meta <- ngs.fitContrastsWithAllMethods(
        counts = counts,
        X = X, ## type = type,
        samples = samples,
        genes = NULL, ##genes=genes,
        methods = methods,
        design = design,
        contr.matrix = contr.matrix,
        prune.samples = prune.samples,
        prior.cpm = PRIOR.CPM,  ## prior count regularization
        ## quantile.normalize = TRUE,  ## only for logCPM???
        remove.batch = FALSE,  ## we do explicit batch correction instead
        conform.output = TRUE,
        do.filter = FALSE,
        correct.AveExpr = TRUE,
        custom = NULL, custom.name = NULL
    )

    message("[compute.testGenesSingleOmics] 13 : fitting done!")
    
    names(gx.meta)
    names(gx.meta$outputs)
    print(gx.meta$timings)
    
    ##--------------------------------------------------------------------------------
    ## set default matrices
    ##--------------------------------------------------------------------------------
    
    rownames(gx.meta$timings) <- paste0("[test.genes]",rownames(gx.meta$timings))
    pgx$timings <- rbind(pgx$timings, gx.meta$timings)
    gx.meta$timings <- NULL
    gx.meta$X <- NULL
    ##pgx$genes = pgx$genes[rownames(pgx$X),]
    ##pgx$Y = pgx$samples[colnames(pgx$X),]
    pgx$model.parameters <- model.parameters
    pgx$gx.meta <- gx.meta
    ## pgx$X = gx.meta$X
    pgx$X <- X  ## replace with filtered
    
    ## remove large outputs... (uncomment if needed!!!)
    if(remove.outputs) {
        pgx$gx.meta$outputs <- NULL
    }

    message("[compute.testGenesSingleOmics] done!")
    
    return(pgx)
}

compute.testGenesMultiOmics <- function(pgx, contr.matrix, max.features=1000, 
                                        test.methods=c("trend.limma","deseq2.wald","edger.qlf"),
                                        use.design=TRUE, prune.samples = FALSE,
                                        remove.outputs=TRUE )
{
    pgx$gx.meta <- NULL
    pgx$model.parameters <- NULL
    pgx$gx.meta$meta <- vector("list",ncol(contr.matrix))
    pgx$X <- c()
    pgx$timings <- c()
    for(j in 1:4) {
        nk <- ncol(contr.matrix)
        pgx$gx.meta$sig.counts[[j]] <- vector("list",nk)
    }

    data.type <- gsub("\\[|\\].*","",rownames(pgx$counts))
    data.types <- unique(data.type)
    data.types
    dt = "cn"
    dt = "gx"
    dt <- data.types[1]
    dt
    for(dt in data.types) {
        
        ## get data block
        pgx1 <- pgx
        jj <- which(data.type == dt)
        pgx1$counts <- pgx1$counts[jj,]
        pgx1$genes  <- pgx1$genes[jj,]
        
        ## determine if datatype are counts or not
        type = "not.counts"
        if(min(pgx1$counts,na.rm=TRUE) >= 0 &&
           max(pgx1$counts,na.rm=TRUE) >= 50 ) {
            type <- "counts"
        }
        dt
        type
        
        ## do test
        pgx1 <- compute.testGenesSingleOmics(
            pgx=pgx1, type=type,
            contr.matrix=contr.matrix,
            max.features=max.features,
            test.methods=test.methods)
        
        ## copy results
        pgx$model.parameters <- pgx1$model.parameters
        names(pgx1$gx.meta)
        for(k in 1:ncol(contr.matrix)) {
            pgx$gx.meta$meta[[k]] <- rbind(pgx$gx.meta$meta[[k]],
                                           pgx1$gx.meta$meta[[k]])
        }
        names(pgx$gx.meta$meta) <- names(pgx1$gx.meta$meta)
        for(j in 1:4) {
            nk <- ncol(contr.matrix)
            for(k in 1:nk) {
                cnt1 <- pgx1$gx.meta$sig.counts[[j]][[k]]
                cnt0 <- pgx$gx.meta$sig.counts[[j]][[k]]
                rownames(cnt1) <- paste0("[",dt,"]",rownames(cnt1))
                pgx$gx.meta$sig.counts[[j]][[k]] <- rbind(cnt0, cnt1)
            }
            names(pgx$gx.meta$sig.counts[[j]]) <- names(pgx1$gx.meta$sig.counts[[j]])            
        }
        names(pgx$gx.meta$sig.counts) <- names(pgx1$gx.meta$sig.counts)
        pgx$timings <- rbind(pgx$timings, pgx1$timings)
        pgx$X <- rbind(pgx$X, pgx1$X)
    }

    gg <- rownames(pgx$counts)
    pgx$X <- pgx$X[match(gg,rownames(pgx$X)),]
    ##pgx$genes <- pgx$genes[match(gg,rownames(pgx$genes)),]
    pgx$model.parameters <- pgx1$model.parameters
    return(pgx)
}

## ---------- clean up ----------------
##contr.matrix <- contr.matrix0  ## RESTORE
##rm(list=setdiff(ls(),SAVE.PARAMS))

