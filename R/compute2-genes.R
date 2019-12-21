## Workflow for differential expression analysis (at gene level)
##
## Input:  contr.matrix to be defined using group as levels
##
##
##
##
##

##SAVE.PARAMS <- ls()
max.features=1000;type="counts";test.methods=c("ttest","trend.limma")

compute.testGenes <- function(ngs, contr.matrix, max.features=1000, type="counts",
                              test.methods=c("trend.limma","deseq2.wald","edger.qlf"))
{
    single.omics <- !any(grepl("\\[",rownames(ngs$counts)))
    single.omics
    data.types <- unique(gsub("\\[|\\].*","",rownames(ngs$counts)))
    ##data.types
    if(single.omics || length(data.types)==1) {
        ## single-omics, no missing values
        cat(">>> computing gene tests for SINGLE-OMICS\n")
        ngs <- compute.testGenesSingleOmics(
            ngs=ngs, type=type,
            contr.matrix=contr.matrix,
            max.features=max.features,
            test.methods=test.methods)
    } else {
        ## multi-omics, missing values allowed
        cat(">>> computing gene tests for MULTI-OMICS\n")
        ngs <- compute.testGenesMultiOmics(
            ngs=ngs,  ## type is inferred
            contr.matrix=contr.matrix,
            max.features=max.features,
            test.methods=test.methods)
    }
    return(ngs)
}

test.methods=c("trend.limma","deseq2.wald","edger.qlf")
test.methods=c("trend.limma","ttest.welch","ttest")
max.features=1000
compute.testGenesMultiOmics <- function(ngs, contr.matrix, max.features=1000, 
                               test.methods=c("trend.limma","deseq2.wald","edger.qlf"))
{
    ngs$gx.meta <- NULL
    ngs$model.parameters <- NULL
    ngs$gx.meta$meta <- vector("list",ncol(contr.matrix))
    ngs$X <- c()
    ngs$timings <- c()
    for(j in 1:4) {
        nk <- ncol(contr.matrix)
        ngs$gx.meta$sig.counts[[j]] <- vector("list",nk)
    }

    data.type <- gsub("\\[|\\].*","",rownames(ngs$counts))
    data.types <- unique(data.type)
    data.types
    dt = "cn"
    dt = "gx"
    dt <- data.types[1]
    dt
    for(dt in data.types) {
        
        ## get data block
        ngs1 <- ngs
        jj <- which(data.type == dt)
        ngs1$counts <- ngs1$counts[jj,]
        ngs1$genes  <- ngs1$genes[jj,]
        
        ## determine if datatype are counts or not
        type = "not.counts"
        if(min(ngs1$counts,na.rm=TRUE) >= 0 &&
           max(ngs1$counts,na.rm=TRUE) >= 50 ) {
            type <- "counts"
        }
        dt
        type
        
        ## do test
        ngs1 <- compute.testGenesSingleOmics(
            ngs=ngs1, type=type,
            contr.matrix=contr.matrix,
            max.features=max.features,
            test.methods=test.methods)
        
        ## copy results
        ngs$model.parameters <- ngs1$model.parameters
        names(ngs1$gx.meta)
        for(k in 1:ncol(contr.matrix)) {
            ngs$gx.meta$meta[[k]] <- rbind(ngs$gx.meta$meta[[k]],
                                           ngs1$gx.meta$meta[[k]])
        }
        names(ngs$gx.meta$meta) <- names(ngs1$gx.meta$meta)
        for(j in 1:4) {
            nk <- ncol(contr.matrix)
            for(k in 1:nk) {
                cnt1 <- ngs1$gx.meta$sig.counts[[j]][[k]]
                cnt0 <- ngs$gx.meta$sig.counts[[j]][[k]]
                rownames(cnt1) <- paste0("[",dt,"]",rownames(cnt1))
                ngs$gx.meta$sig.counts[[j]][[k]] <- rbind(cnt0, cnt1)
            }
            names(ngs$gx.meta$sig.counts[[j]]) <- names(ngs1$gx.meta$sig.counts[[j]])            
        }
        names(ngs$gx.meta$sig.counts) <- names(ngs1$gx.meta$sig.counts)
        ngs$timings <- rbind(ngs$timings, ngs1$timings)
        ngs$X <- rbind(ngs$X, ngs1$X)
    }

    gg <- rownames(ngs$counts)
    ngs$X <- ngs$X[match(gg,rownames(ngs$X)),]
    ##ngs$genes <- ngs$genes[match(gg,rownames(ngs$genes)),]
    ngs$model.parameters <- ngs1$model.parameters
    return(ngs)
}

compute.testGenesSingleOmics <- function(ngs, contr.matrix, max.features=1000,
                                type="counts", filter.low = TRUE,
                                test.methods = c("trend.limma","deseq2.wald","edger.qlf"))
{

    ##-----------------------------------------------------------------------------
    ## Check parameters, decide group level
    ##-----------------------------------------------------------------------------    
    if(!("counts" %in% names(ngs))) {
        stop("cannot find counts in ngs object")
    }
    
    group = NULL
    if( all(rownames(contr.matrix) %in% ngs$samples$group)) {
        cat("testing on groups...\n")
        group = as.character(ngs$samples$group)
    }
    if( all(rownames(contr.matrix) %in% ngs$samples$cluster)) {
        cat("testing on clusters...\n")
        group = as.character(ngs$samples$cluster)
    }
    if( all(rownames(contr.matrix) %in% rownames(ngs$samples))) {
        cat("testing on samples...\n")
        group = rownames(ngs$samples)
    }
    
    if(is.null(group)) {
        stop("invalid contrast matrix. could not assign groups")
    }
    
    table(group)
    ##dim(contr.matrix)
    
    ##-----------------------------------------------------------------------------
    ## normalize contrast matrix to zero mean and signed sums to one
    ##-----------------------------------------------------------------------------
    contr.matrix[is.na(contr.matrix)] <- 0
    contr.matrix0 <- contr.matrix  ## SAVE
    
    ## take out any empty comparisons
    contr.matrix <- contr.matrix0[,which(colSums(contr.matrix0!=0)>0),drop=FALSE]
    for(i in 1:ncol(contr.matrix)) {
        m <- contr.matrix[,i]
        m[is.na(m)] <- 0
        contr.matrix[,i] <- 1*(m>0)/sum(m>0) - 1*(m<0)/sum(m<0)
    }
    dim(contr.matrix)
    
    ##-----------------------------------------------------------------------------
    ## create design matrix from defined contrasts (group or clusters)
    ##-----------------------------------------------------------------------------
    
    no.design <- all(group %in% rownames(ngs$samples))  ## sample-wise design
    no.design
    design=NULL
    
    if(no.design) {
        ## SAMPLE-WISE DESIGN
        design=NULL
        exp.matrix <- contr.matrix
    } else {
        ## GROUP DESIGN
        ##group[is.na(group)] <- "_"
        group[which(!group %in% rownames(contr.matrix))] <- "_"
        design <- model.matrix(~ 0 + group )  ## clean design no batch effects...
        colnames(design) <- sub("^group", "", colnames(design))
        rownames(design) <- colnames(ngs$counts)
        design
        
        ## check contrasts for sample sizes (at least 2 in each group) and
        ## remove otherwise
        design <- design[,match(rownames(contr.matrix),colnames(design)),drop=FALSE]
        colnames(design)
        design = design[,rownames(contr.matrix),drop=FALSE]
        exp.matrix = (design %*% contr.matrix)
        keep <- rep(TRUE,ncol(contr.matrix))
        keep = (colSums(exp.matrix > 0) >= 1 & colSums(exp.matrix < 0) >= 1)
        ##keep = ( colSums(exp.matrix > 0) >= 2 & colSums(exp.matrix < 0) >= 2 )
        table(keep)
        contr.matrix = contr.matrix[,keep,drop=FALSE]
        exp.matrix = (design %*% contr.matrix)
    }

    model.parameters <- list(design = design,
                             contr.matrix = contr.matrix, 
                             exp.matrix = exp.matrix)
    
    ##-----------------------------------------------------------------------------
    ## Filter genes
    ##-----------------------------------------------------------------------------    
    ## get *RAW* counts but use filtered probes from cooked
    counts = ngs$counts  ##??
    genes  = ngs$genes
    samples = ngs$samples
    
    ## prefiltering for low-expressed genes (recommended for edgeR and
    ## DEseq2). Require at least in 2 or 1% of total. Specify the
    ## PRIOR CPM amount to regularize the counts and filter genes
    if(type=="counts" && filter.low) {
        PRIOR.CPM = 0.25
        PRIOR.CPM = 1
        AT.LEAST = ceiling(pmax(2,0.01*ncol(counts)))    
        cat("filtering for low-expressed genes: >",PRIOR.CPM,"CPM in >=",
            AT.LEAST,"samples\n")
        keep <- (rowSums( edgeR::cpm(counts) > PRIOR.CPM, na.rm=TRUE) >= AT.LEAST)
        ##keep <- edgeR::filterByExpr(counts)  ## default edgeR filter
        ngs$filtered <- NULL
        ngs$filtered[["low.expressed"]] <-
            paste(rownames(counts)[which(!keep)],collapse=";")
        table(keep)
        counts <- counts[which(keep),,drop=FALSE]
        genes <- genes[which(keep),,drop=FALSE]
        cat("filtering out",sum(!keep),"low-expressed genes\n")
        cat("keeping",sum(keep),"expressed genes\n")
    }
    
    ##-----------------------------------------------------------------------------
    ## Shrink number of genes before testing
    ##-----------------------------------------------------------------------------
    if(is.null(max.features)) max.features <- -1
    if(max.features > 0 && nrow(counts) > max.features) {
        cat("shrinking data matrices: n=",max.features,"\n")
        ##avg.prior.count <- mean(PRIOR.CPM * Matrix::colSums(counts) / 1e6)  ##
        ##logcpm = edgeR::cpm(counts, log=TRUE, prior.count=avg.prior.count)
        if(type=="counts") {
            logcpm <- log2(PRIOR.CPM + edgeR::cpm(counts, log=FALSE))
            sdx <- apply(logcpm,1,sd)
        } else {
            sdx <- apply(counts,1,sd)
        }
        jj <- head( order(-sdx), max.features )  ## how many genes?
        ## always add immune genes??
        if("gene_biotype" %in% colnames(genes)) {
            imm.gene <- grep("^TR_|^IG_",genes$gene_biotype)
            imm.gene <- imm.gene[which(sdx[imm.gene] > 0.001)]
            jj <- unique(c(jj,imm.gene))
        }
        jj0 <- setdiff(1:nrow(counts),jj)
        ##ngs$filtered[["low.variance"]] <- NULL
        ngs$filtered[["low.variance"]] <- paste(rownames(counts)[jj0],collapse=";")
        counts <- counts[jj,]
        genes <- genes[jj,]
    }
    head(genes)
    genes  = genes[,c("gene_name","gene_title")]
    dim(counts)
    
    ##-----------------------------------------------------------------------------
    ## Do the fitting
    ##-----------------------------------------------------------------------------
    ## Select test methods
    ##
    ## all.methods=c("ttest","ttest.welch","voom.limma","trend.limma","notrend.limma",
    ##               "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
    ## methods=c("trend.limma","edger.qlf","deseq2.wald")
    ## if(ncol(counts)>500) methods=c("trend.limma","edger.qlf","edger.lrt")
    ## methods
    ## if(!is.null(USER.GENETEST.METHODS)) methods = USER.GENETEST.METHODS
    ## if(methods[1]=="*") methods = all.methods
    methods <- test.methods
    
    cat(">>> Testing differential expressed genes (DEG) with methods:",methods,"\n")

    ## Run all test methods
    ##
    gx.meta <- ngs.fitContrastsWithAllMethods(
        X=counts, type=type,
        samples=samples, genes=NULL, ##genes=genes,
        methods=methods, design=design,
        contr.matrix=contr.matrix,
        prior.cpm=PRIOR.CPM,  ## prior count regularization
        quantile.normalize=TRUE,  ## only for logCPM
        remove.batch=FALSE,  ## we do explicit batch correction instead
        conform.output=TRUE, do.filter=FALSE,
        custom=NULL, custom.name=NULL )
    
    names(gx.meta)
    names(gx.meta$outputs)
    names(gx.meta$outputs[[1]])
    names(gx.meta$outputs[[1]][[1]])
    
    print(gx.meta$timings)
    
    ##--------------------------------------------------------------------------------
    ## set default matrices
    ##--------------------------------------------------------------------------------
    
    rownames(gx.meta$timings) <- paste0("[test.genes]",rownames(gx.meta$timings))
    ngs$timings <- rbind(ngs$timings, gx.meta$timings)
    ngs$X = gx.meta$X
    gx.meta$timings <- NULL
    gx.meta$X <- NULL
    ##ngs$genes = ngs$genes[rownames(ngs$X),]
    ##ngs$Y = ngs$samples[colnames(ngs$X),]
    ngs$model.parameters <- model.parameters
    ngs$gx.meta <- gx.meta
    
    ## remove large outputs... (uncomment if needed!!!)
    ngs$gx.meta$outputs <- NULL

    return(ngs)
}

## ---------- clean up ----------------
##contr.matrix <- contr.matrix0  ## RESTORE
##rm(list=setdiff(ls(),SAVE.PARAMS))

