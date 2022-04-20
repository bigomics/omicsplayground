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

##--------------------------------------------------------------------------------------------
##----------------------------- COOKING (i.e. preparing data) --------------------------------
##--------------------------------------------------------------------------------------------

##txi = ngs.importSalmon(sf.files)
ngs.cookForEDGER <- function( counts, samples=NULL, genes=NULL, normalization="none",
                             filter=TRUE, prior.cpm=0, remove.batch=TRUE )
{
    

    if(all(c("counts","samples","genes") %in% names(counts))) {
        samples = counts$samples
        genes = counts$genes
        counts = counts$counts
    }
    ##txi.names = c("abundance","counts","length","countsFromAbundance")
    if(!all(colnames(counts)==rownames(samples))) stop("samples do not match")
    if(!all(rownames(counts)==rownames(genes))) stop("genes do not match")
    if(is.null(samples)) stop("need samples specified")
    if(is.null(genes)) stop("need genes specified")

    ##------------------------------------------------------------
    ## Now create an DGEList object (see also tximport Vignette)
    ##------------------------------------------------------------
    
    
    cooked <- edgeR::DGEList( round(counts), group=NULL)  ## we like integer counts...
    cooked$samples$group <- NULL
    cooked$samples = cbind(cooked$samples, samples)
    if(!is.null(genes)) cooked$genes = genes

    ## filter out non-expressed genes (using edgeR standard function)
    if(filter) {
        keep <- edgeR::filterByExpr(cooked)
        table(keep)
        cooked <- cooked[keep,]
    }

    ## normalized for RNA composition (TMM)
    cooked <- edgeR::calcNormFactors(cooked, method="TMM")
    ##cooked <- edgeR::calcNormFactors(cooked, method="TMM")

    ##------------------------------------------------------------------
    ## prior count regularization
    ##------------------------------------------------------------------
    ## based on the above, we add a prior count
    ##PRIOR.CPM = 3  ## MAGIC NUMBER .... FIDDLE....
    if(prior.cpm > 0) {
        cat("adding prior counts at prior.cpm=",prior.cpm,"\n")
        CPM.FACTORS = Matrix::colSums(counts)  / 1e6
        prior.counts = (prior.cpm * CPM.FACTORS)
        summary(prior.counts)
        ##PRIOR.COUNT = 0.5
        cooked$counts = t(t(cooked$counts) + prior.counts)
    }

    if(normalization!="none") {

        ## remove batch-effects with LIMMA. Be sure to include batch in the
        ## model to avoid overfitting.
        has.batch <-  ("batch" %in% colnames(cooked$samples))
        if(has.batch && remove.batch==TRUE) {
            cat("Found 'batch' column in sample table. Correcting for batch effects...\n")
            ##if(batch[1] %in% colnames(cooked$samples)) batch <- cooked$samples[,batch]
            design1 = model.matrix( ~ group, data=cooked$samples)
            batch1 = cooked$samples[,"batch"]
            xnorm <- limma::removeBatchEffect( xnorm, batch=batch1, design=design1)
        } else {
            ##
        }

        ## quantile normalize. be careful, may introduce artifacts (e.g. clipping)
        if(normalization=="quantile") {
            xnorm <- limma::normalizeQuantiles(xnorm)
        }

        ## we need to undo log and normalizations for further analysis???
        cooked$counts <- 2**(xnorm)  ## undo log2
        cooked$samples$norm.factors <- 1  ## we assumed its already used
        cooked$samples$lib.size = round(Matrix::colSums(cooked$counts)) ## update lib.size

    }

    ## stop here???
    ##res = list( cooked=cooked, raw=raw )
    return(cooked)
}


ngs.cookForDESEQ2 <- function(counts, samples, genes, remove.batch=TRUE,
                              test="Wald", prior.cpm=0, filter=TRUE )
{
    

    if(all(c("counts","samples","genes") %in% names(counts))) {
        samples = counts$samples
        genes = counts$genes
        counts = counts$counts
    }

    ##txi.names = c("abundance","counts","length","countsFromAbundance")
    ##if(!all(txi.names %in% names(txi))) stop("txi must be a tximport object!")
    if(!all(colnames(counts)==rownames(samples))) stop("samples do not match")
    if(!all(rownames(counts)==rownames(genes))) stop("genes do not match")
    if(is.null(samples)) stop("need samples specified")
    if(is.null(genes)) stop("need genes specified")

    ##------------------------------------------------------------------
    ## remove genes with zero values (deseq does this like this?)
    ##------------------------------------------------------------------
    if(filter) {
        keep = rowSums(counts) > 0
        keep = (rowSums(counts) > 100) > 3
        keep = rowSums(counts) > 0  ## check original count because of prior.cpm
        table(keep)
        counts <- counts[keep, ]
        genes <- genes[keep,]
    }

    ##------------------------------------------------------------------
    ## prior count regularization
    ##------------------------------------------------------------------
    ##prior.cpm = 0
    ##prior.cpm = 3  ## MAGIC NUMBER .... FIDDLE....
    prior.cpm
    if(!is.null(prior.cpm) && prior.cpm>0 ) {
        cat("adding prior counts at prior.cpm=",prior.cpm,"\n")
        CPM.FACTORS = Matrix::colSums(counts)  / 1e6
        prior.counts = (prior.cpm * CPM.FACTORS)
        prior.counts <- pmax(prior.counts,1)
        summary(prior.counts)
        ##PRIOR.COUNT = 0.5
        counts = t(t(counts) + prior.counts)
    }

    ##------------------------------------------------------------
    ## Now create an DEseq2 object (see also tximport Vignette)
    ##------------------------------------------------------------
    Matrix::head(samples)
    
    if(!("group" %in% colnames(samples))) {
        stop("samples must have 'group' column for DESeq")
    } else {
         cat("found 'group' column in sample table\n")
    }
    counts = round(counts)
    mode(counts) = "integer"
    colnames(samples)
    nbatch <- sum(grepl("batch",colnames(samples)))
    nbatch
    if(remove.batch && nbatch >1 ) {
        ## It is good practice to add the batch parameters in the
        ## final model, even if they have been remove before. But with
        ## multiple batch effects it often creates '(in)dependency'
        ## problems...
        ##
        cat("found multiple 'batch' columns in sample table\n")
        batch = colnames(samples)[grep("^batch",colnames(samples))]
        if(length(batch)>1) {
            ## ONLY ONE BATCH FOR NOW!!!!!
            batch = batch[1]
            cat("WARNING: only correcting one batch implemented\n")
        }
        design.formula = formula(paste( c("~ 0 + group",batch),collapse=" + "))  ## multiple batch
    ##} else if(is.null(design)) {
    ##    design.formula = formula(" ~ 1")
    } else {
        design.formula = formula(" ~ 0 + group")
    }
    cat("using model design: ",as.character(design.formula),"\n")

    rownames.counts <- rownames(counts)
    rownames(counts) <- NULL
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData=counts, design=design.formula, colData=data.frame(samples))
    rownames(counts) <- rownames.counts
    ##dds <- DESeq2::DESeqDataSetFromMatrix(
    ##counts, design = ~ 0 + group + batch, colData = samples[,c("group","batch")])

    ##------------------------------------------------------------------
    ##to collapse TECHNICAL replicates???
    ##------------------------------------------------------------------
    if(FALSE && "replicate" %in% colnames(samples)) {
        ##repl.id <- factor(paste(dds$clone, dds$strain, dds$treatment,sep="_"))
        repl.id <- samples$replicate
        dds <- DESeq2::collapseReplicates(dds, repl.id, dds$sample)
    }

    ##The following code estimates size factors to account for
    ##differences in sequencing depth. So the value are typically
    ##centered around 1. If all the samples have exactly the same
    ##sequencing depth, you expect these numbers to be near 1.
    ##sf <- DESeq2::estimateSizeFactorsForMatrix(DESeq2::counts(dds))  ## where is it used???

    ##Run DESeq : Modeling counts with generic 'group'
    fitType = 'parametric'
    if(prior.cpm!=0) fitType ='local'
    dds <- DESeq2::DESeq(dds, fitType=fitType, test=test)
    DESeq2::resultsNames(dds)# lists the coefficients

    ## we add the gene annotation here (not standard...)
    ##colnames(SummarizedExperiment::rowData(dds))
    SummarizedExperiment::rowData(dds)$genes <- genes   ## does this work??

    return(dds)
}



##--------------------------------------------------------------------------------------------
##------------------------------------ END OF FILE -------------------------------------------
##--------------------------------------------------------------------------------------------
