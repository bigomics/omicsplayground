

source("../scripts/options.R")
MAX.GENES

if(0) {

    counts  = as.matrix(read.csv("../exampledata/counts.csv", row.names=1))
    samples = read.csv("../exampledata/samples.csv", row.names=1, stringsAsFactors=FALSE)
    genes   = read.csv("../exampledata/genes.csv", row.names=1, stringsAsFactors=FALSE)
    contrasts = as.matrix(read.csv("../exampledata/contrasts.csv", row.names=1))

}

pgx.upload <- function(counts, samples, genes, contrasts,
                       ##gx.methods = c("trend.limma","edger.qlf","deseq2.wald"),
                       gx.methods = c("ttest.welch","trend.limma","edger.qlf"),
                       gset.methods = c("fisher","gsva","fgsea"),
                       extra.methods = c("meta.go","deconv","infer","drugs"),
                       progress=NULL)
{

    library(org.Hs.eg.db)

    if(!is.null(progress)) progress$inc(0.01, detail = "creating object")
    
    ##load(file=rda.file, verbose=1)
    ngs <- list()  ## empty object
    ngs$name = "(uploaded)"
    ngs$date = date()
    ngs$datatype = "unknown"
    ngs$description = "uploaded data set"

    if(0) {
        colnames(counts) == rownames(samples)
        rownames(counts) == rownames(genes)
    }
    
    ##-------------------------------------------------------------------
    ## create ngs object
    ##-------------------------------------------------------------------
    ngs$samples = samples
    ngs$counts  = counts
    ngs$genes   = data.frame(genes)
    ngs$contrasts = contrasts

    cat("DBG [pgx-upload] 1: dim(ngs$counts)=",dim(ngs$counts),"\n")
    cat("DBG [pgx-upload] 1: sum.is.na(ngs$counts)=",sum(is.na(ngs$counts)),"\n")
    cat("DBG [pgx-upload] 1: dim(ngs$genes)=",dim(ngs$genes),"\n")
    
    if(!"gene_title" %in% colnames(ngs$genes)) {
        require(org.Hs.eg.db)
        GENE.TITLE  = unlist(as.list(org.Hs.egGENENAME))
        gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
        names(GENE.TITLE) = gene.symbol    
        gene = rownames(genes)
        ngs$genes$gene_title = gene_title = GENE.TITLE[gene]
        ##ngs$genes$chr  = gene_title = GENE.CHR[gene]
        ##ngs$genes$pos  = gene_title = GENE.POS[gene]
    }
    
    cat("DBG [pgx-upload] 2: dim(ngs$genes)=",dim(ngs$genes),"\n")
    cat("DBG [pgx-upload] 2: dim(ngs$counts)=",dim(ngs$counts),"\n")
    
    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    sum(duplicated(ngs$genes$gene_name))
    gg = as.character(ngs$genes$gene_name)
    x1 = apply( ngs$counts, 2, function(x) tapply(x, gg, sum))
    cat("DBG [pgx-upload] 2: dim(x1)=",dim(x1),"\n")

    ngs$genes = ngs$genes[match(rownames(x1), ngs$genes$gene_name),]
    ngs$counts = x1
    rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
    remove(x1)

    cat("DBG [pgx-upload] 3: dim(ngs$genes)=",dim(ngs$genes),"\n")
    cat("DBG [pgx-upload] 3: dim(ngs$counts)=",dim(ngs$counts),"\n")
    cat("DBG [pgx-upload] 3: sum.is.na(ngs$counts)=",sum(is.na(ngs$counts)),"\n")
    
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    if(!is.null(progress)) progress$inc(0.01, detail = "clustering")
    
    ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=3)
    head(ngs$samples)
    table(ngs$samples$cluster)
    

    ##======================================================================
    ##======================================================================
    ##======================================================================


    ## make model matrix for group vs. rest
    ##contr.matrix <- makeClusterContrasts(ngs$samples$cluster)
    contr.matrix <- ngs$contrasts
    if(!"group" %in% colnames(ngs$samples)) {
        ct = apply(ngs$samples,2,function(px) mean(px %in% rownames(ngs$contrasts),na.rm=TRUE))
        group.col = which(ct > 0.95)
        if(length(group.col)>0) {
            grp.var = colnames(ngs$samples)[group.col[1]]
            cat(paste0("INFO [pgx-upload] assigning '",grp.var,"' variable as 'group' column\n"))
            colnames(ngs$samples)[group.col[1]] <- "group"
        } else {
            stop("sample annotation file must have 'group' column\n")
        }
    }

    ##======================================================================
    ##======================================================================
    ##======================================================================

    ngs$timings <- c()
    
    USER.GENETEST.METHODS=c("ttest","ttest.welch","ttest.rank",
                            "voom.limma","trend.limma","notrend.limma",
                            "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
    USER.GENESETTEST.METHODS = c("fisher","gsva","ssgsea","spearman",
                                 "camera", "fry","fgsea") ## no GSEA, too slow...

    if(0) {
        source("../R/compute-genes.R")
        source("../R/compute-genesets.R")
        source("../R/compute-extra.R")
    } else {

        ## new callling methods
        source("../R/compute2-genes.R")
        source("../R/compute2-genesets.R")
        source("../R/compute2-extra.R")

        ## ------------------ gene level tests ---------------------
        if(!is.null(progress)) progress$inc(0.1, detail = "testing genes")

        ngs <- compute.testGenes(
            ngs, contr.matrix, max.features=MAX.GENES,
            test.methods=gx.methods)
        head(ngs$gx.meta$meta[[1]])        
        
        ## ------------------ gene set tests -----------------------
        if(!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")
        
        ngs <- compute.testGenesets(
            ngs, max.features=MAX.GENES,
            test.methods=gset.methods)
        head(ngs$gset.meta$meta[[1]])

        ## ------------------ extra analyses ---------------------
        if(!is.null(progress)) progress$inc(0.3, detail = "computing extra modules")
        ##extra <- c("meta.go","deconv","infer","drugs")
        ##extra <- c("meta.go","infer","drugs")
        ngs <- compute.extra(ngs, extra=extra.methods)

    }

    ngs$timings
    return(ngs)
}








