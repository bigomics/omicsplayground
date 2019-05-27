## library(limma)
## source("../R/gx-heatmap.r")
## source("../R/gx-limma.r")
## source("../R/gx-util.r")
## source("../R/gset-fisher.r")
## source("../R/gset-gsea.r")
## source("../R/gset-meta.r")

SAVE.PARAMS <- ls()

compute.testGenesets <- function(ngs, max.features=1000,
                                 test.methods = c("gsva","camera","fgsea"))
{
    
    ##-----------------------------------------------------------
    ## Load gene sets
    ##-----------------------------------------------------------    
    pp <- rownames(ngs$counts)
    is.mouse = (mean(grepl("[a-z]",gsub(".*:|.*\\]","",pp))) > 0.8)
    is.mouse
    if(is.mouse) {
        cat("Loading MOUSE gene sets...\n")
        load(file=file.path(FILES,"gmt-all-mouse.rda"))
    } else {
        cat("Loading HUMAN gene sets...\n")
        load(file=file.path(FILES,"gmt-all.rda"))
    }
    table(sub(":.*","",names(gmt.all)))
    summary(sapply(gmt.all,length))
    
    ##-----------------------------------------------------------
    ## Filter gene sets
    ##-----------------------------------------------------------
    cat("Filtering gene sets...\n")
    
    require(Matrix)
    require(org.Hs.eg.db)
    ##GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
    ##genes = head(as.character(unlist(as.list(org.Hs.egSYMBOL))),1000)
    genes = unique(as.character(ngs$genes$gene_name))
    gmt.all = mclapply(gmt.all, function(gs) intersect(gs, genes))
    table(sub(":.*","",names(gmt.all)))
    
    ## filter gene sets on size
    gmt.size = sapply(gmt.all,length)
    summary(gmt.size)
    gmt.all = gmt.all[which(gmt.size >= 15 & gmt.size <= 1e3)]
    gmt.size = sapply(gmt.all,length)
    gmt.all = gmt.all[order(-gmt.size)]
    summary(gmt.size)
    gmt.all = gmt.all[!duplicated(names(gmt.all))]
    length(gmt.all)
    table(sub(":.*","",names(gmt.all)))
    
    ##-----------------------------------------------------------
    ## create the full GENE matrix (always collapsed by gene)
    ##-----------------------------------------------------------
    single.omics <- !any(grepl("\\[",rownames(ngs$counts)))
    single.omics
    if(single.omics) {
        X <- ngs$counts
    } else {
        data.type <- gsub("\\[|\\].*","",rownames(ngs$counts))
        jj <- which(data.type %in% c("gx","mrna"))
        length(jj)
        if(length(jj)==0) {
            stop("FATAL. could not find gx/mrna values.")
        }
        X <- ngs$counts[jj,]
    }
    xgenes = as.character(ngs$genes[rownames(X),"gene_name"])
    X <- apply(X, 2, function(x) tapply(x, xgenes, sum))
    dim(X)
    X <- edgeR::cpm(X, log=TRUE )
    X <- X[which(apply(X,1,sd)>0),,drop=FALSE]
    X <- X[!(rownames(X) %in% c(NA,""," ") ),,drop=FALSE]
    X <- limma::normalizeQuantiles(X)
    dim(X)
    
    ##-----------------------------------------------------------
    ## create the GENESETxGENE matrix
    ##-----------------------------------------------------------
    cat("Building gene set matrix...\n")
    ##GMT = sapply( gmt.all, function(s) 1*(rownames(X) %in% s))
    ##GMT = Matrix(GMT, sparse=TRUE)
    genes <- sort(unique(xgenes))
    GMT <- gmt2mat.nocheck(gmt.all[], bg=genes)  ## in gset-gsea.r
    dim(GMT)
    dim(X)
    table(xgenes %in% rownames(GMT))
    table(names(gmt.all) %in% colnames(GMT))
    GMT[1:4,1:4]
    
    ## align GMT to X (or ngs$X??)
    GMT <- GMT[rownames(X),names(gmt.all)]
    summary(Matrix::colSums(GMT))
    dim(GMT)
    
    ##-----------------------------------------------------------
    ## Prioritize gene sets by fast rank-correlation
    ##-----------------------------------------------------------
    
    if(is.null(max.features)) max.features <- 20000
    if(max.features < 0) max.features <- 20000
    max.features
    
    if(max.features > 0) {
        require(limma)
        ## Reduce gene sets by selecting top varying genesets. We use the
        ## very fast sparse rank-correlation for approximate single sample
        ## geneset activation.
        cX <- X - rowMeans(X, na.rm=TRUE)
        gsetX = qlcMatrix::corSparse( GMT[,], apply( cX[,],2,rank) )
        gsetX = limma::normalizeQuantiles(gsetX) ##???
        grp <- ngs$samples$group
        gsetX.bygroup <- t(apply(gsetX,1,function(x) tapply(x,grp,mean)))
        sdx <- apply(gsetX.bygroup,1,sd)
        names(sdx) <- colnames(GMT)
        jj = head(order(-sdx), max.features) 
        must.include <- "hallmark|kegg|^go|^celltype"
        jj = unique( c(jj, grep(must.include,colnames(GMT),ignore.case=TRUE)))
        jj = jj[order(colnames(GMT)[jj])]
        length(jj)
        GMT = GMT[,jj]
        gmt.all = gmt.all[colnames(GMT)]
    }
    dim(GMT)
    
    sum(duplicated(rownames(GMT)))
    dim(X)
    dim(GMT)
    ngs$GMT <- GMT
    ngs$gmt.all <- gmt.all
    
    ##-----------------------------------------------------------
    ## get design and contrast matrix
    ##-----------------------------------------------------------
    design = ngs$model.parameters$design
    contr.matrix = ngs$model.parameters$contr.matrix
    ##contr.matrix
    ##exp.matrix = (design %*% contr.matrix)
    
    all.gset.methods=c("fisher","ssgsea","gsva", "spearman", "camera", "fry",
                       "fgsea","gsea.permPH","gsea.permGS","gseaPR")
    ##test.methods = c("fisher","fgsea")
    ## test.methods = c("fisher","gsva","ssgsea","spearman","camera","fry","fgsea") ## no GSEA
    ## if(!is.null(USER.GENESETTEST.METHODS)) test.methods = USER.GENESETTEST.METHODS
    ## ##if(test.methods[1]=="*") test.methods = all.gset.methods
    test.methods
    
    ##-----------------------------------------------------------
    ## Run methods
    ##-----------------------------------------------------------
    cat(">>> Testing gene sets with methods:",test.methods,"\n")
    test.methods
    
    Y <- ngs$samples
    ##gmt=ngs$gmt.all[1:100]
    gmt=ngs$gmt.all[]
    gc()
    
    gset.meta = gset.fitContrastsWithAllMethods(
        gmt = ngs$gmt.all, X = X, Y = Y, design=design, ## genes=GENES,
        contr.matrix=contr.matrix, methods=test.methods,
        mc.threads=1, mc.cores=NULL, batch.correct=TRUE )
    
    print(gset.meta$timings)
    
    rownames(gset.meta$timings) <- paste0("[test.genesets]",rownames(gset.meta$timings))
    ngs$timings <- rbind(ngs$timings, gset.meta$timings)
    ngs$gset.meta <- gset.meta
    ngs$gsetX = ngs$gset.meta$matrices[["meta"]]  ## META??!
    
    ##-----------------------------------------------------------------------
    ##------------------------ clean up -------------------------------------
    ##-----------------------------------------------------------------------
    
    ## remove large outputs... (uncomment if needed)
    ngs$gset.meta$outputs <- NULL
    
    remove(X)
    remove(Y)

    return(ngs)
}


##rm(list=setdiff(ls(),SAVE.PARAMS))

