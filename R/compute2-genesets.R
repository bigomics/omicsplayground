##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


##max.features=8000;lib.dir=FILES;test.methods = c("gsva","camera","fgsea")
compute.testGenesets <- function(ngs, max.features=1000, lib.dir="../lib",
                                 test.methods = c("gsva","camera","fgsea"),
                                 remove.outputs=TRUE )
{
    ## Rewritten 24.12.2019. Now much faster but needs gset-sparseG-XL
    ## precomputed.
    ##
    ##
    if(!"X" %in% names(ngs)) {
        stop("[compute.testGenesets] FATAL : object must have normalized matrix X")
    }
    
    ##-----------------------------------------------------------
    ## Load huge geneset matrix
    ##-----------------------------------------------------------    
    G <- t(readRDS(file.path(lib.dir,"gset-sparseG-XL.rds")))
    dim(G)
    
    ##-----------------------------------------------------------
    ## Filter genes
    ##-----------------------------------------------------------

    ## filter genes only in dataset
    require(Matrix)
    require(org.Hs.eg.db)
    ##GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
    ##genes = head(as.character(unlist(as.list(org.Hs.egSYMBOL))),1000)
    genes = unique(as.character(ngs$genes$gene_name))
    genes <- toupper(genes)  ## handle mouse genes...
    G <- G[rownames(G) %in% genes,]
    dim(G)

    ##-----------------------------------------------------------
    ## Filter gene sets
    ##-----------------------------------------------------------

    ## filter gene sets on size
    cat("Filtering gene sets on size...\n")
    gmt.size = colSums(G!=0)
    summary(gmt.size)
    size.ok <- (gmt.size >= 15 & gmt.size <= 1000 )
    G <- G[, which(size.ok)]
    dim(G)   
    table(sub(":.*","",colnames(G)))
    
    ##-----------------------------------------------------------
    ## create the full GENE matrix (always collapsed by gene)
    ##-----------------------------------------------------------

    single.omics <- !any(grepl("\\[",rownames(ngs$counts)))
    single.omics
    if(single.omics) {
        ## normalized matrix
        X <- ngs$X
    } else {
        data.type <- gsub("\\[|\\].*","",rownames(ngs$counts))
        jj <- which(data.type %in% c("gx","mrna"))
        length(jj)
        if(length(jj)==0) {
            stop("FATAL. could not find gx/mrna values.")
        }
        X <- ngs$X[jj,]        
    }
    
    ##-----------------------------------------------------------
    ## create the GENESETxGENE matrix
    ##-----------------------------------------------------------
    cat("Matching gene set matrix...\n")
    gg <- toupper(rownames(X)) ## accomodate for mouse...
    ii <- intersect(gg,rownames(G))
    G <- G[ii,]
    xx <- setdiff(gg,rownames(G))
    matX <- Matrix(0, nrow=length(xx), ncol=ncol(G), sparse=TRUE)
    rownames(matX) <- xx
    colnames(matX) <- colnames(G)
    G <- rbind(G, matX)
    G <- G[match(gg,rownames(G)),]
    rownames(G) <- rownames(X) ## original name (e.g. mouse)
    dim(G)
    dim(X)
    
    ##-----------------------------------------------------------
    ## Prioritize gene sets by fast rank-correlation
    ##-----------------------------------------------------------    
    if(is.null(max.features)) max.features <- 20000
    if(max.features < 0) max.features <- 20000
    max.features
    
    if(max.features > 0) {
        cat("Reducing gene set matrix...\n")
        require(limma)
        ## Reduce gene sets by selecting top varying genesets. We use the
        ## very fast sparse rank-correlation for approximate single sample
        ## geneset activation.
        cX <- X - rowMeans(X, na.rm=TRUE)  ## center!
        gsetX = qlcMatrix::corSparse( G[,], apply( cX[,],2,rank) )
        ## gsetX = limma::normalizeQuantiles(gsetX) ##???
        ##grp <- ngs$samples$group
        grp <- ngs$model.parameters$group
        if(!is.null(grp)) {
            gsetX.bygroup <- t(apply(gsetX,1,function(x) tapply(x,grp,mean)))
            sdx <- apply(gsetX.bygroup,1,sd)
        } else {
            sdx <- apply(gsetX,1,sd)
        }
        names(sdx) <- colnames(G)
        jj = head(order(-sdx), max.features) 
        must.include <- "hallmark|kegg|^go|^celltype"
        jj = unique( c(jj, grep(must.include,colnames(G),ignore.case=TRUE)))
        jj = jj[order(colnames(G)[jj])]
        length(jj)
        G = G[,jj,drop=FALSE]
    }
    dim(G)
    
    
    ##-----------------------------------------------------------
    ## get design and contrast matrix
    ##-----------------------------------------------------------
    design = ngs$model.parameters$design
    contr.matrix = ngs$model.parameters$contr.matrix
    ##contr.matrix
    ##exp.matrix = (design %*% contr.matrix)
    
    ALL.GSET.METHODS=c("fisher","ssgsea","gsva", "spearman", "camera", "fry",
                       "fgsea","gsea.permPH","gsea.permGS","gseaPR")
    test.methods
    
    ##-----------------------------------------------------------
    ## Run methods
    ##-----------------------------------------------------------
    cat(">>> Testing gene sets with methods:",test.methods,"\n")
    test.methods

    ## convert to gene list
    gmt <- lapply(apply(G!=0, 2, which),names)    
    Y <- ngs$samples
    gc()
    
    gset.meta = gset.fitContrastsWithAllMethods(
        gmt = gmt, X = X, Y = Y, G = G,
        design = design, ## genes=GENES,
        contr.matrix = contr.matrix, methods = test.methods,
        mc.threads = 1, mc.cores = NULL, batch.correct = TRUE
    )
    
    rownames(gset.meta$timings) <- paste("[test.genesets]",rownames(gset.meta$timings))
    print(gset.meta$timings)

    ngs$timings <- rbind(ngs$timings, gset.meta$timings)
    ngs$gset.meta <- gset.meta

    names(ngs$gset.meta$matrices)
    ##ngs$gsetX = ngs$gset.meta$matrices[["fc"]]  ## META or average FC??!
    ngs$gsetX = ngs$gset.meta$matrices[["meta"]]  ## META or average FC??!
    ngs$GMT <- G[,rownames(ngs$gsetX)]
    
    ##-----------------------------------------------------------------------
    ##------------------------ clean up -------------------------------------
    ##-----------------------------------------------------------------------
    
    ## remove large outputs... (uncomment if needed)
    if(remove.outputs) {
        ngs$gset.meta$outputs <- NULL
    }
    
    remove(X)
    remove(Y)
    remove(G)
    remove(gmt)

    return(ngs)
}



##rm(list=setdiff(ls(),SAVE.PARAMS))

