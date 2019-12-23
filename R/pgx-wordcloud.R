##---------------------------------------------------------------
##------------- Functions for WordCloud ------------------------
##---------------------------------------------------------------

pgx.calculateWordFreq <- function(ngs, progress=NULL, pg.unit=1) {
    
    if(!is.null(progress)) progress$set(message = "WordCloud", value = 0)
    
    ## get gset meta foldchange-matrix
    S <- sapply( ngs$gset.meta$meta, function(x) x$meta.fx)
    rownames(S) <- rownames(ngs$gset.meta$meta[[1]])
    ##S <- S[order(-apply(S,1,sd)),]
    S <- S[order(-rowMeans(S**2)),]
    
    ## exclude down, GSE gene sets??????
    S <- S[grep("dn|down|^gse",rownames(S),ignore.case=TRUE,invert=TRUE),]
    
    if(!is.null(progress)) progress$inc(0.2*pg.unit, detail="calculating word frequencies")
    
    ## Determine top most frequent terms
    sname <- gsub(".*:","",tolower(rownames(S)))
    sname <- gsub("b.cell","bcell",sname)
    sname <- gsub("t.cell","tcell",sname)
    words <- strsplit(sname,split="[-_ ]")    
    names(words) <- rownames(S)
    terms <- names(sort(-table(unlist(words))))
    stopwords = strsplit(
        "with int like strand tconv pid lee and or mouse dn small big human homo sapiens mus musculus drug early late hsa gse culture in the a g for line r up down events anti large targets tissue vitro process cells ctrl regulation processing common pathway months days white pre post process all mice from", split=" ")[[1]]
    terms <- terms[which(!terms %in% stopwords)]
    terms <- terms[sapply(terms,nchar)>2]
    terms <- grep("[0-9]|^\\(",terms,invert=TRUE,value=TRUE)        
    length(terms)
    terms <- head(terms,1000)
    
    ## Calculate incidence matrix
    words2 <- lapply(words, function(w) intersect(w, terms))
    words2 <- words2[sapply(words2,length)>0]
    idx <- lapply(1:length(words2), function(i) cbind(i,match(words2[[i]],terms)))
    idx <- do.call(rbind, idx)
    library(qlcMatrix)
    W <- sparseMatrix(idx[,1], idx[,2], x=1)
    dim(W)
    rownames(W) = names(words2)
    colnames(W) = terms
    
    ## filter on minimal size and maximum ratio
    nn <- Matrix::colSums(W,na.rm=TRUE)
    nr <- nn / nrow(W)
    W <- W[,which(nn >= 10 & nr <= 0.10)]
    dim(W)
    
    ## align geneset expression matrix
    S <- S[rownames(W),,drop=FALSE]
    
    if(!is.null(progress)) progress$inc(0.3*pg.unit, detail="computing GSEA")        
    
    ## compute for average contrast
    require(fgsea)
    rms.FC <- Matrix::rowMeans(S**2)**0.5
    rms.FC <- rms.FC + 0.01*rnorm(length(rms.FC))
    gmt <- apply(W,2,function(x) names(which(x!=0)))
    res <- fgsea( gmt, rms.FC, nperm=1000 )
    res$leadingEdge <- sapply(res$leadingEdge,paste,collapse="//")
    ## res$leadingEdge <- NULL
    colnames(res)[1] <- "word"
    
    ## --------- only significant and positive
    ##res <- res[(res$padj < 0.20 & res$NES>0),]
    res <- res[(res$padj < 1 & res$NES>0),]
    res <- res[order(-abs(res$NES)),]
    dim(res)
    
    ## now compute significant terms for all contrasts
    all.gsea <- list()
    for(i in 1:ncol(S)) {
        fc <- as.vector(S[,i])
        names(fc) <- rownames(W)
        fc <- fc + 0.01*rnorm(length(fc))
        gmt1 <- gmt[as.character(res$word)]
        res1 <- fgsea( gmt1, fc, nperm=1000 )
        res1$leadingEdge <- sapply(res1$leadingEdge,paste,collapse="//")
        ## res$leadingEdge <- NULL
            colnames(res1)[1] <- "word"
        all.gsea[[colnames(S)[i]]] <- res1
    }
    all.gsea[["rms.FC"]] <- res
    
    if(!is.null(progress)) progress$inc(0.25*pg.unit, detail="clustering")
    
    require(Rtsne)
    library(umap)
    pos1 = Rtsne(as.matrix(t(W)),perplexity=10,check_duplicates=FALSE)$Y
    pos2 = umap(as.matrix(t(W)))$layout
    rownames(pos1) = rownames(pos2) = colnames(W)
    colnames(pos1) = colnames(pos2) = c("x","y")
    pos1 = pos1[res$word,]
    pos2 = pos2[res$word,]
    
    all.res = list(gsea=all.gsea, S=S, W=W, tsne=pos1, umap=pos2)
    all.res
}
