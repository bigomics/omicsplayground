ngs.detectOrganism <- function(ngs) {
    lowcase.ratio <- mean(grepl("[a-z]",substring(rownames(ngs$counts),2,100)))
    c("human","mouse")[1 + 1*(lowcase.ratio>0.5)]
}

ngs.matchFeatures <- function(ngs, genes) {
    jj <- match(toupper(genes), toupper(ngs$genes$gene_name))
    rownames(ngs$genes)[jj]
}

ngs.collapseByGeneSLOW <- function(ngs) {
    ##sum(duplicated(ngs$genes$gene_name))
    x1 = apply(ngs$counts, 2, function(x) tapply(x, ngs$genes$gene_name, sum))
    x1 <- x1[!(rownames(x1) %in% c(NA,"","NA")),,drop=FALSE]
    ngs$genes = ngs$genes[match(rownames(x1),ngs$genes$gene_name),]
    ngs$counts = x1
    rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
    return(ngs)
}


ngs.collapseByGene <- function(ngs) {
    ##sum(duplicated(ngs$genes$gene_name))
    gene <- as.character(ngs$genes$gene_name)
    p1 <- names(which(table(gene)==1))
    p2 <- names(which(table(gene)>1))
    length(p2)
    if(length(p2)==0) {
        gene <- as.character(ngs$genes$gene_name)
        rownames(ngs$genes)  <- gene
        rownames(ngs$counts) <- gene
        return(ngs)
    }
    j1 <- which(gene %in% p1)
    j2 <- which(gene %in% p2)
    x1 <- ngs$counts[j1,,drop=FALSE]
    rownames(x1) <- ngs$genes$gene_name[j1]
    x2 <- apply(ngs$counts[j2,,drop=FALSE], 2, function(x) tapply(x,gene[j2],sum))
    length(p2)
    if(length(p2)==1) {
        x2 <- matrix(x2,nrow=1)
        rownames(x2) <- p2
        colnames(x2) <- colnames(x1)
    }
    x1 <- rbind(x1,x2)
    x1 <- x1[!(rownames(x1) %in% c(NA,"","NA")),,drop=FALSE]
    x1 <- x1[order(rownames(x1)),,drop=FALSE]
    ngs$genes = ngs$genes[match(rownames(x1),ngs$genes$gene_name),]
    ngs$counts = x1
    rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
    return(ngs)
}

