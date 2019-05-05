
ngs.collapseByGene <- function(ngs) {
    ##sum(duplicated(ngs$genes$gene_name))
    x1 = apply(ngs$counts, 2, function(x) tapply(x, ngs$genes$gene_name, sum))
    x1 <- x1[!(rownames(x1) %in% c(NA,"","NA")),,drop=FALSE]
    ngs$genes = ngs$genes[match(rownames(x1),ngs$genes$gene_name),]
    ngs$counts = x1
    rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
    return(ngs)
}

