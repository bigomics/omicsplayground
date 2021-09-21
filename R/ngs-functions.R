##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

ngs.getGeneAnnotation <- function(genes)
{
    require(org.Hs.eg.db)
    require(org.Mm.eg.db)
    
    hs.genes <- unique(unlist(as.list(org.Hs.egSYMBOL)))
    mm.genes <- unique(unlist(as.list(org.Mm.egSYMBOL)))

    ## if multiple genes, take first
    genes0 <- genes ## save
    genes <- sapply(genes, function(s) strsplit(s,split="[;,\\|]")[[1]][1])
    
    is.human <- mean(genes %in% hs.genes) > mean(genes %in% mm.genes)
    is.mouse <- mean(genes %in% hs.genes) < mean(genes %in% mm.genes)
    is.human
    is.mouse

    txlen=SYMBOL=CHRLOC=MAP=NULL
    if(is.human) {
        message("detected organism: human")
        GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
        SYMBOL = unlist(as.list(org.Hs.egSYMBOL))
        names(GENE.TITLE) = SYMBOL
        CHRLOC = as.list(org.Hs.egCHRLOC)
        CHRLOC <- CHRLOC[match(names(SYMBOL),names(CHRLOC))]        
        MAP = as.list(org.Hs.egMAP)
        MAP = MAP[match(names(SYMBOL),names(MAP))]                
        MAP = sapply(MAP,paste,collapse="|")
        names(CHRLOC) = SYMBOL
        names(MAP) = SYMBOL
        
        ##BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")        
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
        tx <- GenomicFeatures::transcriptLengths(txdb)
        TXLEN <- tapply( tx$tx_len, tx$gene_id, mean, na.rm=TRUE)
        TXLEN <- TXLEN[match(names(SYMBOL),names(TXLEN))]
        names(TXLEN) <- SYMBOL        
        Matrix::head(TXLEN)

        ## get gene biotype        
        daf <- GenomicFeatures::transcripts(
                                    EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                    columns = c("gene_name", "gene_biotype"),
                                    return.type="DataFrame")
        GENE.BIOTYPE = daf$gene_biotype
        names(GENE.BIOTYPE) = daf$gene_name
        
    }
    if(is.mouse) {
        message("detected organism: mouse")
        GENE.TITLE = unlist(as.list(org.Mm.egGENENAME))
        SYMBOL = unlist(as.list(org.Mm.egSYMBOL))
        names(GENE.TITLE) = SYMBOL
        CHRLOC = as.list(org.Mm.egCHRLOC)
        CHRLOC <- CHRLOC[match(names(SYMBOL),names(CHRLOC))]
        names(CHRLOC) = SYMBOL
        MAP <- NULL ## no map for mouse???

        ##BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")       
                       
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
        tx <- GenomicFeatures::transcriptLengths(txdb)
        TXLEN <- tapply( tx$tx_len, tx$gene_id, mean, na.rm=TRUE)
        TXLEN <- TXLEN[match(names(SYMBOL),names(TXLEN))]
        names(TXLEN) <- SYMBOL        
        Matrix::head(TXLEN)
        
        daf <- GenomicFeatures::transcripts(
                                    EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                                    columns = c("gene_name", "gene_biotype"),
                                    return.type="DataFrame")
        GENE.BIOTYPE = daf$gene_biotype
        names(GENE.BIOTYPE) = daf$gene_name
        
    }
    remove(txdb)
    remove(tx)
    
    gene_title=gene_biotype=chrom=loc=map=txlen=NULL    
    gene_title <- GENE.TITLE[genes]
    chrloc0 <- CHRLOC[genes]
    loc <- sapply(chrloc0, "[", 1) ## take only first entry !!!
    loc[sapply(loc,is.null)] <- NA
    loc <- abs(as.integer(unlist(loc)))
    chrom <- sapply(chrloc0, function(s) names(s)[1])
    chrom[sapply(chrom,is.null)] <- NA
    chrom <- as.vector(unlist(chrom))
    txlen <- round(TXLEN[genes])
    gene_biotype <- GENE.BIOTYPE[genes]
    
    map <- paste0("chr",chrom)
    map <- sub("chrNA",NA,map)
    if(!is.null(MAP)) map <- MAP[genes]
    
    ## get protein info
    ## fill me    
    annot = data.frame( gene_name = genes,
                       gene_title = gene_title,
                       gene_biotype = gene_biotype,
                       chr=chrom,
                       pos = as.integer(loc), 
                       tx_len = as.integer(txlen),
                       map=map )
    ##genes = apply(genes,2,as.character)
    Matrix::head(annot)
    ##annot[is.na(annot)] <- ""  ## replace NA with empty string
    
    rownames(annot) = genes
    annot
}

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

