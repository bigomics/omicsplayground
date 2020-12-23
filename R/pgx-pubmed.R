##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

require(parallel)
NCORE <- detectCores(all.tests = FALSE, logical = TRUE)/2
NCORE
BLUERED <- colorRampPalette(
    rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0",
          "#92C5DE", "#4393C3", "#2166AC", "#053061")))

if(0) {
    GENERIF.MATRIX <- readRDS(file=file.path(FILES,"geneRIF-matrix.rds"))
    gene="Socs3";context=c("inflamm","cancer")
    gene="Socs3";context="."
}

pmid.getGeneContext <- function(gene, keyword)
{
    library(data.table)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    
    gene1 <- c(gene,sub("([0-9])","-\\1",gene))
    ##gene1 <- paste0(paste0("^",gene1),"|[\\(, ]",gene1,"[\\)-,\\( ]|",gene1,"$")
    gene1 <- paste0("^",gene1,"$|^",gene1,"[-]")
    gene1 <- paste(gene1,collapse="|")
    gene1
    
    if(gene %in% keys(org.Hs.egALIAS2EG)) {
        gname <- get(get(gene, org.Hs.egALIAS2EG),org.Hs.egGENENAME)
        gname <- gsub("[, -]",".",gname)
        gene1 <- paste0(gene1,"|",gname)
    } else if(gene %in% keys(org.Mm.egALIAS2EG)) {
        gname <- get(get(gene, org.Mm.egALIAS2EG),org.Mm.egGENENAME)
        gname <- gsub("[, -]",".",gname)
        gene1 <- paste0(gene1,"|",gname)
    }
    gene1
    
    rif.words <- colnames(GENERIF.MATRIX)
    ii <- grepl(gene1, rif.words, ignore.case=TRUE)
    match0 <- rowSums(GENERIF.MATRIX[,ii,drop=FALSE]) >0
    match1 <- rep(1, length(match0))
    if(length(keyword)>0) {
        i=1
        for(i in 1:length(keyword)) {
            jj <- grepl(keyword[i], rif.words, ignore.case=TRUE)
            match1 <- match1 & (rowSums(GENERIF.MATRIX[,jj,drop=FALSE])>0)
        }
    }
    sel <- ((match0 * match1) >0)
    table(sel)
    rif.hits <- rownames(GENERIF.MATRIX)[sel]
    rif.hits <- rif.hits[!duplicated(rif.hits)]
    ## rif.hits

    ## calculate P-value for this keyword
    A <- table(gene=match0, keyword=match1)
    A
    pv <- NA
    if(nrow(A)==2 && ncol(A)==2) {
        pv <- fisher.test(A, alternative="greater")$p.value
    }
    pv    
    
    context1 <- NULL
    if(1) {
        match2 <- ((match0 * match1)>0)
        m0 <- colSums(match2 * (GENERIF.MATRIX!=0))
        m1 <- colSums(GENERIF.MATRIX!=0)
        pp <- corpora::fisher.pval(m0, sum(match2)+1, m1, nrow(GENERIF.MATRIX)+1, alternative="greater")
        pp <- sort(pp)
        qq <- p.adjust(pp)
        qq <- sort(qq)
        context1 <- head(qq[qq < 1],100)
        head(context1,20)
    }

    list(rifs=rif.hits, table=A, p.value=pv, context=context1)
}

pmid.getPubMedContext <- function(gene, context) {
    ##install.packages("RISmed")
    require(RISmed)
    require(org.Hs.eg.db)
    res <- EUtilsSummary(
        paste0(gene,"[sym] AND ",context),
        type="esearch", db="pubmed", datetype='pdat',
        mindate=2000, maxdate=2099, retmax=1000)
    QueryCount(res)
    if(QueryCount(res)==0) {
        return(NULL)
    }
    gene1 <- c(gene,sub("([0-9])","-\\1",gene),sub("([0-9])","[ ]\\1",gene))
    gene1
    gene1 <- paste(gene1,collapse="|")
    gene1
    ##gname <- get(get("Hmox1", org.Mm.egALIAS2EG),org.Mm.egGENENAME)
    gname <- get(get(toupper(gene), org.Hs.egALIAS2EG),org.Hs.egGENENAME)
    gene2 <- paste0(gene1,"|",gsub("[ -]",".",gname))
    extractRIF <- function(a) {
        s <- strsplit(a,split="[.;:]")[[1]]
        hit <- grepl(gene2,s,ignore.case=TRUE) & grepl(context,s,ignore.case=TRUE)
        if(!any(hit)) return(NULL)
        s <- paste(s[hit], collapse=".")
        sub("^[ ]*","",s)
    }        
    fetch <- EUtilsGet(res)
    tt <- ArticleTitle(fetch)
    aa <- AbstractText(fetch)
    tt2 <- paste(tt,".",aa)
    pp <- PMID(fetch)
    rif <- lapply(tt2, extractRIF)
    rif[sapply(rif,is.null)] <- NA
    rif <- unlist(rif)
    rif[!is.na(rif)] <- paste0(rif," (PMID:",pp,")")[!is.na(rif)]
    list(pmid=pp, title=tt, rif=rif)
}


pmid.buildMatrix <- function() {
    require(org.Hs.eg.db)
    pmid   <- as.list(org.Hs.egPMID2EG)
    symbol <- as.list(org.Hs.egSYMBOL)
    eg <- names(symbol)
    symbol <- sapply(symbol,"[",1)
    names(symbol) <- as.character(eg)
    ngene <- sapply(pmid,length)
    pmid <- pmid[which(ngene <=10)]
    
    ## collapse duplicates
    pmid.gg <- sapply(pmid,paste,collapse=",")
    idx <- tapply(names(pmid), pmid.gg, function(x) x)
    pmid <- pmid[sapply(idx,"[",1)]
    idx <- lapply(idx, function(x) paste0("PMID:",x))
    names(pmid) <- sapply(idx,paste,collapse=",")
    
    ## build PMID2SYMBOL matrix
    require(parallel)
    require(Matrix)
    require(qlcMatrix)
    idx0 <- mclapply(1:length(pmid), function(i) cbind(i, which(eg %in% pmid[[i]])),
                     mc.cores=NCORE)
    idx <- do.call(rbind,idx0)
    dim(idx)
    P <- sparseMatrix( i=idx[,1], j=idx[,2], x=rep(1,nrow(idx)),
                      dims=c(length(idx0), length(eg)) )
    rownames(P) <- names(pmid)
    colnames(P) <- symbol
    P <- P[,which(Matrix::colSums(P)>0)]
    dim(P)
    return(P)
}

pmid.buildGraph <- function(P) {
    require(org.Hs.eg.db)
    library(igraph)
    ##P <- readRDS(file="PMID2SYMBOL_sparsematrix.rds")
    dim(P)
    P <- P[which(Matrix::rowSums(P) <= 10),]
    P <- P[which(Matrix::rowSums(P) >= 2),]
    P <- P[,which(Matrix::colSums(P)>0)]
    ##P <- head(P,4000)
    ##P <- P[sample(nrow(P),5000),]
    P[1:10,1:10]
    dim(P)
    
    ## create graph from overlap
    M <- P[,] %*% t(P[,])
    dim(M)
    diag(M) <- 0
    object.size(M)
    gr <- graph_from_adjacency_matrix(M, mode="undirected",
                                      diag=FALSE, weighted=TRUE)
    V(gr)$name <- rownames(M)
    gr <- subgraph.edges(gr, which(E(gr)$weight>0))
    ##gr
    ##saveRDS(gr, file="PMID2SYMBOL_xgraph_01.rds")
    
    P <- P[V(gr)$name,]
    pmids <- mclapply(V(gr)$name,function(x) gsub("PMID:","",strsplit(x,split=",")[[1]]))
    nref <- sapply(pmids,length)
    ##vgenes <- apply(P[V(gr)$name,],1,function(x) paste(names(which(x!=0)),collapse=","))
    vgenes <- mclapply(1:nrow(P),function(i) names(which(P[i,]!=0)), mc.cores=NCORE)
    vgenes2 <- unlist(sapply(vgenes,paste,collapse=","))
    V(gr)$size <- nref
    V(gr)$genes <- vgenes
    V(gr)$pmid  <- pmids
    return(gr)    
}


pmid.annotateEdges <- function(gr) {
    require(igraph)
    ee <- get.edges(gr, E(gr))
    dim(ee)
    g1 <- V(gr)[ee[,1]]$genes
    g2 <- V(gr)[ee[,2]]$genes
    shared.genes <- mapply(intersect, g1, g2)
    nshared <- sapply(shared.genes, length)
    E(gr)$genes  <- shared.genes
    E(gr)$weight <- nshared
    gr
}

pmid.extractGene <- function(gr, gene, nmin=3) {
    require(igraph)
    jj <- c()
    for(g in gene) {
        j1 <- which(sapply(V(gr)$genes, function(s) (gene %in% s)))
        ##jj <- c(jj, grep(g, V(gr)$genes))
        jj <- c(jj, j1)
    }
    jj <- unique(jj)
    ngene1 <- sapply(V(gr)$genes[jj],length)
    gr1 <- induced_subgraph(gr, jj)
    if(verbose>0) cat("annotating edges...\n")
    gr1 <- pmid.annotateEdges(gr1)
    ##nshared <- sapply(E(gr1)$genes,length)
    nshared <- unlist(mclapply(E(gr1)$genes,length,mc.cores=NCORE))
    ee <- which(!( nshared==1 & sapply(E(gr1)$genes,"[",1) %in% gene))
    gr1 <- subgraph.edges(gr1, ee, delete.vertices=TRUE)
    cmp <- components(gr1)
    jj <- which(cmp$membership %in% which(cmp$csize >= nmin))
    gr1 <- induced_subgraph(gr1, jj)
    gr1
}

pubmedlink <- function(s) {
    paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed/",s,
           "' target='_blank'>PMID:",s,"</a>")
}
