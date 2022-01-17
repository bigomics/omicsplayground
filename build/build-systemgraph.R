##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##install.packages("bigalgebra")
##install.packages("bigmemory")
##install.packages("biganalytics")
##install.packages("bigpca")

source("./R/gx-util.r")

library(Rtsne)
library(qlcMatrix)
library(irlba)
library(sparsesvd)

require(org.Hs.eg.db)
dir("~")

##----------------------------------------------------------------------
## From geneset GMT file, create large sparse 0/1 matrix.
##----------------------------------------------------------------------

load(file="./lib/gmt-all.rda",verbose=1)
G <- build.createSparseGenesetMatrix(gmt.all)
dim(G)
save(G, file="gset-sparseG-XL2.rda")

load(file="gset-sparseG-XL2.rda", verbose=1)
dim(G)
max(G)


build.computeGeneGenesetTSNE <- function(X, G) {
    ##
    ## Compute t-SNE for all genes and gene sets
    ##
    ##
    ##
    
    ##load(file="./files/gset-sparseG-XL.rda",verbose=1)
    ##head(G)[,1:5]
    dim(G)
    cat("dim(G)=",dim(G),"\n")

    ## compute SD (both directions)
    sd1 <- apply(G, 1, sd, na.rm=TRUE)
    sd2 <- apply(G, 2, sd, na.rm=TRUE)
        
    cat("computing SVD...\n")
    require(sparsesvd)
    ##G1 = 1*(G!=0)
    ##system.time( res <- irlba(G[1:1000,], nv=200) )
    ##system.time( res <- big.PCA(G[,], nv=20) )
    jj <- head(order(-sd1),1000)
    system.time( res <- sparsesvd(G[,jj], rank=100) )
    G2 <- res$u %*% diag(res$d)
    dim(G2)
    cat("computing t-SNE for gene sets...\n")
    system.time( tsne <- Rtsne( as.matrix(G2), check_duplicates=FALSE,
                               num_threads=4))
    ##plot(tsne$Y, pch=20, cex=0.5)

    rownames(tsne$Y) <- rownames(G2)
    dim(tsne$Y)
    ## write.csv( tsne$Y, file="tsne-all-genesets-XL-r8k.csv")
    write.csv( tsne$Y, file="gset-sparseG-XL2-tsneGS.csv")
    remove(G2)

    ## -------------------- compute t-SNE ----------------------------

    cat("computing t-SNE for genes...\n")
    G2 <- res$v %*% diag(res$d)
    dim(G2)
    tsne <- Rtsne.multicore( as.matrix(G2), check_duplicates=FALSE, num_threads=24)
    dim(tsne$Y)
    rownames(tsne$Y) <- colnames(G)
    dim(tsne$Y)
    ## write.csv( tsne$Y, file="tsne-all-genes-XL-r8k.csv")
    write.csv( tsne$Y, file="gset-sparseG-XL2-tsneGenes.csv")
    
    ##pos=read.csv(file="tsne-all-genes-XL.csv",row.names=1)
    ##dim(pos)

}


if(0) {

    ##
    ## Create a generic multilayer SNN igraph object from t-SNE
    ## positions for gene and gene set layer.
    ##
    ##
    ##

    library(scran)
    library(igraph)
    require(gplots)
    
    G <- readRDS(file="./lib/gset-sparseG-XL.rds")
    tsne_gsets <- read.csv(file="./lib/tsne-all-genesets-XL.csv",row.names=1)
    tsne_genes <- read.csv(file="./lib/tsne-all-genes-XL.csv",row.names=1)
    ##tsne_gsets <- read.csv(file="tsne-all-genesets-wXL.csv",row.names=1)
    ##tsne_genes <- read.csv(file="tsne-all-genes-wXL.csv",row.names=1)
    dim(tsne_gsets)
    dim(tsne_genes)

    uscale <- function(x) (x-min(x))/(1e-8+max(x)-min(x))-0.5

    ## conform
    ##pos1 <- tsne_gsets[sample(nrow(tsne_gsets),1234),]
    ##pos2 <- tsne_genes[sample(nrow(tsne_genes),1000),]
    pos1 <- tsne_gsets
    pos2 <- tsne_genes
    pos1 <- as.matrix(pos1)
    pos2 <- as.matrix(pos2)
    ss <- intersect( rownames(pos1), rownames(G))
    gg <- intersect( rownames(pos2), colnames(G))
    G1 <- G[ss,gg]
    pos1 <- pos1[ss,]
    pos2 <- pos2[gg,]

    ## make graph object
    KNN=4
    KNN=10
    KNN=20
    ##KNN=100
    g1 <- buildSNNGraph(t(pos1), k=KNN)
    g2 <- buildSNNGraph(t(pos2), k=KNN)
    V(g1)$name <- paste0("{geneset}",rownames(pos1))
    V(g2)$name <- paste0("{gene}",rownames(pos2))

    if(0) {
        V(g2)$label <- NA
        g2$layout <- as.matrix(pos2)
        plot(g2, vertex.size=0.3, edge.width=0.3, layout=pos2)

    }

    ## combine layers, add interlayer edges
    gr <- graph.union(g1, g2)
    gr
    head(V(gr)$name)
    ee <- Matrix::which(G1!=0, arr.ind=TRUE)
    ee.wt <- G1[ee]
    vname <- sub(".*}","",V(gr)$name)
    head(vname)
    j0 <- match( rownames(G1)[ee[,1]], vname)
    j1 <- match( colnames(G1)[ee[,2]], vname)
    sum(is.na(j0))
    sum(is.na(j1))
    gr <- add_edges(gr, rbind(j0,j1), weight=ee.wt )

    edge_attr_names(gr)
    ww <- cbind(E(gr)$weight, E(gr)$weight_1, E(gr)$weight_2)
    ##ww1 <- edge.attributes(gr, index=E(gr))
    E(gr)$weight <- apply(ww,1,max,na.rm=TRUE)
    gr <- remove.edge.attribute(gr, "weight_1")
    gr <- remove.edge.attribute(gr, "weight_2")
    V(gr)$level <- c("gene","geneset")[1 + 1*grepl("geneset\\}",V(gr)$name) ]
    table( V(gr)$level )

    ## match t-SNE positions to new graph
    pos12 <- rbind( cbind(pos1,V3=2), cbind(pos2,V3=1))
    jj <- match(vname, c(rownames(pos1), rownames(pos2)))
    pos3 <- as.matrix(pos12[jj,])
    rownames(pos3) <- V(gr)$name
    gr$layout <- pos3

    ##save(gr, file="pgx-graph-geneXgset-XL.rda")
    saveRDS(gr, file="./files/pgx-graph-geneXgset-XL-snn20.rds")
    ##g1 <- loadRDS(file="pgx-graph-geneXgset.rds")

    if(0) {
        library(igraph)
        gr1 <- readRDS(file="./files/pgx-graph-geneXgset-XL.rds")
        gr2 <- readRDS(file="./files/pgx-graph-geneXgset-XL-snn20.rds")
        gr3 <- readRDS(file="./files/pgx-graph-geneXgset-wXL-snn20.rds")
        gr4 <- readRDS(file="./files/pgx-graph-geneXgset.rds")

        table( V(gr3)$level)
        table( V(gr4)$level)

        load("./pgx/geiger2018-arginine-4k.pgx",verbose=1)
        k=3
        names(ngs$gx.meta$meta)
        fx1 <- unclass(ngs$gx.meta$meta[[k]]$fc)[,"trend.limma"]
        fx2 <- unclass(ngs$gset.meta$meta[[k]]$fc)[,"gsva"]
        names(fx1) <- paste0("{gene}",names(fx1))
        names(fx2) <- paste0("{geneset}",names(fx2))
        fx <- c(fx1,fx2)

        source("./R/pgx-graph.R")
        par(mfrow=c(4,1), mar=c(1,1,1,1))
        pgx.plotDualProjection(gr1, fx=fx)
        pgx.plotDualProjection(gr2, fx=fx)
        pgx.plotDualProjection(gr3, fx=fx)
        pgx.plotDualProjection(gr4, fx=fx)




    }

}




build.createSparseGenesetMatrix <- function(gmt.all) {
    
    ## ----------- Get all official gene symbols
    library(org.Hs.eg.db)
    ##library(org.Mm.eg.db)
    symbol <- as.list(org.Hs.egSYMBOL)
    known.symbols <- sort(unique(unlist(symbol)))
    ##mm.symbol <- as.list( org.Mm.egSYMBOL)
    ##mm.symbol <- sort(unique(unlist(mm.symbol)))
    ##tail(mm.symbol,100)
    ##known.symbols <- unique(c(known.symbols, toupper(mm.symbol)))

    ##------------- filter by size
    gmt.size <- sapply(gmt.all,length)
    summary(gmt.size)
    gmt.all <- gmt.all[which(gmt.size >= 15 & gmt.size <= 1000)]
    ##gmt.all <- gmt.all[1:2000]

    ## ------------- filter genes by minimum frequency and chrom
    ##symbol = unlist(as.list(org.Hs.egSYMBOL))
    ##refseq = unlist(as.list(org.Hs.egREFSEQ))
    genes.table <- table(unlist(gmt.all))
    summary(as.integer(genes.table))
    length(gmt.all)
    ##genes <- names(which( genes.table >= 10 & genes.table <= 1000  ))
    genes <- names(which( genes.table >= 10 ))
    genes <- genes[grep("^LOC|RIK$",genes,invert=TRUE)]
    genes <- intersect(genes, known.symbols)
    annot <- ngs.getGeneAnnotation(genes)
    genes <- genes[ !is.na(annot$chr) ]
    length(genes)

    ## Filter genesets with permitted genes (official and min.sharing)
    gmt.all <- mclapply(gmt.all, function(s) intersect(s,genes))
    gmt.size <- sapply(gmt.all,length)
    summary(gmt.size)
    gmt.all <- gmt.all[which(gmt.size >= 15 & gmt.size <= 500)]
    length(gmt.all)

    ## build huge sparsematrix gene x genesets
    require(parallel)
    genes <- sort(genes)
    idx.j <- mclapply(gmt.all[], function(s) match(s,genes))
    idx.i <- lapply(1:length(gmt.all), function(i) rep(i,length(idx.j[[i]])))
    ii <- unlist(idx.i)
    jj <- unlist(idx.j)
    length(ii)
    length(jj)
    G <- sparseMatrix( i=ii, j=jj, x=rep(1,length(ii)),
                      dims = c(length(gmt.all), length(genes)) )
    dim(G)
    
    colnames(G) = genes
    rownames(G) = names(gmt.all)

    return(G)
}
