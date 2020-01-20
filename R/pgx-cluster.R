

##reduce.sd=1000;reduce.pca=100;methods=c("pca","tsne","umap");dims=c(2,3)
pgx.clusterBigMatrix <- function(X, methods=c("pca","tsne","umap"), dims=c(2,3),
                                 reduce.sd = 1000, reduce.pca = 100 )
{
    require(irlba)
    require(Rtsne)
    require(umap)

    methods <- intersect(methods, c("pca","tsne","umap"))
    if(length(methods)==0) methods <- "pca"
    
    ## Pre-reduce dimensions by SD
    if(reduce.sd>0) {
        cat("pre-reducing to ",reduce.sd," max SD features...\n")                
        sdx <- apply(X,1,sd,na.rm=TRUE)
        rX <- X[head(order(-sdx),reduce.sd),]
    } else {
        rX <- X
    }
        
    ## scale and augment if few samples
    rX[is.na(rX)] <- 0
    rX <- rX - rowMeans(rX) ## do??
    rX <- scale(rX) ## columm normalization
    if(ncol(rX)<=6) rX <- cbind(rX,rX,rX,rX,rX,rX)
    rX <- rX + 1e-2*matrix(rnorm(length(rX)),nrow(rX),ncol(rX))
    dim(rX)
    
    ## Further pre-reduce dimensions using SVD
    res.svd <- NULL
    if(reduce.pca>0) {
        reduce.pca <- max(3,min(c(reduce.pca,dim(X)-1)))
        reduce.pca
        cat("pre-reducing to ",reduce.pca," PCA dimenstions...\n")        
        res.svd <- irlba(rX, nv=reduce.pca)
        rX <- t(res.svd$v) * res.svd$d
    }
    dim(rX)
    
    all.pos <- list()

    if("pca" %in% methods && 2 %in% dims) {
        cat("calculating PCA 2D/3D...\n")
        require(irlba)
        if(is.null(res.svd)) {
            res.svd <- irlba(rX, nv=3)
        }
        pos <- res.svd$v[,1:2]
        pos <- pos[1:ncol(X),]  ## if augmented
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("PC-",c("x","y"))                
        all.pos[["pca2d"]] <- pos
    }

    if("pca" %in% methods && 3 %in% dims) {
        if(is.null(res.svd)) {
            res.svd <- irlba(rX, nv=3)
        }
        pos <- res.svd$v[,1:3]        
        pos <- pos[1:ncol(X),]  ## if augmented
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("PC-",c("x","y","z"))                
        all.pos[["pca3d"]] <- pos
    }
        
    if("tsne" %in% methods && 2 %in% dims) {
        cat("calculating t-SNE 2D...\n")
        require(Rtsne)
        perplexity <- pmax(min(ncol(X)/4,30),2)
        perplexity        
        pos <- Rtsne( t(rX[,]), dims = 2,
                    is_distance = FALSE, check_duplicates=FALSE,                     
                    perplexity = perplexity, num_threads=0)$Y
        pos <- pos[1:ncol(X),]  ## if augmented
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("tSNE-",c("x","y"))                
        all.pos[["tsne2d"]] <- pos
    }
    
    if("tsne" %in% methods && 3 %in% dims) {
        cat("calculating t-SNE 3D...\n")
        require(Rtsne)
        perplexity <- pmax(min(ncol(X)/4,30),2)
        perplexity        
        pos <- Rtsne( t(rX[,]), dims = 3,
                    is_distance = FALSE, check_duplicates=FALSE,
                    perplexity = perplexity, num_threads=0)$Y
        pos <- pos[1:ncol(X),]  ## if augmented        
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("tSNE-",c("x","y","z"))        
        all.pos[["tsne3d"]] <- pos
    }
    
    if("umap" %in% methods && 2 %in% dims) {
        cat("calculating UMAP 2D...\n")
        require(umap)
        pos <- umap( t(rX[,]) )$layout
        ## rownames(pos) <- colnames(X)[1:nrow(pos)]
        pos <- pos[1:ncol(X),]  ## if augmented
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("UMAP-",c("x","y"))        
        all.pos[["umap2d"]] <- pos
    }
    
    if("umap" %in% methods && 3 %in% dims) {
        cat("calculating UMAP 3D...\n")
        require(umap)
        custom.config <- umap.defaults
        custom.config$n_components = 3        
        pos <- umap( t(rX[,]), custom.config )$layout
        dim(pos)
        pos <- pos[1:ncol(X),]  ## if augmented
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("UMAP-",c("x","y","z"))
        all.pos[["umap3d"]] <- pos
    }

    return(all.pos)
}


##skipifexists=0;perplexity=NULL;sv.rank=-1;prefix="C";kclust=1;ntop=1000;fromX=FALSE;prior.counts=1;mean.center=TRUE;method="tsne";determine.clusters=1;dims=c(2,3);find.clusters=TRUE;row.center=TRUE;row.scale=FALSE;clust.detect="louvain"
pgx.clusterSamples <- function(ngs, skipifexists=FALSE, perplexity=NULL,
                               ntop=1000, sv.rank=-1, prefix="C", 
                               fromX=FALSE, is.logx=FALSE,
                               kclust=1, prior.counts=NULL, 
                               dims=c(2,3), find.clusters=TRUE,
                               clust.detect = c("louvain","hclust"),
                               row.center=TRUE, row.scale=FALSE,
                               method=c("tsne","umap","pca") )
{

    clust.detect <- clust.detect[1]
    
    sX <- ngs$counts
    if(fromX) sX <- 2**ngs$X
    res <- NULL
    res <- pgx.clusterSamplesFromMatrix(
        sX, perplexity=perplexity, is.logx=FALSE,
        ntop=ntop, sv.rank=sv.rank, prefix=prefix,         
        kclust=kclust, prior.counts=prior.counts, 
        dims=dims, find.clusters=find.clusters,
        clust.detect = clust.detect,
        row.center=row.center, row.scale=row.scale,
        method=method)

    if(!is.null(res$pos2d)) ngs$tsne2d <- res$pos2d
    if(!is.null(res$pos3d)) ngs$tsne3d <- res$pos3d
    if(!is.null(res$idx)) {
        if(class(ngs$samples)=="data.frame") {
            ngs$samples$cluster <- as.character(res$idx)
        } else {
            ## matrix??
            if("cluster" %in% colnames(ngs$samples)) {
                ngs$samples[,"cluster"] <- as.character(res$idx)
            } else {
                ngs$samples <- cbind(ngs$samples, cluster=as.character(res$idx))
            }
        }        
    }
    
    return(ngs)
}

pgx.clusterSamplesFromMatrix <- function(counts, perplexity=NULL,
                                         ntop=1000, sv.rank=-1, prefix="C", 
                                         fromX=FALSE, is.logx=FALSE, 
                                         prior.counts=NULL, dims=c(2,3),
                                         row.center=TRUE, row.scale=FALSE,
                                         find.clusters=TRUE, kclust=1,
                                         clust.detect = c("louvain","hclust"),
                                         method=c("tsne","umap","pca") )
{
    require(Rtsne)
    require(irlba)
    ##set.seed(0)
    method <- method[1]
    clust.detect <- clust.detect[1]
    
    sX <- counts
    if(is.logx) sX <- 2**sX
    
    if(is.null(prior.counts)) {
        qq <- quantile(as.vector(sX[sX>0]),probs=0.50) ## at 50%
        qq
        prior.counts <- qq
    }
    sX <- log2(prior.counts + sX)
    sX <- limma::normalizeQuantiles(sX)  ## in linear space
    if(row.center) sX <- sX - rowMeans(sX,na.rm=TRUE)
    sX = head( sX[order(-apply(sX,1,sd)),], ntop)
    ## sX = t(scale(t(sX),scale=TRUE))  ## really? or just centering?
    ##sX = t(scale(t(sX),scale=FALSE))  ## really? or just centering?
    if(row.scale) sX <- (sX / apply(sX,1,sd,na.rm=TRUE))
    
    dim(sX)
    ## some randomization is sometimes necessary if the data is 'too
    ## clean' and clusters become lines..
    ##sX = sX + 0.001*matrix(rnorm(length(sX)),nrow(sX),ncol(sX))

    ## ------------ find t-SNE clusters
    if(is.null(perplexity)) {
        perplexity = max(1,min(30, round((ncol(sX)-1)/4)))
        perplexity
    }

    ##sv.rank=20
    if(sv.rank > 0) {
        cat("performing tSNE on reduced PCA k=",sv.rank,"\n")
        svd <- irlba(sX, nv=sv.rank)
        sv <- svd$v %*% diag(svd$d[1:ncol(svd$v)])
        rownames(sv) <- colnames(sX)
        sX <- t(sv)
    }

    pos2=pos3=NULL
    if(method=="umap") {
        require(umap)
        if(2 %in% dims) {
            pos2 = umap(
                t(sX),
                n_neighbours = perplexity,
                n_components = 2,
                metric = "euclidean"
            )$layout
            colnames(pos2) <- c("umap_1","umap_2")
        }
        if(3 %in% dims) {
            pos3 = umap(
                t(sX),
                n_neighbours = perplexity,
                n_components = 3,
                metric = "euclidean"
            )$layout
            colnames(pos3) <- c("umap_1","umap_2","umap_3")
        }
    } else if(method=="tsne") {
        if(2 %in% dims) {
            pos2 = Rtsne( t(sX), dim=2, perplexity=perplexity,
                         check_duplicates=FALSE, num_threads=99)$Y
            colnames(pos2) <- c("tnse_1","tnse_2")
        }
        if(3 %in% dims) {
            pos3 = Rtsne( t(sX), dim=3, perplexity=perplexity,
                         check_duplicates=FALSE, num_threads=99)$Y
            colnames(pos3) <- c("tnse_1","tnse_2","tnse_3")
        }
    } else if(method=="pca") {
        sv.rank
        svd <- irlba(sX, nv=3)
        if(2 %in% dims) {
            pos2 = svd$v[,1:2]
            colnames(pos2) <- c("pca_1","pca_2")
        }
        if(3 %in% dims) {
            pos3 = svd$v[,1:3]
            colnames(pos3) <- c("pca_1","pca_2","pca_3")
        }
    }
    if(!is.null(pos2)) {
        rownames(pos2) = colnames(sX)
    }
    if(!is.null(pos3)) {
        rownames(pos3) = colnames(sX)
    }
    
    ## ------------ find t-SNE clusters from graph
    idx = NULL
    if(find.clusters && clust.detect=="hclust") {
        cat("Finding clusters using hclust...\n")
        hc <- hclust(dist(t(sX)))
        idx2 <- cutree(hc,2)
        idx3 <- cutree(hc,3)
        idx4 <- cutree(hc,4)
        idx5 <- cutree(hc,5)
        idx <- cutree(hc,kclust)

    } else if(find.clusters && clust.detect=="louvain") {
        cat("Finding clusters using Louvain...\n")
        require(igraph)
        if(!is.null(pos2)) pos <- pos2
        if(!is.null(pos3)) pos <- pos3
        dist = as.dist(dist(scale(pos))) ## use 3D distance??
        
        ##dist = dist + 0.1*mean(dist)
        gr = graph_from_adjacency_matrix(1.0/dist, diag=FALSE, mode="undirected")
        ## should we iteratively cluster???
        hc <- hclustGraph(gr)  ##
        dim(hc)
        idx <- hc[,min(kclust,ncol(hc))]
        ##idx <- cluster_louvain(gr)$membership
        table(idx)
        
        ## ------------ zap small clusters to "0"
        sort(table(idx))
        min.size <- pmax(3, 0.01*length(idx))
        min.size
        small.clusters <- names(which(table(idx) < min.size))
        idx[ which(idx %in% small.clusters)] <- "0"
        sort(table(idx))        
    }

    ## rename levels with largest cluster first
    idx <- factor(idx, levels=names(sort(-table(idx))))
    levels(idx) <- paste0(prefix,1:length(levels(idx)))
    table(idx)
    cat("Found",length(unique(idx)),"clusters...\n")

    res <- list(pos2d=pos2, pos3d=pos3, idx=idx)
    return(res)
}
