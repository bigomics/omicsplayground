##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##
##

if(0) {
    reduce.sd=1000;reduce.pca=50;center.rows=TRUE;scale.rows=FALSE;umap.pkg="uwot"
    methods=c("tsne");perplexity=30;dims=2;rank.tf=FALSE 
    ##methods=c("pca","tsne","umap");dims=c(2,3);umap.pkg="uwot"
    
}

pgx.clusterSamples2 <- function(pgx, methods=c("pca","tsne","umap"), dims=c(2,3),
                                reduce.sd=1000, reduce.pca=50, perplexity=30,
                                rank.tf=FALSE, center.rows=TRUE, scale.rows=FALSE,
                                X=NULL, umap.pkg="uwot", replace.orig=TRUE )
{
    if(!is.null(X)) {
        message("using provided X matrix...")
    } else if(!is.null(pgx$X)) {
        message("using pgx$X matrix...")
        X <- pgx$X
    } else {
        message("using logCPM(pgx$counts)...")        
        ## X <- log2(1 + pgx$counts)
        X <- logCPM(pgx$counts, total=NULL)
    }    
    dim(X)
    ## X <- limma::normalizeQuantiles(X)  ##??    
    sdx <- apply(X,1,sd)
    X <- X[head(order(-sdx),reduce.sd),]
    if(center.rows) X <- X  - rowMeans(X)
    if(scale.rows)  X <- X / (1e-6+apply(X,1,sd))
    if(rank.tf) X <- scale(apply(X,2,rank))  ## works nicely
    
    dim(X)
    clust.pos <- pgx.clusterBigMatrix(
        X, methods=methods, dims=dims, perplexity=perplexity,
        reduce.sd=-1, reduce.pca=reduce.pca,
        umap.pkg=umap.pkg )
    names(clust.pos)
    clust.index <- paste0("c",clust.pos$membership)
    
    if(replace.orig) {
        message("[pgx.clusterSamples2] update tsne2d/tsne3d and 'cluster' pheno...")
        pgx$samples$cluster <- clust.index
        pgx$tsne2d <- clust.pos[["tsne2d"]]
        pgx$tsne3d <- clust.pos[["tsne3d"]]
    } else {
        message("[pgx.clusterSamples2] skipping tsne2d/tsne3d update...")
    }

    pgx$cluster <- NULL
    pgx$cluster$pos  <- clust.pos
    pgx$cluster$index <- clust.index

    pgx
}



##skipifexists=0;perplexity=NULL;sv.rank=-1;prefix="C";kclust=1;ntop=1000;fromX=FALSE;mean.center=TRUE;method="tsne";determine.clusters=1;dims=c(2,3);find.clusters=TRUE;row.center=TRUE;row.scale=FALSE;clust.detect="louvain";npca=50
pgx.clusterSamples <- function(pgx, X=NULL, skipifexists=FALSE, perplexity=30,
                               ntop=1000, npca=50, prefix="C", kclust=1, 
                               dims=c(2,3), find.clusters=TRUE,
                               clust.detect = c("louvain","hclust"),
                               row.center=TRUE, row.scale=FALSE,
                               method=c("tsne","umap","pca") )
{
    clust.detect <- clust.detect[1]    
    if(!is.null(X)) {
        message("using provided X matrix...")
        ## X <- X
    } else if(is.null(X) && !is.null(pgx$X)) {
        message("using pgx$X matrix...")
        X <- pgx$X
    } else if(is.null(X) && is.null(pgx$X)) {
        message("using pgx$counts matrix...")
        X <- logCPM(pgx$counts, total=NULL, prior=1)
    } else {
        stop("[pgx.clusterSamples] FATAL ERROR")
    }

    res <- NULL
    res <- pgx.clusterMatrix(
        X, perplexity=perplexity, dims=dims, 
        ntop=ntop, npca=npca, prefix=prefix,         
        kclust=kclust, find.clusters=find.clusters,
        clust.detect = clust.detect,
        row.center=row.center, row.scale=row.scale,
        method=method)

    if(!is.null(res$pos2d)) pgx$tsne2d <- res$pos2d
    if(!is.null(res$pos3d)) pgx$tsne3d <- res$pos3d
    if(!is.null(res$idx)) {
        if(class(pgx$samples)=="data.frame") {
            pgx$samples$cluster <- as.character(res$idx)
        } else {
            ## matrix??
            if("cluster" %in% colnames(pgx$samples)) {
                pgx$samples[,"cluster"] <- as.character(res$idx)
            } else {
                pgx$samples <- cbind(pgx$samples, cluster=as.character(res$idx))
            }
        }        
    }
    
    return(pgx)
}

pgx.FindClusters <- function(X, method=c("kmeans","hclust","louvain","meta"),
                             top.sd=1000, npca=50 )
{    
    if(0) {
        top.sd=1000;npca=50
    }
    message("[FindClusters] called...")

    km.sizes <- c(2,3,4,5,7,10,15,20,25,50,100)
    km.sizes <- km.sizes[km.sizes < ncol(X)]
    km.sizes
    if(length(method)==1 && method[1]=="meta") {
        method <- c("kmeans","hclust","louvain","meta")
    }
    
    ## reduce dimensions
    X <- head(X[order(apply(X,1,sd)),],top.sd)
    X <- t(scale(t(X))) ## scale features??
    dim(X)
    if(nrow(X) > npca) {
        npca <- min(npca, dim(X)-1)
        suppressMessages( suppressWarnings( 
            out <- irlba::irlba(X, nv=npca)
        ))
        X <- t(out$v)
    }
    
    index <- list()
    
    ## perform K-means
    if("kmeans" %in% method) {
        message("perform K-means...")
        km <- lapply(km.sizes, function(k) kmeans(t(X), k, iter.max=10))
        km.idx <- do.call(cbind,lapply(km,function(r) r$cluster))
        dim(km.idx)
        colnames(km.idx) <- paste0("kmeans.",km.sizes)    
        index[["kmeans"]] <- km.idx
    }
    
    ## perform hclust (on positions)
    if("hclust" %in% method) {
        message("perform hclust...")
        hc <- fastcluster::hclust( dist(t(X)), method="ward.D")
        hc.idx <- lapply(km.sizes, function(k) cutree(hc,k))
        hc.idx <- do.call(cbind, hc.idx)
        colnames(hc.idx) <- paste0("hclust.",km.sizes)    
        index[["hclust"]] <- hc.idx
    }

    ## perform Louvain clustering
    if("louvain" %in% method) {
        message("perform Louvain clustering...")
        require(igraph)
        library(scran)
        ##dist = as.dist(dist(t(X))) ## use 3D distance??        
        ##gr = graph_from_adjacency_matrix(1.0/dist, diag=FALSE, mode="undirected")
        gr <- scran::buildSNNGraph(X)
        ##idx <- cluster_louvain(gr)$membership
        ##gr.idx <- cluster_louvain(gr)$membership
        gr.idx <- hclustGraph(gr,k=3)  ## iterative cluster until level3
        rownames(gr.idx) <- rownames(X)
        nc <- apply(gr.idx,2,function(x) length(unique(x)))
        colnames(gr.idx) <- paste0("louvain.",nc)
        dim(gr.idx)        
        index[["louvain"]] <- gr.idx
    }

    ## find meta-index
    if("meta" %in% method && length(index)>1) {
        message("perform meta clustering...")
        K <- do.call(cbind, index)
        dim(K)
        k.rows <- split(K,row(K)) 
        d1 <- outer(k.rows, k.rows, Vectorize(function(x,y) sum(x!=y)))
        rownames(d1) <- colnames(d1) <- rownames(K)
        hc <- hclust(as.dist(d1))
        meta.idx <- do.call(cbind, lapply(km.sizes,function(k) cutree(hc,k)))
        colnames(meta.idx) <- paste0("meta.",km.sizes)
        rownames(meta.idx) <- rownames(K)
        dim(meta.idx)
        index[["meta"]] <- meta.idx
    }

    ## sort cluster index from big to small clusters
    relevelBig2small <- function(idx) as.integer(factor(idx, levels=names(sort(-table(idx)))))
    for(i in 1:length(index)) {
        index[[i]] <- apply(index[[i]],2,relevelBig2small)
    }
    
    return(index)
}

##reduce.sd=1000;reduce.pca=50;methods=c("pca","tsne","umap");dims=c(2,3);umap.pkg="uwot";center.features=TRUE;scale.features=FALSE;perplexity=30
pgx.clusterBigMatrix <- function(X, methods=c("pca","tsne","umap"), dims=c(2,3),
                                 reduce.sd = 1000, reduce.kmeans = 1000, reduce.pca = 50,
                                 center.features=TRUE, scale.features=FALSE,
                                 perplexity=30, umap.pkg="uwot")
{
    require(irlba)
    require(Rtsne)
    require(umap)

    methods <- intersect(methods, c("pca","tsne","umap"))
    if(length(methods)==0) methods <- "pca"
    
    ## Reduce dimensions by SD
    rX <- X
    if(reduce.sd>0 && nrow(rX)>reduce.sd) {
        sdx <- apply(rX,1,sd,na.rm=TRUE)
        is.constant <- all( abs(sdx-mean(sdx,na.rm=TRUE)) < 1e-8 )    
        if(is.constant) {
            cat("WARNING:: SD is constant. Skipping SD reduction...\n")                
        } else {
            cat("Reducing to ",reduce.sd," max SD features...\n")                
            rX <- rX[head(order(-sdx),reduce.sd),]
        }
    }
    if(FALSE && reduce.kmeans>0 && nrow(rX)>reduce.kmeans) {
        cat("Reducing to ",reduce.kmeans," max K-means features...\n")                
        jj <- 1:ncol(rX)
        jj <- head(order(-apply(rX,2,sd,na.rm=TRUE)), reduce.kmeans)
        km <- kmeans(rX[,jj], reduce.kmeans, iter.max=50)
        table(km$cluster)
        rX <- do.call(rbind,tapply(1:nrow(rX),km$cluster,function(i) colMeans(rX[i,,drop=FALSE])))
        dim(rX)
    }
        
    ## scale and augment if few samples
    ## rX <- scale(rX) ## columm normalization??
    if(center.features) {
        rX <- rX - rowMeans(rX,na.rm=TRUE) ## do??
    }
    if(scale.features) {
        rX <- rX / apply(rX,1,sd,na.rm=TRUE)
    }

    ## impute on row median
    rX <- imputeMedian(rX)
    
    if(ncol(rX)<=6) rX <- cbind(rX,rX,rX,rX,rX,rX)
    rX <- rX + 1e-3*matrix(rnorm(length(rX)),nrow(rX),ncol(rX))
    dim(rX)
    
    ## Further pre-reduce dimensions using SVD
    res.svd <- NULL
    if(reduce.pca>0) {
        reduce.pca <- max(3,min(c(reduce.pca,dim(rX)-1)))
        reduce.pca
        cat("Reducing to ",reduce.pca," PCA dimenstions...\n")        
        suppressMessages( suppressWarnings( 
            res.svd <- irlba(rX, nv=reduce.pca)
        ))
        rX <- t(res.svd$v) * res.svd$d
    }
    dim(rX)
    cat("dim(rX)=",dim(rX),"\n")
    
    all.pos <- list()

    if("pca" %in% methods && 2 %in% dims) {
        cat("calculating PCA 2D/3D...\n")
        require(irlba)
        if(is.null(res.svd)) {
            suppressMessages( suppressWarnings( 
                res.svd <- irlba(rX, nv=3)
            ))
        }
        pos <- res.svd$v[,1:2]
        pos <- pos[1:ncol(X),]  ## if augmented
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("PC-",c("x","y"))                
        all.pos[["pca2d"]] <- pos
    }

    if("pca" %in% methods && 3 %in% dims) {
        if(is.null(res.svd)) {
            suppressMessages( suppressWarnings(             
                res.svd <- irlba(rX, nv=3)
            ))
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
        perplexity <- pmax(min(ncol(rX)/4,perplexity),2)
        perplexity        
        pos <- Rtsne( t(rX), dims = 2,
                     ## pca = TRUE, partial_pca = TRUE,
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
        perplexity <- pmax(min(ncol(X)/4,perplexity),2)
        perplexity        
        pos <- Rtsne( t(rX[,]), dims = 3,
                     ## pca = TRUE, partial_pca = TRUE,                     
                     is_distance = FALSE, check_duplicates=FALSE,
                     perplexity = perplexity, num_threads=0)$Y
        pos <- pos[1:ncol(X),]  ## if augmented        
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("tSNE-",c("x","y","z"))        
        all.pos[["tsne3d"]] <- pos
    }
    
    if("umap" %in% methods && 2 %in% dims) {
        cat("calculating UMAP 2D...\n")
        if(umap.pkg=="uwot") {
            require(uwot)
            nb = pmax(min(ncol(X)/4,perplexity),2)
            pos <- uwot::umap( t(rX[,]), n_components=2, n_neighbors=nb)
        } else {
            require(umap)
            custom.config <- umap.defaults
            custom.config$n_components = 2
            custom.config$n_neighbors = pmax(min(ncol(X)/4,perplexity),2)
            pos <- umap::umap( t(rX[,]), custom.config )$layout
        } 
            
        ## rownames(pos) <- colnames(X)[1:nrow(pos)]
        pos <- pos[1:ncol(X),]  ## if augmented
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("UMAP-",c("x","y"))        
        all.pos[["umap2d"]] <- pos
    }
    
    if("umap" %in% methods && 3 %in% dims) {
        cat("calculating UMAP 3D...\n")
        if(umap.pkg=="uwot") {
            require(uwot)
            nb = pmax(min(ncol(X)/4,perplexity),2)
            pos <- uwot::umap( t(rX[,]), n_components=3, n_neighbors=nb )
        } else {
            require(umap)
            custom.config <- umap.defaults
            custom.config$n_components = 3
            custom.config$n_neighbors = pmax(min(ncol(X)/4,perplexity),2)
            pos <- umap::umap( t(rX[,]), custom.config )$layout
        }
        dim(pos)
        pos <- pos[1:ncol(X),]  ## if augmented
        rownames(pos) <- colnames(X)
        colnames(pos) <- paste0("UMAP-",c("x","y","z"))
        all.pos[["umap3d"]] <- pos
    }

    if(1) {
        cat("calculating Louvain memberships...\n")
        library(igraph)
        d <- min(50,dim(rX)-1)
        d
        gr <- scran::buildSNNGraph(rX, d=d)
        idx <- cluster_louvain(gr)$membership
        table(idx)
        idx <- idx[1:ncol(X)]  ## not augmented!!!!
        all.pos$membership <- idx
    }
    
    return(all.pos)
}

pgx.clusterMatrix <- function(X, perplexity=30, dims=c(2,3),
                              ntop=1000, npca=50, prefix="c",
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
    ## X <- limma::normalizeQuantiles(X)  ## in log space
    X = head( X[order(-apply(X,1,sd)),], ntop)
    ## X = t(scale(t(X),scale=TRUE))  ## really? or just centering?
    ## X = t(scale(t(X),scale=FALSE))  ## really? or just centering?
    if(row.center) X <- X - rowMeans(X,na.rm=TRUE)
    if(row.scale)  X <- (X / apply(X,1,sd,na.rm=TRUE))
    
    dim(X)
    ## some randomization is sometimes necessary if the data is 'too
    ## clean' and clusters become lines..
    ##X = X + 0.001*matrix(rnorm(length(X)),nrow(X),ncol(X))

    
    ## ------------ find t-SNE clusters
    ##max.perplexity <-  max(1,min(30, round((ncol(X)-1)/4)))
    max.perplexity <-  max(1,round((ncol(X)-1)/4))
    if(is.null(perplexity)) {
        perplexity  <- max.perplexity
    }
    if(perplexity > max.perplexity) {
        message("[pgx.clusterMatrix] perplexity too large, decreasing perplexity to ",max.perplexity)
        perplexity  <- max.perplexity
    }
    perplexity    

    ##sv.rank=20
    if(npca > 0) {
        npca <- min(npca, dim(X)-1)
        message("performing tSNE on reduced PCA k=",npca)
        suppressMessages( suppressWarnings(
            svd <- irlba(X, nv=npca)
        ))
        sv <- svd$v %*% diag(svd$d[1:ncol(svd$v)])
        rownames(sv) <- colnames(X)
        X <- t(sv)
    }

    pos2=pos3=NULL
    if(method=="umap") {
        require(uwot)
        if(2 %in% dims) {
            pos2 = uwot::umap(
                t(X),
                n_neighbors = perplexity,
                n_components = 2,
                metric = "euclidean"
            )
            colnames(pos2) <- c("umap_1","umap_2")
        }
        if(3 %in% dims) {
            pos3 = uwot::umap(
                t(X),
                n_neighbors = perplexity,
                n_components = 3,
                metric = "euclidean"
            )
            colnames(pos3) <- c("umap_1","umap_2","umap_3")
        }
    } else if(method=="tsne") {
        require(Rtsne)
        if(2 %in% dims) {
            pos2 = Rtsne( t(X), dim=2, perplexity=perplexity,
                         ##pca = TRUE, partial_pca = TRUE,                         
                         check_duplicates=FALSE, num_threads=0)$Y
            colnames(pos2) <- c("tsne_1","tsne_2")
        }
        if(3 %in% dims) {
            pos3 = Rtsne( t(X), dim=3, perplexity=perplexity,
                         ## pca = TRUE, partial_pca = TRUE,
                         check_duplicates=FALSE, num_threads=0)$Y
            colnames(pos3) <- c("tsne_1","tsne_2","tsne_3")
        }
    } else if(method=="pca") {
        require(irlba)
        suppressMessages( suppressWarnings( 
            svd <- irlba(X, nv=3)
        ))
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
        rownames(pos2) = colnames(X)
    }
    if(!is.null(pos3)) {
        rownames(pos3) = colnames(X)
    }
    
    ## ------------ find t-SNE clusters from graph
    idx = NULL
    if(find.clusters && clust.detect=="hclust") {
        cat("Finding clusters using hclust...\n")
        hc <- hclust(dist(t(X)))
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
        hc <- hclustGraph(gr,k=kclust)  ##
        dim(hc)
        idx <- hc[,min(kclust,ncol(hc))]
        ##idx <- cluster_louvain(gr)$membership
        table(idx)        
    }

    if(!is.null(idx)) {
        ## ------------ zap small clusters to "0"
        sort(table(idx))
        min.size <- pmax(3, 0.01*length(idx))
        min.size
        small.clusters <- names(which(table(idx) < min.size))
        idx[which(idx %in% small.clusters)] <- "0"
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


##is.logx=FALSE;ntop=1000
pgx.clusterMatrix.SAVE <- function(X, perplexity=30,
                                   ntop=1000, npca=50, prefix="c", 
                                   dims=c(2,3),
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

    ## X <- limma::normalizeQuantiles(X)  ## in log space
    X = head( X[order(-apply(X,1,sd)),], ntop)
    if(row.center) X <- X - rowMeans(X,na.rm=TRUE)
    if(row.scale)  X <- (X / (1e-4+apply(X,1,sd,na.rm=TRUE)))
    
    dim(X)
    ## some randomization is sometimes necessary if the data is 'too
    ## clean' and clusters become lines..
    ##X = X + 0.001*matrix(rnorm(length(X)),nrow(X),ncol(X))
    
    ## ------------ find t-SNE clusters
    ##max.perplexity <-  max(1,min(30, round((ncol(X)-1)/4)))
    max.perplexity <-  max(1,round((ncol(X)-1)/4))
    if(is.null(perplexity)) {
        perplexity  <- max.perplexity
    }
    if(perplexity > max.perplexity) {
        message("[pgx.clusterMatrix] perplexity too large, decreasing perplexity to ",max.perplexity)
        perplexity  <- max.perplexity
    }

    if(npca > 0) {
        npca <- min(npca, ncol(X))        
        message("performing tSNE on reduced PCA k=",npca)
        suppressMessages( suppressWarnings(
            svd <- irlba(X, nv=npca)
        ))
        sv <- svd$v %*% diag(svd$d[1:ncol(svd$v)])
        rownames(sv) <- colnames(X)
        X <- t(sv)
    }

    pos2=pos3=NULL
    if(method=="umap") {
        require(uwot)
        if(2 %in% dims) {
            pos2 = uwot::umap(
                t(X),
                n_neighbors = perplexity,
                n_components = 2,
                metric = "euclidean"
            )
            colnames(pos2) <- c("umap_1","umap_2")
        }
        if(3 %in% dims) {
            pos3 = uwot::umap(
                t(X),
                n_neighbors = perplexity,
                n_components = 3,
                metric = "euclidean"
            )
            colnames(pos3) <- c("umap_1","umap_2","umap_3")
        }
    } else if(method=="tsne") {
        require(Rtsne)
        if(2 %in% dims) {
            pos2 = Rtsne( t(X), dim=2, perplexity=perplexity,
                         ##pca = TRUE, partial_pca = TRUE,                         
                         check_duplicates=FALSE, num_threads=0)$Y
            colnames(pos2) <- c("tsne_1","tsne_2")
        }
        if(3 %in% dims) {
            pos3 = Rtsne( t(X), dim=3, perplexity=perplexity,
                         ## pca = TRUE, partial_pca = TRUE,
                         check_duplicates=FALSE, num_threads=0)$Y
            colnames(pos3) <- c("tsne_1","tsne_2","tsne_3")
        }
    } else if(method=="pca") {
        require(irlba)
        suppressMessages( suppressWarnings( 
            svd <- irlba(X, nv=3)
        ))
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
        rownames(pos2) = colnames(X)
    }
    if(!is.null(pos3)) {
        rownames(pos3) = colnames(X)
    }
    
    ## ------------ find t-SNE clusters from graph
    idx = NULL
    if(find.clusters && clust.detect=="hclust") {
        cat("Finding clusters using hclust...\n")
        hc <- hclust(dist(t(X)))
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
        hc <- hclustGraph(gr,k=kclust)  ##
        dim(hc)
        idx <- hc[,min(kclust,ncol(hc))]
        ##idx <- cluster_louvain(gr)$membership
        table(idx)        
    }

    if(!is.null(idx)) {
        ## ------------ zap small clusters to "0"
        sort(table(idx))
        min.size <- pmax(3, 0.01*length(idx))
        min.size
        small.clusters <- names(which(table(idx) < min.size))
        idx[which(idx %in% small.clusters)] <- "0"
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
