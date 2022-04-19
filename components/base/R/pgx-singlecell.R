##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


##install.packages("base2grob")
##install.packages("plotrix")
##install.packages("multipanelfigure")
##outs="~/IRB/Gonzalez_10x_combined/outs/"

seurat2pgx <- function(obj, do.cluster=FALSE) {
    ## Convert a Seurat object to a minimal PGX object.
    ##
    ##
    message("[createPGX.10X] creating PGX object...")        
    ## pgx <- pgx.createPGX(
    ##     counts = counts,  samples = pheno, contrasts = ct$contr.matrix,
    ##     do.cluster = FALSE, batch.correct = TRUE, 
    ##     is.logx = FALSE)
    pgx <- list()
    pgx$name   <- "SeuratProject"
    pgx$description  <- "Seurat object converted using seurat2pgx"
    pgx$date   <- Sys.Date()
    pgx$datatype  <- "scRNA-seq"
    pgx$counts <- obj[["RNA"]]@counts
    pgx$X <- obj[["RNA"]]@data
    pgx$samples <- obj@meta.data
    
    pgx$genes <- ngs.getGeneAnnotation(genes=rownames(pgx$counts))  
    rownames(pgx$genes) <- rownames(pgx$counts)
    ##pgx$genes <- cbind(pgx$genes, obj[["RNA"]]@meta.features)
    
    Matrix::head(pgx$samples)
    Matrix::head(pgx$genes)
    if(do.cluster) {
        message("[seurat2pgx] clustering samples")
        pgx <- pgx.clusterSamples2(
            pgx, dims=c(2,3), methods=c("pca","tsne","umap")
        )
        names(pgx$cluster$pos)
    }   

    ## copy clustering from Seurat
    if("tsne" %in% names(obj@reductions)) {
        pos <- obj@reductions[["tsne"]]@cell.embeddings[,1:2]
        pgx$cluster$pos[["tsne2d"]] <- pos
        pgx$tsne2d <- pos
        pgx$tsne3d <- cbind(pos,0)
    }
    if("umap" %in% names(obj@reductions)) {
        pgx$cluster$pos[["umap2d"]] <- obj@reductions[["umap"]]@cell.embeddings[,1:2]
    }    
    pgx$samples$cluster <- obj@meta.data[,"seurat_clusters"]
    names(pgx)
    pgx
}


pgx.createPGX.10X.DEPRECATED <- function(outs, ncells=2000, aggr.file="aggregation.csv")
{    
    ## Read 10X count matrix and aggregation file from outs folder and
    ## create PGX file.
    ##
    message("[createPGX.10X] reading 10X data file...")
    counts.10x <- Seurat::Read10X(file.path(outs,"filtered_feature_bc_matrix"))
    dim(counts.10x)

    message("[createPGX.10X] reading aggregation meta file...")
    aggr <- read.csv(file.path(outs,aggr.file), comment.char="#")
    aggr$molecule_h5 <- NULL
    Matrix::head(aggr)
    idx <- as.integer(sub(".*-","",colnames(counts.10x)))
    pheno <- aggr[idx,]
    rownames(pheno) <- colnames(counts.10x)

    ## QC filter cells
    message("[createPGX.10X] filtering outliers...")    
    counts <- pgx.scFilterOutliers(counts.10x, a=2)
    pheno <- pheno[colnames(counts),]    
    dim(counts)
    dim(pheno)

    ## Reduce number of cells (approx. equalizing libraries)
    message("[createPGX.10X] reducing cells...")    
    out <- pgx.reduceCells(counts, method="subsample",
                           ncells=ncells, group.id=pheno$library_id)
    counts <- out$counts
    pheno <- pheno[colnames(counts),]
    table(pheno$library_id)
    
    ## pre-filter expressed genes
    message("[createPGX.10X] filtering features...")    
    sel <- (rowMeans(counts>=3) > 0.01)
    counts <- counts[sel,]
    dim(counts)

    ## Infer celltype
    message("[createPGX.10X] inferring cell types...")        
    pheno$cell.type <- pgx.inferCellType(counts, low.th=0.03, add.unknown=TRUE)

    ## for the moment just celltype comparison
    ct <- makeDirectContrasts(pheno[,"cell.type",drop=FALSE], ref="others")
    names(ct)
    pheno$group <- ct$group
    #pheno$group <- NULL

    message("[createPGX.10X] creating PGX object...")        
    pgx <- pgx.createPGX(
        counts = counts,  samples = pheno, contrasts = ct$contr.matrix,
        ## batch.correct = FALSE,
        is.logx = FALSE)
    Matrix::head(pgx$genes)
    names(pgx)
    
    ## We perform single-cell integration/batch correction only to
    ## improve clustering!! Note that count/X matrix in PGX is not
    ## replaced with the integrated matrix (not reliable...).
    colnames(pheno)
    if("batch" %in% colnames(pheno)) {

        message("[createPGX.10X] performing batch integration using MNN...")                
        bX <- logCPM(counts, total=1e4)
        bb <- pgx$samples[,"batch"]
        bX <- pgx.scBatchIntegrate(bX, batch=bb, method="MNN")        
        ## pgx <- pgx.clusterSamples(pgx, X=bX, method="tsne", dims=2)
        pgx <- pgx.clusterSamples2(pgx, X=bX, methods=c("pca","tsne","umap"), dims=c(2,3))
        pgx$tsne2d <- pgx$cluster$pos[["tsne2d"]]
        pgx$tsne3d <- pgx$cluster$pos[["tsne3d"]]

        ## update default clusters from tSNE/UMAP
        xpos <- pgx$cluster$pos[["umap3d"]]
        ## xpos <- pgx$cluster$pos[["tsne2d"]]
        ## xpos <- cbind(pgx$cluster$pos[["tsne3d"]],pgx$cluster$pos[["umap3d"]])
        ## xpos <- do.call(cbind, pgx$cluster$pos)
        idx <- pgx.FindClusters(xpos, method="louvain")[[1]][,1]
        table(idx)
        ##pgx$samples$cluster <- pgx$cluster$index[["kmeans"]][,"kmeans.10"]
        pgx$samples$cluster <- paste0("C",idx)        
    }
    message("[createPGX.10X] done!")                
    pgx
}

pgx.scQualityControlPlots.NOTFINISHED <- function(counts) {
    
    mt.genes <- grep("^MT-",rownames(counts),ignore.case=TRUE,value=TRUE)
    rb.genes <- grep("^RP[SL]",rownames(counts),ignore.case=TRUE,value=TRUE)
    percent.mito <- Matrix::colSums(counts[mt.genes,]) / Matrix::colSums(counts)*100
    percent.ribo <- Matrix::colSums(counts[rb.genes,]) / Matrix::colSums(counts)*100
    nfeature <- Matrix::colSums(counts>0)
    ncounts  <- Matrix::colSums(counts)
        
    T1 <- table(sample.id, qc.ok)
    barplot( t(T1), ylab="number of cells")
    legend("topright", legend=c("qc=TRUE","qc=FALSE"),
           fill=c("grey80","grey30"), cex=0.9)

    ##P <- data.frame(nfeature, ncounts, percent.mito, percent.ribo)
    P <- data.frame(percent.mito, percent.ribo) 
    P1 <- apply(P, 2, function(x) tapply(x, sample.id, mean))
    barplot(t(P1), ylab="percentage (%)")
    legend("topright", legend=c("mito","ribo"),fill=c("grey80","grey30"),cex=0.9)
    
    P2 <- apply(P, 2, function(x) tapply(x, cell.type, mean))
    barplot(t(P2), ylab="percentage (%)")
    legend("topright", legend=c("mito","ribo"),fill=c("grey80","grey30"),cex=0.9)
}



##zth=0.25;ref=NULL;samples=NULL;is.count=TRUE
pgx.scTestDifferentialExpression <- function(counts, y, is.count=TRUE, samples=NULL,
                                             ref=NULL, zth=0.25, cpm.total=1e4,
                                             do.nzero=TRUE)
{
    X1 <- NULL
    if(FALSE && is.null(cpm.total)) {
        cpm.total <- mean(Matrix::colSums(counts > 0)) ## mean number of expressing genes
        message("setting column sums to total = ",round(cpm.total,2))
    }

    if(is.count && class(counts)=="dgCMatrix") {
        message("input matrix is counts (sparseMatrix)")
        ## normalize to CPM
        X1 <- counts
        X1@x <- cpm.total * X1@x / rep.int(Matrix::colSums(X1), diff(X1@p))  ## fast divide by columns sum
        ## compute log-expression
        X1@x <- log2(1 + X1@x)
        ##X1 <- limma::normalizeQuantiles(X1)
        zth <- log2(1 + zth)
    } else if(is.count) {
        message("input matrix is counts (Matrix)")
        ## normalize to CPM
        X1 <- cpm.total * t(t(counts) / (1e-8 + Matrix::colSums(counts)) )
        X1 <- log2(1 + X1)
        ##X1 <- limma::normalizeQuantiles(X1)
        zth <- log2(1 + zth)
    } else {
        message("input matrix is normalized log-expression")
        X1 <- counts
        counts <- pmax(2**X1-1,0) ## estimate counts
    }
    
    if(length(setdiff(unique(y),NA))!=2) {
        stop("phenotype must have exactly two groups")
    }
    y <- factor(y)
    if(!is.null(ref)) y <- relevel(y, ref=ref)
    y0 <- levels(y)[1]
    y1 <- levels(y)[2]
    y0
    y1
    
    ## non-zero matrix (all zeros treated as NA)
    Z1 <- X1
    Z1[Z1 <= zth] <- NA  ## notice threshold!
    
    ## t-test
    
    m1 <- matrixTests::row_t_welch(as.matrix(X1)[,y==y1], as.matrix(X1)[,y==y0]) ## all
    m2 <- matrixTests::row_t_welch(as.matrix(Z1)[,y==y1], as.matrix(Z1)[,y==y0]) ## non-zero
    
    
    pct.x   <- m2$obs.x / m1$obs.x
    pct.y   <- m2$obs.y / m1$obs.y
    pct.tot <- m2$obs.tot / m1$obs.tot
    pct.diff <- (pct.x - pct.y)

    ## both k1 and k2 cannot be zero...
    hack.err <- which(m2$obs.x==0 & m2$obs.y==0)
    if(length(hack.err)>0) m2$obs.x[hack.err] <- 1
    ##pv <- corpora::fisher.pval( m2$obs.x, m1$obs.x, m2$obs.y, m1$obs.y, alternative="less")
    ##pv <- corpora::fisher.pval( m2$obs.x, m1$obs.x, m2$obs.y, m1$obs.y, alternative="greater")        
    pv <- corpora::fisher.pval( m2$obs.x, m1$obs.x, m2$obs.y, m1$obs.y)
    ##pv <- corpora::chisq.pval( m2$obs.x, m1$obs.x, m2$obs.y, m1$obs.y)
    Matrix::head(pv)
    if(length(hack.err)>0) m2$obs.x[hack.err] <- 0    
    pct <- cbind(pct.x, pct.y, pct.diff, pct.tot, pct.pvalue=pv)
    
    kk <- c("obs.x","obs.y","obs.tot","mean.x","mean.y","mean.diff","pvalue")
    m1 <- m1[,kk]
    m2 <- m2[,kk]
    colnames(m2) <- paste0("nzero.",sub("^nzero[.]","",colnames(m2)))

    b1 <- NULL
    if(!is.null(samples)) {
        ## geometric average (in log-space) by sample
        S1 <- tapply(1:ncol(X1),samples,function(i) rowMeans(X1[,i,drop=FALSE]))
        S1 <- do.call(cbind,S1)
        yb <- tapply(y, samples, function(y1) names(which.max(table(y1))))
        b1 <- matrixTests::row_t_welch(S1[,yb==y1], S1[,yb==y0]) ## all
        b1 <- b1[,kk]
        colnames(b1) <- paste0("sample.",colnames(b1))
    }
   
    b2 <- NULL
    if(!is.null(samples)) {
        ## sum counts then to average (in log-space) by sample
        C1 <- pmax(2**X1-1,0)
        S2 <- tapply(1:ncol(C1),samples,function(i) rowMeans(C1[,i,drop=FALSE]))
        S2 <- do.call(cbind,S2)
        ##S2 <- t(t(S2) / Matrix::colSums(S2)) * cpm.total
        S2 <- log2(1 + S2) 
        yb <- tapply(y, samples, function(y1) names(which.max(table(y1))))
        b2 <- matrixTests::row_t_welch(S2[,yb==y1], S2[,yb==y0]) ## all
        b2 <- b2[,kk]
        colnames(b2) <- paste0("sample2.",colnames(b2))
        remove(C1)
    }

    ## Group fold-change as used by Seurat: difference of
    ## log-average. Averages are done in linear count space. Note: no
    ## p-values are computed (not possible).
    grpx <- log2(1+rowMeans(X1[,y==y1]))
    grpy <- log2(1+rowMeans(X1[,y==y0]))
    fc0   <- grpx - grpy
    nz.grpx <- log2(1+rowMeans(Z1[,y==y1],na.rm=TRUE))
    nz.grpy <- log2(1+rowMeans(Z1[,y==y0],na.rm=TRUE))
    nz.fc0 <- nz.grpx - nz.grpy
    grp <- cbind( x = grpx, y = grpy, diff = fc0, 
                 nzero.x = nz.grpx, nzero.y = nz.grpy,
                 nzero.diff = nz.fc0)
    colnames(grp) <- paste0("group.",colnames(grp))
    rownames(grp) <- rownames(X1)

    ## Bulk-like fold-change: sum up all counts in one group,
    ## normalize with CPM, then do difference in the log. Note: no
    ## p-values are computed (not possible).
    bulk.x <- rowSums(counts[,y==y1])
    bulk.y <- rowSums(counts[,y==y0])
    bulk.x <- bulk.x / sum(bulk.x) * cpm.total
    bulk.y <- bulk.y / sum(bulk.y) * cpm.total
    bulk.fc   <- bulk.x - bulk.y
    grp2 <- cbind( x = bulk.x, y = bulk.y, diff = bulk.fc)
    colnames(grp2) <- paste0("bulk.",colnames(grp2))
    rownames(grp2) <- rownames(counts)
    
    ##df <- cbind(m1, m2, pct, b1, b2)
    df <- cbind(m1, pct, grp, grp2)
    if(do.nzero) df <- cbind(df, m2)
    if(!is.null(samples)) df <- cbind(df, b1, b2)

    
    P <- df[,grep("pvalue",colnames(df))]
    P[is.na(P)] <- 1
    meta.p <- apply(P,1,function(p) metap::sumlog(p)$p)
    ##meta.p <- apply(P,1,function(p) metap::sumz(p)$p)
    meta.q <- p.adjust(meta.p)
    df$meta.pvalue <- meta.p
    df$meta.qvalue <- meta.q

    ## simplify column names
    colnames(df) <- sub("mean.diff","diff",colnames(df))
    
    df
}

pgx.reduceCells <- function(counts, method, ncells, pheno=NULL, group.id=NULL) {
    dim(counts)

    if(ncol(counts) > ncells) {    
        if(method=="pool") {
            message(">> Pooling cells...")
            ## sample.id <- sub(".*-","",colnames(counts))
            ## group.id <- paste0(sample.id,"-",cell.type)
            pc <- pgx.poolCells(counts, ncells, groups=group.id)
            counts <- pc$counts
            if(!is.null(pheno)) {
                pheno1 <- apply(pheno, 2, function(x) 
                    tapply(x, pc$cluster.id, function(aa) names(which.max(table(aa))))
                    )
                pheno <- pheno1[colnames(counts),]
                rownames(pheno) <- colnames(counts)
            }
        } else if(method=="subsample" && is.null(group.id)) {
            message(">> Subsampling cells...")
            sel <- sample(ncol(counts), ncells)
            counts <- counts[,sel]
            if(!is.null(pheno)) {
                pheno <- pheno[sel,]
                rownames(pheno) <- colnames(counts)
            }
        } else if(method=="subsample" && !is.null(group.id)) {
            message(">> Subsampling cells in groups...")
            table(group.id)        
            n1 <- round(ncells / length(unique(group.id)))
            sel <- tapply(1:ncol(counts), group.id, function(i) Matrix::head(sample(i),n1))
            sel <- unlist(sel)
            counts <- counts[,sel]
            if(!is.null(pheno)) {
                pheno <- pheno[sel,]
                rownames(pheno) <- colnames(counts)
            }
        } else {
            stop("FATAL ERROR:: unknown method",method)
        }
    }    
    list(counts=data.matrix(counts), pheno=data.frame(pheno))
}


##ncells=1000;groups=NULL;verbose=TRUE;clust.method="umap";prior=1
##counts=counts1;ncells=N;groups=ctype1
pgx.poolCells <- function(counts, ncells, groups=NULL, stats="sum",
                          clust.method="umap", prior=1, X=NULL,
                          meta=NULL, verbose=TRUE) {
    
    if(is.null(groups)) {
        groups <- rep("grp1",ncol(counts))
    } 
    groups <- as.integer(factor(as.character(groups)))
    table(groups)
    
    pool.counts <- c()
    cluster.id <- c()

    if(is.null(X)) {
        ## compute log-expression matrix from counts
        X <- counts
        X <- 1e6 * t(t(X) / (1e-8 + Matrix::colSums(X)))  ## CPM
        nz <- Matrix::which(X!=0, arr.ind=TRUE)
        X[nz] <- log2(1 + X[nz])        
    }
    ##X <- limma::normalizeQuantiles(X)
    ##X <- X - Matrix::rowMeans(X,na.rm=TRUE)
    
    ##k=1000;topsd=1000;nv=20;method="tsne";prior=1
    clusterX <- function(X1, k, method, topsd=1000, nv=50) {
        if(ncol(X1) <= k) {
            cluster <- paste0("c",1:ncol(X1))
            names(cluster) <- colnames(X1)
            return(cluster)
        }
        avgx <- Matrix::rowMeans(X1)
        ##sdx <- apply(X1,1,sd)
        sdx <- (Matrix::rowMeans(X1**2) - avgx**2)**0.5
        wt <- sdx * avgx
        X1 <- X1[head(order(-wt),topsd),,drop=FALSE]
        ##X1 <- limma::normalizeQuantiles(X1)        
        X1 <- X1 - Matrix::rowMeans(X1)  ## center features
        nv <- min(nv,ncol(X1)-1)
        sv <- irlba::irlba(X1,nv=nv)
        ##V <- sv$v %*% diag(sv$d) ## ???
        V <- sv$v %*% diag(sv$d**0.5) ## ???
        V <- t(scale(t(V)))  ## rowscale
        cluster <- NULL
        if(method=="tsne") {
            px <- min(nrow(V)/4,30)
            pos <- Rtsne::Rtsne(V, perplexity=px)$Y
            cluster <- kmeans(pos,k,iter.max=100)$cluster
        } else if(method=="umap") {
            pos <- uwot::umap(V)  ## umap is tighter and faster than t-SNE
            cluster <- kmeans(pos,k,iter.max=100)$cluster
        } else if(method=="kmeans") {
            cluster <- kmeans(V,k,iter.max=100)$cluster
        } else if(method=="hclust") {
            cluster <- cutree(fastcluster::hclust(dist(V)),k)
            ## cluster <- cutree(fastcluster::hclust(dist(t(X1))),k)
            ##cluster  <- cutree(fastcluster::hclust(as.dist(1-cor(as.matrix(X1)))),k)
        } else {
            stop("ERROR:: unknown clustering method")
        }
        ##cluster <- kmeans(po s,k)$cluster
        cluster <- paste0("c",cluster)
        names(cluster) <- colnames(X1)        
        cluster
    }

    ngroup <- length(unique(groups))
    if(verbose) {
        message("[pgx.poolCells] ",ngroup," groups defined")
        message("[pgx.poolCells] clustering method = ",clust.method)        
    }

    cluster.id <- rep(NA,ncol(counts))
    names(cluster.id) <- colnames(counts)
    all.groups <- sort(unique(groups))
    g <- groups[1]
    for(g in all.groups) {
        if(verbose) message("[pgx.poolCells] clustering group",g," ...")
        ## do quick clustering 
        sel <- which(groups == g)
        X1 <- X[,sel,drop=FALSE]
        dim(X1)
        k <- ceiling(ncells/ngroup) ## equalize between groups
        k
        ngroup
        cluster <- NULL
        if(ncol(X1) <= k) {
            cluster <- paste0("g",g,"-c",1:ncol(X1))
        } else {
            cluster <- clusterX(X1, k=k, method=clust.method)
            table(cluster)
            if(!is.null(groups)) cluster <- paste0("g",g,"-",cluster)
        }        
        names(cluster) <- colnames(X1)        
        cluster.id[colnames(X1)] <- cluster
    }

    if(verbose) message("[pgx.poolCells] pooling cells...")
    cluster.id <- cluster.id[colnames(counts)]
    table(cluster.id)            
    ##idx <- tapply(1:ncol(counts), cluster.id, list)
    ##pool.counts <- parallel::mclapply(idx, function(ii) Matrix::rowSums(counts[,ii,drop=FALSE]))   
    ##pool.counts <- tapply(1:ncol(counts), cluster.id, function(ii)
    ## Matrix::rowSums(counts[,ii,drop=FALSE]))
    if(stats=="mean") {
        pool.counts <- tapply(1:ncol(counts), cluster.id, function(ii)
            Matrix::rowMeans(counts[,ii,drop=FALSE]))
    }
    if(stats=="sum") {
        pool.counts <- tapply(1:ncol(counts), cluster.id, function(ii)
            Matrix::rowSums(counts[,ii,drop=FALSE]))
    }
    pool.counts <- do.call(cbind, pool.counts)

    new.meta <- NULL
    if(!is.null(meta)) {
        max.class <- function(x) names(which.max(table(x)))        
        new.meta <- apply( meta, 2, function(x)
            tapply(x, cluster.id, max.class))
        new.meta <- data.frame(new.meta, check.names=FALSE)
    }
    
    res <- list()
    res$cluster.id <- cluster.id
    res$counts <- pool.counts
    res$meta <- new.meta
    return(res)
}

pgx.scBatchIntegrate <- function(X, batch,
                               method=c("ComBat","limma","CCA","MNN","Harmony","liger"))
{
    ## method <- method[1]
    res <- list()
    ##pos <- list()
    if("ComBat" %in% method) {
        message("[pgx.scBatchIntegrate] single-cell batch correction using ComBat...")
        
        ## ComBat correction
        res[["ComBat"]] <- sva::ComBat(X, batch=batch, par.prior=TRUE)
    }
    if("limma" %in% method) {
        message("[pgx.scBatchIntegrate] single-cell batch correction using LIMMA...")
        
        ## LIMMA correction
        res[["limma"]] <- limma::removeBatchEffect(X, batch=batch)
    }
    if("CCA" %in% method) {
        message("[pgx.scBatchIntegrate] single-cell batch correction using CCA (Seurat)...")
        ## Seurat CCA correction
        counts <- pmax(2**X - 1,0)
        try( res[["CCA"]] <- pgx.SeuratBatchIntegrate(counts, batch=batch) )
    }
    if("MNN" %in% method) {
        message("[pgx.scBatchIntegrate] single-cell batch correction using MNN...")        
        ## MNN correction        
        try( mnn <- mnnCorrect(X, batch=batch, cos.norm.in=TRUE, cos.norm.out=FALSE))
        res[["MNN"]] <- MultiAssayExperiment::assays(mnn)[["corrected"]]
    }
    if("Harmony" %in% method) {
        message("[pgx.scBatchIntegrate] single-cell batch correction using Harmony...")                
        
        ## Harmony corrections
        nv <- min(floor(dim(X)*0.8),30)
        out <- irlba::irlba(X, nu=nv, nv=nv)
        V <- t(out$v)
        dim(V)
        meta_data <- data.frame(batch=batch)
        try( hm <- harmony::HarmonyMatrix(
                                V, meta_data, 'batch',
                                do_pca = FALSE, npcs=nv,
                                return_object=TRUE)
            )
        dim(hm$Z_corr) ## corrected PCA embeddings
        hX <- (out$u %*% diag(out$d) %*% hm$Z_corr)
        dimnames(hX) <- dimnames(X)
        res[["Harmony"]] <- hX
    }
    if("liger" %in% method) {
        message("[pgx.scBatchIntegrate] single-cell batch correction using LIGER...")
        
        xlist <- tapply(1:ncol(X), batch, function(i) pmax(2**X[,i]-1,0))
        liger <- rliger::createLiger(xlist, take.gene.union=TRUE)
        liger <- rliger::normalize(liger)
        ##liger <- liger::selectGenes(liger, var.thresh=1e-8, alpha.thresh=0.999, do.plot=TRUE)
        ##liger <- liger::selectGenes(liger, num.genes=100, tol=1e-2, do.plot=TRUE)        
        liger@var.genes <- Matrix::head(rownames(X)[order(-apply(X,1,sd))],100)
        liger <- rliger::scaleNotCenter(liger)        
        vg <- liger@var.genes
        vg
        xdim <- sapply(xlist,ncol)
        xdim
        k=15
        k <- round(min(30,length(vg)/3,median(xdim/2)))
        k
        ## OFTEN GIVES ERROR!!!!!
        liger <- try( rliger::optimizeALS(liger, k=k) )
        ##liger <- fliger::optimizeALS(liger, k=k, lambda=5, nrep=3)
        ##liger  <- fliger::online_iNMF(liger, k=k, miniBatch_size=1000, max.epochs = 5)
        if(class(liger)=="try-error") {
            ## res[["LIGER"]] <- NULL
        } else {
            liger <- rliger::quantile_norm(liger)
            ##liger <- rliger::quantile_norm(liger, knn_k=round(k/2), min_cells=k, dims.use=1:k)
            ##liger <- louvainCluster(liger, resolution = 0.25)
            ##liger <- runUMAP(liger, distance='cosine', n_neighbors=30, min_dist=0.3)
            cX <- t(liger@H.norm %*% liger@W)
            dim(cX)
            cat("[pgx.scBatchIntegrate] WARNING:: LIGER returns smaller matrix")
            res[["liger"]] <- cX
        }
    }
    if(length(res)==1) {
        res <- res[[1]]
    }
    res
}

## qc.filter=FALSE;nanchors=-1;sct=FALSE
pgx.SeuratBatchIntegrate <- function(counts, batch, qc.filter=FALSE,
                                     nanchors=-1, sct=FALSE)
{    
    ##
    ## From Seurat vignette: Integration/batch correction using
    ## CCA. Note there is no QC filtering for samples on ribo/mito
    ## content. You need to do that before.
    ##
    
    dim(counts)
    nbatch <- length(unique(batch))
    message("[pgx.SeuratBatchIntegrate] Processing ",nbatch," batches...")        
    obj.list <- list()
    i=1
    b=batch[1]
    batches <- unique(batch)
    batches
    for(i in 1:length(batches)) {
        sel <- which(batch==batches[i])
        ##sel <- Matrix::head(sel,200)
        counts1 <- counts[,sel]
        obj <- Seurat::CreateSeuratObject(counts1)
        if(sct) {
            obj <- Seurat::SCTransform(obj, vars.to.regress = NULL, verbose = FALSE)
            ##obj <- Seurat::SCTransform(obj)
        } else {
            obj <- Seurat::NormalizeData(obj, normalization.method="LogNormalize",
                                 scale.factor=10000, verbose=FALSE)
            obj <- Seurat::FindVariableFeatures(obj, selection.method="vst", verbose=FALSE)
        }
        obj$batch <- b
        obj.list[[i]] <- obj
    }
    
    anchor.features=NULL
    sct
    if(sct) {
        ## See: https://satijalab.org/seurat/v3.0/integration.html
        nfeatures = nrows(counts)
        if(nanchors>0) nfeatures = nanchors 
        anchor.features <- Seurat::SelectIntegrationFeatures(
            object.list = obj.list, nfeatures=nanchors)
        obj.list <- Seurat::PrepSCTIntegration(
            object.list = obj.list,
            anchor.features = anchor.features, 
            verbose = FALSE)        
    } else {        
        ##anchor.features = rownames(counts) ## really?
        anchor.features = nrow(counts) ## really?        
        if(nanchors>0) anchor.features=nanchors ## really?
    }
    
    message("[pgx.SeuratBatchIntegrate] Finding anchors...")
    options(future.globals.maxSize = 8*1024^3) ## set to 8GB
    ##NUM.CC=10;k.filter=10
    NUM.CC = max(min(20, min(table(batch))-1),1)
    NUM.CC
    bdims <- sapply(obj.list,ncol)
    kmax <- max(min(bdims)-1,1)
    kmax
    normalization.method = ifelse(sct,"SCT","LogNormalize")
    message("[pgx.SeuratBatchIntegrate] NUM.CC = ",NUM.CC)
    message("[pgx.SeuratBatchIntegrate] normalization.method = ",normalization.method)
    ##anchors <- Seurat::FindIntegrationAnchors(obj.list, dims = 1:NUM.CC)
    
    anchors <- Seurat::FindIntegrationAnchors(
        obj.list, dims = 1:NUM.CC,
        k.filter = min(200,kmax),
        k.anchor = min(5,kmax),
        k.score = min(30,kmax),
        anchor.features = anchor.features,
        normalization.method = normalization.method,
        verbose = FALSE)
    
    message("[pgx.SeuratBatchIntegrate] Integrating data...")
    len.anchors <- length(anchors@anchor.features)
    message("[pgx.SeuratBatchIntegrate] number of anchors = ",len.anchors)    

    integrated <- Seurat::IntegrateData(anchorset = anchors)
    
    integrated <- Seurat::IntegrateData(
        anchorset = anchors,
        k.weight = min(100,kmax),  ##  troublesome...
        dims = 1:NUM.CC,
        normalization.method = normalization.method,        
        verbose = FALSE)
    dim(integrated)
    dim(counts)
    
    key <- ifelse(sct, "SCT", "integrated")
    mat.integrated <- as.matrix(integrated[[key]]@data)
    mat.integrated <- exp(mat.integrated) ## natural log!!
    dim(mat.integrated)
    mat.integrated <- mat.integrated[,colnames(counts)]
    
    ## set previously zero counts to zero again
    zc <- Matrix::which(counts[rownames(mat.integrated),]==0,arr.ind=TRUE)
    mat.integrated[zc] <- 0

    ## set missing to zero...
    ##mat.integrated[is.nan(mat.integrated)] <- 0
    mat.integrated[is.na(mat.integrated)] <- 0
    
    if(!nrow(mat.integrated)==nrow(counts)) {
        cat("WARNING:: number of rows of integrated matrix has changed!")
    }
    
    return(mat.integrated)    
}


pgx.scFilterOutliers <- function(counts, a=2.5, plot=FALSE)
{
    ## QC filter of (single) cells like Seurat
    ## See https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

    ##--------------------------------------------------------------
    ## calculate percentages
    ##--------------------------------------------------------------        
    mt.genes <- grep("^MT-",rownames(counts),ignore.case=TRUE,value=TRUE)
    rb.genes <- grep("^RP[SL]",rownames(counts),ignore.case=TRUE,value=TRUE)
    percent.mito <- Matrix::colSums(counts[mt.genes,])/Matrix::colSums(counts)*100
    percent.ribo <- Matrix::colSums(counts[rb.genes,])/Matrix::colSums(counts)*100
    nfeature <- Matrix::colSums(counts>0)
    ncounts  <- Matrix::colSums(counts)
    
    if(plot) {
        ##a=2.5
        log.nfeature <- log10(1 + nfeature)
        log.ncounts <- log10(1 + ncounts)
        log.mito <- log10(1 + percent.mito)
        log.ribo <- log10(1 + percent.ribo)
        nfeature.th <- mean(log.nfeature) + c(-a,a) * sd(log.nfeature)
        ncounts.th <- mean(log.ncounts) + c(-a,a) * sd(log.ncounts)    
        mito.th <- mean(log.mito) + c(-a,a) * sd(log.mito)    
        ribo.th <- mean(log.ribo) + c(-a,a) * sd(log.ribo)
        mito.th
        ribo.th
        
        par(mfrow=c(2,2))
        hist(log.nfeature, breaks=100)
        abline(v=nfeature.th, col="red")
        hist(log.ncounts, breaks=100)
        abline(v=ncounts.th, col="red")
        
        hist(log.mito, breaks=100)
        abline(v=mito.th, col="red")    
        hist(log.ribo, breaks=100)
        abline(v=ribo.th, col="red")

    }

    selectInlier <- function(x, a=2.5) {
        xmin <- mean(x) - a * sd(x)
        xmin <- max(xmin, 0.01*mean(x))
        xmax <- mean(x) + a * sd(x)
        x > xmin & x < xmax
    }

    selectInlier <- function(x, a=2.5) {
        x <- log10(1 + x)
        xmin <- mean(x) - a * sd(x)
        xmin <- max(xmin, 0.01*mean(x))
        xmax <- mean(x) + a * sd(x)
        (x > xmin & x < xmax)
    }

    ## sel <- nfeature > 200 & nfeature < 7500 & percent.mito < 10 & percent.ribo < 50
    ## sel <- nfeature < nfeature.th & ncounts < ncounts.th &
    ## percent.mito < mito.th & percent.ribo < ribo.th
    sel <- selectInlier(nfeature,a) &
        selectInlier(ncounts,a) &
        selectInlier(percent.mito,a) &
        selectInlier(percent.ribo,a)
    table(sel)

    counts <- counts[,sel]
    counts
}


pgx.createSeuratObject <- function(counts, aggr.csv=NULL,
                                   project="SeuratProject", max.cells=2000 )
{
    
    obj <- Seurat::CreateSeuratObject(counts, min.cells=5, project=project)    

    if(!is.null(aggr.csv)) {        
        aggr <- read.csv(aggr.csv)
        sample.idx <- as.integer(sub(".*-","",colnames(counts)))
        pheno <- aggr[sample.idx,c("library_id","phenotype")]
        ##colnames(pheno) <- c("sample","phenotype")
        rownames(pheno) <- colnames(counts)

        
        obj@meta.data <- cbind(obj@meta.data, meta.data)
    }

    Seurat::DefaultAssay(obj) <- "RNA"
    obj@project.name <- project
    
    ## QC filtering of cells
    message("performing QC filtering")        
    obj$percent.mito <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-|^Mt-")
    obj$percent.ribo <- Seurat::PercentageFeatureSet(obj, pattern = "^RP[LS]|^Rp[ls]")
    summary(obj$nFeature_RNA)
    obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mito < 20)

    ## Equalizing libraries
    if("library_id" %in% colnames(obj@meta.data)) {
        table(obj$library_id)
        ncells <- median(table(obj$library_id))
        message("library_id parameter found. Equalizing cells to ",ncells)        
        sel <- unlist(tapply(1:ncol(obj), obj$library_id, head, ncells))
        table(obj$library_id[sel])
    }
    
    ##obj <- subset(obj, downsample=4000)
    colnames(obj@meta.data)    
    if("batch" %in% colnames(obj@meta.data)) {
        message("batch parameter found. integrating batches using MNN.")        
        split <- Seurat::SplitObject(obj, split.by = "batch")
        for(i in 1:length(split)) {
            split[[i]] <- Seurat::NormalizeData(split[[i]])
            split[[i]] <- Seurat::FindVariableFeatures(split[[i]])
        }        
        anchors <- Seurat::FindIntegrationAnchors(
            split, dims = 1:30, verbose = FALSE)        
        genes <- rownames(counts)
        integrated <- Seurat::IntegrateData(
            anchorset = anchors, features.to.integrate = genes,
            dims = 1:30, verbose=FALSE)
        dim(integrated)
        obj <- integrated
        Seurat::DefaultAssay(obj) <- "integrated"
    } else {
        obj <- Seurat::NormalizeData(obj)
        obj <- Seurat::FindVariableFeatures(obj)
        Seurat::DefaultAssay(obj) <- "RNA"
    } 

    if(ncol(obj)>max.cells) {
        message("Subsampling cells to ",max.cells)        
        obj <- subset(obj, cells = sample(Seurat::Cells(obj),max.cells))
    }
    
    ## Dimensionality reductions
    message("Calculating dimensionality reductions...")        
    obj <- Seurat::FindVariableFeatures(obj)   
    ##obj <- Seurat::ScaleData(obj, vars.to.regress = c("percent.mito","percent.ribo"))
    obj <- Seurat::ScaleData(obj, vars.to.regress = c("percent.mito"))
    obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:30, dim.embed=2)
    obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:30, dim.embed=2)    
    obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = 1:30)
    obj <- Seurat::FindClusters(obj)

    ## Pre-calculate markers
    if(0) {
        ## cannot store the results somewhere????
        Seurat::Idents(obj) <- "seurat_clusters"
        markers <- Seurat::FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, thresh.use=0.25)
        attributes(obj)$markers <- list()
        attributes(obj)$markers[["seurat_clusters"]] <- markers
    }
    message("Finished!")        
    obj
}


pgx.createSeurateFigures <- function(obj)
{    
    if(0) {
        obj <- pgx.createSeuratObject(pgx)        
        source(file.path(RDIR,"pgx-include.R"))
        pgx.SampleClusterPanel(pgx, pheno="cluster", pheno2=c("cell.type","batch","phenotype"))
        pgx.SampleClusterPanel(pgx, pheno="cell.type", pheno2=c("cluster","batch","phenotype"))
    }
    
    caption1 <- paste("Project:",obj@project.name,"   Date:",Sys.Date())
    caption1
    fig <- list()
    
    ##----------------------------------------------------------------------
    ## QC 
    ##----------------------------------------------------------------------
    Seurat::Idents(obj) <- "library_id"  ## important!!
    Seurat::Idents(obj) <- "phenotype"  ## important!!
    Seurat::Idents(obj) <- "batch"  ## important!!
    Seurat::Idents(obj) <- "orig.indent"  ## important!!
    vplot <- Seurat::VlnPlot(
        obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"),
        pt.size=0.15, ncol = 4)
    vplot <- vplot & ggplot2::xlab(NULL)
    if(0) {
        table(obj$library_id)
        vplot
    }
    
    q1 <- Seurat::FeaturePlot(obj, features = "nCount_RNA")
    q2 <- Seurat::FeaturePlot(obj, features = "nFeature_RNA")
    q3 <- Seurat::FeaturePlot(obj, features = "percent.mito")
    q4 <- Seurat::FeaturePlot(obj, features = "percent.ribo")
    
    plot1 <- Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                            pt.size=0.3)
    plot2 <- Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mito",
                            pt.size=0.3)
    plot3 <- Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.ribo",
                            pt.size=0.3)
        
    theme0 <- ggplot2::theme(
        plot.title = ggplot2::element_text( size=12, face='bold'),
        ##legend.key.size = grid::unit(0.45,"lines")            
        legend.title = ggplot2::element_text(size=11),
        legend.text = ggplot2::element_text(size=9),
        legend.key.size = grid::unit(0.55,"lines"),            
        axis.text.x = ggplot2::element_text(size=10),
        axis.text.y = ggplot2::element_text(size=10),
        axis.title.x = ggplot2::element_text(size=11),
        axis.title.y = ggplot2::element_text(size=11)
    )
    ##(vplot / qq / pp) & theme0
        
    qq <- (q1 | q2 | q3 | q4)
    qq <- qq & ggplot2::guides(color = ggplot2::guide_colourbar(barwidth= 0.5, barheight=3))
    
    pp <- (plot1 | plot2 | plot3) & ggplot2::ggtitle(NULL)
    pp <- pp & ggplot2::theme(legend.key.size = grid::unit(0.45,"lines"))
    ##pp
    
    v1 <- cowplot::plot_grid(vplot & theme0, ncol=1)
    v2 <- cowplot::plot_grid(qq & theme0, ncol=1)
    v3 <- cowplot::plot_grid(pp & theme0, ncol=1) 

    fig1 <- (v1 / v2 / v3)
    fig1 <- fig1 +
        patchwork::plot_annotation(
            title = 'Seurat QC plots',
            subtitle = 'These plots show the QC of your experiment.',
            ##caption = 'made with Omics Playground',
            caption = caption1,
            theme = ggplot2::theme(plot.title = ggplot2::element_text(size=16, face='bold')),
            tag_levels = 'a')
    ##fig1
    fig[["QC"]] <- fig1

    ##----------------------------------------------------------------------
    ## Variable features & PCA
    ##----------------------------------------------------------------------

    vg  <- obj[["RNA"]]@var.features
    hvg <- Matrix::head(Seurat::VariableFeatures(obj), 20)    
    v1 <- Seurat::VariableFeaturePlot(obj) %>% Seurat::LabelPoints(points=hvg, repel=TRUE)
    v1 <- v1 & ggplot2::theme(legend.position = c(0.02, 0.95) )
    
    ##JackStrawPlot(obj, dims=1:15)
    r1 <- Seurat::ElbowPlot(obj) + ggplot2::ggtitle("PCA elbow plot")
    r2 <- base2grob::base2grob(~{
        par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,1,0))
        Seurat::DimHeatmap(obj, dims = 1:9, ncol=3, cells = 500,
                   nfeatures=16, balanced = TRUE)
    })
    
    v1 <- v1 + ggplot2::ggtitle("Variable features") + theme0 +
        ggplot2::theme(plot.margin = ggplot2::margin(0,1,0.5,0, "cm"))
    
    r1 <- r1 + ggplot2::ggtitle("PCA elbow plot") + theme0 +
        ggplot2::theme(plot.margin = ggplot2::margin(0,0.5,0.5,1, "cm"))
        
    fig3 <- (v1 | r1) / r2 + patchwork::plot_layout(design="A\nB\nB") &
        ggplot2::theme( plot.title = ggplot2::element_text(hjust=0)) 

    fig3 <- fig3 +
        patchwork::plot_annotation(
            title = 'Seurat PCA plots',
            subtitle = 'These plots show the PCA of your experiment.',
            ## caption = 'made with Omics Playground',
            caption = caption1,
            theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face='bold')),
            tag_levels = 'a')

    ##fig3    
    fig[["PCA"]] <- fig3
        
    ##----------------------------------------------------------------------
    ## Cluster markers
    ##----------------------------------------------------------------------
    Seurat::Idents(obj) <- "seurat_clusters"
    
    markers <- Seurat::FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, thresh.use=0.25)
    ##markers <- markers[order(markers$avg_logFC*log(markers$p_val)),]
    markers$score <- with(markers, -log(p_val) * avg_logFC * pct.1 * (1 - pct.2))
    
    ##markers <-  attributes(obj)$markers[["seurat_clusters"]]
    top1 <- markers %>% plotly::group_by(cluster) %>% dplyr::top_n(1, -p_val)
    top1 <- markers %>% plotly::group_by(cluster) %>% dplyr::top_n(1, score)
    top1$gene
    vplot2 <- Seurat::VlnPlot( obj, features = top1$gene, pt.size=0.15, ncol = 5)
    vplot2 <- vplot2 & ggplot2::xlab(NULL)
    vplot2 <- vplot2 & ggplot2::ylab("expression")
    vplot2 <- vplot2 & ggplot2::theme(plot.margin = ggplot2::margin(1,3,1,3, "mm"))
    vplot2 <- vplot2 & ggplot2::theme(axis.text.x=element_text(angle=0))

    ntop <- floor(65 / length(unique(obj$seurat_clusters)))
    ##top <- markers %>% plotly::group_by(cluster) %>% dplyr::top_n(ntop, avg_logFC)
    top <- markers %>% plotly::group_by(cluster) %>% dplyr::top_n(ntop, score)
    top
    h1 <- Seurat::DoHeatmap(obj, features=top$gene, angle=0, size=4) + Seurat::NoLegend()
    
    vplot2  <- cowplot::plot_grid(vplot2 & theme0, ncol=1)
    h1x  <- cowplot::plot_grid(h1, ncol=1)
    fig4 <- (vplot2) / h1x + patchwork::plot_layout(design="A\nB\nB")

    fig4 <- fig4 +  patchwork::plot_annotation(
                        title = 'Seurat cluster markers',
                        subtitle = 'These plots show the clusters of your experiment.',
                        ## caption = 'made with Omics Playground',
                        caption = caption1,
                        theme = ggplot2::theme(plot.title = ggplot2::element_text(size=16, face='bold')),
                        tag_levels = 'a')
    ##fig4

    fig[["cluster.markers"]] <-fig4
    
    ##----------------------------------------------------------------------
    ## Celltype assignment
    ##----------------------------------------------------------------------

    CANONICAL.MARKERS = list(
        "B_cells" = c("Ms4a1","Cd79a","Cd79b", "Fcmr","Ebr1"),
        "CD4 T cells" = c("Il17r","Ccr7","Cd3e","Cd3d","Cd3g","Cd4"),
        "CD8 T cells" = c("Cd8a","Cd8b"),
        "NK_cells" = c("Gnly","Nkg7","Gzma", "Klrb1c", "Klrk1", "Klra4"),
        "Dendritic_cells" = c("Fcer1a","Cst3","Siglech","Fscn1","Ccl22"),
        "Macrophages" = c("C1qa","C1qb","C1qc","Lyz2"),
        "Monocytes" = c("Cd14","Lyz","Ly6c2","Fcgr3a","Ms4a7"), 
        ##"B16.melanoma" = c("Mlana", "Dct", "Tyrp1", "Mt1", "Mt2", "Pmel", "Pgk1"),
        "Platelet" = "Ppbp"
    )

    ct1 <- pgx.inferCellType(obj[["RNA"]]@counts, add.unknown=FALSE, low.th=0.01)
    table(ct1)
    ##ct1 <- pgx.inferCellType(obj[["RNA"]]$counts, add.unknown=FALSE)
    tapply(ct1, obj$seurat_clusters, function(x) table(x))
    ct1x <- tapply(ct1, obj$seurat_clusters, function(x) names(which.max(table(x))))    
    ct1x
    
    obj$cell.type <- ct1x[as.character(obj$seurat_clusters)]
    table(obj$cell.type)
    
    d1 <- Seurat::DimPlot(obj, group.by = "seurat_clusters", label=TRUE) +
        ggplot2::ggtitle("Seurat clusters") + Seurat::NoLegend()
    d4 <- Seurat::DimPlot(obj, group.by = "cell.type", label=TRUE) +
        ggplot2::ggtitle("Cell type")  + Seurat::NoLegend()

    sel <- unlist(tapply(1:ncol(obj), obj$cell.type, head, 300))
    markers$cell.type <- ct1x[as.character(markers$cluster)]
    ##markers <- markers[order(markers$cell.type,-markers$avg_logFC),]
    ntop <- floor(45 / length(unique(ct1x)))
    top <- markers %>% plotly::group_by(cell.type) %>% dplyr::top_n(ntop, score)
    Matrix::head(top)
    h2 <- Seurat::DoHeatmap(obj[,sel], features=top$gene, group.by="cell.type",
                    hjust=0.5, angle=0, size=4) + Seurat::NoLegend()
    h2 <- h2 + ggplot2::theme(plot.margin = ggplot2::margin(0,0,0,0,"mm"))
    ##h2x  <- cowplot::plot_grid(h2, ncol=1)
    
    fig5 <- h2 / (d1 + theme0 | d4 + theme0) + patchwork::plot_layout(design="A\nA\nB")
    ## cowplot::plot_grid(d1, d2, d3, d4, ncol=2)
    
    fig5 <- fig5 +  patchwork::plot_annotation(
                        title = 'Seurat cell type identification',
                        subtitle = 'Assignment of cell type identity to clusters.',
                        ## caption = 'made with Omics Playground',
                        caption = caption1,
                        theme = ggplot2::theme(
                            plot.title = ggplot2::element_text(size = 16, face='bold')),
                        tag_levels = 'a')
    ##fig5
    fig[["cell.type"]] <- fig5
    
    ##----------------------------------------------------------------------
    ## Phenotype dimension plots
    ##----------------------------------------------------------------------

    Matrix::head(obj@meta.data)
    not.ph <- grep("ident|nCount|nFeat|mito|[.]mt|ribo|_snn|^percent",colnames(obj@meta.data),value=TRUE)
    ph <- setdiff(colnames(obj@meta.data),not.ph)
    ph
    dd <- list()
    for(p in ph) {
        dd[[p]] <- Seurat::DimPlot(obj, group.by=p, label=TRUE) +
            ggplot2::ggtitle(p) + Seurat::NoLegend() + theme0
    }
    ##d4 <- Seurat::DimPlot(obj, group.by = "cell.type") + ggplot2::ggtitle("cell.type")
    ##fig2 <- cowplot::plot_grid(plotlist=dd, nrow=2, ncol=3)
    fig2 <- patchwork::wrap_plots(dd)

    fig2 <- fig2 +  patchwork::plot_annotation(
                        title = 'Seurat Phenotype plots',
                        subtitle = 'Distribution of phenotypes on the t-SNE/UMAP.',
                        ## caption = 'made with Omics Playground',
                        caption = caption1,
                        theme = ggplot2::theme(plot.title = ggplot2::element_text(size=16, face='bold')),
                        tag_levels = 'a')
    ##fig2

    fig[["phenotypes"]] <- fig2

    
    ##----------------------------------------------------------------------
    ## Gene families
    ##----------------------------------------------------------------------
    if(0) {
        sort(table(sub(":.*","",colnames(pgx$GMT))))
        
        G1 <- pgx$GMT[,grep("HGNC",colnames(pgx$GMT),value=TRUE)]
        hgnc.list <- lapply(apply(G1!=0,2,which),names)
        percent <- lapply(lapply(apply(G1!=0,2,which),names),
                          function(gg) Matrix::colSums(obj[["RNA"]]@counts[gg,]/Matrix::colSums(obj[["RNA"]]@counts)))
        percent <- round(100*do.call(cbind, percent),digits=2)
        colnames(percent) <- gsub("[ ]",".",sub("HGNC:","percent.",colnames(percent)))
        obj@meta.data <- obj@meta.data[,!colnames(obj@meta.data) %in% colnames(percent)]
        obj@meta.data <- cbind(obj@meta.data, percent)
        dd <- list()
        for(p in colnames(percent)) {
            dd[[p]] <- Seurat::FeaturePlot(obj, features=p)
        }
        patchwork::wrap_plots(dd, ncol=3) & theme0

        ## markers with LIMMA (multifactorial)
        y <- obj$cell.type
        design <- model.matrix(~ 1 + y)
        res <- limma::eBayes(limma::lmFit(obj[["RNA"]]@data, design))
        top <- limma::topTable(res, number=1000)
        top <- top[order(top$P.Value),]
        Matrix::head(top,20)


    }
    
    ##----------------------------------------------------------------------
    ## Specific DE markers
    ##----------------------------------------------------------------------

    if(0) {
        Seurat::Idents(obj) <- "phenotype"
        Matrix::head(obj@meta.data)
        stats <- Seurat::FindMarkers(
            obj, ident.1 = "treated", ident.2 = "untreated",
            logfc.threshold = 0.1, ##min.pct = 0.1,
            verbose = FALSE)
        Matrix::head(stats)
        top.deg <- Matrix::head(rownames(stats),9)
        
        vp1 <- Seurat::VlnPlot(obj, features = top.deg, split.plot=FALSE,
                       split.by = "phenotype", group.by = "cell.type",
                       pt.size=0.15, ncol = 3, log=TRUE)
        vp1  <- vp1 & ggplot2::xlab("")
        vp1  <- vp1 & ggplot2::ylab("expression")
        vp1  <- vp1 & theme0
        vp1
        
        ggplot2::theme(plot.margin = ggplot2::margin(0,0,0,0, "cm"))
        vplot2  <- vplot2 & theme0
        ##vplot2
        
        h1 <- Seurat::DoHeatmap(obj, features=top$gene, angle=0, size=4) + Seurat::NoLegend()
        fig4 <- (vplot2) / h1
    }
   
    ##----------------------------------------------------------------------
    ## Return all figures
    ##----------------------------------------------------------------------
    return(fig)
}


if(0) {

    ##
    ## From Seurat vignette: standard QC filtering
    ##
    ##
    
    obj <- Seurat::CreateSeuratObject(counts)    
    obj <- Seurat::AddMetaData(obj, sample.id, col.name = "Sample.id")    
    obj
    slotNames(obj[["RNA"]])
    table(Seurat::Idents(obj))    
    
    Matrix::head(obj@meta.data)
    mt.genes <- rownames(obj)[grep("^MT-",rownames(obj),ignore.case=TRUE)]
    C <- Seurat::GetAssayData(object = obj, slot = "counts")    
    percent.mito <- Matrix::colSums(C[mt.genes,])/Matrix::colSums(C)*100
    hist(percent.mito, breaks=100)
    obj <- Seurat::AddMetaData(obj, percent.mito, col.name = "percent.mito")
    
    rb.genes <- rownames(obj)[grep("^RP[SL]",rownames(obj),ignore.case=TRUE)]
    percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
    obj <- Seurat::AddMetaData(obj, percent.ribo, col.name = "percent.ribo")
    
    Seurat::VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mito","percent.ribo"),
            group.by = "Sample.id",pt.size = 0.1, ncol=4) + Seurat::NoLegend()
    Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    Seurat::FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.mito")
    Seurat::FeatureScatter(obj, feature1 = "percent.ribo", feature2 = "nFeature_RNA")
    
    ## QC select cells
    obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 &
                           percent.mito < 10 & percent.ribo < 50)
    dim(obj)
    
    ## total count normalization
    hist(Matrix::colSums(exp(obj[["RNA"]]@data[,1:1000])),breaks=100)
    obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)    
    hist(Matrix::colSums(exp(obj[["RNA"]]@data[,1:1000])),breaks=100)
    
    ## Find highly variable top 2000 genes
    obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    Seurat::VariableFeaturePlot(obj)
    
    ## scale (=standardize????) really??? you will loose fold-change!
    all.genes <- rownames(obj)
    obj <- Seurat::ScaleData(obj, features=all.genes)
    
    obj <- Seurat::RunPCA(obj, features=VariableFeatures(object=obj))
    Seurat::ElbowPlot(obj)    
}


if(0) {

    ##
    ## From Seurat vignette: Batch correction
    ##
    

    batch = c("BioReplicate1", "BioReplicate2", "BioReplicate3","BioReplicate4",
              "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6",
              "Tattoo1","Tattoo2","Tattoo3","Tattoo4")
    treatment = c("control","control","control","control",
                  "treated","treated","treated","treated","control","control",
                  "neg_control","neg_control","neg_control","neg_control")
    
    obj.list <- list()
    i=1
    for(i in 1:length(outputs)) {
        data.10x <- Seurat::Read10X(data.dir=file.path(outputs[i],"/outs/filtered_feature_bc_matrix"))
        dim(data.10x)
        ## data.10x <- data.10x[,1:1000]  ## just subsample FTM...
        celseq <- Seurat::CreateSeuratObject(data.10x, min.cells=5)
        ##celseq <- FilterCells(celseq, subset.names = "nGene", low.thresholds = 800)
        celseq <- subset(celseq, subset = nFeature_RNA > 500)
        ##-----
        celseq[["percent.mt"]] <- Seurat::PercentageFeatureSet(celseq, pattern = "^mt-")
        celseq <- subset(celseq, subset = percent.mt < 5)
        ##head(celseq@meta.data, 5)
        ##VlnPlot(celseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.3)
        ##-----
        ##celseq <- Seurat::NormalizeData(celseq)
        celseq <- Seurat::NormalizeData(object = celseq, normalization.method = "LogNormalize",
                                scale.factor = 10000)
        celseq <- Seurat::FindVariableFeatures(celseq, selection.method="vst",
                                       do.plot = F, display.progress = F)
        ##celseq$batch <- gsub(".*/|_scRNA","",outputs[i])
        celseq$batch <- batch[i]
        celseq$treatment <- treatment[i]
        obj.list[[i]] <- celseq
    }
    
    NUM.CC = 20
    anchors    <- Seurat::FindIntegrationAnchors(obj.list, dims = 1:NUM.CC,
                                         ##anchor.features=9999,
                                         verbose=FALSE)
    integrated <- Seurat::IntegrateData(anchorset = anchors,
                                dims = 1:NUM.CC,
                                verbose=FALSE)
    dim(integrated)
    Seurat::DefaultAssay(integrated) <- "integrated"
    dim(integrated)
    
    ## obj <- Seurat::CreateSeuratObject(counts)    
    ## obj <- Seurat::AddMetaData(obj, sample.id, col.name = "Sample.id")    
    ## obj
    ## slotNames(obj[["RNA"]])
    ## table(Seurat::Idents(obj))    

}

