##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================================
##============================== GRAPH METHODS ======================================
##===================================================================================


graph_from_knn <- function(pos, k=10)
{
    ## return: edge weight are distances
    if(is.null(rownames(pos))) stop("pos must have rownames")
    
    
    
    if(ncol(pos)>3 || NCOL(pos)==1)  {
        stop("positions must be 2 or 3 columns\n")
    }
    ## use fast KNN package
    
    res = FNN::get.knn(pos, k=k)
    idx  = res$nn.index
    xval = res$nn.dist
    xval = as.vector(unlist(xval))
    sp.idx  = do.call(rbind, lapply(1:nrow(idx), function(i) cbind(i, idx[i,])))
    ##xval = exp(- xval / mean(xval) )
    ##xval = 1 / (1 + xval / mean(xval) )
    ##cat("note: edge weights are distances\n")
    sp = Matrix::sparseMatrix( i=sp.idx[,1], j=sp.idx[,2], x=xval, dims=c(nrow(pos),nrow(pos))  )
    sp = (sp + t(sp))/2
    rownames(sp) <- colnames(sp) <- rownames(pos)
    g = igraph::graph_from_adjacency_matrix(sp, mode="undirected", diag=FALSE, weighted=TRUE)
    g$layout = pos
    return(g)
}

graph_from_pos.DEPRECATED <- function(pos, trh) {
    
    if(ncol(pos)>3 || NCOL(pos)==1)  {
        stop("positions must be 2 or 3 columns\n")
    }
    d <- apply( pos, 1, function(x) Matrix::colSums((t(pos)-x)**2))  ## huge??
    w <- 1/(1 + d/mean(d))
    g = igraph::graph_from_adjacency_matrix(w,mode="undirected",weighted=TRUE,diag=FALSE)
    g = igraph::subgraph.edges(g, which(igraph::E(g)$weight > trh))
    g$layout = pos
    return(g)
}

## Calculate edge values from X
calc.edge.similarity <- function(ee, X, nk=4000, mc.cores=1)
{
    if(mc.cores>1) {
        return(calc.edge.similarityMC(ee, X, nk=nk, mc.cores=mc.cores))
    } else {
        return(calc.edge.similarityKFOLD(ee, X, nk=nk))
    }
}

calc.edge.similarityMC <- function(ee, X, nk=5000, mc.cores=4) {
    ##nk = 50000
    kfold = (nrow(ee) %/% nk) + 1
    kfold
    idx = list()
    i=1
    for(i in 1:kfold) {
        j0 <- (i-1)*nk + 1
        j1 <- min(nrow(ee), (i-1)*nk + nk)
        if(j0>nrow(ee)) next
        idx[[i]] <- j0:j1
    }
    cx <- rep(NA, nrow(ee))
    
    compute.cx <- function(ii) {
        ## NEED RETHINK!!! can be optimized bit more
        cx0 = rep(NA, length(ii))
        which.x = which(!is.na(X[ee[ii,1],1]) & !is.na(X[ee[ii,2],1]) )
        jj = ii[which.x]
        x1 = X[ee[jj,1],,drop=FALSE]
        x2 = X[ee[jj,2],,drop=FALSE]
        n1 = Matrix::rowSums( x1**2 * (!is.na(x2)), na.rm=TRUE )**0.5
        n2 = Matrix::rowSums( x2**2 * (!is.na(x1)), na.rm=TRUE )**0.5        
        n1 = (1e-20 + n1)
        n2 = (1e-20 + n2)
        xx = x1 * x2
        ## sparseMatrix does not propagate NA properly!!
        cx1 = Matrix::rowSums( xx, na.rm=TRUE) / n1 / n2
        ## any NA means missing...
        ##nx = which( Matrix::rowMeans(is.na(x1)) > 0 |
        ##            Matrix::rowMeans(is.na(x2)) > 0 )
        ##if(length(nx)>0) { cx0[nx] <- NA }
        cx0[which.x] = cx1
        cx0
    }
    res = parallel::mclapply( idx, compute.cx, mc.cores=mc.cores)
    cx = as.numeric(unlist(res))
    return(cx)
}
## Calculate edge values from X
calc.edge.similarityKFOLD <- function(ee, X, nk=4000)
{
    calc.similarity <- function(ee, X) {
        cx0 = rep(NA, nrow(ee))
        jj = which(!is.na(X[ee[,1],1]) & !is.na(X[ee[,2],1]))  ## check NA
        x1 = X[ee[jj,1],,drop=FALSE]
        x2 = X[ee[jj,2],,drop=FALSE]
        ##x1 = as.matrix(x1)
        ##x2 = as.matrix(x2)
        n1 = Matrix::rowSums( x1**2 * (!is.na(x2)), na.rm=TRUE )**0.5
        n2 = Matrix::rowSums( x2**2 * (!is.na(x1)), na.rm=TRUE )**0.5        
        n1 = (1e-20 + n1)
        n2 = (1e-20 + n2)
        xx = x1 * x2
        ## sparseMatrix does not propagate NA properly!!
        cx1 = Matrix::rowSums( xx, na.rm=TRUE) / n1 / n2
        ##nx = which( Matrix::rowMeans(is.na(x1))==1 |
        ##            Matrix::rowMeans(is.na(x2))==1 )
        ##if(length(nx)>0) { cx0[nx] <- NA }
        cx0[jj] = cx1
        return(cx0)
    }
    kfold = (nrow(ee) %/% nk) + 1
    kfold
    cx <- rep(NA, nrow(ee))
    i=1
    for(i in 1:kfold) {
        j0 <- (i-1)*nk + 1
        j1 <- min(nrow(ee), (i-1)*nk + nk)
        if(j0>nrow(ee)) next
        jj <- j0:j1
        ee1 = ee[jj,]
        cx0 = calc.similarity(ee1, X)
        cx[jj] = cx0
    }
    return(cx)
}

kcomp <- function(g,k=3,minsize=0) {
    ## Return  K-largest components
    cmp = igraph::components(g)
    sort(cmp$csize)
    cmp.size = cmp$csize
    names(cmp.size) = 1:length(cmp.size)
    cmp.size = cmp.size[which(cmp.size >= minsize)]
    cmp.top = names(cmp.size)
    if(!is.null(k)) cmp.top = Matrix::head(names(sort(-cmp.size)),k) 
    vtx = which( cmp$membership %in% cmp.top)
    igraph::induced_subgraph(g, igraph::V(g)[vtx] )
}

##g=G
cutGraph0 <- function(g, k=3, cut=FALSE)
{
    ## Cluster graph and cut crossing edges if requested.
    clust = igraph::cluster_louvain(g)
    if(cut) {
        g <- igraph::delete_edges(g, igraph::E(g)[crossing(clust, g)])
    }
    cmp = clust$membership
    cmp.size = table(cmp)
    cmp.top = Matrix::head(names(cmp.size)[order(-cmp.size)],k)
    igraph::induced_subgraph(g, which( cmp %in% cmp.top ) )
}

##g=G
graph.cut_crossings <- function(g, idx, max.wt=9999)
{
    ## Cut graph given indices of membership 
    ee = igraph::get.edgelist(g)
    dim(ee)
    jj = which( idx[ee[,1]] != idx[ee[,2]] & igraph::E(g)$weight < max.wt)
    length(jj)
    if(length(jj)>0) g = igraph::delete.edges(g, igraph::E(g)[jj] )
    return(g)    
}

itercluster_louvain <- function(g, n=3) {
    i=1
    idx = rep(1, length(igraph::V(g)))
    K = c()
    for(i in 1:n) {
        k = max(idx)
        newidx = idx
        for(i in 1:k) {
            ii = which(idx==i)
            g1 = igraph::induced_subgraph(g,ii)
            newidx[ii] = paste(i,":",cluster_louvain(g1)$membership)
        }
        levels = names(sort(table(newidx),decreasing=TRUE))
        idx = as.integer(factor(newidx, levels=levels))
        K = cbind(K, idx)
    }
    rownames(K) = igraph::V(g)$name
    colnames(K) = NULL
    table(idx)
    return(K)
}

##g=G;n=3;k=10
cutGraph <- function(g, n=2, k=5, max.wt=9999)
{
    ## Cluster graph and cut crossing edges if requested.
    idx = itercluster_louvain(g, n=n)
    g = graph.cut_crossings(g, idx, max.wt=9999)
    clust = igraph::cluster_louvain(g)
    cmp = clust$membership
    cmp.size = table(cmp)
    cmp.top = as.integer(Matrix::head(names(cmp.size)[order(-cmp.size)],k))
    igraph::induced_subgraph(g, which( cmp %in% cmp.top ) )
}



##k=2;g=G;mc.cores=16
##k=NULL
hclust_graph <- function(g, k=NULL, mc.cores=2)
{
    ## Hierarchical clustering of graph using iterative Louvain
    ## clustering on different levels. If k=NULL iterates until
    ## convergences.
    ##    
    
    idx = rep(1, length(igraph::V(g)))
    K = c()
    maxiter=100
    if(!is.null(k)) maxiter=k
    iter=1
    ok=1
    idx.len = -1
    while( iter <= maxiter && ok ) {
        old.len = idx.len
        newidx0 = newidx = idx
        i=idx[1]
        if(mc.cores>1 && length(unique(idx))>1) {
            idx.list = tapply(1:length(idx),idx,list)
            mc.cores            
            system.time( newidx0 <- parallel::mclapply(idx.list, function(ii) {
                subg = igraph::induced_subgraph(g, ii)
                subi = igraph::cluster_louvain(subg)$membership
                return(subi)
            }, mc.cores=mc.cores) )
            newidx0 = lapply(1:length(newidx0), function(i) paste0(i,"-",newidx0[[i]]))
            newidx0 = as.vector(unlist(newidx0))
            newidx = rep(NA,length(idx))
            newidx[as.vector(unlist(idx.list))] = newidx0
        } else {
            for(i in unique(idx)) {
                ii = which(idx==i)
                subg = igraph::induced_subgraph(g, ii)
                subi = igraph::cluster_louvain(subg)$membership
                newidx[ii] = paste(i,subi,sep="-")
            }
        }
        vv = names(sort(table(newidx),decreasing=TRUE))
        idx = as.integer(factor(newidx, levels=vv))
        K = cbind(K, idx)
        idx.len = length(table(idx))
        ok = (idx.len > old.len)
        iter = iter+1
    }
    rownames(K) = igraph::V(g)$name
    if(!ok && is.null(k)) K = K[,1:(ncol(K)-1)]
    dim(K)
    ##K = K[,1:(ncol(K)-1)]
    colnames(K) <- NULL
    return(K)
}
 

##method="hclust";g=G;k=10;grouped=NULL;nlouvain=3
decompose_graph <- function(g, grouped=NULL, method="louvain", k=NULL,
                            cex1=2, cex2=10, nlouvain=1)
{
    if(is.null(g$layout)) stop("graph must have layout")    
    get.idx <- function(g, method, k) {
        if(method=="louvain") {
            idx <- rep(1, length(igraph::V(g)))
            for(k in 1:nlouvain) {
                idx_new <- rep(0, length(idx))
                i=1
                for(i in unique(idx)) {
                    jj = which(idx == i)
                    g1 = igraph::induced_subgraph(g, igraph::V(g)[jj] )
                    idx1 = igraph::cluster_louvain(g1, weights=E(g1)$weight)$membership
                    idx_new[jj] = max(idx_new) + idx1
                }
                idx = idx_new
                cat("louvain: n=",length(unique(idx)),"\n")
            } ## end of nlouvain iterations
        } else if(method=="hclust") {
            if(is.null(k)) stop("hclust needs k")
            ##D = as.matrix( 1.0 / (1e-8 + g[,]) - 1)
            D = as.matrix( 1.0 - g[,] )
            ##hc = fastcluster::hclust(as.dist(D), method="ward.D2")
            hc = fastcluster::hclust(as.dist(D), method="average")
            k1 = max(min(k, nrow(D)/2),1)
            idx = cutree(hc, k1)
            table(idx)
        } else if(method=="nclust") {
            
            if(is.null(k)) stop("nclust needs k")
            D = as.matrix( 1.0 - g[,] )
            ##hc = nclust(as.matrix(g[,]), link="ward")
            ##hc = nclust( as.matrix(g[,]), link="average")
            hc = nclust( as.dist(D), link="average")
            k1 = max(min(k,nrow(D)/2),1)
            idx = cutree(hc, k1)
            table(idx)
        } else {
            stop("fatal. unknown method",method,"\n")
        }
        return(idx)
    }

    idx = NULL
    group.idx = NULL
    if(!is.null(grouped)) {
        grouped[is.na(grouped)] = "(missing)"
        groups = unique(grouped)
        table(grouped)
        ng = length(igraph::V(g))
        idx <- rep(0, ng)
        group.idx <- list()
        i=1
        for(i in 1:length(groups)) {
            jj <- which(grouped==groups[i])
            sub.idx = rep(1,length(jj))
            length(jj)
            if(length(jj)>3) {
                subg = igraph::induced_subgraph(g, jj )
                k0 = max(ceiling(k * (length(jj)/ng)),2)
                k0
                sub.idx = get.idx(subg, method=method, k=k0)
                table(sub.idx)
            }
            table(sub.idx)
            sub.idx = sub.idx + max(idx)
            idx[jj] = sub.idx
            group.idx[[i]] <- unique(sub.idx)
        }
    } else {
        idx = get.idx(g, method=method, k=k)
    }
    table(idx)
    
    ## resort indices
    idx = as.integer(factor(idx, levels=order(-table(idx))))    
    names(idx) = igraph::V(g)$name
    table(idx)

    ## adjacency matrix
    M = as.matrix(g[,])
    diag(M) <- 1
    M[1:4,1:4]    
    
    ## downsample adjacency matrix
    ngrp = length(unique(idx))
    ngrp
    subM = tapply(1:ncol(M), idx, function(ii) rowMeans(M[,ii,drop=FALSE],na.rm=TRUE))
    subM = lapply(subM, function(x) tapply(x, idx, mean))
    subM = matrix(unlist(subM),ncol=ngrp, nrow=ngrp)
    sum(is.na(subM))
    rownames(subM) <- colnames(subM) <- paste0("cluster.",1:nrow(subM))
    M[1:4,1:4]
    subM[1:4,1:4]
    
    ## create downsamples graph
    subG = igraph::graph_from_adjacency_matrix(
        subM, weighted=TRUE, diag=FALSE, mode="undirected")
    igraph::V(subG)$name <- rownames(subM)
    subG$members <- tapply( idx, idx, names)
    names(subG$members) <- rownames(subM)
    
    ##avg.pos = apply(g$layout,2,function(x) tapply(x,idx,median))
    avg.pos = apply(g$layout[names(idx),],2,function(x) tapply(x,idx,mean))
    rownames(avg.pos) <- rownames(subM)
    subG$layout = avg.pos
    
    subgraphs <- list()
    for(i in 1:length(igraph::V(subG))) {
        vv = subG$members[[i]]
        child <- igraph::induced_subgraph(g, vv)
        child$layout <- g$layout[vv,]
        subgraphs[[i]] <- child
    }

    group.idx2 <- rep(1, length(subgraphs))
    if(!is.null(grouped)) {
        group.idx2 <- sapply(1:length(subgraphs), function(i)
            which(sapply(group.idx,function(x) (i %in% x))))
    }
    
    subG = igraph::simplify(subG)
    subG$group.idx = group.idx2
    subG$membership = idx
    subG$subgraphs = subgraphs
    return(subG)
}


##method="hclust";g=G;k=10;grouped=NULL;nlouvain=3
downsample_graph.DEPRECATED <- function(g, idx=NULL, grouped=NULL, merge.op="max")
{
    if(is.null(g$layout)) stop("graph must have layout")    

    if(is.null(idx)) {
        hc = hclust_graph(g, k=NULL)
        idx = hc[,ncol(hc)]
    }
    
    if(!is.null(grouped)) {
        grouped[is.na(grouped)] = "_other"
        table(grouped)
        grouped.idx <- idx
        groups = unique(grouped)
        i=1
        for(i in 1:length(groups)) {
            jj <- which(grouped==groups[i])
            grouped.idx[jj] = paste0(i,":",idx[jj])
        }
        idx = grouped.idx
    } else {
        idx = as.character(idx)
    }
    table(idx)
    
    ## resort indices
    idx.levels = names(sort(-table(idx)))
    idx = factor(idx, levels=idx.levels)
    table(idx)

    ## get adjacency matrix
    adjM = g[,]
    diag(adjM) <- 0
    adjM[1:4,1:4]    
    
    ## downsample by taking mean of merge blocks in adjM
    
    ngrp = length(unique(idx))
    ngrp
    if(merge.op=="mean") {
        subM = tapply(1:ncol(adjM), idx, function(ii) rowMeans(adjM[,ii,drop=FALSE],na.rm=TRUE))
        subM = parallel::mclapply(subM, function(x) tapply(x, idx, mean, na.rm=TRUE))
    } else if(merge.op=="max") {
        subM = tapply(1:ncol(adjM), idx, function(ii) apply(adjM[,ii,drop=FALSE],1,max,na.rm=TRUE))
        subM = parallel::mclapply(subM, function(x) tapply(x, idx, max, na.rm=TRUE))
    } else {
        stop("unknown merge.op",merge.up,"\n")
    }

    subM.idxnames = names(subM)
    subM = matrix(unlist(subM),ncol=ngrp, nrow=ngrp)
    sum(is.na(subM))
    rownames(subM) <- colnames(subM) <- subM.idxnames
    adjM[1:4,1:4]
    dim(subM)
    subM[1:4,1:4]

    ## diagonal??
    if(0) {
        for(i in 1:length(subM.idxnames)) {
            ii = which(idx==subM.idxnames[i])
            subM[i,i] = min(adjM[ii,ii],na.rm=TRUE)
        }
    }

    ## create downsampled graph
    subG = igraph::graph_from_adjacency_matrix(
        subM, weighted=TRUE, diag=FALSE, mode="undirected")
    igraph::V(subG)$name <- rownames(subM)
    members <- tapply( igraph::V(g)$name, idx, list)
    names(members) <- rownames(subM)
    igraph::V(subG)$members = members
    ## subG$members = members
    ## igraph::V(subG)$size = sapply(members, length)
    
    ##avg.pos = apply(g$layout,2,function(x) tapply(x,idx,median))
    avg.pos = apply(g$layout[igraph::V(g)$name,],2,function(x) tapply(x,idx,mean))
    rownames(avg.pos) <- rownames(subM)
    subG$layout = avg.pos
            
    subG = igraph::simplify(subG)
    names(idx) = igraph::V(g)$name
    subG$membership = idx
    return(subG)
}



##===================================================================================
##============================== END OF FILE ========================================
##===================================================================================
