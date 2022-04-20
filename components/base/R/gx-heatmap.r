##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##



##bluered <- function(n=64) gplots::colorpanel(n,"blue","grey90","red")
##bluered <- function(n=64) gplots::colorpanel(n,"dodgerblue4","grey90","indianred3")
##bluered <- function(n=64) gplots::colorpanel(n,"royalblue3","grey90","indianred3")
BLUERED64 <- colorRampPalette(c("royalblue3","grey90","indianred3"))(64) ## as in pgx-init.R

gx.markermap <- function(X, splitx, n=5, ...)
{
    F <- tapply(1:ncol(X), splitx, function(i) rowMeans(X[,i,drop=FALSE]))
    F <- do.call(cbind,F)
    F <- F - rowMeans(F)
    topg <- apply(F,2,function(x) Matrix::head(order(-x),n))
    topg <- apply(topg, 2, function(i) rownames(F)[i])
    gg <- as.vector(topg)
    gg.idx <- as.vector(sapply(colnames(topg),rep,n))
    ##gx.splitmap(X[gg,], splitx=splitx, split=gg.idx, nmax=999, show_colnames=FALSE, softmax=TRUE, scale='none')
    gx.splitmap(X[gg,], splitx=splitx, split=gg.idx, ...)
}

gx.PCAheatmap <- function(X, nv=5, ngenes=10, ...) {
    if(class(X)=="list" && "X" %in% names(X)) {
        X <- X$X - rowMeans(X$X)
    }
    X <- X - rowMeans(X,na.rm=TRUE)
    res <- irlba::irlba(X, nv=nv)
    gg <- apply(res$u,2,function(x) Matrix::head(rownames(X)[order(-abs(x))],ngenes))
    sx <- paste0("PC.",as.vector(sapply(1:nv,rep,10)))
    gx.splitmap( X[gg,], split=sx, ... )
}

gx.PCAcomponents <- function(X, nv=20, ngenes) {
    if(class(X)=="list" && "X" %in% names(X))  X <- X$X
    res <- irlba::irlba(X - rowMeans(X), nv=nv)
    ##par(mfrow=c(5,4), mar=c(1,1,2,8))
    for(i in 1:nv) {
        gg <- Matrix::head(rownames(X)[order(-abs(res$u[,i]))],ngenes)
        X1 <- X[gg,]
        X1 <- (X1 - rowMeans(X1)) / (1e-4 + apply(X1,1,sd)) ## scale??
        colnames(X1) <- NULL
        gx.imagemap(X1, main=paste0("PC",i), cex=0.8)
    }
}

gx.imagemap <- function(X, main="", cex=1, clust=TRUE)
{
    if(clust) {
        ii <- fastcluster::hclust(dist(X))$order
        jj <- fastcluster::hclust(dist(t(X)))$order
        X <- X[ii,jj]
    }
    Matrix::image( 1:ncol(X), 1:nrow(X), t(X), xaxt="n", yaxt="n",
          xlab="", ylab="")
    mtext(colnames(X), side=1, at=1:ncol(X), las=3,
          cex=0.75*cex, line=0.5 )
    mtext(rownames(X), side=4, at=1:nrow(X), las=1,
          cex=0.85*cex, line=0.5 )
    title(main=main, cex=cex, line=0.5)
}

gx.splitmap <- function(gx, split=5, splitx=NULL,
                        clust.method = "ward.D2",
                        ##dist.method = "pearson",
                        dist.method = "euclidean",                       
                        col.dist.method = "euclidean",
                        plot.method = "heatmap.2",
                        ## col=bluered(64),
                        scale="row", softmax=0, order.groups="clust", symm.scale=FALSE,
                        ## Rowv = NA, Colv = NA,
                        cluster_rows=TRUE, cluster_columns=TRUE, sort_columns=NULL,
                        col.annot=NULL, row.annot=NULL, annot.ht=3,
                        nmax=1000, cmax=NULL, main=" ", verbose=1, denoise = 0,
                        cexRow=1, cexCol=1, mar=c(5,5,5,5), rownames_width = 25,
                        rowlab.maxlen = 20,
                        title_cex=1.2, column_title_rot=0, column_names_rot=90,
                        show_legend=TRUE, show_key=TRUE, zlim=NULL,
                        show_rownames=nmax, lab.len=80, key.offset=c(0.05,1.01),
                        show_colnames=NULL, use.nclust=FALSE)
{
    
    
    
    ComplexHeatmap::ht_global_opt(fast_hclust = TRUE)
    ##
    par(xpd=FALSE)
    
    if(0) {
        split=5; splitx=NULL;
        clust.method="ward.D2";
        dist.method="euclidean";
        col.dist.method="euclidean";
        plot.method="heatmap.2";
        ## col=bluered(64);
        scale="row"; softmax=0; order.groups="clust"; symm.scale=FALSE;
        ## Rowv = NA; Colv = NA;
        cluster_rows=TRUE; cluster_columns=TRUE; sort_columns=NULL;
        col.annot=NULL; row.annot=NULL; annot.ht=3;
        nmax=1000; cmax=NULL; main=" "; verbose=1; denoise = 0;
        cexRow=1; cexCol=1; mar=c(5,5,5,5); rownames_width = 25;
        rowlab.maxlen = 20;
        title_cex=1.2; column_title_rot=0; column_names_rot=90;
        show_legend=TRUE; show_key=TRUE; zlim=NULL;
        show_rownames=nmax; lab.len=80; key.offset=c(0.05,1.01);
        show_colnames=NULL; use.nclust=FALSE
        
    }
    
    if(verbose>1) cat("input.dim.gx=",dim(gx),"\n")
    if(is.null(show_colnames)) {
        show_colnames <- ifelse(ncol(gx)>100, FALSE, TRUE)
    }
    ## give unique name if duplicated
    if(sum(duplicated(rownames(gx)))>0) {
        rownames(gx) <- tagDuplicates(rownames(gx))
        if(!is.null(row.annot)) rownames(row.annot) <- rownames(gx)
    }
    if(!is.null(split) && length(split)==1 && split==1) split <- NULL
    if(!is.null(splitx) && length(splitx)==1 && splitx==1) splitx <- NULL
    if(is.null(main)) main = '  '
    if(length(mar)==1) mar <- rep(mar[1],4) ## old style
    if(length(mar)==2) mar <- c(mar[1],5,5,mar[2]) ## old style


    ## Any scaling before selection??? This should be done before
    ## filtering on SD.
    
    ##-------------------------------------------------------------
    ## scaling options
    ##-------------------------------------------------------------
    ## if(any(c("col","cols","both") %in% scale)) gx = scale(gx)
    ## if(any(c("row","rows","both") %in% scale)) gx = t(scale(t(gx)))
    ## if("col.center" %in% scale) gx = scale(gx,center=TRUE,scale=FALSE)
    ## if("row.center" %in% scale) gx = gx - rowMeans(gx,na.rm=TRUE)
    if("col" %in% scale || "both" %in% scale) {
        ##gx = scale(gx)
        tgx <- t(gx) - colMeans(gx,na.rm=TRUE)
        gx <- t(tgx / (1e-4 + apply(gx, 2, sd,na.rm=TRUE))) ## small EPS maintains SD order!
        remove(tgx)
    }
    if("row" %in% scale || "both" %in% scale) {
        ##gx = t(scale(t(gx)))
        gx <- gx - rowMeans(gx,na.rm=TRUE)
        gx <- gx / (1e-4 + apply(gx, 1, sd,na.rm=TRUE))  ## small EPS maintains SD order!
    }
    if("col.center" %in% scale) {
        gx  <- t( t(gx) - colMeans(gx, na.rm=TRUE))
    }
    if("row.center" %in% scale) {
        gx  <- gx - rowMeans(gx, na.rm=TRUE)
    }    
    if("row.bmc" %in% scale && is.null(splitx)) {
        gx  <- gx - rowMeans(gx, na.rm=TRUE)
    }
    if("row.bmc" %in% scale && !is.null(splitx)) {
        if(class(splitx)=="numeric" && length(splitx)==1 ) {
            ##
            
            ii <- Matrix::head(order(-apply(gx,1,sd,na.rm=TRUE)),1000) ## NEED RETHINK!
            system.time( hc <- fastcluster::hclust(as.dist(1 - stats::cor(gx[ii,])), method="ward.D2" ))
            ##system.time( hc <- nclust( as.dist(1 - stats::cor(gx)), link="ward" ))
            splitx = paste0("cluster",cutree(hc, splitx))
            names(splitx) <- colnames(gx)
        }
        if(class(splitx)=="character" && length(splitx)==1 &&
           !is.null(col.annot) && splitx %in% colnames(col.annot) ) {
            splitx = as.character(col.annot[,splitx])
            names(splitx) <- colnames(gx)            
        }        
        for(ux in unique(splitx)) {
            jj <- which(splitx == ux)
            gx[,jj] = gx[,jj,drop=FALSE] - rowMeans(gx[,jj,drop=FALSE],na.rm=TRUE)
        }
    }
      
    ##-------------------------------------------------------------
    ## Take top SD features 
    ##-------------------------------------------------------------
    ##
    if(!is.null(splitx) && length(splitx)==ncol(gx)) {
        names(splitx) <- colnames(gx)        
    }

    jj1 <- 1:nrow(gx)
    jj2 <- 1:ncol(gx)
    if(!is.null(nmax) && nmax < nrow(gx) && nmax>0 ) {
        ##cat("restricting to top",nmax,"maxSD genes\n")
        jj1 <- Matrix::head(order(-apply(gx,1,sd,na.rm=TRUE) ),nmax)
    }
    if(!is.null(cmax) && cmax < ncol(gx) && cmax>0 ) {
        ##cat("restricting to top",cmax,"maxSD samples\n")
        jj2 <- Matrix::head(order(-apply(gx,2,sd,na.rm=TRUE)),cmax)
    }
    gx <- gx[jj1,jj2]
    if(!is.null(row.annot)) {
        row.annot <- row.annot[jj1,,drop=FALSE]
    }
    if(!is.null(col.annot)) {
        col.annot <- col.annot[,which(colMeans(is.na(col.annot))<1),drop=FALSE]
        col.annot <- col.annot[jj2,,drop=FALSE]
        ##col.annot <- type.convert(col.annot)
    }

    ##--------------------------------------------
    ## split columns
    ##--------------------------------------------
    do.splitx = !is.null(splitx)
    do.splitx
    idx2 = NULL
    if(do.splitx && class(splitx)=="numeric" && length(splitx)==1 ) {
        ##
        
        system.time( hc <- fastcluster::hclust( as.dist(1 - stats::cor(gx)), method="ward.D2" ))
        ##system.time( hc <- nclust( as.dist(1 - stats::cor(gx)), link="ward" ))
        idx2 = paste0("cluster",cutree(hc, splitx))
    }
    if(do.splitx && class(splitx)=="character" && length(splitx)==1 &&
       !is.null(col.annot) && splitx %in% colnames(col.annot) ) {
        idx2 = as.character(col.annot[,splitx])
    }
    if(do.splitx && length(splitx)==ncol(gx) ) {
        idx2 = as.character(splitx[colnames(gx)])
    }
    if(!is.null(idx2)) {
        grp <- tapply( colnames(gx), idx2, list)
        ngrp = length(grp)
    } else {
        grp <- list(colnames(gx))
        names(grp)[1] = main  ## title
        ngrp = 1
    }
    remove(idx2)
    names(grp) <- paste0(names(grp),ifelse(names(grp)=='',' ',''))  ## avoid empty names

    ##--------------------------------------------
    ## split rows
    ##--------------------------------------------
    do.split = !is.null(split)
    do.split
    split.idx = NULL
    if(do.split && class(split)=="numeric" && length(split)==1 ) {
        cor.gx <- stats::cor(t(gx),use="pairwise")
        cor.gx[is.na(cor.gx)] <- 0
        hc = fastcluster::hclust( as.dist(1 - cor.gx), method="ward.D2" )
        split.idx = paste0("group",cutree(hc, split))
    }
    if(do.split && class(split)=="character" && length(split)==1 &&
       !is.null(row.annot) && split %in% colnames(row.annot) ) {
        split.idx = row.annot[,split]
    }
    if(do.split && length(split)>1 ) {
        split.idx = split[jj1]
    }

    if(!is.null(row.annot)) cat("[gx.splitmap] 3: duplicated annot.rownames =",
                                sum(duplicated(rownames(row.annot))),"\n")    
    
    if(!is.null(split.idx)) {
        row.annot = cbind(row.annot, cluster = split.idx)
        row.annot = data.frame(row.annot, check.names=FALSE, stringsAsFactors=FALSE)
        colnames(row.annot) = tolower(colnames(row.annot))
        rownames(row.annot) = rownames(gx)
    }
    
    ##--------------------------------------------
    ## PCA denoising 
    ##--------------------------------------------
    if(denoise>0) {
        sv <- irlba::irlba(gx, nv=denoise)
        dn <- dimnames(gx)
        gx <- sv$u %*% diag(sv$d) %*% t(sv$v)
        dimnames(gx) <- dn
    }

    ##--------------------------------------------
    ## column sorting (for what???)
    ##--------------------------------------------
    if(!cluster_columns && !is.null(sort_columns)) {
        y <- col.annot[,sort_columns,drop=FALSE]
        y1 <- as.integer(factor(y))
        jj <- fastcluster::hclust(dist(t(gx)),method="ward.D2")$order
        jj <- jj[order(1e6*y1[jj] + 1:length(jj))]
        gx <- gx[,jj]
        col.annot <- col.annot[jj,]
    }
    
    
    ##-------------------------------------------------------------------------------
    ## column  HeatmapAnnotation objects
    ##-------------------------------------------------------------------------------
    col.ha = vector("list",ngrp)
    if(!is.null(col.annot)) {
        
        col.annot <- col.annot[,which(colMeans(is.na(col.annot))<1),drop=FALSE]
        col.pars = lapply(data.frame(col.annot),function(x) sort(unique(x[!is.na(x)])))
        ##npar = apply(col.annot,2,function(x) length(unique(x[!is.na(x)])))
        npar <- sapply(col.pars, length)
        par.type <- sapply(col.annot,class)        
        col.colors = list()
        i=1
        for(i in 1:length(npar)) {
            prm = colnames(col.annot)[i]
            klrs = rev(grey.colors(npar[i],start=0.4,end=0.85))
            if(npar[i]==1) klrs = "#E6E6E6"
            names(klrs) <- col.pars[[i]]
            x <- col.annot[,i]
            is.number <- (par.type[i] %in% c("integer","numeric"))
            if(npar[i]>3 & !is.number) {
                klrs = rep(RColorBrewer::brewer.pal(8,"Set2"),99)[1:npar[i]]
                names(klrs) <- col.pars[[i]]
                klrs = klrs[!is.na(names(klrs))]
            }
            if(any(is.na(x))) klrs = c(klrs, "NA"="grey90")
            if(is.number) {
                x <- as.numeric(as.character(x))
                rr <- range(x,na.rm=TRUE)
                rr <- rr + c(-1,1)*1e-8
                klrs <- circlize::colorRamp2(rr, c("grey85", "grey40"))
            }            
            ##if(npar[i]==2) klrs = rep(RColorBrewer::brewer.pal(2,"Paired"),99)[1:npar[i]]
            ##names(klrs) = sort(unique(col.annot[,i]))
            col.colors[[prm]] = klrs
        }

        i=1
        for(i in 1:ngrp) {
            jj = grp[[i]]
            ##ap <- list(labels_gp=gpar(fontsize=6*cexRow))
            ap <- list( title_gp = grid::gpar(fontsize=3.5*annot.ht),
                       labels_gp = grid::gpar(fontsize=3.1*annot.ht),
                       grid_width = grid::unit(1*annot.ht, "mm"),
                       grid_height = grid::unit(1*annot.ht, "mm"))
            aa <- rep(list(ap), ncol(col.annot))
            names(aa) <- colnames(col.annot)
            col.ha[[i]] = ComplexHeatmap::HeatmapAnnotation(
                df = col.annot[jj,,drop=FALSE],
                col = col.colors,  na_col='#FCFCFC',
                ##annotation_height = grid::unit(annot.ht, "mm"),
                simple_anno_size = grid::unit(0.85*annot.ht,"mm"),  ## BioC 3.8!!
                show_annotation_name = (i==ngrp),
                show_legend = show_legend & (npar <= 20),
                annotation_name_gp = grid::gpar(fontsize=3.1*annot.ht),
                annotation_legend_param = aa
            )
        }
    }
    
    ## row annotation bars
    row.ha = NULL
    if(!is.null(row.annot)) {
                
        npar = apply(row.annot,2,function(x) length(setdiff(x,c(NA,"NA"))))        
        row.colors = list()
        i=1
        for(i in 1:length(npar)) {
            prm = colnames(row.annot)[i]
            x = row.annot[,i]
            klrs = rev(grey.colors(npar[i],start=0.3, end=0.8))
            if(npar[i]==1) klrs = "#E6E6E6"
            if(npar[i]>0) klrs = rep(RColorBrewer::brewer.pal(8,"Set2"),99)[1:npar[i]]
            ##if(npar[i]==2) klrs = rep(RColorBrewer::brewer.pal(2,"Paired"),99)[1:npar[i]]
            names(klrs) = sort(setdiff(unique(x),NA))
            if(any(is.na(x))) klrs = c(klrs, "NA"="grey90")
            row.colors[[prm]] = klrs
        }
        ##row.ha = ComplexHeatmap::HeatmapAnnotation(
        ##row.ha = ComplexHeatmap::Heatmap(
        
        row.ha = ComplexHeatmap::rowAnnotation(
            df = row.annot,
            col = row.colors,
            ##show_annotation_name = TRUE,
            show_annotation_name = show_colnames,
            show_legend = FALSE,
            ##annotation_name_gp = grid::gpar(fontsize=9*cexRow),
            annotation_name_gp = grid::gpar(fontsize=3.3*annot.ht),
            simple_anno_size = grid::unit(annot.ht,"mm"),  ## BioC 3.8!!
            width = grid::unit(annot.ht*ncol(row.annot),"mm")
        )
    }
    
    ##-------------------------------------------------------------
    ## Plotting methods
    ##-------------------------------------------------------------    

    if(softmax) {
        gx <- tanh(0.5 * gx / sd(gx))
    }
    
    ## ------------- colorscale options
    col_scale <- NULL
    if(!is.null(zlim)) {
        if(length(zlim)==3) {
            zz <- c(-zlim[1], zlim[2], zlim[3])
        } else {
            zmean = mean(gx[,],na.rm=TRUE)
            zz <- c(zlim[1], zmean, zlim[2])
            if(zlim[1]<0) zz <- c(zlim[1], 0, zlim[2])
        }
        col_scale = circlize::colorRamp2(zz, c("royalblue3","grey90","indianred3"))
    } else if(symm.scale) {
        colmax = 1
        colmax = max(abs(gx[,]),na.rm=TRUE)
        col_scale = circlize::colorRamp2(c(-colmax, 0, colmax),
                                         c("royalblue3","grey90","indianred3"))
    } else {
        colmin = min(gx[,],na.rm=TRUE)
        colmax = max(gx[,],na.rm=TRUE)
        colmean = mean(gx[,],na.rm=TRUE)
        col_scale = circlize::colorRamp2(c(colmin, colmean, colmax),
                                         c("royalblue3","grey90","indianred3"))
    }
    
    ## ------------- cluster blocks
    grp.order <- 1:ngrp
    if(!is.null(order.groups) && ngrp>1 && order.groups[1] == "clust") {
        ## Reorder cluster indices based on similarity clustering
        mx <- do.call(cbind, lapply(grp, function(i) rowMeans(gx[,i,drop=FALSE])))
        ##mx <- Matrix::head(mx[order(-apply(mx,1,sd)),],100)
        mx <- t(scale(t(mx)))
        grp.order <- fastcluster::hclust(dist(t(mx)))$order
    }
    if(!is.null(order.groups) && ngrp>1 && length(order.groups) == ngrp) {
        grp.order <- match(order.groups, names(grp))
    }

    ## ------------- draw heatmap
    ##split.idx = NULL
    ##if(!is.null(split) && split>0) split.idx = row.annot[rownames(gx),split,drop=FALSE]
    hmap = NULL
    for(i in grp.order) {
        jj <- grp[[i]]

        coldistfun1 <- function(x) dist(x)
        rowdistfun1 <- function(x,y) 1 - stats::cor(x, y)
        gx0 <- gx[,jj,drop=FALSE]
        grp.title = names(grp)[i]
        ##grp.title = shortstring(names(grp)[i],15)

        gx0.colnames <- NULL
        if(show_colnames) {
            gx0.colnames <- ComplexHeatmap::HeatmapAnnotation(
                text = ComplexHeatmap::anno_text( colnames(gx0), rot=column_names_rot,
                                 gp = grid::gpar(fontsize=11*cexCol),                                 
                                 location = grid::unit(1, "npc"), just = "right"),                
                annotation_height = ComplexHeatmap::max_text_width(
                    colnames(gx0), grid::gpar(fontsize=11*cexCol)) * sin((column_names_rot/180)*pi)
            )
        }
        
        hmap = hmap + ComplexHeatmap::Heatmap( gx0,
                              ##col = col,  ## from input
                              col = col_scale,  ## from input
                              cluster_rows=cluster_rows,
                              cluster_columns=cluster_columns,
                              ##cluster_columns = as.dendrogram(h1),
                              ##cluster_rows = as.dendrogram(h2),
                              clustering_distance_rows = dist.method,
                              ##clustering_distance_columns = "euclidean",
                              clustering_distance_columns = col.dist.method,
                              clustering_method_rows = "ward.D2",
                              ##clustering_method_rows = "average",
                              clustering_method_columns = "ward.D2",
                              show_heatmap_legend = FALSE,
                              split = split.idx,
                              row_title_gp = grid::gpar(fontsize = 11),  ## cluster titles
                              row_names_gp = grid::gpar(fontsize = 10*cexRow),
                              top_annotation = col.ha[[i]],
                              ## top_annotation_height = grid::unit(annot.ht*ncol(col.annot), "mm"),
                              column_title = grp.title,
                              name = names(grp)[i],
                              column_title_rot = column_title_rot,
                              column_title_gp = grid::gpar(fontsize = 11*title_cex, fontface='plain'),
                              ##column_title_gp = grid::gpar(fontsize = 11*title_cex, fontface='bold')
                              show_row_names=FALSE,
                              ##show_column_names=show_colnames,
                              show_column_names = FALSE,
                              bottom_annotation = gx0.colnames,
                              column_names_gp = grid::gpar(fontsize = 11*cexCol)
                              )
    }
    
    rownames.ha <- NULL
    if(FALSE && show_rownames < nrow(gx) && show_rownames>0 ) {
        ## Show rownames with linked lines
        ## !!!!!!! BUGGY!!!!!!!!!!!!!
        if(1 && !is.null(split.idx) && length(unique(split.idx))>1) {
            nshow=10
            nshow <- show_rownames / length(unique(split.idx))
            nshow <- max(nshow,10)
            subidx <- tapply(1:length(split.idx), split.idx, function(ii)
                ii[head(order(-apply(gx[ii,,drop=FALSE],1,sd,na.rm=TRUE)),nshow)] )
            subset <- unlist(subidx)
        } else {
            subset = Matrix::head(order(-apply(gx,1,sd,na.rm=TRUE)), show_rownames )
        }
        lab = rownames(gx)[subset]
        lab <- substring(lab,1,lab.len)
        rownames.ha = ComplexHeatmap::rowAnnotation(
            ##link = row_anno_link(at = subset, labels=lab,
            link = ComplexHeatmap::anno_mark(at = subset, labels=lab,
                             link_width = grid::unit(0.8,"cm"),
                             labels_gp = grid::gpar(fontsize = 10*cexRow)),
            width = grid::unit(0.8,"cm") + 0.8*max_text_width(lab))
    }
   
    if(1 && is.null(rownames.ha) && show_rownames>0 ) {
        ## empty matrix just for rownames on the far right
        empty.mat = matrix(nrow = nrow(gx), ncol = 0)
        rowlab = substring(trimws(rownames(gx)),1,rowlab.maxlen)
        ## rowlab = paste0(rowlab,'............')
        rownames(empty.mat) = rowlab
        if(!is.null(col.annot)) rowlab <- c(rowlab, colnames(col.annot))
        textwidth <- ComplexHeatmap::max_text_width(
            paste0(rowlab,"XXXXXXXXXX"),
            gp = grid::gpar(fontsize = 10*cexRow))        
        textwidth <- grid::unit(rownames_width,'mm')        
        rownames.ha = ComplexHeatmap::Heatmap(
            empty.mat,
            row_names_max_width = textwidth,
            row_names_gp = grid::gpar(fontsize = 10*cexRow),
            show_row_names = show_rownames>0 )
    }

    if(!is.null(row.ha)) {
        ##message("[gx.splitmap] rendering row annotation")        
        hmap = hmap + row.ha
    }
    if(!is.null(rownames.ha) && show_rownames) {
        ##message("[gx.splitmap] rendering rownames")
        hmap = hmap + rownames.ha
    }
    
    ##draw(hmap, annotation_legend_side="right")
    ComplexHeatmap::draw(hmap, annotation_legend_side="right",
         padding = grid::unit(mar, "mm"), gap = grid::unit(1.0,"mm") )
    ##hmap

    if(show_key && !is.null(col_scale)) {
        brk <- attr(col_scale,"breaks")
        lgd = ComplexHeatmap::Legend(col_fun = col_scale, at=brk[c(1,3)], labels=round(brk[c(1,3)],2),
                     border = "transparent", labels_gp = grid::gpar(fontsize = 8), 
                     legend_width = grid::unit(0.08, "npc"), legend_height = grid::unit(0.01, "npc"),
                     title = "\n", direction = "horizontal")
        ##key.offset=c(0.05,1.01)
        ComplexHeatmap::draw(lgd, x=grid::unit(key.offset[1], "npc"), y=grid::unit(key.offset[2], "npc"), just=c("left", "top"))
    }
    
    ##res <- c()
    ##res$col.clust <- h1
    ##res$row.clust <- h2
    ##invisible(res)
}

##gx=X
gx.heatmap <- function(gx, values=NULL,
                       clust.method="ward.D2",
                       dist.method="pearson",
                       col.dist.method="euclidean",
                       plot.method="heatmap.2",
                       col = BLUERED64, softmax=FALSE,
                       ##col = gplots::colorpanel(64,"blue","grey90","red"),
                       scale="row", verbose=1, symm=FALSE,
                       ## Rowv = NA, Colv = NA,
                       col.annot=NULL, row.annot=NULL, annot.ht=1,
                       nmax=1000, cmax=NULL, show_colnames=TRUE,
                       indent.names=FALSE,
                       ... )
{

    if(0) {
        clust.method="ward.D2"; dist.method="pearson"; col.dist.method="euclidean";
        plot.method="heatmap.2";symm=FALSE;values=NULL;softmax=FALSE;values=NULL
        col=colorpanel(64,"blue","grey90","red"); scale="row";verbose=3;show_colnames=TRUE
        ## Rowv = NA, Colv = NA
        col.annot=NULL; row.annot=NULL; nmax=1000; cmax=NULL; indent.names=FALSE
    }

    
    ##
    par(xpd=FALSE)

    if(verbose>1) cat("input.dim.gx=",dim(gx),"\n")
    
    ##-------------------------------------------------------------
    ## scaling options
    ##-------------------------------------------------------------    
    if(verbose>1) cat("dim.gx=",dim(gx),"\n")
    sym0=FALSE
    if(sum(gx<0,na.rm=TRUE)>0 || scale %in% c("row","col")) sym0=TRUE
    ## scaling options
    if("col" %in% scale || "both" %in% scale) {
        ##gx = scale(gx)
        tgx <- t(gx) - colMeans(gx,na.rm=TRUE)
        gx <- t(tgx / (1e-4 + apply(gx, 2, sd))) ## small eps maintains SD order!
        remove(tgx)
    }
    if("row" %in% scale || "both" %in% scale) {
        ##gx = t(scale(t(gx)))
        gx <- gx - rowMeans(gx,na.rm=TRUE)
        gx <- gx / (1e-4 + apply(gx, 1, sd)) ## small eps maintains SD order!
    }
    if("col.center" %in% scale) {
        gx  <- t( t(gx) - colMeans(gx, na.rm=TRUE))
    }
    if("row.center" %in% scale) {
        gx  <- gx - rowMeans(gx,na.rm=TRUE)
    }

    ##-------------------------------------------------------------
    ## scaling options
    ##-------------------------------------------------------------    
    jj1 <- 1:nrow(gx)
    jj2 <- 1:ncol(gx)
    if(!is.null(nmax) && nmax < nrow(gx) && nmax>0 ) {
        ##cat("restricting to top",nmax,"maxSD genes\n")
        jj1 <- Matrix::head(order(-apply(gx,1,sd,na.rm=TRUE) ),nmax)
    }
    if(!is.null(cmax) && cmax < ncol(gx) && cmax>0 ) {
        ##cat("restricting to top",cmax,"maxSD samples\n")
        jj2 <- Matrix::head(order(-apply(gx,2,sd,na.rm=TRUE)),cmax)
    }
    if(symm && ncol(gx)==nrow(gx)) {
        jj2 <- jj1
    }
    gx <- gx[jj1,jj2]
    if(!is.null(col.annot)) col.annot <- col.annot[jj2,,drop=FALSE]
    if(!is.null(row.annot)) row.annot <- row.annot[jj1,,drop=FALSE]
            
    ##------------------------------------------------------------
    ## Distance calculation
    ##------------------------------------------------------------
    fillNA <- function(x) {
        nx <- x
        nx[is.na(nx)] <- 0
        nx <- x + is.na(x)*rowMeans(x,na.rm=TRUE)
        nx
    }

    ## sample dimension (column)
    h1 <- NULL
    if(!is.null(clust.method) ) {
        if(col.dist.method=="pearson") {
            suppressWarnings( cx <- stats::cor(gx,use="pairwise.complete.obs") )
            cx[is.na(cx)] <- 0 ## really???
            d1 <- as.dist(1 - cx)
        } else if (col.dist.method %in% c("multi.dist","multidist")) {
            gx0 <- t(gx)
            if(any(is.na(gx0))) gx0 <- fillNA(gx0)
            d1 <- multi.dist(gx0)
        } else  {
            ## euclidean
            gx0 <- t(gx)
            if(any(is.na(gx0))) gx0 <- fillNA(gx0)
            d1 <- dist(gx0, method=col.dist.method)
        }
        h1 <- fastcluster::hclust(d1, method=clust.method)
        ##h1 <- as.dendrogram(h1)
    }

    ## gene dimension (rows)
    h2 <- NULL
    if(symm && ncol(gx)==nrow(gx)) {
        h2 <- h1
    } else if(!is.null(clust.method) ) {
        if(dist.method=="pearson") {
            suppressWarnings( cx <- stats::cor(t(gx),use="pairwise.complete.obs") )
            cx[is.na(cx)] <- 0
            d2 <- as.dist(1 - cx)
        } else if (dist.method %in% c("multi.dist","multidist")) {
            gx0 <- gx
            if(any(is.na(gx0))) gx0 <- fillNA(gx0)
            d2 <- multi.dist(gx0)
        } else  {
            gx0 <- gx
            if(any(is.na(gx0))) gx0 <- fillNA(gx0)
            d2 <- dist(gx0, method=dist.method)
        }
        h2 <- fastcluster::hclust( d2 , method=clust.method)
        ##h2 <- as.dendrogram(h2)
    }
    dd <- c("both","row","column","none")[1 + 1*is.null(h1) + 2*is.null(h2)]

    ##------------------------------------------------------------
    ## annotation bars
    ##------------------------------------------------------------    
    ## ident names??
    if(indent.names > 0) {
        nn.sp <- sapply(floor((1:nrow(gx))/indent.names),function(n) paste(rep(" ",n),collapse=""))
        rownames(gx) <- paste(nn.sp,rownames(gx))
    }

    cc0 <- NA
    cc1 <- NA
    if(!is.null(col.annot)) {
        plot.method="heatmap.3"
        ##ry <- apply(col.annot,2,rank,na.last="keep")
        col.annot = col.annot[,which(colMeans(is.na(col.annot))<1),drop=FALSE]
        col.annot = as.data.frame(col.annot)
        aa <- col.annot
        is.num = (sapply(aa,class)=="numeric")
        if(sum(is.num)) {
            for(j in which(is.num)) aa[,j] = cut(aa[,j],breaks=16)
        }
        ry <- sapply(aa,function(x) as.integer(as.factor(x)))
        ry <- t(ry)
        ry <- ry - apply(ry,1,min,na.rm=TRUE)
        ry <- round( (ry / (1e-8+apply(ry,1,max,na.rm=TRUE))  )*7 )
        ##cc0 <- matrix(c("white","black")[ry+1],nrow=nrow(ry),ncol=ncol(ry))
        ry = t(apply(ry, 1, function(x) as.integer(as.factor(x))))
        ##cc0 <- matrix([ry],nrow=nrow(ry),ncol=ncol(ry))
        cc0 = t(apply(ry,1,function(x) rev(grey.colors(max(x,na.rm=TRUE)))[x]))
        jj = which(apply(cc0,1,function(x) length(unique(x))) > 3 & !is.num)
        if(length(jj)) {
            
            klrs = rep(RColorBrewer::brewer.pal(8,"Set2"),99)
            klrs.mat <- matrix(klrs[ry[jj,,drop=FALSE]+1], nrow=length(jj))
            cc0[jj,] <- klrs.mat
        }
        ##cc0 = cc0[which(rowMeans(!is.na(cc0))>0),,drop=FALSE ]
        cc0 <- t(cc0)
        colnames(cc0) <- colnames(col.annot)
        rownames(cc0) <- rownames(col.annot)
    }
    if(!is.null(row.annot)) {
        plot.method="heatmap.3"
        ##ry <- t(apply(row.annot,2,rank,na.last="keep"))
        ry <- sapply(as.data.frame(row.annot),function(x) as.integer(as.factor(x)))
        ry <- t(ry)
        ry <- ry - apply(ry,1,min,na.rm=TRUE)
        ry <- round( (ry / (1e-8+apply(ry,1,max,na.rm=TRUE))  )*7 )
        ry = t(apply(ry, 1, function(x) as.integer(as.factor(x))))
        ##cc0 <- matrix([ry],nrow=nrow(ry),ncol=ncol(ry))
        cc1 = t(apply(ry,1,function(x) rev(grey.colors(max(x,na.rm=TRUE)))[x]))
        jj = which(apply(cc1,1,function(x) length(unique(x))) > 3)
        if(length(jj)) {
            
            klrs = rep(RColorBrewer::brewer.pal(8,"Set2"),99)
            cc1[jj,,drop=FALSE] <- matrix(klrs[ry[jj,,drop=FALSE]+1], nrow=length(jj))
        }
        rownames(cc1) <- colnames(row.annot)
        colnames(cc1) <- rownames(row.annot)
    }
    
    ##-------------------------------------------------------------
    ## Plotting methods
    ##-------------------------------------------------------------    
    if(!is.null(values)) gx <- values[rownames(gx),colnames(gx)]
    if(softmax) gx <- tanh(0.5* gx / sd(gx))  ## just for plotting...
    ##if(all(gx>=0)) col <- Matrix::tail(col,length(col)/2)
    ##if(all(gx<=0)) col <- Matrix::head(col,length(col)/2)
    side.height <- 0.1 * annot.ht * NCOL(cc0)

    if(!show_colnames) {
        colnames(gx) <- rep("",ncol(gx))
    }

    symm
    if(0 && symm && !is.null(cc0)) {
        cc1 <- t(cc0)
    }
    if(!is.null(cc0) && !is.na(cc0)) cc0 <- cc0[colnames(gx),,drop=FALSE]
    if(!is.null(cc1) && !is.na(cc1)) cc1 <- cc1[,rownames(gx),drop=FALSE]

    ## draw heatmap    
    if(plot.method=="heatmap.3" && !is.na(cc0) && !is.na(cc1) ) {
        if(verbose>1) cat("plotting with heatmap.3 + both ColSideColors\n")
        if(is.null(h1) && is.null(h2)) {
            heatmap.3(gx, Colv=NULL, Rowv=NULL,
                      dendrogram=dd, col=col, scale="none",
                      symkey=sym0, symbreaks=sym0, trace="none",
                      side.height.fraction=side.height,
                      ColSideColors=cc0,
                      RowSideColors=cc1,
                      ...)
        } else {
            heatmap.3(gx, Colv=as.dendrogram(h1), Rowv=as.dendrogram(h2),
                      dendrogram = dd, col = col, scale="none",
                      symkey = sym0, symbreaks = sym0, trace="none",
                      side.height.fraction = side.height,
                      ColSideColors = cc0,
                      ##RowSideColors = matrix(cc0[,1],nrow=1) )
                      RowSideColors = cc1,
                      ...)
        }
    } else if(plot.method=="heatmap.3" && !is.na(cc0) && is.na(cc1) ) {
        if(verbose>1) cat("plotting with heatmap.3 + ColSideColors\n")
        if(is.null(h1) && is.null(h2)) {
            heatmap.3(gx, Colv=NULL, Rowv=NULL,
                      dendrogram=dd, col=col, scale="none",
                      symkey=sym0, symbreaks=sym0, trace="none",
                      side.height.fraction=side.height,
                      ColSideColors=cc0, ## RowSideColors=cc1,
                      ...)
        } else {
            heatmap.3(gx, Colv=as.dendrogram(h1), Rowv=as.dendrogram(h2),
                      dendrogram=dd, col=col, scale="none",
                      symkey=sym0, symbreaks=sym0, trace="none",
                      side.height.fraction=side.height,
                      ColSideColors=cc0, ## RowSideColors=cc1,
                      ...)
        }
    } else if(plot.method=="heatmap.3" && is.na(cc0) && !is.na(cc1) ) {
        if(verbose>1) cat("plotting with heatmap.3 + RowSideColors\n")
        if(is.null(h1) && is.null(h2)) {
            heatmap.3(gx, Colv=NULL, Rowv=NULL,
                      dendrogram=dd, col=col, scale="none",
                      symkey=sym0, symbreaks=sym0, trace="none",
                      side.height.fraction=side.height,
                      ## ColSideColors=cc0,
                      RowSideColors=cc1,
                      ...)
        } else {
            heatmap.3(gx, Colv=as.dendrogram(h1), Rowv=as.dendrogram(h2),
                      dendrogram=dd, col=col, scale=scale,
                      symkey=sym0, symbreaks=sym0, trace="none",
                      side.height.fraction=side.height,
                      ##ColSideColors=cc0,
                      RowSideColors=cc1,
                      ...)
        }
    } else if(plot.method=="heatmap.3" && is.na(cc0) && is.na(cc1) ) {
        if(verbose>1) cat("plotting with heatmap.3 no ColSideColors\n")
        if(is.null(h1) && is.null(h2)) {
            heatmap.3(gx, Colv=NULL, Rowv=NULL,
                      dendrogram=dd, col=col, scale="none",
                      symkey=sym0, symbreaks=sym0, trace="none",
                      side.height.fraction=side.height,
                      ...)
        } else {
            heatmap.3(gx, Colv=as.dendrogram(h1), Rowv=as.dendrogram(h2),
                      dendrogram=dd, col=col, scale="none",
                      symkey=sym0, symbreaks=sym0, trace="none",
                      side.height.fraction=side.height,
                      ...)
        }
    
    }  else {
        if(verbose>1) cat("plotting with heatmap.2\n")
        if(is.null(h1) && is.null(h2))  {
            gplots::heatmap.2(gx, Colv=NULL, Rowv=NULL,
                      dendrogram=dd, col=col, scale="none",
                      ##side.height.fraction=side.height,
                      ##symkey=sym, symbreaks=sym, trace="none")
                      symkey=sym0, symbreaks=sym0, trace="none", ...)
        } else {
            gplots::heatmap.2(gx, Colv=as.dendrogram(h1), Rowv=as.dendrogram(h2),
                      dendrogram=dd, col=col, scale="none",
                      ##side.height.fraction=side.height,
                      ##symkey=sym, symbreaks=sym, trace="none")
                      symkey=sym0, symbreaks=sym0, trace="none", ...)
        }
    }

    res <- c()
    res$col.clust <- h1
    res$row.clust <- h2
    invisible(res)
}

##nc=8;nr=8;na=6;p=1;q=0.9;nca=24;method="pearson";na.fill=FALSE;labrow=rownames(x);labcol=colnames(x);nrlab=nrow(x);nclab=ncol(x)
clustermap <- function(x, nc=6, nr=6, na=4, q=0.80, p=2,
                       method="multidist",
                       nca=24, col.annot=NULL, row.annot=NULL, plot=TRUE,
                       na.fill=FALSE, nrlab=nrow(x), nclab=ncol(x),
                       labrow=rownames(x), labcol=colnames(x),
                       ...)
{
    
    
    ## non-linear transformation
    ## x <- abs(x)**p * sign(x)
    ##x[is.na(x)] <- 0

    if(method=="pearson") {
        ##d1 <- as.dist(1 - stats::cor(t(x),use="pairwise")) ## rows
        d1 <- as.dist( stats::cor(t(x),use="pairwise")) ## rows:=drugs
        d2 <- as.dist( stats::cor(x,use="pairwise")) ## cols:=samples
        d1[is.na(d1)] <- mean(d1,na.rm=TRUE)
        d2[is.na(d2)] <- mean(d2,na.rm=TRUE)
        d1 <- (max(d1) - d1)**(p/2)
        d2 <- (max(d2) - d2)**(p/2)
        h1 <- fastcluster::hclust( d1, method="ward.D2" )
        h2 <- fastcluster::hclust( d2, method="ward.D2" )
    } else if(method=="minkowski") {
        d1 <- dist( x, method="minkowski",p=p)
        d2 <- dist( t(x), method="minkowski",p=p)
        d1[is.na(d1)] <- mean(d1,na.rm=TRUE)
        d2[is.na(d2)] <- mean(d2,na.rm=TRUE)
        h1 <- fastcluster::hclust( d1, method="ward.D2" )
        h2 <- fastcluster::hclust( d2, method="ward.D2" )
    } else if(method %in% c("multi.dist","multidist")) {
        d1 <- multi.dist( x )
        d2 <- multi.dist( t(x) )
        d1[is.na(d1)] <- mean(d1,na.rm=TRUE)
        d2[is.na(d2)] <- mean(d2,na.rm=TRUE)
        h1 <- fastcluster::hclust( d1, method="ward.D2" )
        h2 <- fastcluster::hclust( d2, method="ward.D2" )
    } else {
        stop("unknown distance method: ",method)
    }

    ##my.col=colorpanel(64,"blue","grey90","red")
    ##my.col=bluered()
    my.col = BLUERED64
    c1 <- cutree(h1,nr)
    c2 <- cutree(h2,nc)
    kxmap <- function(x, c1, c2, q) {
        nr <- length(unique(c1))
        nc <- length(unique(c2))
        mx <- matrix(0,nr,nc)
        for(i in 1:nr) {
            for(j in 1:nc) {
                xx <- x[which(c1==i),which(c2==j)]
                mx[i,j] <- mean(xx, na.rm=TRUE)
            }
        }
        kx <- x*0
        for(i in 1:nrow(kx)) {
            for(j in 1:ncol(kx)) {
                kx[i,j] <- mx[c1[i],c2[j]]
            }
        }
        x0 <- x
        if(na.fill==TRUE) {
            x0[is.na(x)] <- kx[is.na(x)]  ## fill with cluster mean??
        }
        kx <- q*kx + (1-q)*x0
        return(kx)
    }
    kx <- kxmap(x, c1, c2, q=q)
    dim(kx)

    ## order annotation
    order.annot <- function(A, kx, c1, n) {
        mkx <- apply( kx, 2, function(x) tapply(x,c1,mean))
        rho1 <- stats::cor(t(A),t(mkx),use="pairwise")
        r <- apply(-abs(rho1),2,order)
        jj <- unique(as.vector(t(r)))
        A <- A[jj,,drop=FALSE]
        A <- A[!duplicated(A),,drop=FALSE]
        Matrix::head(A,n)
    }

    ## create annotation side bars
    cc0 <- matrix("white",nrow=ncol(kx),ncol=1)
    cc1 <- matrix("white",nrow=1,ncol=nrow(kx))
    rownames(cc0) <- 1:ncol(kx)
    colnames(cc0) <- " "
    rownames(cc1) <- " "
    colnames(cc1) <- 1:nrow(kx)
    if(!is.null(col.annot)) {
        col.annot0 <- col.annot
        col.annot <- order.annot(col.annot, kx, c1, n=nca)
        ax <- kxmap(col.annot, 1:nrow(col.annot),c2,q=q)
        col.annot <- (col.annot - rowMeans(col.annot,na.rm=TRUE)) /
            (1e-8+apply(col.annot,1,sd,na.rm=TRUE))
        d3 <- dist(col.annot)
        d3 <- as.dist(1 - stats::cor(t(col.annot),use="pairwise"))
        d3[is.na(d3)] <- mean(d3,na.rm=TRUE) ## impute missing...
        h3 <- fastcluster::hclust(d3,method="ward.D2")
        ##h3 <- corclust(ax)
        na0 <- min(na,dim(col.annot))
        col.annot <- kxmap(col.annot, cutree(h3,na0), c2, q=q)
        col.annot <- col.annot[h3$order,,drop=FALSE]
        ## Matrix::image( col.annot[,h2$order], col=my.col)
        ry <- t(apply(col.annot,1,rank,na.last="keep"))
        ry <- ry - apply(ry,1,min,na.rm=TRUE)
        ry <- (ry / (apply(ry,1,max,na.rm=TRUE) + 1e-8))*31
        ##cc0 <- matrix(c("white","black")[ry+1],nrow=nrow(ry),ncol=ncol(ry))
        cc0 <- matrix(rev(grey.colors(32))[ry+1],nrow=nrow(ry),ncol=ncol(ry))
        ##cc0 <- matrix(heat.colors(32)[ry+1],nrow=nrow(ry),ncol=ncol(ry))
        colnames(cc0) <- colnames(col.annot)
        rownames(cc0) <- rownames(col.annot)
        cc0 <- t(cc0)
    }
    if(FALSE && !is.null(row.annot)) {
        jj <- fastcluster::hclust(dist(t(row.annot)))$order
        row.annot <- row.annot[,jj]
        ry <- apply(row.annot,2,rank,na.last="keep")
        ry <- t(t(ry) - apply(ry,2,min,na.rm=TRUE))
        ry <- t(t(ry) / (apply(ry,2,max,na.rm=TRUE) + 1e-8)) * 15
        ##cc0 <- matrix(c("white","black")[ry+1],nrow=nrow(ry),ncol=ncol(ry))
        cc1 <- matrix(rev(grey.colors(16))[ry+1],nrow=nrow(ry),ncol=ncol(ry))
        dim(cc1)
        colnames(cc1) <- NULL
        rownames(cc1) <- NULL
        colnames(cc1) <- colnames(row.annot)
        ##rownames(cc1) <- rownames(row.annot)
        cc1 <- t(cc1)
    }

    if(plot==TRUE) {
        j1 <- nrow(x) - which(diff(c1[h1$order])!=0)
        j2 <- which(diff(c2[h2$order])!=0)
        nh <- max(ncol(cc0),nrow(cc1))**0.5
        nh
        dim(kx)
        dim(cc0)
        dim(cc1)
        ##labrow=rownames(x);labcol=colnames(x)
        ##nrlab=nclab=2
        k1 <- ((1:nrow(x) %% max(1,floor(nrow(x)/nrlab)))==0)
        k2 <- ((1:ncol(x) %% max(1,floor(ncol(x)/nclab)))==0)
        if(sum(!k1)>0) {
            labrow[h1$order[!k1]] <- ""
            labrow[h1$order[k1]]  <- paste(labrow[h1$order[k1]],"+ . . .")
        }
        if(sum(!k2)>0) {
            labcol[h2$order[!k2]] <- ""
            labcol[h2$order[k2]] <- paste(". . . +",labcol[h2$order[k2]])
        }
        heatmap.3(kx, Colv=as.dendrogram(h2), Rowv=as.dendrogram(h1),
                  col=my.col, labRow=labrow, labCol=labcol,
                  ColSideColors=cc0, ...)
        ##side.height.fraction=0.1*nh )
        ##dev.off()
    }

    ## return order clustermap
    kx <- kx[h1$order, h2$order]
    invisible(kx)
}

##n=m=8
frozenmap <- function(x, m=8, n=8, ...) {
    ##n=100;k=3
    x.downsample <- function(x, n=100, k=1, dist.method="pearson") {
        if(n >= nrow(x)) return(x)
        if(dist.method=="multi.dist") {
            d <- multi.dist(x)
        } else if(dist.method=="pearson") {
            d <- as.dist( 1 - stats::cor(t(x),use="pairwise") )
            d[is.na(d)] <- mean(d,na.rm=TRUE)
        } else {
            d <- dist(x)
        }
        h <- fastcluster::hclust(d)
        idx <- cutree(h, k=n)
        ordx <- -apply(h$merge,1,min)
        ordx <- ordx[ordx>0]
        sdx <- apply(x,1,sd,na.rm=TRUE)
        midx  <- function(j) j[which.min(match(j,ordx))]
        maxsd  <- function(j) j[which.max(sdx[j])]
        maxsd2  <- function(j,k) Matrix::head(j[order(-sdx[j])],k)
        ii <- sapply(1:n, function(i) which(idx==i))
        ##jj <- sapply(1:n, function(i) midx(which(idx==i)))
        ##jj <- sapply(1:n, function(i) maxsd(which(idx==i)))
        jj <- lapply(1:n, function(i) maxsd2(which(idx==i),k))
        nn <- (table(idx) > sapply(jj,length))
        hx <- t(sapply(jj, function(j) colMeans(x[j,,drop=FALSE],na.rm=TRUE)))
        ##hx <- t(sapply(ii, function(j) colMeans(x[j,,drop=FALSE],na.rm=TRUE)))
        hx[is.nan(hx)] <- NA
        rn <- sapply(jj, function(j) paste(rownames(x)[j],collapse=" | "))
        rownames(hx) <- paste(rn,c("","+ ...")[1+nn])
        hx
    }

    cx <- x.downsample(x,n=m)
    cx <- t( x.downsample(t(cx), n=n) )
    res <- gx.heatmap(cx, ...)
    cx <- cx[ res$row.clust$order, res$col.clust$order]
    invisible(cx)
}

multi.dist <- function(x, p=4, method=c("pearson","euclidean","manhattan"))
{
    ##method=c("pearson","euclidean","manhattan")
    mm <- matrix(Inf, nrow(x), nrow(x))
    sx <- t(scale(t(x)))
    rx <- t(apply(x,1,rank,na.last="keep"))
    xx <- list(x, sx, rx)  ## parameter transform
    for(this.x in xx) {
        for(m in method) {
            if(m=="pearson") {
                d1 <- as.dist( 1 - stats::cor(t(this.x),use="pairwise"))
            } else {
                d1 <- dist( this.x, method=m, p=p)
            }
            m1 <- matrix(rank(as.vector(as.matrix(d1)),na.last="keep"),nrow=nrow(x))
            m1 <- m1 / max(m1,na.rm=TRUE)
            diag(m1) <- 0
            mm <- pmin(mm, m1, na.rm=TRUE)
        }
    }
    ## mm[is.na(mm)] <- mean(mm,na.rm=TRUE)
    colnames(mm) <- rownames(x)
    rownames(mm) <- rownames(x)
    as.dist(mm)
}


##########################################################################
## from http://www.biostars.org/p/18211/
## https://gist.github.com/nachocab/3853004
##########################################################################

## EXAMPLE USAGE
if(0) {
  ## example of colsidecolors rowsidecolors (single column, single row)
  mat <- matrix(1:100, byrow=T, nrow=10)
  column_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
  column_annotation <- as.matrix(column_annotation)
  colnames(column_annotation) <- c("Variable X")

  row_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
  row_annotation <- as.matrix(t(row_annotation))
  rownames(row_annotation) <- c("Variable Y")

  heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)

  ## multiple column and row
  mat <- matrix(1:100, byrow=T, nrow=10)
  column_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), ncol=2)
  colnames(column_annotation) <- c("Variable X1", "Variable X2")

  row_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), nrow=2)
  rownames(row_annotation) <- c("Variable Y1", "Variable Y2")

  heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)
}

# CODE (see https://gist.github.com/nachocab/3853004)

heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){

    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }

    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dendrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)

        if (!missing(ColSideColors)) {
           #if (!is.matrix(ColSideColors))
           #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
             lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
        }

        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }

    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    plotly::layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                Matrix::image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            Matrix::image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(colnames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }

    if (!missing(ColSideColors)) {

        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            Matrix::image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            Matrix::image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }

    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    Matrix::image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load 
        mmat <- ifelse(is.na(x), 1, NA)
        Matrix::image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }

        z <- seq(min.raw, max.raw, length = length(col))
        Matrix::image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}
