##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

########################################################################
## Plotting functions
########################################################################

if(0) {

    fc <- pgx.getMetaMatrix(ngs, level='geneset')$fc
    x <- Matrix::head( fc[order(-rowMeans(fc**2)),], 60 )
    
    par(mar=c(20,4,4,2), mfrow=c(1,1))
    barplot( t(x), beside=FALSE, las=3)
    ##par(mgp=c(2,1,0))
    pgx.stackedBarplot(x, ylab="cumulative logFC", cex.names=0.001, srt=60, adj=1) 

    par(mar=c(4,0,4,2), mfrow=c(1,2)); frame()
    pgx.stackedBarplot(Matrix::head(x,40), xlab="cumulative logFC", hz=TRUE, cex.names=0.8) 
    
    x=zx0;dim=2;method="pca"
}

heatmapWithAnnot <- function(F, anno.type=c('boxplot','barplot'),
                             bar.height=NULL, map.height=NULL,
                             row_fontsize=9, column_fontsize=9,
                             inset = c(-0.025,-0.1),
                             mar = c(0,0,0,0), legend=TRUE,
                             ...
                               )
{
    col1 <- 1:ncol(F)
    anno.type = anno.type[1]
    if(!is.null(map.height)) map.height <- grid::unit(map.height,'cm')
    if(!is.null(bar.height)) bar.height <- grid::unit(bar.height,'cm')    

    if(anno.type=='barplot') {
        ha <- ComplexHeatmap::HeatmapAnnotation(
            up = ComplexHeatmap::anno_barplot(pmax(F,0.01),gp=gpar(fill=col1)),
            down = ComplexHeatmap::anno_barplot(pmin(F,-0.01),gp=gpar(fill=col1)),
            annotation_height = c(0.5,0.5)*bar.height,
            annotation_label = c("", "")
        )
    } else {
        ha <- ComplexHeatmap::HeatmapAnnotation(
            up = ComplexHeatmap::anno_boxplot(
                t(F),
                gp=gpar(fill='grey80'),
                outline=TRUE, size=unit(1, "mm"), 
                box_width = 0.75),            
            annotation_height = bar.height,
            annotation_label = c("")
        )
    }

    ht <- ComplexHeatmap::Heatmap(
        t(F), name = "logFC",
        top_annotation = ha,
        row_names_gp = grid::gpar(fontsize=row_fontsize),
        column_names_gp = grid::gpar(fontsize=column_fontsize),
        height = map.height,
        heatmap_legend_param = list(
            legend_height = grid::unit(8, "mm")
        )
      , ...
    )
    
    ComplexHeatmap::draw(ht, padding = grid::unit(mar, "cm"),
         heatmap_legend_side = "left")

    if(anno.type=='boxplot') {
        ComplexHeatmap::decorate_annotation("up", { 
            grid::grid.text("logFC",
                      x = grid::unit(0,"npc") - grid::unit(7,"mm"),
                      y = grid::unit(0.5,"npc") - grid::unit(0,"mm"),
                      gp = grid::gpar(fontsize=10),
                      rot=90, hjust = 'center') 
        })
    }
    if(anno.type=='barplot') {
        ComplexHeatmap::decorate_annotation("up", { 
            grid::grid.text("cumulative logFC",
                      x = grid::unit(0,"npc") - grid::unit(7,"mm"),
                      y = grid::unit(0,"npc") - grid::unit(0,"mm"),
                      gp = grid::gpar(fontsize=10),
                      rot=90, hjust = 'center') 
        })
    }
    if(anno.type=='barplot' && legend) {
        legend("topright",colnames(F), fill=col1,
               cex=0.63, y.intersp=0.8,
               inset = inset, xpd=TRUE)
    }
        
}

repelwords <- function (x, y, words, cex = 1, rotate90 = FALSE,
                        xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
                        tstep=0.1, rstep = 0.1, maxiter=2000, ...) 
{
    ## From wordcloud::wordlayout
    tails <- "g|j|p|q|y"
    n <- length(words)
    sdx <- sd(x, na.rm = TRUE)
    sdy <- sd(y, na.rm = TRUE)
    ## rstep <- rstep * min(diff(range(x)),diff(range(y)))     ## IK
    if (sdx == 0) 
        sdx <- 1
    if (sdy == 0) 
        sdy <- 1
    if (length(cex) == 1) 
        cex <- rep(cex, n)
    if (length(rotate90) == 1) 
        rotate90 <- rep(rotate90, n)
    boxes <- list()
    ##theta <- runif(1, 0, 2*pi)
    theta <- 0
    for (i in 1:length(words)) {
        rotWord <- rotate90[i]
        r <- 0
        ## theta <- runif(1, 0, 2*pi)
        ## message("[repelwords] theta0 = ",theta)
        x1 <- xo <- x[i]
        y1 <- yo <- y[i]
        ##wid <- strwidth(words[i], cex = cex[i], ...)
        ##ht  <- strheight(words[i], cex = cex[i], ...)
        wid <- strwidth(words[i], cex = cex[i])
        ht  <- strheight(words[i], cex = cex[i])        
        if (grepl(tails, words[i])) 
            ht <- ht + ht * 0.2
        if (rotWord) {
            tmp <- ht
            ht <- wid
            wid <- tmp
        }
        isOverlapped <- TRUE
        iter=1
        ## maxiter=10000
        while(isOverlapped && iter<maxiter) {
            if(!wordcloud:::is_overlap(x1 - 0.5 * wid, y1 - 0.5 * ht, wid, 
                                       ht, boxes) &&
                x1 - 0.5 * wid > xlim[1] &&
                y1 - 0.5 * ht > ylim[1]  &&
                x1 + 0.5 * wid < xlim[2] &&
                y1 + 0.5 * ht < ylim[2])
            {
                boxes[[length(boxes) + 1]] <-
                    c(x1 - 0.5 * wid, y1 - 0.5 * ht, wid, ht)
                isOverlapped <- FALSE
            } else {
                theta <- theta + tstep
                r <- r + rstep * tstep / (2*pi)
                ##r <- r + rstep * theta / (2*pi)                
                x1 <- xo + sdx * r * cos(theta)
                y1 <- yo + sdy * r * sin(theta)
            }
            iter <- iter + 1
        }
        if(isOverlapped && iter==maxiter) {
            message("[repelwords] warning maximum iterations reached: iter = ",iter)
            boxes[[length(boxes) + 1]] <-
                c(x1 - 0.5 * wid, y1 - 0.5 * ht, wid, ht)
        } else {
            ## message("[repelwords] iter = ",iter)
        }
    }
    result <- do.call(rbind, boxes)

    message("[repelwords] nrow(results) = ",nrow(result))
    message("[repelwords] length(words) = ",length(words))
    
    colnames(result) <- c("x", "y", "width", "ht")
    rownames(result) <- words
    result
}

##=================================================================================
## PGX level plotting API
##=================================================================================


##level="geneset";contrasts=NULL
pgx.ActivationMatrix <- function(pgx, features=NULL, contrasts=NULL,
                                 n = 50, qsize=TRUE,
                                 cex=1, cex.row=1, cex.col=1, srt=90,
                                 flip=FALSE, level="geneset", 
                                 clust.x=TRUE, clust.y=TRUE, plotlib="base")
{
    

    if(level=="geneset") {
        out <- pgx.getMetaFoldChangeMatrix(pgx, what="meta", level="geneset")
    } else {
        out <- pgx.getMetaFoldChangeMatrix(pgx, what="meta", level="gene")
    }
    F <- out$fc
    Q <- out$qv
    if(!is.null(contrasts)) {
        ct <- intersect(contrasts,colnames(F))
        F <- F[,ct,drop=FALSE]
        Q <- Q[,ct,drop=FALSE]
    }
    if(!is.null(features)) {
        features <- intersect(features,rownames(F))
        F <- F[features,,drop=FALSE]
        Q <- Q[features,,drop=FALSE]
    }

    F1 <- Matrix::head(F[order(-apply(F,1,sd)),],n=n)
    F1 <- Matrix::head(F[order(-rowMeans(F**2)),],n=n)
    ##F1 <- F[grep("HALLMARK",rownames(F)),]
    dim(F1)

    ## cluster
    ii <- 1:nrow(F1)
    jj <- 1:ncol(F1)
    if(clust.y) ii <- fastcluster::hclust(dist(F1))$order
    if(clust.x) jj <- fastcluster::hclust(dist(t(F1)))$order
    F1 <- F1[ii,jj]
    Q1 <- Q[rownames(F1),colnames(F1)]

    cpal <- rev(RColorBrewer::brewer.pal(11,"RdYlBu"))
    cpal <- colorRampPalette(c("blue4","lightskyblue1","lightyellow","rosybrown1","red4"))(64)
    cpal <- colorspace::diverge_hcl(64, c=60, l=c(30,100), power=1)
    
    p <- NULL
    if(plotlib=="base") {
        
        F2 <- F1
        if(qsize) F2  <- F2 * Q1
        if(flip)
            F2 <- t(F2)
        colnames.f2 <- colnames(F2)
        colnames(F2) <- rep("",ncol(F2))
        corrplot::corrplot(
            F2, ## p.mat=Q1,
            is.corr=FALSE, col = cpal,
            tl.col = "grey0", tl.cex = 0.8*cex.row,
            tl.pos = "lt", ## tl.srt=45, tl.offset=1,
            cl.cex = 0.6,
            mar = c(10, 0, 1, 0)
        )
        ##mtext(1:ncol(F1), side=1, line=-5, at=1:ncol(F1),
        ##srt=30, cex=0.8)
        text(1:ncol(F2), 0, colnames.f2,
             srt=srt, xpd=TRUE, adj=1, cex=0.8*cex.col )

    }

    df <- NULL
    if(plotlib %in% c("ggplot","plotly")) {
        df <- reshape2::melt(F1)
        df$size <- df$value
        tt.size = "value"
        if(qsize) {
            df$size <- -log10(0.0001+Q[cbind(df$Var1, df$Var2)])
            tt.size <- "-log10q"
        }
        p <- ggplot2::ggplot(ggplot2::aes(x=Var2, y=Var1, color=value, size=size), data=df) +
            ggplot2::theme_minimal() + 
            ggplot2::geom_point() +
            ggplot2::scale_size(range = c(0.1,5*cex)) +       
            ggplot2::scale_color_gradientn( colors=cpal ) +
            ##scale_color_gradientn( colors=cpal, breaks=zz, labels=c(zz[1],zz[2])) +
            ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
            ggplot2::labs(size=tt.size, color="value") +
            ggplot2::scale_x_discrete(guide=guide_axis(angle=srt)) +
            ggplot2::theme(plot.margin = ggplot2::margin(5,0,0,10),
                  ## legend.title = ggplot2::element_blank(),
                  legend.text = ggplot2::element_text(size=9),
                  legend.key.size = grid::unit(10, "pt"),
                  legend.key.height = grid::unit(12, "pt")
                  ) 


        if(!qsize) {
            p <- p + ggplot2::guides(size=FALSE) 
        }
        if(flip) {
            p <- p + ggplot2::coord_flip() +
                ggplot2::scale_x_discrete(guide=guide_axis(angle=0)) +                
                ggplot2::scale_y_discrete(guide=guide_axis(angle=srt))
        }
    }
    if(plotlib == "plotly") {
        p <- p + ggplot2::scale_size(range = c(0.1,3*cex)) 
        p <- plotly::ggplotly(p) %>%
            plotly::layout(xaxis = list(tickangle=-srt, side='bottom') )
    }
    p
}

pgx.scatterPlot <- function(pgx, pheno=NULL, gene=NULL, 
                            contrast=NULL, level="gene", 
                            plotlib="base", pos=NULL, ...)
{
    ## Scatter t-SNE plot samples (or genes) colored on phenotype,
    ## gene expression, geneset expresssion or (correlation with)
    ## contrast.
    ##
    if(level=="geneset" && is.null(pos))  {
        stop("FATAL:: geneset scatter needs position vector")
    }

    if(level=="geneset")  {
        X <- pgx$gsetX
    } else {
        X <- pgx$X
    }
    
    plt <- NULL
    ## args <- list(...)
    if(!is.null(gene)) {
        if(is.null(pos) && level=="gene") pos <- pgx$tsne2d
        var <- X[gene,rownames(pos)]
        title <- gene
        plt <- pgx.scatterPlotXY(
            pos, var, plotlib=plotlib, ## title=title,
            xlab=colnames(pos)[1], ylab=colnames(pos)[2],
            type="numeric", ...)
    }
    if(!is.null(pheno)) {
        if(is.null(pos) && level=="gene") pos <- pgx$tsne2d            
        var <- pgx$samples[rownames(pos),pheno]
        title <- pheno
        plt <- pgx.scatterPlotXY(
            pos, var, plotlib=plotlib, ## title=title,
            xlab=colnames(pos)[1], ylab=colnames(pos)[2],
            ...)
    }
    if(!is.null(contrast)) {
        if(is.null(pos) && 'cluster.genes' %in% names(pgx)) {
            pos <- pgx$cluster.genes$pos[[1]][,1:2]
        }
        if(is.null(pos)) {
            stop('must supply positions for gene contrasts')
        }        
        if(level=="gene") {
            var <- pgx$gx.meta$meta[[contrast]]$meta.fx
            names(var) <- rownames(pgx$gx.meta$meta[[contrast]])
            var <- var[rownames(pos)]            
            tooltip <- pgx$genes[rownames(pos),"gene_title"]
        }
        if(level=="geneset") {
            var <- pgx$gset.meta$meta[[contrast]]$meta.fx
            names(var) <- rownames(pgx$gset.meta$meta[[contrast]])            
            var <- var[rownames(pos)]
            tooltip <- rownames(pos)
        }
        if(0) {
            pgx.scatterPlotXY.BASE(pos, var)
            pgx.scatterPlotXY.GGPLOT(pos, var)
            pgx.scatterPlotXY.PLOTLY(pos, var)
        }

        plt <- pgx.scatterPlotXY(
            pos, var, plotlib=plotlib, ## title=contrast,
            xlab=colnames(pos)[1], ylab=colnames(pos)[2],
            tooltip=tooltip,
            ...)
    }
    
    if(plotlib=="base") return(NULL)
    plt
}

plot.SPLOM <- function(F, F2=NULL, hilight=NULL, cex=0.5, cex.axis=1, cex.space=0.2)    
{

    if(is.null(F2)) F2 <- F    
    symm <- all(colnames(F) == colnames(F2))
    dim(F)
    gg <- intersect(rownames(F),rownames(F2))
    F <- F[gg,,drop=FALSE]
    F2 <- F2[gg,,drop=FALSE]
    
    x0 <- range(as.vector(apply(F,2,quantile, probs=c(0.001,0.999))))
    x1 <- range(as.vector(apply(F2,2,quantile, probs=c(0.001,0.999))))    
    x0 <- range(as.vector(F))
    x1 <- range(as.vector(F2))
    x0 <- x0 + c(-1,1)*diff(x0)*0.05
    x1 <- x1 + c(-1,1)*diff(x1)*0.05    
    
    nr=3;nc=4
    nr = ncol(F2)
    nc = ncol(F)
    par(mfrow=c(nr,nc), mar=c(1,1,1,1)*cex.space, oma=c(4,4,0,0),
        bty='o', xpd=FALSE)
    i=1
    for(j in ncol(F2):1) {        
        for(i in 1:ncol(F)) {
            plot( F[,i], F2[,j], pch=20, col='grey70',
                 cex=cex, cex.axis=0.95*cex.axis, cex.lab=cex.axis,
                 xaxt = ifelse( j==1, 's','n'),
                 yaxt = ifelse( i==1, 's','n'),
                 xlab = "", ylab = "",
                 ##xlab = ifelse( j==nr,colnames(F)[i],''),
                 ##ylab = ifelse( i==1,colnames(F2)[j],''),
                 xlim=x0, ylim=x1)
            abline(v=0, h=0, lty=3, lwd=0.6)
            if(j==1) mtext( colnames(F)[i], 1, line=2.7, cex=0.9*cex.axis)
            if(i==1) mtext( colnames(F2)[j], 2, line=2.7, cex=0.9*cex.axis)            
            ## title(colnames(F)[i], cex.main=1, font.main=1, line=0.3)
            if(!is.null(hilight) && length(hilight)>0) {
                hilight <- intersect(hilight, rownames(F))
                ii <- match(hilight, rownames(F))
                points( F[ii,i], F2[ii,j], pch=20, col='red3' )
                text( F[ii,i], F2[ii,j], labels=hilight, cex=0.85, pos=3, col='black')
            }
            
            ## write correlation value
            rho <- stats::cor(F[,i], F2[,j], use='pairwise')
            rr <- paste("r =",round(rho,digits=3))
            legend("topleft", legend=rr, bty='n', cex=1)
            
        }
    }
}

pgx.plotFoldchangeSPLOM <- function(pgx, contrasts=NULL) {

}

pgx.SankeyFromMatrixList.PLOTLY <- function(matlist, contrast=NULL)
{    
    ##------------------------------------------------------------------------
    ## Prepare matrices 
    ##------------------------------------------------------------------------    
    X <- list()
    for(i in 1:length(matlist)) {
        X[[i]] <- matlist[[i]] - rowMeans(matlist[[i]])
        X[[i]] <- X[[i]] / apply(X[[i]], 1, sd)
    }

    ## Counts cross-table between matrices
    M <- list()
    i=1
    for(i in 1:(length(X)-1)) {
        if(0) {
            k1 <- rownames(X[[i]])[max.col(t(X[[i]]))] ## only max???
            k2 <- rownames(X[[i+1]])[max.col(t(X[[i+1]]))]  ## only max???
            ##M[[i]] <- log(1+table(k1,k2))
            mm <- table(k1,k2)
            ii <- match(rownames(X[[i]]),rownames(mm))
            jj <- match(rownames(X[[i+1]]),colnames(mm))
            mm <- mm[ii,jj]
            rownames(mm) <- rownames(X[[i]])
            colnames(mm) <- rownames(X[[i+1]])
        } else {
            mm <- pmax(X[[i]],0) %*% t(pmax(X[[i+1]],0))
            mm <- mm**4
            mm <- mm / mean(mm)
        }
        M[[i]] <- mm
    }
    
    ## Correlation
    R <- list()
    for(i in 1:(length(X)-1)) {
        r1 <- stats::cor( t(X[[i]]), t(X[[i+1]]) )
        R[[i]] <- pmax(r1,0)
    }

    ## Edge value (i.e. capacity) : rho * contrast/FC
    cty.mode=1
    F = R
    if(!is.null(contrast)) {
        fc <- lapply(X, function(m) stats::cor(t(m),contrast)[,1])
        i=1
        for(i in 1:length(R)) {
            if(cty.mode==1) node.wt <- outer(pmax(fc[[i]],0), pmax(fc[[i+1]],0))
            if(cty.mode==3) node.wt <- abs(outer(fc[[i]], fc[[i+1]]))
            if(cty.mode==2) node.wt <- pmax(outer(fc[[i]], fc[[i+1]]),0)
            ww <- R[[i]] * node.wt
            ww <- ww / max(ww)  ## normalize??
            F[[i]] <- ww
        }
    }
    fill0 <- FALSE
    fill0 <- !is.null(contrast)
    pgx.SankeyFromMRF.PLOTLY(M=M, R=R, F=F, fill=fill0, labels=NULL)

}

pgx.SankeyFromMRF.PLOTLY <- function(M, R, F, fill=TRUE, labels=NULL)
{    

    rho2graph <- function(A, min.rho=0) {
        idx <- which(A > min.rho, arr.ind=TRUE)
        ee <- cbind( rownames(A)[idx[,1]], colnames(A)[idx[,2]] )
        gr <- igraph::graph_from_edgelist(ee, directed=TRUE)
        gr    
    }
        
    gr.list <- list()
    i=1
    for(i in 1:length(M)) {
        gr <- rho2graph(M[[i]], min.rho=0)
        ee <- igraph::get.edgelist(gr)
        pct <- M[[i]] / sum(M[[i]],na.rm=TRUE) * 100
        igraph::E(gr)$count <- pct[ee]
        igraph::E(gr)$rho   <- R[[i]][ee]
        igraph::E(gr)$weight <- F[[i]][ee]
        gr.list[[i]] <- gr
    }

    ## merge graphs
    gr <- gr.list[[1]]
    if(length(gr.list)>1) {
        i=2
        for(i in 2:length(gr.list)) {
            gr <- igraph::union(gr, gr.list[[i]])
            igraph::E(gr)$count <- rowSums(cbind(igraph::E(gr)$count_1,E(gr)$count_2),na.rm=TRUE)
            igraph::E(gr)$weight <- rowSums(cbind(igraph::E(gr)$weight_1,E(gr)$weight_2),na.rm=TRUE)
            igraph::E(gr)$rho <- rowSums(cbind(igraph::E(gr)$rho_1,E(gr)$rho_2),na.rm=TRUE)
        }
        gr <- igraph::delete_edge_attr(gr, "weight_1")
        gr <- igraph::delete_edge_attr(gr, "weight_2")
        gr <- igraph::delete_edge_attr(gr, "rho_1")
        gr <- igraph::delete_edge_attr(gr, "rho_2")
        gr <- igraph::delete_edge_attr(gr, "count_1")
        gr <- igraph::delete_edge_attr(gr, "count_2")
    }
    
    matnames <- c( list(rownames(M[[1]])), lapply(M,colnames))
    vlevel <- sapply(igraph::V(gr)$name, grep, matnames)
    table(vlevel)
        
    ## create Sankey plot
    nv <- length(igraph::V(gr))
    col1 <- rep(RColorBrewer::brewer.pal(12,"Set3"),100)[1:nv]
    ee <- igraph::get.edges(gr,E(gr))-1
    
    ee.label <- paste("rho=",round(igraph::E(gr)$rho,2),
                      "<br>count=",round(igraph::E(gr)$count,2),"%",
                      "<br>weight=",round(igraph::E(gr)$weight,2))
    
    nodes <- data.frame(label = igraph::V(gr)$name, color=col1)
    nodes$info  <- paste(igraph::V(gr)$name, unlist(labels)[igraph::V(gr)$name])
    nodes$label <- sub(".*[:]","",nodes$label)
    nodes$x <- (vlevel-1) / max(vlevel-1)

    if(fill) {
        wt <- igraph::E(gr)$weight
        ev2 <- 0.05 + 0.55 * (wt / max(wt))
        col2 <- paste("rgba(80,80,120,",ev2,")")
    } else {
        col2 <- paste("rgba(80,80,120,0.2)")
    }

    links <- data.frame(
        source = ee[,1],
        target = ee[,2],
        count = igraph::E(gr)$count,
        label = ee.label,
        weight = igraph::E(gr)$weight,
        ##color = "rgba(80,80,120,0.15)"
        color = col2
    )

    fig <- plotly::plot_ly(
        type = "sankey",
        domain = list(
            x =  c(0,1),
            y =  c(0,1)
        ),
        orientation = "h",
        valueformat = ".2f",
        valuesuffix = "",
        arrangement = "snap",

        node = list(
            label = nodes$label,
            ## label = nodes$info,
            ## customdata = nodes$info,
            ## hovertemplate = 'Info: %{customdata}',
            x     = nodes$x,
            y     = 0.01 * (1:length(nodes$x)),
            color = nodes$color,
            pad   = 15,
            thickness = 15,
            line = list(
                color = "black",
                width = 0.5
            )
        ),
        
        link = list(
            source = links$source,
            target = links$target,
            value =  links$count,
            label =  links$label,
            ##line = list( color = "black", width = 0.5)
            color =  links$color
        )
    )
    fig
    
}

pgx.SankeyFromPhenotypes.PLOTLY <- function(pgx, phenotypes, mat=NULL,
                                            fill=NULL, nmin=1, title="")
{    
    prefix1="from";prefix2="to"
    table2graph <- function(tab, prefix1, prefix2) {
        rownames(tab) <- paste0(prefix1,":",rownames(tab))
        colnames(tab) <- paste0(prefix2,":",colnames(tab))
        n <- nrow(tab) + ncol(tab)
        A <- matrix(0, nrow=n, ncol=n)
        rownames(A) <- colnames(A) <- c(rownames(tab),colnames(tab))    
        A[rownames(tab),colnames(tab)] <- tab
        gr <- igraph::graph_from_adjacency_matrix(A, mode="directed", diag=FALSE, weighted=TRUE)
        ## igraph::E(gr)$weight <- A[cbind(igraph::ends(gr,E(gr)))]
        gr    
    }

    if(is.null(mat)) {
        phenotypes <- intersect(phenotypes, colnames(pgx$samples))
        A <- pgx$samples[,phenotypes]
    }
    if(!is.null(mat)) {
        A <- mat
        phenotypes <- colnames(mat)
    }
        
    gr.list <- list()
    for(i in 1:(ncol(A)-1)) {
        p1 <- A[,phenotypes[i]]
        p2 <- A[,phenotypes[i+1]]
        P <- table(p1, p2)
        P[which(P<nmin,arr.ind=TRUE)] <- 0
        gr.list[[i]] <- table2graph(P,phenotypes[i],phenotypes[i+1])
    }

    gr <- gr.list[[1]]
    if(length(gr.list)>1) {
        i=2
        for(i in 2:length(gr.list)) {
            gr <- igraph::union(gr, gr.list[[i]])
            igraph::E(gr)$weight <- rowSums(cbind(igraph::E(gr)$weight_1,E(gr)$weight_2),na.rm=TRUE)
        }
        gr <- igraph::delete_edge_attr(gr, "weight_1")
        gr <- igraph::delete_edge_attr(gr, "weight_2")
    }
    
    ## create Sankey plot
    nv <- length(igraph::V(gr))
    col1 <- rep(RColorBrewer::brewer.pal(12,"Set3"),100)[1:nv]
    ee <- igraph::ends(gr, igraph::E(gr))
    ee <- apply(ee, 2, function(v) match(v, igraph::V(gr)$name)-1)
    
    nodes <- data.frame(label = igraph::V(gr)$name, color=col1)
    nodes$label <- sub(".*[:]","",nodes$label)
    links <- data.frame(
        source = ee[,1],
        target = ee[,2],
        value = igraph::E(gr)$weight,
        label = ""
    )

    fig <- plotly::plot_ly(
        type = "sankey",
        domain = list(
            x =  c(0,1),
            y =  c(0,1)
        ),
        orientation = "h",
        valueformat = ".0f",
        valuesuffix = "cells",
        
        node = list(
            label = nodes$label,
            color = nodes$color,
            pad = 15,
            thickness = 15,
            line = list(
                color = "black",
                width = 0.5
            )
        ),
        
        link = list(
            source = links$source,
            target = links$target,
            value =  links$value,
            label =  links$label,
            ##line = list( color = "black", width = 0.5)
            color = "rgba(80,80,120,0.15)"
        )
    )
    fig
}


## fill=NULL;title="";nmin=5
pgx.SankeyFromPhenotypes.GGPLOT <- function(pgx, phenotypes, mat=NULL, fill=NULL,
                                            sort=FALSE, nmin=1, title="")
{
    

    if(is.null(mat)) {
        phenotypes <- intersect(phenotypes,colnames(pgx$samples))
        phenotypes <- Matrix::head(phenotypes,5)
        pp <- intersect(c(phenotypes,fill), colnames(pgx$samples))   
        A <- pgx$samples[,pp] 
    }
    if(!is.null(mat)) {
        A <- mat
        phenotypes <- colnames(mat)
    }
    
    pt <- apply(A,1,paste,collapse=":")
    N <- table(pt)
    df <- do.call(rbind, strsplit(names(N),split=":"))
    colnames(df) <- phenotypes
    zap.tiny <- (N < nmin)
    if(sum(zap.tiny)>0) {
        n0 <- sum(N[zap.tiny])
        df <- df[!zap.tiny,]
        N  <- N[!zap.tiny]
        if(0) {
            df <- rbind(df, "other")
            N <- c(N,n0)
        }
    }
    df <- data.frame(df)
    df$Frequency <- N
    Matrix::head(df)
    
    if(sort) {
        relevelBig2small <- function(idx) factor(idx, levels=names(sort(table(idx))))
        df[,phenotypes] <- data.frame(lapply(df[,phenotypes], relevelBig2small))
    }
    
    v <- sapply(phenotypes, function(s) {a=sym(s);enquo(a)})    
    if(length(phenotypes)==2) {
        p <- ggplot2::ggplot(df, ggplot2::aes(y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]]))
    }
    if(length(phenotypes)==3) {
        p <- ggplot2::ggplot(df, ggplot2::aes(y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]],
                            axis3 = !!v[[3]] ))        
    }
    if(length(phenotypes)==4) {
        p <- ggplot2::ggplot(df, ggplot2::aes(y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]],
                            axis3 = !!v[[3]], axis4 = !!v[[4]] ))        
    }
    if(length(phenotypes)==5) {
        p <- ggplot2::ggplot(df, ggplot2::aes(y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]],
                            axis3 = !!v[[3]], axis4 = !!v[[4]], axis5 = !!v[[5]] ))        
    }
    if(length(phenotypes)==6) {
        p <- ggplot2::ggplot(df, ggplot2::aes(y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]],
                            axis3 = !!v[[3]], axis4 = !!v[[4]], axis5 = !!v[[5]],
                            axis6 = !!v[[6]] ))        
    }
    if(!is.null(fill)) {
        a <- ggplot2::sym(fill)
        vf <- ggplot2::enquo(a)
        p <- p + geom_alluvium(ggplot2::aes(fill=!!vf, color=!!vf), width=1/12, alpha=0.4) 
    } else {
        p <- p + geom_alluvium(width = 1/12, alpha=0.4) 
    }
    p <- p + 
        geom_stratum(width = 1/3, fill = "grey", color = "grey50") +
        ggplot2::geom_label(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum))) +
        ggplot2::scale_x_discrete(limits = phenotypes, expand = c(.05, .05)) +
        ggplot2::scale_fill_brewer(type = "qual", palette = "Set1", direction=-1) +
        ggplot2::scale_color_brewer(type = "qual", palette = "Set1", direction=-1) +
        ggplot2::ggtitle(title) + ggplot2::theme_minimal() +
        ggplot2::theme(
            ##plot.title = ggplot2::element_text(hjust=0, vjust=3),
            axis.text.x = ggplot2::element_text(size=13, vjust=+7)
        )
    p
}

plot_grid.sharedAxisLabels <- function(plotList, nrow) {
    np <- length(plotList)
    np
    ncol <- ceiling(np / nrow)
    ncol
    ann.y <- which((0:(np-1) %% ncol)!=0)
    ann.x <- which((0:(np-1) %/% ncol)!=(nrow-1))
    for(i in ann.y) {
        plotList[[i]] <- plotList[[i]] + ggplot2::ylab("")
    }
    for(i in ann.x) {
        plotList[[i]] <- plotList[[i]] + ggplot2::xlab("")
    }    
    cowplot::plot_grid(plotlist=plotList, nrow=nrow, labels=NA)
}

pgx.plotContrast <- function(pgx, contrast=NULL, type='scatter',
                             set.par=TRUE, par.sq=FALSE, ...)
{
    if(0) {
        contrast= colnames(pgx$model.parameters$exp.matrix)
        contrast
    }
    if(is.null(contrast)) {
        contrast <- colnames(pgx$model.parameters$exp.matrix)
    }
    
    if(set.par) {
        nc <- ceiling(sqrt(length(contrast)))
        nr <- ceiling(length(contrast)/nc)
        if(par.sq) nr <- nc
        par(mfrow=c(nr,nc))
    }
    
    plist <- list()
    for(i in 1:length(contrast)) {
        ct <- contrast[i]

        if(type == 'volcano') {
            p <- pgx.Volcano(pgx, contrast=ct, ...)        
        } else if(type == 'MA') {
            p <- pgx.plotMA(pgx, contrast=ct, ...)        
        } else if(type == 'UMAP') {
            p <- pgx.plotGeneUMAP(pgx, contrast=ct, set.par=FALSE, ...)        
        } else {
            ## scatter
            p <- pgx.contrastScatter(pgx, contrast=ct, ...)        
        }
        plist[[i]] <- p
    }

    if(length(plist)==1) plist <- plist[[1]]
    plist
}

pgx.Volcano <- function(pgx, contrast, level="gene", methods='meta',
                        psig=0.05, fc=1, cex=1, cex.lab=1, ntop=20,
                        p.min=NULL, fc.max=NULL, hilight=NULL, ## hilight2=NULL, 
                        cpal=c("grey60","red3"), title=NULL,
                        plotlib="base")
{
    if(is.integer(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]
    res=NULL
    if(level=="gene") {
        res <- pgx.getMetaMatrix(pgx, methods, level="gene")
    } else if(level=="geneset") {
        res <- pgx.getMetaMatrix(pgx, methods, level="geneset")        
    } else {
        stop("FATAL:: invalid level=",level)
    }
    f <- res$fc[,contrast]
    q <- res$qv[,contrast]    
    
    sig <- (q <= psig & abs(f) >= fc)
    xy <- cbind(fc=f, y=-log10(q))
    ##rownames(xy) <- rownames(pgx$gx.meta$meta[[contrast]])
    ##names(sig) <- rownames(xy)
    ##cpal <- c("grey70","red3")
    if(is.null(hilight)) {
        wt <- rowSums(scale(xy,center=FALSE)**2)
        hilight <- rownames(xy)[order(-wt)]
        hilight <- intersect(hilight, names(sig[sig==TRUE]))
    }
    hilight <- Matrix::head(hilight,ntop) ## label

    xlim = ylim = NULL
    if(!is.null(fc.max)) {
        xlim = c(-1.1,1.1)*fc.max
    }
    if(!is.null(p.min)) {
        ylim = c(0,-log10(p.min))
    }
    
    if(is.null(title)) title = contrast
    p <- pgx.scatterPlotXY(
        xy, var=sig, type="factor", title=title,
        xlab = "differential expression (log2FC)",
        ylab = "significance (-log10q)",
        hilight = hilight, ## hilight2 = hilight2,
        cex = cex, cex.lab = cex.lab, cex.title = 1.0,
        xlim = xlim, ylim = ylim,
        legend = FALSE, col=cpal, opacity=1,
        plotlib = plotlib)    
    ##ggplotly(p)
    p
}

pgx.plotMA <- function(pgx, contrast, level="gene", psig=0.05, fc=1,
                       cex=1, cex.lab=0.8, hilight=NULL, ntop=20,
                       plotlib="base")
{
    if(is.integer(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]    

    if(level=="gene") {
        f <- pgx$gx.meta$meta[[contrast]]$meta.fx
        q <- pgx$gx.meta$meta[[contrast]]$meta.q
        gg <- rownames(pgx$gx.meta$meta[[contrast]])
        names(f)=names(q)=gg
        m  <- rowMeans(pgx$X[gg,])
    } else if(level=="geneset") {
        f <- pgx$gset.meta$meta[[contrast]]$meta.fx
        q <- pgx$gset.meta$meta[[contrast]]$meta.q
        gg <- rownames(pgx$gset.meta$meta[[contrast]])
        names(f)=names(q)=gg        
        m  <- rowMeans(pgx$gsetX[gg,])
    } else {
        stop("FATAL:: invalid level=",level)
    }

    sig <- (q <= psig & abs(f) >= fc)
    table(sig)
    xy <- cbind(x=m, y=f)    
    ##names(sig) <- rownames(xy) <- gg
    cpal <- c("grey60","red3")
    if(is.null(hilight)) {
        wt <- rowSums(scale(cbind(xy,m),center=FALSE)**2)
        hilight <- rownames(xy)[order(-wt)]
        hilight <- intersect(hilight, names(sig[sig==TRUE]))
    }
    hilight <- Matrix::head(hilight,ntop)
    
    p <- pgx.scatterPlotXY(
        xy, var=sig, type="factor", title=contrast,
        xlab = "average expression  (log2)",
        ylab = "differential expression  (log2FC)",
        hilight = hilight, cex = 0.9*cex,
        cex.lab = cex.lab, cex.title = 1.0,
        legend = FALSE, col=c("grey70","red3"), opacity=1,
        plotlib = plotlib)    
    ##ggplotly(p)
    p
}

pgx.contrastScatter <- function(pgx, contrast, hilight=NULL,
                                cex=1, cex.lab=0.8, 
                                psig=0.05, fc=1, level="gene",
                                ntop=20, dir=0, plotlib="base")
{
    if(0) {
        contrast= colnames(pgx$model.parameters$exp.matrix)[1]
        contrast
    }
    if(is.numeric(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]
    exp.matrix <- pgx$model.parameters$exp.matrix
    ct <- exp.matrix[,contrast]
    ii <- which(ct < 0)
    jj <- which(ct > 0)
    if(level=="gene") {
        x0 <- rowMeans(pgx$X[,ii,drop=FALSE])
        x1 <- rowMeans(pgx$X[,jj,drop=FALSE])
        xy <- cbind(x0,x1)
        gg <- rownames(pgx$gx.meta$meta[[contrast]])
        fx <- pgx$gx.meta$meta[[contrast]]$meta.fx
        q <- pgx$gx.meta$meta[[contrast]]$meta.q
        names(fx)=names(q)=gg
    } else if(level=="gene") {
        x0 <- rowMeans(pgx$gsetX[,ii,drop=FALSE])
        x1 <- rowMeans(pgx$gsetX[,jj,drop=FALSE])
        xy <- cbind(x0,x1)
        gg <- rownames(pgx$gset.meta$meta[[contrast]])
        fx <- pgx$gset.meta$meta[[contrast]]$meta.fx
        q <- pgx$gset.meta$meta[[contrast]]$meta.q
        names(fx)=names(q)=gg
    } else {
        stop("FATAL:: invalid level=",level)
    }

    colnames(xy) <- c(""," ")
    grp1 <- gsub(".*[:]|_vs_.*","",contrast)
    grp0 <- gsub(".*_vs_","",contrast)
    xlab <- paste("expression",grp0,"  (logCPM)")
    ylab <- paste("expression",grp1,"  (logCPM)")
    
    if(is.null(hilight)) {
        top.gg <- c(Matrix::head(names(sort(fx)),ntop/2),
                    Matrix::tail(names(sort(fx)),ntop/2))
        if(dir>0) {
            top.gg <- Matrix::tail(names(sort(fx)),ntop)
        }
        if(dir<0) {
            top.gg <- Matrix::head(names(sort(fx)),ntop)
        }
        hilight <- top.gg
    }
    hilight <- Matrix::head(hilight,ntop)    

    sig <- 1*(q < psig & abs(fx) > fc)
    table(sig)
    names(sig) <- gg
    ## hilight <- hilight[which(sig[hilight]==1)]

    tt <- contrast
    ## tt <- sub(".*:","",tt)
    pgx.scatterPlotXY(
        xy, var=sig, type="factor", title=tt,
        xlab = xlab, ylab=ylab, 
        hilight = hilight, cex = 0.9*cex,
        cex.lab = cex.lab, cex.title = 1.0,
        legend = FALSE, col=c("grey70","red3"), opacity=1,
        plotlib = plotlib)    
    ##plotlib="base"
}

pgx.plotGeneUMAP <- function(pgx, contrast=NULL, value=NULL,
                             pos=NULL, ntop=20, cex=1, cex.lab=0.8,
                             hilight=NULL, title=NULL, zfix=FALSE,
                             set.par=TRUE, par.sq = FALSE,
                             level="gene", plotlib="ggplot")
{
    if(!is.null(contrast)) {
        if(is.numeric(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]
        res=NULL
        if(level=="gene") {
            res <- pgx.getMetaMatrix(pgx, level="gene")
        } else if(level=="geneset") {
            res <- pgx.getMetaMatrix(pgx, level="geneset")        
        } else {
            stop("FATAL:: invalid level=",level)
        }
        F <- res$fc[,contrast,drop=FALSE]
        ##q <- res$qv[,contrast]    
        ##sig <- (q <= psig & abs(f) >= fc)
    }
    dim(F)
    
    if(!is.null(value)) {
        F <- cbind(value)
    }
        
    if(!is.null(pos)) {
        xy <- pos
    } else if("cluster.genes" %in% names(pgx)) {
        xy <- pgx$cluster.genes$pos[['umap2d']]
    } else {
        message("[WARNING] no cluster.genes in object. please supply positions")
        return(NULL)
    }

    F <- F[match(rownames(xy),rownames(F)),,drop=FALSE]
    ## f = scale(rank(abs(f))**2 * sign(f),center=FALSE)[,1]
    rownames(F) <- rownames(xy)
    
    if(set.par) {
        nc <- ceiling(sqrt(ncol(F)))
        nr <- ceiling(ncol(F)/nc)
        if(par.sq) nr <- nc
        par(mfrow=c(nr,nc))
    }
    
    ## z-scale
    zlim = NULL
    if(zfix) {
        zlim <- quantile(F,probs=c(0.01,0.99),na.rm=TRUE)
        zlim <- quantile(F,probs=c(0.002,0.998),na.rm=TRUE)        
        zlim
    }
    
    plist <- list()
    i=1
    for(i in 1:ncol(F)) {

        title1 = contrast[i]
        if(is.null(title)) title1 = paste(title, contrast[i])

        f1 <- F[,i]

        hilight1 <- hilight
        if(is.null(hilight)) {
            hilight1 <- names(sort(-abs(f1)))
            ##hilight <- intersect(hilight, names(sig[sig==TRUE]))
        }
        hilight1 <- Matrix::head(hilight1,ntop) ## label
        opacity = ifelse(length(hilight1)>0,0.66,1)
        
        p1 <- pgx.scatterPlotXY(
            xy, var=f1, type="numeric", 
            xlab = "UMAP-x  (genes)",
            ylab = "UMAP-y  (genes)",
            hilight = hilight1,
            zlim = zlim, zsym = TRUE, softmax=1,
            ## hilight2 = hilight2,
            cex = cex, cex.lab = cex.lab,
            title = title1, cex.title = 1.0,
            ## xlim = xlim, ylim = ylim,
            legend = TRUE,
            ## col=cpal,
            opacity = 0.5,
            plotlib = plotlib)

        plist[[i]] <- p1
    }
    if(plotlib == "base") return()
    ##ggplotly(p)
    if(length(plist)==1) plist <- plist[[1]]    
    return(plist)
}

pgx.plotExpression <- function(pgx, probe, comp, logscale=TRUE,
                               level="gene", grouped=FALSE, srt=NULL, cex=1,
                               collapse.others=TRUE, showothers=TRUE,
                               max.points = 200, group.names=NULL, 
                               main=NULL, xlab=NULL, ylab=NULL, names=TRUE)
{
    if(0) {
        logscale=TRUE;level="gene";grouped=TRUE;srt=90;collapse.others=1;
        max.points=-1;main=NULL;xlab=NULL;ylab=NULL;names=TRUE;group.names=NULL
        probe = rownames(pgx$X)[1];comp=1
    }

    if(is.null(probe)) return(NULL)
    if(is.na(probe)) return(NULL)
    
    if(level=="gene" && !probe %in% rownames(pgx$X)) {
        frame() ## emtpy image
        return(NULL)
    }
    if(level=="geneset" && !probe %in% rownames(pgx$gsetX)) {
        frame() ## emtpy image
        return(NULL)
    }
    
    ## ------------- determine groups
    expmat  <- pgx$model.parameters$exp.matrix
    cntrmat <- pgx$model.parameters$contr.matrix
    expmat  <- expmat[rownames(pgx$samples),,drop=FALSE]
    
    if(class(comp)=="numeric") comp <- colnames(expmat)[comp]
    if(!is.null(group.names) && length(group.names)!=2) stop("group.names must be length=2")
    if(is.null(main)) main <- probe
    comp    

    if(!is.null(group.names) ) {
        message("[pgx.plotExpression] using group names from argument")
    }
    
    ## if a named contrast table is available it is safer
    if(is.null(group.names) && "contrasts" %in% names(pgx)  ) {
        message("[pgx.plotExpression] parsing group names from pgx$contrast labels")
        contr.labels <- pgx$contrasts[,comp]
        contr.idx    <- expmat[,comp]
        group1 <- names(which.max(table(contr.labels[contr.idx > 0])))
        group0 <- names(which.max(table(contr.labels[contr.idx < 0])))
        group.names <- c(group0, group1)
    }
    group.names
    
    ## Otherwise we guess from the contrast title but this is dangerous
    if(is.null(group.names) && grepl("_vs_|_VS_",comp) ) {
        message("[pgx.plotExpression] parsing group names contrast name (warning!)")        
        comp1 <- sub(".*:","",comp)  ## remove prefix
        group.names = strsplit(comp1,split="_vs_|_VS_")[[1]]
        group.names = rev(group.names)  ## first is main group
        ## !!!!!!!! NEED RETHINK !!!!!!!! which one is reference
        ## class???  determine if notation is A_vs_B or B_vs_A (some
        ## people do different) and in sample-based contrasts we don't
        ## know group names.
        if(!is.POSvsNEG(pgx)) {
            ## A_vs_B or B_vs_A notation !!!
            message("[pgx.plotExpression] WARNING! A_vs_B reversed?")            
            group.names <- rev(group.names) ## reversed!!
        } 
        group.names <- gsub("@.*","",group.names)  ## strip postfix
    }
    group.names

    ## label rotation
    if(is.null(srt)) {
        srt = 0
        if(max(nchar(group.names)) >= 7) srt = 45
    }
    
    ## create groups
    ct <- expmat[,comp]
    names(ct) <- rownames(expmat)
    samples <- rownames(expmat)[which(ct!=0)]
    grp0.name=grp1.name=NULL
    if(!is.null(group.names)) {
        grp0.name <- group.names[1]
        grp1.name <- group.names[2]
    } else {
        grp1.name <- comp
        grp0.name <- "REF"
    }

    xgroup <- c("other",grp0.name,grp1.name)[1 + 1*(ct<0) + 2*(ct>0)]       
    table(xgroup)
    
    ## currently cast to character... :(
    names(xgroup) <- rownames(pgx$samples)
    jj <- which(!(xgroup %in% xgroup[samples]))
    jj
    if(length(jj)>0 && collapse.others) {
        xgroup <- as.character(xgroup)
        xgroup[jj] <- "other"
    }
    names(xgroup) <- rownames(expmat)
    table(xgroup)
    class(xgroup)

    if(class(xgroup)=="character") {
        xgroup <- as.character(xgroup)
        ##levels0 <- xgroup[!duplicated(xgroup)]  ## original order??
        levels0 <- group.names
        if("other" %in% xgroup) levels0 <- c(levels0,"other")
        xgroup <- factor(xgroup, levels=levels0)
    }
    table(xgroup)    

    ## ------------- set color of samples    
    ngrp <- length(group.names)
    grp.klr = rep(RColorBrewer::brewer.pal(12,"Paired"),99)[1:ngrp]
    names(grp.klr) <- group.names
    ## if(is.null(comp) && grouped) grp.klr <- rep("grey60",length(grp.klr))
    if(any(grepl("other",xgroup))) {
        grp.klr <- c("other"="grey85", grp.klr)
    }
    grp.klr

    ## -------------- get expression value
    if(level=="geneset") {
        gx <- pgx$gsetX[probe,rownames(pgx$samples)]
    } else {
        gx <- pgx$X[probe,rownames(pgx$samples)]
    }
    if(!logscale) {
        gx <- 2**(gx)
    }
    
    ## -------------- remove others
    if(showothers==FALSE && any(grepl("other",xgroup)) ) {
        jj <- grep("other",xgroup,invert=TRUE)
        xgroup <- factor( xgroup[jj], levels=setdiff(levels(xgroup),"other"))
        gx <- gx[jj]
    }
    
    ## -------------- plot grouped or ungrouped
    if(is.null(main)) main <- probe
    ##if(ncol(X) <= 20) {
    if(!grouped) {
        nx = length(gx)
        if(is.null(ylab)) {
            ylab = "expression (log2CPM)"
            if(!logscale) ylab = "expression (CPM)"
        }
        klr = grp.klr[as.character(xgroup)]
        klr[is.na(klr)] <- "grey90"
        gx.min = 0
        if(min(gx)<0) gx.min <- min(gx)
        ylim <- c(gx.min,1.3*max(gx))
        bx = barplot( gx[], col=klr[], ylim=ylim,
                     ## offset = 0, ylim=c(gx.min,max(gx)),
                     las=3, ylab=ylab, names.arg=NA, border = NA)
        if(nx<20 && names==TRUE) {
            y0 <- min(ylim) - 0.03*diff(ylim)
            text( bx[,1], y0, names(gx)[], adj=1, srt=srt,
                 xpd=TRUE, cex=ifelse(nx>10,0.7,0.9)*cex )
        }
        title(main, cex.main=1.0, line=0)        
    } else {
        if(is.null(ylab)) {
            ylab = "expression (log2CPM)"
            if(!logscale) ylab = "expression (CPM)"
        }
        bee.cex = c(0.3,0.1,0.05)[cut(length(gx),c(0,100,500,99999))]
        ##grp.klr1 <- grp.klr[levels(xgroup)]
        xlevels <- levels(xgroup)
        grp.klr1 <- grp.klr[as.character(xlevels)]
        grp.klr1[is.na(grp.klr1)] <- "grey90"
        names(grp.klr1) <- as.character(xlevels)        
        ##col=grp.klr1;las=3;names.cex=cex;        
        gx.b3plot( gx, xgroup, ##ylim=c(0,max(gx)*1.2),
                  first = xlevels[1],
                  col = grp.klr1, ylab=ylab, bee.cex=bee.cex,
                  max.points=max.points, xlab=xlab, names=names,
                  ## sig.stars=TRUE, max.stars=5,
                  las=3, names.cex=cex, srt=srt)
        title(main, cex.main=1.0, line=0)
    }

}

pgx.plotOmicsNetwork <- function(pgx, gene=NULL, reduced=NULL, levels=c("gene","geneset"),
                                 contrast=NULL, layout=NULL, colorcluster=FALSE,
                                 hilight=NULL )
{
    
    
     ## for col2hex
    ## return(NULL)

    has.graph <- all(c("omicsnet","omicsnet.reduced") %in% names(pgx))
    has.graph
    if(!has.graph) return(NULL)

    ##gr <- pgx$genes_tsne_graph$graph
    gr <- pgx$omicsnet
    if(is.null(gr)) return(NULL)
    if(is.null(reduced) && is.null(gene)) gr <- pgx$omicsnet.reduced

    if(!is.null(gene)) {
        gene0 = paste0("{gene}",gene)
        gene0
        k <- which(igraph::V(gr)$name %in% c(gene, gene0))
        k
        nb <- names(igraph::neighbors(gr, igraph::V(gr)[k]))
        vv <- unique(c(gene0, nb))
        gr <- igraph::induced_subgraph(gr, vv)
    }
    gr
    gr <- igraph::induced_subgraph(gr, which(igraph::V(gr)$level %in% levels))
    if(is.null(gr)) return(NULL)

    ## ------------- get fold-change for node color and size ------------------
    fc0 <- gr$foldchange[igraph::V(gr)$name,,drop=FALSE]
    if(is.null(contrast)) {
        fc <- rowMeans(fc0**2,na.rm=TRUE)**0.5
    } else {
        if(!(contrast %in% colnames(fc0))) stop("unknown contrast")
        fc <- fc0[,contrast]
    }
    fc[is.na(fc)] <- 0
    fc <- fc / max(abs(fc),na.rm=TRUE)
    fc.cex <- (0.01+abs(fc))**0.66

    ## defaults graph parameters
    vlabel <- igraph::V(gr)$name
    if("label" %in% igraph::vertex_attr_names(gr)) vlabel <- igraph::V(gr)$label
    vlabel0 = vlabel
    vklr = c("blue3","grey70","red3")[2+sign(fc)]
    ##do.colorcluster = FALSE
    if(colorcluster) {
        vklr <- rep(rainbow(16),99)[igraph::V(gr)$cluster]
    }
    lab.cex = 1

    ee <- igraph::get.edgelist(gr)
    ee <- igraph::get.edges(gr, igraph::E(gr))
    Matrix::head(ee)
    ##ew = rep(3,nrow(ee))
    ew = 1 + 5 * sqrt(fc.cex[ee[,1]]*fc.cex[ee[,2]])
    ew = 1 + 5 * abs(igraph::E(gr)$weight)
    vsel <- rep(0,length(fc))
    esel <- rep(0,nrow(ee))

    ## ------------------ highlight selection with labels
    if(!is.null(hilight) && length(hilight)) {
        ##sel <- pgx$collections[[10]]
        sel <- hilight
        ##mm  <- pgx$genes_tsne_graph$members
        mm <- igraph::V(gr)$name
        if(!is.null(gr$members)) {
            mm <- gr$members[igraph::V(gr)$name]
        }
        mm <- sapply(mm, function(s) sub(".*\\}","",s))
        vlabel <- sapply( mm, function(x) intersect(x,sel) )
        vlabel <- sapply( vlabel, paste, collapse="\n")
        Matrix::head(vlabel,10)
        table(vlabel!="")

        sel <- which(vlabel!="")
        length(sel)
        if(length(sel)>0) {
            vsel[sel] <- 1
            lab.cex[sel] <- 1 + 18*(fc.cex[sel]/max(fc.cex[sel],na.rm=TRUE))
        }

        jj <- which( vsel[ee[,1]]==1 | vsel[ee[,2]]==1 )
        esel[jj] <- 1
        ew[jj] <- 2.4 * ew[jj]

        nnb <- unique(unlist(sapply(sel, neighbors, graph=gr)))
        is.nb <- (1:length(sel) %in% nnb)
        vlabel[which(vsel==0)] <- NA
        ##vklr[which(vsel==0 & !is.nb)] <- "grey70"
        vklr[which(vsel==0)] <- "grey60"
        lab.cex[which(vsel==0)] <- 1
    }


    vklr = substring(gplots::col2hex(vklr),1,7)
    names(vklr) <- igraph::V(gr)$name
    igraph::V(gr)$label <- vlabel ## filtered labels
    igraph::V(gr)$title <- gsub("\n","<br>",vlabel0)  ## tooltip has complete names
    ##V(gr)$label.cex <- 0.12 * lab.cex
    igraph::V(gr)$size <- 40*fc.cex
    igraph::V(gr)$color <- paste0(vklr,ifelse(vsel==1,"99","55"))
    ##if(input$gr_labels==FALSE) igraph::V(gr)$label <- NA
    ##E(gr)$color <- paste0(vklr[ee[,1]],ifelse(esel==1,"33","22"))
    igraph::E(gr)$color <- paste0(vklr[ee[,1]],ifelse(esel==1,"99","55"))
    igraph::E(gr)$width <- 1 * (2 + 5 * (ew / max(ew)))

    gr

    ##pos <- as.matrix(gr$layout)
    if(!is.null(layout)) {
        layout.fun <- match.fun(layout)
        tmp.gr = gr
        igraph::E(tmp.gr)$weight = abs(igraph::E(tmp.gr)$weight)**0.2
        ##E(tmp.gr)$weight = 1
        pos <- layout.fun(tmp.gr)
        remove(tmp.gr)
        rownames(pos) <- igraph::V(gr)$name
    }

    visdata <- visNetwork::toVisNetworkData(gr, idToLabel=FALSE)
    pos <- pos[igraph::V(gr)$name,]
    pos[,2] <- -pos[,2]
    dim(pos)


    ## ------------------ plot using visNetwork (zoomable) -----------------
    graph <- visNetwork::visNetwork(nodes = visdata$nodes, edges = visdata$edges,
                        height="1200px", width="1600px") %>%
        visNetwork::visNodes(font=list(size=14))  %>%
        visNetwork::visEdges(hidden=FALSE, width=2, color=list(opacity=0.9))  %>%
        visNetwork::visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
        ##visHierarchicalLayout(direction = "LR") %>%
        ##visInteraction(hideEdgesOnDrag = TRUE) %>%
        visNetwork::visIgraphLayout(layout="layout.norm", layoutMatrix=pos)
    graph

}

## gene1="TOPORS";gene2="KCNN4";col="black";cex=1;k=11;samples=NULL;cex.names=1;lab.unit=NULL
pgx.cytoPlot <- function(pgx, gene1, gene2, cex=1, col="grey60",
                         lab.unit=NULL, cex.names=1, samples=NULL, k=11)
{
    
    

    ## some pretty colors
    ##k <- 9
    my.cols <- rev(RColorBrewer::brewer.pal(k, "RdYlBu"))
    ##gene1 = "CD8A"
    ##gene2 = "CD4"
    if(is.null(samples))
        samples <- colnames(pgx$X)
    samples <- intersect(samples, colnames(pgx$X))
    x1 <- pgx$X[gene1,samples]
    x2 <- pgx$X[gene2,samples]
    x1 <- x1 + 1e-3*rnorm(length(x1))
    x2 <- x2 + 1e-3*rnorm(length(x2))
    names(x1) <- samples
    names(x2) <- samples
    m1 <- mean(x1)
    m2 <- mean(x2)

    ## select samples in different quadrants
    j1 <- samples[which(x1 < m1 & x2 > m2)]
    j2 <- samples[which(x1 > m1 & x2 < m2)]
    j3 <- samples[which(x1 > m1 & x2 > m2)]
    j4 <- samples[which(x1 < m1 & x2 < m2)]
    if(0) {
        j1 <- c(j1, sample(samples,5))
        j2 <- c(j2, sample(samples,5))
        j3 <- c(j3, sample(samples,5))
        j4 <- c(j4, sample(samples,5))
    }
    z1=z2=z3=z4=NULL
    if(length(j1)>1) z1 <- MASS::kde2d( x1[j1], x2[j1], n=50)
    if(length(j2)>1) z2 <- MASS::kde2d( x1[j2], x2[j2], n=50)
    if(length(j3)>1) z3 <- MASS::kde2d( x1[j3], x2[j3], n=50)
    if(length(j4)>1) z4 <- MASS::kde2d( x1[j4], x2[j4], n=50)
    ##z0 <- MASS::kde2d( x1[], x2[], n=50)

    ##par(mfrow=c(1,1))
    xlab1 <- paste(gene1, lab.unit, collapse="  ")
    ylab1 <- paste(gene2, lab.unit, collapse="  ")
    plot(x1, x2, xlab=xlab1, ylab=ylab1, col=col, pch=19, cex=cex)
    abline(h=mean(x1), v=mean(x2), lwd=1, lty=2)
    ##z0$z <- log2(z0$z)
    ##contour(z0, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)

    if(length(j1)>10) contour(z1, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
    if(length(j2)>10) contour(z2, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
    if(length(j3)>10) contour(z3, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
    if(length(j4)>10) contour(z4, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)

    N = length(x1)
    d1 = 0.02*max(x1)
    d2 = 0.04*max(x2)
    legend( "topright", paste(round(100*length(j3)/N,2),"%"), cex=1.2, col="gray50",
           bty="n", xpd=TRUE)
    legend( "topleft", paste(round(100*length(j1)/N,2),"%"), cex=1.2, col="gray50",
           bty="n", inset=c(-0.05,0), xpd=TRUE )
    legend( "bottomright", paste(round(100*length(j2)/N,2),"%"), cex=1.2, col="gray50",
           bty="n", xpd=TRUE)
    legend( "bottomleft", paste(round(100*length(j4)/N,2),"%"), cex=1.2, col="gray50",
           bty="n", inset=c(-0.05,0), xpd=TRUE )

    if(!is.null(pgx$deconv)) {
        inferred.celltype <- pgx$deconv[[1]][["meta"]]
        dim(inferred.celltype)
        ##cex.names=2
        lab1 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j1,,drop=FALSE]))),3)
        pos1 <- apply(cbind(x1, x2)[j1,,drop=FALSE],2,median)
        text( pos1[1], pos1[2], paste(lab1,collapse="\n"),cex=0.9*cex.names, xpd=TRUE)

        lab2 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j2,,drop=FALSE]))),3)
        pos2 <- apply(cbind(x1, x2)[j2,,drop=FALSE],2,median)
        text( pos2[1], pos2[2], paste(lab2,collapse="\n"),cex=0.9*cex.names, xpd=TRUE)

        lab3 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j3,,drop=FALSE]))),3)
        pos3 <- apply(cbind(x1, x2)[j3,,drop=FALSE],2,median)
        text( pos3[1], pos3[2], paste(lab3,collapse="\n"),cex=0.9*cex.names, xpd=TRUE)

    }

}

pgx.plotPhenotypeMatrix <- function(annot)
{
    ## ------- set colors
    colors0 = rep("Set2",ncol(annot))    
    names(colors0) = colnames(annot)
    
    grid_params <- iheatmapr::setup_colorbar_grid(
        nrows = 3, 
        y_length = 0.32, 
        y_spacing = 0.39, 
        x_spacing = 0.17,
        x_start = 1.1, 
        y_start = 0.96)

    ## maximize plot area
    mar <- list(l = 150, r = 0, b = 5, t = 0, pad = 3)

    
    plt <- NULL
    empty.X <- matrix(NA,nrow=1,ncol=nrow(annot))
    colnames(empty.X) <- rownames(annot)
    plt <- iheatmapr::main_heatmap(
        empty.X,
        name = "phenotype data",
        show_colorbar = FALSE,
        colorbar_grid = grid_params,
        layout = list(margin = mar)
    )

    ## ------- annot need to be factor
    annotF <- data.frame(as.list(annot),stringsAsFactors=TRUE)
    rownames(annotF) = rownames(annot)
    cvar <- which(sapply(annotF,is.factor))
    cvar
    annotX <- annotF
    annotX[,cvar] <- data.frame(lapply(annotF[,cvar],as.integer))
    annotX[,] <- data.frame(lapply(annotX[,],as.numeric))    
    hc <- fastcluster::hclust(dist(annotX))
    plt <- plt %>% iheatmapr::add_col_dendro(hc, size = 20*ncol(annot)/6 )

    col_annot_height = 11
    plt <- plt %>%
        add_col_annotation(
            ##annotF[colnames(X),,drop=FALSE],
            annotF[,],
            size = col_annot_height, buffer = 0.005, side="bottom",
            colors = colors0, show_title=TRUE)    
    colcex = 1
    if(nrow(annot)<100 && colcex>0) {
        plt <- plt %>% iheatmapr::add_col_labels(
                           side="bottom", size=20*ncol(annot)/6,
                           font=list(size=11*colcex))
    }
    
    return(plt)
}

##annot=ngs$Y
annot.ht=5;cluster.samples=TRUE
pgx.plotPhenotypeMatrix0 <- function(annot, annot.ht=5, cluster.samples=TRUE)
{

    cvar <- pgx.getCategoricalPhenotypes(
        annot, min.ncat=2, max.ncat=10, remove.dup=FALSE)
    cvar

    ##if(length(unique(annot$group)) > 0.33*nrow(annot)) {
    ##    cvar <- setdiff(cvar,"group")
    ##}

    fvar <- pgx.getNumericalPhenotypes(annot)
    fvar
    
    annot.cvar <- annot[,cvar,drop=FALSE]
    annot.fvar <- annot[,fvar,drop=FALSE]
    
    annot.df <- cbind( annot.fvar, annot.cvar )
    colnames(annot.df) <- paste(colnames(annot.df),"        ")

    if(cluster.samples) {
        annotx <- expandAnnotationMatrix(annot.df)
        dim(annot.df)
        dim(annotx)
        hc <- fastcluster::hclust(dist(annotx))  ## cluster samples
        annot.df <- annot.df[ hc$order, ]
    }

    npar <- apply(annot.df,2,function(x) length(setdiff(unique(x),NA)))
    isnum <- c( rep(1,ncol(annot.fvar)), rep(0,ncol(annot.cvar)))
    is.binary <- apply(annot.df,2,function(x) length(setdiff(unique(x),NA))==2)
    is.binary <- apply(annot.df,2,function(x) all(x %in% c(0,1,NA,TRUE,FALSE,"T","F","NA")))
    
    ## set colorscale for each annotation parameter
    ann.colors = list()
    for(i in 1:length(npar)) {
        prm = colnames(annot.df)[i]
        klrs = rev(grey.colors(npar[i],start=0.4,end=0.85))  ## continous scale
        if(npar[i]==1) klrs = "#E6E6E6"
        if(npar[i]>3 & !isnum[i]) klrs = rep(RColorBrewer::brewer.pal(8,"Set2"),99)[1:npar[i]]
        ##if(!is.binary[i] & !isnum[i]) klrs = rep(RColorBrewer::brewer.pal(8,"Set2"),99)[1:npar[i]]        
        ##if(npar[i]==2) klrs = rep(RColorBrewer::brewer.pal(2,"Paired"),99)[1:npar[i]]
        names(klrs) = sort(unique(annot.df[,i]))
        klrs = klrs[!is.na(names(klrs))]
        ann.colors[[prm]] = klrs
    }

    show.legend <- (!is.binary & npar <= 10)
    
    ha = ComplexHeatmap::HeatmapAnnotation(
        df = annot.df,
        ##df2 = annot.fvar[,,drop=FALSE],        
        col = ann.colors,
        na_col = "grey99",
        ##annotation_height = grid::unit(annot.ht, "mm"),
        simple_anno_size = grid::unit(annot.ht,"mm"),  ## BioC 3.8!!
        ##show_annotation_name = (i==ngrp),
        ##annotation_name_gp = grid::gpar(fontsize=3.3*annot.ht),
        ##annotation_legend_param = aa
        show_legend = show.legend
    )
    
    show_colnames = TRUE
    show_colnames = (nrow(annot.df)<100)
    nullmat <- matrix(0,0,nrow(annot.df))
    colnames(nullmat) <- rownames(annot.df)
    colnames(nullmat) <- paste(colnames(nullmat),"     ")
    h <- ComplexHeatmap::Heatmap( nullmat,
                 top_annotation = ha,
                 show_column_names = show_colnames
                 )

    ## some space between heatmap and annotation
    nullmat2 <- matrix(0,0,nrow(annot.df)*0.15)
    h2 <- ComplexHeatmap::Heatmap( nullmat2)
    
    ComplexHeatmap::draw(h + h2,
         padding = grid::unit(c(1,10,1,10),"mm"),
         adjust_annotation_extension = TRUE
         )
    
}

##splitx="cell.family";top.mode="specific";row_annot_width=0.03;ntop=50
pgx.splitHeatmap <- function(ngs, splitx=NULL, top.mode="specific",
                             annot.pheno = NULL, row_annot_width=0.03,
                             scale="row.center", ntop=50, colors=NULL,
                             rowcex=rowcex, colcex=colcex)
{    
    
    
    X0 <- log2(1+ngs$counts)    
    X0[is.na(X0)] <- 0

    splitx
    if(!is.null(splitx) && splitx[1] %in% colnames(ngs$samples)) {
        splitx <- ngs$samples[,splitx]
    } else {
        top.mode = "sd"
        splitx <- NULL
    }

    top.mode
    if(top.mode == "pca") {
        cX <- X0 - rowMeans(X0,na.rm=TRUE)
        dr = svd(cX, nu=5)$u
        ntop1 = ceiling(ntop/ncol(dr))
        jj = as.vector(apply(dr,2,function(x) Matrix::head(order(-abs(x)),ntop1)))
        ##jj = as.vector(apply(dr,2,function(x) Matrix::head(order(-x),10)))
        X1 <- X0[jj,]
        idx <- paste0("PC",as.vector(mapply(rep,1:ncol(dr),ntop1)))
    } else if(top.mode == "specific") {
        grpX <- tapply(colnames(X0), splitx,function(k) rowMeans(X0[,k,drop=FALSE],na.rm=TRUE))
        grpX <- do.call(cbind, grpX)
        cat("dim(grpX)=",dim(grpX),"\n")
        cat("ntop=",ntop,"\n")
        ntop1 = ceiling(ntop/ncol(grpX))
        grpX <- grpX - rowMeans(grpX,na.rm=TRUE)  ## relative to others
        jj = as.vector(apply(grpX,2,function(x) Matrix::head(order(-x),ntop1)))
        X1 <- X0[jj,]
        idx <- paste0("M",as.vector(mapply(rep,1:ncol(grpX),ntop1)))
    } else {
        X1 <- Matrix::head(X0[order(-apply(X0,1,sd)),],ntop)
        hc = fastcluster::hclust( as.dist(1 - stats::cor(t(X1),use="pairwise")), method="ward.D2" )
        idx = paste0("S",cutree(hc, 5))
    }
    table(idx)

    ## ----- Get valid phenotype variables
    if(is.null(annot.pheno)) {
        annot.pheno <- grep("^sample|id$|replicate|ratio|year|month|day",
                            tolower(colnames(ngs$samples)), invert=TRUE)
        annot.pheno <- pgx.getCategoricalPhenotypes(ngs$samples, max.ncat=12, min.ncat=2)
    }
    annot.pheno <- intersect(annot.pheno, colnames(ngs$samples))
    Y <- ngs$samples[,annot.pheno,drop=FALSE]
    Y <- data.frame(apply(Y,2,as.character))
    rownames(Y) <- rownames(ngs$samples)
    colnames(Y)

    sampletips = colnames(X1)
    genetips   = rownames(X1)
    
    ## ----- call plotting function
    plt <- pgx.splitHeatmapFromMatrix(
        X=X1, annot=Y,
        xtips=genetips, ytips=sampletips,
        idx=idx, splitx=splitx,
        row_annot_width=row_annot_width,
        scale=scale, colors=colors,
        rowcex=rowcex, colcex=colcex)
    return(plt)
}

##=================================================================================
## Lower level R level plotting functions
##=================================================================================

plot.gghist <- function(x) {
    
    
    
    dplyr::as_tibble(x) %>%
        tidyr::pivot_longer(
            cols=everything(),
            names_to="cell",
            values_to="expression"
        ) %>%
        ggplot2::ggplot(ggplot2::aes(x=expression, group=cell)) +
        ggplot2::geom_density() +
        ggplot2::coord_cartesian(ylim=c(0,0.6), xlim=c(0,3))
}

plotly2ggplot <- function (plot, width=NULL, height=NULL, scale=1, hjust=0, vjust=0) 
{
    
    
    
    tmpfile <- tempfile()
    tmpfile <- "tmp.png"
    unlink(tmpfile)
    unlink("tmp.png")
    unlink("tmp_1.png")
    try(plotly::orca(plot, file=tmpfile, width=width, height=height, scale=2))
    img <- png::readPNG("tmp_1.png")
    img.grob <- grid::rasterGrob(img, interpolate=TRUE)
    ymin <- xmin <- 1 - scale
    xmax <- ymax <- scale
    message("converting to ggplot...")
    gg <- ggplot2::ggplot(data.frame(x = 0:1, y = 0:1), ggplot2::aes_(x = ~x, y = ~y)) + 
        ggplot2::geom_blank() + ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 
        0)) + ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 
        0)) + ggplot2::annotation_custom(img.grob, xmin = xmin + 
        hjust, xmax = xmax + hjust, ymin = ymin + vjust, ymax = ymax + 
        vjust) + ggplot2::theme_void()
    return(gg)
}

##cex=1;main=NULL;xlab=NULL;ylab=NULL;yth=1;tooltips=NULL;cex.text=1
gsea.enplotly <- function(fc, gset, cex=1, main=NULL, xlab=NULL,
                          ylab=NULL, yth=1, tooltips=NULL, cex.text=1)
{
    if(is.null(xlab))
        xlab <- "Rank in ordered dataset"
    if(is.null(ylab))
        ylab <- "Ranked list metric"
    if(is.null(main))
        main <- "Enrichment plot"

    ## compute running metrix
    fc <- sort(fc,decreasing=TRUE)
    ##in.gset <- 1*(names(fc) %in% gset)
    ##rnk.trace <- cumsum( c(-1/sum(in.gset==0), 1/sum(in.gset))[1+in.gset])
    ##rnk.trace <- 0.66 * rnk.trace / max(abs(rnk.trace)) * max(abs(fc))

    ## weighted cumulative random walk
    x0 <- 1*(names(fc) %in% gset)
    x0 <- x0 * abs(fc)
    n0 <- sum(!(names(fc) %in% gset))
    n1 <- sum(names(fc) %in% gset)
    r0 <- cumsum(x0==0) / sum(x0==0)
    r1 <- cumsum(x0) / (1e-4+sum(x0))
    ##mx <- quantile(abs(fc),probs=1)
    rnk.trace <- (r1 - r0)
    rnk.trace <- rnk.trace / max(abs(rnk.trace)) * 0.8

    qq <- quantile(fc,probs=c(0.005,0.995),na.rm=TRUE)
    qq <- range(fc)
    ##qq <- c(min(fc,na.rm=TRUE),max(fc,na.rm=TRUE)) * 0.8
    y1 <- qq[2]
    y0 <- qq[1]
    dy <- 0.2*(y1-y0)
    if(max(rnk.trace) >= abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y1)
    if(max(rnk.trace) < abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y0)
    
    cc <- sign(fc)*rank(abs(fc))
    df <- data.frame( x=rank(-fc), y=fc, trace=rnk.trace, cc=cc)
    
    ## downsample
    ii <- which(rownames(df) %in% gset)
    ii <- unique(c(ii,seq(1,nrow(df),round(nrow(df)/1000))))
    df <- df[ii,]
    df <- df[order(-df$y),]
    dim(df)
    
    ##dy <- 0.05*diff(range(fc))
    ##y0 <- min(fc)
    cpal <- colorspace::diverge_hcl(64, c=60, l=c(30,100), power=1)
    cpal <- colorspace::diverge_hcl(64)
    
    ## colorbar segments
    db <-  nrow(df)/11
    bb <- round(seq(1, nrow(df), db))
    cbar.x    <- df$x[bb]
    cbar.xend <- df$x[c(bb[-1],nrow(df))]

    fx   <- df$y[bb + round(0.5*db)]
    ##r <- mean(rnk[kk[c(i,i+1)]])
    r1 <- (fx/max(abs(fx),na.rm=TRUE))
    r1 <- abs(r1)**0.66 * sign(r1)
    suppressWarnings(
        cc <- gplots::colorpanel(21,"royalblue3","grey90","indianred3")
    )
    irnk <- 1 + round((length(cc)-1)*(1 + r1)*0.5)
    ##cat("irnk=",irnk,"\n")    
    cbar <- data.frame( x=cbar.x, xend=cbar.xend, color=cc[irnk])

    cex.title=1
    df$text = rownames(df)
    jj <- which(rownames(df) %in% gset)    
    tooltips2 = rownames(df)[jj]
    if(!is.null(tooltips)) {
        sel.tt <- match(rownames(df)[jj],names(tooltips))
        tooltips2 <- paste0('<b>',tooltips2,'</b><br>',tooltips[sel.tt])
    }
    
    ii = seq(1,nrow(df),round(nrow(df)/200))
    fig <- plotly::plot_ly() %>%
        plotly::add_trace(
            x = ~df$x[ii], y = ~df$y[ii],
            type = 'scatter', fill = 'tozeroy',
            fillcolor = '#BBBBBB',
            hoverinfo = 'skip',            
            mode = 'none' )  %>%
        plotly::add_trace(
            x = ~df$x, y = ~df$trace,
            type = 'scatter', mode = 'lines',
            hoverinfo = 'skip',
            line = list(
                color = '#00EE00',
                width = 4
            ))  %>%
        plotly::add_trace(
            ## orange points of genes
            x = ~df$x[jj], y = ~df$y[jj],
            type = 'scatter', ## fill = 'toself', ## mode = 'none',
            marker = list( color = '#FF8C00', size=cex*6 ),
            text = tooltips2,
            hoveron = 'points', hoverinfo='text') %>%
        plotly::add_segments(
            x = df$x[jj], xend = df$x[jj],
            y = y0 - 0.98*dy, yend = y0,
            type = 'scatter', mode='lines',
            line = list(color = '#444444', width=1.5*cex),
            text = rownames(df)[jj],
            hoveron = 'points', hoverinfo='text'
        )

    ## colorbar/color scale
    for(i in 1:nrow(cbar)) {
        fig <- fig %>%
            plotly::add_segments(
                x = cbar$x[i], xend = cbar$xend[i],
                y = y0 - 0.9*dy, yend = y0 - 0.9*dy,
                ##type = 'scatter', mode='lines',
                line = list(color=cbar$color[i], width=32)
        )
    }

    if(!is.null(yth) && yth>=0 ) {
        ## static labels
        jj2 = jj[which( abs(df$y[jj]) >= yth )]
        if(length(jj2)) {
            fig <- fig %>%
                plotly::add_annotations(
                    x = df$x[jj2],
                    y = df$y[jj2],
                    text = rownames(df)[jj2],
                    xref = 'x',
                    yref = 'y',
                    xanchor = 'left',
                    xshift = 4,
                    yshift = 3,
                    showarrow = FALSE
                )
        }
    }
    
    fig <- fig %>%
        plotly::layout(
            font = list(size=12*cex.text),
            title = list(text=main, y=0.99),
            xaxis = list(title = xlab, gridwidth=0.3),
            yaxis = list(title = ylab, gridwidth=0.3,
                         range=c(y0-1.1*dy,y1) )) %>% 
        plotly::config(toImageButtonOptions = list(format = "svg")) %>% 
        plotly::hide_legend()
    
    fig
    
}

ggenplot <- function(fc, gset, cex=1, main=NULL, xlab=NULL, ylab=NULL)
{
    if(is.null(xlab))
        xlab <- "Rank in ordered dataset"
    if(is.null(ylab))
        ylab <- "Ranked list metric"
    if(is.null(main))
        main <- "Enrichment plot"

    ## compute running metrix
    fc <- sort(fc,decreasing=TRUE)
    ##in.gset <- 1*(names(fc) %in% gset)
    ##rnk.trace <- cumsum( c(-1/sum(in.gset==0), 1/sum(in.gset))[1+in.gset])
    ##rnk.trace <- 0.66 * rnk.trace / max(abs(rnk.trace)) * max(abs(fc))

    ## weighted cumulative random walk
    x0 <- 1*(names(fc) %in% gset)
    x0 <- x0 * abs(fc)
    n0 <- sum(!(names(fc) %in% gset))
    n1 <- sum(names(fc) %in% gset)
    r0 <- cumsum(x0==0) / sum(x0==0)
    r1 <- cumsum(x0) / (1e-4+sum(x0))
    ##mx <- quantile(abs(fc),probs=1)
    rnk.trace <- (r1 - r0)
    rnk.trace <- rnk.trace / max(abs(rnk.trace)) * 0.8

    qq <- quantile(fc,probs=c(0.005,0.995),na.rm=TRUE)
    qq <- range(fc)
    ##qq <- c(min(fc,na.rm=TRUE),max(fc,na.rm=TRUE)) * 0.8
    y1 <- qq[2]
    y0 <- qq[1]
    dy <- 0.2*(y1-y0)
    if(max(rnk.trace) >= abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y1)
    if(max(rnk.trace) < abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y0)
    
    cc <- sign(fc)*rank(abs(fc))
    df <- data.frame(rank=rank(-fc),fc=fc,run=rnk.trace,cc)
    jj <- which(names(fc) %in% gset)
    ##dy <- 0.05*diff(range(fc))
    ##y0 <- min(fc)
    cpal <- colorspace::diverge_hcl(64, c=60, l=c(30,100), power=1)
    ii <- 1 + 32 + range(round(32*(fc / max(abs(fc)))))
    cpal <- colorspace::diverge_hcl(65)[ii[1]:ii[2]]
    
    cex.title=1

    ggplot2::ggplot(data=df, ggplot2::aes(x=rank, y=fc, color=cc)) +
        ggplot2::geom_segment(ggplot2::aes(x=rank, y=0, xend=rank, yend=fc), color="grey70") +
        ggplot2::geom_line(ggplot2::aes(x=rank, y=run), color="green", size=0.8) +
        ggplot2::geom_segment(data=df[jj,], color="black", size=0.3*cex,
                     ggplot2::aes(x=rank, y=(y0-dy), xend=rank, yend=(y0))) +
        ggplot2::geom_segment(ggplot2::aes(x=rank, y=y0-1*dy, xend=rank, yend=y0-0.6*dy)) +
        ggplot2::scale_color_gradientn( colors=cpal ) +
        ggplot2::ylim(y0-1.1*dy,y1) + 
        ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(main) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position='none',
            plot.title = ggplot2::element_text(size=9*cex.title),
            axis.text.x = ggplot2::element_text(size=8, vjust=+5),
            axis.text.y = ggplot2::element_text(size=8, hjust=0),
            axis.title.x = ggplot2::element_text(size=9, vjust=+5.5),
            axis.title.y = ggplot2::element_text(size=9, hjust=+0.5)
        ) 

}

plot.ggsplom <- function(F, F2=NULL, title_cex=2, no.axes=FALSE, ...)
{
    if(is.null(F2)) {
        F2 <- F
    }
    lim0 <- range(F)
    lim1 <- range(F2)
    symm <- all(colnames(F)==colnames(F2))

    gg <- intersect(rownames(F),rownames(F2))
    F  <- F[gg,]
    F2 <- F2[gg,]    

    plt <- list()
    k=1
    for(i in 1:ncol(F)) {
        for(j in 1:ncol(F2)) {
            if(i==j && symm) {
                p1 <- plot.ggscatter(0,0,cex=0) +
                    ggplot2::theme_void() +
                    ggplot2::geom_text(x=0,y=0,label=colnames(F)[i],size=title_cex)
            } else {
                p1 <- plot.ggscatter(F[,i], F[,j] ) +
                ##p1 <- plot.ggscatter(F[,i], F2[,j], ... ) +
                    ggplot2::xlim(lim0) + ggplot2::ylim(lim1) +
                    ggplot2::xlab(colnames(F)[i]) +
                    ggplot2::ylab(colnames(F2)[j]) 
            }
            plt[[k]] <- p1
            k <- k + 1
        }
    }

    blankx <- ggplot2::theme(
        plot.margin = ggplot2::margin(0,0,0,0),
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
    )
    blanky <- ggplot2::theme(
        plot.margin = ggplot2::margin(0,0,0,0),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
    )

    np <- length(plt)
    nc <- ncol(F)
    nr <- ncol(F2)    
    border.y <- ((0:(np-1) %% nc)==0)
    border.x <- ((0:(np-1) %/% nc)==(nr-1))
    for(i in 1:np) {
        if(no.axes || !border.y[i]) plt[[i]] <- plt[[i]] + blanky
        if(no.axes || !border.x[i]) plt[[i]] <- plt[[i]] + blankx
    }
    cowplot::plot_grid(
                 plotlist=plt, labels=NA,
                 nrow=nr, ncol=nc, 
                 rel_widths = c(1.2,rep(1,nc-1)),
                 rel_heights = c(rep(1,nr-1),1.2) )
    
}

##col=NULL;pch=20;xlab="";ylab="";cex=1;opacity=1;main=""
plot.ggscatterFILL <- function(x, y=NULL, col=NULL, shape=NULL,
                          main=NULL, cex=1,
                          pch=20, legend=TRUE,  xlab=NULL, ylab=NULL,
                          legend.ysp=0.8, cex.legend=1,
                          barscale=0.3, opacity=1, gamma=1) {
    if(is.null(y) && NCOL(x)==2) {
        xlab <- colnames(x)[1]
        ylab <- colnames(x)[2]
        y <- x[,2]
        x <- x[,1]
    }
    df <- data.frame(x=x, y=y)
    if(!is.null(col)) df$col <- col
    if(!is.null(shape)) df$shape <- shape
    Matrix::head(df)    
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color=col, shape=shape)) +
        ggplot2::geom_point( shape=pch, alpha=opacity, size=2.0*cex ) +
        ggplot2::ggtitle(main) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    p
    if(!is.null(col)) {
        cpal <- rev(RColorBrewer::brewer.pal(11,"RdYlBu"))
        if(0 && opacity<1) {
            cpal <- add_opacity(cpal, opacity**0.33)
        }
        zr <- range(col)
        zz <- round(c(zr[1], zr[2]),digits=2)
        cgamma <- seq(0,1,1/(length(cpal)-1))**(1/gamma)
        p <- p +
            ggplot2::guides(colour = ggplot2::guide_colourbar(
                       barwidth= 1.5*barscale, barheight= 8*barscale),
                   shape = ggplot2::guide_legend(override.aes = list(size=2.5*cex.legend))
                   ) +
            ggplot2::scale_color_gradientn(
                colors=cpal, breaks=zz, values = cgamma, 
                labels=c(zz[1],zz[2]) ) +
            ggplot2::expand_limits( color = zr + c(-0.01,0.01))
    }
    if(legend) {
        p <- p + ggplot2::theme(
                     legend.justification = c(0,0),
                     legend.position=c(0.01,0.01)
                 )
    } else {
        p <- p + ggplot2::theme(legend.position="none")
    }
    p <- p + ggplot2::theme(legend.title = ggplot2::element_blank())          
    p
}

plot.ggscatter <- function(x, y=NULL, col=NULL, main=NULL, 
                      cex=1, col.scale=NULL, shape=NULL, pch=20, 
                      legend=TRUE, legend.ysp=0.8, cex.legend=1,
                      legend.pos = "right",
                      xlab=NULL, ylab=NULL, base_size=12)
{
    if(is.null(y) && NCOL(x)==2) {
        y <- x[,2]
        x <- x[,1]
    }
    if(is.null(col.scale)) col.scale <- rep(1:9,99)
    legend.justification <- NULL
    legend.position <- legend.pos
    if(legend.position=="bottomright") {
        legend.justification <- c(1,0)
        legend.position <- c(1.0, 0.0)
    } else if(legend.position=="topleft") {
        legend.justification <- c(0,1)
        legend.position <- c(.0, 1.0)
    } else if(legend.position=="topright") {
        legend.justification <- c(1,1)
        legend.position <- c(1.0, 1.0)
    } else if(legend.position=="bottomleft") {
        legend.justification <- c(0,0)
        legend.position <- c(.0, 0.0)
    }

    df <- data.frame(x=x, y=y)
    if(!is.null(col)) df$col <- col
    if(!is.null(shape)) df$shape <- shape
    Matrix::head(df)
    is.factor <- class(type.convert(as.character(col)))=="factor"
    if(is.factor) {
        p <- ggplot2::ggplot(df, ggplot2::aes(y=y, x=x, color=col, shape=shape)) +
            ggplot2::geom_point(size = 2.0*cex) +
            ggplot2::scale_color_manual( values=col.scale, name='' ) +
            ggplot2::ggtitle(main) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    } else {
        p <- ggplot2::ggplot(df, ggplot2::aes(y=y, x=x, color=col, shape=shape)) +
            ggplot2::geom_point(size = 2.0*cex) +
            ggplot2::ggtitle(main) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
        if(!is.null(col.scale)) {
            p <- p + ggplot2::scale_color_gradient( low=col.scale[1], high=col.scale[2] ) 
        }
    }
    p
    
    if(legend) {
        p <- p +  ggplot2::theme(
                      legend.title = ggplot2::element_blank(),
                      legend.text = ggplot2::element_text(size=9*cex.legend),
                      legend.key.size = grid::unit(legend.ysp*0.8,"lines"),
                      legend.key = ggplot2::element_rect(color="transparent", fill=scales::alpha("white",0.5)),
                      legend.justification = legend.justification,
                      legend.position = legend.position,
                      legend.margin = ggplot2::margin(1, 3, 2, 1),
                      legend.box.just = "right",
                      legend.box.background = ggplot2::element_rect(color="#888888", size=0.25),
                      legend.box.margin = ggplot2::margin(1,2,1,1)
                  ) +
            ggplot2::guides(
                color = ggplot2::guide_legend(override.aes = list(size=2.0*cex.legend)),
                shape = ggplot2::guide_legend(override.aes = list(size=2.0*cex.legend)) )
        
    } else {
        p <- p + ggplot2::theme(legend.position="none")
    }
    ## p <- p + ggplot2::theme_minimal(base_size=base_size)    
    p
}

##srt=30;xlab=ylab="";n.dodge=1;cex=1;ylab=xlab=""
plot.ggviolin <- function(x, y, group=NULL, main="", ylim=NULL, add.dots=TRUE,
                     col="#AAAAAA", cex=1, xlab="", ylab="y", srt=0,
                     pdodge=1.5, n.dodge=1, base_size=13)
{

    df <- data.frame(y=y, x=x, group="")
    if(!is.null(group)) {
        df$group <- group
        col <- RColorBrewer::brewer.pal(8,"Set3")
    }
    if(is.null(ylim)) ylim <- range(y)
    
    p <- ggplot2::ggplot(df,aes(y=y, x=x, fill=group)) +
        ggplot2::ggtitle(main) +
        ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + 
        ggplot2::ylim( ylim[1], ylim[2] ) +
        ggplot2::scale_fill_manual(values=col) +
        ggplot2::scale_x_discrete(guide=guide_axis(angle=srt)) +
        ## ggplot2::scale_x_discrete(guide=guide_axis(n.dodge=n.dodge)) +
        ggplot2::geom_violin(trim=TRUE, position=position_dodge(pdodge)) +        
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=srt, vjust=0)) +
        ggplot2::theme_minimal(base_size=base_size)
    if(is.null(group)) {
        p <- p + ggplot2::theme(legend.position="none") 
    }
    if(add.dots && is.null(group) ) {
        p <- p +
            ggplot2::geom_jitter(shape=20, size=1.2*cex, position=position_jitter(0.07))
        ## geom_beeswarm(priority='density',cex=0.9*cex) +
    }
    p
}



##group.name="group";xlab="x";ylab="y";srt=0;main=NULL;base_size=14
plot.ggbarplot <- function(mat, xlab="x", ylab="y", srt=0, main=NULL, 
                      las=NULL, col=NULL, beside=FALSE,
                      legend.pos = c(0.016,1), legend.cex = 1,
                      bar_width=0.7, base_size=12, group.name="group")
{
    require(ggplot2)
    
    if(NCOL(mat)==1) mat <- rbind(mat)
    mat <- mat[nrow(mat):1,,drop=FALSE]
    df <- reshape2::melt(t(mat), value.name = "value")
    colnames(df)[1:2] <- c("x","y")
    ##df$y <- paste0(df$y," ")
    ##df$x <- paste0(" ",df$x," ")
    df$y <- factor(df$y, levels=rownames(mat))
    df$x <- factor(df$x, levels=colnames(mat))
    if(!is.null(las) && las==3) srt <- 90

    cpal <- rev(grey.colors(nrow(mat)))
    ## cpal <- grey.colors(nrow(mat))
    if(nrow(mat)==1) cpal <- "grey70"
    if(!is.null(col)) cpal <- rep(col,99)
    posmode <- ifelse(beside, "dodge", "stack")
    
    p <- ggplot(df, aes(x=x, y=value, fill=y)) + 
        ##geom_bar(stat="identity", aes(fill="transparent"), size=1) +
        geom_bar(stat="identity", color="black", size=0.3,
                 width=bar_width, position=posmode) + 
        xlab(xlab) + ylab(ylab) + labs(fill=group.name) +
        ggtitle(main) +
        ##scale_fill_grey(start=0.8,end=0.2) +
        scale_fill_manual(values=cpal) +
        theme_classic(base_size=base_size) +
        scale_x_discrete(guide=guide_axis(angle=srt)) +
        ## scale_x_discrete(guide=guide_axis(n.dodge=n.dodge)) + 
        theme(
            axis.text.x = element_text(angle=srt, vjust=0),
            axis.title.x = element_text(size=10),
            axis.title.y = element_text(size=10)
        ) 
    
    p <- p + theme(
                 legend.title = element_blank(), 
                 legend.justification = legend.pos,
                 legend.text = element_text(size=9*legend.cex),
                 legend.position = legend.pos,
                 legend.key.size = grid::unit(9*legend.cex, "pt"),
                 legend.key.height = grid::unit(7*legend.cex, "pt"))

    if(nrow(mat)==1) {
        p <- p + theme(legend.position = 'none')
    }
    p
}

##jitter=TRUE;vcol="lightcyan2";maxbee=200;xth=0
pgx.violinPlot <- function(x, y, group=NULL, xlab='', ylab='',
                           srt=0, cex.lab=1, cex.main=1.1,
                           jitter=0.015, vcol="grey85",
                           plotlib='base', maxbee=NULL, ...)
{
    fig = NULL
    if(plotlib=='base') {
        
        
        if(is.null(maxbee)) maxbee <- 100*length(setdiff(unique(y),NA))
        if(jitter>0) x <- x + jitter*diff(range(x))*rnorm(length(x))
        ## vioplot::vioplot(x ~ y, col=vcol, plotCentre="line", xaxt='n' )
        ## vioplot::vioplot(x ~ y, col=vcol, plotCentre="line", xaxt='t' )        
        vioplot::vioplot(x ~ y, col=vcol, ## cex.main=cex.main,
                plotCentre="line", xaxt='n', xlab=xlab, ylab=ylab,
                ... )
        if(maxbee>0) {
            ii <- Matrix::head(sample(length(x)), maxbee)
            beeswarm::beeswarm(x[ii] ~ y[ii], add=TRUE, pch=20, cex=0.6, col="grey10")
        }
        yy <- sort(unique(y))
        text(x = 1:length(yy),
             y = par("usr")[3] - 0.03*diff(range(x)),
             labels = yy,
             xpd = NA,
             srt = srt,
             adj = ifelse(srt==0, 0.5, 0.965),
             cex = cex.lab)        
    }
    if(plotlib=='ggplot') {
        d1 <- !is.null(maxbee) && maxbee>0
        fig <- plot.ggviolin( y, x, add.dots=d1, xlab=xlab, ylab=ylab, srt=srt,
                             main=main, cex.main=cex.main)        
    }
    if(plotlib=='plotly') {
        df = data.frame(x=x, y=y, group=0)
        if(!is.null(group)) df$group <- group
        fig <- df %>%
            plotly::plot_ly(
                x = ~y,
                y = ~x,
                ## split = ~group,
                type = 'violin',
                box = list(
                    visible = TRUE
                ),
                meanline = list(
                    visible = TRUE
                )
            ) %>%
            plotly::layout(
                font = list(size=16*cex.lab),
                xaxis = list(
                    title = xlab
                ),
                yaxis = list(
                    title = ylab,
                    zeroline = FALSE
                )
            )
    }
    fig
}

pgx.scatterPlotXY <- function(..., plotlib="base") {
    if(plotlib=="plotly") {
        pgx.scatterPlotXY.PLOTLY(...)
    } else if(plotlib=="ggplot") {
        pgx.scatterPlotXY.GGPLOT(...)
    } else if(plotlib=="scatterD3") {
        pgx.scatterPlotXY.D3(...)
    } else {
        pgx.scatterPlotXY.BASE(...)
    }
}
 
pgx.scatterPlotXY.BASE <- function(pos, var=NULL, type=NULL, col=NULL, title="",
                                   zlim=NULL, zlog=FALSE, zsym=FALSE, softmax=FALSE, pch=20,
                                   cex=NULL, cex.lab=1, cex.title=1.2, cex.legend=1,
                                   zoom=1, legend=TRUE, bty='o', legend.ysp=0.85,
                                   legend.pos = 'bottomleft', lab.pos=NULL, repel=TRUE,
                                   xlab = NULL, ylab=NULL, xlim=NULL, ylim=NULL, dlim=0.05,
                                   hilight2=hilight, hilight.cex = NULL, lab.xpd=TRUE,
                                   hilight=NULL, hilight.col=NULL, hilight.lwd=0.8,
                                   label.clusters=FALSE, cex.clust=1.5,
                                   tstep=0.1, rstep=0.1,
                                   tooltip=NULL, theme=NULL, set.par=TRUE, 
                                   axt='s', xaxs=TRUE, yaxs=TRUE, 
                                   labels=NULL, label.type=NULL, opacity=1)    
{
    
    if(0) {
        var=NULL; type=NULL; col=NULL; title="";
        zlim=NULL; zlog=FALSE; zsym=FALSE; softmax=FALSE; pch=20;
        cex=NULL; cex.lab=1; cex.title=1.2; cex.legend=1;
        zoom=1; legend=TRUE; bty='o'; legend.ysp=0.85;
        legend.pos = 'bottomleft'; lab.pos=NULL; repel=TRUE;
        xlab = NULL; ylab=NULL; xlim=NULL; ylim=NULL; dlim=0.05;
        hilight2=hilight; hilight.cex = NULL; lab.xpd=TRUE;
        hilight=NULL; hilight.col=NULL; hilight.lwd=0.8;
        label.clusters=FALSE; cex.clust=1.5; rstep=0.1;
        tooltip=NULL; theme=NULL; set.par=TRUE; 
        axt='s'; labels=NULL; label.type=NULL; opacity=1
    }

## automatically set pointsize of dots
    if(is.null(cex)) {
        nr <- nrow(pos)
        i <- as.integer(cut(nr,breaks=c(0,100,500,1000,5000,Inf)))
        cex  <- c(2,1.4,1,0.7,0.4)[i]
    }
    if(is.null(hilight.cex)) hilight.cex = cex
    if(!is.null(var) && !is.null(ncol(var))) {
        var <- array(var[,1], dimnames=list(rownames(var)))
    }
    if(!is.null(var) && length(var)==nrow(pos) && is.null(names(var)) ) {
        names(var) <- rownames(pos)
    }
    if(is.null(var)) {
        var <- as.character(rep("_",nrow(pos)))
        names(var) <- rownames(pos)
    }
    if(is.null(type)) {
        type <- c("numeric","factor")[1 + class(var) %in% c("factor","character")]
    }
    if(is.null(colnames(pos))) {
        colnames(pos) <- c("x","y")
    }
    var <- var[match(rownames(pos),names(var))]
    
    ## normalize pos
    xlim0 <- range(pos[,1])
    ylim0 <- range(pos[,2])
    if(zoom!=1) {
        cx <- mean(range(pos[,1]))
        cy <- mean(range(pos[,2]))
        dx <- diff(range(pos[,1]))
        dy <- diff(range(pos[,2]))
        xlim0 <- cx + 0.5 * c(-1,1.05) * dx / zoom
        ylim0 <- cy + 0.5 * c(-1,1.05) * dy / zoom
    }

    if(length(dlim)==1) dlim <- rep(dlim,2)
    xlim0[1] <- xlim0[1] - dlim[1]*diff(xlim0)
    xlim0[2] <- xlim0[2] + dlim[1]*diff(xlim0)
    ylim0[1] <- ylim0[1] - dlim[2]*diff(ylim0)
    ylim0[2] <- ylim0[2] + dlim[2]*diff(ylim0)    

    if(is.null(xlab)) xlab <- colnames(pos)[1]
    if(is.null(ylab)) ylab <- colnames(pos)[2]

    if(set.par) {
        par.save <- par()
        ##par(mar=c(0.4,0.3,1.7,0.4), mgp=c(1.2,0.4,0),
        par(mar=c(2.6,2.8,1.9,0.5), mgp=c(1.5,0.4,0),        
            cex.axis=1, cex.lab=1, tcl=-0.3, las=1)
    }

    if(is.null(labels))
        labels <- rownames(pos)
    
    ## Plot the discrete variables
    if(type=="factor") {
        
        z1 <- factor(var)
        nz <- length(levels(z1)) 
        if(is.null(col) && nz>2) {
            col1 <- c(RColorBrewer::brewer.pal(8,"Set1"),
                      RColorBrewer::brewer.pal(8,"Set2"),
                      RColorBrewer::brewer.pal(12,"Set3"))    
        } else if(is.null(col) && nz==2) {
            col1 <- rev(grey.colors(2, end=0.8))
            col1 <- c("#AAAAAA55","#555555FF")
            col1 <- c("#00008855","#AA0000FF")
            col1 <- c("#CCCCCC55","#AA0000FF")
        } else if(is.null(col) && nz==1) {
            col1 <- c("#22222255")
        } else {
            col1 <- col
        }
        col1 <- Matrix::head(rep(col1,99),nz)
        pt.col <- col1[z1]
        pt.col[is.na(pt.col)] <- "#DDDDDD33"
        pt.col0 = pt.col
        if(opacity<1) {
            pt.col <- add_opacity(pt.col, opacity)
            col1 <- add_opacity(col1, opacity**0.33)
        }

        jj <- order(-table(pt.col)[pt.col]) ## plot less frequent points last...            
        if(length(cex)>1) cex <- cex[jj]
        plot(pos[jj,,drop=FALSE],
             col = pt.col[jj], pch=20, cex=cex,
             xlim = xlim0, ylim = ylim0,
             xlab = xlab, ylab = ylab,
             ##xaxt = axt, yaxt = axt,
             axes = FALSE,
             bty = 'n' )
        if(bty!='n') graphics::box(lwd=0.8, bty=bty, col="black")
        graphics::Axis(side=1, labels=xaxs)
        graphics::Axis(side=2, labels=yaxs)
        graphics::grid(lwd=0.8)
        
        ## label cluster
        if(label.clusters) {
            mpos <- apply( pos,2,function(x) tapply(x,z1,median))
            mlab <- rownames(mpos)
            ##if(!is.null(labels)) mlab <- labels[rownames(mpos)]
            text(mpos, labels=mlab, cex=cex.clust, lheight=0.8)
        }
                
        ## discrete parameter legend
        nlev <- length(levels(z1))
        if(legend.pos!='none' && legend && nlev < 30) {
            cex1 <- ifelse(length(levels(z1))>=8,0.85,1)
            cex1 <- ifelse(length(levels(z1))>=15,0.75,cex1)
            legend(legend.pos, levels(z1), bty=bty, fill=col1,
                   ## cex=cex1, y.intersp=0.75, x.intersp=0.7
                   cex = 0.85*cex1*cex.legend,
                   inset = c(0.01,0.02),
                   y.intersp=legend.ysp, x.intersp=0.75
                   )
        }
    }
    
    ## Plot continous variable    
    if(type=="numeric") {
        ##if(NCOL(cvar)==1) cvar <- cbind(cvar=cvar)
        z <- var
        
        if(is.null(zlim)) {
            ## global zlim
            zlim <- range(z, na.rm=TRUE)
            ##zlim <- quantile(z,probs=c(0.01,0.99),na.rm=TRUE)
            ##zlim <- quantile(z,probs=c(0.002,0.998),na.rm=TRUE)        
        }
        
        if(zsym) {
            zlim <- c(-1,1)*max(abs(zlim),na.rm=TRUE)
        }
        
        ## z1 is normalized [0;1] for coloring        
        z1 <- (z - min(zlim)) / diff(zlim)
        z1 <- pmin(pmax(z1,0),1) ## clip
        if(softmax) {
            z1 <- 0.5*(tanh(4*(z1-0.5))+1)
        }
        
        ##-------------- set colors
        ##cpal <- rev(viridis::viridis(11))
        cpal <- colorspace::diverge_hcl(64, c=60, l=c(30,100), power=1)
        cpal <- colorRampPalette(c("#313695","#FFFFDF","#A50026"))(64)
        cpal <- colorRampPalette(c("#3136B5","#FFFFDF","#B50026"))(64)                        
        cpal <- sapply(1:length(cpal),function(i) add_opacity(cpal[i],0.2+0.8*abs(i-32.5)/32))        
        ##cpal <- rev(RColorBrewer::brewer.pal(11,"RdYlBu"))
        pt.col <- cpal[1+ceiling(z1*(length(cpal)-1))]
        ##pt.col <- viridis::viridis(16)[1+round(15*z1)]
        ##pt.col <- matlab.like(21)[1+round(20*z1)]        
        pt.col0 = pt.col
        if(opacity<1) {
            pt.col <- add_opacity(pt.col, opacity)
            cpal <- add_opacity(cpal, opacity**0.33)                
        }        
        pt.col[is.na(pt.col)] <- "#DDDDDD33"  ## missing values color

        jj <- 1:nrow(pos)
        jj <- order(abs(z),na.last=FALSE) ## higher values last??        
        if(length(cex)>1) cex <- cex[jj]
        plot(pos[jj,],
             col=pt.col[jj],
             pch=20, cex=cex,
             xlim = xlim0, ylim=ylim0, 
             xlab = xlab, ylab = ylab,
             ##xaxt=axt, yaxt=axt,
             axes = FALSE, ## later
             bty = 'n' )
        if(bty!='n') graphics::box(lwd=0.8, bty=bty, col="black")
        graphics::Axis(side=1, labels=xaxs)
        graphics::Axis(side=2, labels=yaxs)
        graphics::grid(lwd=0.8)

        ## colorscale bar
        if(legend.pos!='none' && legend) {
            zr <- range(z,na.rm=TRUE)
            if(!is.null(zlim)) zr <- zlim
            if(zsym) zr <- c(-1,1)*max(abs(zr),na.rm=TRUE)
            if(zlog) zr <- round(10**zr-1)  ##???
            zr <- 0.01 * c(ceiling(100*zr[1]), floor(100*zr[2])) ## round
            legend(legend.pos, cex=0.8*cex.legend, ## text.width=2,
                   y.intersp=0.18, x.intersp=0.5, border=NA, bty=bty,
                   inset = c(0.01,0.02),
                   fill = rev(cpal[seq(from=1,to=length(cpal),length.out=11)]),
                   legend = c(zr[2],rep(NA,9),zr[1]), )
        }

    }
    
    if(!is.null(hilight) && length(hilight)>0) {
        jj <- which(rownames(pos) %in% hilight)
        if(length(jj)) {
            hcol1 = hilight.col
            if(is.null(hcol1)) hcol1 <- pt.col0[jj]
            hcex <- hilight.cex
            if(length(hcex)>1) hcex <- hcex[jj]
            points(pos[jj,,drop=FALSE], pch=20, col=hcol1, cex=1.05*hcex)
            points(pos[jj,,drop=FALSE], pch=1, lwd=hilight.lwd, cex=0.85*hcex)
        }
    }

    if(!is.null(hilight2) && length(hilight2)>0) {        
        jj <- which(rownames(pos) %in% hilight2)
        if(length(jj)) {
            if(length(cex.lab)==1) cex.lab <- rep(cex.lab,nrow(pos))
            df <- data.frame(x=pos[jj,1], y=pos[jj,2], z=labels[jj], cex=cex.lab[jj])
            if(is.null(lab.pos)) {
                ## repelling text using wordcloud package
                df2 <- rbind(df, df)
                df2$z[1:nrow(df)] <- "x"
                xlim1 = ylim1 = c(-Inf, Inf)
                if(!lab.xpd) {
                    xlim1 = xlim0
                    ylim1 = ylim0
                }
                cex.lab2 <- c(rep(1,nrow(df)),df$cex)
                if(repel) {
                    dbg("[pgx.scatterPlotXY.BASE] 6: repelling labels...")
                    ##nc <- wordcloud::wordlayout(df2$x, df2$y, df2$z, 
                    nc <- repelwords(df2$x, df2$y, df2$z, 
                                     xlim = xlim1, ylim = ylim1,
                                     tstep = tstep, rstep = rstep,
                                     cex = cex.lab2)
                    nc <- Matrix::tail(nc,nrow(df))                    
                    lab.pos <- data.frame(x=nc[,1]+.5*nc[,3], y=nc[,2]+.5*nc[,4])
                } else {
                    ##dbg("[pgx.scatterPlotXY.BASE] 6: fixed labels...")
                    lab.pos <- data.frame(x=pos[jj,1], y=pos[jj,2])
                }
                rownames(lab.pos) <- rownames(pos)[jj]
            } else {
                dbg("[pgx.scatterPlotXY.BASE] 6: using user labels positions...")
                lab.pos <- lab.pos[match(rownames(df),rownames(lab.pos)),]
            }
            segments(df$x, df$y, lab.pos$x, lab.pos$y, col='#222222AA', lwd=0.85)
            text(lab.pos$x, lab.pos$y, labels=df$z, cex=0.7*df$cex)
            
        }
    }
    
    ## parameter name
    if(!is.null(title) && title!="") {
        ##legend("topleft", legend=title, cex=cex.title, bty=bty, 
        ##       inset=c(0,0), x.intersp=0, y.intersp=0, xpd=NA)
        ##legend("topleft", legend=title, cex=cex.title, bty='n', ## inset=c(-0.07,-0.07),
        ##       y.intersp=0.10, x.intersp=0.0, xpd=TRUE)
        mtext(title, 3, adj=0, padj=-0.35, cex=0.9*cex.title)
    }
    
    if(set.par) {
        ## suppressWarnings(par(par.save))
    }
    out <- list()
    out$lab.pos <- lab.pos
    invisible(out)
}

pgx.scatterPlotXY.GGPLOT <- function(pos, var=NULL, type=NULL, col=NULL, cex=NULL,
                                     cex.lab=0.8, cex.title=1.2, cex.clust=1.5, cex.legend=1,
                                     zoom=1, legend=TRUE, bty='n', hilight=NULL, 
                                     zlim=NULL, zlog=FALSE, softmax=FALSE, zsym=FALSE, 
                                     xlab = NULL, ylab=NULL, xlim=NULL, ylim=NULL,
                                     hilight2=hilight, hilight.col='black',
                                     hilight.lwd=0.8, hilight.cex=NULL,
                                     opacity=1, label.clusters=FALSE, labels=NULL,
                                     legend.ysp=0.85, legend.pos = "bottomleft",
                                     tooltip=NULL, theme=NULL, set.par=TRUE,
                                     label.type=c("text","box"), 
                                     title=NULL, nrows=NULL,  barscale=0.8 )    
{
    require(ggplot2)

    if(0) {
        cex.lab=0.8; cex.title=1.2; cex.clust=1.5; cex.legend=1;
        zoom=1; legend=TRUE; bty='n'; hilight=NULL; 
        zlim=NULL; zlog=FALSE; softmax=FALSE; zsym=FALSE; 
        xlab = NULL; ylab=NULL; xlim=NULL; ylim=NULL;
        hilight2=hilight; hilight.col='black';
        hilight.lwd=0.8; hilight.cex=NULL;
        opacity=1; label.clusters=FALSE; labels=NULL;
        legend.ysp=0.85; legend.pos = "bottomleft";
        tooltip=NULL; theme=NULL; set.par=TRUE;
        label.type="text"
    }

    if(!is.null(var) && !is.null(ncol(var))) {
        var <- var[,1]
    }
    if(!is.null(var) && length(var)==nrow(pos) && is.null(names(var)) ) {
        names(var) <- rownames(pos)
    }
    if(is.null(var)) {
        var <- rep("_",nrow(pos))
        names(var) <- rownames(pos)
    }
    var <- var[match(rownames(pos),names(var))]    
    names(var) <- rownames(pos)
    
    if(is.null(type)) {
        type <- c("numeric","factor")[1 + class(var) %in% c("factor","character")]
    }
    label.type <- label.type[1]
    if(is.null(colnames(pos))) {
        colnames(pos) <- c("x","y")
    }
    
    ## automatically set pointsize of dots
    if(is.null(cex)) {
        nr <- nrow(pos)
        i <- as.integer(cut(nr,breaks=c(0,100,500,1000,5000,Inf)))
        cex  <- c(2,1.4,1,0.7,0.4)[i]
    }
    if(is.null(hilight.cex))
        hilight.cex = cex
    
    ## normalize pos
    if(is.null(xlim)) xlim <- range(pos[,1])
    if(is.null(ylim)) ylim <- range(pos[,2])
    if(zoom!=1) {
        cx <- mean(range(pos[,1]))
        cy <- mean(range(pos[,2]))
        dx <- diff(range(pos[,1]))
        dy <- diff(range(pos[,2]))
        xlim <- cx + 0.5 * c(-1,1.05) * dx / zoom
        ylim <- cy + 0.5 * c(-1,1.05) * dy / zoom
    }
    
    ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
    )

    legend.justification <- NULL
    legend.position <- legend.pos
    if(legend.pos=="bottomright") {
        legend.justification <- c(1,0)
        legend.position <- c(1.0, 0.0)
    } else if(legend.pos=="topleft") {
        legend.justification <- c(0,1)
        legend.position <- c(.0, 1.0)
    } else if(legend.pos=="topright") {
        legend.justification <- c(1,1)
        legend.position <- c(1.0, 1.0)
    } else if(legend.pos=="bottomleft") {
        legend.justification <- c(0,0)
        legend.position <- c(0.0, 0.0)
    }

    if(is.null(xlab)) xlab <- colnames(pos)[1]
    if(is.null(ylab)) ylab <- colnames(pos)[2]
    
    plt <- NULL
    ## Plot the discrete variables
    if(type=="factor") {
        
        z1 <- var
        if(class(z1)!='factor') z1 <- factor(as.character(z1))
        nz <- length(levels(z1))
        col1 <- NULL
        if(!is.null(col)) {
            col1 <- col
        } else if(nz==2) {
            col1 <- rev(grey.colors(2, end=0.8))
            col1 <- c("#AAAAAA55","#555555FF")
            col1 <- c("#00008855","#AA0000FF")
            col1 <- c("#CCCCCC55","#AA0000FF")
        } else if(nz==1) {
            col1 <- c("#22222255")
        } else {
            col1 <- c(RColorBrewer::brewer.pal(8,"Set1"),
                      RColorBrewer::brewer.pal(8,"Set2"),
                      RColorBrewer::brewer.pal(12,"Set3"))    
        }
        col1 <- Matrix::head(rep(col1,99),nz)
        pt.col <- col1[z1]
        pt.col[is.na(pt.col)] <- "#DDDDDD11"
        if(opacity<1) {
            pt.col <- add_opacity(pt.col, opacity**0.33)
            col1 <- add_opacity(col1, opacity**0.33)
        }
        
        ##tooltip <- paste0(pgx$genes$gene_name,"<br>",pgx$genes$gene_name)
        tt <- ifelse(!is.null(title),title,"var")
        tt  <- "value"
        tooltip <- rownames(pos)
        tooltip <- paste(rownames(pos),"<br>",tt,"=",z1)
        label1 <- rownames(pos)
        if(!is.null(labels) && length(labels)==nrow(pos)) label1 <- labels

        df <- data.frame(x=pos[,1], y=pos[,2], name=rownames(pos),
                         ## variable=z1,
                         text=tooltip, label=label1)
        jj <- order(-table(pt.col)[pt.col]) ## plot less frequent points last...            
        df <- df[jj,]
        pt.col <- pt.col[jj]
        plt <- ggplot(df, aes(x,y), legend=legend) +
            geom_point( shape=20, alpha=opacity, size=1.8*cex, col=pt.col ) +
            scale_color_manual( values=col1, name=title, na.value="#DDDDDD44") 
        
        ## label cluster
        if(label.clusters) {
            
            mpos <- apply(pos,2,function(x) tapply(x,z1,median))
            ## text(med.pos, labels=rownames(med.pos),cex=1.6)
            mlab <- rownames(mpos)
            ##if(!is.null(labels)) mlab <- labels[rownames(mpos)]
            ## plt <- plt + annotate(
            ##                  geom="text",x=mpos[,1],  y=mpos[,2],
            ##                  label=mlab, size=4.5*cex.clust,
            ##                  lineheight=0.7)
            df1 <- data.frame( x=mpos[,1], y=mpos[,2], name=rownames(mpos))
            if(label.type=="text")  labelFUN <- ggrepel::geom_text_repel
            if(label.type=="box") labelFUN <- ggrepel::geom_label_repel
            plt <- plt +
                labelFUN(
                ##geom_label_repel(
                    data = df1,
                    aes(x=x, y=y, label=name),
                    size = 3.0 * cex.clust,
                    color = "black",
                    label.size = 0.10,
                    max.overlaps = 99,
                    fill = scales::alpha(c("white"),0.7),
                    segment.color = "grey30",
                    segment.size = 0.2,
                    box.padding = grid::unit(0.4, "lines"),
                    point.padding = grid::unit(0.0, "lines")
                )
            ##plt
        }
        
        
        nlev <- length(levels(z1))
        if(legend && nlev <= 10) {
            ## plt <- plt + theme(legend.position = "bottom") 
            plt <- plt +
                theme(
                    legend.title = element_blank(),
                    legend.text = element_text(size=9*cex.legend),
                    legend.key.size = grid::unit(legend.ysp*0.8*cex.legend,"lines"),
                    legend.key.height = grid::unit(legend.ysp*0.8*cex.legend,"lines"),
                    legend.key = element_rect(color="transparent", fill=scales::alpha("white",0.0)),
                    legend.justification = legend.justification,
                    legend.position = legend.position,
                    legend.background = element_rect(fill=scales::alpha("white",0.5)),
                    legend.margin = margin(0,4,4,4),
                    legend.box.just = "right",
                    legend.box.background = element_rect(color="#888888", size=0.25),
                    legend.box.margin = margin(0.8,1,1,1)
                ) +
                guides(color = guide_legend(override.aes=list(size=2.8*cex.legend)))                    
        } else {
            plt <- plt + theme(legend.position = "none") 
        }
        
    }
    
    ## Plot the continous variables    
    if(type=="numeric") {
        z <- as.numeric(var)
        z1 <- NULL
        ## if(!is.null(zlim)) {
        ##     z1 <- (z - min(zlim)) / diff(zlim)
        ##     z1 <- pmin(pmax(z1,0),1) ## clip
        ## } else {
        ##     if(zsym) {
        ##         dz <- max(abs(z),na.rm=TRUE)
        ##         z1 <- (z + dz) / (2*dz)
        ##     } else {
        ##         z1 <- (z - min(z,na.rm=TRUE)) / diff(range(z,na.rm=TRUE))
        ##     }
        ## }
        ## ## z1 is [0:1] scaled variable
        ## if(softmax) z1 <- 0.5*(tanh(4*(z1-0.5))+1)
        cpal <- rev(viridis::viridis(11))
        cpal <- rev(RColorBrewer::brewer.pal(11,"RdYlBu"))
        if(opacity<1) {
            cpal <- add_opacity(cpal, opacity**0.33)
        }
        
        ##tooltip <- paste0(pgx$genes$gene_name,"<br>",pgx$genes$gene_title)
        tt <- ifelse(!is.null(title),title,"var")
        tt  <- "value"
        tooltip <- paste(rownames(pos),"<br>",tt,"=",round(z,digits=4))
        label1 <- rownames(pos)
        if(!is.null(labels) && length(labels)==nrow(pos)) label1 <- labels
        df <- data.frame( x=pos[,1], y=pos[,2], name=rownames(pos), 
                         variable=z, text=tooltip, label=label1 )
        df <- df[order(abs(z),na.last=FALSE),] ## strongest last??

        zr <- range(z)
        if(!is.null(zlim)) zr <- zlim
        if(zsym && min(zr,na.rm=TRUE)<0 ) zr <- c(-1,1)*max(abs(zr),na.rm=TRUE)
        zz <- round(c(zr[1], zr[2]),digits=2)

        plt <- ggplot(df, aes(x, y, color=variable)) +
            geom_point(shape=20, alpha=opacity, size=1.8*cex ) +
            scale_color_gradientn(colors=cpal, breaks=zz,
                                  labels=c(zz[1],zz[2]),
                                  na.value="#DDDDDD44") +
            expand_limits(color = zr + c(-0.01,0.01)) 
        
        ## colorscale bar
        if(legend) {
            xmax <- round(max(z,na.rm=TRUE),2)
            xmin <- round(min(z,na.rm=TRUE),2)
            plt <- plt +
                guides(colour = guide_colourbar(
                           barwidth= 0.5*barscale, barheight= 2.2*barscale )) +                
                theme(legend.title = element_blank(),
                      legend.text = element_text(size=9*cex.legend),
                      legend.justification = legend.justification,
                      legend.position = legend.position,
                      ##legend.justification = c(0,0), 
                      ##legend.position = c(0.22, 0.02),
                      ##legend.background = element_rect(fill=scales::alpha("white",0.5)),
                      legend.background = element_blank(),
                      legend.key = element_blank())             
            ##legend("bottomleft", cex=0.8, ## text.width=2,
            ##       y.intersp=0.25, x.intersp=0.5, border=NA, bty=bty,
            ##       fill=rev(cpal), legend=c(xmax,rep(NA,9),xmin), )
        } else {
            plt <- plt + theme(legend.position="none") 
        }
    }    
    
    if(!is.null(hilight)) {
        sel <- which(df$name %in% hilight)
        dbg("[pgx.scatterPlotXY.GGPLOT] sel = ",sel)
        if(is.null(hilight.col)) hilight.col <- 'transparent'
        plt <- plt +
            geom_point(
                ##data = subset(df, name %in% hilight),
                data = df[sel,],
                mapping = aes(x,y),                 
                size = 2.0 * hilight.cex,
                shape = 21,
                stroke = 0.5 * hilight.lwd,
                fill = hilight.col,                  
                color = 'black'
            )

    }

    if(!is.null(hilight2)) {
                    
        if(label.type=="text")  labelFUN <- ggrepel::geom_text_repel
        if(label.type=="box") labelFUN <- ggrepel::geom_label_repel            
        ##geom_text_repel(
        plt <- plt + labelFUN(
                         data = subset(df, name %in% hilight2),
                         aes(label = label),
                         size = 3.0*cex.lab,
                         color = "black",
                         label.size = 0.08,
                         max.overlaps = 99,
                         fill = scales::alpha(c("white"),0.6),
                         segment.color = "grey30",
                         segment.size = 0.5,                    
                         box.padding = grid::unit(0.25, "lines"),
                         point.padding = grid::unit(0.2, "lines")
                     )
    }

    
    if(!is.null(theme)) {
        plt <- plt + theme
    } else {
        plt <- plt + theme_bw()
    }

    ## additional theme
    plt <- plt + 
        xlim(xlim[1], xlim[2]) +
        ylim(ylim[1], ylim[2]) + 
        xlab(xlab) +  ylab(ylab) + ggtitle(title) +
        theme(
            plot.title = element_text(size=11*cex.title, hjust=0, vjust=+0),
            axis.text.x = element_text(size=7),
            axis.text.y = element_text(size=7),            
            axis.title.x = element_text(size=9, vjust=+2),
            axis.title.y = element_text(size=9, vjust=+0),            
            panel.border = element_rect(fill=NA, color="grey20", size=0.15)
        ) 
    plt
}

pgx.scatterPlotXY.PLOTLY <- function(pos,
                                     var=NULL, type=NULL, col=NULL,
                                     cex=NULL, cex.lab=0.8, cex.title=1.2,
                                     cex.clust=1.5, cex.legend = 1,
                                     xlab = NULL, ylab = NULL, xlim=NULL, ylim=NULL,
                                     axis=TRUE, zoom=1, legend=TRUE, bty='n',
                                     hilight=NULL, hilight2=hilight, hilight.col=NULL,
                                     hilight.cex=NULL, hilight.lwd=0.8,
                                     zlim=NULL, zlog=FALSE, zsym=FALSE, softmax=FALSE, 
                                     opacity=1, label.clusters=FALSE,
                                     labels=NULL, label.type=NULL, 
                                     tooltip=NULL, theme=NULL, set.par=TRUE,
                                     title="", nrows=NULL, source=NULL, key=NULL,
                                     displayModeBar=FALSE)    
{
    
    
    

    if(0) {
        var=NULL; type=NULL; col=NULL;
        cex=NULL; cex.lab=0.8; cex.title=1.2;
        cex.clust=1.5; cex.legend = 1;
        xlab = NULL; ylab = NULL; xlim=NULL; ylim=NULL;
        axis=TRUE; zoom=1; legend=TRUE; bty='n';
        hilight=NULL; hilight2=hilight; hilight.col=NULL;
        hilight.cex=NULL; hilight.lwd=0.8;
        zlim=NULL; zlog=FALSE; zsym=FALSE; softmax=FALSE; 
        opacity=1; label.clusters=FALSE;
        labels=NULL; label.type=NULL; 
        tooltip=NULL; theme=NULL; set.par=TRUE;
        title=""; nrows=NULL; source=NULL; key=NULL;
        displayModeBar=FALSE    
    }

    if(!is.null(var) && !is.null(ncol(var))) {
        var <- var[,1]
    }
    if(!is.null(var) && length(var)==nrow(pos) && is.null(names(var)) ) {
        names(var) <- rownames(pos)
    }
    if(is.null(var)) {
        var <- rep("_",nrow(pos))
        names(var) <- rownames(pos)
    }
    var <- var[match(rownames(pos),names(var))]
    names(var) <- rownames(pos)
    
    if(is.null(type)) {
        type <- c("numeric","factor")[1 + class(var) %in% c("factor","character")]
    }
    
    ## automatically set pointsize of dots
    if(is.null(cex)) {
        nr <- nrow(pos)
        i <- as.integer(cut(nr,breaks=c(0,100,500,1000,5000,Inf)))
        cex  <- c(2,1.4,1,0.7,0.4)[i]
    }
    if(is.null(colnames(pos))) {
        colnames(pos) <- c("x","y")
    }
    if(is.null(hilight.cex)) hilight.cex = cex
    
    ## normalize pos
    ## normalize pos
    xlim0 <- range(pos[,1])
    ylim0 <- range(pos[,2])
    if(zoom!=1) {
        cx <- mean(range(pos[,1]))
        cy <- mean(range(pos[,2]))
        dx <- diff(range(pos[,1]))
        dy <- diff(range(pos[,2]))
        xlim0 <- cx + 0.5 * c(-1,1.05) * dx / zoom
        ylim0 <- cy + 0.5 * c(-1,1.05) * dy / zoom
    }
    
    ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
    )
    if(is.null(xlab)) xlab <- colnames(pos)[1]
    if(is.null(ylab)) ylab <- colnames(pos)[2]

    ##tooltip <- paste0(pgx$genes$gene_name,"<br>",pgx$genes$gene_name)
    ##gg.title <- pgx$genes[rownames(pos),"gene_title"]
    ##tooltip <- paste(rownames(pos),"<br>", title,"=",z1,"<br>",gg.title)
    if(type=='numeric') var <- round(var, digits=4)
    tooltip1=NULL
    tooltip1 <- paste0(
        rownames(pos),
        "<br>value = ",var,
        "<br>x = ",round(pos[,1],2),"; y = ",round(pos[,2],2)
    )
    if(!is.null(tooltip)) {
        tooltip1 <- paste0(tooltip1,"<br>",tooltip)
    }
    label1 <- rownames(pos)
    if(!is.null(labels) && labels[1]!=FALSE) label1 <- labels
    label1 <- gsub('[)(]','',label1)  ## HTML labels do not like it...
    label1 <- gsub('[_]',' ',label1)  ## HTML labels do not like it...
    label1 <- sapply(label1, function(s) paste(strwrap(s,50),collapse='<br>'))
        
    plt <- NULL
    ## Plot the discrete variables
    if(type=="factor") {
        
        z1 <- factor(var)
        nz <- length(levels(z1))
        cpal <- NULL
        if(!is.null(col)) {
            cpal <- col
        } else if(nz==2) {
            cpal <- rev(grey.colors(2, end=0.8))
            cpal <- c("#AAAAAA55","#555555FF")
            cpal <- c("#00008855","#AA0000FF")
            cpal <- c("#CCCCCC55","#AA0000FF")
        } else if(nz==1) {
            cpal <- c("#22222255")
        } else {
            cpal <- c(RColorBrewer::brewer.pal(8,"Set1"),
                      RColorBrewer::brewer.pal(8,"Set2"),
                      RColorBrewer::brewer.pal(12,"Set3"))    
        }
        cpal <- Matrix::head(rep(cpal,99),nz)

        ##pt.col <- cpal[z1]
        ##pt.col[is.na(pt.col)] <- "#DDDDDD55"
        if(opacity<1) {
            ##pt.col <- add_opacity(pt.col, opacity)
            cpal <- add_opacity(cpal, opacity**0.33)
        }
        
        ## prepare dataframe
        df <- data.frame(
            x = pos[,1], y = pos[,2], name = rownames(pos),
            value = z1, text = tooltip1, label = label1)

        ## plot less frequent points first... (NEED RETHINK)           
        jj <- order(-table(z1)[z1]) 
        df <- df[jj,]
              
    }
    
    ## Plot the continous variables    
    if(type=="numeric") {
        z <- as.numeric(var)
        z1 <- NULL
        ## if(!is.null(zlim)) {
        ##     z1 <- (z - min(zlim)) / diff(zlim)
        ##     z1 <- pmin(pmax(z1,0),1) ## clip
        ## } else {
        ##     z1 <- (z - min(z,na.rm=TRUE)) / diff(range(z,na.rm=TRUE))
        ## }
        ## if(softmax) z1 <- 0.5*(tanh(4*(z1-0.5))+1)
        cpal <- rev(viridis::viridis(11))
        cpal <- rev(RColorBrewer::brewer.pal(11,"RdYlBu"))
        if(opacity<1) {
            cpal <- add_opacity(cpal, opacity**0.33)
        }

        ## df <- data.frame( x=pos[,1], y=pos[,2], variable=z, text=tooltip1 )
        df <- data.frame(x=pos[,1], y=pos[,2], name=rownames(pos),
                         value=z, text=tooltip1, label=label1)
        
        ## plot low values first
        df <- df[order(abs(z)),,drop=FALSE]

        cmin = min(df$value, na.rm=TRUE)
        cmax = max(df$value, na.rm=TRUE)
        if(zsym) {
            cmax = max(abs(df$value),na.rm=TRUE)
            cmin = -cmax
        }
    }

    ## remove NA??
    ##df <- df[!is.na(df$value),]
    
    ##---------------- call PLOTLY -----------
    if(is.null(source)) source <- paste0(sample(LETTERS,10),collapse='')
    ##plt <- plotly::plot_ly(df,
    plt <- plotly::plot_ly(
        ##data = df[jj,,drop=FALSE],
        source = source,
        ##key = rownames(df),
        showlegend = FALSE                    
    ) 

    any(is.na(df$value))
    if(any(is.na(df$value))) {
        jj <- which(is.na(df$value))
        plt <- plt %>%            
            plotly::add_markers(
                data = df[jj,,drop=FALSE],
                x = ~x, y = ~y,
                ##colors = ~value,
                colors = cpal,
                text = ~text, hoverinfo='text', 
                marker = list(size=7*cex, opacity=opacity, color='#DDDDDD44'),            
                showlegend = FALSE,
                key = ~label,
                mode = "markers",
                type = "scattergl"
            )        
    }

    jj <- which(!is.na(df$value))
    plt <- plt %>%
        plotly::add_markers(
            data = df[jj,,drop=FALSE],
            x = ~x, y = ~y,
            color = ~value, colors = cpal,
            text = ~text, hoverinfo='text', 
            marker = list(size=7*cex, opacity=opacity),            
            showlegend = FALSE,
            key = ~label,
            mode = "markers",
            type = "scattergl")   

    if(!is.null(hilight)) {

        jj <- which(rownames(df) %in% hilight)
        col1 = 'transparent'
        if(!is.null(hilight.col)) col1 <- hilight.col
        df1 <- df[jj,]
        plt <- plt %>%
            ##add_trace(
            plotly::add_markers(            
                data = df[jj,],
                x = ~x,  y = ~y,
                color = ~value, colors = cpal,
                color = NULL, 
                key = ~name,
                mode = 'markers', type = 'scattergl', ## color=NULL,
                text = ~text, hoverinfo='text',
                ## showlegend = FALSE,
                marker = list(
                    color = col1,
                    size = 5*hilight.cex,
                    showlegend = FALSE,                    
                    showscale = FALSE,                    
                    line = list(
                        color = '#000000',
                        width = 1.0*hilight.lwd)
                )
            ) 
        
    }

    if(!is.null(hilight2)) {
        jj <- which(rownames(df) %in% hilight2)            
        message("[pgx.scatterPlotXY.PLOTLY] hilight2 points : len.jj = ",length(jj))
        plt <- plt %>%            
                plotly::add_annotations( 
                    data = df[jj,,drop=FALSE],
                    x = ~x,  y = ~y, text = ~label,                   
                    ##textposition = "top center",
                    yanchor = "bottom",
                    xanchor = "center",  ## left,center,right
                    showarrow = FALSE,
                    showlegend = FALSE,
                    font = list(size = 12*cex.lab),
                    xref = "x", yref = "y")
    }
    ## label cluster
    if(label.clusters) {
        mpos <- apply(pos,2,function(x) tapply(x,z1,median))
        mlab <- rownames(mpos)
        ##if(!is.null(labels)) mlab <- labels[rownames(mpos)]
        plt <- plt %>%
            plotly::add_annotations(
                x = mpos[,1],  y = mpos[,2], 
                showarrow = FALSE, text = mlab,
                font = list(size = 15*cex.clust),
                xref = "x", yref = "y")
    }

    if(!legend) {
        plt <- plt %>%
            plotly::hide_colorbar() %>%
            plotly::hide_legend() %>%
            plotly::layout(showlegend=FALSE)
    }

    if(legend && type=='numeric') {
        plt <- plt %>%
            plotly::colorbar(limits = c(cmin, cmax), len=0.15, thickness=9,
                     x=0.01, y=0.1, title='', tickfont = list(size=9) )             
    }
    
    ## add legend and title
    plt <- plt %>%
        plotly::layout(
            showlegend = legend,
            xaxis = list(title=xlab, titlefont=list(size=12)),
            yaxis = list(title=ylab, titlefont=list(size=12)),
            margin = list(l=5, r=5, b=25, t=25, pad=3)
        )

    if(!is.null(title) && title!="") {
        plt <- plt %>%
            plotly::layout(
                annotations = list(text=title, font = list(size=14*cex.title),
                                   xref="paper", yref="paper",
                                   yanchor = "bottom", xanchor = "left",
                                   align = "right", x=0, y=1 , showarrow = FALSE )
            )
    }
    
    if(axis==FALSE) {
        plt <- plt %>% plotly::layout(xaxis = ax, yaxis = ax)
    }
    plt
}

pgx.scatterPlotXY.D3 <- function(pos, var=NULL, type=NULL, col=NULL, cex=1,
                                  cex.lab=0.8, cex.title=1.2, cex.clust=1.5, cex.legend=1,
                                  zoom=1, legend=TRUE, bty='n', hilight=NULL, 
                                  zlim=NULL, zlog=FALSE, softmax=FALSE,
                                  xlab = NULL, ylab=NULL, hilight2=hilight,
                                  opacity=1, label.clusters=FALSE, labels=NULL,
                                  legend.ysp=0.85, legend.pos = "bottomleft",
                                  tooltip=NULL, theme=NULL, set.par=TRUE,
                                  title=NULL, nrows=NULL,  barscale=0.8 )    
{
    
    
    

    if(0) {
        type=NULL;col=NULL;cex=NULL
        cex.lab=0.8;cex.title=1.2;cex.clust=1.5;cex.legend=1
        zoom=1;legend=TRUE;bty='n';hilight=NULL
        zlim=NULL;zlog=FALSE;softmax=FALSE
        xlab = NULL;ylab=NULL
        opacity=1;label.clusters=FALSE;labels=NULL
        legend.ysp=0.85;legend.pos = "bottomleft"
        title=NULL;nrows=NULL; barscale=0.8
    }
    if(is.null(colnames(pos))) {
        colnames(pos) <- c("x","y")
    }
    if(is.null(xlab))
        xlab <- colnames(pos)[1]
    if(is.null(ylab))
        ylab <- colnames(pos)[2]

    df <- data.frame(x=pos[,1], y=pos[,2], z=var, names=rownames(pos))
    Matrix::head(df)
    if(!is.null(var)) {
        plt <- scatterD3::scatterD3(
            data = df, x = x, y = y, ## lab = names,
            col_var = z, ## symbol_var = am,
            xlab = xlab, ylab = ylab,
            point_size = 32 * cex,
            ##symbol_lab = "Manual transmission",
            legend_width = 70,            
            col_lab = "value")
    } else {
        plt <- scatterD3::scatterD3(data = df, x = x, y = y, ## lab = names,
                         xlab = xlab, ylab = ylab,
                         point_size = 32 * cex,
                         ##symbol_lab = "Manual transmission",
                         col_lab = "value")
    }
    plt
}

file="tmp.pdf";width=height=8
plotWidget.PLEASECHECK <- function(plt,file,width=8,height=8) {
    HTMLFILE <- paste0(tempfile(),"_plotwidget.html")
    HTMLFILE
    htmlwidgets::saveWidget(plt, HTMLFILE)
    res <- 1
    if(grepl("png$",file)) res=100
    webshot::webshot(HTMLFILE, file=file, vwidth=width*res,vheight=height*res)
}

pgx.plotSampleClustering <- function(x, dim=2, 
                                     method=c("tsne","umap","pca"),
                                     ntop=1000, ...)
{
    method = method[1]    
    clust <- pgx.clusterMatrix(
        x, is.logx=TRUE, perplexity=NULL,
        ntop=ntop, sv.rank=-1, dims=dim,
        row.center=TRUE, row.scale=FALSE,
        find.clusters=FALSE, kclust=1,
        prefix="C", clust.detect = "louvain",
        method = method )

    ##col1 <- as.integer(clust$idx)
    plot( clust$pos2d, ... )

}

pgx.stackedBarplot <- function(x, hz=FALSE, srt=NULL, cex.text=0.9, ...)
{
    ##x <- x[order(rowMeans(x,na.rm=TRUE)),]    
    ##barplot( t(x), beside=FALSE, las=3)
    x.pos <- pmax(x,0)
    x.neg <- pmin(x,0)
    y0 <- max(abs(rowSums(x,na.rm=TRUE)))
    y0 <- max(rowSums(pmax(x,0),na.rm=TRUE),
              rowSums(pmax(-x,0),na.rm=TRUE))

    rownames(x.neg) <- NULL
    p <- NULL
    if(hz==TRUE) {
        ##p <- barplot( t(x.pos), horiz=TRUE, beside=FALSE, las=1, xlim=c(-1,1)*y0 )
        ##barplot( t(x.neg), horiz=TRUE, beside=FALSE, las=1, add=TRUE )
        
        p <- barplot( t(x.pos), horiz=TRUE, beside=FALSE, las=1, xlim=c(-1,1)*y0, ... )
        barplot( t(x.neg), horiz=TRUE, beside=FALSE, las=1, add=TRUE, ... )


    } else {
        p <- barplot( t(x.pos), beside=FALSE, las=3, ylim=c(-1.1,1.1)*y0, ... )
        barplot( t(x.neg), beside=FALSE, las=3, add=TRUE, ... )

        if(!is.null(srt)) {
            text(p, par("usr")[3], labels=rownames(x), srt=srt, adj=1, xpd=TRUE, cex=cex.text)
        }

    }
        
}

## for plotly
darkmode <- function(p, dim=2) {
    font.par <- list(
        color = "#AAA"
    )
    axis.par <- list(
        color='#AAA',
        linecolor='#AAA',
        tickcolor='#AAA',
        zerolinecolor='#AAA'
    )
    p <- plotly::layout(p,
                plot_bgcolor = "rgb(10,20,40)",
                paper_bgcolor = "rgb(10,20,40)",
                xaxis = axis.par,
                yaxis = axis.par,
                ## zaxis = axis.par,
                title = font.par
                )
    if(dim==3) {
        p <- plotly::layout(p,
                    zaxis = axis.par
                    )
    }
    return(p)
}

myplot_ly <- function(..., theme="default") {
    ## 'Themed' plotly 
    ##
    ##
    if(theme=="default") {
        p <- plotly::plot_ly(...)

    } else if(theme=="dark") {
        font.par <- list(
            color = "#FFF"
        )
        axis.par <- list(
            color='#FFF',
            linecolor='#FFF'
        )
        p <- plotly::plot_ly(...) %>%
            plotly::layout(
                plot_bgcolor = "rgb(10,30,50)",
                paper_bgcolor = "rgb(10,30,50)",
                xaxis = axis.par,
                yaxis = axis.par,
                title = font.par
        )
    }
    return(p)
}


##lfc=1;psig=0.05;showlegend=FALSE;xlab=ylab="";group.names=c("group1","group2")
##showlegend=TRUE;highlight=NULL;marker.size=5;label=NULL;marker.type="scatter";displayModeBar=TRUE
plotlyMA <- function(x, y, names, source="plot1",
                     group.names=c("group1","group2"),
                     xlab = "average expression (log2.CPM)",
                     ylab = "effect size (log2.FC)",
                     lfc=1, psig=0.05, showlegend=TRUE, highlight=NULL,
                     marker.size = 5, label = NULL, label.cex = 1,
                     marker.type="scatter", displayModeBar=TRUE )
{

    

    if(is.null(highlight)) highlight <- names
    i0 <- which(!names %in% highlight)
    i1 <- which(names %in% highlight)
    
    p <- plotly::plot_ly(
        type = marker.type, mode = 'markers'
        ##type='scattergl', mode='markers',
        ##source=source, key=1:length(x)
    )
    
    p <- p %>%
        plotly::event_register('plotly_hover') %>%
        plotly::event_register('plotly_selected')

    if(length(i0)) {
        p <- p %>%
            plotly::add_trace(
                x = x[i0],
                y = y[i0],
                text = names[i0],
                marker = list(
                    size = marker.size,
                    color = '#ccc'
                ),
                showlegend = showlegend
            )
    }
    
    if(length(i1)) {
        p <- p %>%
            plotly::add_trace(
                x = x[i1],
                y = y[i1],
                text = names[i1],
                marker = list(
                    size = marker.size,
                    color = '#1f77b4'
                ),
                showlegend = showlegend
            )
    }
    
    if(!is.null(label) && length(label)>0) {
        i2 = which(names %in% label)
        p <- p %>%
            plotly::add_annotations(
                x = x[i2],
                y = y[i2],
                text = names[i2],
                font = list(
                    size = 12 * label.cex,
                    color = '#1f77b4'
                ),
                showarrow = FALSE,
                yanchor = "bottom",
                yshift = 2,
                textposition = 'top'
            )        
    }

    x1 = 1.05*max(x)
    yy = 1.05*max(abs(y))
    abline1 = list(type='line', y0= -lfc, y1= -lfc, x0=0, x1=x1,
                   line=list(dash='dot', width=1, color="grey"))
    abline2 = list(type='line', y0= +lfc, y1= +lfc, x0=0, x1=x1,
                   line=list(dash='dot', width=1, color="grey"))
        
    xrange <- c(0,1)*max(abs(x))*1.05
    xrange <- range(x,na.rm=TRUE)
    yrange <- c(-1,1)*max(abs(y),na.rm=TRUE)*1.05        
    xaxis = list( title = xlab, range = xrange, gridwidth=0.2, showgrid=FALSE )
    yaxis = list( title = ylab, range = yrange, gridwidth=0.2, showgrid=FALSE )    

    p <- p %>%
        plotly::layout(
            shapes = list(abline1,abline2),
            xaxis = xaxis, yaxis = yaxis,
            showlegend = showlegend,
            hovermode='closest',
            dragmode= 'select') %>%
        plotly::config(displayModeBar = displayModeBar)
    
    xann <- c(1,1)*0.95
    yann <- c(0,0.99)
    ann.text <- paste("UP in", group.names[c(2,1)])
    
    p <- p %>%
        plotly::add_annotations(
            x = xann,
            y = yann,
            text = ann.text,
            font = list(size=10),
            xanchor = c('right','right'),
            align = c('right','right'),
            showarrow = FALSE,
            xref = 'paper',
            yref = 'paper',
            borderpad = 3, 
            bordercolor = 'black',
            borderwidth = 0.6)

    p
}

##lfc=1;psig=0.05;showlegend=FALSE;xlab=ylab="";group.names=c("group1","group2");highlight=NULL
plotlyVolcano <- function(x, y, names, source="plot1", group.names=c("group1","group2"),
                          xlab="effect size (logFC)", ylab="significance (-log10p)",
                          lfc=1, psig=0.05, showlegend=TRUE, highlight=NULL,
                          marker.size = 5, label=NULL, label.cex = 1,
                          marker.type='scatter', displayModeBar=TRUE )
{

    

    if(is.null(highlight)) highlight <- names
    i0 <- which(!names %in% highlight)
    i1 <- which(names %in% highlight)
    
    p <- plotly::plot_ly(
        type = marker.type, mode='markers'
        ##source=source, key=1:length(x)
    )
    
    p <- p %>%
        plotly::event_register('plotly_hover') %>%
        plotly::event_register('plotly_selected')

    if(length(i0)) {
        p <- p %>%
            plotly::add_trace(
                x = x[i0],
                y = y[i0],
                text = names[i0],
                marker = list(
                    size = marker.size,
                    color = '#ccc'
                ),
                showlegend = showlegend
            )
    }
    
    if(length(i1)) {
        p <- p %>%
            plotly::add_trace(
                x = x[i1],
                y = y[i1],
                text = names[i1],
                marker = list(
                    size = 5,
                    color = '#1f77b4'
                ),
                showlegend = showlegend
            )
    }
    
    if(!is.null(label) && length(label)>0) {
        i2 = which(names %in% label)
        p <- p %>%
            plotly::add_annotations(
                x = x[i2],
                y = y[i2],
                text = names[i2],
                font = list(
                    size = 12*label.cex,
                    color = '#1f77b4'
                ),
                showarrow = FALSE,
                yanchor = "bottom",
                yshift = 2,
                textposition = 'top'
            )
        
    }

    y0 = -log10(psig)
    y1 = 1.05*max(y)
    xx = 1.05*max(abs(x))
    abline1 = list(type='line', x0= -lfc, x1= -lfc, y0=0, y1=y1,
                   line=list(dash='dot', width=1, color="grey"))
    abline2 = list(type='line', x0= +lfc, x1= +lfc, y0=0, y1=y1,
                   line=list(dash='dot', width=1, color="grey"))
    abline3 = list(type='line', x0= -xx, x1= +xx, y0=y0, y1=y0,
                   line=list(dash='dot', width=1, color="grey"))
    
    max.absx <- max(max(abs(x),na.rm=TRUE),lfc*1.2)
    max.absy <- max(max(abs(y),na.rm=TRUE),y0*1.2)
    xrange <- c(-1,1)*max.absx*1.05
    if(min(x)>=0) xrange <- c(0,1)*max.absx*1.05
    yrange <- c(0,1)*max.absy*1.05
    ## if(min(y)>=0) yrange <- c(0,1)*max.absy*1.05
        
    xaxis = list( title = xlab, range = xrange, gridwidth=0.2, showgrid=FALSE )
    yaxis  = list( title = ylab, range = yrange, gridwidth=0.2, showgrid=FALSE )    
    p <- p %>%
        plotly::layout(
            shapes = list(abline1,abline2,abline3),
            xaxis = xaxis, yaxis = yaxis,
            showlegend = showlegend,
            hovermode='closest',
            dragmode= 'select') %>%
        plotly::config(displayModeBar = displayModeBar)
    
    xann <- c(0.01,0.99)
    yann <- c(1,1)*1.04
    ann.text <- paste("UP in", group.names[c(2,1)])
    
    p <- p %>%
        plotly::add_annotations(
            x = xann,
            y = yann,
            text = ann.text,
            font = list(size = 10),
            xanchor = c('left','right'),
            align = c('left','right'),
            showarrow = FALSE,
            xref = 'paper',
            yref = 'paper',
            borderpad = 3, 
            bordercolor = 'black',
            borderwidth = 0.6)

    p
}

corclust <- function(x) {
    dd <- as.dist(1 - stats::cor(t(x),use="pairwise"))
    hc <- fastcluster::hclust(dd, method="ward.D2" )
    hc
}

## Override add_col_annotation to be able to suppress titles
##
##
require(iheatmapr)
setMethod(add_col_annotation,
          c(p = "Iheatmap"),
          function(p,
                   annotation,
                   colors = NULL,
                   side = c("top","bottom"),
                   size = 0.05,
                   buffer = 0.015,
                   inner_buffer = buffer / 2,
                   show_title = TRUE,
                   layout = list()){
            
            side <- match.arg(side)
            # Convert to data.frame
            x <- as.data.frame(annotation)
            
            for (i in seq_len(ncol(x))){
              if (is.character(x[,i]) || is.factor(x[,i]) || is.logical(x[,i])){
                if (!is.null(colors) && colnames(x)[i] %in% names(colors)){
                  tmp_colors <- colors[[colnames(x)[i]]]
                } else{
                  tmp_colors <- iheatmapr:::pick_discrete_colors(as.factor(x[,i]), p)
                }
                p <- iheatmapr::add_col_groups(p, 
                                    x[,i],
                                    name = colnames(x)[i],
                                    title = colnames(x)[i],
                                    colors = tmp_colors,
                                    side = side,
                                    size = size,
                                    buffer = if (i == 1)
                                      buffer else inner_buffer,
                                    layout = layout,
                                    show_title = show_title)
              } else if (is.numeric(x[,i])){
                if (!is.null(colors) && colnames(x)[i] %in% names(colors)){
                  tmp_colors <- colors[[colnames(x)[i]]]
                } else{
                  tmp_colors <- pick_continuous_colors(zmid = 0, 
                                                       zmin = min(x[,i]),
                                                       zmax = max(x[,i]), p)
                }
                p <- iheatmapr::add_col_signal(p, 
                                    x[,i],
                                    name = colnames(x)[i],
                                    colors = tmp_colors,
                                    side = side,
                                    size = size,
                                    buffer = if (i == 1)
                                      buffer else inner_buffer,
                                    layout = layout,
                                    show_title = show_title)
              } else{
                stop("Input should be character, factor, logical, or numeric")
              }
            }
            return(p)
          })


##X=head(ngs$X,100);annot=ngs$samples
##row_annot_width=0.03;colors=NULL;label_size=11;scale="row.center"
##xtips=NULL;ytips=NULL;lmar=60

pgx.splitHeatmapFromMatrix <- function(X, annot, idx=NULL, splitx=NULL, 
                                       xtips=NULL, ytips=NULL, row_clust=TRUE,
                                       row_annot_width=0.03, scale="row.center",
                                       colors=NULL, lmar=60,
                                       rowcex=1, colcex=1)
{
    
    ## constants
    col_annot_height = 0.021
    if(!is.null(idx)) idx = as.character(idx)
    if(!is.null(splitx)) splitx = as.character(splitx)
    
    ## --------- defaults
    if(is.null(xtips)) xtips = colnames(X)
    if(is.null(ytips)) ytips = rownames(X)
    if(is.null(names(xtips))) names(xtips) = colnames(X)
    if(is.null(names(ytips))) names(ytips) = rownames(X) 
    
    ## --------- scaling
    if("row.center" %in% scale) X <- X - rowMeans(X, na.rm=TRUE)
    if("row" %in% scale) X <- t(scale(t(X)))
    if("col.center" %in% scale) X <- t(t(X) - colMeans(X, na.rm=TRUE))
    if("col" %in% scale) X <- scale(X)
    if("row.bmc" %in% scale) {
        if(is.null(splitx)) {
            X <- X - rowMeans(X, na.rm=TRUE)
        } else {
            for(ux in unique(splitx)) {
                jj <- which(splitx==ux)
                X[,jj] = X[,jj,drop=FALSE] - rowMeans(X[,jj,drop=FALSE],na.rm=TRUE)
            }
        }        
    }
    
    ## ------ split Y-axis (genes) by factor
    hc.order <- function(x) {
        suppressWarnings( dd <- as.dist(1 - stats::cor(t(x),use="pairwise")) )
        ## cat("DBG >pgx-plotting:pgx.splitHeatmapX> sum.is.na.dd=",sum(is.na(dd)),"\n")
        if(sum(is.na(dd))) dd[is.na(dd)] <- 1
        hc <- fastcluster::hclust(dd, method="ward.D2" )
        rownames(x)[hc$order]
    }
    if(!is.null(idx)) {        
        if(row_clust) {
            kk  <- tapply(rownames(X),idx,function(k) c(hc.order(X[k,]),"   "))
        } else {
            kk  <- tapply(rownames(X),idx,function(k) c(k,"   "))
        }
        idx <- tapply(idx,idx,function(k) c(k,NA))
        idx = as.vector(unlist(idx))
        kk <- unlist(kk)
        kk  <- kk[1:(length(kk)-1)] ## remove trailing spacer
        idx <- idx[1:(length(idx)-1)]
        X <- rbind(X,"   " = 0)[kk,]

        ## invert
        X <- X[nrow(X):1,]
        idx <- rev(idx)
        
    } else {
        if(row_clust) {
            kk <- hc.order(X[,])
            X <- X[kk,]
        }
    }
    
    ## ------ split X-axis by some group factor
    if(!is.null(splitx)) {
        xx <- tapply( colnames(X), splitx, function(i) X[,i,drop=FALSE] )
    } else {
        xx <- list("Samples"=X)
    }
    length(xx)
    
    ## ------- set colors
    colors0 = rep("Set2",ncol(annot))    
    names(colors0) = colnames(annot)
    
    if(!is.null(colors) && any(names(colors) %in% names(colors0))) {
        for(v in intersect(names(colors), names(colors0))) colors0[[v]] <- colors[[v]]
    }

    ## ------- annot need to be factor
    annotF <- data.frame(as.list(annot),stringsAsFactors=TRUE)
    rownames(annotF) = rownames(annot)
    
    grid_params <- iheatmapr::setup_colorbar_grid(
        nrows = 5, 
        y_length = 0.16, 
        y_spacing = 0.18, 
        x_spacing = 0.2,
        x_start = 1.1, 
        y_start = 0.85)

    ## maximize plot area
    mar <- list(l = 50, r = 50, b = 100, t = 100, pad = 4)
    mar <- list(l = lmar, r = 0, b = 5, t = 0, pad = 3)

    ex <- ncol(X)/ncol(xx[[1]])  ## expansion factor
    hc <- fastcluster::hclust(as.dist(1 - stats::cor(xx[[1]],use="pairwise")))
    dd <- 0.08 * ex ## left dendrogram width

    hovertemplate = "Row: %{y} <br>Column: %{x}<br>Value: %{z}<extra> </extra>"
    tooltip = iheatmapr::setup_tooltip_options(
        prepend_row = "Row: ", prepend_col = "Column: ", 
        prepend_value = "Value: ") 

    x1 <- xx[[1]]
    ##x1 <- x1[nrow(x1):1,]
    plt <- iheatmapr::main_heatmap(
        x1, name = "expression",
        colorbar_grid = grid_params,
        x = xtips[colnames(x1)],
        y = ytips[rownames(x1)],
        tooltip = tooltip,
        ## hovertemplate = hovertemplate,
        layout = list(margin = mar) ) %>%        
        ##add_col_clustering() %>%
        ##add_row_clustering() %>%
        ##add_row_clustering(method="groups", groups=idx) %>%
        iheatmapr::add_col_dendro(hc, size = 0.06) %>%
        ##add_row_dendro(hr, size = dd) %>%
        iheatmapr::add_row_title("Genes") %>%
        iheatmapr::add_col_title(names(xx)[1], side="top") %>%
        add_col_annotation(
            size = col_annot_height, buffer = 0.005, side="bottom",
            colors = colors0, show_title=TRUE,
            annotF[colnames(xx[[1]]),,drop=FALSE])
    ##plt
    length(xx)
    dim(X)
    if(ncol(X)<100 && colcex>0) {
        plt <- plt %>% iheatmapr::add_col_labels(
                           side="bottom", size=0.15*ex,
                           font=list(size=11*colcex))
    }
    
    if(length(xx)>1) {
        
        sizes = sapply(xx, ncol) / ncol(xx[[1]])
        sizes
        i=5
        i=2
        for(i in 2:length(xx)) {
            
            x1 <- xx[[i]]
            if(ncol(x1)==1) x1 <- cbind(x1,x1)
            ##x1 <- x1[nrow(x1):1,]
            hc <- fastcluster::hclust(as.dist(1 - stats::cor(x1,use="pairwise")))

            plt <- plt %>%
                iheatmapr::add_main_heatmap(
                    x1, , name = "expression",
                    x = xtips[colnames(x1)],
                    y = ytips[rownames(x1)],
                    tooltip = tooltip,
                    ## hovertemplate = hovertemplate,
                    size=sizes[i], buffer = 0.007*ex) %>%
                ##add_col_clustering() %>%
                iheatmapr::add_col_dendro(hc, size = 0.06) %>%
                iheatmapr::add_col_title(names(xx)[i], side="top") %>%
                add_col_annotation(
                    size = col_annot_height, buffer = 0.005, side="bottom",
                    colors = colors0, show_title=FALSE,
                    data.frame(annotF[colnames(x1),,drop=FALSE]))
            
            if(ncol(X)<100 && colcex>0) {
                plt <- plt %>%
                    iheatmapr::add_col_labels(side="bottom",
                                   size=0.15*ex,
                                   font=list(size=11*colcex))
            }
        }

    }

    ## ----------- row annotation (i.e. gene groups)
    if(!is.null(idx) && !is.null(row_annot_width) && row_annot_width>0 ) {
        plt <- iheatmapr::add_row_annotation(
            plt, data.frame("gene.group"=idx),
            size = row_annot_width*ex)
    }

    ## ----------- add gene/geneset names
    if(rowcex > 0) {
        
        gnames <- rownames(X)
        gnames <- gsub("[&].*[;]","",gnames) ## HTML special garbage...
        gnames <- gsub("^.*:","",gnames) ## strip prefix
        gnames <- shortstring(gnames,25)  ## shorten
        gnames <- sub("   ","-",gnames)  ## necessary otherwise strange error
        ##gnames <- sub("   ","-------------",gnames)  ## necessary otherwise strange error
        
        maxlen <- max(sapply(gnames,nchar))
        w = ifelse(maxlen >= 20, 0.45, 0.20)
        s1 <- ifelse(maxlen >= 20, 9, 11)*rowcex
        
        plt <- iheatmapr::add_row_labels(
            plt, side="right", ticktext = gnames,
            size=w*ex, font=list(size=s1) ) 
    }

    
    return(plt)
}
