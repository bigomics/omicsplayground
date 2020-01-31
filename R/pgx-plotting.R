########################################################################
## Plotting functions
########################################################################


##df=ngs$samples
pgx.testPhenoCorrelation <- function(df, plot=TRUE, cex=1)
{

    require(corrplot)
    
    cl <- sapply(df,class)
    cvar <- which(cl %in% c("numeric","integer"))
    dvar <- which(cl %in% c("factor","character"))    
    dc <- df[,cvar,drop=FALSE]
    dd <- df[,dvar,drop=FALSE]

    dim(dc)
    dim(dd)
    
    ## discrete vs discreate -> Fisher test
    fisher.P <- NULL
    if(ncol(dd)) {
        fisher.P <- matrix(NA,ncol(dd),ncol(dd))
        i=1;j=2
        for(i in 1:(ncol(dd)-1)) {
            for(j in (i+1):ncol(dd)) {
                tb <- table(dd[,i], dd[,j])
                fisher.P[i,j] <- fisher.test(tb, simulate.p.value=TRUE)$p.value
            }
        }
        rownames(fisher.P) <- colnames(dd)
        colnames(fisher.P) <- colnames(dd)
    }
        
    ## discrete vs continuous -> ANOVA or Kruskal-Wallace
    kruskal.P <- NULL
    if(ncol(dc)>0) {
        kruskal.P <- matrix(NA,ncol(dd),ncol(dc))
        i=1;j=2
        for(i in 1:ncol(dd)) {
            for(j in 1:ncol(dc)) {
                kruskal.P[i,j] <- kruskal.test(dc[,j], dd[,i])$p.value
            }
        }
        rownames(kruskal.P) <- colnames(dd)
        colnames(kruskal.P) <- colnames(dc)
    }
    
    ## continuous vs continuous -> correlation test
    cor.P <- NULL
    if(ncol(dc)>1) {
        cor.P <- matrix(NA,ncol(dc),ncol(dc))
        i=1;j=2
        for(i in 1:(ncol(dc)-1)) {
            for(j in (i+1):ncol(dc)) {
                cor.P[i,j] <- cor.test(dc[,i], dc[,j])$p.value
            }
        }
        rownames(cor.P) <- colnames(dc)
        colnames(cor.P) <- colnames(dc)
    }
    
    P <- matrix(NA,ncol(df),ncol(df))
    rownames(P) <- colnames(P) <- colnames(df)

    if(!is.null(fisher.P)) {
        ii <- match(rownames(fisher.P),rownames(P))
        jj <- match(colnames(fisher.P),colnames(P))
        P[ii,jj] <- fisher.P
    }

    if(!is.null(kruskal.P)) {
        ii <- match(rownames(kruskal.P),rownames(P))
        jj <- match(colnames(kruskal.P),colnames(P))
        P[ii,jj] <- kruskal.P
    }

    if(!is.null(cor.P)) {
        ii <- match(rownames(cor.P),rownames(P))
        jj <- match(colnames(cor.P),colnames(P))
        P[ii,jj] <- cor.P
    }
    
    ij <- which(!is.na(P),arr.ind=TRUE)
    qv <- p.adjust(P[ij], method="BH")
    Q <- P
    Q[ij] <- qv
    
    P[is.na(P)] <- 0
    P <- (P + t(P))/2
    Q[is.na(Q)] <- 0
    Q <- (Q + t(Q))/2
    
    if(plot==TRUE) {

        require(corrplot)    
        logP <- -log10(P+1e-8)
        logQ <- -log10(Q+1e-8)
        diag(logQ) <- 0
        ##par(oma=c(0,0,0,1))
        corrplot( logQ, is.corr=FALSE, type="upper",
                 mar = c(0,0,0,2),
                 p.mat = Q, sig.level = 0.05, ##insig = "blank",
                 tl.cex = cex, tl.col="black", tl.offset = 1,
                 cl.align.text = "l", cl.offset = 0.25, cl.cex = 0.7, 
                 pch.col = "grey50",
                 order="hclust")

    }

    return(list(P=P, Q=Q))
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
    p <- layout(p,
                plot_bgcolor = "rgb(10,20,40)",
                paper_bgcolor = "rgb(10,20,40)",
                xaxis = axis.par,
                yaxis = axis.par,
                ## zaxis = axis.par,
                title = font.par
                )
    if(dim==3) {
        p <- layout(p,
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
            layout(
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
plotlyMA <- function(x, y, names, source="plot1",
                     group.names=c("group1","group2"),
                     xlab = "average expression (log2.CPM)",
                     ylab = "effect size (log2.FC)",
                     lfc=1, psig=0.05, showlegend=TRUE, highlight=NULL,
                     marker.size = 5, label=NULL, displayModeBar=TRUE )
{

    require(plotly)

    if(is.null(highlight)) highlight <- names
    i0 <- which(!names %in% highlight)
    i1 <- which(names %in% highlight)
    
    p <- plot_ly(
        type='scatter', mode='markers'
        ##type='scattergl', mode='markers',
        ##source=source, key=1:length(x)
    )
    
    p <- p %>%
        event_register('plotly_hover') %>%
        event_register('plotly_selected')

    if(length(i0)) {
        p <- p %>%
            add_trace(
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
            add_trace(
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
    
    if(!is.null(label)) {
        i2 = which(names %in% label)
        p <- p %>%
            add_annotations(
                x = x[i2],
                y = y[i2],
                text = names[i2],
                font = list(
                    size = 12,
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
    yrange <- c(-1,1)*max(abs(y))*1.05        
    xaxis = list( title = xlab, xrange = xrange )
    yaxis = list( title = ylab, yrange = yrange )    

    p <- p %>%
        layout(
            shapes = list(abline1,abline2),
            xaxis = xaxis, yaxis = yaxis,
            showlegend = showlegend,
            hovermode='closest',
            dragmode= 'select') %>%
        config(displayModeBar = displayModeBar)
    
    xann <- c(1,1)*0.95
    yann <- c(0,0.99)
    ann.text <- paste("UP in", group.names[c(2,1)])
    
    p <- p %>%
        add_annotations(
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
                          marker.size = 5, label=NULL, displayModeBar=TRUE )
{

    require(plotly)

    if(is.null(highlight)) highlight <- names
    i0 <- which(!names %in% highlight)
    i1 <- which(names %in% highlight)
    
    p <- plot_ly(
        type='scatter', mode='markers'
        ##type='scattergl', mode='markers',
        ##source=source, key=1:length(x)
    )
    
    p <- p %>%
        event_register('plotly_hover') %>%
        event_register('plotly_selected')

    if(length(i0)) {
        p <- p %>%
            add_trace(
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
            add_trace(
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
    
    if(!is.null(label)) {
        i2 = which(names %in% label)
        p <- p %>%
            add_annotations(
                x = x[i2],
                y = y[i2],
                text = names[i2],
                font = list(
                    size = 12,
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
    
    
    xrange <- c(-1,1)*max(abs(x))*1.05
    if(min(x)>=0) xrange <- c(0,1)*max(abs(x))*1.05
    yrange <- c(-1,1)*max(abs(y))*1.05
    if(min(y)>=0) yrange <- c(0,1)*max(abs(y))*1.05
        
    xaxis = list( title = xlab, xrange = xrange )
    yaxis  = list( title = ylab, yrange = yrange )    
    p <- p %>%
        layout(
            shapes = list(abline1,abline2,abline3),
            xaxis = xaxis, yaxis = yaxis,
            showlegend = showlegend,
            hovermode='closest',
            dragmode= 'select') %>%
        config(displayModeBar = displayModeBar)
    
    xann <- c(0.01,0.99)
    yann <- c(1,1)*1.04
    ann.text <- paste("UP in", group.names[c(2,1)])
    
    p <- p %>%
        add_annotations(
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
    dd <- as.dist(1 - cor(t(x),use="pairwise"))
    hc <- fastcluster::hclust(dd, method="ward.D2" )
    hc
}

## Override add_col_annotation to be able to suppress titles
##
##
library(iheatmapr)
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
                p <- add_col_groups(p, 
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
                p <- add_col_signal(p, 
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


##splitx="cell.family";top.mode="specific";row_annot_width=0.03;ntop=50
pgx.splitHeatmap <- function(ngs, splitx=NULL, top.mode="specific",
                             annot.pheno = NULL, row_annot_width=0.03,
                             scale="row.center", ntop=50, colors=NULL,
                             rowcex=rowcex, colcex=colcex)
{    
    require(iheatmapr)
    
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
        jj = as.vector(apply(dr,2,function(x) head(order(-abs(x)),ntop1)))
        ##jj = as.vector(apply(dr,2,function(x) head(order(-x),10)))
        X1 <- X0[jj,]
        idx <- paste0("PC",as.vector(mapply(rep,1:ncol(dr),ntop1)))
    } else if(top.mode == "specific") {
        grpX <- tapply(colnames(X0), splitx,function(k) rowMeans(X0[,k,drop=FALSE],na.rm=TRUE))
        grpX <- do.call(cbind, grpX)
        cat("dim(grpX)=",dim(grpX),"\n")
        cat("ntop=",ntop,"\n")
        ntop1 = ceiling(ntop/ncol(grpX))
        grpX <- grpX - rowMeans(grpX,na.rm=TRUE)  ## relative to others
        jj = as.vector(apply(grpX,2,function(x) head(order(-x),ntop1)))
        X1 <- X0[jj,]
        idx <- paste0("M",as.vector(mapply(rep,1:ncol(grpX),ntop1)))
    } else {
        X1 <- head(X0[order(-apply(X0,1,sd)),],ntop)
        hc = fastcluster::hclust( as.dist(1 - cor(t(X1),use="pairwise")), method="ward.D2" )
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

##X=head(ngs$X,100);annot=ngs$samples
##row_annot_width=0.03;colors=NULL;label_size=11;scale="row.center"
xtips=NULL;ytips=NULL;lmar=60

pgx.splitHeatmapFromMatrix <- function(X, annot, idx=NULL, splitx=NULL, 
                                       xtips=NULL, ytips=NULL, row_clust=TRUE,
                                       row_annot_width=0.03, scale="row.center",
                                       colors=NULL, lmar=60,
                                       rowcex=rowcex, colcex=colcex)
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
    
    ## ------ split Y-axis (genes) by factor
    if(!is.null(idx)) {        
        hc.order <- function(x) {
            suppressWarnings( dd <- as.dist(1 - cor(t(x),use="pairwise")) )
            ## cat("DBG >pgx-plotting:pgx.splitHeatmapX> sum.is.na.dd=",sum(is.na(dd)),"\n")
            if(sum(is.na(dd))) dd[is.na(dd)] <- 1
            hc <- fastcluster::hclust(dd, method="ward.D2" )
            rownames(x)[hc$order]
        }
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
    
    grid_params <- setup_colorbar_grid(
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
    hc <- hclust(as.dist(1 - cor(xx[[1]],use="pairwise")))
    dd <- 0.08 * ex ## left dendrogram width

    hovertemplate = "Row: %{y} <br>Column: %{x}<br>Value: %{z}<extra> </extra>"
    tooltip = setup_tooltip_options(
        prepend_row = "Row: ", prepend_col = "Column: ", 
        prepend_value = "Value: ") 

    x1 <- xx[[1]]
    ##x1 <- x1[nrow(x1):1,]
    plt <- main_heatmap(
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
        add_col_dendro(hc, size = 0.06) %>%
        ##add_row_dendro(hr, size = dd) %>%
        add_row_title("Genes") %>%
        add_col_title(names(xx)[1], side="top") %>%
        add_col_annotation(
            size = col_annot_height, buffer = 0.005, side="bottom",
            colors = colors0, show_title=TRUE,
            annotF[colnames(xx[[1]]),,drop=FALSE])
    ##plt
    length(xx)
    dim(X)
    if(ncol(X)<100 && colcex>0) {
        plt <- plt %>% add_col_labels(
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
            hc <- hclust(as.dist(1 - cor(x1,use="pairwise")))

            plt <- plt %>%
                add_main_heatmap(
                    x1, , name = "expression",
                    x = xtips[colnames(x1)],
                    y = ytips[rownames(x1)],
                    tooltip = tooltip,
                    ## hovertemplate = hovertemplate,
                    size=sizes[i], buffer = 0.007*ex) %>%
                ##add_col_clustering() %>%
                add_col_dendro(hc, size = 0.06) %>%
                add_col_title(names(xx)[i], side="top") %>%
                add_col_annotation(
                    size = col_annot_height, buffer = 0.005, side="bottom",
                    colors = colors0, show_title=FALSE,
                    data.frame(annotF[colnames(x1),,drop=FALSE]))
            
            if(ncol(X)<100 && colcex>0) {
                plt <- plt %>%
                    add_col_labels(side="bottom",
                                   size=0.15*ex,
                                   font=list(size=11*colcex))
            }
        }

    }

    ## ----------- row annotation (i.e. gene groups)
    if(!is.null(idx) && !is.null(row_annot_width) && row_annot_width>0 ) {
        plt <- add_row_annotation(
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
        
        plt <- add_row_labels(
            plt, side="right", ticktext = gnames,
            size=w*ex, font=list(size=s1) ) 
    }

    
    return(plt)
}

##annot=ngs$Y
annot.ht=4;cluster.samples=TRUE
pgx.plotPhenotypeMatrix0 <- function(annot, annot.ht=4, cluster.samples=TRUE)
{

    cvar <- pgx.getCategoricalPhenotypes(
        annot, min.ncat=2, max.ncat=10, remove.dup=FALSE)
    cvar
    if(length(unique(annot$group)) > 0.33*nrow(annot)) {
        cvar <- setdiff(cvar,"group")
    }

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
        hc <- hclust(dist(annotx))  ## cluster samples
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
        if(npar[i]>3 & !isnum[i]) klrs = rep(brewer.pal(8,"Set2"),99)[1:npar[i]]
        ##if(!is.binary[i] & !isnum[i]) klrs = rep(brewer.pal(8,"Set2"),99)[1:npar[i]]        
        ##if(npar[i]==2) klrs = rep(brewer.pal(2,"Paired"),99)[1:npar[i]]
        names(klrs) = sort(unique(annot.df[,i]))
        klrs = klrs[!is.na(names(klrs))]
        ann.colors[[prm]] = klrs
    }

    show.legend <- (!is.binary & npar <= 10)
    
    ha = HeatmapAnnotation(
        df = annot.df,
        ##df2 = annot.fvar[,,drop=FALSE],        
        col = ann.colors,
        na_col = "grey99",
        ##annotation_height = unit(annot.ht, "mm"),
        simple_anno_size = unit(annot.ht,"mm"),  ## BioC 3.8!!
        ##show_annotation_name = (i==ngrp),
        ##annotation_name_gp = gpar(fontsize=3.3*annot.ht),
        ##annotation_legend_param = aa
        show_legend = show.legend
    )
    
    show_colnames = TRUE
    show_colnames = (nrow(annot.df)<100)
    nullmat <- matrix(0,0,nrow(annot.df))
    colnames(nullmat) <- rownames(annot.df)
    colnames(nullmat) <- paste(colnames(nullmat),"     ")
    h <- Heatmap( nullmat,
                 top_annotation = ha,
                 show_column_names = show_colnames
                 )

    ## some space between heatmap and annotation
    nullmat2 <- matrix(0,0,nrow(annot.df)*0.15)
    h2 <- Heatmap( nullmat2)
    
    draw(h + h2,
         padding = unit(c(1,10,1,10),"mm"),
         adjust_annotation_extension = TRUE
         )
    
}


pgx.plotPhenotypeMatrix <- function(annot)
{
    ## ------- set colors
    colors0 = rep("Set2",ncol(annot))    
    names(colors0) = colnames(annot)
    
    grid_params <- setup_colorbar_grid(
        nrows = 3, 
        y_length = 0.32, 
        y_spacing = 0.39, 
        x_spacing = 0.17,
        x_start = 1.1, 
        y_start = 0.96)

    ## maximize plot area
    mar <- list(l = 150, r = 0, b = 5, t = 0, pad = 3)

    require(iheatmapr)
    plt <- NULL
    empty.X <- matrix(NA,nrow=1,ncol=nrow(annot))
    colnames(empty.X) <- rownames(annot)
    plt <- main_heatmap(
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
    hc <- hclust(dist(annotX))
    plt <- plt %>% add_col_dendro(hc, size = 20*ncol(annot)/6 )

    col_annot_height = 10
    plt <- plt %>%
        add_col_annotation(
            ##annotF[colnames(X),,drop=FALSE],
            annotF[,],
            size = col_annot_height, buffer = 0.005, side="bottom",
            colors = colors0, show_title=TRUE)    
    colcex = 1
    if(nrow(annot)<100 && colcex>0) {
        plt <- plt %>% add_col_labels(
                           side="bottom", size=20*ncol(annot)/6,
                           font=list(size=11*colcex))
    }
    
    return(plt)
}


pgx.plotGeneExpression <- function(ngs, probe, comp=NULL, logscale=TRUE,
                                   level="gene", grouped=FALSE, srt=0,
                                   collapse.others=TRUE, max.points=-1,
                                   group.names=NULL,
                                   main=NULL, xlab=NULL, ylab=NULL, names=TRUE )
{
    if(0) {
        comp=NULL;
        logscale=TRUE;level="gene";grouped=TRUE;srt=90;collapse.others=1;
        max.points=-1;main=NULL;xlab=NULL;ylab=NULL;names=TRUE;group.names=NULL
        comp=1;group.names=c("CTRL","treated")
        probe=ngs$genes$gene_name[1]
    }

    if(is.null(probe)) return(NULL)
    if(is.na(probe)) return(NULL)

    if(level=="gene" && !probe %in% rownames(ngs$X)) {
        frame() ## emtpy image
        return(NULL)
    }
    if(level=="geneset" && !probe %in% rownames(ngs$gsetX)) {
        frame() ## emtpy image
        return(NULL)
    }
    
    ## ------------- determine groups
    expmat  <- ngs$model.parameters$exp.matrix
    cntrmat <- ngs$model.parameters$contr.matrix
    expmat  <- expmat[rownames(ngs$samples),,drop=FALSE]
    
    if(class(comp)=="numeric") comp <- colnames(expmat)[comp]
    if(!is.null(group.names) && length(group.names)!=2) stop("group.names must be length=2")
    if(is.null(main)) main <- probe
    comp    

    if(grepl("_vs_|_VS_",comp) && is.null(group.names) ) {
        comp1 <- sub(".*:","",comp)  ## remove prefix
        group.names = strsplit(comp1,split="_vs_|_VS_")[[1]]
        ## determine if notation is A_vs_B or B_vs_A 
        if(is.POSvsNEG(ngs)) {
            ## A_vs_B or B_vs_A notation !!!
            group.names <- rev(group.names) ## reversed!!
        } 

    }
    
    if(!is.null(comp)) {
        ct <- expmat[,comp]
        names(ct) <- rownames(expmat)
        ct
        samples <- rownames(expmat)[which(ct!=0)]
        samples
        grp0.name=grp1.name=NULL
        if(!is.null(group.names)) {
            grp0.name <- group.names[1]
            grp1.name <- group.names[2]
        } else {
            grp1.name <- comp
            grp0.name <- "REF"
        }
        xgroup <- c("other",grp1.name,grp0.name)[1 + 1*(ct>0) + 2*(ct<0)]       
    } else  {
        samples <- rownames(expmat)
        xgroup <- ngs$samples$group
    }
    
    ## currently cast to character... :(
    names(xgroup) <- rownames(ngs$samples)
    table(xgroup)
    jj <- which(!(xgroup %in% xgroup[samples]))
    jj
    if(length(jj)>0 && collapse.others) {
        xgroup <- as.character(xgroup)
        xgroup[jj] <- "other"
    }
    names(xgroup) <- rownames(expmat)
    table(xgroup)

    if(class(xgroup)=="character") {
        xgroup <- as.character(xgroup)
        levels0 <- xgroup[!duplicated(xgroup)]
        if("other" %in% levels0) levels0 <- c(levels0[levels0!="other"],"other")
        xgroup <- factor(xgroup, levels=levels0)
    }

    
    ## ------------- set color of samples
    require(RColorBrewer)
    grps <- sort(setdiff(unique(xgroup),"other"))
    grps
    ngrp <- length(grps)
    grp.klr = rep(brewer.pal(12,"Paired"),99)[1:ngrp]
    names(grp.klr) <- grps
    ## if(is.null(comp) && grouped) grp.klr <- rep("grey60",length(grp.klr))
    if(any(grepl("other",xgroup))) {
        grp.klr <- c("other"="grey85", grp.klr)
    }
    grp.klr

    ## -------------- get expression value
    if(level=="geneset") {
        gx <- ngs$gsetX[probe,rownames(ngs$samples)]
    } else {
        gx <- ngs$X[probe,rownames(ngs$samples)]
    }
    if(!logscale) gx <- 2**(gx)

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
            y0 <- min(ylim) - 0.08*diff(ylim)
            text( bx[,1], y0, names(gx)[], adj=1, srt=srt,
                 xpd=TRUE, cex=ifelse(nx>10,0.6,0.9) )
        }
        title(main, cex.main=1.0)

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
        grp.klr1
        gx.b3plot( gx, xgroup, ##ylim=c(0,max(gx)*1.2),
                  col = grp.klr1, ylab=ylab, bee.cex=bee.cex,
                  max.points=max.points, xlab=xlab, names=names,
                  ## sig.stars=TRUE, max.stars=5,
                  las=3, cex.names=0.75, srt=srt)
        title(main, cex.main=1.0)
    }

}

pgx.plotOmicsNetwork <- function(ngs, gene=NULL, reduced=NULL, levels=c("gene","geneset"),
                                 contrast=NULL, layout=NULL, colorcluster=FALSE,
                                 hilight=NULL )
{
    require(igraph)
    require(visNetwork)
    require(gplots) ## for col2hex
    ## return(NULL)

    has.graph <- all(c("omicsnet","omicsnet.reduced") %in% names(ngs))
    has.graph
    if(!has.graph) return(NULL)

    ##gr <- ngs$genes_tsne_graph$graph
    gr <- ngs$omicsnet
    if(is.null(gr)) return(NULL)
    if(is.null(reduced) && is.null(gene)) gr <- ngs$omicsnet.reduced

    if(!is.null(gene)) {
        gene0 = paste0("{gene}",gene)
        gene0
        k <- which(V(gr)$name %in% c(gene, gene0))
        k
        nb <- names(neighbors(gr, V(gr)[k]))
        vv <- unique(c(gene0, nb))
        gr <- induced_subgraph(gr, vv)
    }
    gr
    gr <- induced_subgraph(gr, which(V(gr)$level %in% levels))
    if(is.null(gr)) return(NULL)

    ## ------------- get fold-change for node color and size ------------------
    fc0 <- gr$foldchange[V(gr)$name,,drop=FALSE]
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
    vlabel <- V(gr)$name
    if("label" %in% vertex_attr_names(gr)) vlabel <- V(gr)$label
    vlabel0 = vlabel
    vklr = c("blue3","grey70","red3")[2+sign(fc)]
    ##do.colorcluster = FALSE
    if(colorcluster) {
        vklr <- rep(rainbow(16),99)[V(gr)$cluster]
    }
    lab.cex = 1

    ee <- get.edgelist(gr)
    ee <- get.edges(gr, E(gr))
    head(ee)
    ##ew = rep(3,nrow(ee))
    ew = 1 + 5 * sqrt(fc.cex[ee[,1]]*fc.cex[ee[,2]])
    ew = 1 + 5 * abs(E(gr)$weight)
    vsel <- rep(0,length(fc))
    esel <- rep(0,nrow(ee))

    ## ------------------ highlight selection with labels
    if(!is.null(hilight)) {
        ##sel <- ngs$collections[[10]]
        sel <- hilight
        ##mm  <- ngs$genes_tsne_graph$members
        mm <- V(gr)$name
        if(!is.null(gr$members)) {
            mm <- gr$members[V(gr)$name]
        }
        mm <- sapply(mm, function(s) sub(".*\\}","",s))
        vlabel <- sapply( mm, function(x) intersect(x,sel) )
        vlabel <- sapply( vlabel, paste, collapse="\n")
        head(vlabel,10)
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


    vklr = substring(col2hex(vklr),1,7)
    names(vklr) <- V(gr)$name
    V(gr)$label <- vlabel ## filtered labels
    V(gr)$title <- gsub("\n","<br>",vlabel0)  ## tooltip has complete names
    ##V(gr)$label.cex <- 0.12 * lab.cex
    V(gr)$size <- 40*fc.cex
    V(gr)$color <- paste0(vklr,ifelse(vsel==1,"99","55"))
    ##if(input$gr_labels==FALSE) V(gr)$label <- NA
    ##E(gr)$color <- paste0(vklr[ee[,1]],ifelse(esel==1,"33","22"))
    E(gr)$color <- paste0(vklr[ee[,1]],ifelse(esel==1,"99","55"))
    E(gr)$width <- 1 * (2 + 5 * (ew / max(ew)))

    gr

    ##pos <- as.matrix(gr$layout)
    if(!is.null(layout)) {
        layout.fun <- match.fun(layout)
        tmp.gr = gr
        E(tmp.gr)$weight = abs(E(tmp.gr)$weight)**0.2
        ##E(tmp.gr)$weight = 1
        pos <- layout.fun(tmp.gr)
        remove(tmp.gr)
        rownames(pos) <- V(gr)$name
    }

    visdata <- toVisNetworkData(gr, idToLabel=FALSE)
    pos <- pos[V(gr)$name,]
    pos[,2] <- -pos[,2]
    dim(pos)


    ## ------------------ plot using visNetwork (zoomable) -----------------
    graph <- visNetwork(nodes = visdata$nodes, edges = visdata$edges,
                        height="1200px", width="1600px") %>%
        visNodes(font=list(size=14))  %>%
        visEdges(hidden=FALSE, width=2, color=list(opacity=0.9))  %>%
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
        ##visHierarchicalLayout(direction = "LR") %>%
        ##visInteraction(hideEdgesOnDrag = TRUE) %>%
        visIgraphLayout(layout="layout.norm", layoutMatrix=pos)
    graph

}



##gene1="CD4";gene2="CD8A";col="black";cex=1;k=11
pgx.cytoPlot <- function(ngs, gene1, gene2, cex=1, col="grey60",
                         cex.names=1, samples=NULL, k=11)
{
    library(MASS)
    library(RColorBrewer)

    ## some pretty colors
    ##k <- 9
    my.cols <- rev(brewer.pal(k, "RdYlBu"))
    ##gene1 = "CD8A"
    ##gene2 = "CD4"
    if(is.null(samples))
        samples <- colnames(ngs$X)
    samples <- intersect(samples, colnames(ngs$X))
    x1 <- ngs$X[gene1,samples]
    x2 <- ngs$X[gene2,samples]
    names(x1) <- samples
    names(x2) <- samples
    m1 <- mean(x1)
    m2 <- mean(x2)

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
    z1 <- kde2d( x1[j1], x2[j1], n=50)
    z2 <- kde2d( x1[j2], x2[j2], n=50)
    z3 <- kde2d( x1[j3], x2[j3], n=50)
    z4 <- kde2d( x1[j4], x2[j4], n=50)
    ##z0 <- kde2d( x1[], x2[], n=50)

    ##par(mfrow=c(1,1))
    plot(x1, x2, xlab=gene1, ylab=gene2, col=col, pch=19, cex=cex)
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

    if(!is.null(ngs$deconv)) {
        inferred.celltype <- ngs$deconv[[1]][["meta"]]
        dim(inferred.celltype)
        ##cex.names=2
        lab1 <- names(sort(-colSums(inferred.celltype[j1,])))[1:3]
        pos1 <- apply(cbind(x1, x2)[j1,],2,median)
        text( pos1[1], pos1[2], paste(lab1,collapse="\n"),cex=cex.names, xpd=TRUE)

        lab2 <- names(sort(-colSums(inferred.celltype[j2,])))[1:3]
        pos2 <- apply(cbind(x1, x2)[j2,],2,median)
        text( pos2[1], pos2[2], paste(lab2,collapse="\n"),cex=cex.names, xpd=TRUE)

        lab3 <- names(sort(-colSums(inferred.celltype[j3,])))[1:3]
        pos3 <- apply(cbind(x1, x2)[j3,],2,median)
        text( pos3[1], pos3[2], paste(lab3,collapse="\n"),cex=cex.names, xpd=TRUE)

    }

}
