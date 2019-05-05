pgx.plotGeneExpression <- function(ngs, probe, comp=NULL, logscale=TRUE,
                                   level="gene", grouped=FALSE, srt=0,
                                   collapse.others=TRUE, max.points=-1,
                                   group.names=NULL,
                                   main=NULL, xlab=NULL, ylab=NULL, names=TRUE )
{
    if(0) {
        comp=NULL;logscale=TRUE;level="gene";grouped=TRUE;srt=90;collapse.others=1;
        max.points=-1;main=NULL;xlab=NULL;ylab=NULL;names=TRUE;group.names=NULL
        grouped=0;collapse.others=0;
        comp=6;group.names=c("CTRL","treated")
        probe="TOPORS"
        probe="KCNN4"
    }
    
    ## ------------- determine groups
    expmat  <- ngs$model.parameters$exp.matrix
    cntrmat <- ngs$model.parameters$contr.matrix
    expmat <- expmat[rownames(ngs$samples),,drop=FALSE]
    
    if(class(comp)=="numeric") comp <- colnames(expmat)[comp]
    if(!is.null(group.names) && length(group.names)!=2) stop("group.names must be length=2")
    if(is.null(main)) main <- probe
    comp
    
    if(!is.null(comp)) {
        ct <- expmat[,comp]
        ct
        samples <- names(which(ct!=0))
        samples
        grp0.name=grp1.name=NULL
        if(grepl("_vs_|_VS_",comp) && is.null(group.names) ) {
            comp1 <- sub(".*:","",comp)  ## remove prefix
            comp1 <- strsplit(comp1,split="_vs_|_VS_")[[1]]
            grp1.name <- comp1[1]
            grp0.name <- comp1[2]
        } else if(!is.null(group.names)) {
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
                     las=3, ylab=ylab, names.arg=NA)
        if(nx<200 && names==TRUE) {
            y0 <- min(ylim) - 0.08*diff(ylim)
            text( bx[,1], y0, names(gx)[], adj=1, srt=srt,
                 xpd=TRUE, cex=ifelse(nx>30,0.6,0.9) )
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
