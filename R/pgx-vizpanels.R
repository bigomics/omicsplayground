##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

library(ggplot2)
library(patchwork)

## just to list functions in this file
viz.ClusterMarkers <- function(pgx){}
viz.PhenoMaps <- function(pgx){}
viz.PhenoStats <- function(pgx){}
viz.PhenoStatsBy <- function(pgx, by.pheno){}
viz.Expression <- function(pgx, pheno, genes){}
viz.GeneSetEnrichment <- function(pgx, genesets, contrast) {}
viz.Contrasts <- function(pgx, contrasts){}
viz.MitoRiboQC <- function(pgx){}
viz.VHVLusage <- function(pgx) {}
viz.BatchCorrection <- function(pgx, cX){}
viz.showFigure <- function(fig){}
viz._showPlotly <- function(fig){}
viz._showShiny  <- function(fig){}
viz._showGGplot  <- function(fig){}

if(0) {

    RDIR = "../R"
    FILES = "../lib"
    PGX.DIR = "../data"
    source("../R/pgx-include.R")
    ##source("options.R")
    FILES
    source("pgx-include.R")
    load("../data/geiger2016-arginine.pgx")
    viz.Contrasts(ngs, contrasts=NULL)
    
    
}


##cX2=NULL;pos0=pos1=pos2=NULL;npca=3;cex=1;nmax=40;main=c("not-corrected", "corrected","corrected2");pca.heatmap=FALSE
viz.CompareDatasets <- function(pgx1, pgx2, nmax=50, cex=1, 
                                main=c("heatmap1", "heatmap2"),
                                title=NULL, subtitle=NULL, caption=NULL)
{    
    ##-------------------------------------------
    ## Heatmaps
    ##-------------------------------------------    
    xlist <- list(X0=X0, X1=X1, X2=X2)
    hlist <- list()        
    for(i in 1:length(xlist)) {
        hlist[[i]] <- grid.grabExpr(
            gx.splitmap(
                xlist[[i]], main=main[i],
                col.annot=pheno, softmax=TRUE,
                show_legend=FALSE, scale="row", split=NULL,
                nmax = nmax, show_rownames = 1, 
                title_cex = 1.1, cexRow=0.7, cexCol=0.78,
                annot.ht=2.5, mar=c(4,1,1,10)*1,
                    key.offset=c(0.05,1.03),
                rownames_width = 10,
                show_colnames = ifelse(ncol(X0)<25,1,0)
            )
        )
    }
    
        
    ##-------------------------------------------
    ## Arrange plots
    ##-------------------------------------------
    np <- length(xlist)
    np
    lab1 <- letters[1:np]
    lab2 <- letters[(np+1):(2*np)]
    lab3 <- letters[(2*np+1):(3*np)]
    plt1 <- cowplot::plot_grid(plotlist=hlist, nrow=1, labels=lab1, label_size=15)
    plt2 <- cowplot::plot_grid(plotlist=flist, nrow=1, labels=lab2, label_size=15)
    plt3 <- cowplot::plot_grid(plotlist=plist, ncol=1, labels=lab3, label_size=15)
    fig  <- cowplot::plot_grid(
                         cowplot::plot_grid(plt1, plt2, ncol=1, rel_heights=c(1.95,1)),
                         plt3, ncol=2, rel_widths=c(np,1),
                         labels=c("","")
                     )

    ##if(is.null(title)) title <- "Batch effect analysis"
    ##if(is.null(subtitle)) subtitle = "The plots show possible batch effects your experiments."
    ##if(is.null(caption)) caption <- paste0("Project: ",pgx$name)
    viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption, tag=FALSE)    
    
}

viz.FoldChangeHeatmap <- function(pgx, comparisons=NULL, hilight=NULL,
                                  ntop=10, plot.diag=NULL, cex=1, nrow=NULL,
                                  title=NULL, subtitle=NULL, caption=NULL)                              
{
    ## 
    ##
    ##

    if(is.null(comparisons)) {
        ct <- names(pgx$gx.meta$meta)
    }
    
    out <- pgx.getMetaFoldChangeMatrix(pgx)
    F <- out$fc
    Q <- out$qv

    ii <- head(order(-rowMeans(F**2)),60)
    F1 <- F[ii,]
    F1 <- F1[order(rowMeans(F1)),]
    
    hm <- Heatmap( t(F1),
                  cluster_rows = TRUE,
                  cluster_columns = FALSE)

    p1 <- ggbarplot( t(F1), las=3,
                    legend.pos=c(1.2,1.05), legend.cex=0.9,
                    xlab="", ylab="cumulative fold-change") +
        theme(
            plot.margin = unit(c(0.3,5.6,0,0), "cm")
        )
    p2 <- grid.grabExpr(draw(hm))
    
    ##--------------------------------------------------
    ## Arrange
    ##--------------------------------------------------
    fig <- cowplot::plot_grid(p1, p2, nrow=nrow, rel_heights=c(1,1.8)) 
    fig
    
    if(is.null(title)) title = "Fold Change heatmap"
    if(is.null(subtitle)) subtitle = "Clustering of fold-changes signatures."    
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)

    viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption)
    
}

##pos="tsne2d";pheno1="cluster";pheno2="cell.type"
viz.ClusterMarkers <- function(pgx, pheno1, pheno2, n=NULL, pos="tsne2d",
                                  genes=NULL, subsample=300, scale='none', legend = TRUE,
                               theme=NULL, title=NULL, subtitle=NULL, caption=NULL)
{
    ##posx <- pgx$cluster$pos[[pos]]
    posx <- pgx$tsne2d
    X <- as.matrix(pgx$X)
    if(!is.null(genes)) {
        X <- X[intersect(genes,rownames(X)),]
    }
    if(length(legend)==1) legend <- c(legend,legend)
    
    ph1 <- as.character(pgx$samples[,pheno1])
    lg1 <- legend[1]
    p1 <- pgx.scatterPlotXY(
        posx, var=ph1, label.clusters=!lg1, legend=lg1,
        cex.clust=1.1, title=pheno1, plotlib="ggplot")
    if(!is.null(theme)) { p1 <- p1 + theme }
    
    if(!is.null(pheno2)) {
        ph2 <- pgx$samples[,pheno2]
        lg2 <- legend[2]
        p2 <- pgx.scatterPlotXY(
            posx, var=ph2, plotlib="ggplot",
            label.clusters=!lg2, legend=lg2,
            title=pheno2, cex.clust=1.1)
        if(!is.null(theme)) { p2 <- p2 + theme }        
    } else {
        p2 <- plot_spacer()
    }
    
    if(is.null(n)) n <- floor(65 / length(unique(ph1)))
    sel <- tapply(1:ncol(X), ph1, function(i) head(sample(i),subsample))
    sel <- unlist(sel)
    p3 <- grid.grabExpr(
        gx.markermap(
            X[,sel], splitx=ph1[sel], n=n,
            scale=scale, softmax=TRUE, title_cex=1,
            show_rownames=100, cexRow=0.8)
    )
    fig <- ((p1 / p2) | p3) + plot_layout(widths=c(1,2.5))
    
    if(is.null(title)) title = "Cluster markers"
    if(is.null(subtitle)) subtitle = "The heatmap shows the top markers genes for each phenotype."
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)
    viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption)

}


##pos="tsne2d";phenotypes=NULL;legend.pos="bottomleft"
viz.PhenoMaps <- function(pgx, phenotypes=NULL, pos=NULL, cex=1,
                             legend.pos = "bottomleft", theme=NULL,
                             title=NULL, subtitle=NULL, caption=NULL)
{
    posx <- pgx$tsne2d
    if(!is.null(pos) && is.character(pos) && pos %in% names(pgx$cluster$pos)) {
        posx <- pgx$cluster$pos[[pos]]
    }
    if(!is.null(pos) && is.matrix(pos)) {
        posx <- pos
    }

    if(is.null(phenotypes)) {
        phenotypes <- colnames(pgx$samples)
    } else {
        phenotypes <- intersect(phenotypes, colnames(pgx$samples))
    }

    plt <- list()
    for(ph in phenotypes) {
        y <- pgx$samples[,ph]
        p <- pgx.scatterPlotXY(
            posx, var=y, title=ph, cex=cex, legend.pos=legend.pos,
            plotlib="ggplot", theme=theme)
        ## if(!is.null(theme)) { p <- p + theme }
        p <- p + theme( plot.margin = ggplot2::margin(0,8,0,8) )
        plt[[ph]] <- p
    }

    if(is.null(title)) title = "Phenotype Maps"
    if(is.null(subtitle)) subtitle = "The plot show the phenotypes of your experiments."    
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)    

    ##fig <- patchwork::wrap_plots(plt) 
    ##viz.showFigure(fig=fig, title=title, subtitle=subtitle, caption=caption)
    viz.showFigure(fig=plt, title=title, subtitle=subtitle, caption=caption)
}

##dtype="bar";plotlib="ggplot";title=NULL;subtitle=NULL;caption=NULL;ctype="density";theme=NULL
viz.PhenoStats <- function(pgx, phenotypes=NULL,
                           dtype="bar", ctype="density", plot.lib="ggplot",
                           theme=NULL, title="*", subtitle=NULL, caption=NULL)
{
    if(is.null(phenotypes)) {
        phenotypes <- colnames(pgx$samples)
    } else {
        phenotypes <- intersect(phenotypes, colnames(pgx$samples))
    }

    df <- data.frame(pgx$samples[,phenotypes,drop=FALSE])
    sel <- which(colSums(!is.na(df))>1)
    phenotypes <- phenotypes[sel]
    phenotypes

    cat("[viz.PhenoStats] phenotypes = ",phenotypes,"\n")
    ptype <- sapply(df[,phenotypes,drop=FALSE],class)
    names(ptype) <- phenotypes
    ptype
    phenotypes
    df <- df[,phenotypes,drop=FALSE]

    i=1
    plotlist <- list()
    for(i in 1:length(phenotypes)) {
        i
        isfactor = TRUE
        isfactor = ptype[i] %in% c("factor","character","logical")
        isfactor
        ##plt <- NULL
        ## violin plot or barplot
        p <- phenotypes[i]
        p
        if(isfactor) {
            ##-----------------------------------------
            ## for factor/discrete variables
            ##-----------------------------------------
            ## get expression per group
            if(dtype=="bar") {
                plt <- ggplot(df, aes_string(x=p, fill=p)) +
                    geom_bar(width=0.85, position='dodge')
            } else if(dtype=="pie") {
                plt <- ggplot(df, aes_string(x=factor(1), fill=p)) +
                    geom_bar(width=1) + coord_polar("y")
            } else {
                stop("FATAL:: unknown dtype=",dtype)
            }
            if(!is.null(theme)) {
                plt <- plt + theme
            }
            plt <- plt + xlab("") + ggtitle(p) +
                ## scale_fill_brewer(palette="Dark2") + 
                theme(
                    legend.position='none',
                    legend.title = element_blank()
                )            
            if(length(unique(df[,p]))>10) {
                plt <- plt + theme(legend.position='none')
            }
            plotlist[[i]] <- plt
        } else {
            ##-----------------------------------------
            ## for continuous variables
            ##-----------------------------------------
            df1 <- data.frame(x=df[,p])
            if(ctype=="density") {
                plt <- ggplot(data=df1, aes(x=x)) + 
                    geom_density(fill="grey50", alpha=0.35, size=0.8) +
                    labs(title=paste("Density for ",p), x=p, y="density") 
            } else if(ctype=="histogram") {
                plt <- ggplot(data=df1, aes(x=x)) + 
                    geom_histogram(col="grey30", fill="grey70", size=0.4, alpha=0.8) + 
                    labs(title=paste("Histogram for ",p), x=p, y="count") 
            } else {
                stop("FATAL:: unknown ctype=",ctype)
            }
            if(!is.null(theme)) {
                plt <- plt + theme
            }
            ##plt
            plotlist[[i]] <- plt
        }
        ##plotlist[[p]] <- plt
    }

    if(!is.null(title) && title=="*") {
        title = "Phenotype distribution"
        subtitle = "The plot show the statistics of your phenotypes."    
        caption <- paste0("Project: ",pgx$name)
    }
    
    ##fig <- patchwork::wrap_plots(plotlist)    
    ##viz.showFigure(fig, title, subtitle, caption)
    viz.showFigure(fig=plotlist, title=title,
                   subtitle=subtitle, caption=caption)
}

viz.PhenoStatsBy <- function(pgx, by.pheno, phenotypes=NULL,
                             pct=FALSE, srt=0, xlab=by.pheno,
                             dtype="bar", ctype = "density",
                             title=NULL, subtitle=NULL, caption=NULL)
{
    ##
    ##
    ##
    
    if(0) {
        by.pheno="dlbcl.type";pct=TRUE;phenotypes=NULL;srt=45;
        xlab=title=subtitle=caption="?"
        dtype="bar";ctype="box"
    }

    if(is.null(phenotypes))
        phenotypes <- colnames(pgx$samples)    
    phenotypes <- intersect(phenotypes, colnames(pgx$samples))
    phenotypes <- setdiff(phenotypes, by.pheno)

    x <- pgx$samples[,by.pheno]
    p <- phenotypes[1]    
    Y <- data.frame(pgx$samples[,phenotypes,drop=FALSE])
    Y <- tidy.dataframe(Y)
    
    ptype <- sapply(Y, class)
    names(ptype) <- phenotypes
    nlev <- apply(Y,2,function(y) length(unique(y[!is.na(y)])))
    ptype[nlev<=5 & ptype %in% c("integer","numeric")] <- "discrete"
    ptype
    
    plotlist <- list()
    for(p in phenotypes) {
        isfactor = TRUE
        isfactor = ptype[p] %in% c("factor","character","logical","discrete")
        isfactor
        if(isfactor) {
            ##-----------------------------------------
            ## for factor/discrete variables
            ##-----------------------------------------
            y <- as.factor(Y[,p])
            F <- table(x,y)
            if(pct) {
                ##F <- F / rowSums(F)  ## normalize to 100% ??
                F <- t(t(F) / colSums(F))  ## normalize to 100% ??
                F <- round(100*F,2)
                xlab0 <- "percentage  (%)"
            }
            rownames(F) <- paste0(" ",rownames(F)," ")
            plt <- ggbarplot((F), base_size=12,
                             legend.pos = "right",
                             ylab = "counts  (N)", srt=srt,
                             xlab = "") 
            ## scale_x_discrete(guide=guide_axis(n.dodge=2)) +                        
            if(length(unique(y))>10) {
                plt <- plt + theme(legend.position='none')
            }
        } else {
            ##-----------------------------------------
            ## for numerical variables
            ##-----------------------------------------
            y = as.numeric(Y[,p])
            df <- data.frame(x=x, y=y)       
            if(ctype=="box") {
                suppressWarnings(
                    plt <- ggplot(df, aes(x=x, y=y, fill=x)) +
                        geom_boxplot() +
                        theme(
                            ##legend.position = 'none'
                            legend.position = 'right',
                            legend.title = element_blank()
                        ) +
                        xlab("") + ylab(p) +
                        labs(fill=element_blank())
                )
            } else if(ctype=="violin") {
                suppressWarnings(
                    plt <- ggviolin(
                        df$x, df$y, group=df$group, base_size=13,
                        srt=srt, pdodge=0.6, cex=cex, main=title) +
                        xlab("") + ylab(p)
                )
            } else if(ctype=="density") {
                plt <- ggplot(data=df, aes(x=y, color=x, fill=x)) + 
                    geom_density(alpha=0.15, size=1, position="stack") +
                    guides(fill=FALSE) +
                    labs(title=paste("Density for ",p),
                         x=p, y="density", color=by.pheno) 
            } else if(ctype=="histogram") {
                plt <- ggplot(data=df, aes(x=y, color=x, fill=x)) + 
                    geom_histogram(alpha=0.25, size=0.5, position="stack") +
                    guides(fill=FALSE) +
                    labs(title=paste("Histogram for ",p),
                         x=p, y="counts", color=by.pheno) 
            } else {
                stop("FATAL:: unknown ctype=",ctype)
            }
            
            plt <- plt + 
                theme_classic(base_size=12) 
            ## theme(legend.position = 'none')
            plt
        }

        cpal <- rep(brewer.pal(12, "Set3"),99)
        cpal <- rep(brewer.pal(12, "Paired"),99)
        suppressWarnings( suppressMessages(
            plt <- plt + ggtitle(p) +
                scale_color_manual(values=cpal) +
                scale_fill_manual(values=cpal) +
                theme(
                    plot.margin = ggplot2::margin(2,10,5,10),
                    legend.key.size = unit(9, "pt"),
                    legend.key.height = unit(9, "pt")
                )
        ))

        plotlist[[p]] <- plt        
    }

    if(is.null(title)) title = paste("Phenotype statistics by",by.pheno)
    if(is.null(subtitle)) subtitle = "The plots show the phenotypes of your experiments."    
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)
    ##fig <- patchwork::wrap_plots(plotlist)    
    ##viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption)
    viz.showFigure(fig=plotlist, title=title, subtitle=subtitle, caption=caption)    
}

##pheno="dlbcl.type";contrast="ABC_vs_GCB";nrow=3;ngenes=24;cex=1;srt=25
viz.Expression <- function(pgx, pheno, contrast, genes=NULL, 
                           ngenes=24, nrow=NULL, cex=1, srt=25,
                           title=NULL, subtitle=NULL, caption=NULL)
{
    ## get cluster markers (only works if vs_others is defined!!)
    require(ggpubr)       
    if(is.null(genes) && !is.null(contrast)) {
        px.markers <- pgx.getMarkerGenes(pgx, n=ngenes, dir=0)
        genes <- px.markers[[contrast]]
    } else if(!is.null(genes) && is.null(contrast)) {
        genes <- head(genes, ngenes)
    } else {
        sdx <- apply(pgx$X, 1, sd)
        genes <- head(rownames(pgx$X)[order(-sdx)],ngenes)
    }
    genes <- intersect(genes, rownames(pgx$X))
    genes <- head(genes,ngenes)    

    ## phenoplot
    pos <- pgx$tsne2d
    df <- data.frame( x=pos[,1], y=pos[,2],
                     pheno = pgx$samples[,pheno] )
    p1 <- ggpubr::ggscatter(df, "x", "y", color="pheno",
                            size=1.2*cex ) +
        xlab(colnames(pos)[1]) + ylab(colnames(pos)[2]) +
        theme(legend.title = element_blank())
    ##p1
    
    ## violin plot or barplot
    ## get expression per group
    df <- melt(pgx$X[genes,])
    cpal <- "grey70"
    if(!is.null(pheno)) {
        df$pheno <- pgx$samples[df$Var2,pheno]
        np <- length(unique(df$pheno))
        cpal <-  rep(c("#00AFBB", "#E7B800"),99)[1:np]
        cpal <-  rep(brewer.pal(12,"Paired"),99)[1:np]
        cpal <-  rep(brewer.pal(12,"Set3"),99)[1:np]
        cpal <-  rep(brewer.pal(8,"Dark2"),99)[1:np]
    }
    p2 <- ggpubr::ggbarplot(
                      df, x = "Var1", y = "value", fill = "pheno",
                      add = "mean_sd", palette = cpal,
                      position = position_dodge(width=0.8)) +
        scale_x_discrete(guide=guide_axis(angle=srt)) +
        xlab(NULL) + ylab("expression   (log2)") +  
        theme(
            axis.text.x = element_text(angle=srt, vjust=0, size=10),
            legend.justification = c(1,1),
            legend.position = c(1,1),
            legend.key.size = unit(9, "pt"),
            legend.key.height = unit(9, "pt")
        )
    ##p2

    ## make marker grid
    gene_plots <- list()
    g="Cd8a"
    zlim <- range(pgx$X[genes,])
    ## zlim <- NULL
    for(i in 1:length(genes)) {        
        ##m1 <- pgx.plotSampleProjection (
        m1 <- pgx.scatterPlot(
            pgx, gene=genes[i], title=genes[i],
            ##label.clusters=TRUE, ## labels=markers,
            cex=cex, cex.clust=1, cex.title=0.9, 
            legend = ifelse(i==1,TRUE,FALSE),
            cex.legend=0.7, legend.ysp=1.8, 
            zlim=zlim, opacity=0.8, plotlib="ggplot") 
        m1 <- m1 + theme(
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank()
                   )        
        gene_plots[[i]] <- m1
    }
    
    require(grid)
    require(cowplot)
    p2 <- p2 + theme(plot.margin = ggplot2::margin(5,5,5,25))
    row1 <- cowplot::plot_grid(p1, p2, nrow=1, rel_widths=c(1,4), labels=c("a","b"))
    row2 <- cowplot::plot_grid(plotlist=gene_plots, nrow=nrow) +
        theme(plot.margin = ggplot2::margin(5,5,5,15))
    fig <- cowplot::plot_grid(row1, row2, ncol=1, rel_heights=c(1,2), labels=c("","c"))

    if(is.null(title)) title = "Gene Expression"
    if(is.null(subtitle)) subtitle = "The plots show the gene expression of your experiments."    
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)

    class(fig)
    viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption)    
}


type="volcano";ntop=20
viz.GeneSetEnrichment <- function(pgx, genesets, contrast, pos=NULL,
                                  fill=NULL, ntop=20, type="volcano",
                                  plots.only=FALSE, plotlib="ggplot",
                                  strip=NULL,
                                  title=NULL, subtitle=NULL, caption=NULL)
{

    if(!is.null(genesets)) {
        F <- pgx.getMarkerGenes(pgx, n=ntop, dir=0, sym=TRUE)
        markers  <- F[[contrast]]
        gs.genes <- lapply(apply(pgx$GMT[,genesets]!=0,2,which),names)        
        top.genes <- gs.genes
        for(i in 1:length(gs.genes)) {
            top.genes[[i]] <- intersect(markers, gs.genes[[i]])
        }
    } else {
        out <- pgx.getTopGeneSets(pgx, n=15, ng=ntop, dir=0, sym=TRUE)
        gs.genes  <- out$genes[[contrast]]
        top.genes <- out$top.genes[[contrast]]
        names(top.genes)
    }
    
    fc <- pgx$gx.meta$meta[[contrast]]$meta.fx
    mq <- pgx$gx.meta$meta[[contrast]]$meta.q
    names(fc) <- names(mq) <- rownames(pgx$gx.meta$meta[[contrast]])
    
    pos1 <- pos
    if(type=="volcano" && is.null(pos)) {
        pos1 <- cbind(fc, -log10(mq))
        rownames(pos1) <- rownames(pgx$gx.meta$meta[[contrast]])
        colnames(pos1) <- c("difference (log2FC)","significance (-log10q)")
    }
    if(type=="rank" && is.null(pos)) {
        pos1 <- cbind(rank(-fc), fc)
        rownames(pos1) <- rownames(pgx$gx.meta$meta[[contrast]])
        colnames(pos1) <- c("rank","difference (log2FC)")
    }

    if(is.null(fill))
        fill <- fc    

    plist <- list()
    i=1
    for(i in 1:length(top.genes)) {
        tt <- names(top.genes)[i]
        if(!is.null(strip)) tt <- gsub(strip,"",tt)
        if(type=="enplot") {
            p1 <- ggenplot(fc, gs.genes[[i]], main=tt )
        } else {
            p1 <- pgx.scatterPlotXY(
                pos=pos1, var=fill, 
                cex=1, cex.title=0.8,
                hilight = gs.genes[[i]],
                ##hilight2 = NULL,
                hilight2 = top.genes[[i]],
                title = tt,
                plotlib=plotlib ## , set.par=FALSE
            )
        }
        p1
        plist[[i]] <- p1
    }

    if(plots.only) return(plist)
    if(plotlib=="base") return(NULL)
    
    if(is.null(title)) {
        if(!is.null(pos)) {
            title = paste0("Geneset Enrichment (custom)")
        } else {
            title = paste0("Geneset Enrichment (",type,")")
        }
    }
    if(is.null(subtitle)) subtitle = "The plots show the geneset enrichment of your experiments."
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)
    viz.showFigure(plist, title=title, subtitle=subtitle, caption=caption)    
}

viz.Contrasts <- function(pgx, contrasts=NULL, ntop=10, dir=1, pos=NULL,
                          psig=0.05, fc=0.20,
                          cex=1, type = c("pair","MA","volcano"),
                          level="gene", filt=NULL, label.type="box",
                          strip=NULL, plots.only = FALSE, 
                          title=NULL, subtitle=NULL, caption=NULL)
{
    ## 
    ##
    ##

    type <- type[1]
    if(!is.null(pos)) type <- "custom"
    if(level=="geneset" && is.null(pos)) {
        stop("FATAL:: geneset level requires position vector (for now)")
    }
    
    available.ct <- unique(names(pgx$gx.meta$meta))
    available.ct
    if(is.null(contrasts)) {
        contrasts <- available.ct
    }
    contrasts <- intersect(contrasts, available.ct)    
    if(length(contrasts)==0) {
        cat("ERROR:: No valid contrasts. Available contrasts: ",
            paste(available.ct,collapse=" "),"\n")
        return(NULL)
    }

    ##--------------------------------------------------
    ## Scatter plots
    ##--------------------------------------------------
    type
    plist <- list()
    i=1
    for(i in 1:length(contrasts)) {
        ct = contrasts[i]
        labels=gg=NULL
        if(level=="gene") {
            gg <- pgx.getMarkerGenes(
                pgx, dir=dir, n=ntop, filt=filt)[[ct]]
            labels <- rownames(pgx$X)
        } 
        if(level=="geneset") {
            out <- pgx.getTopGeneSets(
                pgx, dir=dir, n=ntop, filt=filt)
            gg <- out$gsets[[ct]]
            labels <- rownames(pgx$gsetX)            
        }
        if(!is.null(strip)) labels <- gsub(strip,"",labels)
        if(type=="pair") {
            grpx <- gsub(".*_vs_|@.*","",ct)
            grpy <- gsub(".*:|_vs_.*","",ct)
            ylab1 <- paste0("expression ",grpy," (logCPM)")
            xlab1 <- paste0("expression ",grpx," (logCPM)")
            p1 <- pgx.plotContrast(
                pgx, ct, plotlib="ggplot", level=level,
                psig=psig, fc=fc,
                hilight=gg, dir=dir, ntop=ntop, cex=cex) 
            p1 <- p1 + ## theme_classic(base_size=12) +
                xlab(xlab1) + ylab(ylab1)
        } else if(type=="MA") {
            p1 <- pgx.plotMA(
                pgx, ct, cex=cex, ntop=ntop, level=level,
                psig=psig, fc=fc,
                hilight=gg, plotlib="ggplot")
        } else if(type=="custom") {
            p1 <- pgx.scatterPlot(
                pgx, pos=pos, contrast=ct, hilight=gg,
                level=level, cex.lab=0.6, cex=cex, 
                labels=labels, label.type=label.type,
                title=NULL, plotlib="ggplot")
        } else {
            ## if(type=="volcano") {
            p1 <- pgx.Volcano(
                pgx, ct, cex=cex, ntop=ntop, level=level,
                psig=psig, fc=fc,
                hilight=gg, plotlib="ggplot")            
        }
        p1 <- p1 + ## theme_classic(base_size=12) +
            ggtitle(ct) +
            theme(
                plot.margin = ggplot2::margin(5,10,5,10),
                legend.position = "none",            
                axis.title=element_text(size=11)
            )
        plist[[ct]] <- p1
    }
    ##p1
    
    if(plots.only) return(plist)

    ##--------------------------------------------------
    ## Arrange
    ##--------------------------------------------------    
    if(is.null(title)) title = "Contrasts"
    if(is.null(subtitle)) subtitle = "The plots show the contrasts."    
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)

    viz.showFigure(plist, title=title, subtitle=subtitle, caption=caption)
    
}
 
##title=NULL;subtitle=NULL;caption=NULL
##hilight=NULL;ntop=20
viz.FoldChangePairs <- function(pgx, comparisons=NULL, hilight=NULL,
                                ntop=10, plot.diag=NULL, cex=1, nrow=NULL,
                                title=NULL, subtitle=NULL, caption=NULL)                              
{
    ## 
    ##
    ##

    ##--------------------------------------------------
    ## Main scatter plots
    ##--------------------------------------------------

    if(is.null(comparisons)) {
        ct <- names(pgx$gx.meta$meta)
        comparisons <- t(combn(head(ct,5),m=2))  ## restrict
        comparisons <- head(comparisons, 24)
    }

    out <- pgx.getMetaFoldChangeMatrix(pgx)
    F <- out$fc
    Q <- out$qv

    plotlist <- list()
    i=1
    for(i in 1:nrow(comparisons)) {
        ij <- comparisons[i,]
        pos <- F[,ij]
        qq  <- Q[,ij]
        top.gg <- hilight
        if(is.null(hilight)) {
            df <- sort(rowSums(pos**2))
            top.gg <- names(tail(df,ntop))
        }
        sig <- c("","top")[1 + 1*(rownames(pos) %in% top.gg)]
        sig <- c("","sig1","sig2")[1 + rowSums(qq < 0.05)]        
        cpal <- c("grey70","darkorange","red2")
        p2 <- pgx.scatterPlotXY(
            pos, var=sig, plotlib="ggplot", cex=1.0*cex,
            hilight=top.gg, legend=FALSE) +
            scale_fill_manual(values=cpal) +
            scale_color_manual(values=cpal)
        p2 <- p2 +
            geom_hline(yintercept=0, size=0.3, color="grey70") +
            geom_vline(xintercept=0, size=0.3, color="grey70") 

        rho <- cor(pos[,1],pos[,2])
        plot.diag1 <- ((is.null(plot.diag) && rho>0.3) || plot.diag==TRUE)
        if(plot.diag1) {
            p2 <- p2 + 
                geom_abline(slope=1, intercept=0, linetype='dotted', size=0.4) +
                geom_abline(slope=1, intercept=+1, linetype='dashed', size=0.4) +
                geom_abline(slope=1, intercept=-1, linetype='dashed', size=0.4)             
        }
        plotlist[[i]] <- p2
    }
    ##scat1 <- patchwork::wrap_plots(plotlist)
    
    ##--------------------------------------------------
    ## Arrange
    ##--------------------------------------------------
    fig <- cowplot::plot_grid(plotlist=plotlist, nrow=nrow) 
    ##fig
    
    if(is.null(title)) title = "Fold Change comparisons"
    if(is.null(subtitle)) subtitle = "Compare pairs of fold-change signatures."    
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)

    viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption)
    
}

viz.MitoRiboQC <- function(pgx, group, srt=0, pos="tsne2d",
                           title=NULL, subtitle=NULL, caption=NULL)
{
    nFeature_RNA <- colSums(pgx$counts>0)
    nCounts_RNA  <- colSums(pgx$counts)
    sel.mito <- grep("^mt-",rownames(pgx$counts),ignore.case=TRUE)
    sel.ribo <- grep("^rp[ls]",rownames(pgx$counts),ignore.case=TRUE)
    percent.mito <- colSums(pgx$counts[sel.mito,]) / nCounts_RNA * 100
    percent.ribo <- colSums(pgx$counts[sel.ribo,]) / nCounts_RNA * 100

    y <- pgx$samples[,group]
    v1 <- ggviolin(y, nCounts_RNA, group=NULL, srt=srt,
                   main="ncounts", ylab="counts")
    v2 <- ggviolin(y, nFeature_RNA, group=NULL, srt=srt,
                   main="nfeature", ylab="number")
    v3 <- ggviolin(y, percent.mito, group=NULL, srt=srt,
                   main="percent.mito", ylab="percentage")
    v4 <- ggviolin(y, percent.ribo, group=NULL, srt=srt,
                   main="percent.ribo", ylab="percentage")
    vv <- patchwork::wrap_plots(v1, v2, v3, v4, nrow=1)

    x <- pgx$tsne2d
    s1 <- ggscatterFILL(x, col=nCounts_RNA, barscale=0.5,  main="ncounts", gamma=0.7)
    s2 <- ggscatterFILL(x, col=nFeature_RNA, barscale=0.5, main="nfeature", gamma=0.7)
    s3 <- ggscatterFILL(x, col=percent.mito, barscale=0.5, main="percent.mito", gamma=1)
    s4 <- ggscatterFILL(x, col=percent.ribo, barscale=0.5, main="percent.ribo", gamma=2)
    ss <- patchwork::wrap_plots(s1, s2, s3 ,s4, nrow=1) 
    ##ss <- ss & theme_minimal()
    
    fig <- (vv / ss)
    
    if(is.null(title)) title = "Experiment QC"
    if(is.null(subtitle)) subtitle = "The plot show the quality of your experiments."
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)
    viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption)
}

viz.NormalizeCounts <- function(pgx, methods=NULL, post.qn=FALSE,
                                title=NULL, subtitle=NULL, caption=NULL)
{

    plotDensity <- function(xx, main) {
        if(nrow(xx)>1000) xx <- xx[sample(1:nrow(xx),1000),,drop=FALSE]
        if(ncol(xx)>100)  xx <- xx[,sample(1:ncol(xx),100)]
        dc <- reshape2::melt(xx)
        ##dc$value[dc$value==0] <- NA
        dc <- dc[dc$value>0,,drop=FALSE]
        dc$Var2 <- paste0("column.",dc$Var2)  ## colnames
        tt2 <- paste(nrow(counts),"x",ncol(counts))
        xx1 <- xx
        xx1[xx1==0] <- NA  ## non-zero
        avgx <- mean(xx1,na.rm=TRUE)
        sdx  <- sd(xx1,na.rm=TRUE)
        ggplot(dc, aes(x=value, color=Var2)) +
            geom_density() + xlab("log2(counts+1)") +
            theme( legend.position = "none") +
            ggtitle(main) +
            geom_vline( xintercept = avgx,
                       linetype="dashed",
                       color="grey40", size=0.3) +
            geom_vline( xintercept = avgx + sdx*c(-2,2),
                       linetype="dotted",
                       color="grey50", size=0.6)            

    }
    
    counts <- as.matrix(pgx$counts)
    ## counts <- t(t(counts) * 1e6*runif(ncol(counts)))
                
    xlist <- list()
    ##NORMALIZATION.METHODS <- c("none","mean","scale","NC","CPM","TMM","RLE","quantile")
    NORMALIZATION.METHODS <- c("none","scale","quantile","CPM","TMM","RLE")

    if(is.null(methods))
        methods <- NORMALIZATION.METHODS
    methods <- intersect(methods, NORMALIZATION.METHODS)
    for(m in methods) {
        xlist[[m]] <- pgx.countNormalization(counts, m)
    }

    if(post.qn) {
        doQN <- function(x) limma::normalizeQuantiles(x)
        doQN <- function(x) pmax(0.01*2**limma::normalizeQuantiles(log2(100*x+1))-1,0)
        for(i in 1:length(xlist)) {
            xlist[[i]] <- doQN( xlist[[i]] )
        }
    }
    
    plotlist <- list()
    for(i in 1:length(xlist)) {
        x1 <- log2(1 + xlist[[i]])
        plotlist[[i]] <- plotDensity(x1, names(xlist)[i])
    }    

    ## ss <- patchwork::wrap_plots(plotlist, nrow=2, ncol=4)
    fig <- patchwork::wrap_plots(plotlist)
    fig
    
    ## if(is.null(title)) title = "Experiment QC"
    ## if(is.null(subtitle)) subtitle = "The plot show the quality of your experiments."
    ## if(is.null(caption)) caption <- paste0("Project: ",pgx$name)
    viz.showFigure(plotlist, title=title, subtitle=subtitle, caption=caption)
    viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption)
}

##by.pheno="isotype";ng=60;nmin=8
viz.VHVLusage <- function(pgx, by.pheno="isotype", ng=30, nmin=1,
                          title=NULL, subtitle=NULL, caption=NULL)
{
    ##
    ##
    ##

    grep("^IG",rownames(pgx$counts),value=TRUE,ignore.case=TRUE)
    X <- pgx$counts
    ##X <- logCPM(pgx$counts, total=NULL)
    X <- X[grep("^IG",rownames(X)),]

    ## subtract background
    ##cX <- pmax(X - apply(X, 1, quantile, probs=0.5),0)
    ##cX <- pmax(t(t(X) - apply(X,2,median,na.rm=TRUE)), 0)
    ##X <- pmax(X - apply(X,1,median,na.rm=TRUE), 0)

    igX  <- X[grep("^IG[HKL]",rownames(X),ignore.case=TRUE),]
    igH  <- X[grep("^IGH[ADGEM]$|^IGH[ADGEM][1-9]",rownames(X),ignore.case=TRUE),]
    VH   <- X[grep("^IGHV",rownames(X),ignore.case=TRUE),]
    VL   <- X[grep("^IGLV|^IGKV",rownames(X),ignore.case=TRUE,value=TRUE),]
    isotype  <- sub("IGH","Ig",rownames(igH)[max.col(t(igH))])
    table(isotype)
    if(by.pheno %in% colnames(pgx$samples)) {
        y <- pgx$samples[,by.pheno]
    } else if(by.pheno=="isotype") {
        y <- isotype
    } else {
        stop("ERROR:: unknown phenotype=",by.pheno)
    }
    
    VH.avg <- tapply(1:ncol(VH), y, function(i) rowMeans(VH[,i,drop=FALSE]))
    VL.avg <- tapply(1:ncol(VL), y, function(i) rowMeans(VL[,i,drop=FALSE]))
    VH.avg <- do.call(cbind, VH.avg)
    VL.avg <- do.call(cbind, VL.avg)
    ## VH.avg <- VH.avg - rowMeans(VH.avg)
    ## VH.avg <- VL.avg - rowMeans(VL.avg)
    VH.avg <- VH.avg[order(-rowSums(abs(VH.avg))),]
    VL.avg <- VL.avg[order(-rowSums(abs(VL.avg))),]        

    B1 <- ggbarplot( t(head(VH.avg,ng)),
                    col=brewer.pal(12,"Set3"), las=3,
                    beside=FALSE, xlab="", ylab="cumulative log-expression",
                    main=paste("VH usage by",by.pheno)) 
    B2 <- ggbarplot( t(head(VL.avg,ng)), col=brewer.pal(12,"Set3"), las=3,
                    beside=FALSE, xlab="", ylab="cumulative log-expression",
                    main=paste("VL usage by",by.pheno)) 

    
    pgx$samples$isotype  <- sub("IGH","Ig",rownames(igH)[max.col(t(igH))])
    pgx$samples$isotype0 <- substring(pgx$samples$isotype,1,3)
    pgx$samples$VH  <- rownames(VH)[max.col(t(VH))]
    pgx$samples$VL  <- rownames(VL)[max.col(t(VL))]
    head(pgx$samples)
    
    T1 <- unclass(table(pgx$samples$VL, y))
    T2 <- unclass(table(pgx$samples$VH, y))
    T1 <- head(T1[order(-rowSums(T1)),], ng)
    T2 <- head(T2[order(-rowSums(T2)),], ng)        

    VH.avg <- VH.avg[order(-rowSums(abs(VH.avg))),]
    VL.avg <- VL.avg[order(-rowSums(abs(VL.avg))),]        
    VH.avg <- VH.avg[order(-apply(VH.avg,1,sd)),]
    VL.avg <- VL.avg[order(-apply(VL.avg,1,sd)),]
    T1 <- head(VL.avg,ng)
    T2 <- head(VH.avg,ng)
    
    H1 <- grid.grabExpr(
        gx.splitmap(t(T1), split=NULL, cexCol=0.85, scale="none",
                    main=paste("VL usage by",by.pheno))
    )
    H2 <- grid.grabExpr(
        gx.splitmap(t(T2), split=NULL, cexCol=0.85,
                    main=paste("VH usage by",by.pheno))
    )
    
    if(1) {
        mm <- pgx.phenoMatrix(pgx, by.pheno)
        matlist <- list(VL, VH, mm)            
        S1 <- pgx.SankeyFromMatrixList.PLOTLY(matlist)
        ##S1 <- pgx.SankeyFromPhenotypes.GGPLOT(
        ##    pgx, phenotypes=c(by.pheno,"VL","VH"), fill=by.pheno,
        ##    nmin=nmin, title="VH/VL pairing")
        S1.grob <- NULL
        S1.grob <- plotly2ggplot(S1, width=600, height=800, scale=0.97)
        class(S1.grob)
    }

    ## arrange
    c1  <- cowplot::plot_grid(H2, H1, B1, B2, nrow=2, labels=c("a","b","c","d"))
    fig <- cowplot::plot_grid(c1, S1.grob, ncol=2, labels=c("","e"), rel_widths=c(1.5,1))
    fig
    
    if(is.null(title)) title  <- paste("VH/VL usage by",by.pheno)
    if(is.null(subtitle)) subtitle  <- "The plot show the VH/VL usage of your experiments."    
    if(is.null(caption)) caption <- paste0("Project: ",pgx$name)
    viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption, tag=FALSE)

}

if(0) {
    X0=cX=cX2=NULL;main1="not-corrected";pca.heatmap=FALSE;nmax=40;
    pos1=pos2=NULL;npca=3
    main=c("not-corrected", "batch-corrected", "batch-corrected2")
}


viz.GeneFamilies <- function(pgx, by.pheno=NULL, gset=NULL, ntop=20, srt=0,
                             sort="avg", mode="pct", lab.cex=1)
{
    if(is.null(gset)) {
        fam <- gsub("[0-9]","",names(tail(sort(table(substring(rownames(pgx$X),1,3))),100)))   
        fam <- setdiff(c(fam,"S100"),c("S"))
        fam <- sort(unique(fam))
        gset <- lapply(fam, function(f) grep(paste0("^",f),rownames(pgx$counts),value=TRUE))
        names(gset) <- fam
    }

    genes <- rownames(pgx$counts)
    gset <- lapply(gset, function(f) intersect(f,genes))
    ylab <- NULL
    if(mode=="pct") {
        fx <- t(sapply(gset, function(gg) colSums(pgx$counts[gg,])))
        rownames(fx) <- names(gset)
        fx <- t(t(fx) / colSums(pgx$counts)) * 100
        ylab <- "relative abundance  (%)"
    } else if(mode=="avg") {
        fx <- t(sapply(gset, function(gg) colMeans(pgx$X[gg,])))
        rownames(fx) <- names(gset)
        ## fx <- t(t(fx) / colSums(pgx$counts)) * 100
        ## fx <- fx - mean(fx)
        ## fx <- fx - rowMeans(fx)
        ylab <- "average expression  (logCPM)"
    } else {
        stop("FATAL ERROR:: invalid mode", mode)
    }
    dim(fx)

    if(!is.null(by.pheno) && sort=="diff") {
        y <- pgx$samples[,by.pheno]
        mx <- tapply(1:ncol(fx), y, function(i) rowMeans(fx[,i,drop=FALSE]))
        mx <- do.call(cbind,mx)
        fx <- fx[order(-apply(mx,1,function(x) diff(range(x)))),]    
    } else {
        fx <- fx[order(-rowSums(fx)),]    
    }


    ## barplot(log2(1+rowMeans(head(fx,50))), las=3, ylim=c() )
    df <- reshape2::melt(head(fx,ntop))
    ##by.pheno = "Chemotherapy"
    df$pheno <- ""
    if(!is.null(by.pheno)) {
        df$pheno <- pgx$samples[df$Var2,by.pheno]
    }
    colnames(df) <- c("gene family","sample","value","pheno")
    dim(df)
    
    nlev <- length(unique(df$pheno))
    nlev
    xoff <- 0
    if(1) {
        xoff <- min(df$value)
        df$value <- df$value - xoff
    }
    xbreaks <- seq(min(df$value),max(df$value),length.out=7)
    xbreaks
    
    if(nlev<=5) {
        plt <- ggpubr::ggbarplot(
                          df, x = "gene family", y = "value", fill = "pheno",
                          add = "mean_se", ## palette = c("#00AFBB", "#E7B800"),
                          position = position_dodge(width=0.8)) +
            ##scale_x_discrete(guide=guide_axis(angle=srt)) +
            xlab("gene family") + ylab(ylab) + 
            labs(fill = by.pheno, color = by.pheno) +
            scale_x_discrete(guide=guide_axis(angle=srt)) +
            ## xlim(c(min(df$value), max(df$value))) +
            scale_y_continuous(breaks=xbreaks, labels=round(xbreaks+xoff,2)) +
            theme(
                axis.text.x = element_text(angle=srt, vjust=0.2, size=10*lab.cex),
                ## axis.text.y = element_text(angle=srt, vjust=0.2, size=10*lab.cex),
                axis.text.y = element_text(size=10*lab.cex),
                legend.title = element_text(by.pheno),
                legend.justification = c(1,1),
                legend.position = c(1,1),
                legend.key.size = unit(9, "pt"),
                legend.key.height = unit(9, "pt"),
                ## scale_y_continuous(expand = c(0, 0)),                
                plot.margin = ggplot2::margin(5,5,20,5)
            )
        
        if(nlev==1) {
            plt <- plt + theme(legend.position = 'none')
        }
    } else {
        plt <- ggpubr::ggbarplot(
                           df,
                           x = "pheno", y = "value", fill = "gene family",
                           add = "mean",
                           ## palette = c("#00AFBB", "#E7B800"),
                           position = position_stack()) +
            scale_y_continuous(breaks=xbreaks, labels=round(xbreaks + xoff,1)) +            
            xlab(by.pheno) + ylab(ylab)            
    }
    if(nlev==1) {
        plt <- plt + ggtitle("Gene families")
    } else {
        plt <- plt + ggtitle(paste("Gene families by",by.pheno))
    }
    plt
}

viz.BatchCorrection <- function(pgx, cX, cX2=NULL, phenotype, stat="F", 
                                pca.heatmap=FALSE, nmax=40, cex=1,
                                pos1=NULL, pos2=NULL, npca=3, pheno=NULL,
                                main=c("not-corrected", "corrected","corrected2"),
                                title=NULL, subtitle=NULL, caption=NULL)
{
    X0 <- pgx$X  ## ???
    if(is.null(pheno))  pheno <- pgx$samples
    pos0 <- pgx$tsne2d
    if(0) {
        phenotype='Treatment'
        viz.BatchCorrectionMatrix(X0=X0, pheno=pheno, cX=cX, cX2=cX2,
                                  phenotype=phenotype)
    }
    viz.BatchCorrectionMatrix(
        X0=X0, pheno=pheno, cX=cX, cX2=cX2,
        phenotype=phenotype, stat=stat, 
        pca.heatmap=pca.heatmap, nmax=nmax, cex=cex,
        pos0=pos0, pos1=pos1, pos2=pos2, npca=npca,
        main=main, title=title, subtitle=subtitle, caption=caption)    
}


viz.BatchCorrectionMatrix <- function(X0, pheno, cX, cX2=NULL, phenotype, stat="F", 
                                      pca.heatmap=FALSE, nmax=40, cex=1, 
                                      pos0=NULL, pos1=NULL, pos2=NULL, npca=3,
                                      main=c("not-corrected", "corrected","corrected2"),
                                      title=NULL, subtitle=NULL, caption=NULL)
{
    ##X0 <<- pgx$X  ## ???
    ##X0 <- pgx$X  ## ???
    ## cX <<- cX  ## ???
    ## cX2 <<- cX2  ## ???
    ##pheno <- pgx$samples
    if(is.null(phenotype))
        phenotype <- colnames(pheno)[1]
    phenotype <- head(intersect(phenotype,colnames(pheno)),4) ## max 4
    phenotype

    if(1) {
        is.na1=is.na2=is.na3=0
        is.na1 <- rowMeans(is.na(X0)) > 0.9
        if(!is.null(cX)) is.na2 <- rowMeans(is.na(cX)) > 0.9
        if(!is.null(cX2)) is.na3 <- rowMeans(is.na(cX2)) > 0.9
        sel <- which(!is.na1 & !is.na2 & !is.na3)
        X0 <- X0[sel,,drop=FALSE]
        cX <- cX[sel,,drop=FALSE]
        if(!is.null(cX2)) cX2 <- cX2[sel,,drop=FALSE]
    }
    
    ##X1 <- head(X1[order(-apply(X1,1,sd)),],50)
    xlist <- list(X0, cX, cX2)
    pos <- list()
    pos[[1]] <- pos0
    pos[[2]] <- pos1
    pos[[3]] <- pos2
    if(is.null(pos0)) {
        pos[[1]] <- pgx.clusterBigMatrix(xlist[[1]], method="tsne", dims=2)[[1]]
    }
    if(is.null(pos1) && !is.null(cX)) {
        pos[[2]] <- pgx.clusterBigMatrix(xlist[[2]], method="tsne", dims=2)[[1]]
    }
    if(is.null(pos2) && !is.null(cX2)) {
        pos[[3]] <- pgx.clusterBigMatrix(xlist[[3]], method="tsne", dims=2)[[1]]
    }
    pos <- pos[!sapply(xlist,is.null)]
    xlist <- xlist[!sapply(xlist,is.null)]
    length(xlist)
    
    ##-------------------------------------------
    ## Heatmaps
    ##-------------------------------------------    
    hlist <- list()
    if(pca.heatmap) {

        for(i in 1:length(xlist)) {
            hlist[[i]] <- grid.grabExpr(
                ##gx.splitmap(
                gx.PCAheatmap( 
                    X=xlist[[i]], main=main[i],
                    nv=5, ngenes=10, ## split=NULL,               
                    col.annot=pheno, softmax=TRUE,
                    show_legend=FALSE, scale="row", 
                    nmax = nmax, show_rownames = 40, 
                    title_cex = 1.1, cexRow=0.7, cexCol=0.78,
                    annot.ht=2.5, mar=c(4,1,1,10)*1,
                    key.offset=c(0.05,1.03),
                    show_colnames = ifelse(ncol(X0)<25,1,0)
                )
            )
        }
        
    } else {
        for(i in 1:length(xlist)) {
            hlist[[i]] <- grid.grabExpr(
                gx.splitmap(
                    xlist[[i]], main=main[i],
                    col.annot=pheno, softmax=TRUE,
                    show_legend=FALSE, scale="row", split=NULL,
                    nmax = nmax, show_rownames = 1, 
                    title_cex = 1.1, cexRow=0.7, cexCol=0.78,
                    annot.ht=2.5, mar=c(4,1,1,10)*1,
                    key.offset=c(0.05,1.03),
                    rownames_width = 10,
                    show_colnames = ifelse(ncol(X0)<25,1,0)
                )
            )
        }
    }
    
    ##-------------------------------------------
    ## PCA variance plots (F-statistics)
    ##-------------------------------------------    
    flist <- list()
    pheno1 <- pheno[,order(colnames(pheno)),drop=FALSE]
    i=1
    for(i in 1:length(xlist)) {
        f1 <- pgx.PC_correlation(
            xlist[[i]], pheno1, nv=npca, stat="F", plot=TRUE,
            main = paste0("PC variance (",main[i],")"))
        f1 <- f1 + theme(plot.margin = ggplot2::margin(4,4,0,4,"mm"),
                         legend.justification = c(0,1),
                         legend.position = c(0.01,1) )        
        flist[[i]] <- f1
    }
    
    ##-------------------------------------------
    ## Scatterplots
    ##-------------------------------------------    
    plist <- list()
    y1 <- pheno[,phenotype[1]]
    y2 <- NULL
    if(length(phenotype)>1) y2 <- pheno[,phenotype[2]]
    for(i in 1:length(xlist)) {
        lg = 'right'
        if(length(unique(y1))>20) lg = 'none'
        plist[[i]] <- ggscatter(
            pos[[i]], col=y1, shape=y2,
            cex=0.7*cex) +
            ## theme(legend.position="top") +
            ## theme_classic() +
            xlab("tSNE-1") + ylab("tSNE-2") + 
            ggtitle(paste0(phenotype[1]," (",main[i],")")) +
            theme(
                legend.position = lg,
                plot.margin = ggplot2::margin(2,2,0,2,"mm")
            )
    }
        
    ##-------------------------------------------
    ## Arrange plots
    ##-------------------------------------------
    np <- length(xlist)
    np
    lab1 <- letters[1:np]
    lab2 <- letters[(np+1):(2*np)]
    lab3 <- letters[(2*np+1):(3*np)]
    plt1 <- cowplot::plot_grid(plotlist=hlist, nrow=1, labels=lab1, label_size=15)
    plt2 <- cowplot::plot_grid(plotlist=flist, nrow=1, labels=lab2, label_size=15)
    plt3 <- cowplot::plot_grid(plotlist=plist, ncol=1, labels=lab3, label_size=15)
    fig  <- cowplot::plot_grid(
                         cowplot::plot_grid(plt1, plt2, ncol=1, rel_heights=c(1.95,1)),
                         plt3, ncol=2, rel_widths=c(np,1),
                         labels=c("","")
                     )

    ##if(is.null(title)) title <- "Batch effect analysis"
    ##if(is.null(subtitle)) subtitle = "The plots show possible batch effects your experiments."
    ##if(is.null(caption)) caption <- paste0("Project: ",pgx$name)
    viz.showFigure(fig, title=title, subtitle=subtitle, caption=caption, tag=FALSE)    
    
}


viz.System <- function(pgx, contrast, umap, gs.umap)
{

    ## determine clusters
    cl1 <- pgx.FindClusters(t(umap), method="kmeans")[[1]][,"kmeans.10"]
    cl2 <- pgx.FindClusters(t(gs.umap), method="kmeans")[[1]][,"kmeans.5"]
    ##cl3 <- pgx.FindClusters(t(pgx$tsne2d), method="kmeans")[[1]][,"kmeans.4"]
    cl1 <- paste0("G",cl1)
    cl2 <- paste0("S",cl2)
    ##cl3 <- paste0("C",cl3)
    cl3 <- pgx$samples$cluster

    ## contrast plots
    plotlib="ggplot"
    cluster.plots <- function(contrast, plotlib="ggplot") {
        
        p1 <- pgx.scatterPlotXY(
            umap, var=factor(cl1), label.clusters=TRUE, legend=FALSE,
            cex=1, cex.clust=0.85, plotlib=plotlib, title="gene clusters")
        p2 <- pgx.scatterPlotXY(
            gs.umap, var=factor(cl2), label.clusters=TRUE, legend=FALSE,
            cex=1, cex.clust=0.85, plotlib=plotlib, title="geneset clusters")
        p3 <- pgx.scatterPlotXY(
            pgx$tsne2d, var=factor(cl3), label.clusters=TRUE, legend=FALSE,
            cex=1, cex.clust=0.85, plotlib=plotlib, title="sample clusters")    

        ##p4 <- viz.Contrasts(pgx, contrasts=contrast, pos=umap, plots.only=TRUE)[[1]]
        ##p5 <- viz.Contrasts(pgx, contrasts=contrast, level="geneset",
        ##                    filt="HALLMARK", strip=".*HALLMARK_",
        ##                    ntop=5, plots.only=TRUE, pos=gs.umap)[[1]]
        ##p6 <- pgx.scatterPlot(pgx, pheno="celltype", cex.legend=0.7,
        ##                      plotlib=plotlib, title="celltype")
        p4 <- pgx.scatterPlot(
            pgx, contrast=contrast, pos=umap, plotlib=plotlib, title="gene")
        p5 <- pgx.scatterPlot(
            pgx, contrast=contrast, pos=gs.umap, level="geneset",
            plotlib=plotlib, title="geneset")
        p6 <- pgx.scatterPlot(
            pgx, contrast=contrast, plotlib=plotlib, title="sample")

        plot_grid(p1, p2, p3, p4, p5, p6, nrow=2)        
        ##p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(nrow=2)
    }
    
    ##----------------------------------------------
    ## Setup multi-layer network
    ##----------------------------------------------    
    X1 <- pgx$X - rowMeans(pgx$X)
    X2 <- pgx$gsetX - rowMeans(pgx$gsetX)
    X1 <- X1 / apply(X1,1,sd) ## scale SD -> 1
    X2 <- X2 / apply(X2,1,sd) ## scale SD -> 1

    M1 <- matGroupMeans(X1, group=cl1, dir=2)
    M2 <- matGroupMeans(X2, group=cl2, dir=2)
    M1 <- M1 - rowMeans(M1)
    M2 <- M2 - rowMeans(M2)
    M1 <- M1 / apply(M1,1,sd)
    M2 <- M2 / apply(M2,1,sd)
    
    y1 <- pgx$samples$cluster
    P1 <- t(model.matrix(~ 0 + y1))
    rownames(P1) <- sub("^y1","",rownames(P1))
    colnames(P1) <- rownames(pgx$samples)

    sel2 <- grep("celltype|cell.type",colnames(pgx$samples))[1]
    y2 <- pgx$samples[,sel2]
    P2 <- t(model.matrix(~ 0 + y2))
    rownames(P2) <- sub("^y2","",rownames(P2))
    colnames(P2) <- rownames(pgx$samples)
    
    M1.top <- tapply( rownames(X1), cl1, function(i) head(names(sort(rowMeans(-X1[i,])))))
    M2.top <- tapply( rownames(X2), cl2, function(i)
        head(grep("HALLMARK",names(sort(rowMeans(-X2[i,]))),value=TRUE)) )
    M2.top <- lapply(M2.top, function(s) gsub(".*:|.*HALLMARK_","",s))
    M1.top <- as.character(sapply(M1.top,paste,collapse="\n"))
    M2.top <- as.character(sapply(M2.top,paste,collapse="\n"))
    names(M1.top) = rownames(M1)
    names(M2.top) = rownames(M2)
    P1.top <- array(rownames(P1), dimnames=list(rownames(P1)))
    P2.top <- array(rownames(P2), dimnames=list(rownames(P2)))
    
    ##------------------------------------------------------------------------
    ## Prepare matrices 
    ##------------------------------------------------------------------------

    xlist <- list(t(X1), t(X2), t(P1), t(P2))
    matlist <- list(t(M1), t(M2), t(P1), t(P2))
    ##classlist <- list(cl1, cl2, y1, y2)
    toplist <- list( M1.top, M2.top, P1.top, P2.top)

    fc.plots <- function(contrast, dir) {
        ##dir <- as.numeric(dir)
        y0 <- pgx$model.parameters$exp.matrix[,contrast]
        ##rho <- lapply(xlist, function(x) cor(x,y0)[,1])
        fc1 <- rowMeans(X1[,y0>0]) - rowMeans(X1[,y0<0])
        fc2 <- rowMeans(X2[,y0>0]) - rowMeans(X2[,y0<0])
        fc3 <- rowMeans(P1[,y0>0]) - rowMeans(P1[,y0<0])
        fc4 <- rowMeans(P2[,y0>0]) - rowMeans(P2[,y0<0])
        
        fc2 <- fc2[grep("HALLMARK",names(fc2))]
        names(fc2) <- gsub(".*HALLMARK_","",names(fc2))    
        if(dir=="both") {
            i1 <- head(order(-abs(fc1)),20)
            i2 <- head(order(-abs(fc2)),20)
            i3 <- head(order(-abs(fc3)),20)
            i4 <- head(order(-abs(fc4)),20)
        } else {
            xdir <- ifelse(dir=="up",1,-1)
            i1 <- head(order(-xdir*fc1),20)
            i2 <- head(order(-xdir*fc2),20)
            i3 <- head(order(-xdir*fc3),20)
            i4 <- head(order(-xdir*fc4),20)
        }
        f1 <- ggbarplot(sort(fc1[i1]),col="grey80",main="gene",ylab="logFC",xlab="") 
        f2 <- ggbarplot(sort(fc2[i2]),col="grey80",main="geneset",ylab="logFC",xlab="") 
        f3 <- ggbarplot(sort(fc3[i3]),col="grey80",main="cluster",ylab="logFC",xlab="") 
        f4 <- ggbarplot(sort(fc4[i4]),col="grey80",main="cell.type",ylab="logFC",xlab="") 
        
        fc.plots <- f1 + f2 + f3 + f4 + plot_layout(ncol=4) &
            coord_flip() & ## xlab("") & ylab("logFC") &
            theme_minimal() & theme(legend.position="none") &
            theme(
                axis.title.x = element_text(size=9),
                axis.title.y = element_text(size=9)
            ) 
            
        fc.plots
    }
    
    ##------------------------------------------------------------------------
    ## Prepare multi-layer system
    ##------------------------------------------------------------------------

    ## Counts cross-table between matrices
    M <- list()
    for(i in 1:(length(matlist)-1)) {
        k1 <- colnames(matlist[[i]])[max.col(matlist[[i]])]
        k2 <- colnames(matlist[[i+1]])[max.col(matlist[[i+1]])]        
        ##M[[i]] <- log(1+table(k1,k2))
        M[[i]] <- table(k1,k2)
    }
    
    ## Correlation
    R <- list()
    for(i in 1:(length(matlist)-1)) {
        r1 <- cor( matlist[[i]], matlist[[i+1]] )
        R[[i]] <- pmax(r1,0)
    }

    ## Edge value (i.e. capacity) : rho * contrast/FC
    cty.mode=1
    ##wt.contrast=FALSE
    plot.sk <- function(contrast, fill=TRUE) {
        F = R
        cty <- sign(pgx$model.parameters$exp.matrix[,contrast])        
        fc <- lapply(matlist, function(m) cor(m,cty)[,1])
        i=1
        for(i in 1:length(R)) {
            if(cty.mode==1) node.wt <- outer(pmax(fc[[i]],0), pmax(fc[[i+1]],0))
            if(cty.mode==3) node.wt <- abs(outer(fc[[i]], fc[[i+1]]))
            if(cty.mode==2) node.wt <- pmax(outer(fc[[i]], fc[[i+1]]),0)
            ww <- R[[i]] * node.wt
            ww <- ww / max(ww)  ## normalize??
            F[[i]] <- ww
        }
        pgx.SankeyFromMRF.PLOTLY(
            M=M, R=R, F=F, fill=fill, labels=toplist)
    }

    ##plot.sk(contrast, fill=TRUE)
    
    ## ---------------------------------------------------------
    ## ---------------------------------------------------------
    ## ---------------------------------------------------------

    all.ct <- colnames(pgx$model.parameters$exp.matrix)

    ## list of plot objects
    ## plots <- list(sk,p1,p2,p3,p4,p5,p6)
    params <- NULL
    ##plots <- list(plot.sk,p1,p2,p3,p4,p5,p6,fc.plots)
    plots <- list(plot.sk,cluster.plots,fc.plots)
    shared.params <- list(contrast=all.ct)
    params <- list(
        list(contrast="*", fill=c(TRUE,FALSE)),
        ##list(), list(), list(), list(), list(), list(),
        list(contrast="*"),
        list(contrast="*", dir=c("up","down","both"))
    )
    ##types <- c("plotly","ggplot","ggplot","ggplot","plotly","ggplot","ggplot","ggplot")
    types <- c("plotly","ggplot","ggplot")
    widths  <- c(6,6,12)
    heights <- c(2,2,1.5)*200
    
    viz._showShiny(plots, params=params, shared.params=shared.params,
                   types=types, widths=widths, heights=heights)
    
}


##=============================================================================================
##=============================================================================================
##=============================================================================================
##plots=plotList

## just a shortcut
viz <- function(...) viz.showFigure(...)  
vizShiny <- function(fig,...) viz._showShiny(fig, ...)

widths=heights=NULL;title=subtitle=caption="Hello";tag=FALSE;plotlib=NULL;sort=FALSE;height=NULL;ncol=NULL

viz.showFigure <- function(fig, widths=NULL, heights=NULL, page.height=NULL,
                           ncol=NULL, plotlib=NULL, sort=FALSE, scale=1, 
                           title=NULL, subtitle=NULL, caption=NULL, tag=FALSE)
{
    require(cowplot)    
    is.plot <- any(class(fig[[1]]) %in% c("ggplot","plotly"))
    is.plot
    if(class(fig)=="list" && is.plot && (is.null(plotlib) || plotlib!="shiny"))
    {
        cat("viz.showFigure:: plotlist detected\n")
        plotlib0 <- intersect(c("ggplot","plotly"),class(fig[[1]]))
        cat("viz.showFigure:: plotlib0=",plotlib0,"\n")
        ##fig <- autoarrange_plots(            
        ##    fig, widths=widths, heights=heights, page.height=page.height, sort=sort,
        ##    tag=tag, ncol=ncol, plotlib=plotlib0 )
    } else {
        ## single figure
        plotlib0 <- intersect(c("ggplot","plotly"),class(fig))        
    }
    if(is.null(plotlib))
        plotlib <- plotlib0
    cat("viz.showFigure:: plotlib=",plotlib,"\n")
    
    z <- NULL
    if(plotlib=="plotly") {
        z <- viz._showPlotly(
            fig, title=title, subtitle=subtitle, caption=caption, tag=tag)
    } else if(plotlib=="shiny") {
        z <- viz._showShiny(
            fig, widths=widths, heights=heights,  ## page.height=page.height, 
            title=title, subtitle=subtitle, caption=caption)
    } else if(plotlib=="ggplot") {
        z <- viz._showGGplot(
            fig,
            ##cowplot::plot_grid(fig, scale=scale),  ## allow scaling
            title=title, subtitle=subtitle, caption=caption, tag=FALSE) ## tagging in layout
    } else {
        stop("FATAL:: unknown plotlib=",plotlib)
    }
    z
}

widths=heights=3;page.height=700;sort=tag=FALSE;plotlib="ggplot";ncol=NULL
autoarrange_plots <- function(plotlist, widths=NULL, heights=NULL, ncol=NULL, 
                              sort=FALSE, tag=FALSE, page.height=NULL, plotlib="ggplot")
{
    if(!is.null(widths)) {
        autoarrange_plots.SIZED(
            plotlist, widths=widths, heights=heights, ncol=ncol,
            page.height=page.height, sort=sort, tag=tag,
            plotlib=plotlib)
    } else {
        autoarrange_plots.UNIFORM(
            plotlist, ncol=ncol, page.height=page.height, tag=tag,
            plotlib=plotlib) 
    }        
}

n=8;ncol=12;widths=3;heights=1
autolayout.GREEDY <- function(n, ncol=12, widths=3, heights=1)
{
    M <- matrix("#",ncol,10)
    widths  <- head(rep(widths,n),n)
    heights <- head(rep(heights,n),n)    
    k=2
    k=1
    for(k in 1:n) {
        w <- widths[k]
        h <- heights[k]
        i=1;ok=FALSE
        free <- which(M=="#",arr.ind=TRUE)
        while(!ok && i<=nrow(free)) {
            ii <- free[i,1]:(free[i,1]+w-1)
            jj <- free[i,2]:(free[i,2]+h-1)
            if(max(jj) <= ncol(M) &&
               max(ii) <= nrow(M) &&
               sum(M[ii,jj]!="#")==0) {
                ok <- TRUE
                M[ii,jj] <- LETTERS[k]
            }
            i <- i+1
        }
    }
    M <- t(M[,colMeans(M=='#')<1])
    M
}

autoarrange_plots.UNIFORM <- function(plotlist, ncol=NULL, page.height=NULL,
                                      tag=FALSE, plotlib="ggplot")
{
    cat("autoarrange_plots.UNIFORM\n")
    np <- length(plotlist)
    np
    if(is.null(ncol)) {
        nrow <- floor(sqrt(np))
        nrow <- floor(sqrt(np+1))
        nrow <- pmax(floor(sqrt(np+1)*0.85),1)
        ncol <- ceiling(np/nrow)
    } else {
        nrow <- ceiling(np/ncol)
    }
    tags <- rep("",np)
    if(tag) tags <- letters[1:np]
    idx <- lapply(1:nrow, function(i) ((i-1)*ncol+1):(i*ncol))
    idx <- lapply(idx, function(i) i[which(i<=np)])
    idx
    fig <- NULL
    if(plotlib=="ggplot") {
        subrows <- lapply(idx, function(i)
            do.call(cowplot::plot_grid, list(nrow=1, plotlist=plotlist[i], labels=tags[i])))
        fig <- cowplot::plot_grid(plotlist=subrows, nrow=length(subrows))
    } else if(plotlib=="plotly") {
        subrows <- lapply(idx, function(i) do.call(subplot, plotlist[i]))
        fig <- subplot(subrows, nrows=length(subrows), margin=0.04)
    } else if(plotlib %in% c("shiny","shiny-plotly")) {
        plotnames <- names(plotlist)
        plotnames <- paste0("plot_",1:length(plotlist))
        plotnames
        outFUN <- "plotOutput"
        if(is.null(page.height)) page.height <- 720
        ht <- page.height / nrow
        if(plotlib=="shiny-plotly") outFUN <- "plotlyOutput"
        ##subrows <- lapply(idx, function(i) paste0(outFUN,"('",plotnames[i],"')"))
        subrows <- lapply(idx, function(i) paste0(outFUN,"('",plotnames[i],"',height=",ht,")"))

        ## add shinydashboard BOX
        ##subrows <- sapply(subrows, function(s) paste0("shinydashboard::box(",s,",width=NULL,title='BOX')"))
        subrows <- sapply(subrows, function(s) paste(s,collapse=","))
        subrows <- sapply(subrows, function(s)
            paste0("fillRow(",s,",height=",ht,")")) ## height important!!
        ##subrows <- sapply(subrows, function(s) paste0("fillRow(",s,")"))
        subrows <- lapply(subrows, function(s) eval(parse(text=s)))
        fig <- tagList(subrows)
    } else {
        stop("FATAL:: unknown plotlib",plotlib)
    }
    fig
}


viz._showShiny <- function(plots, params=NULL, types=NULL, shared.params=NULL,
                           tag=FALSE, widths = 3, heights = 3, 
                           title="TITLE", subtitle="Subtitle", caption="Caption"
                           ##title="", subtitle="", caption=""
                           )
{
    ## ---------------------------------------------------------------------
    ## Show a shiny page provided a list of plots. It automatically
    ## detects the type of the plot objects.
    ## ---------------------------------------------------------------------    
    
    require(shiny)
    if(class(plots)!="list") plots <- list(plots)    
    if(!is.null(params) && length(params)!=length(plots)) {
        stop("ERROR:: length of params must be equal to number of plots")
    }
    if(is.null(names(plots))) {
        names(plots) <- LETTERS[1:length(plots)]
        if(!is.null(params)) names(params) <- names(plots)
    }

    lapply(plots,class)
    if(any(sapply(plots,class)=="function") && is.null(types)) {
        stop("ERROR:: you must supply types if using functions")
    }
    
    ## guess plotlib
    if(is.null(types)) {
        types <- rep(NA,length(plots))
        i=1
        for(i in 1:length(plots)) {
            types[i] <- intersect(class(plots[[i]]),
                                  c("ggplot","plotly","scatterD3","function"))
        }
        types
    }
    
    server <- function(input, output, session) {        

        message("viz._showShiny:: server : 1")

        ## This create array of plot on the output object with names
        ## plot_1, plot_2, plot_3, etc ...
        i=1
        lapply(1:length(plots), function(i) {
            plotFUN <- NULL
            if(class(plots[[i]])[[1]]=="function") {
                np <- length(params[[i]])
                if(np>0) {
                    plotFUN <- reactive({
                        ## req(input$input_1_1)
                        inputVALUES <- list()
                        for(j in 1:np) {
                            if(params[[i]][j]=="*") {
                                ## get input from shared parameters
                                k <- match(names(params[[i]])[j],names(shared.params))
                                inputVALUES[[j]] <- eval(parse(text=paste0("input$input_0_",k)))
                            } else {
                                inputVALUES[[j]] <- eval(parse(text=paste0("input$input_",i,"_",j)))
                            }
                        }
                        do.call(plots[[i]], as.list(inputVALUES))
                    })
                } else {
                    ##plotFUN <- reactive({ plots[[1]]()})
                    ##plotFUN <- function(){ plots[[i]] }
                    plotFUN = plots[[i]]
                }
            } else {
                ##plotFUN <- reactive({ plots[[i]] })
                plotFUN <- function(){ plots[[i]] }
            }
            pclass <- types[i]
            ##pclass = "ggplot"
            r <- NULL
            if(any(pclass=="ggplot")) r <- renderPlot({plotFUN()},res=90)
            if(any(pclass=="plotly")) r <- renderPlotly(plotFUN())
            if(any(pclass=="scatterD3"))  r <- renderScatterD3(plotFUN())
            plotname <- paste0("plot_",i)
            output[[plotname]] <- r
        })
        
        message("viz._showShiny:: server : 2")
        
        ## This create array of input widgets on the input object with
        ## names input_1, input_2, input_3, etc ...
        inputWIDGETS <- vector("list",length(plots))
        for(i in 1:length(params)) {
            inputPAR <- params[[i]]
            if(length(inputPAR)) {
                ww <- list()
                for(j in 1:length(inputPAR)) {
                    if(inputPAR[[j]][1]!="*") {
                        iname <- paste0("input_",i,"_",j)
                        ititle <- names(inputPAR)[j]
                        npar <- length(inputPAR[[j]])
                        if(npar<=10) {
                            ww[[j]] <- radioButtons(iname, ititle, inputPAR[[j]], inline=TRUE)
                        } else {
                            ww[[j]] <- selectInput(iname, ititle, inputPAR[[j]])
                        }
                    } else {
                        ww[[j]] <- NULL
                    }
                }
                inputWIDGETS[[i]] <- ww
            }
        }

        ## Any shared parameters
        sharedWIDGET <- NULL
        if(!is.null(shared.params) && length(shared.params)) {
            sharedWIDGET <- list()
            for(j in 1:length(shared.params)) {
                iname <- paste0("input_0_",j)
                ititle <- names(shared.params)[j]
                npar <- length(shared.params[[j]])
                if(npar<=10) {
                    sharedWIDGET[[j]] <- radioButtons(iname, ititle, shared.params[[j]], inline=TRUE)
                } else {
                    sharedWIDGET[[j]] <- selectInput(iname, ititle, shared.params[[j]])
                }
            }
        }
        
        message("viz._showShiny:: server : 3")
        plot_output <- function(i, ht) {
            pclass <- types[[i]]
            ## pclass = "ggplot"
            plotname <- paste0("plot_",i)
            if(any(pclass=="ggplot")) {
                return(plotOutput(plotname, height=ht))
            }
            if(any(pclass=="plotly")) {
                return(plotlyOutput(plotname, height=ht))
            }
            if(any(pclass=="scatterD3")) {
                return(scatterD3Output(plotname, height=ht))
            }
            return(NULL)
        }
        
        output$plotsUI <- renderUI({
            np <- length(plots)
            wd <- rep(widths,np)[1:np]
            ht <- rep(heights,np)[1:np]
            outs <- lapply(1:np, function(i) {
                title <- names(plots)[i]
                title <- HTML(paste0("<h4>",title,"</h4>"))                
                column(width=wd[i], title, plot_output(i,ht=100*ht[i]))
            })
            outs <- do.call(tagList, outs)
            outs <- jqui_sortable(div(outs)) ## make the boxes sortable..
            outs
        })        

        output$parametersUI <- renderUI({
            if(is.null(params)) return(NULL)
            require(shinydashboardPlus)
            require(shinyWidgets)
            useShinydashboardPlus()
            i=1
            make.accordion <- function(i) {
                if(i==0 && !is.null(sharedWIDGET)) {
                    ww <- sharedWIDGET
                    acc.name <- "shared"
                    acc.id <- "accordion_0"
                } else if(i==0 && is.null(sharedWIDGET)) {
                    return(NULL)
                } else {
                    ww <- inputWIDGETS[[i]]
                    acc.name <- names(plots)[i]
                    acc.id <- paste0("accordion_",i)
                }
                if(is.null(ww) || length(ww)==0) return(NULL)
                items <- list()
                for(j in 1:length(ww)) {
                    if(!is.null(ww[[j]])) {
                        a <- paste0("input_",i,"_",j)
                        items[[j]] <- column(width=12, ww[[j]]) ## width inside sidebar
                    } else {
                        items[[j]]  <- NULL
                    }
                }
                items <- items[!sapply(items,is.null)]
                length(items)
                acc <- accordion(
                    inputId=acc.id,
                    accordionItem( title=acc.name,
                                  fluidRow(jqui_sortable(div(items))),
                                  collapsed=FALSE))
                ##acc <- accordion(inputId=a, do.call(tagList,items))
                ##acc <- do.call(accordion, c(do.call(tagList,items),inputId=a))
                ##acc <- do.call(accordion, do.call(tagList,c(items,inputId=a)))
                column(width=12, acc)
            }
            aa <- lapply(0:length(params), function(i) make.accordion(i))
            ##eval(call("accordion", list(tags, inputId="accordion1")))
            fluidRow(jqui_sortable(div(do.call(tagList,aa))))
        })      
        outputOptions(output, "plotsUI", suspendWhenHidden=FALSE) ## important!!!
        outputOptions(output, "parametersUI", suspendWhenHidden=FALSE) ## important!!!        
        observeEvent(input$done,{stopApp()})

        message("viz._showShiny:: server : done")

    }

    ##--------------------------------------------------------------------
    ## UI
    ##--------------------------------------------------------------------
    require(miniUI)

    full.ui <- fluidPage(
        div(titlePanel(title), style="color: darkblue;"),
        HTML(subtitle,"<br><br>"),
        uiOutput("plotsUI"),
        div(HTML("<br>",caption), style="text-align:right; font-size:small;")
    )    

    sidebar.ui <- fluidPage(
        div(titlePanel(title), style="color: darkblue;"),
        HTML(subtitle,"<br><br>"),
        jqui_sortable(
            sidebarLayout(
                sidebarPanel(width=3, NULL, uiOutput("parametersUI")),
                mainPanel(uiOutput("plotsUI"),width=9)
            )
        ),
        div(HTML("<br>",caption), style="text-align:right; font-size:small;")
    )    

    gadget.ui <- miniPage(
        gadgetTitleBar(title, left=NULL),
        miniContentPanel(sidebar.ui)
    )    

    ##runGadget(mini.ui, server)
    if(!is.null(params)) {
        runGadget(sidebar.ui, server)
        ##shinyApp(sidebar.ui, server)
    } else {
        runGadget(full.ui, server)
        ##shinyApp(full.ui, server)
    }
}

viz._showGGplot <- function(fig, title="", subtitle="", caption="", tag=FALSE)
{
    require(patchwork)
    class(fig)
    if(class(fig)[1]=="list" && "ggplot" %in% class(fig[[1]]) ) {
        fig <- patchwork::wrap_plots(fig)
    }
    tag1 <- NULL
    if(tag) tag1="a"
    fig <- fig + plot_annotation(
                     title = title,
                     subtitle = subtitle,
                     caption = caption,
                     tag_levels = tag1,
                     theme = theme(
                         plot.title = element_text(
                             size=18, face='bold', color="royalblue4")
                     ))
    fig <- fig & theme(plot.tag = element_text(face='bold'))    
    if(0) {
        logo <- magick::image_read("~/Playground/logo/bigomics-logo-blue.png")
        grid::grid.raster(logo, x = 0.99, y = 0.99, just = c('right', 'top'),
                          width = unit(1.2, 'cm'))
    }
    fig
}

viz._showPlotly <- function(fig, title="", subtitle="", caption="", tag=FALSE)
{
    require(plotly)
    require(patchwork)    

    if(class(fig)=="list" &&
       ( "ggplot" %in% class(fig[[1]]) ||
         "plotly" %in% class(fig[[1]]) ))
    {
        nr <- floor(sqrt(length(fig)+1))
        fig <- subplot(fig, nrows=nr, margin=0.06)
    }

    header <- plotly_empty()
    caption <- plotly_empty()        

    if(!is.null(subtitle)) {
        header <- header %>%
            add_annotations(
                x = 0.0, y = 1,
                yshift = -3, ## xshift = 12, 
                ##align= "right",
                xanchor = "left",
                xref = "paper", yref = "paper",
                text = subtitle,
                showarrow = FALSE
            )
    }

    if(!is.null(caption)) {    
        caption <- caption %>% 
            add_annotations(
                x= 1.0, y=0.0, align= "right",
                yshift = -20, ## xshift = 80, 
                xref = "paper", yref = "paper",
                text = caption,
                showarrow = FALSE
            )
    }

    if(is.null(title) && is.null(subtitle) && is.null(caption)) {
        plt <- fig
    } else {
        ht <- c(0.08,0.90,0.02)
        plt <- subplot(header, fig, caption, nrows=3, heights=ht, margin=0.0) %>%
            layout(showlegend = FALSE)
    }

    if(!is.null(title)) {
        btitle <- paste0("<b>",title,"</b>")
        plt <- plt %>%
            layout(title = list(text=btitle, x=0.001, y=0.99, font=list(size=24))) 
    }
    plt <- plt %>%
        config(toImageButtonOptions = list(format='svg', height=900, width=1500))
    plt

}


##====================================================================================
##======================== END OF FILE ===============================================
##====================================================================================




