##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Single cell plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
singlecell_plot_markersplot_ui <- function(id,
                                          label='',
                                          height,
                                          width,
                                          parent){
  ns <- shiny::NS(id)

  markersplot.opts = shiny::tagList(
    withTooltip(shiny::selectInput(parent("mrk_level"),"Level:", choices=c("gene","geneset")),
                "Specify the level of the marker analysis: gene or gene set level.",
                placement="top", options = list(container = "body")),
    withTooltip(shiny::selectInput(parent("mrk_features"),"Feature set:", choices=NULL,
                                   multiple=FALSE),
                "Select a particular functional group for the analysis.",
                placement="top", options = list(container = "body")),
    withTooltip(shiny::textInput(parent("mrk_search"),"Filter:"),
                "Filter markers by a specific keywords.",
                placement="top", options = list(container = "body")),
    withTooltip(shiny::radioButtons(parent("mrk_sortby"),"Sort by:",
                                    choices=c("intensity","name"), inline=TRUE),
                "Sort by name or intensity.", placement="top",
                options = list(container = "body"))
  )

  markersplot_info = "The Markers section produces for the top marker genes, a t-SNE with samples colored in red when the gene is overexpressed in corresponding samples. The top genes (N=36) with the highest standard deviation are plotted. <p>In the plotting options, users can also restrict the marker analysis by selecting a particular functional group in which genes are divided into 89 groups, such as chemokines, transcription factors, genes involved in immune checkpoint inhibition, and so on."


  PlotModuleUI(ns("plot"),
               label = label,
               info.text = markersplot_info,
               options = markersplot.opts,
               download.fmt=c("png","pdf","csv"),
               height = height,
               width = width)
}

#' Single cell plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @export
singlecell_plot_markersplot_server <- function(id,
                                               inputData,
                                               pfGetClusterPositions,
                                               mrk_level,
                                               mrk_features,
                                               mrk_search,
                                               mrk_sortby,
                                               watermark = FALSE){
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    markers.plotFUNC <- shiny::reactive({
      ##if(!input$tsne.all) return(NULL)

      ngs <- inputData()
      shiny::req(ngs)

      mrk_level <- mrk_level()
      mrk_features <- mrk_features()
      mrk_search <- mrk_search()
      mrk_sortby <- mrk_sortby()

      clust.pos <- pfGetClusterPositions()
      if(is.null(clust.pos)) return(NULL)
      ##pos <- ngs$tsne2d
      pos <- clust.pos

      ##markers <- ngs$families[["CD family"]]
      if(is.null(mrk_features)) return(NULL)
      if(mrk_features=="") return(NULL)

      term = ""
      if(mrk_level=="gene") {
        markers <- ngs$families[["Transcription factors (ChEA)"]]
        if(mrk_search!="") {
          term = mrk_search
          jj <- grep(term, ngs$genes$gene_name, ignore.case=TRUE )
          markers <- ngs$genes$gene_name[jj]
          term = paste("filter:",term)
        } else if(mrk_features %in% names(ngs$families)) {
          markers <- ngs$families[[mrk_features]]
          term = mrk_features
        } else {
          markers <- ngs$genes$gene_name
        }
        ##markers <- intersect(markers, rownames(ngs$X))
        markers <- intersect(toupper(markers),toupper(ngs$genes$gene_name))
        jj <- match(markers,toupper(ngs$genes$gene_name))
        pmarkers <- intersect(rownames(ngs$genes)[jj],rownames(ngs$X))
        gx <- ngs$X[pmarkers,rownames(pos),drop=FALSE]

      } else if(mrk_level=="geneset") {
        ##markers <- ngs$families[["Immune checkpoint (custom)"]]
        markers <- COLLECTIONS[[1]]
        if(is.null(mrk_features)) return(NULL)
        ft <- mrk_features
        if(mrk_search=="" && ft %in% names(COLLECTIONS)) {
          markers <- COLLECTIONS[[mrk_features]]
          markers <- intersect(markers, rownames(ngs$gsetX))
          term = mrk_features
        } else if(mrk_search!="") {
          term = mrk_search
          jj <- grep(term, rownames(ngs$gsetX), ignore.case=TRUE )
          markers <- rownames(ngs$gsetX)[jj]
          term = paste("filter:",term)
        } else {
          markers <- rownames(ngs$gsetX)
        }
        gx <- ngs$gsetX[markers,rownames(pos),drop=FALSE]
      } else {
        cat("fatal error")
        return(NULL)
      }

      if(!"group" %in% names(ngs$model.parameters)) {
        stop("[markers.plotFUNC] FATAL: no group in model.parameters")
      }

      ## prioritize gene with large variance (groupwise)
      ##grp <- as.character(ngs$samples[rownames(pos),"group"])
      grp <- ngs$model.parameters$group[rownames(pos)]
      zx <- t(apply(gx,1,function(x) tapply(x,grp,mean)))
      gx <- gx[order(-apply(zx,1,sd)),,drop=FALSE]
      gx <- gx - min(gx,na.rm=TRUE) + 0.01 ## subtract background??
      rownames(gx) = sub(".*:","",rownames(gx))

      ##gx <- tanh(gx/sd(gx) ) ## softmax
      cex1 = 1.0
      cex1 <- 0.8*c(2.2,1.1,0.6,0.3)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]
      klrpal <- colorRampPalette(c("grey90", "grey80", "grey70", "grey60","red4", "red3"))(16)
      klrpal = colorRampPalette(c("grey90", "grey60", "red3"))(16)
      klrpal = paste0(gplots::col2hex(klrpal),"66")

      NP=25
      if(mrk_level=="gene") NP=36
      top.gx = head(gx,NP)  ## match number of plot below!
      if(mrk_sortby=="name") {
        top.gx = top.gx[order(rownames(top.gx)),,drop=FALSE]
      } else {
        top.gx = top.gx[order(-rowMeans(top.gx)),,drop=FALSE]
      }
      top.gx = pmax(top.gx,0)
      ##top.gx <- tanh(top.gx/mean(top.gx))

      plevel="gene"
      plevel <- mrk_level

      par(mfrow=c(1,1)*sqrt(NP), mar=c(0,0.2,0.5,0.2)*0.6, oma=c(1,1,1,1)*0.5)
      par(mfrow=c(1,1)*sqrt(NP), mar=c(0,0.2,0.5,0.2)*0.6, oma=c(1,1,6,1)*0.5)
      ##par(mfrow=c(6,6), mar=c(0,0.2,0.5,0.2), oma=c(1,1,1,1)*0.5)
      i=1
      for(i in 1:min(NP,nrow(top.gx))) {

        colvar = pmax(top.gx[i,],0)
        colvar = 1+round(15*(colvar/(0.7*max(colvar)+0.3*max(top.gx))))
        klr0 = klrpal[colvar]

        ii <- order(colvar)
        ##ii <- sample(nrow(pos))
        base::plot( pos[ii,], pch=19, cex=cex1, col=klr0[ii],
                    xlim=1.1*range(pos[,1]), ylim=1.1*range(pos[,2]),
                    fg = gray(0.8), bty = "o",
                    xaxt='n', yaxt='n', xlab="tSNE1", ylab="tSNE2")

        if(plevel=="gene") {
          gene <- sub(".*:","",rownames(top.gx)[i])
          ##title( gene, cex.main=1.0, line=0.3, col="grey40", font.main=1)
          legend( "topleft", legend=gene, bg="#AAAAAA88",
                  cex=1.2, text.font=1, y.intersp=0.8,
                  bty="n", inset=c(-0.05,-0.0) )
        } else {
          gset <- sub(".*:","",rownames(top.gx)[i])
          gset1 <- breakstring(substring(gset,1,80),24,force=TRUE)
          gset1 <- tolower(gset1)
          ##title( gset1, cex.main=0.9, line=0.4, col="grey40", font.main=1)
          legend( "topleft", legend=gset1, cex=0.95, bg="#AAAAAA88",
                  text.font=2, y.intersp=0.8, bty="n",
                  inset=c(-0.05,-0.0) )
        }
      }
      mtext(term, outer=TRUE, cex=1.0, line=0.6)

    })

    PlotModuleServer(
      "plot",
      func = markers.plotFUNC,
      res = c(85,90),
      pdf.width = 10, pdf.height = 10,
      add.watermark = watermark
    )
  })## end of moduleServer
}