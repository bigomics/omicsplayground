##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
expression_plot_volcano_ui <- function(id,
                                       label='',
                                       height=c(600, 800)) {
  ns <- shiny::NS(id)
  info_text = "A volcano plot of genes for the selected comparison under the <code>Contrast</code> settings plotting fold-change versus significance on the x and y axes, respectively."

  PlotModuleUI(ns("pltmod"),
               title = "Volcano plot",
               label = label,
               plotlib = "plotly",
               ##outputFunc = plotly::plotlyOutput,
               ##outputFunc2 = plotly::plotlyOutput,
               info.text = info_text,
               options = NULL,
               download.fmt=c("png","pdf","csv"),
               width = c("auto","100%"),
               height = height

}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#' @param pgx
#' @param pgx_fdr
#' @param pgx_contrast
#' @param pgx_lfc
#' @param pgx_features
#' @param res
#' @param watermark
#'
#'
#' @export
expression_plot_volcano_server <- function(id,
                                           pgx,
                                           pgx_fdr = 0.1,
                                           pgx_contrast,
                                           pgx_lfc = 1.0,
                                           pgx_features,
                                           res,
                                           watermark = FALSE
                                           )
{
  moduleServer(id, function(input, output, session, watermark) {

    dbg("[plots_volcano.PLOTLY] reacted")

    #calculate required inputs for plotting
    serverSideComputation <- function(pgx,
                                      pgx_fdr,
                                      pgx_contrast,
                                      pgx_lfc,
                                      pgx_features,
                                      res,
                                      fam.genes){
      shiny::req(pgx)
      comp1 = gx_contrast
      alertDataLoaded(session,pgx)
      fdr = as.numeric(pgx_fdr)
      res = fullDiffExprTable()
      lfc = as.numeric(pgx_lfc)
      fam.genes = res$gene_name
      ##fam.genes = unique(unlist(pgx$families[input$gx_features]))

      if(is.null(res)) return(NULL)
      if(length(comp1)==0) return(NULL)
      if(is.null(pgx_features)) return(NULL)
      if(pgx_features!="<all>") {
        ##gset <- GSETS[input$gx_features]
        gset <- getGSETS(pgx_features)
        fam.genes = unique(unlist(gset))
      }

      jj <- match(toupper(fam.genes),toupper(res$gene_name))
      sel.genes <- res$gene_name[setdiff(jj,NA)]

      fc.genes = as.character(res[,grep("^gene$|gene_name",colnames(res))])
      qval = res[,grep("adj.P.Val|meta.q|qval|padj",colnames(res))[1]]
      qval = pmax(qval, 1e-20)
      x = res[,grep("logFC|meta.fx|fc",colnames(res))[1]]
      y <- -log10(qval + 1e-12)

      sig.genes = fc.genes[which(qval <= fdr & abs(x) > lfc )]
      sel.genes = intersect(sig.genes, sel.genes)
      scaled.x <- scale(x,center=FALSE)
      scaled.y <- scale(y,center=FALSE)
      impt <- function(g) {
        j = match(g, fc.genes)
        x1 = scaled.x[j]
        y1 = scaled.y[j]
        x = sign(x1)*(x1**2 + 0.25*y1**2)
        names(x)=g
        x
      }

      sel1 = genetable$rows_selected()
      df1 = filteredDiffExprTable()

      sel2 = gsettable$rows_selected()
      df2 <- gx_related_genesets()

      lab.cex = 1
      gene.selected <- !is.null(sel1) && !is.null(df1)
      gset.selected <- !is.null(sel2) && !is.null(df2)
      if(gene.selected && !gset.selected) {
        lab.genes = rownames(df1)[sel1]
        sel.genes = lab.genes
        lab.cex = 1.3
      } else if(gene.selected && gset.selected) {
        gs <- rownames(df2)[sel2]
        dbg("[plots_volcano.PLOTLY] gs = ",gs)
        ##gset <- GSETS[[gs]]
        gset <- unlist(getGSETS(gs))
        sel.genes = intersect(sel.genes, gset)
        lab.genes = c( head(sel.genes[order(impt(sel.genes))],10),
                       head(sel.genes[order(-impt(sel.genes))],10) )
        lab.cex = 1
      } else {
        lab.genes = c( head(sel.genes[order(impt(sel.genes))],10),
                       head(sel.genes[order(-impt(sel.genes))],10) )
        lab.cex = 1
      }
      xlim = c(-1,1)*max(abs(x),na.rm=TRUE)
      ylim = c(0, max(12, 1.1*max(-log10(qval),na.rm=TRUE)))

      ## par(mfrow=c(1,1), mar=c(4,3,1,1.5), mgp=c(2,0.8,0), oma=c(0,0,0,0))

      return(list(x,y,fc.genes, sel.genes,lab.genes,lab.cex,fdr,lfc))

    }

    #reactive function listening for changes in input
    plot_data <- shiny::reactive({
      serverSideComputation(pgx,
                            pgx_fdr,
                            pgx_contrast,
                            pgx_lfc,
                            pgx_features,
                            res,
                            fam.genes)
      })


    plotly.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      df <- pd
      par(mfrow=c(1,1), mar=c(4,3,1,1.5), mgp=c(2,0.8,0), oma=c(0,0,0,0))
      plt <- plotlyVolcano(
        x=x, y=y, names=fc.genes,
        source = "plot1", marker.type = "scattergl",
        highlight = sel.genes,
        label = lab.genes, label.cex = lab.cex,
        group.names = c("group1","group0"),
        ##xlim=xlim, ylim=ylim, ## hi.col="#222222",
        ##use.fdr=TRUE,
        psig = fdr, lfc = lfc,
        xlab = "effect size (log2FC)",
        ylab = "significance (-log10q)",
        marker.size = 4,
        displayModeBar = FALSE,
        showlegend = FALSE
        ) %>% plotly::layout( margin = list(b=65))
      plt
    }

    modal_plotly.RENDER <- function() {
      plotly.RENDER() %>%
        plotly::layout(
          ## showlegend = TRUE,
          font = list(
            size = 16
          )
        )
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data,   ##  *** downloadable data as CSV
      res = c(80,170),                ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
    }## end of moduleServer

    dbg("[plots_volcano.PLOTLY] done!")
}