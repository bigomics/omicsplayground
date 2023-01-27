##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

enrichment_plot_geneplot_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "An expression barplot per sample group for the gene that is selected from the genes Table <code>II</code>. Samples can be ungrouped in the barplot by selecting <code>ungroup samples</code> from the plot <i>Settings</i>."

  options <- shiny::tagList(
    withTooltip( shiny::checkboxInput(ns('gs_ungroup2'),'ungroup samples',FALSE),
                 "Ungroup samples in the plot", placement="top",
                 options = list(container = "body")))

  PlotModuleUI(
    ns("plot"),
    title = "Expression geneplot",
    label = "c",
    info.text = info_text,
    options = options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_geneplot_server <- function(id,
                                            inputData,
                                            gs_contrast,
                                            gene_selected,
                                            subplot.MAR,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    subplot_geneplot.RENDER <- shiny::reactive({

      par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0,0.4) )
      par(mar=subplot.MAR)

      ngs <- inputData()
      shiny::req(ngs)

      comp0 = colnames(ngs$model.parameters$contr.matrix)[1]
      comp0 = gs_contrast()

      has.design <- !is.null(ngs$model.parameters$design)
      collapse.others <- ifelse(has.design, FALSE, TRUE)
      ##collapse.others=TRUE

      sel  = gene_selected()
      if(is.null(sel) || is.na(sel) || length(sel)==0) {
        frame()
      } else {
        probe = sel$probe
        gene = sel$gene
        if(length(probe)>1) {
          probe <- grep("\\[gx|\\[mrna",probe,value=TRUE)
        }
        ngrp <- length(unique(ngs$samples$group))
        grouped=TRUE
        grouped <- !input$gs_ungroup2
        srt <- ifelse(!grouped || ngrp>4, 30, 0)
        if(!grouped && ncol(ngs$X) > 15) srt <- 60
        pgx.plotExpression(
          ngs, probe, comp=comp0, logscale=TRUE, level="gene",
          collapse.others=collapse.others, grouped=grouped,
          srt=srt, main="")
        title(gene, cex.main=0.9)
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = subplot_geneplot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(78,100),
      add.watermark = watermark
    )
  })
}
