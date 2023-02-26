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
singlecell_plot_cytoplot_ui <- function(id,
                                        label = "",
                                        height,
                                        width,
                                        parent) {
  ns <- shiny::NS(id)

  cyto.opts <- shiny::tagList(
    withTooltip(shiny::selectInput(parent("cytovar1"), label = "x-axis:", choices = NULL, multiple = FALSE),
      "Select your prefered gene on the x-axis.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(parent("cytovar2"), label = "y-axis:", choices = NULL, multiple = FALSE),
      "Choose your prefered gene on the y-axis.",
      placement = "top", options = list(container = "body")
    )
  )

  cytoModule_info <- "For each combination of gene pairs, the platform can generate a cytometry-like plot of samples under the Cytoplot tab. The aim of this feature is to observe the distribution of samples in relation to the selected gene pairs. For instance, when applied to single-cell sequencing data from immunological cells, it can mimic flow cytometry analysis and distinguish T helper cells from the other T cells by selecting the CD4 and CD8 gene combination."


  PlotModuleUI(ns("plot"),
    label = label,
    info.text = cytoModule_info,
    options = cyto.opts,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width
  )
}

#' Single cell plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @export
singlecell_plot_cytoplot_server <- function(id,
                                            inputData,
                                            pfGetClusterPositions,
                                            samplefilter, # input$samplefilter
                                            cytovar1,
                                            cytovar2,
                                            selectSamplesFromSelectedLevels,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    cyto.plotFUNC <- shiny::reactive({
      ## if(!input$tsne.all) return(NULL)

      ngs <- inputData()
      ## if(is.null(ngs)) return(NULL)
      shiny::req(ngs)

      cytovar1 <- cytovar1()
      cytovar2 <- cytovar2()
      samplefilter <- samplefilter()

      if (is.null(cytovar1)) {
        return(NULL)
      }
      if (is.null(cytovar2)) {
        return(NULL)
      }
      if (cytovar1 == "") {
        return(NULL)
      }
      if (cytovar2 == "") {
        return(NULL)
      }

      kk <- selectSamplesFromSelectedLevels(ngs$Y, samplefilter)
      gene1 <- cytovar1
      gene2 <- cytovar2
      ## if(gene1 == gene2) return(NULL)
      par(mfrow = c(1, 1), mar = c(10, 5, 4, 1))
      pgx.cytoPlot(ngs, gene1, gene2,
        samples = kk, cex = 0.8,
        col = "grey60", cex.names = 1, lab.unit = "(log2CPM)"
      )
    })

    PlotModuleServer(
      "plot",
      func = cyto.plotFUNC,
      res = c(80, 80),
      pdf.width = 6, pdf.height = 8,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
