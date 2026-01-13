##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
singlecell_plot_cytoplot_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width,
  parent
) {
  ns <- shiny::NS(id)

  cyto.opts <- shiny::tagList(
    withTooltip(shiny::selectInput(parent("cytovar1"), label = "x-axis:", choices = NULL, multiple = FALSE),
      "Select your prefered gene on the x-axis.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(parent("cytovar2"), label = "y-axis:", choices = NULL, multiple = FALSE),
      "Choose your prefered gene on the y-axis.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::sliderInput(parent("nbins"), label = "nbins:", min = 0, max = 30, value = 10, step = 5),
      "Select the maximum number of bins for histogram distribution.",
      placement = "top", options = list(container = "body")
    ),
  )

  PlotModuleUI(
    id = ns("plotmodule"),
    plotlib = "plotly",
    info.text = info.text,
    title = title,
    caption = caption,
    options = cyto.opts,
    download.fmt = c("png", "pdf", "svg"),
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
                                            pgx,
                                            pfGetClusterPositions,
                                            samplefilter, # input$samplefilter
                                            cytovar1,
                                            cytovar2,
                                            nbins,
                                            selectSamplesFromSelectedLevels,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({})

    ##    cyto.plotFUNC <- shiny::reactive({
    cyto.plotFUNC <- function() {
      shiny::req(pgx$X)
      shiny::req(cytovar1(), cytovar2(), cytovar1() != "", cytovar2() != "")
      cytovar1 <- cytovar1()
      cytovar2 <- cytovar2()
      samplefilter <- samplefilter()
      nbins <- nbins()
      kk <- playbase::selectSamplesFromSelectedLevels(pgx$Y, samplefilter)

      playbase::plotlyCytoplot(pgx,
        gene1 = cytovar1,
        gene2 = cytovar2,
        nbinsx = nbins,
        nbinsy = nbins,
        marker.size = 8,
        samples = kk,
        lab.unit = "  (log2CPM)",
        contour.coloring = "none"
      )

      ##      plotly::plot_ly(z = ~volcano, type = "contour")
    } ## )

    PlotModuleServer(
      "plotmodule",
      func = cyto.plotFUNC,
      ##      func2 = plotly_modal.RENDER,
      plotlib = "plotly",
      res = c(90, 130), ## resolution of plots
      pdf.width = 8,
      pdf.height = 8,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
