##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics Sagl. All rights reserved.
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
compare_plot_compare1_ui <- function(id,
                                     height,
                                     width) {
  ns <- shiny::NS(id)
  info_text <- "<b>Expression scatter plot.</b> Scatter plots of gene expression."

  PlotModuleUI(ns("plot"),
    title = "Dataset 1",
    plotlib = "base",
    label = "a",
    info.text = NULL,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
compare_plot_compare1_server <- function(id,
                                         inputData,
                                         input.contrast1,
                                         hilightgenes,
                                         createPlot,
                                         plottype,
                                         dataset2,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      inputData()
    })

    scatter1.RENDER <- shiny::reactive({
      ngs1 <- plot_data()
      ngs2 <- dataset2()
      all.ct <- names(ngs1$gx.meta$meta)
      ct1 <- input.contrast1()
      shiny::req(ct1)
      if (!all(ct1 %in% all.ct)) {
        return(NULL)
      }
      higenes <- hilightgenes()
      cex.lab <- 1.0
      ntop <- 9999
      type <- plottype()

      if (length(higenes) <= 3) cex.lab <- 1.3
      createPlot(ngs1, ngs1, ngs2, ct1, type, cex.lab, higenes, ntop)
    })

    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = scatter1.RENDER,
      csvFunc = plot_data,
      res = c(90, 110), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
